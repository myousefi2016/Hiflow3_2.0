// Copyright (C) 2011-2017 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for more
// details.
//
// You should have received a copy of the European Union Public Licence (EUPL) v1.2
// along with HiFlow3.  If not, see <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

/// @author Dimitar Lukarski

#include "lpreconditioner_mc.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "lvector.h"
#include "lvector_cpu.h"

#include "lmatrix.h"
#include "lmatrix_csr_cpu.h"

#include "lmp_log.h"

#include <math.h>

namespace hiflow
{
    namespace la
    {

        const int DEBUG_LEVEL = 1;

        template <typename ValueType>
        lPreconditioner_MultiColoring<ValueType>::lPreconditioner_MultiColoring ( )
        {
            this->Operator_ = NULL;
            this->Operator_perm_ = NULL;

            this->x_ = NULL;
            this->x_const_ = NULL;
            this->rhs_ = NULL;

            this->LU_ = NULL;
            this->LU_cpu_ = NULL;

            this->aux_mat_analysis_cpu_ = NULL;

            this->ncolors_ = 0;

            this->flag_mc_ = true;
            this->flag_dropoff_mc_ = true;
            this->flag_ls_ = false;

            this->omega_relax_ = 1.0;
        }

        template <typename ValueType>
        lPreconditioner_MultiColoring<ValueType>::~lPreconditioner_MultiColoring ( )
        {
            this->Clear ( );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::Init ( )
        {
            LOG_ERROR ( "lPreconditioner_MultiColoring<ValueType>::Init() - wrong initialization call" );
            exit ( -1 );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::SetupPermVectors ( hiflow::la::lVector<ValueType> *x,
                                                                          hiflow::la::lVector<ValueType> *rhs )
        {
            assert ( x != NULL );
            assert ( rhs != NULL );

            this->x_ = x;
            this->rhs_ = rhs;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::SetupVector ( const hiflow::la::lVector<ValueType> *x )
        {
            assert ( x != NULL );

            this->x_const_ = x;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::SetRelax ( const ValueType omega_relax )
        {
            assert ( omega_relax > 0.0 );
            assert ( omega_relax < 2.0 );

            this->omega_relax_ = omega_relax;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::Clear ( )
        {
            if ( this->ncolors_ > 0 )
            {

                for ( int i = 0; i<this->ncolors_; ++i )
                {
                    // diag vectors
                    delete this->Dv_[i];
                    delete this->iDv_[i];

                    // output chunks
                    delete this->output_chunks_[i];
                }

                for ( int i = 0; i < ( this->ncolors_ )*( this->ncolors_ - 1 ) / 2; ++i )
                {
                    // L, R matrices
                    delete this->L_[i];
                    delete this->R_[i];
                }

                // free the permut vector and the color array
                delete [] this->color_sizes_;
                delete [] this->permut_index_;
            }

            this->Operator_ = NULL;
            this->Operator_perm_ = NULL;

            this->x_ = NULL;
            this->x_const_ = NULL;
            this->rhs_ = NULL;

            this->LU_ = NULL;
            this->LU_cpu_ = NULL;

            this->aux_mat_analysis_cpu_ = NULL;

            this->ncolors_ = 0;

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::SetupPermOperator ( hiflow::la::lMatrix<ValueType> &op )
        {
            assert ( &op != NULL );
            this->Operator_perm_ = &op;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::Build ( )
        {
            this->Analyse ( );

            // Multi-Coloring Permutation
            if ( this->flag_mc_ )
            {
                // on local level
                this->permute ( );

                // on global level
                // one should use get_permut() to permut the matrix
            }

            //  this->Build_LU();

            this->factorize ( );
            //  this->LU_->WriteFile("LU_original.mtx");

            if ( this->flag_dropoff_mc_ )
                this->dropoffs ( );

            // Level-Scheduling Permutation
            if ( this->flag_ls_ )
            {

                // Level-Scheduling / Analyse
                this->levelscheduling ( );

                // on local level
                this->permute_all ( );

                // on global level
                // one should use get_permut() to permut the matrix
                // drop again (should drop 0 elements)
                this->dropoffs ( );

            }

            this->Prepare_LU ( );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::Analyse ( )
        {
            assert ( this->Operator_ != NULL );

            this->build_aux_matrix ( );

            if ( this->flag_mc_ )
            {
                this->multicoloring ( );
            }

            this->allocate_LU ( );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::Prepare_LU ( )
        {
            this->extract_mat ( );
            this->allocate_vec ( );
            this->convert_Dmat2vec ( );
            this->scale_D ( );
            this->check_for_dropoffs ( );

            //    this->Operator_->WriteFile("op.mtx");
            //    this->LU_->WriteFile("LU.mtx");

            this->cleanup ( );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::Build_LU ( )
        {
            this->factorize ( );

            if ( this->flag_dropoff_mc_ )
                this->dropoffs ( );

            // Level-Scheduling Permutation
            if ( this->flag_ls_ )
            {

                // Level-Scheduling / Analyse
                this->levelscheduling ( );

                // on local level
                this->permute_all ( );

                // on global level
                // one should use get_permut() to permut the matrix

                // drop again (should drop 0 elements)
                this->dropoffs ( );
            }

            this->extract_mat ( );
            this->allocate_vec ( );
            this->convert_Dmat2vec ( );
            this->scale_D ( );
            this->check_for_dropoffs ( );
            this->cleanup ( );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::levelscheduling ( )
        {
            if ( this->flag_mc_ )
            {
                // delete the old colors
                this->ncolors_ = 0;
                delete [] this->color_sizes_;
                delete [] this->permut_index_;
            }

            LU_cpu_->Levelscheduling ( this->ncolors_,
                                       &( this->color_sizes_ ),
                                       &( this->permut_index_ ) );

            LOG_DEBUG ( 1, "lprecond mc number of levels: " << this->ncolors_ );

            for ( int k = 0; k<this->ncolors_; ++k )
                LOG_DEBUG ( 2, "lprecond mc level=" << k << " ;number of elements in this level= " << this->color_sizes_[k] );

            char data_info[255];
            sprintf ( data_info, " / levels = %d", this->ncolors_ );
            this->precond_name_.append ( data_info );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::multicoloring ( )
        {

            aux_mat_analysis_cpu_->Multicoloring ( this->ncolors_,
                                                   &( this->color_sizes_ ),
                                                   &( this->permut_index_ ) );

            LOG_DEBUG ( 1, "lprecond mc number of colours: " << this->ncolors_ );

            for ( int k = 0; k<this->ncolors_; ++k )
                LOG_DEBUG ( 2, "lprecond mc color=" << k << " ;number of elements in this color= " << this->color_sizes_[k] );

            char data_info[255];
            sprintf ( data_info, " / colors = %d", this->ncolors_ );
            this->precond_name_.append ( data_info );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::permute_all ( )
        {
            // Permute the operator

            assert ( this->Operator_perm_ != NULL );

            CPU_CSR_lMatrix<ValueType> *mat_cpu = new CPUsimple_CSR_lMatrix<ValueType>;

            // Permute the Operator
            mat_cpu->CloneFrom ( *this->Operator_perm_ );
            mat_cpu->Reorder ( this->permut_index_ );

            this->Operator_perm_->Clear ( );
            this->Operator_perm_->CloneFrom ( *mat_cpu );

            // Permute the LU and LU_cpu
            mat_cpu->CloneFrom ( *this->LU_cpu_ );
            mat_cpu->Reorder ( this->permut_index_ );

            this->LU_cpu_->Clear ( );
            this->LU_cpu_->CloneFrom ( *mat_cpu );

            //  this->LU_cpu_->WriteFile("LU_ls.mtx");

            this->LU_->Clear ( );
            this->LU_->CloneFrom ( *mat_cpu );

            delete mat_cpu;

            // make x and rhs a cpu copy
            lVector<ValueType> *x_cpu = new CPUsimple_lVector<ValueType>;
            lVector<ValueType> *rhs_cpu = new CPUsimple_lVector<ValueType>;

            x_cpu ->CloneFrom ( *this->x_ );
            rhs_cpu->CloneFrom ( *this->rhs_ );

            // permute the x and rhs
            x_cpu->Reorder ( this->permut_index_ );
            rhs_cpu->Reorder ( this->permut_index_ );

            // copy back
            this->x_ ->CloneFrom ( *x_cpu );
            this->rhs_->CloneFrom ( *rhs_cpu );

            delete x_cpu;
            delete rhs_cpu;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::permute ( )
        {
            // Permute the operator

            assert ( this->Operator_perm_ != NULL );

            lMatrix<ValueType> *mat_cpu = new CPUsimple_CSR_lMatrix<ValueType>;
            mat_cpu->CloneFrom ( *this->Operator_perm_ );

            mat_cpu->Reorder ( this->permut_index_ );

            this->Operator_perm_->Clear ( );
            this->Operator_perm_->CloneFrom ( *mat_cpu );
            delete mat_cpu;

            // again *LU_ = *Operator_
            LU_cpu_->CloneFrom ( *this->Operator_perm_ );

            LU_->Clear ( );
            LU_->CloneFrom ( *LU_cpu_ );

            // make x and rhs a cpu copy
            lVector<ValueType> *x_cpu = new CPUsimple_lVector<ValueType>;
            lVector<ValueType> *rhs_cpu = new CPUsimple_lVector<ValueType>;

            x_cpu ->CloneFrom ( *this->x_ );
            rhs_cpu->CloneFrom ( *this->rhs_ );

            // permute the x and rhs
            x_cpu->Reorder ( this->permut_index_ );
            rhs_cpu->Reorder ( this->permut_index_ );

            // copy back
            this->x_ ->CloneFrom ( *x_cpu );
            this->rhs_->CloneFrom ( *rhs_cpu );

            delete x_cpu;
            delete rhs_cpu;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::PermuteBack ( hiflow::la::lMatrix<ValueType> *op ) const
        {
            assert ( op != NULL );
            assert ( op->get_num_row ( ) == this->Operator_->get_num_row ( ) );
            assert ( op->get_num_col ( ) == this->Operator_->get_num_col ( ) );

            lMatrix<ValueType> *mat_cpu = new CPUsimple_CSR_lMatrix<ValueType>;
            mat_cpu->CloneFrom ( *op );

            int *ind = new int[op->get_num_row ( )];
            assert ( ind != NULL );

            // make back permutation
            for ( int i = 0; i < op->get_num_row ( ); ++i )
                ind[ this->permut_index_[i] ] = i;

            mat_cpu->Reorder ( ind );

            delete [] ind;

            op->Clear ( );
            op->CloneFrom ( *mat_cpu );

            delete mat_cpu;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::PermuteBack ( hiflow::la::lVector<ValueType> *vec ) const
        {
            assert ( vec != NULL );
            assert ( vec->get_size ( ) == this->x_->get_size ( ) );

            lVector<ValueType> *vec_cpu = new CPUsimple_lVector<ValueType>;
            vec_cpu->CloneFrom ( *vec );

            int *ind = new int[vec->get_size ( )];
            assert ( ind != NULL );

            // make back permutation
            for ( int i = 0; i < vec->get_size ( ); ++i )
                ind[ this->permut_index_[i] ] = i;

            //  vec_cpu->Reorder(this->permut_index_);
            vec_cpu->Reorder ( ind );

            delete [] ind;

            vec->Clear ( );
            vec->CloneFrom ( *vec_cpu );

            delete vec_cpu;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::reInitOperator ( const hiflow::la::lMatrix<ValueType> &op )
        {

            assert ( &op != NULL );
            this->Operator_ = &op;

            this->Operator_perm_ = NULL;

            if ( this->LU_ != NULL )
                delete this->LU_;

            if ( this->LU_cpu_ != NULL )
                delete this->LU_cpu_;

            // re allocate LU_ and LU_cpu_
            this->allocate_LU ( );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::extract_mat ( )
        {
            assert ( this->ncolors_ > 0 );

            int offset = 0;
            int offset_col = 0;
            int L_index = 0, R_index = 0;

            // allocate D,L,R array of matrices
            this->D_ = new lMatrix<ValueType>* [this->ncolors_];
            this->L_ = new lMatrix<ValueType>* [this->ncolors_ * ( this->ncolors_ - 1 ) / 2];
            this->R_ = new lMatrix<ValueType>* [this->ncolors_ * ( this->ncolors_ - 1 ) / 2];

            // Scaling Offdiaognal entries
            this->LU_cpu_->CloneFrom ( *this->LU_ );
            this->LU_cpu_->ScaleOffdiag ( this->omega_relax_ );
            this->LU_->CloneFrom ( *this->LU_cpu_ );

            for ( int i = 0; i<this->ncolors_; ++i )
            {

                if ( i >= 1 )
                    offset += this->color_sizes_[i - 1];

                //    D matrices
                this->D_[i] = this->LU_->extract_submatrix ( offset,
                                                             offset,
                                                             offset + this->color_sizes_[i],
                                                             offset + this->color_sizes_[i] );
                //        this->D_[i]->print();

                //        this->D_[i]->compress_me();

                offset_col = 0;
                for ( int j = 0; j<this->ncolors_; ++j )
                {

                    if ( j >= 1 )
                        offset_col += this->color_sizes_[j - 1];

                    // R matrices
                    if ( j > i )
                    {

                        LOG_DEBUG ( 7, "lPreconditioner_MultiColoring<ValueType>::extract_mat() Extracting R LU ("
                                    << this->ncolors_ * ( this->ncolors_ - 1 ) / 2 << "," << R_index << ") with sizes "
                                    << offset << "," << offset_col << " x "
                                    << offset + this->color_sizes_[i] << "," << offset_col + this->color_sizes_[j] );

                        this->R_[R_index] = this->LU_->extract_submatrix ( offset,
                                                                           offset_col,
                                                                           offset + this->color_sizes_[i],
                                                                           offset_col + this->color_sizes_[j] );

                        //        this->R_[R_index]->print();
                        //        this->R_[R_index]->compress_me();

                        ++R_index;
                    }

                    // L matrices
                    if ( j < i )
                    {

                        LOG_DEBUG ( 7, "lPreconditioner_MultiColoring<ValueType>::extract_mat() Extracting L LU ("
                                    << this->ncolors_ * ( this->ncolors_ - 1 ) / 2 << "," << L_index << ") with sizes "
                                    << offset << "," << offset_col << " x "
                                    << offset + this->color_sizes_[i] << "," << offset_col + this->color_sizes_[j] );

                        this->L_[L_index] = this->LU_->extract_submatrix ( offset,
                                                                           offset_col,
                                                                           offset + this->color_sizes_[i],
                                                                           offset_col + this->color_sizes_[j] );
                        //        this->L_[L_index]->print();
                        //        this->L_[L_index]->compress_me();

                        ++L_index;
                    }

                }
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::allocate_vec ( )
        {
            this->output_chunks_ = new lVector<ValueType>* [this->ncolors_];

            int offset = 0;
            for ( int i = 0; i<this->ncolors_; ++i )
            {

                if ( i >= 1 )
                    offset += this->color_sizes_[i - 1];

                // extract the sub vectors
                if ( this->x_const_ != NULL )
                {
                    this->output_chunks_[i] = this->x_const_->extract_subvector ( offset,
                                                                                  offset + this->color_sizes_[i] );
                }
                else
                {

                    assert ( this->x_ != NULL );
                    this->output_chunks_[i] = this->x_->extract_subvector ( offset,
                                                                            offset + this->color_sizes_[i] );
                }
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::convert_Dmat2vec ( )
        {
            this-> Dv_ = new lVector<ValueType>* [this->ncolors_];
            this->iDv_ = new lVector<ValueType>* [this->ncolors_];

            int offset = 0;
            for ( int i = 0; i<this->ncolors_; ++i )
            {

                if ( i >= 1 )
                    offset += this->color_sizes_[i - 1];

                // the sub diag vectors
                this->Dv_[i] = this->output_chunks_[i]->CloneWithoutContent ( );
                //    this->Dv_[i]->CloneFrom(*this->output_chunks_[i]);

                this->iDv_[i] = this->output_chunks_[i]->CloneWithoutContent ( );
                //    this->iDv_[i]->CloneFrom(*this->output_chunks_[i]);

                // extract them
                this->LU_cpu_->extract_diagelements ( offset,
                                                      offset + this->color_sizes_[i],
                                                      this->Dv_[i] );

                this->LU_cpu_->extract_invdiagelements ( offset,
                                                         offset + this->color_sizes_[i],
                                                         this->iDv_[i] );

            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::check_for_dropoffs ( )
        {

            int offset = 0;

            for ( int i = 0; i<this->ncolors_; ++i )
            {

                if ( i >= 1 )
                    offset += this->color_sizes_[i - 1];

                if ( ( this->D_[i]->get_nnz ( ) != this->D_[i]->get_num_row ( ) ) || ( this->D_[i]->get_nnz ( ) != this->D_[i]->get_num_col ( ) ) )
                {
                    LOG_INFO ( "lprecond mc", " D Matrix " << i << " has wrong nnz!!!!!!!!!!!!!" );
                    LOG_INFO ( "lprecond mc", "Extracting D LU (" << this->ncolors_ << "/" << i << ") with sizes "
                               << offset << "," << offset << " x "
                               << offset + this->color_sizes_[i] << "," << offset + this->color_sizes_[i] );

                    this->D_[i]->print ( );
                    this->D_[i]->WriteFile ( "/home/dlukarski/matrices/D_wrong.mtx" );
                }

            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::dropoffs ( )
        {

            CPU_CSR_lMatrix<ValueType> *mat_cpu = new CPUsimple_CSR_lMatrix<ValueType>;

            // Permute the Operator
            mat_cpu->CloneFrom ( *this->LU_cpu_ );

            // Simple MC dropoffs
            int ndel = 0;

            for ( int i = 0; i < mat_cpu->get_num_row ( ); ++i )
                for ( int j = mat_cpu->matrix.row[i]; j < mat_cpu->matrix.row[i + 1]; ++j )
                {

                    int offset = 0;
                    for ( int ic = 0; ic<this->ncolors_; ++ic )
                    {

                        if ( ic >= 1 )
                            offset += this->color_sizes_[ic - 1];

                        if ( ( i >= offset ) && ( i < offset + color_sizes_[ic] ) )
                            if ( ( mat_cpu->matrix.col[j] >= offset ) && ( mat_cpu->matrix.col[j] < offset + color_sizes_[ic] ) )
                                if ( i != mat_cpu->matrix.col[j] )
                                {
                                    mat_cpu->matrix.val[j] = 0.0;
                                    ++ndel;
                                }

                    }
                }

            LOG_DEBUG ( 2, "lprecond mc (drop off) Number of deleted entries = " << ndel );

            if ( ndel > 0 )
                mat_cpu->compress_me ( );

            this->LU_cpu_->Clear ( );
            this->LU_cpu_->CloneFrom ( *mat_cpu );

            this->LU_->Clear ( );
            this->LU_->CloneFrom ( *mat_cpu );

            delete mat_cpu;

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::cleanup ( )
        {
            // delete the diagonal matrices
            for ( int i = 0; i<this->ncolors_; ++i )
            {
                delete this->D_[i];
            }

            // delete the LU cpu matrix
            if ( this->LU_cpu_ != this->aux_mat_analysis_cpu_ )
                delete this->LU_cpu_;

            // delete LU
            delete this->LU_;

            // delete the aux cpu matrix
            delete aux_mat_analysis_cpu_;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                              hiflow::la::lVector<ValueType> *output )
        {

            this->splitvectors ( input );

            this->forwardstep ( );
            this->diagonalstep ( );
            this->backwardstep ( );

            this->concatenatevectors ( output );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::splitvectors ( const hiflow::la::lVector<ValueType> &input )
        {

            int offset = 0;

            for ( int i = 0; i<this->ncolors_; ++i )
            {

                if ( i >= 1 )
                    offset += this->color_sizes_[i - 1];

                // ouput = input
                // split ouput into chunks
                this->output_chunks_[i]->partial_replace_subvector ( 0,
                                                                     offset,
                                                                     this->color_sizes_[i],
                                                                     input );
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring<ValueType>::concatenatevectors ( hiflow::la::lVector<ValueType> *output )
        {

            int offset = 0;

            for ( int i = 0; i<this->ncolors_; ++i )
            {

                if ( i >= 1 )
                    offset += this->color_sizes_[i - 1];

                // output = ouput_chunks
                output->partial_replace_subvector ( offset,
                                                    0,
                                                    this->color_sizes_[i],
                                                    *this->output_chunks_[i] );
            }

        }

        template <typename ValueType>
        int* lPreconditioner_MultiColoring<ValueType>::get_permut ( )
        {
            return this->permut_index_;
        }

        template class lPreconditioner_MultiColoring<double>;
        template class lPreconditioner_MultiColoring<float>;

    } // namespace hiflow
} // namespace la
