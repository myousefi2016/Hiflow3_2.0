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

#include "lpreconditioner.h"
#include "lpreconditioner_mc.h"
#include "lpreconditioner_mc_ilup.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "lvector.h"
#include "lvector_cpu.h"

#include "lmatrix.h"
#include "lmatrix_csr_cpu.h"

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        lPreconditioner_MultiColoring_ILUp<ValueType>::lPreconditioner_MultiColoring_ILUp ( )
        {
            this->ilu_p_ = 0;
            this->mat_power_ = 1;
            this->precond_name_ = "Multi-Coloring / ILU ";

            this->flag_mc_ = true;
            this->flag_dropoff_mc_ = true;
            this->flag_ls_ = false;

            this->ext_power_pat_ = false;
            this->ext_ilup_pat_ = false;
        }

        template <typename ValueType>
        lPreconditioner_MultiColoring_ILUp<ValueType>::~lPreconditioner_MultiColoring_ILUp ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::Init ( )
        {
            // use default data
            // ILU(0,1) MC, DROPOFF, NO LS
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::set_ext_power_pattern ( const lMatrix<ValueType> &mat )
        {
            assert ( &mat != NULL );

            this->aux_mat_analysis_cpu_ = new CPUsimple_CSR_lMatrix<ValueType>;
            this->aux_mat_analysis_cpu_->CloneFrom ( mat );
            this->ext_power_pat_ = true;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::set_ext_ilup_pattern ( const lMatrix<ValueType> &mat )
        {
            assert ( &mat != NULL );

            this->ext_ilup_pat_mat_ = new CPUsimple_CSR_lMatrix<ValueType>;
            this->ext_ilup_pat_mat_->CloneFrom ( mat );
            this->ext_ilup_pat_ = true;
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::Init ( const int ilu_p,
                                                                   const int mat_power,
                                                                   const bool mc,
                                                                   const bool dropoff_mc,
                                                                   const bool ls )
        {
            assert ( this->ilu_p_ >= 0 );
            this->ilu_p_ = ilu_p;

            assert ( this->mat_power_ >= 1 );
            this->mat_power_ = mat_power;

            this->flag_mc_ = mc;
            this->flag_dropoff_mc_ = dropoff_mc;
            this->flag_ls_ = ls;

            char data_info[255];
            sprintf ( data_info, "(%d,%d)", ilu_p, mat_power );
            this->precond_name_.append ( data_info );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::build_aux_matrix ( )
        {
            assert ( this->Operator_ != NULL );

            if ( this->ext_power_pat_ )
            {

                // do nothing
                // this->aux_mat_analysis_cpu_ is set

                assert ( this->aux_mat_analysis_cpu_ != NULL );

            }
            else
            {

                if ( this->mat_power_ == 1 )
                {

                    this->aux_mat_analysis_cpu_ = new CPUsimple_CSR_lMatrix<ValueType>;
                    this->aux_mat_analysis_cpu_->CloneFrom ( *this->Operator_ );

                }
                else
                {

                    CPUsimple_CSR_lMatrix<ValueType> *matrix_cpu = new CPUsimple_CSR_lMatrix<ValueType>;
                    matrix_cpu->CloneFrom ( *this->Operator_ );

                    this->aux_mat_analysis_cpu_ = matrix_cpu->MatrixSupSPower ( this->mat_power_ );

                    delete matrix_cpu;
                }
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::allocate_LU ( )
        {
            this->LU_cpu_ = new CPUsimple_CSR_lMatrix<ValueType>;
            this->LU_cpu_->CloneFrom ( *this->Operator_ );

            // make LU the same as Operator_
            this->LU_ = this->Operator_->CloneWithoutContent ( );

            // copy the ILU(p) matrix to LU_
            this->LU_->CloneFrom ( *this->LU_cpu_ );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::factorize ( )
        {

            // do ILU(p) factorization on the cpu
            if ( this->ilu_p_ == 0 )
            {

                this->LU_cpu_->ilu0 ( );
                //    this->LU_cpu_->ilup(0);

            }
            else
            {

                if ( this->ext_ilup_pat_ )
                {
                    // use external pattern
                    this->LU_cpu_->ilup ( *this->ext_ilup_pat_mat_, this->ilu_p_ );
                }
                else
                {
                    // buld inside the pattern
                    this->LU_cpu_->ilup ( this->ilu_p_ );
                }

                //    this->LU_cpu_->ilusp(this->ilu_p_, this->ncolors_, this->color_sizes_, this->permut_index_);
                //    this->LU_cpu_->WriteFile("LU.mtx");
            }

            // make LU the same as Operator_perm_
            this->LU_ = this->Operator_->CloneWithoutContent ( );

            // copy the ILU(p) matrix to LU_
            this->LU_->CloneFrom ( *this->LU_cpu_ );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::scale_D ( )
        {

            for ( int i = 0; i<this->ncolors_; ++i )
            {
                this->iDv_[i]->Scale ( -1.0 );
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::forwardstep ( )
        {

            int l_index = 0; // forward accessing

            for ( int i = 0; i<this->ncolors_; ++i )
            {

                // this->output_chunks_[i]:=-this->output_chunks_[i];
                this->output_chunks_[i]->Scale ( -1.0 );

                // this->output_chunks_[i] += this->L_*this->output_chunks_[j]
                for ( int j = 0; j < i; ++j )
                {
                    this->L_[l_index]->VectorMultAdd ( *this->output_chunks_[j], this->output_chunks_[i] );
                    ++l_index;
                }

                this->output_chunks_[i]->Scale ( -1.0 );
                //    this->output_chunks_[i]->ElementWiseMult(*this->iDv_[i]);
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::diagonalstep ( )
        {
            // Diag step
            // z:=-z
            for ( int i = 0; i<this->ncolors_; ++i )
                this->output_chunks_[i]->Scale ( -1.0 );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_ILUp<ValueType>::backwardstep ( )
        {

            int r_index = this->ncolors_ * ( this->ncolors_ - 1 ) / 2 - 1; // backward accessing

            for ( int i = this->ncolors_ - 1; i >= 0; --i )
            {

                // this->output_chunks_[i] += this->R_*this->output_chunks_[j]
                for ( int j = this->ncolors_ - 1; j > i; --j )
                {
                    this->R_[r_index]->VectorMultAdd ( *this->output_chunks_[j], this->output_chunks_[i] );
                    --r_index;
                }

                // this->output_chunks_[i] := -(D^-1)*this->output_chunks_[i];
                //    this->output_chunks_[i]->Scale(-1.0); // in side the iDv_
                this->output_chunks_[i]->ElementWiseMult ( *this->iDv_[i] );
            }

        }

        template class lPreconditioner_MultiColoring_ILUp<double>;
        template class lPreconditioner_MultiColoring_ILUp<float>;

    } // namespace hiflow
} // namespace la
