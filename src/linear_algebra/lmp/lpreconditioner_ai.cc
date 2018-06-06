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

/// @author Dimitar Lukarski, Niels Wegh

#include "lpreconditioner_ai.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "lvector.h"
#include "lvector_cpu.h"

#include "lmatrix.h"
#include "lmatrix_csr_cpu.h"
#include "lmatrix_dense_cpu.h"
#include "lmp_log.h"

#include "solvers/cg.h"

#include <math.h>

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        lPreconditioner_ApproximateInverse<ValueType>::lPreconditioner_ApproximateInverse ( )
        {
            this->Operator_ = NULL;
            this->AI_ = NULL;
        }

        template <typename ValueType>
        lPreconditioner_ApproximateInverse<ValueType>::~lPreconditioner_ApproximateInverse ( )
        {
            this->Clear ( );
        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse<ValueType>::Init ( )
        {
            // do nothing
        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse<ValueType>::Clear ( )
        {
            if ( this->AI_ != NULL )
                delete this->AI_;

            this->AI_ = NULL;
            this->Operator_ = NULL;
        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                                   hiflow::la::lVector<ValueType> *output )
        {
            assert ( this->AI_ != NULL );

            this->AI_->VectorMult ( input, output );
        }

        template <typename ValueType>
        lPreconditioner_ApproximateInverse_FSAI<ValueType>::lPreconditioner_ApproximateInverse_FSAI ( )
        {
            this->precond_name_ = "Approximate Inverse - FSAI";

            // Parameters
            this->solver_max_iter_ = 1000;
            this->solver_rel_eps_ = 1e-10;
            this->solver_abs_eps_ = 1e-10;
            this->matrix_power_ = 1;

            this->mult_tmp_ = NULL;
            this->ext_matrix_pat_ = false;

            this->AI_L_ = NULL;
            this->AI_Lt_ = NULL;
            this->mult_tmp_ = NULL;

        }

        template <typename ValueType>
        lPreconditioner_ApproximateInverse_FSAI<ValueType>::~lPreconditioner_ApproximateInverse_FSAI ( )
        {
            this->mult_tmp_->Clear ( );
        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse_FSAI<ValueType>::Init ( const int solver_max_iter, const ValueType solver_rel_eps, const ValueType solver_abs_eps, const int matrix_power )
        {
            assert ( solver_max_iter > 0 );
            this->solver_max_iter_ = solver_max_iter;

            assert ( solver_rel_eps > 0 );
            this->solver_rel_eps_ = solver_rel_eps;

            assert ( solver_abs_eps > 0 );
            this->solver_abs_eps_ = solver_abs_eps;

            assert ( matrix_power >= 1 );
            this->matrix_power_ = matrix_power;

            char data_info[255];
            this->precond_name_ = "Approximate Inverse - FSAI";
            sprintf ( data_info, "(%d)", matrix_power );
            this->precond_name_.append ( data_info );
        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse_FSAI<ValueType>::set_matrix_power ( const int matrix_power )
        {
            assert ( matrix_power >= 1 );
            this->matrix_power_ = matrix_power;

            char data_info[255];
            this->precond_name_ = "Approximate Inverse - FSAI";
            sprintf ( data_info, "(%d)", matrix_power );
            this->precond_name_.append ( data_info );
        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse_FSAI<ValueType>::set_ext_matrix_pattern ( const hiflow::la::lMatrix<ValueType> &mat )
        {
            assert ( &mat != NULL );

            this->matrix_pat_ = new CPUsimple_CSR_lMatrix<ValueType>;
            this->matrix_pat_->CloneFrom ( mat );

            for ( int i = 0; i<this->matrix_pat_->get_nnz ( ); ++i )
                this->matrix_pat_->matrix.val[i] = 1.0;

            this->ext_matrix_pat_ = true;
        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse_FSAI<ValueType>::Init ( )
        {
            this->solver_max_iter_ = 50000;
            this->solver_rel_eps_ = 1e-10;
            this->solver_abs_eps_ = 1e-10;
        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse_FSAI<ValueType>::SetupVector ( const hiflow::la::lVector<ValueType> *x )
        {
            assert ( x != NULL );

            this->mult_tmp_ = x->CloneWithoutContent ( );
            mult_tmp_->Zeros ( );
        }

        // Source code by Niels Wegh

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse_FSAI<ValueType>::Build ( )
        {
            // For using the FSAI method, we need to solve the system: A G_j = e_j (for (i,j) in P)
            // P is a prescribed triangular pattern of the preconditioner.
            // A is a dense matrix. The values of A correspond to the values of the coefficient matrix and the structure arise from the pattern P.
            // G_j is the j-th column of this pattern (only the nonzeros)
            // e_j is a vector with only entry corresponding to the diagonal is nonzero.

            int dense_size; // the size of the small dense matrix
            int dense_col; // column index of the dense matrix
            int start_col; // starting column index
            int pat_col; // the column index of the pattern
            int insert; // continuous index for inserting the sol_ vector into the preconditioner
            int max_iter; // maximum number of iterations

            bool Stop = false; // stop criterion for loops
            ValueType buffer; // temporary buffer for values

            CPU_CSR_lMatrix<ValueType> *mat = new CPUsimple_CSR_lMatrix<ValueType>; // the orginal coefficient matrix
            CPU_CSR_lMatrix<ValueType> *pat = new CPUsimple_CSR_lMatrix<ValueType>; // the prescribed pattern of the preconditioner

            mat->Init ( this->Operator_->get_nnz ( ),
                        this->Operator_->get_num_row ( ),
                        this->Operator_->get_num_col ( ),
                        "Matrix" );

            mat->CloneFrom ( *this->Operator_ );

            // Build or load the pattern of the preconditioner
            if ( this->ext_matrix_pat_ )
            {

                pat->CloneFrom ( *this->matrix_pat_ );

                delete this->matrix_pat_;
                this->matrix_pat_ = NULL;

            }
            else
            {

                // Ensure symmetric pattern if Dirichlet BC are prescribed

                for ( int i = 0; i < mat->get_num_row ( ); ++i )
                    for ( int row = mat->matrix.row[i]; row < mat->matrix.row[i + 1]; ++row )
                        if ( mat->matrix.val[row] == 0.0 )
                        {

                            for ( int col = mat->matrix.row[ mat->matrix.col[row] ]; col < mat->matrix.row[mat->matrix.col[row] + 1]; ++col )
                                if ( i == mat->matrix.col[col] )
                                {
                                    mat->matrix.val[col] = 0.0;

                                }
                        }

                mat->compress_me ( );

                if ( mat->issymmetric ( ) == false )
                {
                    LOG_ERROR ( "lprecond FSAI - input matrix is not symmetric!!!" );
                    exit ( -1 );
                }

                lMatrix<ValueType> *matrix_cpu_power = new CPUsimple_CSR_lMatrix<ValueType>;

                matrix_cpu_power = mat->MatrixSupSPower ( this->matrix_power_ );
                pat->CloneFrom ( *matrix_cpu_power );

                delete matrix_cpu_power;
            }

            //  pat->CloneFrom(*this->Operator_);

            pat->delete_lower_triangular ( );

            for ( int i = 0; i < mat->get_num_row ( ); ++i )
            {

                //The size of the dense matrix is equal to the nonzero elements in G_j
                dense_size = 0;
                dense_size = pat->matrix.row[i + 1] - pat->matrix.row[i];

                //If only the diagonal entry of G_j is nonzero
                if ( dense_size == 1 )
                {
                    for ( int row = mat->matrix.row[i]; row < mat->matrix.row[i + 1]; ++row )
                    {
                        if ( mat->matrix.col[row] == i )
                        {
                            pat->matrix.val[ pat->matrix.row[i] ] = 1.0 / mat->matrix.val[row];
                            break;
                        }
                    }
                }

                if ( dense_size > 1 )
                {
                    lVector<ValueType> *sol_, *rhs_;
                    lMatrix<ValueType> *matrix_;
                    // matrix_ is the dense matrix (size: dense_size x dense_size )
                    // rhs_ is the right hand side and only the entry corresponding to the diagonal is nonzero (e_j)
                    // sol_ is the solution vector and the entries of sol_ correspond to the nonzero elements of the precondtioner (G_j)
                    matrix_ = new CPUsimple_DENSE_lMatrix<ValueType>( dense_size*dense_size, dense_size, dense_size, "dense Matrix" );
                    sol_ = new CPUsimple_lVector<ValueType>( dense_size, "sol" );
                    rhs_ = new CPUsimple_lVector<ValueType>( dense_size, "rhs" );

                    // Set the initial guess to zero
                    matrix_->Zeros ( );
                    sol_->Zeros ( );
                    rhs_->Zeros ( );

                    // Set the diagonal to 1
                    rhs_->add_value ( dense_size - 1, 1.0 );

                    // Build the symmetric dense matrix by going through all indices (dense_row, dense_col) and filling them with the corresponding values.
                    for ( int dense_row = 0; dense_row < dense_size; ++dense_row )
                    {

                        start_col = pat->matrix.col[ pat->matrix.row[i] + dense_row];
                        dense_col = 0;

                        for ( int mat_col = mat->matrix.row[ start_col ]; mat_col < mat->matrix.row[ start_col + 1 ]; ++mat_col )
                        {

                            pat_col = pat->matrix.col[ pat->matrix.row[i] + dense_col ];

                            if ( mat->matrix.col[mat_col] < pat_col )
                                continue;

                            while ( mat->matrix.col[mat_col] > pat_col )
                            {
                                dense_col += 1;
                                pat_col = pat->matrix.col[ pat->matrix.row[i] + dense_col ];
                                if ( dense_col >= dense_size )
                                {
                                    Stop = true;
                                    break;
                                }
                            }

                            if ( Stop == true )
                            {
                                Stop = false;
                                break;
                            }

                            if ( mat->matrix.col[mat_col] == pat_col )
                            {

                                matrix_->add_value ( dense_row, dense_col, mat->matrix.val[mat_col] );
                                dense_col += 1;

                                if ( dense_col >= dense_size )
                                {
                                    break;
                                }
                            }

                        }
                    }
                    max_iter = this->solver_max_iter_;

                    // Solve the linear equation system and insert the entries of sol_ into the precondtioner
                    cg ( sol_, rhs_, matrix_, this->solver_rel_eps_, this->solver_abs_eps_, max_iter, -1 );

                    insert = 0;

                    for ( int j = pat->matrix.row[i]; j < pat->matrix.row[i + 1]; ++j )
                    {

                        if ( pat->matrix.col[j] > i )
                        {
                            LOG_ERROR ( "lpreconditioner_ai FSAI - insert below the diagonal" );
                            exit ( -1 );
                            break;
                        }
                        sol_->GetValues ( &insert, 1, &buffer );
                        insert += 1;
                        if ( insert > dense_size )
                        {
                            LOG_ERROR ( "lpreconditioner_ai FSAI - too many elements to be inserted" );
                            exit ( -1 );
                        }
                        pat->matrix.val[j] = buffer;

                    }

                    delete sol_;
                    delete rhs_;
                    delete matrix_;

                }
                if ( dense_size < 1 )
                {
                    LOG_ERROR ( "lpreconditioner_ai FSAI - zero column" );
                    exit ( -1 );
                }
            }

            // Diagonal scaling of the preconditioner, in order to satisfy G A G^T = I ( for (i,j) in the prescribed pattern)
            for ( int i = 0; i < pat->get_num_row ( ); ++i )
            {
                buffer = sqrt ( 1.0 / pat->matrix.val[pat->matrix.row[i + 1] - 1] );
                for ( int j = pat->matrix.row[i]; j < pat->matrix.row[i + 1]; ++j )
                {
                    pat->matrix.val[j] = pat->matrix.val[j] * buffer;
                }
            }

            this->AI_L_ = this->Operator_->CloneWithoutContent ( );
            this->AI_L_->CloneFrom ( *pat );
            //  pat->WriteFile("FSAI_L.mtx");

            pat->transpose_me ( );
            this->AI_Lt_ = this->Operator_->CloneWithoutContent ( );
            this->AI_Lt_->CloneFrom ( *pat );
            //  pat->WriteFile("FSAI_Lt.mtx");

            delete mat;
            delete pat;

            //  this->AI_L_->print();
            //  this->AI_Lt_->print();

            char data_info[255];
            sprintf ( data_info, " / NNZ(L+L^T)=%d", this->AI_L_->get_nnz ( ) + this->AI_L_->get_nnz ( ) );
            this->precond_name_.append ( data_info );

        }

        template <typename ValueType>
        void lPreconditioner_ApproximateInverse_FSAI<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                                        hiflow::la::lVector<ValueType> *output )
        {
            assert ( this->AI_L_ != NULL );
            assert ( this->AI_Lt_ != NULL );
            assert ( this->mult_tmp_ != NULL );

            this->AI_L_->VectorMult ( input, this->mult_tmp_ );
            this->AI_Lt_->VectorMult ( *this->mult_tmp_, output );

        }

        template class lPreconditioner_ApproximateInverse<double>;
        template class lPreconditioner_ApproximateInverse<float>;

        template class lPreconditioner_ApproximateInverse_FSAI<double>;
        template class lPreconditioner_ApproximateInverse_FSAI<float>;

    } // namespace hiflow
} // namespace la
