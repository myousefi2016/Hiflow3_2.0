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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-11-26

#include "common/pointers.h"
#include "common/log.h"
#include "cxx98_prettyprint.h"
#include "petsc.h"
#include "petsc_environment.h"
#include "petsc_matrix.h"
#include "petsc_matrix_interface.h"
#include "petsc_vector_interface.h"
#include <vector>
#include <cstdlib>

namespace hiflow
{
    namespace la
    {

        template <class DataType>
        PETScMatrix<DataType>::PETScMatrix ( )
        : comm_ ( MPI_COMM_NULL ),
        my_rank_ ( ),
        nb_procs_ ( ),
        cp_ ( NULL ),
        ilower_ ( ),
        iupper_ ( ),
        size_local_ ( ),
        initialized_ ( false ),
        assembled_ ( false ),
        structure_ ( ),
        row_length_diag_ ( ),
        row_length_offdiag_ ( ),
        ptr_mat_wrapper_ ( new petsc::Mat_wrapper )
        {
        }

        template <class DataType>
        PETScMatrix<DataType>::~PETScMatrix ( )
        {
            this->Clear ( );
        }

        template <class DataType>
        Matrix<DataType>* PETScMatrix<DataType>::Clone ( ) const
        {
            LOG_ERROR ( "PETScMatrix::Clone not yet implemented!!!" );
            exit ( -1 );
            return NULL;
        }

        template <class DataType>
        void PETScMatrix<DataType>::Init ( const MPI_Comm& comm, const LaCouplings& cp )
        {
            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
            }

            assert ( comm != MPI_COMM_NULL );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( comm, &nb_procs_ );
            assert ( info == MPI_SUCCESS );
            assert ( nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( comm, &my_rank_ );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank_ >= 0 );
            assert ( my_rank_ < nb_procs_ );

            info = MPI_Comm_split ( comm, 0, my_rank_, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            this->cp_ = &cp;
            assert ( this->cp_ != NULL );
        }

        template <class DataType>
        void PETScMatrix<DataType>::InitStructure ( const int* rows_diag, const int* cols_diag, const int nnz_diag,
                                                    const int* rows_offdiag, const int* cols_offdiag,
                                                    const int nnz_offdiag )
        {
            // clear possibly existing data
            this->Clear ( );

            // Get information about size of local matrix
            ilower_ = this->cp_->dof_offset ( my_rank_ );
            iupper_ = ilower_ + this->cp_->nb_dofs ( my_rank_ ) - 1;
            size_local_ = this->cp_->nb_dofs ( my_rank_ );

            structure_.resize ( size_local_ );
            row_length_diag_.resize ( size_local_, 0 );
            row_length_offdiag_.resize ( size_local_, 0 );

            // 1. step: fill structure_ and row_length_diag_ with diagonal part
            for ( int i = 0; i < nnz_diag; ++i )
            {
                structure_[rows_diag[i] - ilower_].find_insert ( cols_diag[i] );
                row_length_diag_[rows_diag[i] - ilower_] += 1;
            }

            // 2. step: fill structure_ and row_length_offdiag_ with offdiagonal part
            for ( int i = 0; i < nnz_offdiag; ++i )
            {
                structure_[rows_offdiag[i] - ilower_].find_insert ( cols_offdiag[i] );
                row_length_offdiag_[rows_offdiag[i] - ilower_] += 1;
            }

            // Prepare PETSc MPI interface
            PETScEnvironment::initialize ( );

            // Create PETSC Matrix
            MatCreate ( comm_, &ptr_mat_wrapper_->mat_ );
            MatCreateAIJ ( comm_, cp_->nb_dofs ( my_rank_ ), cp_->nb_dofs ( my_rank_ ), cp_->nb_total_dofs ( ),
                           cp_->nb_total_dofs ( ), 0, vec2ptr ( row_length_diag_ ), 0, vec2ptr ( row_length_offdiag_ ), &ptr_mat_wrapper_->mat_ );

            initialized_ = true;
            assembled_ = true;
        }

        template <class DataType>
        void PETScMatrix<DataType>::Clear ( )
        {
            if ( initialized_ )
            {
                structure_.clear ( );
                row_length_diag_.clear ( );
                row_length_offdiag_.clear ( );
                MatDestroy ( &ptr_mat_wrapper_->mat_ );
                initialized_ = false;
            }
        }

        template <class DataType>
        void PETScMatrix<DataType>::VectorMult ( Vector<DataType>& in, Vector<DataType>* out ) const
        {
            PETScVector<DataType> *hv_in, *hv_out;

            hv_in = dynamic_cast < PETScVector<DataType>* > ( &in );
            hv_out = dynamic_cast < PETScVector<DataType>* > ( out );

            if ( ( hv_in != 0 ) && ( hv_out != 0 ) )
            {
                this->VectorMult ( *hv_in, hv_out );
            }
            else
            {
                if ( hv_in == 0 )
                {
                    LOG_ERROR ( "Called PETScMatrix::VectorMult with incompatible input vector type." );
                }
                if ( hv_out == 0 )
                {
                    LOG_ERROR ( "Called PETScMatrix::VectorMult with incompatible output vector type." );
                }
                exit ( -1 );
            }
        }

        template <class DataType>
        void PETScMatrix<DataType>::VectorMult ( PETScVector<DataType>& in, PETScVector<DataType>* out ) const
        {
            Assembly ( );
            MatMult ( ptr_mat_wrapper_->mat_, in.ptr_vec_wrapper_->vec_, out->ptr_vec_wrapper_->vec_ );
        }

        template <class DataType>
        DataType PETScMatrix<DataType>::GetValue ( int row, int col ) const
        {
            assert ( initialized_ );
            Assembly ( );
            DataType value;
            MatGetValues ( ptr_mat_wrapper_->mat_, 1, &row, 1, &col, &value );
            return value;
        }

        template <class DataType>
        void PETScMatrix<DataType>::GetValues ( const int* row_indices, int num_rows,
                                                const int* col_indices, int num_cols, DataType* values ) const
        {
            assert ( initialized_ );
            Assembly ( );
            MatGetValues ( ptr_mat_wrapper_->mat_, num_rows, row_indices, num_cols, col_indices, values );
        }

        template <class DataType>
        DataType PETScMatrix<DataType>::NormFrobenius ( ) const
        {
            Assembly ( );
            DataType value;
            MatNorm ( ptr_mat_wrapper_->mat_, NORM_FROBENIUS, &value );
            return value;
        }

        template <class DataType>
        void PETScMatrix<DataType>::Add ( int global_row_id, int global_col_id, DataType value )
        {
            this->Add ( &global_row_id, 1, &global_col_id, 1, &value );
        }

        template <class DataType>
        void PETScMatrix<DataType>::Add ( const int* row_indices, int num_rows, const int* col_indices, int num_cols, const DataType* values )
        {
            for ( size_t i = 0; i != num_rows; ++i )
            {
                std::vector<int> col_ind_row;
                std::vector<DataType> val_row;
                col_ind_row.reserve ( num_cols );
                val_row.reserve ( num_cols );

                const int current_row_global = row_indices[i];
                const int current_row_local = current_row_global - ilower_;
                int k = 0;
                const int row_offset = i * num_cols;
                for ( size_t j = 0, j_e = structure_[current_row_local].size ( ); k < num_cols && j != j_e; ++j )
                {
                    if ( col_indices[k] == structure_[current_row_local][j] )
                    {
                        col_ind_row.push_back ( col_indices[k] );
                        val_row.push_back ( values[row_offset + k] );
                        ++k;
                    }
                    else if ( col_indices[k] < structure_[current_row_local][j] )
                    {
                        while ( col_indices[k] < structure_[current_row_local][j] && k < num_cols )
                        {
                            ++k;
                        }
                        if ( k >= num_cols )
                        {
                            break;
                        }
                        if ( col_indices[k] == structure_[current_row_local][j] )
                        {
                            col_ind_row.push_back ( col_indices[k] );
                            val_row.push_back ( values[row_offset + k] );
                            ++k;
                        }
                    }
                }
                int ncols = col_ind_row.size ( );
                if ( ncols > 0 )
                {
                    MatSetValues ( ptr_mat_wrapper_->mat_, num_rows, row_indices, num_cols, col_indices, values, ADD_VALUES );
                }
            }
        }

        template <class DataType>
        void PETScMatrix<DataType>::SetValue ( int row, int col, DataType value )
        {
            this->SetValues ( &row, 1, &col, 1, &value );
        }

        template <class DataType>
        void PETScMatrix<DataType>::SetValues ( const int* row_indices, int num_rows,
                                                const int* col_indices, int num_cols, const DataType* values )
        {
            assert ( initialized_ );
            for ( size_t i = 0; i != num_rows; ++i )
            {
                std::vector<int> col_ind_row;
                std::vector<DataType> val_row;
                col_ind_row.reserve ( num_cols );
                val_row.reserve ( num_cols );

                const int current_row_global = row_indices[i];
                const int current_row_local = current_row_global - ilower_;
                int k = 0;
                const int row_offset = i * num_cols;
                for ( size_t j = 0, j_e = structure_[current_row_local].size ( ); k < num_cols && j != j_e; ++j )
                {
                    if ( col_indices[k] == structure_[current_row_local][j] )
                    {
                        col_ind_row.push_back ( col_indices[k] );
                        val_row.push_back ( values[row_offset + k] );
                        ++k;
                    }
                    else if ( col_indices[k] < structure_[current_row_local][j] )
                    {
                        while ( col_indices[k] < structure_[current_row_local][j] && k < num_cols )
                        {
                            ++k;
                        }
                        if ( k >= num_cols )
                        {
                            break;
                        }
                        if ( col_indices[k] == structure_[current_row_local][j] )
                        {
                            col_ind_row.push_back ( col_indices[k] );
                            val_row.push_back ( values[row_offset + k] );
                            ++k;
                        }
                    }
                }
                int ncols = col_ind_row.size ( );
                if ( ncols > 0 )
                {
                    MatSetValues ( ptr_mat_wrapper_->mat_, 1, &current_row_global, ncols, vec2ptr ( col_ind_row ), vec2ptr ( val_row ), INSERT_VALUES );
                }
            }
            assembled_ = false;
        }

        template <class DataType>
        void PETScMatrix<DataType>::Zeros ( )
        {
            Assembly ( );
            MatZeroEntries ( ptr_mat_wrapper_->mat_ );
            assembled_ = false;
        }

        template <class DataType>
        void PETScMatrix<DataType>::diagonalize_rows ( const int* row_indices, int num_rows, DataType diagonal_value )
        {
            for ( size_t i = 0; i != num_rows; ++i )
            {
                assert ( row_indices[i] >= ilower_ && row_indices[i] <= iupper_ );
                const int row_index_global = row_indices[i];
                const int row_index_local = row_index_global - ilower_;
                int ncols = structure_[row_index_local].size ( );
                std::vector<DataType> val ( ncols, 0. );

                for ( size_t j = 0, e_j = ncols; j != e_j; ++j )
                {
                    if ( row_index_global == structure_[row_index_local][j] )
                    {
                        val[j] = diagonal_value;
                    }
                }
                //HYPRE_IJMatrixSetValues(A_, 1, &ncols, &row_index_global, &(structure_[row_index_local].front()), vec2ptr(val));
                MatSetValues ( ptr_mat_wrapper_->mat_, 1, &row_index_global, ncols, &( structure_[row_index_local].front ( ) ), vec2ptr ( val ), INSERT_VALUES );
            }
        }

        template <class DataType>
        void PETScMatrix<DataType>::Scale ( DataType alpha )
        {
            assert ( initialized_ );
            Assembly ( );
            MatScale ( ptr_mat_wrapper_->mat_, alpha );
        }

        template <class DataType>
        void PETScMatrix<DataType>::Assembly ( ) const
        {
            if ( assembled_ ) return;
            MatAssemblyBegin ( ptr_mat_wrapper_->mat_, MAT_FINAL_ASSEMBLY );
            MatAssemblyEnd ( ptr_mat_wrapper_->mat_, MAT_FINAL_ASSEMBLY );
            assembled_ = true;
        }

        // template instantiation
        template class PETScMatrix<double>;

    } // namespace la
} // namespace hiflow
