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

/// @author Simon Gawlok

#include "linear_algebra/hypre_matrix.h"
#include "common/pointers.h"
#include "common/log.h"
#include <vector>
#include <cstdlib>

namespace hiflow
{

    namespace la
    {

        template<class DataType>
        HypreMatrix<DataType>::HypreMatrix ( )
        {
            this->initialized_ = false;
            this->comm_ = MPI_COMM_NULL;
            this->cp_row_ = NULL;
            this->cp_col_ = NULL;
        }

        template<class DataType>
        HypreMatrix<DataType>::~HypreMatrix ( )
        {
#ifdef WITH_HYPRE
            this->Clear ( );

            int is_finalized;
            MPI_Finalized ( &is_finalized );
            if ( !is_finalized )
            {
                if ( this->comm_ != MPI_COMM_NULL )
                    MPI_Comm_free ( &this->comm_ );
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        Matrix<DataType>* HypreMatrix<DataType>::Clone ( ) const
        {
            LOG_ERROR ( "HypreMatrix::Clone not yet implemented!!!" );
            exit ( -1 );
            return NULL;
        }

        template<class DataType>
        void HypreMatrix<DataType>::Init ( const MPI_Comm &comm, const LaCouplings &cp )
        {
            // clear possibly existing DataType
            if ( initialized_ )
            {
                this->Clear ( );
            }

            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
            }

            assert ( comm != MPI_COMM_NULL );

            MPI_Comm_dup ( comm, &this->comm_ );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( this->comm_, &nb_procs_ );
            assert ( info == MPI_SUCCESS );
            assert ( nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( this->comm_, &my_rank_ );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank_ >= 0 );
            assert ( my_rank_ < nb_procs_ );

            this->cp_row_ = &cp;
            this->cp_col_ = &cp;
            assert ( this->cp_row_ != NULL );
            assert ( this->cp_col_ != NULL );
            assert ( this->cp_row_->initialized ( ) );
            assert ( this->cp_col_->initialized ( ) );
        }

        template<class DataType>
        void HypreMatrix<DataType>::Init ( const MPI_Comm &comm, const LaCouplings &cp_row, const LaCouplings &cp_col )
        {
            // clear possibly existing DataType
            if ( initialized_ )
            {
                this->Clear ( );
            }

            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
            }

            assert ( comm != MPI_COMM_NULL );

            MPI_Comm_dup ( comm, &this->comm_ );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( this->comm_, &nb_procs_ );
            assert ( info == MPI_SUCCESS );
            assert ( nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( this->comm_, &my_rank_ );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank_ >= 0 );
            assert ( my_rank_ < nb_procs_ );

            this->cp_row_ = &cp_row;
            this->cp_col_ = &cp_col;
            assert ( this->cp_row_ != NULL );
            assert ( this->cp_col_ != NULL );
            assert ( this->cp_row_->initialized ( ) );
            assert ( this->cp_col_->initialized ( ) );
        }

        template<class DataType>
        void HypreMatrix<DataType>::InitStructure ( const int* rows_diag, const int* cols_diag, const int nnz_diag,
                                                    const int* rows_offdiag, const int* cols_offdiag,
                                                    const int nnz_offdiag )
        {
#ifdef WITH_HYPRE

            // Get information about size of local matrix
            ilower_ = this->cp_row_->dof_offset ( my_rank_ );
            iupper_ = ilower_ + this->cp_row_->nb_dofs ( my_rank_ ) - 1;
            jlower_ = this->cp_col_->dof_offset ( my_rank_ );
            jupper_ = jlower_ + this->cp_col_->nb_dofs ( my_rank_ ) - 1;
            size_local_ = this->cp_row_->nb_dofs ( my_rank_ );
            nnz_local_ = nnz_diag + nnz_offdiag;

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

            // Create the HYPRE matrix
            HYPRE_IJMatrixCreate ( comm_, ilower_, iupper_, jlower_, jupper_, &A_ );

            // Use parallel csr format
            HYPRE_IJMatrixSetObjectType ( A_, HYPRE_PARCSR );

            HYPRE_IJMatrixSetPrintLevel ( A_, 100 );

            HYPRE_IJMatrixSetDiagOffdSizes ( A_, vec2ptr ( row_length_diag_ ), vec2ptr ( row_length_offdiag_ ) );

            // Tell HYPRE that no matrix entries need to be communicated to other processors
            HYPRE_IJMatrixSetMaxOffProcElmts ( A_, 0 );

            // Initialize
            HYPRE_IJMatrixInitialize ( A_ );

            // Now initialize exact structure of matrix. To achieve this we set every element to zero
            this->Zeros ( );

            HYPRE_IJMatrixAssemble ( A_ );
            HYPRE_IJMatrixGetObject ( A_, ( void** ) &parcsr_A_ );
            initialized_ = true;

            int m, n;
            HYPRE_ParCSRMatrixGetDims ( parcsr_A_, &m, &n );
            LOG_INFO ( "Global number of rows", m );
            LOG_INFO ( "Global number of columns", n );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::Clear ( )
        {
#ifdef WITH_HYPRE
            structure_.clear ( );
            row_length_diag_.clear ( );
            row_length_offdiag_.clear ( );
            if ( initialized_ )
            {
                HYPRE_IJMatrixDestroy ( A_ );
            }
            initialized_ = false;

#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::print_statistics ( ) const
        {
#ifdef WITH_HYPRE

            int ilower, iupper, jlower, jupper;
            HYPRE_IJMatrixGetLocalRange ( this->A_, &ilower, &iupper, &jlower, &jupper );

            std::vector<int> rows ( iupper - ilower + 1, 0 );
#    pragma clang loop vectorize(enable)
            for ( int i = 0; i < iupper - ilower + 1; ++i )
            {
                rows[i] = ilower + i;
            }

            std::vector<int> nnz_count ( iupper - ilower + 1, 0 );
            HYPRE_IJMatrixGetRowCounts ( A_, iupper - ilower + 1, vec2ptr ( rows ), vec2ptr ( nnz_count ) );

            // print statistics
            for ( int i = 0; i < nb_procs_; ++i )
            {
                MPI_Barrier ( comm_ );
                if ( i == my_rank_ )
                {
                    std::cout << "HypreMatrix on process " << my_rank_ << ":" << std::endl;
                    // print size information
                    std::cout << "\t ilower: " << ilower << std::endl;
                    std::cout << "\t iupper: " << iupper << std::endl;
                    std::cout << "\t jlower: " << jlower << std::endl;
                    std::cout << "\t jupper: " << jupper << std::endl;
                    std::cout << "\t Nonzero elements (row: nnz)" << std::endl;
                    for ( int j = 0; j < iupper - ilower + 1; ++j )
                    {
                        std::cout << "\t\t " << ilower + j << ": " << nnz_count[j] << std::endl;
                    }
                }
                MPI_Barrier ( comm_ );
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::VectorMult ( Vector<DataType>& in, Vector<DataType>* out ) const
        {
            HypreVector<DataType> *hv_in, *hv_out;

            hv_in = dynamic_cast < HypreVector<DataType>* > ( &in );
            hv_out = dynamic_cast < HypreVector<DataType>* > ( out );

            if ( ( hv_in != 0 ) && ( hv_out != 0 ) )
            {
                this->VectorMult ( *hv_in, hv_out );
            }
            else
            {
                if ( hv_in == 0 )
                {
                    LOG_ERROR ( "Called HypreMatrix::VectorMult with incompatible input vector type." );
                }
                if ( hv_out == 0 )
                {
                    LOG_ERROR ( "Called HypreMatrix::VectorMult with incompatible output vector type." );
                }
                exit ( -1 );
            }
        }

        template<class DataType>
        void HypreMatrix<DataType>::VectorMult ( HypreVector<DataType>& in, HypreVector<DataType>* out ) const
        {
#ifdef WITH_HYPRE
            HYPRE_ParCSRMatrixMatvec ( 1., parcsr_A_, *( in.GetParVector ( ) ), 0., *( out->GetParVector ( ) ) );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::VectorMultAdd ( DataType alpha, Vector<DataType>& in, DataType beta, Vector<DataType>* out ) const
        {
            HypreVector<DataType> *hv_in, *hv_out;

            hv_in = dynamic_cast < HypreVector<DataType>* > ( &in );
            hv_out = dynamic_cast < HypreVector<DataType>* > ( out );

            if ( ( hv_in != 0 ) && ( hv_out != 0 ) )
            {
                this->VectorMultAdd ( alpha, *hv_in, beta, hv_out );
            }
            else
            {
                if ( hv_in == 0 )
                {
                    LOG_ERROR ( "Called HypreMatrix::VectorMult with incompatible input vector type." );
                }
                if ( hv_out == 0 )
                {
                    LOG_ERROR ( "Called HypreMatrix::VectorMult with incompatible output vector type." );
                }
                exit ( -1 );
            }
        }

        template<class DataType>
        void HypreMatrix<DataType>::VectorMultAdd ( DataType alpha, HypreVector<DataType>& in, DataType beta, HypreVector<DataType>* out ) const
        {
#ifdef WITH_HYPRE
            HYPRE_ParCSRMatrixMatvec ( alpha, parcsr_A_, *( in.GetParVector ( ) ), beta, *( out->GetParVector ( ) ) );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::GetValues ( const int* row_indices, int num_rows,
                                                const int* col_indices, int num_cols, DataType* values ) const
        {
#ifdef WITH_HYPRE    
            for ( size_t i = 0; i < num_rows; ++i )
            {
                std::vector<int> col_ind_row, col_map_row;
                std::vector<DataType> val_row;
                col_ind_row.reserve ( num_cols );
                col_map_row.reserve ( num_cols );

                int current_row_global = row_indices[i];
                const int current_row_local = current_row_global - ilower_;
                int k = 0;
                const int* const struct_curr_row = &( structure_[current_row_local][0] );
                for ( size_t j = 0, j_e = structure_[current_row_local].size ( ); k < num_cols && j != j_e; ++j )
                {
                    if ( col_indices[k] == struct_curr_row[j] )
                    {
                        col_ind_row.push_back ( col_indices[k] );
                        col_map_row.push_back ( k );
                        ++k;
                    }
                    else if ( col_indices[k] < struct_curr_row[j] )
                    {
                        while ( col_indices[k] < struct_curr_row[j] && k < num_cols )
                        {
                            ++k;
                        }
                        if ( k >= num_cols )
                        {
                            break;
                        }
                        if ( col_indices[k] == struct_curr_row[j] )
                        {
                            col_ind_row.push_back ( col_indices[k] );
                            col_map_row.push_back ( k );
                            ++k;
                        }
                    }
                }
                int ncols = col_ind_row.size ( );
                val_row.resize ( col_ind_row.size ( ), 0. );
                if ( ncols > 0 )
                {
                    HYPRE_IJMatrixGetValues ( A_, 1, &ncols, &current_row_global, vec2ptr ( col_ind_row ), vec2ptr ( val_row ) );

                    const size_t offset = i*num_cols;
                    for ( size_t j = 0; j < ncols; ++j )
                    {
                        values[offset + col_map_row[j]] = val_row[j];
                    }
                }
                col_ind_row.clear ( );
                col_map_row.clear ( );
                val_row.clear ( );
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::Add ( int global_row_id, int global_col_id, DataType value )
        {
            this->Add ( &global_row_id, 1, &global_col_id, 1, &value );
        }

        template<class DataType>
        void HypreMatrix<DataType>::Add ( const int* rows, int num_rows, const int* cols, int num_cols, const DataType* values )
        {
#ifdef WITH_HYPRE
            int nrows_add = 0;
            std::vector<int> ncols_add;
            ncols_add.reserve ( num_rows );
            std::vector<int> rows_add;
            rows_add.reserve ( num_rows );
            std::vector<int> cols_add;
            cols_add.reserve ( num_rows * num_cols );
            std::vector<DataType> vals_add;
            vals_add.reserve ( num_rows * num_cols );

            for ( size_t i = 0; i != num_rows; ++i )
            {
                int cols_in_row = 0;

                const int current_row_global = rows[i];
                const int current_row_local = current_row_global - ilower_;
                int k = 0;
                const int row_offset = i*num_cols;
                const int* const struct_curr_row = &( structure_[current_row_local][0] );
                for ( size_t j = 0, j_e = structure_[current_row_local].size ( ); k < num_cols && j != j_e; ++j )
                {
                    if ( cols[k] == struct_curr_row[j] )
                    {
                        cols_add.push_back ( cols[k] );
                        vals_add.push_back ( values[row_offset + k] );
                        ++cols_in_row;
                        ++k;
                    }
                    else if ( cols[k] < struct_curr_row[j] )
                    {
                        while ( cols[k] < struct_curr_row[j] && k < num_cols )
                        {
                            ++k;
                        }
                        if ( k >= num_cols )
                        {
                            break;
                        }
                        if ( cols[k] == struct_curr_row[j] )
                        {
                            cols_add.push_back ( cols[k] );
                            vals_add.push_back ( values[row_offset + k] );
                            ++cols_in_row;
                            ++k;
                        }
                    }
                }
                ncols_add.push_back ( cols_in_row );
                rows_add.push_back ( current_row_global );
                ++nrows_add;
            }
            HYPRE_IJMatrixAddToValues ( A_, nrows_add, vec2ptr ( ncols_add ), vec2ptr ( rows_add ), vec2ptr ( cols_add ), vec2ptr ( vals_add ) );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::SetValue ( int row, int col, DataType value )
        {
            this->SetValues ( &row, 1, &col, 1, &value );
        }

        template<class DataType>
        void HypreMatrix<DataType>::SetValues ( const int* row_indices, int num_rows, const int* col_indices, int num_cols, const DataType* values )
        {
#ifdef WITH_HYPRE
            int nrows_add = 0;
            std::vector<int> ncols_add;
            ncols_add.reserve ( num_rows );
            std::vector<int> rows_add;
            rows_add.reserve ( num_rows );
            std::vector<int> cols_add;
            cols_add.reserve ( num_rows * num_cols );
            std::vector<DataType> vals_add;
            vals_add.reserve ( num_rows * num_cols );

            for ( size_t i = 0; i != num_rows; ++i )
            {
                int cols_in_row = 0;

                const int current_row_global = row_indices[i];
                const int current_row_local = current_row_global - ilower_;
                int k = 0;
                const int row_offset = i*num_cols;
                const int* const struct_curr_row = &( structure_[current_row_local][0] );
                for ( size_t j = 0, j_e = structure_[current_row_local].size ( ); k < num_cols && j != j_e; ++j )
                {
                    if ( col_indices[k] == struct_curr_row[j] )
                    {
                        cols_add.push_back ( col_indices[k] );
                        vals_add.push_back ( values[row_offset + k] );
                        ++cols_in_row;
                        ++k;
                    }
                    else if ( col_indices[k] < struct_curr_row[j] )
                    {
                        while ( col_indices[k] < struct_curr_row[j] && k < num_cols )
                        {
                            ++k;
                        }
                        if ( k >= num_cols )
                        {
                            break;
                        }
                        if ( col_indices[k] == struct_curr_row[j] )
                        {
                            cols_add.push_back ( col_indices[k] );
                            vals_add.push_back ( values[row_offset + k] );
                            ++cols_in_row;
                            ++k;
                        }
                    }
                }
                ncols_add.push_back ( cols_in_row );
                rows_add.push_back ( current_row_global );
                ++nrows_add;
            }
            HYPRE_IJMatrixSetValues ( A_, nrows_add, vec2ptr ( ncols_add ), vec2ptr ( rows_add ), vec2ptr ( cols_add ), vec2ptr ( vals_add ) );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::Zeros ( )
        {
#ifdef WITH_HYPRE

            int nrows_add = structure_.size ( );
            std::vector<int> ncols_add;
            ncols_add.reserve ( structure_.size ( ) );
            std::vector<int> rows_add;
            rows_add.reserve ( structure_.size ( ) );
            std::vector<int> cols_add;
            cols_add.reserve ( nnz_local_ );
            std::vector<DataType> vals_add ( nnz_local_, 0. );

            for ( size_t i = 0, e_i = structure_.size ( ); i != e_i; ++i )
            {
                cols_add.insert ( cols_add.end ( ), structure_[i].begin ( ), structure_[i].end ( ) );

                ncols_add.push_back ( structure_[i].size ( ) );
                rows_add.push_back ( ilower_ + i );
            }
            HYPRE_IJMatrixSetValues ( A_, nrows_add, vec2ptr ( ncols_add ), vec2ptr ( rows_add ), vec2ptr ( cols_add ), vec2ptr ( vals_add ) );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::diagonalize_rows ( const int* row_indices, int num_rows, DataType diagonal_value )
        {
#ifdef WITH_HYPRE
            for ( size_t i = 0; i != num_rows; ++i )
            {
                assert ( row_indices[i] >= ilower_ && row_indices[i] <= iupper_ );
                const int row_index_global = row_indices[i];
                const int row_index_local = row_index_global - ilower_;
                int ncols = structure_[row_index_local].size ( );
                std::vector<DataType> val ( ncols, 0. );

                int pos = -1;
                structure_[row_index_local].find ( row_index_global, &pos );
                assert ( pos >= 0 );
                val[pos] = diagonal_value;

                HYPRE_IJMatrixSetValues ( A_, 1, &ncols, &row_index_global, &( structure_[row_index_local].front ( ) ), vec2ptr ( val ) );

                val.clear ( );
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreMatrix<DataType>::Scale ( const DataType alpha )
        {
#ifdef WITH_HYPRE
            for ( size_t i = 0, e_i = structure_.size ( ); i != e_i; ++i )
            {
                int ncols = structure_[i].size ( );
                if ( ncols > 0 )
                {
                    int global_row_index = ilower_ + i;
                    std::vector<DataType> val ( ncols, 0. );
                    HYPRE_IJMatrixGetValues ( A_, 1, &ncols, &global_row_index, &( structure_[i].front ( ) ), vec2ptr ( val ) );

                    for ( size_t j = 0; j != ncols; ++j )
                    {
                        val[j] *= alpha;
                    }
                    HYPRE_IJMatrixSetValues ( A_, 1, &ncols, &global_row_index, &( structure_[i].front ( ) ), vec2ptr ( val ) );

                    val.clear ( );
                }
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

#ifdef WITH_HYPRE

        template<class DataType>
        HYPRE_ParCSRMatrix* HypreMatrix<DataType>::GetParCSRMatrix ( )
        {
            HYPRE_IJMatrixGetObject ( A_, ( void** ) &parcsr_A_ );
            return &parcsr_A_;
        }

        template<class DataType>
        const HYPRE_ParCSRMatrix* HypreMatrix<DataType>::GetParCSRMatrix ( ) const
        {
            return &parcsr_A_;
        }
#endif

        template class HypreMatrix<double>;
    } // namespace la
} // namespace hiflow
