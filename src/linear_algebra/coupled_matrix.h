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

/// @author Chandramowli Subramanian, Nico Trost, Dimitar Lukarski, Martin Wlotzka

#ifndef HIFLOW_LINEARALGEBRA_COUPLED_MATRIX_H_
#    define HIFLOW_LINEARALGEBRA_COUPLED_MATRIX_H_

#    include <cstdio>
#    include <map>

#    include "linear_algebra/matrix.h"
#    include "linear_algebra/la_couplings.h"
#    include "linear_algebra/coupled_matrix_factory.h"
#    include "linear_algebra/lmp/la_global.h"
#    include "linear_algebra/lmp/lmatrix.h"
#    include "common/property_tree.h"
#    include "linear_algebra/lmp/platform_management.h"

#    include "mpi.h"

namespace hiflow
{
    namespace la
    {

        template<class DataType> class lMatrix;
        template<class DataType> class CoupledVector;

        /// @brief Distributed matrix.
        ///
        /// The entries are stored on local processes, the memory is allocated dynamically.
        /// Diagonal block and offdiagonal block are stored seperately.

        template<class DataType>
        class CoupledMatrix : public Matrix<DataType>
        {
          public:
            /// Standard constructor
            CoupledMatrix ( );
            /// Destructor
            virtual ~CoupledMatrix ( );

            /// Inits an empty (non-square) matrix.
            /// @param comm MPI communicator
            /// @param row_cp LaCouplings describing the distribution of the rows
            /// @param col_cp LaCouplings describing the distribution of the columns
            /// @param plat System platform (see lmp/la_global.h)
            /// @param impl Implementation (see lmp/la_global.h)
            /// @param format Matrix format (see lmp/la_global.h)
            void Init ( const MPI_Comm& comm, const LaCouplings& row_cp, const LaCouplings& col_cp,
                        PLATFORM plat, IMPLEMENTATION impl, MATRIX_FORMAT format );
            void Init ( const MPI_Comm& comm, const LaCouplings& row_cp, const LaCouplings& col_cp,
                        PLATFORM plat, IMPLEMENTATION impl, MATRIX_FORMAT format,
                        const SYSTEM& my_system );

            /// Inits an empty square matrix.
            /// @param comm MPI communicator
            /// @param cp LaCouplings describing the distribution of the rows and columns
            /// @param plat System platform (see lmp/la_global.h)
            /// @param impl Implementation (see lmp/la_global.h)
            /// @param format Matrix format (see lmp/la_global.h)
            void Init ( const MPI_Comm& comm, const LaCouplings& cp,
                        PLATFORM plat, IMPLEMENTATION impl, MATRIX_FORMAT format );
            void Init ( const MPI_Comm& comm, const LaCouplings& cp,
                        PLATFORM plat, IMPLEMENTATION impl, MATRIX_FORMAT format,
                        const SYSTEM& my_system );

            /// 3 new funcitons for Init
            void Init ( const MPI_Comm& comm_hf, const LaCouplings& cp );
            void Init ( const MPI_Comm& comm_hf, const LaCouplings& row_cp, const LaCouplings& col_cp );
            void Init_la_system ( PLATFORM plat, IMPLEMENTATION impl, MATRIX_FORMAT format );
            /// Initializes the structure of the matrix with pairs of indices.
            /// The LaCouplings describing the row and column distribution
            /// need to be initialized. Use rather @em CloneFromWithoutContent
            /// whenever possible.
            /// @param rows_diag Global row indices of diagonal block
            /// @param cols_diag Global column indices of diagonal block
            /// @param nnz_diag Size of @em rows_diag and @em cols_diag arrays
            /// @param rows_offdiag Global row indices of offdiagonal block
            /// @param cols_offdiag Global column indices of offdiagonal block
            /// @param nnz_offdiag Size of @em rows_offdiag and @em cols_offdiag arrays
            void InitStructure ( const int* rows_diag, const int* cols_diag, const int nnz_diag,
                                 const int* rows_offdiag, const int* cols_offdiag,
                                 const int nnz_offdiag );

            /// Compute Matrix-Vector product only with offdiagonal parts of local matrices
            /// @param[in] in input vector
            /// @param[out] out output vector
            void VectorMultOffdiag ( Vector<DataType>& in, Vector<DataType>* out ) const;

            // Member functions inherited from Matrix base class
            void VectorMult ( Vector<DataType>& in, Vector<DataType>* out ) const;
            void VectorMultAdd ( DataType alpha, Vector<DataType>& in, DataType beta, Vector<DataType>* out ) const;
            void diagonalize_rows ( const int* row_indices, int num_rows, DataType diagonal_value );
            void Scale ( const DataType alpha );
            void Add ( const int global_row_id, const int global_col_id, const DataType val );
            void Add ( const int* rows, int num_rows, const int* cols, int num_cols, const DataType* values );

            // Specializations of inherited member functions to kind of parallel implementation
            void VectorMult ( CoupledVector<DataType>& in, CoupledVector<DataType>* out ) const;
            void VectorMultAdd ( DataType alpha, CoupledVector<DataType>& in, DataType beta, CoupledVector<DataType>* out ) const;
            void VectorMultOffdiag ( CoupledVector<DataType>& in, CoupledVector<DataType>* out ) const;

            void GetValues ( const int* row_indices, int num_rows,
                             const int* col_indices, int num_cols, DataType* values ) const;

            virtual void SetValue ( int row, int col, DataType value )
            {
                LOG_INFO ( "CoupledMatrix::SetValue", "not implemented" );
            };

            virtual void SetValues ( const int* row_indices, int num_rows, const int* col_indices, int num_cols, const DataType* values )
            {
                LOG_INFO ( "CoupledMatrix::SetValues", "not implemented" );
            };

            /// Perform matrix-vector multiplication with the submatrix which is
            /// defined by the num_rows x num_cols matrix with the given row and
            /// column indices.
            /// @param row_indices - row indices of submatrix
            /// @param num_rows - number of rows in submatrix
            /// @param col_indices - column indices of submatrix
            /// @param num_cols - number of columns in submatrix
            /// @param in_values - input vector
            /// @param out_values - output_vector
            void VectorMult_submatrix ( const int* row_indices,
                                        int num_rows,
                                        const int* col_indices,
                                        int num_cols,
                                        const DataType* in_values,
                                        DataType* out_values ) const;

            /// Perform matrix-vector multiplication with the submatrix which is
            /// defined by the rows given in row_indices. The multiplication is
            /// performed EXCEPT the column indices specified in row_indices (yes, row_indices!!).
            /// This special multiplication is needed for the Vanka preconditioner.
            /// @param row_indices - row indices of submatrix
            /// @param num_rows - number of rows in submatrix
            /// @param col_indices - column indices of submatrix
            /// @param num_cols - number of columns in submatrix
            /// @param in_values - input vector
            /// @param out_values - output_vector
            void VectorMult_submatrix_vanka ( const int* row_indices,
                                              int num_rows,
                                              CoupledVector<DataType> &in,
                                              DataType* out_values ) const;

            /// Compress matrix by deleting zero elements
            void Compress ( void );

            /// Extracts CSR structure of matrix (only CPU matrices with CSR format)
            /// @param ia Indices of local row pointers (needs to be allocated)
            /// @param ja Global column indices (needs to be allocated)
            /// @param val Values (needs to be allocated)
            void ExtractCSR ( int* ia, int* ja, DataType* val ) const;

            /// Extracts CSR structure of diagonal block (only CPU matrices with CSR format)
            /// @param ia Indices of local row pointers (needs to be allocated)
            /// @param ja Local column indices (needs to be allocated)
            /// @param val Values (needs to be allocated)
            void ExtractDiagonalCSR ( int* ia, int* ja, DataType* val ) const;

            void ExtractDiagElements ( CoupledVector<DataType>& diag ) const;
            void ExtractInverseDiagElements ( CoupledVector<DataType>& inv_diag ) const;

            /// Sets every element to zero.
            void Zeros ( );

            virtual Matrix<DataType>* Clone ( ) const;

            /// Clones the whole matrix (everything).
            /// @param mat Matrix to copy
            void CloneFrom ( const CoupledMatrix<DataType>& mat );

            /// Copies only the values (no structure, no platform).
            /// Matrix already needs to be initialized.
            void CopyFrom ( const CoupledMatrix<DataType>& mat );

            void CastFrom ( const CoupledMatrix<double>& other );
            void CastFrom ( const CoupledMatrix<float>& other );

            /// Copies only the structure.
            void CopyStructureFrom ( const CoupledMatrix<DataType>& mat );

            /// Copies everything but the values.
            void CloneFromWithoutContent ( const CoupledMatrix<DataType>& mat );

            /// Convert CSR format to COO.
            /// @param mat Matrix to be converted from CSR to COO
            void ConvertFromCSR2COO ( const CoupledMatrix<DataType>& mat );

            /// Convert COO format to CSR.
            /// @param mat Matrix to be converted from COO to CSR
            void ConvertFromCOO2CSR ( const CoupledMatrix<DataType>& mat );

            /// Creates the transpose of a square matrix.
            /// @param other Matrix to create the transpose from.
            void CreateTransposedFrom ( const CoupledMatrix<DataType>& other );

            /// Clears allocated local matrices.
            void Clear ( );

            /// @return Global number of rows
            int nrows_global ( ) const;

            /// @return Global number of columns = global number of rows (square matrix!)
            int ncols_global ( ) const;

            /// Prints information to stream \em out .
            /// @param out Stream for output
            void Print ( std::ostream &out = std::cout ) const;

            // inline functions
            /// @return Local number of rows
            inline int nrows_local ( ) const;

            /// @return Local number of columns
            inline int ncols_local ( ) const;

            /// @return Local number of nonzero elements of diagonal block
            inline int nnz_local_diag ( ) const;

            /// @return Local number of nonzero elements of offdiagonal block
            inline int nnz_local_offdiag ( ) const;

            /// @return Local number of nonzero elements
            inline int nnz_local ( ) const;

            /// @return Global number of nonzero elements
            inline int nnz_global ( ) const;

            /// @return @c True if square matrix, otherwise @c false
            inline bool is_square_matrix ( ) const;

            /// Initiate update

            void begin_update ( )
            {
            }
            /// Finalize update

            void end_update ( )
            {
            }

            // Update matrix entries

            void Update ( )
            {
            }

            /// Global number of rows

            int num_rows_global ( ) const
            {
                return this->nrows_global ( );
            }
            /// Global number of columns

            int num_cols_global ( ) const
            {
                return this->ncols_global ( );
            }
            /// Local number of rows

            int num_rows_local ( ) const
            {
                return this->nrows_local ( );
            }
            /// Local number of columns

            int num_cols_local ( ) const
            {
                return this->ncols_local ( );
            }

            /// @return MPI communicator

            const MPI_Comm& comm ( ) const
            {
                return this->comm_;
            }

            /// @return LaCouplings for the distribution of the rows

            const LaCouplings& row_couplings ( ) const
            {
                return *( this->row_couplings_ );
            }

            /// @return LaCouplings for the distribution of the columns

            const LaCouplings& col_couplings ( ) const
            {
                return *( this->col_couplings_ );
            }

            /// @return Local diagonal block of the matrix

            const lMatrix<DataType>& diagonal ( ) const
            {
                return *( this->diagonal_ );
            }

            /// @return Local diagonal block of the matrix

            lMatrix<DataType>& diagonal ( )
            {
                return *( this->diagonal_ );
            }

            /// @return Local offdiagonal block of the matrix

            const lMatrix<DataType>& offdiagonal ( ) const
            {
                return *( this->offdiagonal_ );
            }

            /// @return Local offdiagonal block of the matrix

            lMatrix<DataType>& offdiagonal ( )
            {
                return *( this->offdiagonal_ );
            }

            /// @return Size of communicator

            int comm_size ( ) const
            {
                return this->comm_size_;
            }

            /// @return Rank of this process

            int my_rank ( ) const
            {
                return this->my_rank_;
            }

            /// @return Global index of the first local row

            int row_ownership_begin ( ) const
            {
                return this->row_ownership_begin_;
            }

            /// @return One more than the global index of the last local row

            int row_ownership_end ( ) const
            {
                return this->row_ownership_end_;
            }

            /// @return Global index of the first local column

            int col_ownership_begin ( ) const
            {
                return this->col_ownership_begin_;
            }

            /// @return One more than the global index of the last local column

            int col_ownership_end ( ) const
            {
                return this->col_ownership_end_;
            }

          private:
            // no implementation of copy constructor or assignement operator
            CoupledMatrix ( const CoupledMatrix<DataType>& );
            CoupledMatrix<DataType>& operator= ( const CoupledMatrix<DataType>& );

            /// Computes ownership range via LaCouplings.
            void ComputeOwnershipRange ( );

            MPI_Comm comm_;
            const LaCouplings* row_couplings_;
            const LaCouplings* col_couplings_;
            lMatrix<DataType>* diagonal_;
            lMatrix<DataType>* offdiagonal_;

            /// Size of MPI communicator.
            int comm_size_;
            /// Rank of this process.
            int my_rank_;
            /// Global index of the first local row.
            int row_ownership_begin_;
            /// One more than the global index of the last local row.
            int row_ownership_end_;
            /// Global index of the first local column.
            int col_ownership_begin_;
            /// One more than the global index of the last local column.
            int col_ownership_end_;

            int nnz_global_;
        };

        template<class DataType>
        inline int CoupledMatrix<DataType>::nrows_local ( ) const
        {
            assert ( this->diagonal_ != NULL );
            return this->diagonal_->get_num_row ( );
        }

        template<class DataType>
        inline int CoupledMatrix<DataType>::ncols_local ( ) const
        {
            assert ( this->diagonal_ != NULL );
            return this->diagonal_->get_num_col ( );
        }

        template<class DataType>
        inline int CoupledMatrix<DataType>::nnz_local_diag ( ) const
        {
            assert ( this->diagonal_ != NULL );
            return this->diagonal ( ).get_nnz ( );
        }

        template<class DataType>
        inline int CoupledMatrix<DataType>::nnz_local_offdiag ( ) const
        {
            assert ( this->offdiagonal_ != NULL );
            return this->offdiagonal ( ).get_nnz ( );
        }

        template<class DataType>
        inline int CoupledMatrix<DataType>::nnz_local ( ) const
        {
            return (this->nnz_local_diag ( ) + this->nnz_local_offdiag ( ) );
        }

        template<class DataType>
        inline int CoupledMatrix<DataType>::nnz_global ( ) const
        {
            return (this->nnz_global_ );
        }

        template<class DataType>
        inline bool CoupledMatrix<DataType>::is_square_matrix ( ) const
        {
            return (this->nrows_global ( ) == this->ncols_global ( ) );
        }

        template<class DataType>
        class CoupledMatrixCreator
        {
          public:

            CoupledMatrix<DataType>* params ( const PropertyTree& c )
            {
                CoupledMatrix<DataType>* newCoupledMatrix = new CoupledMatrix<DataType>( );
                if ( c.contains ( "Platform" ) && c.contains ( "Implementation" ) &&
                     c.contains ( "MatrixFormat" ) )
                {
                    const std::string platform_str = c["Platform"].template get<std::string>( );
                    if ( platform_str == "CPU" )
                    {
                        la_sys_.Platform = CPU;
                    }
                    else if ( platform_str == "GPU" )
                    {
                        la_sys_.Platform = GPU;
                    }
                    else
                    {
                        LOG_ERROR ( "CoupledMatrixCreator::params: No format of this name registered(platform)." );
                        return NULL;
                    }
                    init_platform ( la_sys_ );
                    const std::string impl_str = c["Implementation"].template get<std::string>( );
                    if ( impl_str == "Naive" )
                    {
                        la_impl_ = NAIVE;
                    }
                    else if ( impl_str == "BLAS" )
                    {
                        la_impl_ = BLAS;
                    }
                    else if ( impl_str == "MKL" )
                    {
                        la_impl_ = MKL;
                    }
                    else if ( impl_str == "OPENMP" )
                    {
                        la_impl_ = OPENMP;
                    }
                    else if ( impl_str == "SCALAR" )
                    {
                        la_impl_ = SCALAR;
                    }
                    else if ( impl_str == "SCALAR_TEX" )
                    {
                        la_impl_ = SCALAR_TEX;
                    }
                    else
                    {
                        LOG_ERROR ( "CoupledMatrixCreator::params: No format of this name registered(implementation)." );
                        return NULL;
                    }
                    const std::string matrix_str = c["MatrixFormat"].template get<std::string>( );
                    if ( matrix_str == "CSR" )
                    {
                        la_matrix_format_ = CSR;
                    }
                    else if ( matrix_str == "COO" )
                    {
                        la_matrix_format_ = COO;
                    }
                    else
                    {
                        LOG_ERROR ( "CoupledMatrixCreator::params: No format of this name registered(matrix format)." );
                        return NULL;
                    }
                    newCoupledMatrix->Init_la_system ( la_sys_.Platform, la_impl_, la_matrix_format_ );
                    return newCoupledMatrix;
                }
                return 0;
            }
          private:
            SYSTEM la_sys_;
            IMPLEMENTATION la_impl_;
            MATRIX_FORMAT la_matrix_format_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_COUPLED_MATRIX_H_
