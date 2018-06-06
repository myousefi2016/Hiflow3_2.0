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

#ifndef HIFLOW_LINEARALGEBRA_PETSC_MATRIX_H_
#    define HIFLOW_LINEARALGEBRA_PETSC_MATRIX_H_

#    include "common/log.h"
#    include "common/sorted_array.h"
#    include "config.h"
#    include "linear_algebra/matrix.h"
#    include "linear_algebra/petsc_vector.h"
#    include "mpi.h"

// TODO: Fix order dependency of next inclusion
#    include "common/smart_pointers.h"

namespace hiflow
{
    namespace la
    {

        /// Forwarding PETSc matrix wrapper class
        namespace petsc
        {
            class Mat_wrapper;
        }

        /// Forwarding PETSc linear solver
        template <class LAD>
        class PETScCG;

        /// @brief Wrapper class to PETSc matrix

        template <class DataType>
        class PETScMatrix : public Matrix<DataType>
        {
          public:
            /// Default constructor
            PETScMatrix ( );

            /// Destructor
            virtual ~PETScMatrix ( );

            virtual Matrix<DataType>* Clone ( ) const;

            /// Initialize matrix
            /// @param[in] comm MPI communicator to be used by matrix
            void Init ( const MPI_Comm& comm, const LaCouplings& cp );

            /// Initializes the structure of the matrix with pairs of indices.
            /// @param rows_diag Global row indices of diagonal block
            /// @param cols_diag Global column indices of diagonal block
            /// @param nnz_diag Size of @em rows_diag and @em cols_diag arrays
            /// @param rows_offdiag Global row indices of offdiagonal block
            /// @param cols_offdiag Global column indices of offdiagonal block
            /// @param nnz_offdiag Size of @em rows_offdiag and @em cols_offdiag arrays
            void InitStructure ( const int* rows_diag, const int* cols_diag, const int nnz_diag,
                                 const int* rows_offdiag, const int* cols_offdiag,
                                 const int nnz_offdiag );

            /// Clear all allocated data
            void Clear ( );

            /// Print statistical data
            void print_statistics ( ) const;

            /// out = this * in
            void VectorMult ( Vector<DataType>& in, Vector<DataType>* out ) const;

            void VectorMult ( PETScVector<DataType>& in, PETScVector<DataType>* out ) const;

            /// out = beta * out + alpha * this * in

            void VectorMultAdd ( DataType alpha, Vector<DataType>& in, DataType beta, Vector<DataType>* out ) const
            {
                LOG_ERROR ( "Function not available in PETSc interface" );
                exit ( -1 );
            }

            /// Get value at specified indices
            DataType GetValue ( int row, int col ) const;
            /// Get values at specified indices
            void GetValues ( const int* row_indices, int num_rows,
                             const int* col_indices, int num_cols, DataType* values ) const;

            /// Euclidean length of vector
            DataType NormFrobenius ( ) const;

            // Mutating functions: after calling any of these, a call to
            // begin_update()/end_update() or update() must be made before
            // any other function can be called. It is, however, possible
            // to call the same mutating function several times in a row,
            // without calling update() in between.

            /// Add value to given indices
            void Add ( int global_row_id, int global_col_id, DataType value );

            /// \brief Add submatrix of values at positions (rows x cols).
            /// The row and column numbers are assumed to correspond to global dof ids.
            /// Size of values is assumed to be |rows| x |cols|.
            void Add ( const int* rows, int num_rows, const int* cols, int num_cols, const DataType* values );
            /// Set value at given indices
            void SetValue ( int row, int col, DataType value );
            /// Set submatrix of values
            void SetValues ( const int* row_indices, int num_rows, const int* col_indices, int num_cols, const DataType* values );
            /// Set Matrix to zero
            void Zeros ( );

            /// Sets rows to zero except the diagonal element to alpha.
            /// @param row_indices Global row indices (must be owned by this process)
            /// @param num_rows Size of array @em row_indices
            /// @param diagonal_value Value to be set for diagonal element
            void diagonalize_rows ( const int* row_indices, int num_rows, DataType diagonal_value );

            /// Scale Matrix: this = alpha * this
            /// @param alpha Scaling factor
            void Scale ( const DataType alpha );

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
                int nrows_local = this->num_rows_local ( );
                int nrows_global = 0;
                MPI_Allreduce ( &nrows_local, &nrows_global, 1, MPI_INT, MPI_SUM, this->comm_ );
                return nrows_global;
            }
            /// Global number of columns

            int num_cols_global ( ) const
            {
                int ncols_local = this->num_cols_local ( );
                int ncols_global = 0;
                MPI_Allreduce ( &ncols_local, &ncols_global, 1, MPI_INT, MPI_SUM, this->comm_ );
                return ncols_global;
            }
            /// Local number of rows

            int num_rows_local ( ) const
            {
                return this->size_local_;
            }
            /// Local number of columns

            int num_cols_local ( ) const
            {
                return this->size_local_;
            }

          private:

            /// Friends
            template <class LAD> friend class PETScCG;
            template <class LAD> friend class PETScGeneralKSP;

            /// Final assemble of cached values. Must be called before using the matrix.
            void Assembly ( ) const;

            /// MPI communicator
            MPI_Comm comm_;
            /// Rank of current process
            int my_rank_;
            /// Global number of processes
            int nb_procs_;
            /// Linear algebra couplings describing global dof distribution
            const LaCouplings* cp_;

            /// Global number of first row owned by this process
            int ilower_;
            /// Global number of last row owned by this process
            int iupper_;
            /// Number of rows owned by this process
            int size_local_;

            /// Flag if matrix is initialized
            bool initialized_;

            /// Flag if matrix is assembled
            mutable bool assembled_;

            /// Structure of the matrix
            std::vector<SortedArray<int> > structure_;
            std::vector<int> row_length_diag_;
            std::vector<int> row_length_offdiag_;

            /// Pointer to PETSc matrix object
            hiflow::scoped_ptr<petsc::Mat_wrapper> ptr_mat_wrapper_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_PETSC_MATRIX_H_
