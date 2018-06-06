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

#ifndef HIFLOW_LINEARALGEBRA_HYPRE_MATRIX_H_
#    define HIFLOW_LINEARALGEBRA_HYPRE_MATRIX_H_

#    include <cstdlib>
#    include <iostream>
#    include <mpi.h>
#    include <map>
#    include <set>
#    include "config.h"
#    include "linear_algebra/matrix.h"
#    include "linear_algebra/hypre_vector.h"
#    include "common/log.h"
#    include "common/sorted_array.h"
#    ifdef WITH_HYPRE
#        include "_hypre_utilities.h"
#        include "HYPRE.h"
#        include "HYPRE_parcsr_ls.h"
#        include "HYPRE_parcsr_mv.h"
#    endif

namespace hiflow
{
    namespace la
    {

        /// @author Simon Gawlok

        /// @brief Wrapper to HYPRE matrix

        template<class DataType>
        class HypreMatrix : public Matrix<DataType>
        {
          public:
            /// Standard constructor
            HypreMatrix ( );
            /// Destructor
            ~HypreMatrix ( );

            virtual Matrix<DataType>* Clone ( ) const;

            /// Initialize matrix
            /// @param[in] comm MPI communicator to be used by matrix
            void Init ( const MPI_Comm &comm, const LaCouplings &cp );

            /// Initialize matrix
            void Init ( const MPI_Comm &comm, const LaCouplings &cp_row, const LaCouplings &cp_col );

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

            void VectorMult ( HypreVector<DataType>& in, HypreVector<DataType>* out ) const;

            /// out = beta * out + alpha * this * in
            void VectorMultAdd ( DataType alpha, Vector<DataType>& in, DataType beta, Vector<DataType>* out ) const;

            /// out = beta * out + alpha * this * in
            void VectorMultAdd ( DataType alpha, HypreVector<DataType>& in, DataType beta, HypreVector<DataType>* out ) const;

            /// Get values at specified indices
            void GetValues ( const int* row_indices, int num_rows,
                             const int* col_indices, int num_cols, DataType* values ) const;

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

#    ifdef WITH_HYPRE
            /// Get pointer to HXPRE_ParCSRMatrix objects
            HYPRE_ParCSRMatrix* GetParCSRMatrix ( );

            /// Get pointer to HXPRE_ParCSRMatrix objects
            const HYPRE_ParCSRMatrix* GetParCSRMatrix ( ) const;
#    endif

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
                return this->iupper_ - this->ilower_ + 1;
            }
            /// Local number of columns

            int num_cols_local ( ) const
            {
                return this->jupper_ - this->jlower_ + 1;
            }

            /// Get MPI communicator

            const MPI_Comm& comm ( ) const
            {
                return comm_;
            }

            /// Get MPI communicator

            MPI_Comm& comm ( )
            {
                return comm_;
            }

            /// Get LaCouplings for rows

            const LaCouplings& row_coulings ( ) const
            {
                return *cp_row_;
            }

            /// Get LaCouplings for columns

            const LaCouplings& col_couplings ( ) const
            {
                return *cp_col_;
            }

          private:
            /// MPI communicator
            MPI_Comm comm_;
            /// Rank of current process
            int my_rank_;
            /// Global number of processes
            int nb_procs_;
            /// Linear algebra couplings describing global row dof distribution
            const LaCouplings* cp_row_;
            /// Linear algebra couplings describing global column dof distribution
            const LaCouplings* cp_col_;

            /// Global number of first row owned by this process
            int ilower_;
            /// Global number of last row owned by this process
            int iupper_;
            /// Global number of first column owned by this process
            int jlower_;
            /// Global number of last column owned by this process
            int jupper_;
            /// Number of rows owned by this process
            int size_local_;
            /// Number of non-zero entries on this process
            int nnz_local_;

            /// Flag if vector is already initialized
            bool initialized_;

            /// Structure of the matrix
            std::vector< SortedArray<int> > structure_;
            std::vector<int> row_length_diag_;
            std::vector<int> row_length_offdiag_;
#    ifdef WITH_HYPRE
            HYPRE_IJMatrix A_;
            HYPRE_ParCSRMatrix parcsr_A_;
#    endif

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_HYPRE_MATRIX_H_
