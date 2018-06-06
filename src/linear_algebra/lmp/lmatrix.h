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

#ifndef __LMATRIX_H
#    define __LMATRIX_H

#    include <iostream>
#    include <stdlib.h>

#    include "la_global.h"
#    include "lmatrix_formats.h"
#    include "lvector.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Provides the base class to all local multi-platform matrices
        /// @author Dimitar Lukarski
        ///
        /// lMatrix (local Matrix) - lmpLAtoolbox (Local Multi-Platform Linear Algebra Tool Box)
        /// Provides the base class to all local matrices for different platforms
        /// (like CPU,GPU,etc) with different implementations (sequential,parallel,mkl,cblas,...)

        template <typename ValueType>
        class lMatrix
        {
          public:

            lMatrix ( );
            /// The virtual destructor free the buffers of the matrix
            virtual ~lMatrix ( );

            /// Compress matrix by deleting zero entries
            virtual void Compress ( void );

            /// Return the number of columns in the matrix
            /// @return The number of columns in the matrix
            inline int get_num_col ( void ) const;

            /// Return the number of rows in the matrix
            /// @return The number of rows in the matrix
            inline int get_num_row ( void ) const;

            /// Needed for transpose_me().
            /// Internal COO matrix has no access to protected variables.
            inline void swap_dimensions ( void );

            /// Return the number of non-zero elements in the matrix
            /// @return The number of non-zeros elements in the matrix
            inline int get_nnz ( void ) const;

            /// Return the matrix format
            /// @return The matrix format
            enum MATRIX_FORMAT get_matrix_format ( void ) const;

            /// Return the name of the matrix format
            /// @return The name of the matrix format
            std::string get_matrix_format_name ( void ) const;

            virtual void get_as_coo ( std::vector<ValueType>& vals, std::vector<int>& rows, std::vector<int>& cols ) const;

            /// Return the name of the matrix
            /// @return The name of the matrix
            std::string get_name ( void ) const;

            /// Provides human-readable information for the matrix - matrix name,
            /// nnz, row, col, precision, platform, implementation,
            /// @param out - output stream
            virtual void print ( std::ostream &out = std::cout ) const;

            /// Return the implementation method of the Blas Level 2 routines
            /// @return The implementation method of the Blas Level 2 routines
            enum IMPLEMENTATION get_implementation ( void ) const;

            /// Return the implementation name
            /// @return The implementation name
            std::string get_implementation_name ( void ) const;

            /// Return the platform of the matrix
            /// @return The platform of the matrix
            enum PLATFORM get_platform ( void ) const;

            /// Return the platform name of the matrix
            /// @return The platform name of the matrix
            std::string get_platform_name ( void ) const;

            /// Initialize (i.e. allocate) the matrix
            /// @param init_nnz - nnz
            /// @param init_row - number of row
            /// @param init_col - number of col
            /// @param init_name - the name of the matrix
            virtual void Init ( const int init_nnz,
                                const int init_num_row,
                                const int init_num_col,
                                const std::string init_name ) = 0;

            /// Clear (i.e. free) the matrix
            virtual void Clear ( void ) = 0;

            /// Set a specified row set to zero (the diagonal elements = alpha)
            /// @param index_set - the set of rows
            /// @param size - the size of the set
            /// @param alpha - the scalar value for the diagonal elements
            virtual void ZeroRows ( const int *index_set,
                                    const int size,
                                    const ValueType alpha );

            virtual void ZeroCols ( const int *index_set,
                                    const int size,
                                    const ValueType alpha );

            /// Set the matrix to zero
            virtual void Zeros ( void ) = 0;

            /// Reorder the matrix
            /// @param index - the dest index permutation vector
            virtual void Reorder ( const int *index );

            /// Multi-coloring algorithm
            /// @return ncolors - the number of colors
            /// @return color_sizes - the number of elements with respect to the colors
            /// @return - the permuration vector (dest index)
            virtual void Multicoloring ( int &ncolors, int **color_sizes, int **permut_index ) const;

            // LS
            virtual void Levelscheduling ( int &nlevels, int **level_sizes, int **permut_index ) const;

            /// ILU(0) Builder
            virtual void ilu0 ( void );

            /// ILU(p) Builder
            virtual void ilup ( const int p );

            /// ILU(p) Builder
            virtual void ilup ( const lMatrix<ValueType> &mat, const int p );

            /// Solver for the (I)LU preconditioners
            virtual void ilu_solve ( const lVector<ValueType> &invec,
                                     lVector<ValueType> *outvec ) const;

            /// Scale the matrix with a scalar
            /// @param alpha - the scalar
            virtual void Scale ( const ValueType alpha );

            /// Scale the off diagonal elements with a scalar
            /// @param alpha - the scalar
            virtual void ScaleOffdiag ( const ValueType alpha );

            /// Initialize the structure of the matrix with pairs of indexes
            /// @param rows - row array contains the row index
            /// @param rows - col array contains the col index
            virtual void init_structure ( const int *rows, const int *cols );

            /// Add a value to a specific index pair; check for the
            /// structure is not performed
            /// @param row - the row index
            /// @param col - the col index
            /// @param val - the value
            virtual void add_value ( const int row, const int col, const ValueType val );

            /// Add values to a index pairs rows x cols; check for the
            /// structure is not performed
            /// @param rows - the row indices
            /// @param num_rows - number of row indices
            /// @param cols - the column indices
            /// @param num_cols - number of column indices
            /// @param values - the value s
            virtual void add_values ( const int* rows,
                                      int num_rows,
                                      const int* cols,
                                      int num_cols,
                                      const ValueType* values );

            /// Get a value to a specific index pair; check for the
            /// structure is not performed
            /// @param row - the row index
            /// @param col - the col index
            /// @param val - the value
            virtual void get_value ( const int row, const int col, ValueType *val ) const;

            /// Get values at specified indices and add them to values; check for the
            /// structure is not performed
            /// @param rows - the row indices
            /// @param num_rows - number of row indices
            /// @param cols - the column indices
            /// @param num_cols - number of column indices
            /// @param cols_target - the column indices in the target values
            /// @param num_cols_target - number of columns in target values
            /// @param values - the value s
            virtual void get_add_values ( const int* rows,
                                          int num_rows,
                                          const int* cols,
                                          int num_cols,
                                          const int* cols_target,
                                          int num_cols_target,
                                          ValueType* values ) const;

            /// Perform matrix-vector multiplication with the submatrix which is
            /// defined by the num_rows x num_cols matrix with the given row and
            /// column indices.
            /// @param rows - the row indices
            /// @param num_rows - number of row indices
            /// @param cols - the column indices
            /// @param num_cols - number of column indices
            /// @param cols_input - the column indices in the input vector
            /// @param in_values - input vector
            /// @param out_values - output vector
            virtual void VectorMultAdd_submatrix ( const int* rows,
                                                   int num_rows,
                                                   const int* cols,
                                                   int num_cols,
                                                   const int* cols_input,
                                                   const ValueType* in_values,
                                                   ValueType* out_values ) const;

            /// Perform matrix-vector multiplication with the submatrix which is
            /// defined by the rows given in rows.
            /// This special multiplication is needed for the Vanka preconditioner.
            /// @param rows - the row indices
            /// @param num_rows - number of row indices
            /// @param invec - input vector
            /// @param out_values - output vector
            virtual void VectorMultAdd_submatrix_vanka ( const int* rows,
                                                         int num_rows,
                                                         const hiflow::la::lVector< ValueType > &invec,
                                                         ValueType* out_values ) const;

            /// Transpose the matrix
            virtual void transpose_me ( void );

            /// Check is the matrix symmetric
            virtual bool issymmetric ( void );

            /// Compress the matrix (delete the entries for val=0.0)
            virtual void compress_me ( void );

            /// Delete the diagonal entries of the matrix
            virtual void delete_diagonal ( void );

            /// Delete the off diagonal entries of the matrix
            virtual void delete_offdiagonal ( void );

            /// Delete the lower triangular entries of the matrix
            virtual void delete_lower_triangular ( void );

            /// Delete the stricktly lower triangular entries of the matrix
            virtual void delete_strictly_lower_triangular ( void );

            /// Delete the upper triangular entries of the matrix
            virtual void delete_upper_triangular ( void );

            /// Delete the stricktly upper triangular entries of the matrix
            virtual void delete_strictly_upper_triangular ( void );

            /// Extract a sub matrix. The returned matrix
            /// has the same platform/implementation as the
            /// original one.
            /// @param start_row - beginning row of the block
            /// @param start_col - beginning col of the block
            /// @param end_row - ending row of the block
            /// @param end_col - ending col of the block
            /// @return lVector<ValueType>
            virtual lMatrix<ValueType> *extract_submatrix ( const int start_row, const int start_col,
                                                            const int end_row, const int end_col ) const;

            /// Extract the diagonal elements of the matrix
            /// @parm start_i - begin index
            /// @parm end_i - end index
            /// @return vec - insert the values in this vector (the vector has the allocated)
            virtual void extract_diagelements ( const int start_i, const int end_i, lVector<ValueType> *vec ) const;

            /// Extract the inverse diagonal elements of the matrix
            virtual void extract_invdiagelements ( const int start_i, const int end_i, lVector<ValueType> *vec ) const;

            /// See CopyTo() and CopyFrom()
            virtual lMatrix<ValueType> &operator= ( const lMatrix<ValueType> &mat2 ) = 0;

            /// Copy from and Copy to operators - can be used for different platforms
            virtual void CopyFrom ( const lMatrix<ValueType> &mat2 ) = 0;
            virtual void CopyTo ( lMatrix<ValueType> &mat2 ) const = 0;

            virtual void CastFrom ( const lMatrix<double>& other ) = 0;
            virtual void CastFrom ( const lMatrix<float>& other ) = 0;

            virtual void CastTo ( lMatrix<double>& other ) const = 0;
            virtual void CastTo ( lMatrix<float>& other ) const = 0;

            /// convert a matrix from other type
            virtual void ConvertFrom ( const lMatrix<ValueType> &mat2 ) = 0;

            /// Copy the structure of a matrix
            virtual void CopyStructureFrom ( const lMatrix<ValueType> &mat2 ) = 0;
            virtual void CopyStructureTo ( lMatrix<ValueType> &mat2 ) const = 0;

            /// Clone a matrix - clean(), init(), initstructure(), copyfrom()
            virtual void CloneFrom ( const lMatrix<ValueType> &mat2 );

            /// Clone without content of the class - only initialization of the matrix is done
            /// @return the same class as this one
            virtual lMatrix<ValueType> *CloneWithoutContent ( void ) const;

            /// Synchronize (sync threads in CUDA, OpenCL, etc)
            virtual void Sync ( void ) const;

            /// Read the matrix from file
            /// @param filename - input file name
            virtual void ReadFile ( const char* filename );

            /// Write the matrix to file
            /// @param filename - output file name
            virtual void WriteFile ( const char* filename ) const;

            /// Matrix-Vector Multiplication y=Ax
            /// @param invec - the input vector (x)
            /// @return outvec - the output vector (y)
            virtual void VectorMult ( const lVector<ValueType> &invec,
                                      lVector<ValueType> *outvec ) const = 0;
            /// Add Matrix-Vector Multiplication y=y+Ax
            /// @param invec - the input vector (x)
            /// @return outvec - the output vector (y)
            virtual void VectorMultAdd ( const lVector<ValueType> &invec,
                                         lVector<ValueType> *outvec ) const = 0;

            virtual void VectorMultNoDiag ( const lVector<ValueType> &in,
                                            lVector<ValueType> *out ) const = 0;

            /// Compute the residual res=b-Ax
            /// @param b - the rhs
            /// @param x - the approxmiation
            /// @return res - the resitual = b-Ax
            virtual void Residual ( const lVector<ValueType> &b,
                                    const lVector<ValueType> &x,
                                    lVector<ValueType> *res ) const;

            /// Matrix-Matrix Multiplication C=this*B;
            /// the pointer C should be not allocated!
            /// @param inmat - the input matrix (B)
            /// @return outmat - the output matrix (C)
            virtual lMatrix<ValueType> *MatrixMult ( const lMatrix<ValueType> &inmat ) const;

            /// Compute the matrix pattern of the matrix product $|this|*|B|$
            /// This works only for matrices with symmetric pattern !!!
            /// @param inmat - the input matrix (B)
            /// @return outmat - the output matrix pattern (vals = 1.0)
            virtual lMatrix<ValueType> *MatrixMultSupStructure ( const lMatrix<ValueType> &inmat ) const;

            /// Compute the matrix pattern of $|this|^p$
            /// @param p - interger power > 0
            /// @return outmat - the output matrix pattern (vals = 1.0)
            virtual hiflow::la::lMatrix<ValueType> *MatrixSupSPower ( const int p ) const;

            /// Add a matrix to the current matrix this = this+B;
            /// The Matrix pattern should be the same or
            /// a subset!
            /// @param inmat - the second matrix (B)
            virtual void MatrixAdd ( const lMatrix<ValueType> &inmat );

            /// Gershgorin circles - provides a bound of the spectrum
            /// of a square matrix.
            /// @return lambda_min - the min eigenvalue approxmiation
            /// @return lambda_max - the max eigenvalue approxmiation
            virtual void GershgorinSpectrum ( ValueType *lambda_min, ValueType *lambda_max ) const;

            /// Perform Jacobi Precoditioner (Solve Mz=r)
            /// @param invec - input vector (z)
            /// @return outvec - output vector (r)
            virtual void Pjacobi ( const lVector<ValueType> &invec,
                                   lVector<ValueType> *outvec ) const;

            /// Perform Gauss-Sedel Precoditioner (Solve Mz=r)
            /// @param invec - input vector (z)
            /// @return outvec - output vector (r)
            virtual void Pgauss_seidel ( const lVector<ValueType> &invec,
                                         lVector<ValueType> *outvec ) const;

            /// Perform Symmetric Gauss-Seidel Precoditioner (Solve Mz=r)
            /// @param invec - input vector (z)
            /// @return outvec - output vector (r)
            virtual void Psgauss_seidel ( const lVector<ValueType> &invec,
                                          lVector<ValueType> *outvec ) const;

            /// Perform Block Symmetric Gauss-Seidel Precoditioner (Solve Mz=r)
            /// @param start_i - beginning of the block
            /// @param end_i - ending of the block
            /// @param invec - input vector (z)
            /// @return outvec - output vector (r)
            virtual void BlockPsgauss_seidel ( const lVector<ValueType> &invec,
                                               lVector<ValueType> *outvec,
                                               const int start_i,
                                               const int end_i ) const;

            /// Perform N Blocks Symmetric Gauss-Seidel Precoditioner (Solve Mz=r)
            /// If a parallel implemation is available then the blocks are perform
            /// in a parallel manner.
            /// @param num_blocks - number of blocks (degree of parallelism)
            /// @param invec - input vector (z)
            /// @return outvec - output vector (r)
            virtual void BlocksPsgauss_seidel ( const lVector<ValueType> &invec,
                                                lVector<ValueType> *outvec,
                                                const int num_blocks ) const;

            /// Perform SOR Precoditioner (Solve Mz=r)
            /// @param omega - relaxation parameter
            /// @param invec - input vector (z)
            /// @return outvec - output vector (r)
            virtual void Psor ( const ValueType omega,
                                const lVector<ValueType> &invec,
                                lVector<ValueType> *outvec ) const;

            /// Perform Symmetric SOR Precoditioner (Solve Mz=r)
            /// @param omega - relaxation parameter
            /// @param invec - input vector (z)
            /// @return outvec - output vector (r)
            virtual void Pssor ( const ValueType omega,
                                 const lVector<ValueType> &invec,
                                 lVector<ValueType> *outvec ) const;

          protected:

            int num_col_, num_row_, nnz_;

            enum MATRIX_FORMAT matrix_format_id_; // matrix ID
            enum IMPLEMENTATION implementation_id_; //  implementation ID
            enum PLATFORM platform_id_; // platform ID

            std::string matrix_format_name_;
            std::string name_;
            std::string platform_name_;
            std::string implementation_name_;

            bool symmetric_;

        };

        template <typename ValueType>
        inline int lMatrix<ValueType>::get_num_col ( void ) const
        {
            return this->num_col_;
        }

        template <typename ValueType>
        inline int lMatrix<ValueType>::get_num_row ( void ) const
        {
            return this->num_row_;
        }

        template <typename ValueType>
        inline int lMatrix<ValueType>::get_nnz ( void ) const
        {
            return this->nnz_;
        }

        template <typename ValueType>
        inline void lMatrix<ValueType>::swap_dimensions ( void )
        {
            std::swap ( num_row_, num_col_ );
        }

    } // namespace la
} // namespace hiflow

#endif
