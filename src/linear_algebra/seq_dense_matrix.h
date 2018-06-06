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

/// @author Chandramowli Subramanian, Martin Wlotzka, Simon Gawlok, Volker Lange

#ifndef HIFLOW_LINEARALGEBRA_SEQ_DENSE_MATRIX_H_
#    define HIFLOW_LINEARALGEBRA_SEQ_DENSE_MATRIX_H_

#    include <vector>
#    include <cassert>
#    include <iostream>
#    include <cstdlib>

#    include "lmp/la_global.h"
#    include "lmp/lvector_cpu.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Sequential dense matrix
        ///
        /// Sequential dense m x n matrix. Entries are stored in one array.

        template<class DataType>
        class SeqDenseMatrix
        {
          public:

            /// Standard constructor
            SeqDenseMatrix ( );
            /// Destructor.
            ~SeqDenseMatrix ( );

            /// Clears allocated memory.
            void Clear ( );

            /// Resizes the matrix.
            /// @param m number of rows
            /// @param n number of columns
            void Resize ( int m, int n );

            /// Sets all entries of the matrix to zero.
            void Zeros ( );

            inline const DataType& operator() ( int i, int j ) const;
            inline DataType& operator() ( int i, int j );

            /// @return number of rows

            int nrows ( ) const
            {
                return this->nrows_;
            }
            /// @return number of columns

            int ncols ( ) const
            {
                return this->ncols_;
            }

            /// Matrix-Matrix mulitplication
            /// @param inMat matrix to be multiplied with the matrix
            /// @param outMat resulting matrix. outMat = this * inMat
            void MatrixMult ( const SeqDenseMatrix<DataType> &inMat, SeqDenseMatrix<DataType> &outMat );

            /// @return eigenvalues in ascending order

            std::vector<DataType> get_eigenvalues ( ) const
            {
                return this->eigenvalues_;
            }

            /// @return eigenvectors stored one after another

            std::vector<DataType> get_eigenvectors ( ) const
            {
                return this->eigenvectors_;
            }

            /// Compute eigenvalues and eigenvectors if the matrix is symmetric
            void compute_eigenvalues_and_vectors ( );

            /// Matrix-Vector multiplication
            /// @param invec vector to be multiplied with matrix
            /// @param outvec result vector
            void VectorMult ( const CPU_lVector<DataType> &invec, CPU_lVector<DataType> &outvec );

            /// Matrix-Vector multiplication
            /// @param invec vector to be multiplied with matrix
            /// @param outvec result vector
            void VectorMult ( const std::vector<DataType> &invec, std::vector<DataType> &outvec );

            /// LU factorization
            void Factorize ( );

            /// QR factorization
            /// @param[out] Q Q (orthogonal) matrix of factorization
            /// @param[out] R R (upper triangular) matrix of factorization
            void QRFactorize ( SeqDenseMatrix<DataType> &Q, SeqDenseMatrix<DataType> &R );

            /// Reduction of this matrix to Hessenberg matrix
            /// @param[out] H Hessenberg reduction of this matrix
            /// @param[out] V Matrix storing similarity transformation
            void HessenbergReduction ( SeqDenseMatrix<DataType> &H, SeqDenseMatrix<DataType> &V );

            /// Forward/Backward solve linear system Ax=b; need to run Factorize() or Solve() prior to execution
            /// @param x result vector
            /// @param b right hand side vector
            void ForwardBackward ( const lVector<DataType> &b, lVector<DataType> &x );

            /// Forward/Backward solve linear system Ax=b; need to run Factorize() or Solve() prior to execution
            /// @param x result vector
            /// @param b right hand side vector
            void ForwardBackward ( const std::vector<DataType> &b, std::vector<DataType> &x );

            /// Solve linear system Ax=b
            /// @param x result vector
            /// @param b right hand side vector
            void Solve ( const lVector<DataType> &b, lVector<DataType> &x );

            /// Solve linear system Ax=b
            /// @param x result vector
            /// @param b right hand side vector
            void Solve ( const std::vector<DataType> &b, std::vector<DataType> &x );

            /// Add Matrix entrywise. A_ij := A_ij + B_ij
            /// @param B Matrix of same dimensions, that is to be added.
            void Add ( const SeqDenseMatrix<DataType> &B );

            /// scale Matrix entrywise. A_ij := alpha * A_ij
            /// @param alpha scalar
            void Scale ( const DataType alpha );

            /// Add a scaled Matrix entrywise. A_ij := A_ij + alpha * B_ij
            /// @param B Matrix of same dimensions, that is to be added.
            /// @param alpha scalar for Matrix B
            void Axpy ( const SeqDenseMatrix<DataType> &B, const DataType alpha );

            /// Add smaller Matrix entrywise at position ij. A_(i+k)(j+l) := A_(i+k)(j+l) + B_kl
            /// @param B Matrix of smaller dimensions, that is to be added.
            void Add ( const SeqDenseMatrix<DataType> &B, const int i, const int j );

            /// Insert smaller Matrix at position ij. A_(i+k)(j+l) := B_kl
            /// @param B Matrix of smaller dimensions, that is to be inserted.
            void ins ( const SeqDenseMatrix<DataType> &B, const int i, const int j );

            /// Extract part of Matrix from position ij to kl.
            /// @param B Matrix of smaller dimensions, that is to be inserted.
            void extr ( SeqDenseMatrix<DataType> &B, const int i, const int j, const int k, const int l );

            /// Transposes the Matrix.
            /// @param trans_mat transposed Matrix.
            void transpose_me ( SeqDenseMatrix<DataType> &trans_mat );

            /// Frobenius norm auf matrix
            DataType Frob ( void );

            /// Set block size for LU decomposition
            /// @param size block size

            inline void set_blocksize ( int size )
            {
                this->block_ = size;
            }

            /// Get block size for LU decomposition

            inline int get_blocksize ( )
            {
                return this->block_;
            }

            /// @return Implementation ID

            enum IMPLEMENTATION get_implementation ( void ) const
            {
                return implementation_id_;
            }

            /// @return Platform ID

            enum PLATFORM get_platform ( void ) const
            {
                return platform_id_;
            }

            /// @return Implemenation name

            std::string get_implementation_name ( void ) const
            {
                return implementation_name_;
            }

            /// @return Platform name

            std::string get_platform_name ( void ) const
            {
                return platform_name_;
            }

            /// Set implementation
            /// @param impl Implementation ID
            /// @param plat Platform ID
            void set_implementation_and_platform ( enum IMPLEMENTATION impl, enum PLATFORM plat );

          protected:
            std::vector<DataType> val_; // entries row by row
            int nrows_; // number of rows
            int ncols_; // number of columns

            std::vector<DataType> eigenvalues_; //eigenvalues
            std::vector<DataType> eigenvectors_; //eigenvectors

            std::vector<DataType> lu_; // lu decomposition
            std::vector<int> permutation_; // permutation vector of rows in LU Decomposition
            int block_; // block size for LU decomposition

            std::string platform_name_; // platform name
            std::string implementation_name_; // implementation specification
            enum IMPLEMENTATION implementation_id_; //  implementation ID
            enum PLATFORM platform_id_; // platform ID

          private:
            /// Compute eigenvalues and eigenvectors with BLAS/LAPACK routine
            void compute_eigenvalues_and_vectors_cblas ( );
            /// Compute eigenvalues and eigenvectors with MKL routine
            void compute_eigenvalues_and_vectors_mkl ( );
            /// Compute eigenvalues and eigenvectors with simple routine
            void compute_eigenvalues_and_vectors_simple ( );
            /// Compute eigenvalues and eigenvectors with OpenMP parallelized routine
            void compute_eigenvalues_and_vectors_openmp ( );

            /// Compute Matrix-Matrix product with BLAS/LAPACK routine
            void MatrixMult_cblas ( const SeqDenseMatrix<DataType> *inMat, SeqDenseMatrix<DataType> *outMat );
            /// Compute Matrix-Matrix product with MKL routine
            void MatrixMult_mkl ( const SeqDenseMatrix<DataType> *inMat, SeqDenseMatrix<DataType> *outMat );
            /// Compute Matrix-Matrix product with simple routine
            void MatrixMult_simple ( const SeqDenseMatrix<DataType> *inMat, SeqDenseMatrix<DataType> *outMat );
            /// Compute Matrix-Matrix product with OpenMP routine
            void MatrixMult_openmp ( const SeqDenseMatrix<DataType> *inMat, SeqDenseMatrix<DataType> *outMat );

            /// Compute Matrix-Vector product with BLAS routine
            void VectorMult_cblas ( const DataType *invec, DataType *outvec );
            /// Compute Matrix-Vector product with MKL routine
            void VectorMult_mkl ( const DataType *invec, DataType *outvec );
            /// Compute Matrix-Vector product with simple routine
            void VectorMult_simple ( const DataType *invec, DataType *outvec );
            /// Compute Matrix-Vector product with OpenMP routine
            void VectorMult_openmp ( const DataType *invec, DataType *outvec );

            /// Factorize linear system with BLAS routine
            void Factorize_cblas ( );
            /// Factorize linear system with MKL routine
            void Factorize_mkl ( );
            /// Factorize linear system with simple routine
            void Factorize_simple ( );
            /// Factorize linear system with OpenMP routine
            void Factorize_openmp ( );

            /// ForwardBackward with BLAS routine
            void ForwardBackward_cblas ( DataType *x );
            /// ForwardBackward with MKL routine
            void ForwardBackward_mkl ( DataType *x );
            /// ForwardBackward with simple routine
            void ForwardBackward_simple ( DataType *x );

            /// QR factorization with simple routine
            void QRFactorize_simple ( SeqDenseMatrix<DataType> &Q, SeqDenseMatrix<DataType> &R );
            /// QR factorization with OpenMP routine
            void QRFactorize_openmp ( SeqDenseMatrix<DataType> &Q, SeqDenseMatrix<DataType> &R );
            /// QR factorization with MKL routine
            void QRFactorize_mkl ( SeqDenseMatrix<DataType> &Q, SeqDenseMatrix<DataType> &R );
            /// QR factorization with BLAS routine
            void QRFactorize_cblas ( SeqDenseMatrix<DataType> &Q, SeqDenseMatrix<DataType> &R );

            /// Hessenberg reduction with simple routine
            void HessenbergReduction_simple ( SeqDenseMatrix<DataType> &H, SeqDenseMatrix<DataType> &V );
            /// Hessenberg reduction with OpenMP routine
            void HessenbergReduction_openmp ( SeqDenseMatrix<DataType> &H, SeqDenseMatrix<DataType> &V );
            /// Hessenberg reduction with MKL routine
            void HessenbergReduction_mkl ( SeqDenseMatrix<DataType> &H, SeqDenseMatrix<DataType> &V );
            /// Hessenberg reduction with BLAS routine
            void HessenbergReduction_cblas ( SeqDenseMatrix<DataType> &H, SeqDenseMatrix<DataType> &V );

            /// simple gemm routine
            void gemm_simple ( const int i, const int blockend );
            /// openmp gemm routine
            void gemm_openmp ( const int i, const int blockend );

            template<class T> friend std::ostream& operator<< ( std::ostream& os, const SeqDenseMatrix<T>& m );

            // no implementation of copy constructor or assignement operator
            SeqDenseMatrix ( const SeqDenseMatrix<DataType>& );
            SeqDenseMatrix<DataType>& operator= ( const SeqDenseMatrix& );
        };

        /// Returns the (i,j)-th element of the matrix.
        /// @param i row
        /// @param j column
        /// @return (i,j)-th element

        template<class DataType>
        inline
        const DataType& SeqDenseMatrix<DataType>::operator() ( int i, int j ) const
        {
            assert ( i >= 0 );
            assert ( i < this->nrows ( ) );
            assert ( j >= 0 );
            assert ( j < this->ncols ( ) );
            return this->val_[j + i * ncols_];
        }

        /// Returns the (i,j)-th element of the matrix.
        /// @param i row
        /// @param j column
        /// @return (i,j)-th element

        template<class DataType>
        inline DataType& SeqDenseMatrix<DataType>::operator() ( int i, int j )
        {
            assert ( i >= 0 );
            assert ( i < this->nrows ( ) );
            assert ( j >= 0 );
            assert ( j < this->ncols ( ) );
            return this->val_[j + i * ncols_];
        }

        template<class T>
        std::ostream& operator<< ( std::ostream& os, const SeqDenseMatrix<T>& m )
        {
            os << "[";
            for ( int i = 0; i < m.nrows ( ); ++i )
            {
                for ( int j = 0; j < m.ncols ( ); ++j )
                {
                    os << m ( i, j ) << " ";
                }
                os << "\n";
            }
            os << "]";
            return os;
        }

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_SEQ_DENSE_MATRIX_H_
