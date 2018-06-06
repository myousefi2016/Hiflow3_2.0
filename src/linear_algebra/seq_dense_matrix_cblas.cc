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

#include "config.h"

#include "seq_dense_matrix.h"

#include "lmp/lmp_log.h"

#ifdef WITH_CLAPACK
#    include <lapacke.h>
extern "C"
{
#    include <cblas.h>
}
#else
#    define ERROR LOG_ERROR("no Clapack support");  exit(-1);
#endif

namespace hiflow
{
    namespace la
    {

        template<>
        void SeqDenseMatrix<double>::compute_eigenvalues_and_vectors_cblas ( )
        {
#ifdef WITH_CLAPACK
            //prepare working arguments
            lapack_int info = 0;

            //convert dimensions to appropriate format
            lapack_int n = static_cast < lapack_int > ( ncols_ );

            //define kind of job to do
            char jobz = 'V';
            char joba = 'U';

            //computation
            info = LAPACKE_dsyevd ( LAPACK_ROW_MAJOR, jobz, joba, n, &( eigenvectors_[0] ), n, &( eigenvalues_[0] ) );

            //check correctness
            if ( info < 0 )
            {
                LOG_ERROR ( "argument number " << -info << " had an illegal value" );
                exit ( -1 );
            }
            else if ( info > 0 )
            {
                LOG_ERROR ( "algorithm failed to compute an eigenvalue" );
                exit ( -1 );
            }
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<float>::compute_eigenvalues_and_vectors_cblas ( )
        {
#ifdef WITH_CLAPACK
            //prepare working arguments
            lapack_int info = 0;

            //convert dimensions to appropriate format
            lapack_int n = static_cast < lapack_int > ( ncols_ );

            //define kind of job to do
            char jobz = 'V';
            char joba = 'U';

            //computation
            info = LAPACKE_ssyevd ( LAPACK_ROW_MAJOR, jobz, joba, n, &eigenvectors_[0], n, &eigenvalues_[0] );

            //check correctness
            if ( info < 0 )
            {
                LOG_ERROR ( "argument number " << -info << " had an illegal value" );
                exit ( -1 );
            }
            else if ( info > 0 )
            {
                LOG_ERROR ( "algorithm failed to compute an eigenvalue" );
                exit ( -1 );
            }
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<double>::VectorMult_cblas ( const double *invec, double *outvec )
        {
#ifdef WITH_CLAPACK
            cblas_dgemv ( CblasRowMajor, CblasNoTrans, nrows_, ncols_, 1.0, &val_[0], ncols_, invec, 1, 0.0, outvec, 1 );
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<float>::VectorMult_cblas ( const float *invec, float *outvec )
        {
#ifdef WITH_CLAPACK
            cblas_sgemv ( CblasRowMajor, CblasNoTrans, nrows_, ncols_, 1.0, &val_[0], ncols_, invec, 1, 0.0, outvec, 1 );
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<double>::MatrixMult_cblas ( const SeqDenseMatrix<double> *inMat, SeqDenseMatrix<double> *outMat )
        {
#ifdef WITH_CLAPACK
            cblas_dgemm ( CblasRowMajor, CblasNoTrans, CblasNoTrans, nrows_, outMat->ncols ( ), ncols_, 1.0, &val_[0], ncols_, &( inMat->val_[0] ), ncols_, 0.0, &( outMat->val_[0] ), outMat->ncols ( ) );
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<float>::MatrixMult_cblas ( const SeqDenseMatrix<float> *inMat, SeqDenseMatrix<float> *outMat )
        {
#ifdef WITH_CLAPACK
            cblas_sgemm ( CblasRowMajor, CblasNoTrans, CblasNoTrans, nrows_, outMat->ncols ( ), ncols_, 1.0, &val_[0], ncols_, &( inMat->val_[0] ), ncols_, 0.0, &( outMat->val_[0] ), outMat->ncols ( ) );
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<double>::Factorize_cblas ( )
        {
#ifdef WITH_CLAPACK
            lapack_int n = static_cast < lapack_int > ( nrows_ );

            lapack_int *perm = new lapack_int[nrows_];

            lapack_int info = 0;
            info = LAPACKE_dgetrf ( LAPACK_ROW_MAJOR, n, n, &( lu_[0] ), n, perm );

            //check for correct determination
            if ( info > 0 )
            {
                LOG_ERROR ( "The diagonal element of the triangular factor of A, U(" << info << "," << info << ") is zero, so that A is singular; the solution could not be computed" );
                exit ( -1 );
            }

            // copy perm to permutation_
            for ( int i = 0; i < nrows_; ++i )
            {
                permutation_[i] = static_cast < int > ( perm[i] );
            }

            delete [] perm;

#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<float>::Factorize_cblas ( )
        {
#ifdef WITH_CLAPACK
            lapack_int n = static_cast < lapack_int > ( nrows_ );

            lapack_int *perm = new lapack_int[nrows_];

            lapack_int info = 0;
            info = LAPACKE_sgetrf ( LAPACK_ROW_MAJOR, n, n, &( lu_[0] ), n, perm );

            //check for correct determination
            if ( info > 0 )
            {
                LOG_ERROR ( "The diagonal element of the triangular factor of A, U(" << info << "," << info << ") is zero, so that A is singular; the solution could not be computed" );
                exit ( -1 );
            }

            // copy perm to permutation_
            for ( int i = 0; i < nrows_; ++i )
            {
                permutation_[i] = static_cast < int > ( perm[i] );
            }

            delete [] perm;
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<double>::QRFactorize_cblas ( SeqDenseMatrix<double> &Q, SeqDenseMatrix<double> &R )
        {
#ifdef WITH_CLAPACK
            // copy this to Q
            for ( int i = 0; i < nrows_; ++i )
            {
                const int offset = i * ncols_;
                for ( int j = 0; j < ncols_; ++j )
                {
                    Q ( i, j ) = val_[offset + j];
                }
            }

            // clear R
            R.Zeros ( );

            double *tau = new double[nrows_];
            lapack_int info = 0;
            info = LAPACKE_dgeqrf ( LAPACK_ROW_MAJOR, nrows_, ncols_, &Q ( 0, 0 ), ncols_, tau );

            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "Argument " << -info << " of sgeqrf has an illegal value!" );
                exit ( -1 );
            }

            info = 0;
            info = LAPACKE_dorgqr ( LAPACK_ROW_MAJOR, nrows_, ncols_, ncols_, &Q ( 0, 0 ), ncols_, tau );
            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "Argument " << -info << " of sorgqr has an illegal value!" );
                exit ( -1 );
            }

            // compute R = Q^T * this
            for ( int i = 0; i < ncols_; ++i )
            {
                for ( int k = 0; k < nrows_; ++k )
                {
                    const int k_index = k * ncols_;
                    const double temp = Q ( k, i );
                    for ( int j = i; j < ncols_; ++j )
                    {
                        R ( i, j ) += temp * val_[k_index + j];
                    }
                }
            }

            delete [] tau;
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<float>::QRFactorize_cblas ( SeqDenseMatrix<float> &Q, SeqDenseMatrix<float> &R )
        {
#ifdef WITH_CLAPACK
            // copy this to Q
            for ( int i = 0; i < nrows_; ++i )
            {
                const int offset = i * ncols_;
                for ( int j = 0; j < ncols_; ++j )
                {
                    Q ( i, j ) = val_[offset + j];
                }
            }

            // clear R
            R.Zeros ( );

            float *tau = new float[nrows_];
            lapack_int info = 0;
            info = LAPACKE_sgeqrf ( LAPACK_ROW_MAJOR, nrows_, ncols_, &Q ( 0, 0 ), ncols_, tau );

            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "Argument " << -info << " of sgeqrf has an illegal value!" );
                exit ( -1 );
            }

            info = 0;
            info = LAPACKE_sorgqr ( LAPACK_ROW_MAJOR, nrows_, ncols_, ncols_, &Q ( 0, 0 ), ncols_, tau );
            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "Argument " << -info << " of sorgqr has an illegal value!" );
                exit ( -1 );
            }

            // compute R = Q^T * this
            for ( int i = 0; i < ncols_; ++i )
            {
                for ( int k = 0; k < nrows_; ++k )
                {
                    const int k_index = k * ncols_;
                    const float temp = Q ( k, i );
                    for ( int j = i; j < ncols_; ++j )
                    {
                        R ( i, j ) += temp * val_[k_index + j];
                    }
                }
            }

            delete [] tau;
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<double>::HessenbergReduction_cblas ( SeqDenseMatrix<double> &H, SeqDenseMatrix<double> &V )
        {
#ifdef WITH_CLAPACK

            // copy this to V
            for ( int i = 0; i < nrows_; ++i )
            {
                const int offset = i * ncols_;
                for ( int j = 0; j < ncols_; ++j )
                {
                    V ( i, j ) = val_[offset + j];
                }
            }

            // clear H
            H.Zeros ( );

            double *tau = new double[nrows_];
            lapack_int info = 0;

            info = LAPACKE_dgehrd ( LAPACK_ROW_MAJOR, nrows_, 1, nrows_, &V ( 0, 0 ), nrows_, tau );

            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "Argument " << -info << " of dgehrd has an illegal value!" );
                exit ( -1 );
            }

            // copy Hessenberg part from V to H
            for ( int i = 0; i < nrows_; ++i )
            {
                for ( int j = std::max ( 0, i - 1 ); j < ncols_; ++j )
                {
                    H ( i, j ) = V ( i, j );
                }
            }

            // Compute Matrix V
            info = 0;
            info = LAPACKE_dorghr ( LAPACK_ROW_MAJOR, nrows_, 1, nrows_, &V ( 0, 0 ), nrows_, tau );
            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "Argument " << -info << " of sorgqr has an illegal value!" );
                exit ( -1 );
            }

            delete [] tau;
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<float>::HessenbergReduction_cblas ( SeqDenseMatrix<float> &H, SeqDenseMatrix<float> &V )
        {
#ifdef WITH_CLAPACK

            // copy this to V
            for ( int i = 0; i < nrows_; ++i )
            {
                const int offset = i * ncols_;
                for ( int j = 0; j < ncols_; ++j )
                {
                    V ( i, j ) = val_[offset + j];
                }
            }

            // clear H
            H.Zeros ( );

            float *tau = new float[nrows_];
            lapack_int info = 0;

            info = LAPACKE_sgehrd ( LAPACK_ROW_MAJOR, nrows_, 1, nrows_, &V ( 0, 0 ), nrows_, tau );

            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "Argument " << -info << " of dgehrd has an illegal value!" );
                exit ( -1 );
            }

            // copy Hessenberg part from V to H
            for ( int i = 0; i < nrows_; ++i )
            {
                for ( int j = std::max ( 0, i - 1 ); j < ncols_; ++j )
                {
                    H ( i, j ) = V ( i, j );
                }
            }

            // Compute Matrix V
            info = 0;
            info = LAPACKE_sorghr ( LAPACK_ROW_MAJOR, nrows_, 1, nrows_, &V ( 0, 0 ), nrows_, tau );
            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "Argument " << -info << " of sorgqr has an illegal value!" );
                exit ( -1 );
            }

            delete [] tau;
#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<double>::ForwardBackward_cblas ( double *x )
        {
#ifdef WITH_CLAPACK
            //setup parameters
            lapack_int info = 0;
            lapack_int n = static_cast < lapack_int > ( nrows_ );

            lapack_int nrhs = 1;

            lapack_int *perm = new lapack_int[nrows_];
            for ( int i = 0; i < nrows_; ++i )
            {
                perm[i] = static_cast < lapack_int > ( permutation_[i] );
            }

            //perform solving
            info = LAPACKE_dgetrs ( LAPACK_ROW_MAJOR, 'N', n, nrhs, &( lu_[0] ), n, perm, x, nrhs );

            //check for correct determination
            if ( info < 0 )
            {
                LOG_ERROR ( "The parameter " << -info << ") has an illegal value. The solution could not be computed!" );
                exit ( -1 );
            }

            delete [] perm;

#else
            ERROR;
#endif
        }

        template<>
        void SeqDenseMatrix<float>::ForwardBackward_cblas ( float *x )
        {
#ifdef WITH_CLAPACK
            //setup parameters
            lapack_int info = 0;
            lapack_int n = static_cast < lapack_int > ( nrows_ );

            lapack_int nrhs = 1;

            lapack_int *perm = new lapack_int[nrows_];
            for ( int i = 0; i < nrows_; ++i )
            {
                perm[i] = static_cast < lapack_int > ( permutation_[i] );
            }

            //perform solving
            info = LAPACKE_sgetrs ( LAPACK_ROW_MAJOR, 'N', n, nrhs, &( lu_[0] ), n, perm, x, nrhs );

            //check for correct determination
            if ( info > 0 )
            {
                LOG_ERROR ( "The diagonal element of the triangular factor of A, U(" << info << "," << info << ") is zero, so that A is singular; the solution could not be computed" );
                exit ( -1 );
            }

            delete [] perm;

#else
            ERROR;
#endif
        }

    } // namespace la
} // namespace hiflow
