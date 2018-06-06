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

#include "seq_dense_matrix.h"

#include "common/pointers.h"

#include "lmp/lmp_log.h"
#include "lmp/lmp_mem.h"

#ifdef WITH_OPENMP
#    include <omp.h>
#    define ERROR LOG_ERROR("no OpenMP support");  exit(-1);
#else
#    define ERROR LOG_ERROR("no OpenMP support");  exit(-1);
#endif

namespace hiflow
{

    namespace la
    {

        template<class DataType>

        SeqDenseMatrix<DataType>::SeqDenseMatrix ( )
        {

            this->Clear ( );

            this->platform_name_ = "CPU (x86)";

            this->platform_id_ = CPU;

            this->implementation_name_ = "NAIVE";
            this->implementation_id_ = NAIVE;

#ifdef WITH_OPENMP
            this->implementation_name_ = "OPENMP";
            this->implementation_id_ = OPENMP;
#endif

#ifdef WITH_CLAPACK
            this->implementation_name_ = "BLAS";
            this->implementation_id_ = BLAS;
#endif

#ifdef WITH_MKL
            this->implementation_name_ = "Intel MKL";
            this->implementation_id_ = MKL;
#endif

            this->block_ = 16;

        }

        template<class DataType>

        SeqDenseMatrix<DataType>::~SeqDenseMatrix ( )
        {

            this->Clear ( );

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Clear ( )
        {

            this->val_.clear ( );

            this->nrows_ = 0;

            this->ncols_ = 0;

            eigenvalues_.clear ( );

            eigenvectors_.clear ( );

            lu_.clear ( );

            permutation_.clear ( );

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Resize ( int m, int n )
        {

            assert ( m >= 0 && n >= 0 );

            this->Clear ( );

            if ( m > 0 && n > 0 )
            {

                this->val_.resize ( m * n, static_cast < DataType > ( 0.0 ) );

                this->nrows_ = m;

                this->ncols_ = n;

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Zeros ( )
        {

            const int e = this->val_.size ( );

            this->val_.clear ( );

            this->val_.resize ( e, static_cast < DataType > ( 0.0 ) );

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::set_implementation_and_platform ( enum IMPLEMENTATION impl, enum PLATFORM plat )
        {

            if ( plat == CPU )
            {

                switch ( impl )
                {

                    case NAIVE:

                        platform_id_ = plat;

                        platform_name_ = "CPU (x86)";

                        implementation_id_ = impl;

                        implementation_name_ = "NAIVE";

                        return;

                        break;

                    case OPENMP:

                        platform_id_ = plat;

                        platform_name_ = "CPU (x86)";

                        implementation_id_ = impl;

                        implementation_name_ = "OPENMP";

                        return;

                        break;

                    case MKL:

                        platform_id_ = plat;

                        platform_name_ = "CPU (x86)";

                        implementation_id_ = impl;

                        implementation_name_ = "Intel MKL";

                        return;

                        break;

                    case BLAS:

                        platform_id_ = plat;

                        platform_name_ = "CPU (x86)";

                        implementation_id_ = impl;

                        implementation_name_ = "BLAS";

                        return;

                        break;

                    default:

                        ;

                }

            }

            LOG_ERROR ( "set_implementation_and_platform() incompatibility PLATFORM/IMPLEMENTATION" );

            LOG_ERROR ( " Platform ID=" << plat << " Implementation ID=" << impl );

            exit ( -1 );

            return;

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Scale ( const DataType alpha )
        {

            const int e = this->val_.size ( );

            switch ( implementation_id_ )
            {

                case OPENMP:

#ifdef WITH_OPENMP

#    pragma omp parallel for

                    for ( int i = 0; i < e; ++i )
                    {

                        val_[i] = alpha * val_[i];

                    }

                    return;

                    break;

#endif

                default:

                    for ( int i = 0; i < e; ++i )
                    {

                        val_[i] = alpha * val_[i];

                    }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Add ( const SeqDenseMatrix<DataType> &B )
        {

            assert ( nrows_ == B.nrows ( ) );

            assert ( ncols_ == B.ncols ( ) );

            const DataType *b = &B ( 0, 0 );

            const int e = this->val_.size ( );

            switch ( implementation_id_ )
            {

                case OPENMP:

#ifdef WITH_OPENMP

#    pragma omp parallel for

                    for ( int i = 0; i < e; ++i )
                    {

                        val_[i] += b[i];

                    }

                    return;

                    break;

#endif

                default:

                    for ( int i = 0; i < e; ++i )
                    {

                        val_[i] += b[i];

                    }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Axpy ( const SeqDenseMatrix<DataType> &B, const DataType alpha )
        {

            assert ( nrows_ == B.nrows ( ) );

            assert ( ncols_ == B.ncols ( ) );

            const DataType *b = &B ( 0, 0 );

            const int e = this->val_.size ( );

            switch ( implementation_id_ )
            {

                case OPENMP:

#ifdef WITH_OPENMP

#    pragma omp parallel for

                    for ( int i = 0; i < e; ++i )
                    {

                        val_[i] += alpha * b[i];

                    }

                    return;

                    break;

#endif

                default:

                    for ( int i = 0; i < e; ++i )
                    {

                        val_[i] += alpha * b[i];

                    }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Add ( const SeqDenseMatrix<DataType> &B, const int row, const int col )
        {

            assert ( row >= 0 );

            assert ( col >= 0 );

            assert ( nrows_ >= ( B.nrows ( ) + row ) );

            assert ( ncols_ >= ( B.ncols ( ) + col ) );

            const DataType *b = &B ( 0, 0 );

            switch ( implementation_id_ )
            {

                case OPENMP:

#ifdef WITH_OPENMP

#    pragma omp parallel for

                    for ( int i = 0; i < B.nrows ( ); ++i )
                    {

                        const int b_index = i * B.ncols ( );

                        const int offset = col + ( i + row ) * ncols_;

                        for ( int j = 0; j < B.ncols ( ); ++j )
                        {

                            val_[offset + j] += b[b_index + j];

                        }

                    }

                    return;

                    break;

#endif

                default:

                    for ( int i = 0; i < B.nrows ( ); ++i )
                    {

                        const int b_index = i * B.ncols ( );

                        const int offset = col + ( i + row ) * ncols_;

                        for ( int j = 0; j < B.ncols ( ); ++j )
                        {

                            val_[offset + j] += b[b_index + j];

                        }

                    }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::ins ( const SeqDenseMatrix<DataType> &B, const int row, const int col )
        {

            assert ( row >= 0 );

            assert ( col >= 0 );

            assert ( nrows_ >= ( B.nrows ( ) + row ) );

            assert ( ncols_ >= ( B.ncols ( ) + col ) );

            const DataType *b = &B ( 0, 0 );

            switch ( implementation_id_ )
            {

                case OPENMP:

#ifdef WITH_OPENMP

#    pragma omp parallel for

                    for ( int i = 0; i < B.nrows ( ); ++i )
                    {

                        const int b_index = i * B.ncols ( );

                        const int offset = col + ( i + row ) * ncols_;

                        for ( int j = 0; j < B.ncols ( ); ++j )
                        {

                            val_[offset + j] = b[b_index + j];

                        }

                    }

                    return;

                    break;

#endif

                default:

                    for ( int i = 0; i < B.nrows ( ); ++i )
                    {

                        const int b_index = i * B.ncols ( );

                        const int offset = col + ( i + row ) * ncols_;

                        for ( int j = 0; j < B.ncols ( ); ++j )
                        {

                            val_[offset + j] = b[b_index + j];

                        }

                    }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::extr ( SeqDenseMatrix<DataType> &B, const int row1, const int col1, const int row2, const int col2 )
        {

            assert ( row1 >= 0 );

            assert ( col1 >= 0 );

            assert ( row2 > row1 );

            assert ( col2 > col1 );

            assert ( nrows_ >= row2 );

            assert ( ncols_ >= col2 );

            assert ( B.nrows ( ) == ( row2 - row1 + 1 ) );

            assert ( B.ncols ( ) == ( col2 - col1 + 1 ) );

            const int diffRow = row2 - row1;

            const int diffCol = col2 - col1;

            DataType *b = &B ( 0, 0 );

            switch ( implementation_id_ )
            {

                case OPENMP:

#ifdef WITH_OPENMP

#    pragma omp parallel for

                    for ( int i = 0; i <= diffRow; ++i )
                    {

                        const int b_index = i * B.ncols ( );

                        const int offset = col1 + ( i + row1 ) * ncols_;

                        for ( int j = 0; j <= diffCol; ++j )
                        {

                            b[b_index + j] = val_[offset + j];

                        }

                    }

                    return;

                    break;

#endif

                default:

                    for ( int i = 0; i <= diffRow; ++i )
                    {

                        const int b_index = i * B.ncols ( );

                        const int offset = col1 + ( i + row1 ) * ncols_;

                        for ( int j = 0; j <= diffCol; ++j )
                        {

                            b[b_index + j] = val_[offset + j];

                        }

                    }

            }

        }

        template<class Datatype>

        void SeqDenseMatrix<Datatype>::transpose_me ( SeqDenseMatrix<Datatype> &trans_mat )
        {

            assert ( ncols_ == trans_mat.nrows ( ) );

            assert ( nrows_ == trans_mat.ncols ( ) );

            switch ( implementation_id_ )
            {

                case OPENMP:

#ifdef WITH_OPENMP

#    pragma omp parallel for

                    for ( int i = 0; i < ncols_; ++i )
                    {

                        const int offset = i * nrows_;

                        for ( int j = i + 1; j < nrows_; ++j )
                        {

                            trans_mat.val_[offset + j] = val_[j * ncols_ + i];

                        }

                    }

                    return;

                    break;

#endif

                default:

                    for ( int i = 0; i < ncols_; ++i )
                    {

                        const int offset = i * nrows_;

                        for ( int j = i + 1; j < nrows_; ++j )
                        {

                            trans_mat.val_[offset + j] = val_[j * ncols_ + i];

                        }

                    }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::compute_eigenvalues_and_vectors ( )
        {

            //assert aquare matrix

            assert ( nrows_ == ncols_ );

            //initialize eigenvectors_

            eigenvectors_ = val_;

            //allocate eigenvalue vector

            eigenvalues_.clear ( );

            eigenvalues_.resize ( ncols_ );

            if ( platform_id_ == CPU )
            {

                switch ( implementation_id_ )
                {

                    case NAIVE:

                        compute_eigenvalues_and_vectors_simple ( );

                        return;

                        break;

                    case OPENMP:

                        compute_eigenvalues_and_vectors_openmp ( );

                        return;

                        break;

                    case MKL:

                        compute_eigenvalues_and_vectors_mkl ( );

                        return;

                        break;

                    case BLAS:

                        compute_eigenvalues_and_vectors_cblas ( );

                        return;

                        break;

                    default:

                        ;

                }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::MatrixMult ( const SeqDenseMatrix<DataType> &inMat, SeqDenseMatrix<DataType> &outMat )
        {

            // check for compatible dimensions

            assert ( this->ncols_ == inMat.nrows ( ) );

            assert ( this->nrows_ == outMat.nrows ( ) );

            assert ( inMat.ncols ( ) == outMat.ncols ( ) );

            if ( platform_id_ == CPU )
            {

                switch ( implementation_id_ )
                {

                    case NAIVE:

                        MatrixMult_simple ( &inMat, &outMat );

                        return;

                        break;

                    case OPENMP:

                        MatrixMult_openmp ( &inMat, &outMat );

                        return;

                        break;

                    case MKL:

                        MatrixMult_mkl ( &inMat, &outMat );

                        return;

                        break;

                    case BLAS:

                        MatrixMult_cblas ( &inMat, &outMat );

                        return;

                        break;

                    default:

                        ;

                }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::VectorMult ( const CPU_lVector<DataType> &invec, CPU_lVector<DataType> &outvec )
        {

            // check for compatible dimensions

            assert ( this->ncols_ == invec.get_size ( ) );

            assert ( this->nrows_ == outvec.get_size ( ) );

            if ( platform_id_ == CPU )
            {

                switch ( implementation_id_ )
                {

                    case NAIVE:

                        VectorMult_simple ( invec.buffer, outvec.buffer );

                        return;

                        break;

                    case OPENMP:

                        VectorMult_openmp ( invec.buffer, outvec.buffer );

                        return;

                        break;

                    case MKL:

                        VectorMult_mkl ( invec.buffer, outvec.buffer );

                        return;

                        break;

                    case BLAS:

                        VectorMult_cblas ( invec.buffer, outvec.buffer );

                        return;

                        break;

                    default:

                        ;

                }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::VectorMult ( const std::vector<DataType> &invec, std::vector<DataType> &outvec )
        {

            // check for compatible dimensions

            assert ( this->ncols_ == invec.size ( ) );

            assert ( this->nrows_ == outvec.size ( ) );

            if ( platform_id_ == CPU )
            {

                switch ( implementation_id_ )
                {

                    case NAIVE:

                        VectorMult_simple ( &invec[0], &outvec[0] );

                        return;

                        break;

                    case OPENMP:

                        VectorMult_openmp ( &invec[0], &outvec[0] );

                        return;

                        break;

                    case MKL:

                        VectorMult_mkl ( &invec[0], &outvec[0] );

                        return;

                        break;

                    case BLAS:

                        VectorMult_cblas ( &invec[0], &outvec[0] );

                        return;

                        break;

                    default:

                        ;

                }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::ForwardBackward ( const std::vector<DataType> &b, std::vector<DataType> &x )
        {

            //assert correct dimensions of x and b

            assert ( ncols_ == x.size ( ) );

            assert ( nrows_ == b.size ( ) );

            // copy b to x

            for ( int i = 0; i != ncols_; ++i )
            {

                x[i] = b[i];

            }

            if ( platform_id_ == CPU )
            {

                switch ( implementation_id_ )
                {

                    case NAIVE:

                        ForwardBackward_simple ( &x[0] );

                        return;

                        break;

                    case OPENMP:

                        ForwardBackward_simple ( &x[0] );

                        return;

                        break;

                    case MKL:

                        ForwardBackward_mkl ( &x[0] );

                        return;

                        break;

                    case BLAS:

                        ForwardBackward_cblas ( &x[0] );

                        return;

                        break;

                    default:

                        ;

                }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::ForwardBackward ( const lVector<DataType> &b, lVector<DataType> &x )
        {

            //assert correct dimensions of x and b

            assert ( ncols_ == x.get_size ( ) );

            assert ( nrows_ == b.get_size ( ) );

            // copy b to x

            std::vector<DataType> x_values ( ncols_, 0. );

            std::vector<DataType> b_values ( ncols_, 0. );

            b.GetBlockValues ( 0, ncols_, vec2ptr ( b_values ) );

            ForwardBackward ( b_values, x_values );

            x.SetBlockValues ( 0, ncols_, vec2ptr ( x_values ) );

            return;

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Solve ( const lVector<DataType> &b, lVector<DataType> &x )
        {

            Factorize ( );

            ForwardBackward ( b, x );

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Solve ( const std::vector<DataType> &b, std::vector<DataType> &x )
        {

            Factorize ( );

            ForwardBackward ( b, x );

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::Factorize ( )
        {

            //assert square matrix

            assert ( nrows_ == ncols_ );

            //allocate and prepare lu_ and permutation_

            lu_ = val_;

            permutation_.resize ( ncols_, 0 );

            if ( platform_id_ == CPU )
            {

                switch ( implementation_id_ )
                {

                    case NAIVE:

                        Factorize_simple ( );

                        return;

                        break;

                    case OPENMP:

                        Factorize_openmp ( );

                        return;

                        break;

                    case MKL:

                        Factorize_mkl ( );

                        return;

                        break;

                    case BLAS:

                        Factorize_cblas ( );

                        return;

                        break;

                    default:

                        ;

                }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::QRFactorize ( SeqDenseMatrix<DataType> &Q, SeqDenseMatrix<DataType> &R )
        {

            // assert suitable dimensions of this

            assert ( this->nrows_ >= this->ncols_ );

            // assert correct dimensions of Q

            assert ( Q.nrows ( ) == this->nrows_ );

            assert ( Q.ncols ( ) == this->ncols_ );

            // assert correct dimensions of R

            assert ( R.nrows ( ) == this->ncols_ );

            assert ( R.ncols ( ) == this->ncols_ );

            if ( platform_id_ == CPU )
            {

                switch ( implementation_id_ )
                {

                    case NAIVE:

                        QRFactorize_simple ( Q, R );

                        return;

                        break;

                    case OPENMP:

                        QRFactorize_openmp ( Q, R );

                        return;

                        break;

                    case MKL:

                        QRFactorize_mkl ( Q, R );

                        return;

                        break;

                    case BLAS:

                        QRFactorize_cblas ( Q, R );

                        return;

                        break;

                    default:

                        ;

                }

            }

        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::HessenbergReduction ( SeqDenseMatrix<DataType> &H, SeqDenseMatrix<DataType> &V )
        {
            if ( platform_id_ == CPU )
            {

                switch ( implementation_id_ )
                {

                    case NAIVE:

                        this->HessenbergReduction_simple ( H, V );

                        return;

                        break;

                    case OPENMP:

                        this->HessenbergReduction_openmp ( H, V );

                        return;

                        break;

                    case MKL:

                        this->HessenbergReduction_mkl ( H, V );

                        return;

                        break;

                    case BLAS:

                        this->HessenbergReduction_cblas ( H, V );

                        return;

                        break;

                    default:

                        ;

                }

            }
        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::HessenbergReduction_simple ( SeqDenseMatrix<DataType> &H, SeqDenseMatrix<DataType> &V )
        {
            // Initialize H to this
            for ( int i = 0; i != nrows_; ++i )
            {
                const int offset = i * ncols_;
                for ( int j = 0; j != ncols_; ++j )
                {
                    H ( i, j ) = val_[offset + j];
                }
            }

            /* transform H to Hessenberg form */
            SeqDenseMatrix<DataType> V_temp;
            V_temp.Resize ( nrows_, ncols_ );
            V.Zeros ( );

            for ( int i = 0; i < ncols_; ++i )
            {
                /* construct reflection vector */
                const int size_x = nrows_ - ( i + 1 );
                DataType *x = new DataType[size_x];

                // x = H(i:nrows_, i)
                for ( int j = 0; j < size_x; ++j )
                {
                    x[j] = H ( ( i + 1 ) + j, i );
                }

                //compute norm of x
                DataType norm_x = 0.;
                for ( int j = 0; j < size_x; ++j )
                {
                    norm_x += x[j] * x[j];
                }
                norm_x = std::sqrt ( norm_x );

                DataType tmp = 0.;
                if ( x[0] > 0 )
                {
                    tmp = norm_x;
                }
                else
                {
                    tmp = -norm_x;
                }

                x[0] += tmp;

                // recompute norm of x
                norm_x = 0.;
                for ( int j = 0; j < size_x; ++j )
                {
                    norm_x += x[j] * x[j];
                }
                norm_x = std::sqrt ( norm_x );

                if ( norm_x > 0.0 )
                {
                    // normalize x
                    for ( int j = 0; j < size_x; ++j )
                    {
                        x[j] /= norm_x;
                    }

                    /* apply reflection from the left */
                    // precompute x'*H(k+1:end,k:end)
                    int size_h = ncols_ - i;
                    DataType *h = new DataType[size_h];
                    memsethost ( h, 0, size_h, sizeof (DataType ) );

                    for ( int k = 0; k < size_x; ++k )
                    {
                        const DataType x_temp = x[k];
                        const int i_index = i + 1 + k;
                        for ( int j = 0; j < size_h; ++j )
                        {
                            h[j] += x_temp * H ( i_index, i + j );
                        }
                    }

                    // H((k+1:end,k:end) = H(k+1:end,k:end) - 2*x*(x'*H(k+1:end,k:end))
                    for ( int j = i + 1; j < nrows_; ++j )
                    {
                        DataType x_temp = x[j - ( i + 1 )];
                        for ( int k = i; k < ncols_; ++k )
                        {
                            H ( j, k ) -= 2 * x_temp * h[k - i];
                        }
                    }

                    /* apply refelction from the right */
                    // precompute h = H(1:n, k+1:n)*x
                    delete [] h;
                    h = new DataType[nrows_];
                    memsethost ( h, 0, nrows_, sizeof (DataType ) );
                    for ( int k = 0; k < nrows_; ++k )
                    {
                        for ( int j = i + 1; j < ncols_; ++j )
                        {
                            h[k] += H ( k, j ) * x[j - ( i + 1 )];
                        }
                    }

                    // H(1:end,k+1:end) = H(1:end,k+1:end) - 2*(H(1:end,k+1:end)*x)*x'
                    for ( int j = 0; j < nrows_; ++j )
                    {
                        const DataType h_temp = h[j];
                        for ( int k = i + 1; k < ncols_; ++k )
                        {
                            H ( j, k ) -= 2 * h_temp * x[k - ( i + 1 )];
                        }
                    }

                    // save reflection vector
                    for ( int j = 0; j < size_x; ++j )
                    {
                        V_temp ( i, i + j ) = x[j];
                    }

                    delete [] h;
                }

                delete [] x;
            }

            // Initialize V to identity
            for ( int i = 0; i != nrows_; ++i )
            {
                V ( i, i ) = 1.;
            }

            // Compute EV_matrix_new
            for ( int i = ncols_ - 1; i >= 0; i-- )
            {
                /* apply reflection */
                const int size_x = nrows_ - ( i + 1 );
                DataType *x = new DataType[size_x];
                for ( int j = 0; j < size_x; ++j )
                {
                    x[j] = V_temp ( i, i + j );
                }

                // precompute x'*EV_matrix_new(k+1:end,k+1:end)
                const int size_h = ncols_ - ( i + 1 );
                DataType *h = new DataType[size_h];
                memsethost ( h, 0, size_h, sizeof (DataType ) );

                for ( int k = 0; k < size_x; ++k )
                {
                    const DataType x_temp = x[k];
                    const int i_index = i + 1 + k;
                    for ( int j = 0; j < size_h; ++j )
                    {
                        h[j] += x_temp * V ( i_index, ( i + 1 ) + j );
                    }
                }

                // EV_matrix_new((k+1:end,k+1:end) = EV_matrix_new(k+1:end,k+1:end) - 2*x*(x'*EV_matrix_new(k+1:end,k+1:end))
                for ( int j = i + 1; j < nrows_; ++j )
                {
                    const DataType x_temp = x[j - ( i + 1 )];
                    for ( int k = i + 1; k < ncols_; ++k )
                    {
                        V ( j, k ) -= 2 * x_temp * h[k - ( i + 1 )];
                    }
                }

                delete [] x;
                delete [] h;
            }
        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::HessenbergReduction_openmp ( SeqDenseMatrix<DataType> &H, SeqDenseMatrix<DataType> &V )
        {
#ifdef WITH_OPENMP
            // Initialize H to this
#    pragma omp parallel for
            for ( int i = 0; i != nrows_; ++i )
            {
                const int offset = i * ncols_;
                for ( int j = 0; j != ncols_; ++j )
                {
                    H ( i, j ) = val_[offset + j];
                }
            }

            /* transform H to Hessenberg form */
            SeqDenseMatrix<DataType> V_temp;
            V_temp.Resize ( nrows_, ncols_ );
            V.Zeros ( );

            for ( int i = 0; i < ncols_; ++i )
            {
                /* construct reflection vector */
                const int size_x = nrows_ - ( i + 1 );
                DataType *x = new DataType[size_x];

                // x = H(i:nrows_, i)
#    pragma omp parallel for
                for ( int j = 0; j < size_x; ++j )
                {
                    x[j] = H ( ( i + 1 ) + j, i );
                }

                //compute norm of x
                DataType norm_x = 0.;
#    pragma omp parallel for reduction(+:norm_x) shared(x)
                for ( int j = 0; j < size_x; ++j )
                {
                    norm_x += x[j] * x[j];
                }
                norm_x = std::sqrt ( norm_x );

                DataType tmp = 0.;
                if ( x[0] > 0 )
                {
                    tmp = norm_x;
                }
                else
                {
                    tmp = -norm_x;
                }

                x[0] += tmp;

                // recompute norm of x
                norm_x = 0.;
#    pragma omp parallel for reduction(+:norm_x) shared(x)
                for ( int j = 0; j < size_x; ++j )
                {
                    norm_x += x[j] * x[j];
                }
                norm_x = std::sqrt ( norm_x );

                if ( norm_x > 0.0 )
                {
                    // normalize x
#    pragma omp parallel for
                    for ( int j = 0; j < size_x; ++j )
                    {
                        x[j] /= norm_x;
                    }

                    /* apply reflection from the left */
                    // precompute x'*H(k+1:end,k:end)
                    int size_h = ncols_ - i;
                    DataType *h = new DataType[size_h];
                    memsethost ( h, 0, size_h, sizeof (DataType ) );

                    for ( int k = 0; k < size_x; ++k )
                    {
                        const DataType x_temp = x[k];
                        const int i_index = i + 1 + k;
#    pragma omp parallel for
                        for ( int j = 0; j < size_h; ++j )
                        {
                            h[j] += x_temp * H ( i_index, i + j );
                        }
                    }

                    // H((k+1:end,k:end) = H(k+1:end,k:end) - 2*x*(x'*H(k+1:end,k:end))
#    pragma omp parallel for
                    for ( int j = i + 1; j < nrows_; ++j )
                    {
                        DataType x_temp = x[j - ( i + 1 )];
                        for ( int k = i; k < ncols_; ++k )
                        {
                            H ( j, k ) -= 2 * x_temp * h[k - i];
                        }
                    }

                    /* apply refelction from the right */
                    // precompute h = H(1:n, k+1:n)*x
                    delete [] h;
                    h = new DataType[nrows_];
                    memsethost ( h, 0, nrows_, sizeof (DataType ) );
#    pragma omp parallel for
                    for ( int k = 0; k < nrows_; ++k )
                    {
                        for ( int j = i + 1; j < ncols_; ++j )
                        {
                            h[k] += H ( k, j ) * x[j - ( i + 1 )];
                        }
                    }

                    // H(1:end,k+1:end) = H(1:end,k+1:end) - 2*(H(1:end,k+1:end)*x)*x'
#    pragma omp parallel for
                    for ( int j = 0; j < nrows_; ++j )
                    {
                        const DataType h_temp = h[j];
                        for ( int k = i + 1; k < ncols_; ++k )
                        {
                            H ( j, k ) -= 2 * h_temp * x[k - ( i + 1 )];
                        }
                    }

                    // save reflection vector
#    pragma omp parallel for
                    for ( int j = 0; j < size_x; ++j )
                    {
                        V_temp ( i, i + j ) = x[j];
                    }

                    delete [] h;
                }

                delete [] x;
            }

            // Initialize V to identity
#    pragma omp parallel for
            for ( int i = 0; i != nrows_; ++i )
            {
                V ( i, i ) = 1.;
            }

            // Compute EV_matrix_new
            for ( int i = ncols_ - 1; i >= 0; i-- )
            {
                /* apply reflection */
                const int size_x = nrows_ - ( i + 1 );
                DataType *x = new DataType[size_x];
#    pragma omp parallel for
                for ( int j = 0; j < size_x; ++j )
                {
                    x[j] = V_temp ( i, i + j );
                }

                // precompute x'*EV_matrix_new(k+1:end,k+1:end)
                const int size_h = ncols_ - ( i + 1 );
                DataType *h = new DataType[size_h];
                memsethost ( h, 0, size_h, sizeof (DataType ) );

                for ( int k = 0; k < size_x; ++k )
                {
                    const DataType x_temp = x[k];
                    const int i_index = i + 1 + k;
#    pragma omp parallel for
                    for ( int j = 0; j < size_h; ++j )
                    {
                        h[j] += x_temp * V ( i_index, ( i + 1 ) + j );
                    }
                }

                // EV_matrix_new((k+1:end,k+1:end) = EV_matrix_new(k+1:end,k+1:end) - 2*x*(x'*EV_matrix_new(k+1:end,k+1:end))
#    pragma omp parallel for
                for ( int j = i + 1; j < nrows_; ++j )
                {
                    const DataType x_temp = x[j - ( i + 1 )];
                    for ( int k = i + 1; k < ncols_; ++k )
                    {
                        V ( j, k ) -= 2 * x_temp * h[k - ( i + 1 )];
                    }
                }

                delete [] x;
                delete [] h;
            }
#else
            ERROR;
#endif
        }

        template<class DataType>

        void SeqDenseMatrix<DataType>::ForwardBackward_simple ( DataType *x )
        {

            // result vector for forward solving

            DataType *z = new DataType[nrows_];

            // forward solving

            z[0] = x[permutation_[0]];

            for ( int i = 1; i < nrows_; ++i )
            {

                DataType sum = 0;

                const int index = i*ncols_;

                for ( int j = 0; j < i; ++j )
                {

                    sum += lu_[index + j] * z[j];

                }

                z[i] = x[permutation_[i]] - sum;

            }

            //backward solving

            x[nrows_ - 1] = z[nrows_ - 1] / lu_[( nrows_ - 1 ) * ncols_ + ncols_ - 1];

            for ( int i = nrows_; i--; )
            {

                DataType sum = 0;

                const int index = i*ncols_;

                for ( int j = i + 1; j < ncols_; ++j )
                {

                    sum += lu_[index + j] * x[j];

                }

                x[i] = ( z[i] - sum ) / lu_[index + i];

            }

            delete [] z;

        }

        template<class DataType>

        DataType SeqDenseMatrix<DataType>::Frob ( void )
        {

            DataType frob_norm = 0.0;

            for ( int i = 0; i < val_.size ( ); ++i )
            {

                frob_norm += val_[i] * val_[i];

            }

            return std::sqrt ( frob_norm );

        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::gemm_simple ( int i, int blockend )
        {
            for ( int k = blockend; k < nrows_; ++k )
            {
                const int index = k *ncols_;
                for ( int m = i; m < blockend; ++m )
                {
                    const DataType temp = lu_[index + m];
                    const int m_index = m*ncols_;
                    if ( temp != 0. )
                    {
                        for ( int l = blockend; l < ncols_; ++l )
                        {
                            lu_[index + l] -= ( temp * lu_[m_index + l] );
                        }
                    }
                }
            }
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::gemm_openmp ( const int i, const int blockend )
        {
#ifdef WITH_OPENMP
#    pragma omp parallel for
            for ( int k = blockend; k < nrows_; ++k )
            {
                const int index = k *ncols_;
                for ( int m = i; m < blockend; ++m )
                {
                    const DataType temp = lu_[index + m];
                    const int m_index = m*ncols_;
                    if ( temp != 0. )
                    {
                        for ( int l = blockend; l < ncols_; ++l )
                        {
                            lu_[index + l] -= ( temp * lu_[m_index + l] );
                        }
                    }
                }
            }
#else
            ERROR;
#endif
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::MatrixMult_simple ( const SeqDenseMatrix<DataType> *inMat,
                                                           SeqDenseMatrix<DataType> *outMat )
        {
            outMat->Zeros ( );
            const int out_cols = outMat->ncols ( );
            // loop over all rows of outMat
            for ( int i = 0; i != nrows_; ++i )
            {
                // index of first element of current row in outMat
                const int out_row_index = i * out_cols;

                // index of first element of current row in this
                const int this_row_index = i*ncols_;

                // loop over one row of this and one column of inMat
                for ( int k = 0; k != ncols_; ++k )
                {
                    // current element of this
                    const DataType temp = val_[this_row_index + k];

                    // index of first element of current row in inMat
                    // observe, that inMat.ncols() == outMat.ncols()
                    const int in_row_index = k * out_cols;

                    // loop over all columns of outMat
                    for ( int j = 0; j != out_cols; ++j )
                    {
                        outMat->val_[out_row_index + j] += temp * inMat->val_[in_row_index + j];
                    }
                }
            }
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::MatrixMult_openmp ( const SeqDenseMatrix<DataType> *inMat,
                                                           SeqDenseMatrix<DataType> *outMat )
        {
#ifdef WITH_OPENMP
            outMat->Zeros ( );
            const int out_cols = outMat->ncols ( );
            // loop over all rows of outMat
#    pragma omp parallel for
            for ( int i = 0; i < nrows_; i++ )
            {
                // index of first element of current row in outMat
                const int out_row_index = i * out_cols;

                // index of first element of current row in this
                const int this_row_index = i*ncols_;

                // loop over one row of this and one column of inMat
                for ( int k = 0; k < ncols_; k++ )
                {
                    // current element of this
                    const DataType temp = val_[this_row_index + k];

                    // index of first element of current row in inMat
                    // observe, that inMat.ncols() == outMat.ncols()
                    const int in_row_index = k * out_cols;

                    // loop over all columns of outMat
                    for ( int j = 0; j < out_cols; j++ )
                    {
                        outMat->val_[out_row_index + j] += temp * inMat->val_[in_row_index + j];
                    }
                }
            }
#else
            ERROR;
#endif
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::compute_eigenvalues_and_vectors_simple ( )
        {
            // Initialize data structures
            SeqDenseMatrix<DataType> Z_new, EV_matrix_new;
            Z_new.set_implementation_and_platform ( implementation_id_, platform_id_ );
            EV_matrix_new.set_implementation_and_platform ( implementation_id_, platform_id_ );

            Z_new.Resize ( nrows_, ncols_ );

            EV_matrix_new.Resize ( nrows_, ncols_ );

            this->HessenbergReduction ( Z_new, EV_matrix_new );

            /* Compute eigenvalues and -vectors of Z_new which is similar to this */

            // local data structures
            SeqDenseMatrix<DataType> Q_local, Z_new_local, R_local, Z_local, EV_local, EV_new_local;
            Z_local.set_implementation_and_platform ( implementation_id_, platform_id_ );
            R_local.set_implementation_and_platform ( implementation_id_, platform_id_ );
            EV_local.set_implementation_and_platform ( implementation_id_, platform_id_ );

            // Iteration
            for ( int i = nrows_ - 1; i >= 1; i-- )
            {

                const int len = i + 1;

                // prepare local data structures
                Q_local.Resize ( len, len );
                Z_new_local.Resize ( len, len );
                R_local.Resize ( len, len );
                Z_local.Resize ( len, len );
                EV_local.Resize ( len, len );
                EV_new_local.Resize ( len, len );

                Z_new.extr ( Z_new_local, 0, 0, i, i );
                EV_matrix_new.extr ( EV_new_local, 0, 0, i, i );

                DataType error = 0.;

                // repeat until convergence
                do
                {
                    const DataType sigma = Z_new_local ( i, i );

                    //copy Z_new_local to Z_local
                    for ( int j = 0; j < len; ++j )
                    {
                        for ( int k = 0; k < len; ++k )
                        {
                            Z_local ( j, k ) = Z_new_local ( j, k );
                        }
                    }

                    // modify Z_local
                    for ( int j = 0; j < len; ++j )
                    {
                        Z_local ( j, j ) -= sigma;
                    }

                    // Compute Z_local = Q_local * R_local
                    Z_local.QRFactorize ( Q_local, R_local );

                    // Compute Z_new_local = R * Q
                    R_local.MatrixMult ( Q_local, Z_new_local );

                    // add sigma*I
                    for ( int j = 0; j < len; ++j )
                    {
                        Z_new_local ( j, j ) += sigma;
                    }

                    // Copy EV_new_local to EV_local
                    for ( int j = 0; j < len; ++j )
                    {
                        for ( int k = 0; k < len; ++k )
                        {
                            EV_local ( j, k ) = EV_new_local ( j, k );
                        }
                    }

                    //Compute EV_new_local = EV_local * Q_local
                    EV_local.MatrixMult ( Q_local, EV_new_local );

                    // convergence criterion
                    const DataType diff = Z_new_local ( i, i - 1 );
                    error = std::abs ( diff );
                }
                while ( error > static_cast < DataType > ( 0 ) );

                // Insert small matrices into global ones
                Z_new.ins ( Z_new_local, 0, 0 );
                EV_matrix_new.ins ( EV_new_local, 0, 0 );
            }

            // get eigenvalues
            DataType *ev_temp = new DataType[nrows_];
            for ( int i = 0; i < nrows_; ++i )
            {
                ev_temp[i] = Z_new ( i, i );
            }

            // sort eigenvalues and vectors
            int *perm = new int[nrows_];
            for ( int i = 0; i < nrows_; ++i )
            {
                perm[i] = i;
            }

            //sort eigenvalues in ascending order
            for ( int i = 0; i < nrows_; ++i )
            {
                // find smallest eigenvalue in later indizes
                DataType current_min = ev_temp[perm[i]];
                int current_index = i;
                for ( int j = i; j < nrows_; ++j )
                {
                    const DataType temp = ev_temp[perm[j]];
                    if ( temp < current_min )
                    {
                        current_min = temp;
                        current_index = j;
                    }
                }
                int temp = perm[i];
                perm[i] = perm[current_index];
                perm[current_index] = temp;
            }

            // sort eigenvalues and vectors for output
            for ( int i = 0; i < nrows_; ++i )
            {
                const int perm_index = perm[i];
                // eigenvalues
                eigenvalues_[i] = ev_temp[perm_index];

                // eigenvectors
                for ( int j = 0; j < nrows_; ++j )
                {
                    eigenvectors_[j * ncols_ + i] = EV_matrix_new ( j, perm_index );
                }
            }

            delete [] ev_temp;
            delete [] perm;
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::compute_eigenvalues_and_vectors_openmp ( )
        {
#ifdef WITH_OPENMP
            // Initialize data structures
            SeqDenseMatrix<DataType> Z_new, EV_matrix_new;
            Z_new.set_implementation_and_platform ( implementation_id_, platform_id_ );
            EV_matrix_new.set_implementation_and_platform ( implementation_id_, platform_id_ );

            Z_new.Resize ( nrows_, ncols_ );

            EV_matrix_new.Resize ( nrows_, ncols_ );

            this->HessenbergReduction ( Z_new, EV_matrix_new );

            /* Compute eigenvalues and -vectors of Z_new which is similar to this */

            // local data structures
            SeqDenseMatrix<DataType> Q_local, Z_new_local, R_local, Z_local, EV_local, EV_new_local;
            Z_local.set_implementation_and_platform ( implementation_id_, platform_id_ );
            R_local.set_implementation_and_platform ( implementation_id_, platform_id_ );
            EV_local.set_implementation_and_platform ( implementation_id_, platform_id_ );

            // Iteration
            for ( int i = nrows_ - 1; i >= 1; i-- )
            {

                const int len = i + 1;

                // prepare local data structures
                Q_local.Resize ( len, len );
                Z_new_local.Resize ( len, len );
                R_local.Resize ( len, len );
                Z_local.Resize ( len, len );
                EV_local.Resize ( len, len );
                EV_new_local.Resize ( len, len );

                Z_new.extr ( Z_new_local, 0, 0, i, i );
                EV_matrix_new.extr ( EV_new_local, 0, 0, i, i );

                DataType error = 0.;

                // repeat until convergence
                do
                {
                    const DataType sigma = Z_new_local ( i, i );

                    //copy Z_new_local to Z_local
#    pragma omp parallel for
                    for ( int j = 0; j < len; ++j )
                    {
                        for ( int k = 0; k < len; ++k )
                        {
                            Z_local ( j, k ) = Z_new_local ( j, k );
                        }
                    }

                    // modify Z_local
#    pragma omp parallel for
                    for ( int j = 0; j < len; ++j )
                    {
                        Z_local ( j, j ) -= sigma;
                    }

                    // Compute Z_local = Q_local * R_local
                    Z_local.QRFactorize ( Q_local, R_local );

                    // Compute Z_new_local = R * Q
                    R_local.MatrixMult ( Q_local, Z_new_local );

                    // add sigma*I
#    pragma omp parallel for
                    for ( int j = 0; j < len; ++j )
                    {
                        Z_new_local ( j, j ) += sigma;
                    }

                    // Copy EV_new_local to EV_local
#    pragma omp parallel for
                    for ( int j = 0; j < len; ++j )
                    {
                        for ( int k = 0; k < len; ++k )
                        {
                            EV_local ( j, k ) = EV_new_local ( j, k );
                        }
                    }

                    //Compute EV_new_local = EV_local * Q_local
                    EV_local.MatrixMult ( Q_local, EV_new_local );

                    // convergence criterion
                    const DataType diff = Z_new_local ( i, i - 1 );
                    error = std::abs ( diff );
                }
                while ( error > static_cast < DataType > ( 0 ) );

                // Insert small matrices into global ones
                Z_new.ins ( Z_new_local, 0, 0 );
                EV_matrix_new.ins ( EV_new_local, 0, 0 );
            }

            // get eigenvalues
            DataType *ev_temp = new DataType[nrows_];
#    pragma omp parallel for
            for ( int i = 0; i < nrows_; ++i )
            {
                ev_temp[i] = Z_new ( i, i );
            }

            // sort eigenvalues and vectors
            int *perm = new int[nrows_];
#    pragma omp parallel for
            for ( int i = 0; i < nrows_; ++i )
            {
                perm[i] = i;
            }

            //sort eigenvalues in ascending order
            for ( int i = 0; i < nrows_; ++i )
            {
                // find smallest eigenvalue in later indizes
                DataType current_min = ev_temp[perm[i]];
                int current_index = i;
                for ( int j = i; j < nrows_; ++j )
                {
                    const DataType temp = ev_temp[perm[j]];
                    if ( temp < current_min )
                    {
                        current_min = temp;
                        current_index = j;
                    }
                }
                int temp = perm[i];
                perm[i] = perm[current_index];
                perm[current_index] = temp;
            }

            // sort eigenvalues and vectors for output
            for ( int i = 0; i < nrows_; ++i )
            {
                const int perm_index = perm[i];
                // eigenvalues
                eigenvalues_[i] = ev_temp[perm_index];

                // eigenvectors
                for ( int j = 0; j < nrows_; ++j )
                {
                    eigenvectors_[j * ncols_ + i] = EV_matrix_new ( j, perm_index );
                }
            }

            delete [] ev_temp;
            delete [] perm;
#else
            ERROR;
#endif
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::VectorMult_simple ( const DataType *invec, DataType *outvec )
        {
            for ( int i = 0; i != nrows_; ++i )
            {
                outvec[i] = 0;
                const int index = i * ncols_;
                for ( int j = 0; j != ncols_; ++j )
                {
                    outvec[i] += val_[index + j] * invec[j];
                }
            }
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::VectorMult_openmp ( const DataType *invec, DataType *outvec )
        {
#ifdef WITH_OPENMP

#    pragma omp parallel for shared(invec, outvec)
            for ( int i = 0; i < nrows_; ++i )
            {
                const int index = i * ncols_;
                register DataType res = 0.;
                for ( int j = 0; j < ncols_; ++j )
                {
                    res += val_[index + j] * invec[j];
                }
                outvec[i] = res;
            }
#else
            ERROR;
#endif
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::Factorize_simple ( )
        {
            //initialize permutation_
            for ( int i = 0; i != nrows_; ++i )
            {
                permutation_[i] = i;
            }

            int j = 0;
            // Beginning of gaussian elimination
            for ( int i = 0; i < nrows_; i += block_ )
            {
                // number of row or column, respectively, of the beginning of the block of the delayed updates
                j = i + block_;

                const int blockend = std::min ( nrows_, j );

                //gaussian elimination on the diagonal block
                for ( int k = i; k < blockend; ++k )
                {
                    //find maximum entry in current column
                    int max_index = k;
                    DataType max_value = std::abs ( lu_[k * ncols_ + k] );
                    for ( int l = k + 1; l < blockend; ++l )
                    {
                        DataType current_abs = std::abs ( lu_[l * ncols_ + k] );
                        if ( current_abs > max_value )
                        {
                            max_value = current_abs;
                            max_index = l;
                        }
                    }

                    //check for singular matrix
                    if ( max_value == 0.0 )
                    {
                        LOG_ERROR ( "The diagonal element of the triangular factor of A, U(" << k << "," << k << ") is zero, so that A is singular; the solution could not be computed" );
                        exit ( -1 );
                    }

                    // interchange rows
                    const int old_row_offset = k * ncols_;
                    const int new_row_offset = max_index * ncols_;
                    for ( int i = 0; i != ncols_; ++i )
                    {
                        const DataType temp = lu_[old_row_offset + i];
                        lu_[old_row_offset + i] = lu_[new_row_offset + i];
                        lu_[new_row_offset + i] = temp;
                    }

                    const int temp = permutation_[k];
                    permutation_[k] = permutation_[max_index];
                    permutation_[max_index] = temp;

                    //update diagonal block
                    const int k_index = k * ncols_;
                    const DataType diag = lu_[k_index + k];
                    for ( int l = k + 1; l < blockend; ++l )
                    {
                        const int index = l * ncols_;
                        lu_[index + k] /= diag;
                        const DataType temp = lu_[index + k];
                        for ( int m = k + 1; m < blockend; ++m )
                        {
                            lu_[index + m] -= ( temp * lu_[k_index + m] );
                        }
                    }
                }

                // update L block
                for ( int l = blockend; l < nrows_; ++l )
                {
                    const int index = l * ncols_;
                    for ( int k = i; k < blockend; ++k )
                    {
                        const int k_index = k * ncols_;
                        lu_[index + k] /= lu_[k_index + k];
                        const DataType temp = lu_[index + k];
                        for ( int m = k + 1; m < blockend; ++m )
                        {
                            lu_[index + m] -= ( temp * lu_[k_index + m] );
                        }
                    }
                }

                // update U-block
                for ( int l = i + 1; l < blockend; ++l )
                {
                    const int index = l * ncols_;
                    for ( int m = i; m < l; ++m )
                    {
                        const DataType fac = lu_[index + m];
                        const int index_m = m*ncols_;
                        if ( fac != 0. )
                        {
                            for ( int k = blockend; k < ncols_; ++k )
                            {
                                lu_[index + k] -= ( fac * lu_[index_m + k] );
                            }
                        }
                    }
                }

                gemm_simple ( i, blockend );
            }
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::Factorize_openmp ( )
        {
#ifdef WITH_OPENMP

            //initialize permutation_
#    pragma omp parallel for
            for ( int i = 0; i < nrows_; ++i )
            {
                permutation_[i] = i;
            }

            int j = 0;
            // Beginning of gaussian elimination
            for ( int i = 0; i < nrows_; i += block_ )
            {
                // number of row or column, respectively, of the beginning of the block of the delayed updates
                j = i + block_;

                int blockend = std::min ( nrows_, j );

                //gaussian elimination on the diagonal block
                for ( int k = i; k < blockend; ++k )
                {
                    //find maximum entry in current column
                    int max_index = k;
                    DataType max_value = std::abs ( lu_[k * ncols_ + k] );
#    pragma omp parallel
                    {
                        int loc_ind = k;
                        DataType loc_max = std::abs ( lu_[k * ncols_ + k] );
#    pragma omp for
                        for ( int l = k + 1; l < blockend; ++l )
                        {
                            DataType current_abs = std::abs ( lu_[l * ncols_ + k] );
                            if ( current_abs > max_value )
                            {
                                loc_max = current_abs;
                                loc_ind = l;
                            }
                        }

#    pragma omp critical
                        {
                            if ( loc_max > max_value )
                            {
                                max_value = loc_max;
                                max_index = loc_ind;
                            }
                        }

                    }

                    //check for singular matrix
                    if ( max_value == 0.0 )
                    {
                        LOG_ERROR ( "The diagonal element of the triangular factor of A, U(" << k << "," << k << ") is zero, so that A is singular; the solution could not be computed" );
                        exit ( -1 );
                    }

                    // interchange rows
                    const int old_row_offset = k * ncols_;
                    const int new_row_offset = max_index * ncols_;
#    pragma omp parallel for
                    for ( int i = 0; i < ncols_; ++i )
                    {
                        const DataType temp = lu_[old_row_offset + i];
                        lu_[old_row_offset + i] = lu_[new_row_offset + i];
                        lu_[new_row_offset + i] = temp;
                    }

                    const int temp = permutation_[k];
                    permutation_[k] = permutation_[max_index];
                    permutation_[max_index] = temp;

                    //update diagonal block
                    const int k_index = k * ncols_;
                    const DataType diag = lu_[k_index + k];
#    pragma omp parallel for
                    for ( int l = k + 1; l < blockend; ++l )
                    {
                        const int index = l * ncols_;
                        lu_[index + k] /= diag;
                        const DataType temp = lu_[index + k];
                        for ( int m = k + 1; m < blockend; ++m )
                        {
                            lu_[index + m] -= ( temp * lu_[k_index + m] );
                        }
                    }
                }

                // update L block
#    pragma omp parallel for
                for ( int l = blockend; l < nrows_; ++l )
                {
                    const int index = l * ncols_;
                    for ( int k = i; k < blockend; ++k )
                    {
                        const int k_index = k * ncols_;
                        lu_[index + k] /= lu_[k_index + k];
                        const DataType temp = lu_[index + k];
                        for ( int m = k + 1; m < blockend; ++m )
                        {
                            lu_[index + m] -= ( temp * lu_[k_index + m] );
                        }
                    }
                }

                // update U-block
                for ( int l = i + 1; l < blockend; ++l )
                {
                    const int index = l * ncols_;
                    for ( int m = i; m < l; ++m )
                    {
                        const DataType fac = lu_[index + m];
                        const int index_m = m*ncols_;
                        if ( fac != 0. )
                        {
#    pragma omp parallel for
                            for ( int k = blockend; k < ncols_; ++k )
                            {
                                lu_[index + k] -= ( fac * lu_[index_m + k] );
                            }
                        }
                    }
                }

                gemm_openmp ( i, blockend );
            }
#else
            ERROR;
#endif
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::QRFactorize_simple ( SeqDenseMatrix<DataType> &Q, SeqDenseMatrix<DataType> &R )
        {

            // initialize data structures
            SeqDenseMatrix<DataType> temp;
            temp.set_implementation_and_platform ( implementation_id_, platform_id_ );
            temp.Resize ( nrows_, ncols_ );

            for ( int i = 0; i < nrows_; ++i )
            {
                const int i_index = i*ncols_;
                for ( int j = 0; j < ncols_; ++j )
                {
                    temp ( i, j ) = val_[i_index + j];
                }
            }

            SeqDenseMatrix<DataType> V;
            V.set_implementation_and_platform ( implementation_id_, platform_id_ );
            V.Resize ( ncols_, nrows_ );
            V.Zeros ( );

            Q.Zeros ( );
            for ( int i = 0; i < ncols_; ++i )
            {
                Q ( i, i ) = 1.;
            }

            R.Zeros ( );

            for ( int i = 0; i < ncols_; ++i )
            {
                /* construct reflection vector */
                const int size_x = nrows_ - i;
                DataType *x = new DataType[size_x];

                // x = temp(i:nrows_, i)
                for ( int j = 0; j < size_x; ++j )
                {
                    x[j] = temp ( i + j, i );
                }

                //compute norm of x
                DataType norm_x = 0.;
                for ( int j = 0; j < size_x; ++j )
                {
                    norm_x += x[j] * x[j];
                }
                norm_x = std::sqrt ( norm_x );

                DataType tmp = 0.;
                if ( x[0] > 0 )
                {
                    tmp = norm_x;
                }
                else
                {
                    tmp = -norm_x;
                }

                x[0] += tmp;

                // recompute norm of x
                norm_x = 0.;
                for ( int j = 0; j < size_x; ++j )
                {
                    norm_x += x[j] * x[j];
                }
                norm_x = std::sqrt ( norm_x );

                if ( norm_x > 0.0 )
                {
                    // normalize x
                    for ( int j = 0; j < size_x; ++j )
                    {
                        x[j] /= norm_x;
                    }

                    /* apply reflection */
                    // precompute x'*temp(k:end,k:end)
                    const int size_h = ncols_ - i;
                    DataType *h = new DataType[size_h];
                    memsethost ( h, 0, size_h, sizeof (DataType ) );

                    for ( int k = 0; k < size_x; ++k )
                    {
                        const DataType x_temp = x[k];
                        const int i_index = i + k;
                        for ( int j = 0; j < size_h; ++j )
                        {
                            h[j] += x_temp * temp ( i_index, i + j );
                        }
                    }

                    // tmp(k:end,k:end) = temp(k:end,k:end) - 2*x*(x'*temp(k:end,k:end))
                    for ( int j = i; j < nrows_; ++j )
                    {
                        const DataType x_temp = x[j - i];
                        for ( int k = i; k < ncols_; ++k )
                        {
                            temp ( j, k ) -= 2 * x_temp * h[k - i];
                        }
                    }

                    // save reflection vector
                    for ( int j = 0; j < size_x; ++j )
                    {
                        V ( i, i + j ) = x[j];
                    }

                    delete [] h;
                }

                delete [] x;
            }

            // Compute Q
            for ( int i = ncols_ - 1; i >= 0; --i )
            {
                /* apply reflection */
                const int size_x = nrows_ - i;
                DataType *x = new DataType[size_x];
                for ( int j = i; j < nrows_; ++j )
                {
                    x[j - i] = V ( i, j );
                }
                // precompute x'*Q(k:end,k:end)
                const int size_h = ncols_ - i;
                DataType *h = new DataType[size_h];
                memsethost ( h, 0, size_h, sizeof (DataType ) );

                for ( int k = 0; k < size_x; ++k )
                {
                    const DataType x_temp = x[k];
                    const int i_index = i + k;
                    for ( int j = 0; j < size_h; ++j )
                    {
                        h[j] += x_temp * Q ( i_index, i + j );
                    }
                }

                // Q(k:end,k:end) = Q(k:end,k:end) - 2*x*(x'*Q(k:end,k:end))
                for ( int j = i; j < nrows_; ++j )
                {
                    const DataType x_temp = x[j - i];
                    for ( int k = i; k < ncols_; ++k )
                    {
                        Q ( j, k ) -= 2 * x_temp * h[k - i];
                    }
                }

                delete [] x;
                delete [] h;
            }

            // Compute R = upper triangular of temp
            for ( int i = 0; i < ncols_; ++i )
            {
                if ( temp ( i, i ) >= 0 )
                {
                    for ( int j = i; j < ncols_; ++j )
                    {
                        R ( i, j ) = temp ( i, j );
                    }
                }
                else
                {
                    for ( int j = i; j < ncols_; ++j )
                    {
                        R ( i, j ) = -temp ( i, j );
                    }
                    for ( int j = 0; j < nrows_; ++j )
                    {
                        Q ( j, i ) = -Q ( j, i );
                    }
                }
            }
        }

        template<class DataType>
        void SeqDenseMatrix<DataType>::QRFactorize_openmp ( SeqDenseMatrix<DataType> &Q, SeqDenseMatrix<DataType> &R )
        {
#ifdef WITH_OPENMP

            // initialize data structures
            SeqDenseMatrix<DataType> temp;
            temp.set_implementation_and_platform ( implementation_id_, platform_id_ );
            temp.Resize ( nrows_, ncols_ );

#    pragma omp parallel for
            for ( int i = 0; i < nrows_; ++i )
            {
                const int i_index = i*ncols_;
                for ( int j = 0; j < ncols_; ++j )
                {
                    temp ( i, j ) = val_[i_index + j];
                }
            }

            SeqDenseMatrix<DataType> V;
            V.set_implementation_and_platform ( implementation_id_, platform_id_ );
            V.Resize ( ncols_, nrows_ );
            V.Zeros ( );

            Q.Zeros ( );

#    pragma omp parallel for
            for ( int i = 0; i < ncols_; ++i )
            {
                Q ( i, i ) = 1.;
            }

            R.Zeros ( );

            for ( int i = 0; i < ncols_; ++i )
            {
                /* construct reflection vector */
                const int size_x = nrows_ - i;
                DataType *x = new DataType[size_x];

                // x = temp(i:nrows_, i)
#    pragma omp parallel for shared(x)
                for ( int j = 0; j < size_x; ++j )
                {
                    x[j] = temp ( i + j, i );
                }

                //compute norm of x
                DataType norm_x = 0.;
#    pragma omp parallel for reduction(+:norm_x) shared(x)
                for ( int j = 0; j < size_x; ++j )
                {
                    norm_x += x[j] * x[j];
                }
                norm_x = std::sqrt ( norm_x );

                DataType tmp = 0.;
                if ( x[0] > 0 )
                {
                    tmp = norm_x;
                }
                else
                {
                    tmp = -norm_x;
                }

                x[0] += tmp;

                // recompute norm of x
                norm_x = 0.;
#    pragma omp parallel for reduction(+:norm_x) shared(x)
                for ( int j = 0; j < size_x; ++j )
                {
                    norm_x += x[j] * x[j];
                }
                norm_x = std::sqrt ( norm_x );

                if ( norm_x > 0.0 )
                {
                    // normalize x
#    pragma omp parallel for shared(x)
                    for ( int j = 0; j < size_x; ++j )
                    {
                        x[j] /= norm_x;
                    }

                    /* apply reflection */
                    // precompute x'*temp(k:end,k:end)
                    const int size_h = ncols_ - i;
                    DataType *h = new DataType[size_h];
                    memsethost ( h, 0, size_h, sizeof (DataType ) );

#    pragma omp parallel shared(x,h)
                    {
                        DataType *loc_h = new DataType[size_h];
                        for ( int j = 0; j < size_h; ++j )
                        {
                            loc_h[j] = 0.;
                        }
#    pragma omp for
                        for ( int k = 0; k < size_x; ++k )
                        {
                            const DataType x_temp = x[k];
                            const int i_index = i + k;
                            for ( int j = 0; j < size_h; ++j )
                            {
                                loc_h[j] += x_temp * temp ( i_index, i + j );
                            }
                        }
#    pragma omp critical
                        {
                            for ( int j = 0; j < size_h; ++j )
                            {
                                h[j] += loc_h[j];
                            }
                        }

                        delete [] loc_h;
                    }

                    // tmp(k:end,k:end) = temp(k:end,k:end) - 2*x*(x'*temp(k:end,k:end))
#    pragma omp parallel for shared(x,h)
                    for ( int j = i; j < nrows_; ++j )
                    {
                        const DataType x_temp = x[j - i];
                        for ( int k = i; k < ncols_; ++k )
                        {
                            temp ( j, k ) -= 2 * x_temp * h[k - i];
                        }
                    }

                    // save reflection vector
#    pragma omp parallel for shared(x)
                    for ( int j = 0; j < size_x; ++j )
                    {
                        V ( i, i + j ) = x[j];
                    }

                    delete [] h;
                }

                delete [] x;
            }

            // Compute Q
            for ( int i = ncols_ - 1; i >= 0; --i )
            {
                /* apply reflection */
                const int size_x = nrows_ - i;
                DataType *x = new DataType[size_x];
#    pragma omp parallel for shared(x)
                for ( int j = i; j < nrows_; ++j )
                {
                    x[j - i] = V ( i, j );
                }
                // precompute x'*Q(k:end,k:end)
                const int size_h = ncols_ - i;
                DataType *h = new DataType[size_h];
                memsethost ( h, 0, size_h, sizeof (DataType ) );

#    pragma omp parallel shared(x, h)
                {
                    DataType *loc_h = new DataType[size_h];
                    for ( int j = 0; j < size_h; ++j )
                    {
                        loc_h[j] = 0.;
                    }
#    pragma omp for
                    for ( int k = 0; k < size_x; ++k )
                    {
                        const DataType x_temp = x[k];
                        const int i_index = i + k;
                        for ( int j = 0; j < size_h; ++j )
                        {
                            loc_h[j] += x_temp * Q ( i_index, i + j );
                        }
                    }
#    pragma omp critical
                    {
                        for ( int j = 0; j < size_h; ++j )
                        {
                            h[j] += loc_h[j];
                        }
                    }

                    delete [] loc_h;
                }
                // Q(k:end,k:end) = Q(k:end,k:end) - 2*x*(x'*Q(k:end,k:end))
#    pragma omp parallel for shared(x,h)
                for ( int j = i; j < nrows_; ++j )
                {
                    DataType x_temp = x[j - i];
                    for ( int k = i; k < ncols_; ++k )
                    {
                        Q ( j, k ) -= 2 * x_temp * h[k - i];
                    }
                }

                delete [] x;
                delete [] h;
            }

            // Compute R = upper triangular of temp
#    pragma omp parallel for
            for ( int i = 0; i < ncols_; ++i )
            {
                if ( temp ( i, i ) >= 0 )
                {
                    for ( int j = i; j < ncols_; ++j )
                    {
                        R ( i, j ) = temp ( i, j );
                    }
                }
                else
                {
                    for ( int j = i; j < ncols_; ++j )
                    {
                        R ( i, j ) = -temp ( i, j );
                    }
                    for ( int j = 0; j < nrows_; ++j )
                    {
                        Q ( j, i ) = -Q ( j, i );
                    }
                }
            }

#else
            ERROR
#endif
        }

        /// template instantiation

        template class SeqDenseMatrix<double>;

        template class SeqDenseMatrix<float>;

    } // namespace la

} // namespace hiflow
