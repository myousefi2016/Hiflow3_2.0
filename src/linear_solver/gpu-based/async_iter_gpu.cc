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

#include "async_iter_gpu.h"
#include "linear_algebra/la_descriptor.h"

#ifdef WITH_CUDA
#    include "linear_algebra/lmp/cuda/cuda_async_iter.h"
#endif

#include "tools/mpi_tools.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        AsynchronousIterationGPU<LAD>::AsynchronousIterationGPU ( const MPI_Comm& comm )
        : LinearSolver<LAD>( ),
        A_ ( 0 ), res_vec_ ( 0 ),
        lmat_A_diag_ ( 0 ), lmat_A_offdiag_ ( 0 ),
        lvec_b_ ( 0 ), lvec_x_ ( 0 ), lvec_xg_ ( 0 ),
        b_host ( 0 ), x_host ( 0 ), xg_host ( 0 ),
        block_res_host ( 0 ),
        A_diag_gpu ( 0 ), A_offdiag_gpu ( 0 ),
        b_gpu ( 0 ), x_gpu ( 0 ), xg_gpu ( 0 ),
        inv_D_gpu ( 0 ), D_gpu ( 0 ), block_res_gpu ( 0 ),
#ifdef WITH_CUDA
        stream_ ( 0 ),
#endif
        n_gpus_ ( 0 ), used_gpus_ ( 0 ),
        enable_hyperq_ ( 0 ), enable_p2p_ ( 0 ),
        grid_dim_ ( 0 ), block_dim_ ( 128 ),
        n_global_iter_ ( 0 ), n_block_iter_ ( 0 ), n_inner_iter_ ( 0 ),
        solve_mode_ ( 0 ), max_threads_per_block_ ( 0 ),
        w_ ( 1.0 ),
        comm_ ( comm ), rank_ ( -1 ), nproc_ ( -1 ),
        row_offs_ ( 0 ),
        nnz_diag_ ( 0 ), row_diag_ ( 0 ), col_diag_ ( 0 ),
        nnz_offdiag_ ( 0 ), row_offdiag_ ( 0 ), col_offdiag_ ( 0 ),
        have_offdiag_ ( false ), compute_residual_on_dev_ ( false )
        {
            MPI_Comm_rank ( comm_, &rank_ );
            MPI_Comm_size ( comm_, &nproc_ );
            assert ( rank_ >= 0 );
            assert ( rank_ < nproc_ );

#ifdef WITH_CUDA
            // get number of devices
            cudaGetDeviceCount ( &n_gpus_ );
            if ( n_gpus_ <= 0 )
            {
                LOG_ERROR ( "AsynchronousIterationGPU rank " << rank_ << " couldn't detect any CUDA-capable device" );
                exit ( -1 );
            }
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "AsynchronousIterationGPU",
                           " rank " << rank_ << " detected " << n_gpus_ << " CUDA-capable devices" );
            }

            // get device properties
            cudaDeviceProp prop;
            int have_hyperq = 1;
            max_threads_per_block_ = 1000000;

            for ( int i = 0; i < n_gpus_; ++i )
            {
                cudaSetDevice ( i );
                cudaGetDeviceProperties ( &prop, i );

                max_threads_per_block_ = std::min ( max_threads_per_block_, prop.maxThreadsPerBlock );

                // check for compute capability
                if ( ( prop.major < 3 ) || ( ( prop.major == 3 ) && ( prop.minor < 5 ) ) || ( prop.major > 3 ) )
                {
                    // compute capability is less than 3.5
                    have_hyperq = 0;
                }
            }

            MPI_Allreduce ( &have_hyperq, &enable_hyperq_, 1, MPI_INT, MPI_MIN, comm_ );

            if ( enable_hyperq_ > 0 )
            {
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "AsynchronousIterationGPU",
                               " may use Hyper-Q" );
                }
            }
            else
            {
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "AsynchronousIterationGPU",
                               " Hyper-Q not available" );
                }
            }
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "AsynchronousIterationGPU",
                           " CUDA kernels can use max. " << max_threads_per_block_ << " threads per block" );
            }

#else
            LOG_ERROR ( "AsynchronousIterationGPU :  no CUDA-support" );
            exit ( -1 );
#endif
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::Prepare ( const MatrixType& A,
                                                      const VectorType& b,
                                                      VectorType& x,
                                                      const bool compute_residual_on_dev,
                                                      const int use_gpus,
                                                      const int procs_per_node )
        {
#ifdef WITH_CUDA

            assert ( procs_per_node > 0 );
            assert ( A.nrows_local ( ) == x.size_local ( ) );
            assert ( A.nrows_local ( ) == b.size_local ( ) );
            assert ( A.nrows_global ( ) == x.size_global ( ) );
            assert ( A.nrows_global ( ) == b.size_global ( ) );
            assert ( A.ncols_global ( ) == x.size_global ( ) );

            this->Clear ( );

            A_ = &A;
            lmat_A_diag_ = dynamic_cast < const CPU_CSR_lMatrix<DataType>* > ( &A.diagonal ( ) );
            lmat_A_offdiag_ = dynamic_cast < const CPU_CSR_lMatrix<DataType>* > ( &A.offdiagonal ( ) );
            lvec_b_ = dynamic_cast < const CPU_lVector<DataType>* > ( &b.interior ( ) );
            lvec_x_ = dynamic_cast < CPU_lVector<DataType>* > ( &x.interior ( ) );
            lvec_xg_ = dynamic_cast < CPU_lVector<DataType>* > ( &x.ghost ( ) );

            assert ( A_ != 0 );
            assert ( lmat_A_diag_ != 0 );
            assert ( lmat_A_offdiag_ != 0 );
            assert ( lvec_b_ != 0 );
            assert ( lvec_x_ != 0 );
            assert ( lvec_xg_ != 0 );

            compute_residual_on_dev_ = compute_residual_on_dev;
            used_gpus_ = use_gpus;

            if ( lmat_A_offdiag_->get_nnz ( ) > 0 ) have_offdiag_ = true;
            else have_offdiag_ = false;

            if ( used_gpus_ > n_gpus_ )
            {
                LOG_ERROR ( "Requested " << use_gpus << " devices but only " << n_gpus_ << " available." );
                exit ( -1 );
            }
            if ( A.diagonal ( ).get_matrix_format ( ) != CSR )
            {
                LOG_ERROR ( "Called AsynchronousIterationGPU::Prepare with non-CSR matrix argument." );
                exit ( -1 );
            }
            if ( A.offdiagonal ( ).get_matrix_format ( ) != CSR )
            {
                LOG_ERROR ( "Called AsynchronousIterationGPU::Prepare with non-CSR matrix argument." );
                exit ( -1 );
            }
            if ( have_offdiag_ && ( A.diagonal ( ).get_num_row ( ) != A.offdiagonal ( ).get_num_row ( ) ) )
            {
                LOG_ERROR ( "Number of rows in diagonal (" << A.diagonal ( ).get_num_row ( )
                            << ") and offdiagonal (" << A.offdiagonal ( ).get_num_row ( )
                            << ") part of matrix are not equal." )
                        exit ( -1 );
            }

            if ( ( procs_per_node > used_gpus_ ) && ( enable_hyperq_ == 0 ) )
            {
                LOG_ERROR ( "Requested for " << procs_per_node << " MPI processes to share " << used_gpus_ << " GPUs, but Hyper-Q is not enabled." );
                exit ( -1 );
            }

            // allocate one host memory pointer per GPU
            b_host = new DataType*[used_gpus_];
            x_host = new DataType*[used_gpus_];

            // allocate one GPU memory pointer per GPU
            A_diag_gpu = new DataType*[used_gpus_];
            b_gpu = new DataType*[used_gpus_];
            x_gpu = new DataType*[used_gpus_];
            inv_D_gpu = new DataType*[used_gpus_];

            if ( have_offdiag_ )
            {
                xg_host = new DataType*[used_gpus_];
                A_offdiag_gpu = new DataType*[used_gpus_];
                xg_gpu = new DataType*[used_gpus_];
            }
            else
            {
                xg_host = 0;
                A_offdiag_gpu = 0;
                xg_gpu = 0;
            }

            if ( compute_residual_on_dev_ )
            {
                D_gpu = new DataType*[used_gpus_];
                block_res_gpu = new DataType*[used_gpus_];
                block_res_host = new DataType[used_gpus_];
            }
            else
            {
                D_gpu = 0;
                block_res_gpu = 0;
                block_res_host = 0;
            }

            // create a stream for each GPU
            stream_ = new cudaStream_t[used_gpus_];
            for ( int i = 0; i < used_gpus_; ++i )
            {
                cudaSetDevice ( i );
                cudaStreamCreate ( &( stream_[i] ) );
            }

            // temporary arrays for host to GPU copy
            int tmpsize = std::max ( lmat_A_diag_->get_nnz ( ), lmat_A_offdiag_->get_nnz ( ) );
            int *tmpi = new int[tmpsize];
            DataType *tmpd = new DataType[tmpsize];

            // determine CSR structure for each GPU
            row_offs_ = new int[used_gpus_ + 1];
            row_offs_[0] = 0;
            nnz_diag_ = new int[used_gpus_];
            row_diag_ = new int*[used_gpus_];
            col_diag_ = new int*[used_gpus_];
            if ( have_offdiag_ )
            {
                nnz_offdiag_ = new int[used_gpus_];
                row_offdiag_ = new int*[used_gpus_];
                col_offdiag_ = new int*[used_gpus_];
            }
            else
            {
                nnz_offdiag_ = 0;
                row_offdiag_ = 0;
                col_offdiag_ = 0;
            }

            for ( int i = 0; i < used_gpus_; ++i )
            {
                cudaSetDevice ( i );

                // ==============================================================
                // determine number of rows per gpu, store as offsets
                int nrows = ( A.nrows_local ( ) / used_gpus_ );
                if ( i < ( A.nrows_local ( ) % used_gpus_ ) ) ++nrows;
                row_offs_[i + 1] = row_offs_[i] + nrows;

                // ================================================================
                // CSR row arrays for diagonal and offdiagonal part
                cudaMalloc ( ( void** ) ( &( row_diag_[i] ) ), ( nrows + 1 ) * sizeof (int ) );
                if ( have_offdiag_ )
                    cudaMalloc ( ( void** ) ( &( row_offdiag_[i] ) ), ( nrows + 1 ) * sizeof (int ) );

                memcpy ( tmpi,
                         &( lmat_A_diag_->matrix.row[row_offs_[i]] ),
                         ( nrows + 1 ) * sizeof (int ) );

                for ( int j = 0; j <= nrows; ++j )
                {
                    tmpi[j] -= lmat_A_diag_->matrix.row[row_offs_[i]] + j;
                }

                assert ( tmpi[0] == 0 );

                cudaMemcpy ( row_diag_[i], tmpi, ( nrows + 1 ) * sizeof (int ), cudaMemcpyHostToDevice );

                if ( have_offdiag_ )
                {
                    memcpy ( tmpi,
                             &( lmat_A_offdiag_->matrix.row[row_offs_[i]] ),
                             ( nrows + 1 ) * sizeof (int ) );

                    assert ( tmpi[0] == 0 );
                    cudaMemcpy ( row_offdiag_[i], tmpi, ( nrows + 1 ) * sizeof (int ), cudaMemcpyHostToDevice );
                }

                // ===================================================================
                // count nnz in diag part ignoring diagonal elements
                nnz_diag_[i] = 0;
                for ( int j = row_offs_[i]; j < row_offs_[i + 1]; ++j )
                {
                    for ( int k = lmat_A_diag_->matrix.row[j];
                          k < lmat_A_diag_->matrix.row[j + 1];
                          ++k )
                    {
                        if ( j != lmat_A_diag_->matrix.col[k] ) ++( nnz_diag_[i] );
                    }
                }
                assert ( nnz_diag_[i] ==
                         lmat_A_diag_->matrix.row[row_offs_[i + 1]] - lmat_A_diag_->matrix.row[row_offs_[i]] - nrows );

                // ===================================================================
                // count nnz in offdiag part (anyway no diagonal elements)
                if ( have_offdiag_ )
                {
                    nnz_offdiag_[i] = lmat_A_offdiag_->matrix.row[row_offs_[i + 1]] - lmat_A_offdiag_->matrix.row[row_offs_[i]];
                }

                // ===================================================================
                // allocate memory
                cudaMallocHost ( ( void** ) ( &( b_host[i] ) ), nrows * sizeof (DataType ) );
                cudaMallocHost ( ( void** ) ( &( x_host[i] ) ), x.size_local ( ) * sizeof (DataType ) );

                cudaMalloc ( ( void** ) ( &( col_diag_[i] ) ), ( nnz_diag_[i] ) * sizeof (int ) );
                cudaMalloc ( ( void** ) ( &( A_diag_gpu[i] ) ), ( nnz_diag_[i] ) * sizeof (DataType ) );
                cudaMalloc ( ( void** ) ( &( b_gpu[i] ) ), nrows * sizeof (DataType ) );
                cudaMalloc ( ( void** ) ( &( x_gpu[i] ) ), x.size_local ( ) * sizeof (DataType ) );
                cudaMalloc ( ( void** ) ( &( inv_D_gpu[i] ) ), nrows * sizeof (DataType ) );

                if ( have_offdiag_ )
                {
                    cudaMallocHost ( ( void** ) ( &( xg_host[i] ) ), x.size_local_ghost ( ) * sizeof (DataType ) );
                    cudaMalloc ( ( void** ) ( &( col_offdiag_[i] ) ), ( nnz_offdiag_[i] ) * sizeof (int ) );
                    cudaMalloc ( ( void** ) ( &( A_offdiag_gpu[i] ) ), ( nnz_offdiag_[i] ) * sizeof (DataType ) );
                    cudaMalloc ( ( void** ) ( &( xg_gpu[i] ) ), x.size_local_ghost ( ) * sizeof (DataType ) );
                }

                if ( compute_residual_on_dev_ )
                {
                    cudaMalloc ( ( void** ) ( &( D_gpu[i] ) ), nrows * sizeof (DataType ) );
                    cudaMalloc ( ( void** ) ( &( block_res_gpu[i] ) ), nrows * sizeof (DataType ) );
                    cudaMallocHost ( ( void** ) ( &( block_res_host[i] ) ), sizeof (DataType ) );
                }
                // ===================================================================
                // CSR col arrays for diagonal and offdiagonal parts without diagonal elements
                int ctr = 0;
                int inv_D_ctr = 0;
                for ( int j = row_offs_[i]; j < row_offs_[i + 1]; ++j )
                {
                    for ( int k = lmat_A_diag_->matrix.row[j];
                          k < lmat_A_diag_->matrix.row[j + 1];
                          ++k )
                    {
                        if ( j != lmat_A_diag_->matrix.col[k] )
                        {
                            tmpi[ctr] = lmat_A_diag_->matrix.col[k];
                            ++ctr;
                        }
                        else
                        {
                            tmpd[inv_D_ctr] = static_cast < DataType > ( 1.0 ) / lmat_A_diag_->matrix.val[k];
                            assert ( std::isnormal ( tmpd[inv_D_ctr] ) );
                            ++inv_D_ctr;
                        }
                    }
                }
                assert ( ctr == nnz_diag_[i] );
                assert ( inv_D_ctr == nrows );

                cudaMemcpy ( col_diag_[i], tmpi, ( nnz_diag_[i] ) * sizeof (int ), cudaMemcpyHostToDevice );
                cudaMemcpy ( inv_D_gpu[i], tmpd, nrows * sizeof (DataType ), cudaMemcpyHostToDevice );

                if ( compute_residual_on_dev_ )
                {
                    ctr = 0;
                    for ( int j = row_offs_[i]; j < row_offs_[i + 1]; ++j )
                    {
                        for ( int k = lmat_A_diag_->matrix.row[j];
                              k < lmat_A_diag_->matrix.row[j + 1];
                              ++k )
                        {
                            if ( j == lmat_A_diag_->matrix.col[k] )
                            {
                                tmpd[ctr] = lmat_A_diag_->matrix.val[k];
                                ++ctr;
                            }
                        }
                    }
                    assert ( ctr == nrows );
                    cudaMemcpy ( D_gpu[i], tmpd, nrows * sizeof (DataType ), cudaMemcpyHostToDevice );
                }

                if ( have_offdiag_ )
                {
                    int idx_start = lmat_A_offdiag_->matrix.row[row_offs_[i]];
                    memcpy ( tmpi, &( lmat_A_offdiag_->matrix.col[idx_start] ), ( nnz_offdiag_[i] ) * sizeof (int ) );
                    cudaMemcpy ( col_offdiag_[i], tmpi, ( nnz_offdiag_[i] ) * sizeof (int ), cudaMemcpyHostToDevice );
                }

                // ===================================================================
                // copy diagonal and offdiagonal matrix parts to GPU
                ctr = 0;
                for ( int j = row_offs_[i]; j < row_offs_[i + 1]; ++j )
                {
                    for ( int k = lmat_A_diag_->matrix.row[j];
                          k < lmat_A_diag_->matrix.row[j + 1];
                          ++k )
                    {
                        if ( j != lmat_A_diag_->matrix.col[k] )
                        {
                            tmpd[ctr] = lmat_A_diag_->matrix.val[k];
                            ++ctr;
                        }
                    }
                }

                assert ( ctr == nnz_diag_[i] );
                cudaMemcpy ( A_diag_gpu[i], tmpd, ( nnz_diag_[i] ) * sizeof (DataType ), cudaMemcpyHostToDevice );

                if ( have_offdiag_ )
                {
                    int idx_start = lmat_A_offdiag_->matrix.row[row_offs_[i]];
                    memcpy ( tmpd, &( lmat_A_offdiag_->matrix.val[idx_start] ), ( nnz_offdiag_[i] ) * sizeof (DataType ) );

                    cudaMemcpy ( A_offdiag_gpu[i], tmpd, ( nnz_offdiag_[i] ) * sizeof (DataType ), cudaMemcpyHostToDevice );
                }
            } // END for (int i = 0; i < used_gpus_; ++i)

            delete [] tmpi;
            delete [] tmpd;

#endif
        }

        template <class LAD>
        void AsynchronousIterationGPU<LAD>::CopyHostToDev_x ( )
        {
#ifdef WITH_CUDA
            assert ( lvec_x_ != 0 );

            // copy x and xg to page-locked host memory
            for ( int i = 0; i < used_gpus_; ++i )
            {
                memcpy ( x_host[i], lvec_x_->buffer, lvec_x_->get_size ( ) * sizeof (DataType ) );
                if ( have_offdiag_ )
                {
                    assert ( lvec_xg_ != 0 );
                    memcpy ( xg_host[i], lvec_xg_->buffer, lvec_xg_->get_size ( ) * sizeof (DataType ) );
                }
            }

            // async copy b, x and xg to GPUs
            for ( int i = 0; i < used_gpus_; ++i )
            {
                cudaSetDevice ( i );
                cudaMemcpyAsync ( x_gpu[i], x_host[i], lvec_x_->get_size ( ) * sizeof (DataType ),
                                  cudaMemcpyHostToDevice, stream_[i] );

                if ( have_offdiag_ )
                {
                    cudaMemcpyAsync ( xg_gpu[i], xg_host[i], lvec_xg_->get_size ( ) * sizeof (DataType ),
                                      cudaMemcpyHostToDevice, stream_[i] );
                }
            }
#endif
        }

        template <class LAD>
        void AsynchronousIterationGPU<LAD>::CopyHostToDev_b_x ( )
        {
#ifdef WITH_CUDA
            assert ( lvec_b_ != 0 );
            assert ( lvec_x_ != 0 );

            // copy b, x and xg to page-locked host memory
            for ( int i = 0; i < used_gpus_; ++i )
            {
                memcpy ( b_host[i], &( lvec_b_->buffer[row_offs_[i]] ), ( row_offs_[i + 1] - row_offs_[i] ) * sizeof (DataType ) );
                memcpy ( x_host[i], lvec_x_->buffer, lvec_x_->get_size ( ) * sizeof (DataType ) );
                if ( have_offdiag_ )
                {
                    assert ( lvec_xg_ != 0 );
                    memcpy ( xg_host[i], lvec_xg_->buffer, lvec_xg_->get_size ( ) * sizeof (DataType ) );
                }
            }

            // copy complete local and ghost vector to GPUs
            for ( int i = 0; i < used_gpus_; ++i )
            {
                cudaSetDevice ( i );
                cudaMemcpyAsync ( b_gpu[i], b_host[i], ( row_offs_[i + 1] - row_offs_[i] ) * sizeof (DataType ),
                                  cudaMemcpyHostToDevice, stream_[i] );
                cudaMemcpyAsync ( x_gpu[i], x_host[i], lvec_x_->get_size ( ) * sizeof (DataType ),
                                  cudaMemcpyHostToDevice, stream_[i] );
                if ( have_offdiag_ )
                {
                    cudaMemcpyAsync ( xg_gpu[i], xg_host[i], lvec_xg_->get_size ( ) * sizeof (DataType ),
                                      cudaMemcpyHostToDevice, stream_[i] );
                }
            }
#endif
        }

        template <class LAD>
        void AsynchronousIterationGPU<LAD>::CopyDevToHost_x ( )
        {
#ifdef WITH_CUDA
            assert ( lvec_x_ != 0 );

            // async copy of result x to page-locked host memory
            for ( int i = 0; i < used_gpus_; ++i )
            {
                cudaSetDevice ( i );
                cudaMemcpyAsync ( &( x_host[i][row_offs_[i]] ), &( x_gpu[i][row_offs_[i]] ),
                                  ( row_offs_[i + 1] - row_offs_[i] ) * sizeof (DataType ),
                                  cudaMemcpyDeviceToHost, stream_[i] );
            }

            // copy result x to vector
            for ( int i = 0; i < used_gpus_; ++i )
            {
                cudaSetDevice ( i );
                cudaStreamSynchronize ( stream_[i] );
            }

            for ( int i = 0; i < used_gpus_; ++i )
            {
                memcpy ( &( lvec_x_->buffer[row_offs_[i]] ), &( x_host[i][row_offs_[i]] ),
                         ( row_offs_[i + 1] - row_offs_[i] ) * sizeof (DataType ) );
            }

#endif
        }

        template <class LAD>
        void AsynchronousIterationGPU<LAD>::Run ( VectorType *x )
        {
#ifdef WITH_CUDA
            assert ( x != 0 );

            // do computations
            for ( int ii = 0; ii < n_global_iter_; ++ii )
            {
                for ( int jj = 0; jj < n_block_iter_; ++jj )
                {
                    // launch kernel on each GPU
                    for ( int i = 0; i < used_gpus_; ++i )
                    {
                        cudaSetDevice ( i );

                        assert ( block_dim_ > 0 );
                        const int nrows = row_offs_[i + 1] - row_offs_[i];
                        grid_dim_ = nrows / block_dim_;
                        if ( grid_dim_ * block_dim_ < nrows ) ++grid_dim_;

                        if ( have_offdiag_ )
                        {
                            cuda_async_inner_iter ( A_diag_gpu[i], row_diag_[i], col_diag_[i], nrows, row_offs_[i],
                                                    A_offdiag_gpu[i], row_offdiag_[i], col_offdiag_[i],
                                                    inv_D_gpu[i], b_gpu[i], x_gpu[i], xg_gpu[i], n_inner_iter_,
                                                    grid_dim_, block_dim_, stream_[i] );
                        }
                        else
                        {
                            cuda_async_inner_iter_no_offdiag ( A_diag_gpu[i], row_diag_[i], col_diag_[i], nrows, row_offs_[i],
                                                               inv_D_gpu[i], b_gpu[i], x_gpu[i], n_inner_iter_,
                                                               grid_dim_, block_dim_, stream_[i] );
                        }
                    }

                    this->iter_ += n_inner_iter_;

                    // update across GPUs used by this MPI proc
                    if ( ( used_gpus_ > 1 ) && ( jj < n_block_iter_ - 1 ) ) BlockUpdate ( );
                }
                // update across MPI procs
                if ( ( nproc_ > 1 ) && ( ii < n_global_iter_ - 1 ) ) GlobalUpdate ( x );
            }
#endif
        }

        template <class LAD>
        void AsynchronousIterationGPU<LAD>::RunDamped ( VectorType *x )
        {
#ifdef WITH_CUDA
            assert ( x != 0 );
            assert ( w_ > 0.0 );

            // do computations
            for ( int ii = 0; ii < n_global_iter_; ++ii )
            {
                for ( int jj = 0; jj < n_block_iter_; ++jj )
                {
                    // launch kernel on each GPU
                    for ( int i = 0; i < used_gpus_; ++i )
                    {
                        cudaSetDevice ( i );

                        assert ( block_dim_ > 0 );
                        const int nrows = row_offs_[i + 1] - row_offs_[i];
                        grid_dim_ = nrows / block_dim_;
                        if ( grid_dim_ * block_dim_ < nrows ) ++grid_dim_;

                        if ( have_offdiag_ )
                        {
                            cuda_async_inner_damped_iter ( A_diag_gpu[i], row_diag_[i], col_diag_[i], nrows, row_offs_[i],
                                                           A_offdiag_gpu[i], row_offdiag_[i], col_offdiag_[i],
                                                           inv_D_gpu[i], b_gpu[i], x_gpu[i], xg_gpu[i], w_, n_inner_iter_,
                                                           grid_dim_, block_dim_, stream_[i] );
                        }
                        else
                        {
                            cuda_async_inner_damped_iter_no_offdiag ( A_diag_gpu[i], row_diag_[i], col_diag_[i], nrows, row_offs_[i],
                                                                      inv_D_gpu[i], b_gpu[i], x_gpu[i], w_, n_inner_iter_,
                                                                      grid_dim_, block_dim_, stream_[i] );
                        }
                    }

                    this->iter_ += n_inner_iter_;

                    // update inner
                    if ( jj < n_block_iter_ - 1 ) BlockUpdate ( );
                }
                if ( ii < n_global_iter_ - 1 ) GlobalUpdate ( x );
            }
#endif
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::ComputeResidualOnDev ( )
        {
            DataType local_res = ComputeBlockSquaredResidualOnDev ( );
            DataType res;

            MPI_Allreduce ( &local_res, &res, 1, mpi_data_type<DataType>::get_type ( ), MPI_SUM, comm_ );

            this->res_ = std::sqrt ( res );
        }

        template <class LAD>
        typename LAD::DataType AsynchronousIterationGPU<LAD>::ComputeBlockSquaredResidualOnDev ( )
        {
#ifdef WITH_CUDA
            // launch kernel on each GPU
            for ( int i = 0; i < used_gpus_; ++i )
            {
                cudaSetDevice ( i );

                assert ( block_dim_ > 0 );
                const int nrows = row_offs_[i + 1] - row_offs_[i];
                grid_dim_ = nrows / block_dim_;
                if ( grid_dim_ * block_dim_ < nrows ) ++grid_dim_;

                if ( have_offdiag_ )
                {
                    cuda_compute_block_squared_residual
                            ( A_diag_gpu[i], row_diag_[i], col_diag_[i], nrows, row_offs_[i],
                              A_offdiag_gpu[i], row_offdiag_[i], col_offdiag_[i],
                              D_gpu[i], b_gpu[i], x_gpu[i], xg_gpu[i], block_res_gpu[i],
                              grid_dim_, block_dim_, max_threads_per_block_,
                              stream_[i] );
                }
                else
                {
                    cuda_compute_block_squared_residual_no_offdiag
                            ( A_diag_gpu[i], row_diag_[i], col_diag_[i], nrows, row_offs_[i],
                              D_gpu[i], b_gpu[i], x_gpu[i], block_res_gpu[i],
                              grid_dim_, block_dim_, max_threads_per_block_,
                              stream_[i] );
                }

                cudaMemcpyAsync ( &( block_res_host[i] ), block_res_gpu[i], sizeof (DataType ), cudaMemcpyDeviceToHost, stream_[i] );
            }

            for ( int i = 0; i < used_gpus_; ++i )
            {
                cudaSetDevice ( i );
                cudaStreamSynchronize ( stream_[i] );
            }

            DataType result = 0.0;
            for ( int i = 0; i < used_gpus_; ++i )
            {
                result += block_res_host[i];
            }

            return result;
#else
            return 0.0;
#endif
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::SetNumIter ( const int n_global,
                                                         const int n_block, // = 1
                                                         const int n_inner ) // = 5
        {
            assert ( n_global > 0 );
            assert ( n_block > 0 );
            assert ( n_inner > 0 );

            n_global_iter_ = n_global;
            n_block_iter_ = n_block;
            n_inner_iter_ = n_inner;
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::SetSolveMode ( const int mode )
        {
            assert ( mode >= 0 );
            assert ( mode <= 3 );

            // 0: solve normal (check for residual norm after each run)
            // 1: solve damped (check for residual norm after each damped run)
            // 2: smooth normal (performs one run, no residual norm check)
            // 3: smooth damped (performs one damped run, no residual norm check)

            solve_mode_ = mode;
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::SetDampingParameter ( const DataType w )
        {
            assert ( w > 0 );
            w_ = w;
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::SetCudaBlockSize ( const int block_size )
        {
            assert ( block_size > 0 );

            block_dim_ = block_size;
        }

        template<class LAD>
        LinearSolverState AsynchronousIterationGPU<LAD>::Solve ( const VectorType& b, VectorType *x )
        {
            assert ( block_dim_ > 0 );

            this->iter_ = 0;

            lvec_b_ = dynamic_cast < const CPU_lVector<DataType>* > ( &( b.interior ( ) ) );
            lvec_x_ = dynamic_cast < CPU_lVector<DataType>* > ( &( x->interior ( ) ) );
            lvec_xg_ = dynamic_cast < CPU_lVector<DataType>* > ( &( x->ghost ( ) ) );

            assert ( A_ != 0 );
            assert ( lvec_b_ != 0 );
            assert ( lvec_x_ != 0 );
            assert ( lvec_xg_ != 0 );

            switch ( solve_mode_ )
            {
                case 0:
                    return SolveNormal ( b, x );
                    break;

                case 1:
                    return SolveDamped ( b, x );
                    break;

                case 2:
                    return SmoothNormal ( b, x );
                    break;

                case 3:
                    return SmoothDamped ( b, x );
                    break;

                default:
                    LOG_ERROR ( "No valid solver mode: " << solve_mode_ );
                    return kSolverError;
                    break;
            }
        }

        template<class LAD>
        LinearSolverState AsynchronousIterationGPU<LAD>::SolveNormal ( const VectorType& b, VectorType *x )
        {
#ifdef WITH_CUDA
            this->iter_ = 0;

            x->UpdateGhost ( );
            CopyHostToDev_b_x ( );

            if ( compute_residual_on_dev_ )
            {
                ComputeResidualOnDev ( );
            }
            else
            {
                if ( res_vec_ == 0 )
                {
                    res_vec_ = new VectorType;
                    res_vec_->CloneFromWithoutContent ( b );
                }

                // compute initial residual
                A_->VectorMult ( *x, res_vec_ );
                res_vec_->ScaleAdd ( b, -1.0 );
                this->res_ = res_vec_->Norm2 ( );
            }

            IterateControl::State iterctrl = this->control ( ).Check ( this->iter_, this->res_ );

            while ( iterctrl == IterateControl::kIterate )
            {
                Run ( x );

                if ( compute_residual_on_dev_ )
                {
                    GlobalUpdate ( x );
                    ComputeResidualOnDev ( );
                }
                else
                {
                    CopyDevToHost_x ( );

                    A_->VectorMult ( *x, res_vec_ );
                    res_vec_->ScaleAdd ( b, -1.0 );
                    this->res_ = res_vec_->Norm2 ( );

                    CopyHostToDev_x ( );
                }

                iterctrl = this->control ( ).Check ( this->iter_, this->res_ );
            }
            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "AsynchronousIterationGPU",
                           " normal solve stops after " << this->iter_ << " iterations with residual norm " << this->res_ );
            }

            if ( res_vec_ != 0 ) delete res_vec_;

            return kSolverSuccess;
#else
            return kSolverError;
#endif
        }

        template<class LAD>
        LinearSolverState AsynchronousIterationGPU<LAD>::SolveDamped ( const VectorType& b,
                                                                       VectorType *x )
        {

#ifdef WITH_CUDA
            this->iter_ = 0;

            x->UpdateGhost ( );
            CopyHostToDev_b_x ( );

            if ( compute_residual_on_dev_ )
            {
                ComputeResidualOnDev ( );
            }
            else
            {
                if ( res_vec_ == 0 )
                {
                    res_vec_ = new VectorType;
                    res_vec_->CloneFromWithoutContent ( b );
                }

                // compute initial residual
                A_->VectorMult ( *x, res_vec_ );
                res_vec_->ScaleAdd ( b, -1.0 );
                this->res_ = res_vec_->Norm2 ( );
            }

            IterateControl::State iterctrl = this->control ( ).Check ( this->iter_, this->res_ );

            while ( iterctrl == IterateControl::kIterate )
            {
                RunDamped ( x );

                if ( compute_residual_on_dev_ )
                {
                    GlobalUpdate ( x );
                    ComputeResidualOnDev ( );
                }
                else
                {
                    CopyDevToHost_x ( );

                    A_->VectorMult ( *x, res_vec_ );
                    res_vec_->ScaleAdd ( b, -1.0 );
                    this->res_ = res_vec_->Norm2 ( );

                    CopyHostToDev_x ( );
                }

                iterctrl = this->control ( ).Check ( this->iter_, this->res_ );
            }
            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "AsynchronousIterationGPU",
                           " damped solve stops after " << this->iter_ << " iterations with residual norm " << this->res_ );
            }

            if ( res_vec_ != 0 ) delete res_vec_;

            return kSolverSuccess;
#else
            return kSolverError;
#endif
        }

        template<class LAD>
        LinearSolverState AsynchronousIterationGPU<LAD>::SmoothNormal ( const VectorType& b, VectorType* x )
        {
#ifdef WITH_CUDA
            x->UpdateGhost ( );
            this->CopyHostToDev_b_x ( );
            this->Run ( x );
            this->CopyDevToHost_x ( );

            return kSolverSuccess;
#else
            return kSolverError;
#endif
        }

        template<class LAD>
        LinearSolverState AsynchronousIterationGPU<LAD>::SmoothDamped ( const VectorType& b, VectorType* x )
        {

#ifdef WITH_CUDA  
            x->UpdateGhost ( );
            this->CopyHostToDev_b_x ( );
            this->RunDamped ( x );
            this->CopyDevToHost_x ( );

            return kSolverSuccess;
#else
            return kSolverError;
#endif
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::BlockUpdate ( )
        {
#ifdef WITH_CUDA
            if ( used_gpus_ > 1 )
            {
                if ( enable_p2p_ )
                {
                    LOG_ERROR ( "P2P block update not yet implemented" );
                    exit ( -1 );
                }
                else
                {
                    // copy local vector parts from GPUs to host
                    for ( int i = 0; i < used_gpus_; ++i )
                    {
                        cudaSetDevice ( i );
                        cudaMemcpyAsync ( &( x_host[i][row_offs_[i]] ), &( x_gpu[i][row_offs_[i]] ),
                                          ( row_offs_[i + 1] - row_offs_[i] ) * sizeof (DataType ),
                                          cudaMemcpyDeviceToHost, stream_[i] );
                    }

                    // wait for transfers to complete
                    for ( int i = 0; i < used_gpus_; ++i )
                    {
                        cudaSetDevice ( i );
                        cudaStreamSynchronize ( stream_[i] );
                    }

                    // update inner parts i <=> j of x_host[i] and x_host[j]
                    for ( int i = 0; i < used_gpus_; ++i )
                    {
                        for ( int j = 0; j < used_gpus_; ++j )
                        {
                            if ( i != j )
                            {
                                memcpy ( &( x_host[j][row_offs_[i]] ), &( x_host[i][row_offs_[i]] ),
                                         ( row_offs_[i + 1] - row_offs_[i] ) * sizeof (DataType ) );
                            }
                        }
                    }

                    // copy complete local vector to GPUs
                    for ( int i = 0; i < used_gpus_; ++i )
                    {
                        cudaSetDevice ( i );
                        cudaMemcpyAsync ( x_gpu[i], x_host[i], lvec_x_->get_size ( ) * sizeof (DataType ),
                                          cudaMemcpyHostToDevice, stream_[i] );
                    }
                }
            }
#endif
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::GlobalUpdate ( VectorType *x )
        {
            if ( nproc_ > 1 )
            {
                CopyDevToHost_x ( );
                x->UpdateGhost ( );
                CopyHostToDev_x ( );
            }
            else
            {
                BlockUpdate ( );
            }
        }

        template<class LAD>
        AsynchronousIterationGPU<LAD>::~AsynchronousIterationGPU ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void AsynchronousIterationGPU<LAD>::Clear ( )
        {
#ifdef WITH_CUDA

            if ( stream_ != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    cudaStreamSynchronize ( stream_[i] );
                    cudaStreamDestroy ( stream_[i] );
                }
                delete [] stream_;
            }
            stream_ = 0;

            A_ = 0;
            lmat_A_diag_ = 0;
            lmat_A_offdiag_ = 0;
            lvec_b_ = 0;
            lvec_x_ = 0;
            lvec_xg_ = 0;

            if ( res_vec_ != 0 ) delete res_vec_;
            res_vec_ = 0;

            // free host memory  
            if ( b_host != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( b_host[i] != 0 ) cudaFreeHost ( b_host[i] );
                }
                delete [] b_host;
            }
            b_host = 0;

            if ( x_host != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( x_host[i] != 0 ) cudaFreeHost ( x_host[i] );
                }
                delete [] x_host;
            }
            x_host = 0;

            if ( xg_host != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( xg_host[i] != 0 ) cudaFreeHost ( xg_host[i] );
                }
                delete [] xg_host;
            }
            xg_host = 0;

            if ( block_res_host != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    cudaFreeHost ( &( block_res_host[i] ) );
                }
                delete [] block_res_host;
            }
            block_res_host = 0;

            // free GPU memory
            if ( A_diag_gpu != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( A_diag_gpu[i] != 0 ) cudaFree ( A_diag_gpu[i] );
                }
                delete [] A_diag_gpu;
            }
            A_diag_gpu = 0;

            if ( A_offdiag_gpu != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( A_offdiag_gpu[i] != 0 ) cudaFree ( A_offdiag_gpu[i] );
                }
                delete [] A_offdiag_gpu;
            }
            A_offdiag_gpu = 0;

            if ( b_gpu != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( b_gpu[i] != 0 ) cudaFree ( b_gpu[i] );
                }
                delete [] b_gpu;
            }
            b_gpu = 0;

            if ( x_gpu != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( x_gpu[i] != 0 ) cudaFree ( x_gpu[i] );
                }
                delete [] x_gpu;
            }
            x_gpu = 0;

            if ( xg_gpu != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( xg_gpu[i] != 0 ) cudaFree ( xg_gpu[i] );
                }
                delete [] xg_gpu;
            }
            xg_gpu = 0;

            if ( inv_D_gpu != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( inv_D_gpu[i] != 0 ) cudaFree ( inv_D_gpu[i] );
                }
                delete [] inv_D_gpu;
            }
            inv_D_gpu = 0;

            if ( D_gpu != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( D_gpu[i] != 0 ) cudaFree ( D_gpu[i] );
                }
                delete [] D_gpu;
            }
            D_gpu = 0;

            if ( block_res_gpu != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( block_res_gpu[i] != 0 ) cudaFree ( block_res_gpu[i] );
                }
                delete [] block_res_gpu;
            }
            block_res_gpu = 0;

            if ( row_diag_ != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( row_diag_[i] != 0 ) cudaFree ( row_diag_[i] );
                }
                delete [] row_diag_;
            }
            row_diag_ = 0;

            if ( col_diag_ != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( col_diag_[i] != 0 ) cudaFree ( col_diag_[i] );
                }
                delete [] col_diag_;
            }
            col_diag_ = 0;

            if ( row_offdiag_ != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( row_offdiag_[i] != 0 ) cudaFree ( row_offdiag_[i] );
                }
                delete [] row_offdiag_;
            }
            row_offdiag_ = 0;

            if ( col_offdiag_ != 0 )
            {
                for ( int i = 0; i < used_gpus_; ++i )
                {
                    cudaSetDevice ( i );
                    if ( col_offdiag_[i] != 0 ) cudaFree ( col_offdiag_[i] );
                }
                delete [] col_offdiag_;
            }
            col_offdiag_ = 0;

            if ( nnz_diag_ != 0 ) delete [] nnz_diag_;
            if ( nnz_offdiag_ != 0 ) delete [] nnz_offdiag_;
            if ( row_offs_ != 0 ) delete [] row_offs_;
            row_offs_ = 0;
            nnz_diag_ = 0;
            nnz_offdiag_ = 0;

            used_gpus_ = 0;
#endif
        }

        /// template instantiation
        template class AsynchronousIterationGPU<LADescriptorCoupledD>;
        template class AsynchronousIterationGPU<LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow
