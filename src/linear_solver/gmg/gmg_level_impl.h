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

#include "gmg_level.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            template<class LAD, class ConnectionType>
            void GMGLevel<LAD, ConnectionType>::prepare_grid_transfer_on_device ( LinearSolver<LAD> * smoother )
            {
                if ( this->is_scheduled_to_this_process ( ) )
                {
                    assert ( this->matrix ( ) );
                    assert ( this->restriction_matrix ( ) );
                    assert ( this->interpolation_matrix ( ) );
                    assert ( this->matrix ( )->nrows_local ( ) == this->restriction_matrix ( )->nrows_local ( ) );
                    assert ( this->matrix ( )->nrows_local ( ) == this->interpolation_matrix ( )->nrows_local ( ) );
                    assert ( this->sol ( ) );
                    assert ( this->rhs ( ) );
                    assert ( this->res ( ) );

#ifdef WITH_CUDA
                    int ngpus = 0;
                    cudaGetDeviceCount ( &ngpus );
                    if ( ngpus <= 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): No CUDA-capable device found!" );
                        exit ( -1 );
                    }
                    cudaSetDevice ( 0 );
                    cudaStreamCreate ( &stream_ );

                    const CPU_CSR_lMatrix<DataType> * lmat_A_diag
                            = dynamic_cast < const CPU_CSR_lMatrix<DataType>* > ( &( this->matrix ( )->diagonal ( ) ) );
                    if ( lmat_A_diag == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): System matrix diagonal part is not CPU CSR." );
                        exit ( -1 );
                    }

                    const CPU_CSR_lMatrix<DataType> * lmat_A_offdiag
                            = dynamic_cast < const CPU_CSR_lMatrix<DataType>* > ( &( this->matrix ( )->offdiagonal ( ) ) );
                    if ( lmat_A_offdiag == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): System matrix offdiagonal part is not CPU CSR." );
                        exit ( -1 );
                    }

                    const CPU_CSR_lMatrix<DataType> * lmat_R_diag
                            = dynamic_cast < const CPU_CSR_lMatrix<DataType>* > ( &( this->restriction_matrix ( )->diagonal ( ) ) );
                    if ( lmat_R_diag == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): Restriction matrix diagonal part is not CPU CSR." );
                        exit ( -1 );
                    }

                    const CPU_CSR_lMatrix<DataType> * lmat_R_offdiag
                            = dynamic_cast < const CPU_CSR_lMatrix<DataType>* > ( &( this->restriction_matrix ( )->offdiagonal ( ) ) );
                    if ( lmat_R_offdiag == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): Restriction matrix offdiagonal part is not CPU CSR." );
                        exit ( -1 );
                    }

                    const CPU_CSR_lMatrix<DataType> * lmat_P_diag
                            = dynamic_cast < const CPU_CSR_lMatrix<DataType>* > ( &( this->interpolation_matrix ( )->diagonal ( ) ) );
                    if ( lmat_P_diag == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): Interpolation matrix diagonal part is not CPU CSR." );
                        exit ( -1 );
                    }

                    const CPU_CSR_lMatrix<DataType> * lmat_P_offdiag
                            = dynamic_cast < const CPU_CSR_lMatrix<DataType>* > ( &( this->interpolation_matrix ( )->offdiagonal ( ) ) );
                    if ( lmat_P_offdiag == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): Interpolation matrix offdiagonal part is not CPU CSR." );
                        exit ( -1 );
                    }

                    lvec_sol_ = dynamic_cast < CPU_lVector<DataType>* > ( &( this->sol ( )->interior ( ) ) );
                    if ( lvec_sol_ == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): Solution vector interior part is not on CPU." );
                        exit ( -1 );
                    }

                    lvec_sol_ghost_ = dynamic_cast < CPU_lVector<DataType>* > ( &( this->sol ( )->ghost ( ) ) );
                    if ( lvec_sol_ghost_ == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): Solution vector ghost part is not on CPU." );
                        exit ( -1 );
                    }

                    lvec_rhs_ = dynamic_cast < CPU_lVector<DataType>* > ( &( this->rhs ( )->interior ( ) ) );
                    if ( lvec_rhs_ == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): RHS vector interior part is not on CPU." );
                        exit ( -1 );
                    }

                    lvec_rhs_ghost_ = dynamic_cast < CPU_lVector<DataType>* > ( &( this->rhs ( )->ghost ( ) ) );
                    if ( lvec_rhs_ghost_ == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): RHS vector ghost part is not on CPU." );
                        exit ( -1 );
                    }

                    lvec_res_ = dynamic_cast < CPU_lVector<DataType>* > ( &( this->res ( )->interior ( ) ) );
                    if ( lvec_res_ == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): Residual vector interoir part is not on CPU." );
                        exit ( -1 );
                    }

                    lvec_res_ghost_ = dynamic_cast < CPU_lVector<DataType>* > ( &( this->res ( )->ghost ( ) ) );
                    if ( lvec_res_ghost_ == 0 )
                    {
                        LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): Residual vector ghost part is not on CPU." );
                        exit ( -1 );
                    }

                    async_iter_gpu_ = dynamic_cast < AsynchronousIterationGPU<LAD>* > ( smoother );
                    if ( async_iter_gpu_ != 0 )
                    {
                        if ( async_iter_gpu_->n_gpus ( ) > 1 ) async_iter_gpu_ = 0;
                    }

                    // restriction matrix R
                    nnz_R_diag_ = this->restriction_matrix ( )->nnz_local_diag ( );
                    nrows_ = this->restriction_matrix ( )->nrows_local ( );

                    nnz_R_offdiag_ = this->restriction_matrix ( )->nnz_local_offdiag ( );

                    cudaMalloc ( ( void** ) &R_diag_, nnz_R_diag_ * sizeof (DataType ) );
                    cudaMalloc ( ( void** ) &col_R_diag_, nnz_R_diag_ * sizeof (int ) );
                    cudaMalloc ( ( void** ) &row_R_diag_, ( nrows_ + 1 ) * sizeof (int ) );

                    cudaMemcpy ( R_diag_, lmat_R_diag->matrix.val, nnz_R_diag_ * sizeof (DataType ), cudaMemcpyHostToDevice );
                    cudaMemcpy ( col_R_diag_, lmat_R_diag->matrix.col, nnz_R_diag_ * sizeof (int ), cudaMemcpyHostToDevice );
                    cudaMemcpy ( row_R_diag_, lmat_R_diag->matrix.row, ( nrows_ + 1 ) * sizeof (int ), cudaMemcpyHostToDevice );

                    if ( nnz_R_offdiag_ > 0 )
                    {
                        assert ( nrows_ == lmat_R_offdiag->get_num_row ( ) );

                        cudaMalloc ( ( void** ) &R_offdiag_, nnz_R_offdiag_ * sizeof (DataType ) );
                        cudaMalloc ( ( void** ) &col_R_offdiag_, nnz_R_offdiag_ * sizeof (int ) );
                        cudaMalloc ( ( void** ) &row_R_offdiag_, ( nrows_ + 1 ) * sizeof (int ) );

                        cudaMemcpy ( R_offdiag_, lmat_R_offdiag->matrix.val, nnz_R_offdiag_ * sizeof (DataType ), cudaMemcpyHostToDevice );
                        cudaMemcpy ( col_R_offdiag_, lmat_R_offdiag->matrix.col, nnz_R_offdiag_ * sizeof (int ), cudaMemcpyHostToDevice );
                        cudaMemcpy ( row_R_offdiag_, lmat_R_offdiag->matrix.row, ( nrows_ + 1 ) * sizeof (int ), cudaMemcpyHostToDevice );
                    }

                    // system matrix A
                    nnz_A_diag_ = this->matrix ( )->nnz_local_diag ( );
                    assert ( nrows_ == this->matrix ( )->nrows_local ( ) );

                    nnz_A_offdiag_ = this->matrix ( )->nnz_local_offdiag ( );

                    cudaMalloc ( ( void** ) &A_diag_, nnz_A_diag_ * sizeof (DataType ) );
                    cudaMalloc ( ( void** ) &col_A_diag_, nnz_A_diag_ * sizeof (int ) );
                    cudaMalloc ( ( void** ) &row_A_diag_, ( nrows_ + 1 ) * sizeof (int ) );

                    cudaMemcpy ( A_diag_, lmat_A_diag->matrix.val, nnz_A_diag_ * sizeof (DataType ), cudaMemcpyHostToDevice );
                    cudaMemcpy ( col_A_diag_, lmat_A_diag->matrix.col, nnz_A_diag_ * sizeof (int ), cudaMemcpyHostToDevice );
                    cudaMemcpy ( row_A_diag_, lmat_A_diag->matrix.row, ( nrows_ + 1 ) * sizeof (int ), cudaMemcpyHostToDevice );

                    if ( nnz_A_offdiag_ > 0 )
                    {
                        assert ( nrows_ == lmat_A_offdiag->get_num_row ( ) );

                        cudaMalloc ( ( void** ) &A_offdiag_, nnz_A_offdiag_ * sizeof (DataType ) );
                        cudaMalloc ( ( void** ) &col_A_offdiag_, nnz_A_offdiag_ * sizeof (int ) );
                        cudaMalloc ( ( void** ) &row_A_offdiag_, ( nrows_ + 1 ) * sizeof (int ) );

                        cudaMemcpy ( A_offdiag_, lmat_A_offdiag->matrix.val, nnz_A_offdiag_ * sizeof (DataType ), cudaMemcpyHostToDevice );
                        cudaMemcpy ( col_A_offdiag_, lmat_A_offdiag->matrix.col, nnz_A_offdiag_ * sizeof (int ), cudaMemcpyHostToDevice );
                        cudaMemcpy ( row_A_offdiag_, lmat_A_offdiag->matrix.row, ( nrows_ + 1 ) * sizeof (int ), cudaMemcpyHostToDevice );
                    }

                    // interpolation matrix P
                    nnz_P_diag_ = this->interpolation_matrix ( )->nnz_local_diag ( );
                    assert ( nrows_ == this->interpolation_matrix ( )->nrows_local ( ) );

                    nnz_P_offdiag_ = this->interpolation_matrix ( )->nnz_local_offdiag ( );

                    cudaMalloc ( ( void** ) &P_diag_, nnz_P_diag_ * sizeof (DataType ) );
                    cudaMalloc ( ( void** ) &col_P_diag_, nnz_P_diag_ * sizeof (int ) );
                    cudaMalloc ( ( void** ) &row_P_diag_, ( nrows_ + 1 ) * sizeof (int ) );

                    cudaMemcpy ( P_diag_, lmat_P_diag->matrix.val, nnz_P_diag_ * sizeof (DataType ), cudaMemcpyHostToDevice );
                    cudaMemcpy ( col_P_diag_, lmat_P_diag->matrix.col, nnz_P_diag_ * sizeof (int ), cudaMemcpyHostToDevice );
                    cudaMemcpy ( row_P_diag_, lmat_P_diag->matrix.row, ( nrows_ + 1 ) * sizeof (int ), cudaMemcpyHostToDevice );

                    if ( nnz_P_offdiag_ > 0 )
                    {
                        assert ( nrows_ == lmat_P_offdiag->get_num_row ( ) );

                        cudaMalloc ( ( void** ) &P_offdiag_, nnz_P_offdiag_ * sizeof (DataType ) );
                        cudaMalloc ( ( void** ) &col_P_offdiag_, nnz_P_offdiag_ * sizeof (int ) );
                        cudaMalloc ( ( void** ) &row_P_offdiag_, ( nrows_ + 1 ) * sizeof (int ) );

                        cudaMemcpy ( P_offdiag_, lmat_P_offdiag->matrix.val, nnz_P_offdiag_ * sizeof (DataType ), cudaMemcpyHostToDevice );
                        cudaMemcpy ( col_P_offdiag_, lmat_P_offdiag->matrix.col, nnz_P_offdiag_ * sizeof (int ), cudaMemcpyHostToDevice );
                        cudaMemcpy ( row_P_offdiag_, lmat_P_offdiag->matrix.row, ( nrows_ + 1 ) * sizeof (int ), cudaMemcpyHostToDevice );
                    }

                    // solution vector
                    if ( async_iter_gpu_ != 0 )
                    {
                        // use device and host memory from AsynchronousIterationGPU smoother
                        dev_sol_ = async_iter_gpu_-> dev_x_ptr ( 0 );
                        host_sol_ = async_iter_gpu_->host_x_ptr ( 0 );
                        dev_sol_ghost_ = async_iter_gpu_-> dev_xg_ptr ( 0 );
                        host_sol_ghost_ = async_iter_gpu_->host_xg_ptr ( 0 );
                    }
                    else
                    {

                        cudaMalloc ( ( void** ) &dev_sol_, nrows_ * sizeof (DataType ) );
                        cudaMallocHost ( ( void** ) &host_sol_, nrows_ * sizeof (DataType ) );
                        cudaMalloc ( ( void** ) &dev_sol_ghost_, lvec_sol_ghost_->get_size ( ) * sizeof (DataType ) );
                        cudaMallocHost ( ( void** ) &host_sol_ghost_, lvec_sol_ghost_->get_size ( ) * sizeof (DataType ) );
                    }

                    // rhs vector
                    if ( async_iter_gpu_ != 0 )
                    {
                        // use device and host memory from AsynchronousIterationGPU smoother
                        dev_rhs_ = async_iter_gpu_-> dev_b_ptr ( 0 );
                        host_rhs_ = async_iter_gpu_->host_b_ptr ( 0 );
                    }
                    else
                    {
                        cudaMalloc ( ( void** ) &dev_rhs_, nrows_ * sizeof (DataType ) );
                        cudaMallocHost ( ( void** ) &host_rhs_, nrows_ * sizeof (DataType ) );
                    }
                    cudaMalloc ( ( void** ) &dev_rhs_ghost_, lvec_rhs_ghost_->get_size ( ) * sizeof (DataType ) );
                    cudaMallocHost ( ( void** ) &host_rhs_ghost_, lvec_rhs_ghost_->get_size ( ) * sizeof (DataType ) );

                    // residual vector
                    cudaMalloc ( ( void** ) &dev_res_, nrows_ * sizeof (DataType ) );
                    cudaMallocHost ( ( void** ) &host_res_, nrows_ * sizeof (DataType ) );
                    cudaMalloc ( ( void** ) &dev_res_ghost_, lvec_res_ghost_->get_size ( ) * sizeof (DataType ) );
                    cudaMallocHost ( ( void** ) &host_res_ghost_, lvec_res_ghost_->get_size ( ) * sizeof (DataType ) );
                    cudaMalloc ( ( void** ) &dev_tmp_, nrows_ * sizeof (DataType ) );

                    grid_dim_ = nrows_ / block_dim_;
                    if ( grid_dim_ * block_dim_ < nrows_ ) ++grid_dim_;

                    grid_transfer_on_device_ = true;

#else
                    LOG_ERROR ( "BasicLevel<LAD>::prepare_grid_transfer_on_device(): No CUDA support!" );
#endif
                }
            }

        } // namespace gmg
    } // namespace la
} // namespace hiflow