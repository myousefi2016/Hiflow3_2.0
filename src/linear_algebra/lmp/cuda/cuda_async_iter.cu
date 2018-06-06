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

#include "cuda_async_iter.h"
#include "linear_algebra/lmp/lmp_log.h"

#ifndef NOCUDA
#include "cuda_async_iter_kernels.h"
#endif

template <typename DataType>
void cuda_async_inner_iter(const DataType *A_diag,
                      const int *row_diag,
                      const int *col_diag,
                      const int nrows,
                      const int row_offset,
                      const DataType *A_offdiag,
                      const int *row_offdiag,
                      const int *col_offdiag,
                      const DataType *inv_D,
                      const DataType *b,
                      DataType *x,
                      const DataType *xg,
                      const int n_iter,
                      const int dim_grid,
                      const int dim_block,
                      cudaStream_t stream)
{
#ifndef NOCUDA
  dim3 dimGrid(dim_grid);
  dim3 dimBlock(dim_block);

  kernel_async_inner_iter<<<dimGrid, dimBlock, 0, stream>>>
  (A_diag, row_diag, col_diag, nrows, row_offset,
   A_offdiag, row_offdiag, col_offdiag,
   inv_D, b, x, xg, n_iter);
#endif
}

template <typename DataType>
void cuda_async_inner_iter_no_offdiag(const DataType *A_diag,
                      const int *row_diag,
                      const int *col_diag,
                      const int nrows,
                      const int row_offset,
                      const DataType *inv_D,
                      const DataType *b,
                      DataType *x,
                      const int n_iter,
                      const int dim_grid,
                      const int dim_block,
                      cudaStream_t stream)
{
#ifndef NOCUDA
  dim3 dimGrid(dim_grid);
  dim3 dimBlock(dim_block);

  kernel_async_inner_iter_no_offdiag<<<dimGrid, dimBlock, 0, stream>>>
  (A_diag, row_diag, col_diag, nrows, row_offset,
   inv_D, b, x, n_iter);
#endif
}

template <typename DataType>
void cuda_async_inner_damped_iter(const DataType *A_diag,
                      const int *row_diag,
                      const int *col_diag,
                      const int nrows,
                      const int row_offset,
                      const DataType *A_offdiag,
                      const int *row_offdiag,
                      const int *col_offdiag,
                      const DataType *inv_D,
                      const DataType *b,
                      DataType *x,
                      const DataType *xg,
                      const DataType w,
                      const int n_iter,
                      const int dim_grid,
                      const int dim_block,
                      cudaStream_t stream)
{
#ifndef NOCUDA
  dim3 dimGrid(dim_grid);
  dim3 dimBlock(dim_block);

  kernel_async_inner_damped_iter<<<dimGrid, dimBlock, 0, stream>>>
  (A_diag, row_diag, col_diag, nrows, row_offset,
   A_offdiag, row_offdiag, col_offdiag,
   inv_D, b, x, xg, w, n_iter);
#endif
}

template <typename DataType>
void cuda_async_inner_damped_iter_no_offdiag
     (const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
      const DataType *inv_D, const DataType *b, DataType *x,
      const DataType w, const int n_iter,
      const int dim_grid, const int dim_block, cudaStream_t stream)
{
#ifndef NOCUDA
  dim3 dimGrid(dim_grid);
  dim3 dimBlock(dim_block);

  kernel_async_inner_damped_iter_no_offdiag<<<dimGrid, dimBlock, 0, stream>>>
  (A_diag, row_diag, col_diag, nrows, row_offset,
   inv_D, b, x, w, n_iter);  
#endif
}

template <typename DataType>
void cuda_compute_block_squared_residual
                     (const DataType *A_diag,
                      const int *row_diag,
                      const int *col_diag,
                      const int nrows,
                      const int row_offset,
                      const DataType *A_offdiag,
                      const int *row_offdiag,
                      const int *col_offdiag,
                      const DataType *D,
                      const DataType *b,
                      const DataType *x,
                      const DataType *xg,
                      DataType *block_res,
                      const int dim_grid,
                      const int dim_block,
                      const int max_dim_block,
                      cudaStream_t stream)
{
#ifndef NOCUDA
  dim3 dimGrid(dim_grid);
  dim3 dimBlock(dim_block);

  kernel_compute_block_squared_residual<<<dimGrid, dimBlock, 0, stream>>>
  (A_diag, row_diag, col_diag, nrows, row_offset,
   A_offdiag, row_offdiag, col_offdiag,
   D, b, x, xg, block_res);

  switch (max_dim_block)
  {
    case 4096:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 4096 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 2048:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 2048 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 1024: // Tesla K40
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 1024 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 512:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 512 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 256:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 256 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 128:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 128 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 64:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 64 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 32:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 32 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 16:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 16 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 8:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 8 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 4:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 4 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 2:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 2 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 1:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 1 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    default:
      LOG_ERROR("Thread block size is not a power of 2.");
      exit(-1);
      break;
  }
#endif
}

template <typename DataType>
void cuda_compute_block_squared_residual_no_offdiag
                     (const DataType *A_diag,
                      const int *row_diag,
                      const int *col_diag,
                      const int nrows,
                      const int row_offset,
                      const DataType *D,
                      const DataType *b,
                      const DataType *x,
                      DataType *block_res,
                      const int dim_grid,
                      const int dim_block,
                      const int max_dim_block,
                      cudaStream_t stream)
{
#ifndef NOCUDA
  dim3 dimGrid(dim_grid);
  dim3 dimBlock(dim_block);

  kernel_compute_block_squared_residual_no_offdiag<<<dimGrid, dimBlock, 0, stream>>>
  (A_diag, row_diag, col_diag, nrows, row_offset,
   D, b, x, block_res);

  switch (max_dim_block)
  {
    case 4096:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 4096 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 2048:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 2048 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 1024: // Tesla K40
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 1024 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 512:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 512 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 256:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 256 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 128:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 128 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 64:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 64 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 32:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 32 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 16:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 16 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 8:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 8 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 4:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 4 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 2:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 2 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    case 1:
      kernel_reduce_block_squared_residual
      <<<1, max_dim_block, 1 * sizeof(DataType), stream>>>
      (block_res, nrows);
      break;

    default:
      LOG_ERROR("Thread block size for reduction kernel is not a power of 2.");
      exit(-1);
      break;
  }
#endif
}

// explicit template instantiation
template void cuda_async_inner_iter<double>(const double*,
                                            const int*,
                                            const int*,
                                            const int,
                                            const int,
                                            const double*,
                                            const int*,
                                            const int*,
                                            const double*,
                                            const double*,
                                            double*,
                                            const double*,
                                            const int,
                                            const int,
                                            const int,
                                            cudaStream_t);

template void cuda_async_inner_iter<float>(const float*,
                                           const int*,
                                           const int*,
                                           const int,
                                           const int,
                                           const float*,
                                           const int*,
                                           const int*,
                                           const float*,
                                           const float*,
                                           float*,
                                           const float*,
                                           const int,
                                           const int,
                                           const int,
                                           cudaStream_t);

template void cuda_async_inner_iter_no_offdiag<double>(const double*,
                                            const int*,
                                            const int*,
                                            const int,
                                            const int,
                                            const double*,
                                            const double*,
                                            double*,
                                            const int,
                                            const int,
                                            const int,
                                            cudaStream_t);

template void cuda_async_inner_iter_no_offdiag<float>(const float*,
                                           const int*,
                                           const int*,
                                           const int,
                                           const int,
                                           const float*,
                                           const float*,
                                           float*,
                                           const int,
                                           const int,
                                           const int,
                                           cudaStream_t);

template void cuda_async_inner_damped_iter<double>(const double*,
                                            const int*,
                                            const int*,
                                            const int,
                                            const int,
                                            const double*,
                                            const int*,
                                            const int*,
                                            const double*,
                                            const double*,
                                            double*,
                                            const double*,
                                            const double,
                                            const int,
                                            const int,
                                            const int,
                                            cudaStream_t);

template void cuda_async_inner_damped_iter<float>(const float*,
                                           const int*,
                                           const int*,
                                           const int,
                                           const int,
                                           const float*,
                                           const int*,
                                           const int*,
                                           const float*,
                                           const float*,
                                           float*,
                                           const float*,
                                           const float,
                                           const int,
                                           const int,
                                           const int,
                                           cudaStream_t);

template void cuda_async_inner_damped_iter_no_offdiag<double>(const double*,
                                            const int*,
                                            const int*,
                                            const int,
                                            const int,
                                            const double*,
                                            const double*,
                                            double*,
                                            const double,
                                            const int,
                                            const int,
                                            const int,
                                            cudaStream_t);

template void cuda_async_inner_damped_iter_no_offdiag<float>(const float*,
                                           const int*,
                                           const int*,
                                           const int,
                                           const int,
                                           const float*,
                                           const float*,
                                           float*,
                                           const float,
                                           const int,
                                           const int,
                                           const int,
                                           cudaStream_t);

template void cuda_compute_block_squared_residual<double>(const double*,
                                            const int*,
                                            const int*,
                                            const int,
                                            const int,
                                            const double*,
                                            const int*,
                                            const int*,
                                            const double*,
                                            const double*,
                                            const double*,
                                            const double*,
                                            double*,
                                            const int,
                                            const int,
                                            const int,
                                            cudaStream_t);

template void cuda_compute_block_squared_residual<float>(const float*,
                                           const int*,
                                           const int*,
                                           const int,
                                           const int,
                                           const float*,
                                           const int*,
                                           const int*,
                                           const float*,
                                           const float*,
                                           const float*,
                                           const float*,
                                           float*,
                                           const int,
                                           const int,
                                           const int,
                                           cudaStream_t);

template void cuda_compute_block_squared_residual_no_offdiag<double>(const double*,
                                            const int*,
                                            const int*,
                                            const int,
                                            const int,
                                            const double*,
                                            const double*,
                                            const double*,
                                            double*,
                                            const int,
                                            const int,
                                            const int,
                                            cudaStream_t
                                            );

template void cuda_compute_block_squared_residual_no_offdiag<float>(const float*,
                                           const int*,
                                           const int*,
                                           const int,
                                           const int,
                                           const float*,
                                           const float*,
                                           const float*,
                                           float*,
                                           const int,
                                           const int,
                                           const int,
                                           cudaStream_t
                                           );
