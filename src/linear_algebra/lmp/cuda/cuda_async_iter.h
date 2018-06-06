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

#include <cuda_runtime_api.h>

template <typename DataType>
void cuda_async_inner_iter
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *A_offdiag, const int *row_offdiag, const int *col_offdiag,
  const DataType *inv_D, const DataType *b, DataType *x, const DataType *xg, const int n_iter,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_async_inner_iter_no_offdiag
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *inv_D, const DataType *b, DataType *x, const int n_iter,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_async_inner_damped_iter
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *A_offdiag, const int *row_offdiag, const int *col_offdiag,
  const DataType *inv_D, const DataType *b, DataType *x, const DataType *xg,
  const DataType w, const int n_iter,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_async_inner_damped_iter_no_offdiag
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *inv_D, const DataType *b, DataType *x,
  const DataType w, const int n_iter,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_compute_block_squared_residual
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *A_offdiag, const int *row_offdiag, const int *col_offdiag,
  const DataType *D, const DataType *b, const DataType *x, const DataType *xg, DataType *block_res,
  const int dim_grid, const int dim_block, const int max_dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_compute_block_squared_residual_no_offdiag
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *D, const DataType *b, const DataType *x, DataType *block_res,
  const int dim_grid, const int dim_block, const int max_dim_block, cudaStream_t stream );