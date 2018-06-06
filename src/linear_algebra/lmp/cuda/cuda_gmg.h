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
void cuda_gmg_compute_residual
( const DataType *A_diag, const int *col_A_diag, const int *row_A_diag, const int nrows,
  const DataType *A_offdiag, const int *col_A_offdiag, const int *row_A_offdiag,
  const DataType *sol, const DataType *rhs, DataType *res,
  const DataType *sol_ghost,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_gmg_compute_residual_no_offdiag
( const DataType *A_diag, const int *col_A_diag, const int *row_A_diag, const int nrows,
  const DataType *sol, const DataType *rhs, DataType *res,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_gmg_compute_restriction
( const DataType *R_diag, const int *col_R_diag, const int *row_R_diag, const int nrows,
  const DataType *R_offdiag, const int *col_R_offdiag, const int *row_R_offdiag,
  DataType *res, const DataType *res_ghost, const DataType *tmp,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_gmg_compute_restriction_no_offdiag
( const DataType *R_diag, const int *col_R_diag, const int *row_R_diag, const int nrows,
  DataType *res, const DataType *tmp,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_gmg_compute_updated_solution
( const DataType *P_diag, const int *col_P_diag, const int *row_P_diag, const int nrows_P_diag,
  const DataType *P_offdiag, const int *col_P_offdiag, const int *row_P_offdiag,
  DataType *sol, const DataType *res, const DataType *res_ghost,
  const int dim_grid, const int dim_block, cudaStream_t stream );

template <typename DataType>
void cuda_gmg_compute_updated_solution_no_offdiag
( const DataType *P_diag, const int *col_P_diag, const int *row_P_diag, const int nrows_P_diag,
  DataType *sol, const DataType *res,
  const int dim_grid, const int dim_block, cudaStream_t stream );