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

#include "cuda_gmg.h"
#include "cuda_gmg_kernel.h"

template <typename DataType>
void cuda_gmg_compute_residual
     (const DataType *A_diag,    const int *col_A_diag,    const int *row_A_diag,    const int nrows,
      const DataType *A_offdiag, const int *col_A_offdiag, const int *row_A_offdiag,
      const DataType *sol, const DataType *rhs, DataType *res,
      const DataType *sol_ghost,
      const int dim_grid, const int dim_block, cudaStream_t stream)
{
  dim3 dimGrid (dim_grid);
  dim3 dimBlock(dim_block);

  kernel_gmg_compute_residual<<<dimGrid, dimBlock, 0, stream>>>
     (A_diag, col_A_diag, row_A_diag, nrows,
      A_offdiag, col_A_offdiag, row_A_offdiag,
      sol, rhs, res, sol_ghost);
}

template <typename DataType>
void cuda_gmg_compute_residual_no_offdiag
     (const DataType *A_diag, const int *col_A_diag, const int *row_A_diag, const int nrows,
      const DataType *sol, const DataType *rhs, DataType *res,
      const int dim_grid, const int dim_block, cudaStream_t stream)
{
  dim3 dimGrid (dim_grid);
  dim3 dimBlock(dim_block);

  kernel_gmg_compute_residual_no_offdiag<<<dimGrid, dimBlock, 0, stream>>>
     (A_diag, col_A_diag, row_A_diag, nrows,
      sol, rhs, res);
}

template <typename DataType>
void cuda_gmg_compute_restriction
     (const DataType *R_diag,    const int *col_R_diag,    const int *row_R_diag,    const int nrows,
      const DataType *R_offdiag, const int *col_R_offdiag, const int *row_R_offdiag,
      DataType *res, const DataType *res_ghost, const DataType *tmp,
      const int dim_grid, const int dim_block, cudaStream_t stream)
{
  dim3 dimGrid (dim_grid);
  dim3 dimBlock(dim_block);

  kernel_gmg_compute_restriction<<<dimGrid, dimBlock, 0, stream>>>
  (R_diag,    col_R_diag,    row_R_diag,    nrows,
   R_offdiag, col_R_offdiag, row_R_offdiag,
   res, res_ghost, tmp);
}

template <typename DataType>
void cuda_gmg_compute_restriction_no_offdiag
     (const DataType *R_diag, const int *col_R_diag, const int *row_R_diag, const int nrows,
      DataType *res, const DataType *tmp,
      const int dim_grid, const int dim_block, cudaStream_t stream)
{
  dim3 dimGrid (dim_grid);
  dim3 dimBlock(dim_block);

  kernel_gmg_compute_restriction_no_offdiag<<<dimGrid, dimBlock, 0, stream>>>
  (R_diag, col_R_diag, row_R_diag, nrows,
   res, tmp);
}

template <typename DataType>
void cuda_gmg_compute_updated_solution
     (const DataType *P_diag, const int *col_P_diag, const int *row_P_diag, const int nrows_P_diag,
      const DataType *P_offdiag, const int *col_P_offdiag, const int *row_P_offdiag,
      DataType *sol, const DataType *res, const DataType *res_ghost,
      const int dim_grid, const int dim_block, cudaStream_t stream)
{
  dim3 dimGrid (dim_grid);
  dim3 dimBlock(dim_block);

  kernel_gmg_compute_updated_solution<<<dimGrid, dimBlock, 0, stream>>>
     (P_diag, col_P_diag, row_P_diag, nrows_P_diag,
      P_offdiag, col_P_offdiag, row_P_offdiag,
      sol, res, res_ghost);
}

template <typename DataType>
void cuda_gmg_compute_updated_solution_no_offdiag
     (const DataType *P_diag, const int *col_P_diag, const int *row_P_diag, const int nrows_P_diag,
      DataType *sol, const DataType *res,
      const int dim_grid, const int dim_block, cudaStream_t stream)
{
  dim3 dimGrid (dim_grid);
  dim3 dimBlock(dim_block);

  kernel_gmg_compute_updated_solution_no_offdiag<<<dimGrid, dimBlock, 0, stream>>>
     (P_diag, col_P_diag, row_P_diag, nrows_P_diag,
      sol, res);
}

// explicit template instantiations

template void cuda_gmg_compute_residual<double>
     (const double*, const int*, const int*, const int,
      const double*, const int*, const int*,
      const double*, const double*, double*,
      const double*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_residual<float>
     (const float*, const int*, const int*, const int,
      const float*, const int*, const int*,
      const float*, const float*, float*,
      const float*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_residual_no_offdiag<double>
     (const double*, const int*, const int*, const int,
      const double*, const double*, double*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_residual_no_offdiag<float>
     (const float*, const int*, const int*, const int,
      const float*, const float*, float*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_restriction<double>
     (const double*, const int*, const int*, const int,
      const double*, const int*, const int*,
      double*, const double*, const double*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_restriction<float>
     (const float*, const int*, const int*, const int,
      const float*, const int*, const int*,
      float*, const float*, const float*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_restriction_no_offdiag<double>
     (const double*, const int*, const int*, const int,
      double*, const double*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_restriction_no_offdiag<float>
     (const float*, const int*, const int*, const int,
      float*, const float*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_updated_solution<double>
     (const double*, const int*, const int*, const int,
      const double*, const int*, const int*,
      double*, const double*, const double*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_updated_solution<float>
     (const float*, const int*, const int*, const int,
      const float*, const int*, const int*,
      float*, const float*, const float*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_updated_solution_no_offdiag<double>
     (const double*, const int*, const int*, const int,
      double*, const double*,
      const int, const int, cudaStream_t);

template void cuda_gmg_compute_updated_solution_no_offdiag<float>
     (const float*, const int*, const int*, const int,
      float*, const float*,
      const int, const int, cudaStream_t);