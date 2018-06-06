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

template <typename DataType>
__global__ void kernel_gmg_compute_residual
( const DataType *A_diag, const int *col_A_diag, const int *row_A_diag, const int nrows,
  const DataType *A_offdiag, const int *col_A_offdiag, const int *row_A_offdiag,
  const DataType *sol, const DataType *rhs, DataType *res,
  const DataType *sol_ghost )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j;
    DataType result;

    if ( idx < nrows )
    {
        // inititalise result with 0
        result = rhs[idx];

        // compute A_diag * rhs
        for ( j = row_A_diag[idx]; j != row_A_diag[idx + 1]; ++j )
        {
            result -= A_diag[j] * sol[col_A_diag[j]];
        }

        // add A_offdiag * rhs_ghost
        for ( j = row_A_offdiag[idx]; j != row_A_offdiag[idx + 1]; ++j )
        {
            result -= A_offdiag[j] * sol_ghost[col_A_offdiag[j]];
        }

        res[idx] = result;
    }
}

template <typename DataType>
__global__ void kernel_gmg_compute_residual_no_offdiag
( const DataType *A_diag, const int *col_A_diag, const int *row_A_diag, const int nrows,
  const DataType *sol, const DataType *rhs, DataType *res )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j;
    DataType result;

    if ( idx < nrows )
    {
        // inititalise result with 0
        result = rhs[idx];

        // compute A_diag * rhs
        for ( j = row_A_diag[idx]; j != row_A_diag[idx + 1]; ++j )
        {
            result -= A_diag[j] * sol[col_A_diag[j]];
        }

        res[idx] = result;
    }
}

template <typename DataType>
__global__ void kernel_gmg_compute_restriction
( const DataType *R_diag, const int *col_R_diag, const int *row_R_diag, const int nrows,
  const DataType *R_offdiag, const int *col_R_offdiag, const int *row_R_offdiag,
  DataType *res, const DataType *res_ghost, const DataType *tmp )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j;
    DataType result;

    if ( idx < nrows )
    {
        // inititalise result with 0
        result = 0.0;

        // compute R_diag * tmp
        for ( j = row_R_diag[idx]; j != row_R_diag[idx + 1]; ++j )
        {
            result += R_diag[j] * tmp[col_R_diag[j]];
        }

        // add R_offdiag * rhs_ghost
        for ( j = row_R_offdiag[idx]; j != row_R_offdiag[idx + 1]; ++j )
        {
            result += R_offdiag[j] * res_ghost[col_R_offdiag[j]];
        }

        res[idx] = result;
    }
}

template <typename DataType>
__global__ void kernel_gmg_compute_restriction_no_offdiag
( const DataType *R_diag, const int *col_R_diag, const int *row_R_diag, const int nrows,
  DataType *res, const DataType *tmp )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j;
    DataType result;

    if ( idx < nrows )
    {
        // inititalise result with 0
        result = 0.0;

        // compute R_diag * tmp
        for ( j = row_R_diag[idx]; j != row_R_diag[idx + 1]; ++j )
        {
            result += R_diag[j] * tmp[col_R_diag[j]];
        }

        res[idx] = result;
    }
}

template <typename DataType>
__global__ void kernel_gmg_compute_updated_solution
( const DataType *P_diag, const int *col_P_diag, const int *row_P_diag, const int nrows_P_diag,
  const DataType *P_offdiag, const int *col_P_offdiag, const int *row_P_offdiag,
  DataType *sol, const DataType *res, const DataType *res_ghost )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j;
    DataType result;

    if ( idx < nrows_P_diag )
    {
        // inititalise result
        result = sol[idx];

        // add P_diag * res
        for ( j = row_P_diag[idx]; j != row_P_diag[idx + 1]; ++j )
        {
            result += P_diag[j] * res[col_P_diag[j]];
        }

        // add P_offdiag * res_ghost
        for ( j = row_P_offdiag[idx]; j != row_P_offdiag[idx + 1]; ++j )
        {
            result += P_offdiag[j] * res_ghost[col_P_offdiag[j]];
        }

        sol[idx] = result;
    }
}

template <typename DataType>
__global__ void kernel_gmg_compute_updated_solution_no_offdiag
( const DataType *P_diag, const int *col_P_diag, const int *row_P_diag, const int nrows_P_diag,
  DataType *sol, const DataType *res )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j;
    DataType result;

    if ( idx < nrows_P_diag )
    {
        // inititalise result
        result = sol[idx];

        // add P_diag * res
        for ( j = row_P_diag[idx]; j != row_P_diag[idx + 1]; ++j )
        {
            result += P_diag[j] * res[col_P_diag[j]];
        }

        sol[idx] = result;
    }
}