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
__global__ void kernel_async_inner_iter
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *A_offdiag, const int *row_offdiag, const int *col_offdiag,
  const DataType *inv_D, const DataType *b, DataType *x, const DataType *xg, const int n_iter )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idx_offset = idx + row_offset;
    int j, n;

    DataType c, tmp;

    if ( idx < nrows )
    {
        // precompute c = b - A_offdiag x_ghost
        c = b[idx];
        for ( j = row_offdiag[idx]; j < row_offdiag[idx + 1]; ++j )
        {
            c -= A_offdiag[j] * xg[col_offdiag[j]];
        }

        // inner iterations
        for ( n = 0; n < n_iter; ++n )
        {
            // tmp = b - A_offdiag x_ghost
            tmp = c;

            // tmp = b - Ax
            for ( j = row_diag[idx]; j < row_diag[idx + 1]; ++j )
            {
                tmp -= A_diag[j] * x[col_diag[j]];
            }

            // x = D^-1(b - Ax)
            x[idx_offset] = inv_D[idx] * tmp;
        }
    }
}

template <typename DataType>
__global__ void kernel_async_inner_iter_no_offdiag
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *inv_D, const DataType *b, DataType *x, const int n_iter )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idx_offset = idx + row_offset;
    int j, n;

    DataType tmp;

    if ( idx < nrows )
    {
        // inner iterations
        for ( n = 0; n < n_iter; ++n )
        {
            // tmp = b
            tmp = b[idx];

            // tmp = b - Ax
            for ( j = row_diag[idx]; j < row_diag[idx + 1]; ++j )
            {
                tmp -= A_diag[j] * x[col_diag[j]];
            }

            // x = D^-1(b - Ax)
            x[idx_offset] = inv_D[idx] * tmp;
        }
    }
}

template <typename DataType>
__global__ void kernel_async_inner_damped_iter
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *A_offdiag, const int *row_offdiag, const int *col_offdiag,
  const DataType *inv_D, const DataType *b, DataType *x, const DataType *xg,
  const DataType w, const int n_iter )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idx_offset = idx + row_offset;
    int j, n;

    DataType one_minus_w = ( DataType ) ( 1.0 ) - w;
    DataType c, tmp;

    if ( idx < nrows )
    {
        // precompute c = b - A_offdiag x_ghost
        c = b[idx];
        for ( j = row_offdiag[idx]; j < row_offdiag[idx + 1]; ++j )
        {
            c -= A_offdiag[j] * xg[col_offdiag[j]];
        }

        // inner iterations
        for ( n = 0; n < n_iter; ++n )
        {
            // tmp = b - A_offdiag x_ghost
            tmp = c;

            // tmp = b - (L+U)x
            for ( j = row_diag[idx]; j < row_diag[idx + 1]; ++j )
            {
                tmp -= A_diag[j] * x[col_diag[j]];
            }

            // x = w D^-1 [b - (L+U)x] + (1-w)x
            x[idx_offset] = ( w * inv_D[idx] * tmp ) + ( one_minus_w * x[idx_offset] );
        }
    }
}

template <typename DataType>
__global__ void kernel_async_inner_damped_iter_no_offdiag
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *inv_D, const DataType *b, DataType *x,
  const DataType w, const int n_iter )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idx_offset = idx + row_offset;
    int j, n;

    DataType one_minus_w = ( DataType ) ( 1.0 ) - w;
    DataType tmp;

    if ( idx < nrows )
    {

        // inner iterations
        for ( n = 0; n < n_iter; ++n )
        {
            // tmp = b
            tmp = b[idx];

            // tmp = b - (L+U)x
            for ( j = row_diag[idx]; j < row_diag[idx + 1]; ++j )
            {
                tmp -= A_diag[j] * x[col_diag[j]];
            }

            // x = w D^-1 [b - (L+U)x] + (1-w)x
            x[idx_offset] = ( w * inv_D[idx] * tmp ) + ( one_minus_w * x[idx_offset] );
        }
    }
}

template <typename DataType>
__global__ void kernel_compute_block_squared_residual
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *A_offdiag, const int *row_offdiag, const int *col_offdiag,
  const DataType *D, const DataType *b, const DataType *x, const DataType *xg, DataType *block_res )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idx_offset = idx + row_offset;
    int j;

    if ( idx < nrows )
    {
        // compute c = b - Ax
        block_res[idx] = b[idx];

        for ( j = row_diag[idx]; j < row_diag[idx + 1]; ++j )
        {
            block_res[idx] -= A_diag[j] * x[col_diag[j]];
        }

        for ( j = row_offdiag[idx]; j < row_offdiag[idx + 1]; ++j )
        {
            block_res[idx] -= A_offdiag[j] * xg[col_offdiag[j]];
        }

        block_res[idx] -= D[idx] * x[idx_offset];

        // square res
        block_res[idx] *= block_res[idx];
    }
}

template <typename DataType>
__global__ void kernel_compute_block_squared_residual_no_offdiag
( const DataType *A_diag, const int *row_diag, const int *col_diag, const int nrows, const int row_offset,
  const DataType *D, const DataType *b, const DataType *x, DataType *block_res )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idx_offset = idx + row_offset;
    int j;

    if ( idx < nrows )
    {
        // compute res = b - Ax
        block_res[idx] = b[idx];

        for ( j = row_diag[idx]; j < row_diag[idx + 1]; ++j )
        {
            block_res[idx] -= A_diag[j] * x[col_diag[j]];
        }

        block_res[idx] -= D[idx] * x[idx_offset];

        // square res
        block_res[idx] *= block_res[idx];
    }
}

__device__ void warpReduce ( volatile double *red, int idx )
{
    if ( blockDim.x >= 64 ) red[idx] += red[idx + 32];
    if ( blockDim.x >= 32 ) red[idx] += red[idx + 16];
    if ( blockDim.x >= 16 ) red[idx] += red[idx + 8];
    if ( blockDim.x >= 8 ) red[idx] += red[idx + 4];
    if ( blockDim.x >= 4 ) red[idx] += red[idx + 2];
    if ( blockDim.x >= 2 ) red[idx] += red[idx + 1];
}

__device__ void warpReduce ( volatile float *red, int idx )
{
    if ( blockDim.x >= 64 ) red[idx] += red[idx + 32];
    if ( blockDim.x >= 32 ) red[idx] += red[idx + 16];
    if ( blockDim.x >= 16 ) red[idx] += red[idx + 8];
    if ( blockDim.x >= 8 ) red[idx] += red[idx + 4];
    if ( blockDim.x >= 4 ) red[idx] += red[idx + 2];
    if ( blockDim.x >= 2 ) red[idx] += red[idx + 1];
}

__global__ void kernel_reduce_block_squared_residual ( double *block_res, const int n )
{
    extern __shared__ double reducd[];
    int idx = threadIdx.x;
    int j = idx;

    reducd[idx] = 0.0;

    // accumulate whole array in shared memory of this block
    while ( j < n )
    {
        reducd[idx] += block_res[j];
        j += blockDim.x;
    }
    __syncthreads ( );

    // rely on block size equal to a power of 2
    if ( blockDim.x >= 4096 )
    {
        if ( idx < 2048 ) reducd[idx] += reducd[idx + 2048];
        __syncthreads ( );
    }
    if ( blockDim.x >= 2048 )
    {
        if ( idx < 1024 ) reducd[idx] += reducd[idx + 1024];
        __syncthreads ( );
    }
    if ( blockDim.x >= 1024 )
    {
        if ( idx < 512 ) reducd[idx] += reducd[idx + 512];
        __syncthreads ( );
    }
    if ( blockDim.x >= 512 )
    {
        if ( idx < 256 ) reducd[idx] += reducd[idx + 256];
        __syncthreads ( );
    }
    if ( blockDim.x >= 256 )
    {
        if ( idx < 128 ) reducd[idx] += reducd[idx + 128];
        __syncthreads ( );
    }
    if ( blockDim.x >= 128 )
    {
        if ( idx < 64 ) reducd[idx] += reducd[idx + 64];
        __syncthreads ( );
    }

    if ( idx < 32 ) warpReduce ( reducd, idx );
    if ( idx == 0 ) block_res[blockIdx.x] = reducd[0];
}

__global__ void kernel_reduce_block_squared_residual ( float *block_res, const int n )
{
    extern __shared__ float reducf[];
    int idx = threadIdx.x;
    int j = idx;

    reducf[idx] = 0.0;

    // accumulate whole array in shared memory of this block
    while ( j < n )
    {
        reducf[idx] += block_res[j];
        j += blockDim.x;
    }
    __syncthreads ( );

    // rely on block size equal to a power of 2
    if ( blockDim.x >= 4096 )
    {
        if ( idx < 2048 ) reducf[idx] += reducf[idx + 2048];
        __syncthreads ( );
    }
    if ( blockDim.x >= 2048 )
    {
        if ( idx < 1024 ) reducf[idx] += reducf[idx + 1024];
        __syncthreads ( );
    }
    if ( blockDim.x >= 1024 )
    {
        if ( idx < 512 ) reducf[idx] += reducf[idx + 512];
        __syncthreads ( );
    }
    if ( blockDim.x >= 512 )
    {
        if ( idx < 256 ) reducf[idx] += reducf[idx + 256];
        __syncthreads ( );
    }
    if ( blockDim.x >= 256 )
    {
        if ( idx < 128 ) reducf[idx] += reducf[idx + 128];
        __syncthreads ( );
    }
    if ( blockDim.x >= 128 )
    {
        if ( idx < 64 ) reducf[idx] += reducf[idx + 64];
        __syncthreads ( );
    }

    if ( idx < 32 ) warpReduce ( reducf, idx );
    if ( idx == 0 ) block_res[blockIdx.x] = reducf[0];
}

