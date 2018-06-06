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

/// @author Dimitar Lukarski

#ifndef __CUDA_SPMV_KERNEL_CSR_SCALAR_TEX_H
#    define __CUDA_SPMV_KERNEL_CSR_SCALAR_TEX_H

// extracted from math_functions.h because #include <math_functions.h>
//   caused compile errors
extern __device__ double __hiloint2double ( int, int );

texture<int2, 1> tex_double;
texture<float, 1> tex_float;

void bind2tex ( const double *x )
{
    cudaBindTexture ( NULL, tex_double, x );
}

void unbind2tex ( const double *x )
{
    cudaUnbindTexture ( tex_double );
}

void bind2tex ( const float *x )
{
    cudaBindTexture ( NULL, tex_float, x );
}

void unbind2tex ( const float *x )
{
    cudaUnbindTexture ( tex_float );
}

__inline__ __device__ double read_from_tex ( const int& i, const double *x )
{
    int2 temp = tex1Dfetch ( tex_double, i );
    return __hiloint2double ( temp.y, temp.x );
}

__inline__ __device__ float read_from_tex ( const int& i, const float *x )
{
    return tex1Dfetch ( tex_float, i );
}

template <typename ValueType>
__global__ void kernel_spmv_csr_scalar_tex ( ValueType *val, int *col, int *row, ValueType *x, ValueType *y, int n )
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int j;

    if ( index < n )
    {
        y[index] = 0;
        for ( j = row[index]; j < row[index + 1]; j++ )
        {
            y[index] += val[j] * read_from_tex ( col[j], x );
        }
    }
}

template <typename ValueType>
__global__ void kernel_add_spmv_csr_scalar_tex ( ValueType *val, int *col, int *row, ValueType *x, ValueType *y, int n )
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int j;

    if ( index < n )
    {
        for ( j = row[index]; j < row[index + 1]; j++ )
        {
            y[index] += val[j] * read_from_tex ( col[j], x );
        }
    }

}

template <typename ValueType>
__global__ void kernel_spmvnd_csr_scalar_tex ( const ValueType *val, const int *col, const int *row,
                                               const ValueType *in, ValueType *out, const int n_rows )
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int j;

    if ( idx < n_rows )
    {
        out[idx] = 0.0;
        for ( j = row[idx]; j < row[idx + 1]; ++j )
        {
            if ( idx != col[j] )
            {
                out[idx] += val[j] * in[col[j]];
            }
        }
    }
}

#endif
