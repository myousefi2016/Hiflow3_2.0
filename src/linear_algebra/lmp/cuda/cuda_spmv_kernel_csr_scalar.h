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

#ifndef __CUDA_SPMV_KERNEL_CSR_SCALAR_H
#    define __CUDA_SPMV_KERNEL_CSR_SCALAR_H

template <typename ValueType>
__global__ void kernel_spmv_csr_scalar ( ValueType *val, int *col, int *row, ValueType *x, ValueType *y, int n )
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int j;

    if ( index < n )
    {
        y[index] = 0;
        for ( j = row[index]; j < row[index + 1]; j++ )
        {
            y[index] += val[j] * x[col[j]];
        }
    }
}

template <typename ValueType>
__global__ void kernel_add_spmv_csr_scalar ( ValueType *val, int *col, int *row, ValueType *x, ValueType *y, int n )
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int j;

    if ( index < n )
    {
        for ( j = row[index]; j < row[index + 1]; j++ )
        {
            y[index] += val[j] * x[col[j]];
        }
    }
}

template <typename ValueType>
__global__ void kernel_spmvnd_csr_scalar ( const ValueType *val, const int *col, const int *row,
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

template <typename ValueType>
__device__ void gpu_BlockPsgauss_seidel ( const ValueType *invec, ValueType *outvec,
                                          const int *rows, const int *cols, const ValueType *vals,
                                          const int start_i, const int end_i )
{

    // M = (D+L) D^-1 (D+R)
    // Mz=r

    // (D+L) y = r
    // forward step

    int last_j = 0;
    int i = 0;
    for ( i = start_i; i < end_i; ++i )
    {

        outvec[i] = invec[i];

        for ( int j = rows[i]; j < rows[i + 1]; ++j )
        {

            if ( ( cols[j] >= start_i ) &&
                 ( cols[j] < end_i ) )
            {

                if ( i > cols[j] )
                    outvec[i] -= outvec[ cols[j] ] * vals[j];
                if ( i == cols[j] )
                {
                    last_j = j;
                    break;
                }
            }

        }

        outvec[i] /= vals[last_j];

    }

    // Dy
    for ( int i = start_i; i < end_i; ++i )
        for ( int j = rows[i]; j < rows[i + 1]; ++j )
            if ( i == cols[j] )
            {
                outvec[i] *= vals[j];
                break;
            }

    // backward
    // (D+R)r = Dy
    for ( int i = end_i - 1; i >= start_i; --i )
    {

        for ( int j = rows[i]; j < rows[i + 1]; ++j )
        {

            if ( ( cols[j] >= start_i ) &&
                 ( cols[j] < end_i ) )
            {

                if ( i == cols[j] )
                    last_j = j;
                if ( i < cols[j] )
                    outvec[i] -= outvec[ cols[j] ] * vals[j];
            }

        }
        outvec[i] /= vals[last_j];
    }

}

template <typename ValueType>
__global__ void kernel_BlockPsgauss_seidel ( const ValueType *val, const int *col, const int *row,
                                             const ValueType *invec, ValueType *outvec,
                                             const int max_partitions, const int size_max,
                                             const int step_size )
{
    int np = blockIdx.x * blockDim.x + threadIdx.x;

    if ( np < max_partitions )
    {

        int start_i = np*step_size;
        int end_i = ( np + 1 ) * step_size;

        if ( np == max_partitions - 1 ) end_i = size_max;

        gpu_BlockPsgauss_seidel ( invec, outvec,
                                  row, col, val,
                                  start_i, end_i );
    }
}

#endif
