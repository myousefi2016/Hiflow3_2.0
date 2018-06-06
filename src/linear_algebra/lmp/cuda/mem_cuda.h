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

#ifndef __MEM_CUDA_H
#    define __MEM_CUDA_H

template <typename ValueType>
__global__ void kernel_cudasetvalues ( const int *index, const size_t size, const ValueType *values, ValueType *buffer )
{

    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        buffer[ index[ind] ] = values[ ind ];
    }

}

template <typename ValueType>
__global__ void kernel_cudagetvalues ( const int *index, const size_t size, ValueType *values, const ValueType *buffer )
{

    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        values[ ind ] = buffer[ index[ind] ];
    }

}

template <typename ValueType>
__global__ void kernel_cudasetblockvalues ( const int start_i, const int start_sub_vec, const int size, const ValueType *values, ValueType *buffer )
{

    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        buffer[ ind + start_i ] = values[ ind + start_sub_vec ];
    }

}

template <typename ValueType>
__global__ void kernel_cudaaddblockvalues ( const int start_i, const int start_sub_vec, const int size, const ValueType *values, ValueType *buffer, ValueType weight )
{

    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        buffer[ ind + start_i ] += weight * values[ ind + start_sub_vec ];
    }

}

template <typename ValueType>
__global__ void kernel_cudagetblockvalues ( const int start_i, const int end_i, ValueType *values, const ValueType *buffer )
{

    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < end_i - start_i )
    {
        values[ ind ] = buffer[ ind + start_i ];
    }

}

template <typename ValueType>
__global__ void kernel_cudamultvalues ( const int size, const ValueType *values, ValueType *buffer )
{

    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        buffer[ind] *= values[ ind ];
    }

}

template <typename ValueType>
__global__ void kernel_cudacastfromdouble ( const int size, ValueType *dest, const double *src )
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        dest[ind] = ( ValueType ) ( src[ind] );
    }
}

template <typename ValueType>
__global__ void kernel_cudacastfromfloat ( const int size, ValueType *dest, const float *src )
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        dest[ind] = ( ValueType ) ( src[ind] );
    }
}

template <typename ValueType>
__global__ void kernel_cudacasttodouble ( const int size, double *dest, const ValueType *src )
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        dest[ind] = ( double ) ( src[ind] );
    }
}

template <typename ValueType>
__global__ void kernel_cudacasttofloat ( const int size, float *dest, const ValueType *src )
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;

    if ( ind < size )
    {
        dest[ind] = ( float ) ( src[ind] );
    }
}

template <typename ValueType>
__global__ void kernel_cudaswapdiagelemtorowfront ( ValueType *val, int *col, const int *row, const int n_rows )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i < n_rows )
    {
        for ( int j = row[i]; j < row[i + 1]; ++j )
        {
            if ( i == col[j] )
            {
                // index of first element in row i
                int j0 = row[i];

                // swap first and diagonal element
                ValueType diag = val[j];
                val[j] = val[j0];
                val[j0] = diag;

                // swap index of first and diagonal element
                int jd = col[j];
                col[j] = col[j0];
                col[j0] = jd;

                break;
            }
        }
    }
}

#endif
