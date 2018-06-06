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

// Copyright (C) 2011-2015 Vincent Heuveline

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-06-27

#include <stdio.h>

/**
 * CUDA Kernel Device code
 *
 * Static loop unrolling for the thread within one warp.
 */
template <unsigned int block_size>
__device__ void warpReduce ( volatile double *data, int tid )
{
    if ( block_size >= 64 ) data[tid] += data[tid + 32];
    if ( block_size >= 32 ) data[tid] += data[tid + 16];
    if ( block_size >= 16 ) data[tid] += data[tid + 8];
    if ( block_size >= 8 ) data[tid] += data[tid + 4];
    if ( block_size >= 4 ) data[tid] += data[tid + 2];
    if ( block_size >= 2 ) data[tid] += data[tid + 1];
}

/**
 * CUDA Kernel Device code
 *
 * Computes the dot product of two vectors.
 */
template <unsigned int block_size>
__global__ void dot_kernel ( double *v1, double *v2, double *dot, int size )
{
    int tid = threadIdx.x;
    double sum = 0.0;

    __shared__ double dot_local[block_size];

    for ( int i = tid; i < size; i += block_size )
    {
        sum += v1[i] * v2[i];
    }

    dot_local[tid] = sum;
    __syncthreads ( );

    // Parallel reduction
    if ( block_size >= 512 )
    {
        if ( tid < 256 )
        {
            dot_local[tid] += dot_local[tid + 256];
        }
        __syncthreads ( );
    }
    if ( block_size >= 256 )
    {
        if ( tid < 128 )
        {
            dot_local[tid] += dot_local[tid + 128];
        }
        __syncthreads ( );
    }
    if ( block_size >= 128 )
    {
        if ( tid < 64 )
        {
            dot_local[tid] += dot_local[tid + 64];
        }
        __syncthreads ( );
    }

    // Static loop unrolling for the thread within one warp.
    if ( tid < 32 ) warpReduce<block_size>( dot_local, tid );

    // Copy accumulated local value to global array firstStep
    if ( tid == 0 ) *dot = dot_local[0];
}
