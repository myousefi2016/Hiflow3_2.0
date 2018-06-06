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

/// @author Nico Trost, Benedikt Galler, Dimitar Lukarski

#include <iostream>
#include <string.h>
#include <cmath>

#include "opencl/opencl_utils.h"
#include "opencl/opencl_kernel_mapper.h"

#ifdef WITH_OPENCL

#    include <malloc.h>
#    ifdef __APPLE__
#        include <cl.h>
#    else
#        include <CL/cl.h>
#    endif
#    include "opencl/mem_opencl.h"

void openclsetvalues ( cl_mem vector, cl_mem index_buffer, cl_mem value_buffer, int size, cl_kernel opencl_kernel,
                       cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &vector );
    checkOpenCLErr ( status, "Error: openclsetvalues (vector)" );

    status = clSetKernelArg ( opencl_kernel, 1, sizeof (unsigned int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: openclsetvalues (size)" );

    status = clSetKernelArg ( opencl_kernel, 2, sizeof (cl_mem ), ( void * ) &index_buffer );
    checkOpenCLErr ( status, "Error: openclsetvalues (index_buffer)" );

    status = clSetKernelArg ( opencl_kernel, 3, sizeof (cl_mem ), ( void * ) &value_buffer );
    checkOpenCLErr ( status, "Error: openclsetvalues (value_buffer)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (setvalues) onto command queue (clEnqueueNDRangeKernel)" );

    /* Wait for the kernel call to finish execution */
    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );
    clReleaseEvent ( my_event );
}

void openclgetvalues ( cl_mem vector, cl_mem index_buffer, cl_mem value_buffer, int size, cl_kernel opencl_kernel,
                       cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &vector );
    checkOpenCLErr ( status, "Error: openclgetvalues (vector)" );

    status = clSetKernelArg ( opencl_kernel, 1, sizeof (unsigned int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: openclgetvalues (size)" );

    status = clSetKernelArg ( opencl_kernel, 2, sizeof (cl_mem ), ( void * ) &index_buffer );
    checkOpenCLErr ( status, "Error: openclgetvalues (index_buffer)" );

    status = clSetKernelArg ( opencl_kernel, 3, sizeof (cl_mem ), ( void * ) &value_buffer );
    checkOpenCLErr ( status, "Error: openclgetvalues (value_buffer)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (getvalues) onto command queue (clEnqueueNDRangeKernel)" );

    /* Wait for the kernel call to finish execution */
    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );
    clReleaseEvent ( my_event );
}

void openclsetblockvalues ( cl_mem vector, cl_mem value_buffer, int start_i, int start_sub_vec, int size, cl_kernel opencl_kernel,
                            cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &vector );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (vector)" );
    status = clSetKernelArg ( opencl_kernel, 1, sizeof (unsigned int ), ( void * ) &start_i );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (start_i)" );
    status = clSetKernelArg ( opencl_kernel, 2, sizeof (unsigned int ), ( void * ) &start_sub_vec );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (start_sub_vec)" );
    status = clSetKernelArg ( opencl_kernel, 3, sizeof (unsigned int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (size)" );
    status = clSetKernelArg ( opencl_kernel, 4, sizeof (cl_mem ), ( void * ) &value_buffer );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (value_buffer)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (getblockvalues) onto command queue (clEnqueueNDRangeKernel)" );

    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );

    clReleaseEvent ( my_event );
}

void openclgetblockvalues ( cl_mem vector, cl_mem value_buffer, int start_i, int end_i, cl_kernel opencl_kernel,
                            cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &vector );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (vector)" );
    status = clSetKernelArg ( opencl_kernel, 1, sizeof (unsigned int ), ( void * ) &start_i );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (start_i)" );
    status = clSetKernelArg ( opencl_kernel, 2, sizeof (unsigned int ), ( void * ) &end_i );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (end_i)" );
    status = clSetKernelArg ( opencl_kernel, 3, sizeof (cl_mem ), ( void * ) &value_buffer );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (value_buffer)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (getblockvalues) onto command queue (clEnqueueNDRangeKernel)" );

    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );

    clReleaseEvent ( my_event );
}

void openclzeros ( cl_mem vector, int size, cl_kernel opencl_kernel,
                   cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &vector );
    checkOpenCLErr ( status, "Error: zeros (vector)" );

    status = clSetKernelArg ( opencl_kernel, 1, sizeof (int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: zeros (size)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (zeros) onto command queue (clEnqueueNDRangeKernel)" );

    /* Wait for the kernel call to finish execution */
    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );
    status = clReleaseEvent ( my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (zeros) onto command queue (clEnqueueNDRangeKernel)" );
}

void openclmultvalues ( cl_mem vector, cl_mem value_buffer, int size, cl_kernel opencl_kernel,
                        cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (unsigned int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (size)" );
    status = clSetKernelArg ( opencl_kernel, 1, sizeof (cl_mem ), ( void * ) &vector );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (vector)" );
    status = clSetKernelArg ( opencl_kernel, 2, sizeof (cl_mem ), ( void * ) &value_buffer );
    checkOpenCLErr ( status, "Error: openclgetblockvalues (value_buffer)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (getblockvalues) onto command queue (clEnqueueNDRangeKernel)" );

    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );

    clReleaseEvent ( my_event );
}

template <typename ValueType>
void openclscale ( cl_mem vector, int size, ValueType a, cl_kernel opencl_kernel,
                   cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &vector );
    checkOpenCLErr ( status, "Error: zeros (vector)" );

    status = clSetKernelArg ( opencl_kernel, 1, sizeof (ValueType ), ( void * ) &a );
    checkOpenCLErr ( status, "Error: zeros (a)" );

    status = clSetKernelArg ( opencl_kernel, 2, sizeof (unsigned int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: zeros (size)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (scale) onto command queue (clEnqueueNDRangeKernel)" );

    /* Wait for the kernel call to finish execution */
    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );
    status = clReleaseEvent ( my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (zeros) onto command queue (clEnqueueNDRangeKernel)" );
}

template <typename ValueType>
void openclaxpy ( cl_mem vector1, cl_mem vector2, int size, ValueType a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &vector1 );
    checkOpenCLErr ( status, "Error: zeros (vector1)" );

    status = clSetKernelArg ( opencl_kernel, 1, sizeof (cl_mem ), ( void * ) &vector2 );
    checkOpenCLErr ( status, "Error: zeros (vector2)" );

    status = clSetKernelArg ( opencl_kernel, 2, sizeof (ValueType ), ( void * ) &a );
    checkOpenCLErr ( status, "Error: zeros (a)" );

    status = clSetKernelArg ( opencl_kernel, 3, sizeof (unsigned int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: zeros (size)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (axpy) onto command queue (clEnqueueNDRangeKernel)" );

    /* Wait for the kernel call to finish execution */
    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );
    status = clReleaseEvent ( my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (axpy) onto command queue (clEnqueueNDRangeKernel)" );
}

template <typename ValueType>
void openclscaledadd ( cl_mem vector1, cl_mem vector2, int size, ValueType a, cl_kernel opencl_kernel,
                       cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // set kernel arguments
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &vector1 );
    checkOpenCLErr ( status, "Error: zeros (vector1)" );

    status = clSetKernelArg ( opencl_kernel, 1, sizeof (cl_mem ), ( void * ) &vector2 );
    checkOpenCLErr ( status, "Error: zeros (vector2)" );

    status = clSetKernelArg ( opencl_kernel, 2, sizeof (ValueType ), ( void * ) &a );
    checkOpenCLErr ( status, "Error: zeros (a)" );

    status = clSetKernelArg ( opencl_kernel, 3, sizeof (unsigned int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: zeros (size)" );

    // enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (axpy) onto command queue (clEnqueueNDRangeKernel)" );

    /* Wait for the kernel call to finish execution */
    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );
    status = clReleaseEvent ( my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (scaled add) onto command queue (clEnqueueNDRangeKernel)" );
}

// This software contains source code provided by NVIDIA Corporation.
// nextPow2 by NVIDIA

unsigned int nextPow2 ( unsigned int x )
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

bool isPow2 ( unsigned int x )
{
    return (( x & ( x - 1 ) ) == 0 );
}

// This software contains source code provided by NVIDIA Corporation.
// getNumBlocksAndThreads by NVIDIA

void getNumBlocksAndThreads ( int n, cl_int &blocks, cl_int &threads, cl_int kerneltype, cl_int maxThreads )
{
    // 0 = NVIDIA DEVICE
    if ( kerneltype == 0 )
    {
        threads = ( n < maxThreads * 2 ) ? nextPow2 ( ( n + 1 ) / 2 ) : maxThreads;
        blocks = ( n + ( threads * 2 - 1 ) ) / ( threads * 2 );
    }
    // 1 = ATI DEVICE
    // 2 = CPU DEVICE
    if ( kerneltype == 1 || kerneltype == 2 )
    {
        threads = ( n < maxThreads ) ? nextPow2 ( n ) : maxThreads;
        blocks = ( n + threads - 1 ) / threads;
    }
}

template <typename ValueType>
ValueType opencldot ( cl_mem vector1,
                      cl_mem vector2,
                      int size,
                      cl_context opencl_context,
                      cl_kernel opencl_kernel_dot,
                      cl_kernel opencl_kernel_reduction,
                      cl_command_queue commandQueue,
                      const size_t globalThreads,
                      const size_t localThreads,
                      const size_t globalBlocks,
                      cl_int kerneltype,
                      bool cpuFinalReduction,
                      cl_int cpuFinalThreshold,
                      cl_int maxThreads,
                      const size_t max_work_group_size )
{

    // call cpu function for cpu devices
    if ( kerneltype == 2 )
    {

        // no CPU version so far!
        exit ( -1 );

    }

    cl_int status, i = 0;
    cl_int numBlocks = globalBlocks, numThreads = localThreads;
    cl_bool needReadBack = true;
    cl_kernel finalReductionKernel[10];
    cl_int finalReductionIterations = 0;
    cl_event events[2];
    cl_mem DOTBuffer;
    ValueType gpu_result = 0;
    cl_int nIsPowOf2;
    ValueType *host_buffer = ( ValueType* ) malloc ( numBlocks * sizeof (ValueType ) );

    opencl_memallocdev<ValueType>( DOTBuffer, size, opencl_context );

    // Set appropriate arguments to the kernel
    nIsPowOf2 = isPow2 ( size );

    clSetKernelArg ( opencl_kernel_dot, 0, sizeof (cl_mem ), ( void * ) &vector1 );
    clSetKernelArg ( opencl_kernel_dot, 1, sizeof (cl_mem ), ( void * ) &vector2 );
    clSetKernelArg ( opencl_kernel_dot, 2, sizeof (cl_mem ), ( void * ) &DOTBuffer );
    clSetKernelArg ( opencl_kernel_dot, 3, sizeof (ValueType ) * numThreads, NULL );
    clSetKernelArg ( opencl_kernel_dot, 4, sizeof (cl_int ), &size );
    clSetKernelArg ( opencl_kernel_dot, 5, sizeof (cl_int ), &numThreads );
    clSetKernelArg ( opencl_kernel_dot, 6, sizeof (cl_int ), &nIsPowOf2 );

    if ( !cpuFinalReduction )
    {
        cl_int s = numBlocks;
        cl_int threads = 0, blocks = 0;

        while ( s > cpuFinalThreshold )
        {
            getNumBlocksAndThreads ( s, blocks, threads, kerneltype, maxThreads );
            nIsPowOf2 = isPow2 ( s );
            finalReductionKernel[finalReductionIterations] = opencl_kernel_reduction;
            // Set appropriate arguments to the kernel

            clSetKernelArg ( finalReductionKernel[finalReductionIterations], 0, sizeof (cl_mem ), ( void * ) &DOTBuffer );
            clSetKernelArg ( finalReductionKernel[finalReductionIterations], 1, sizeof (cl_mem ), ( void * ) &DOTBuffer );
            clSetKernelArg ( finalReductionKernel[finalReductionIterations], 2, sizeof (ValueType ) * numThreads, NULL );
            clSetKernelArg ( finalReductionKernel[finalReductionIterations], 3, sizeof (cl_int ), &size );
            clSetKernelArg ( finalReductionKernel[finalReductionIterations], 4, sizeof (cl_int ), &threads );

            //s = (s + (threads*2-1)) / (threads*2);
            s = blocks;
            finalReductionIterations++;
        }
    }

    size_t globalWorkSize[1];
    size_t localWorkSize[1];

    gpu_result = 0;

    // Execute the kernel
    globalWorkSize[0] = numBlocks * numThreads;
    localWorkSize[0] = numThreads;

    // Enqueue a kernel run call
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel_dot, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &events[0] );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (DotProduct) onto command queue (clEnqueueNDRangeKernel)" );

    // Wait for the kernel call to finish execution
    status = clWaitForEvents ( 1, &events[0] );
    checkOpenCLErr ( status, "Error: Waiting for kernel (DotProduct) run to finish (clWaitForEvents)" );
    clReleaseEvent ( events[0] );

    if ( cpuFinalReduction )
    {
        // sum partial sums from each block on CPU
        // copy result from device to host
        clEnqueueReadBuffer ( commandQueue, DOTBuffer, CL_TRUE, 0, numBlocks * sizeof (ValueType ), host_buffer, 0, NULL, NULL );

        int num_sum = 0;
        for ( int i = 0; i < numBlocks; i++ )
        {
            gpu_result += host_buffer[i];
            num_sum++;
        }
        needReadBack = false;
    }
    else
    {
        // Sum partial block sums on GPU
        int s = numBlocks;
        int it = 0;

        while ( s > cpuFinalThreshold )
        {
            int threads = 0, blocks = 0;
            getNumBlocksAndThreads ( s, blocks, threads, kerneltype, maxThreads );

            globalWorkSize[0] = threads * blocks;
            localWorkSize[0] = threads;

            // Enqueue a kernel run call
            status = clEnqueueNDRangeKernel ( commandQueue, finalReductionKernel[it], 1, 0, globalWorkSize, localWorkSize, 0, NULL, &events[0] );
            checkOpenCLErr ( status, "Error: Enqueueing kernel (Reduction) onto command queue (clEnqueueNDRangeKernel)" );

            //s = (s + (threads*2-1)) / (threads*2);
            s = blocks;
            it++;

            status = clWaitForEvents ( 1, &events[0] );
            checkOpenCLErr ( status, "Error: Waiting for kernel (Reduction) run to finish (clWaitForEvents)" );
            clReleaseEvent ( events[0] );
        }

        if ( s > 1 )
        {
            // copy result from device to host
            clEnqueueReadBuffer ( commandQueue, DOTBuffer, CL_TRUE, 0, s * sizeof (ValueType ), host_buffer, 0, NULL, NULL );
            for ( int i = 0; i < s; i++ ) gpu_result += host_buffer[i];
            needReadBack = false;
        }
    }

    if ( needReadBack )
    {
        // copy final sum from device to host
        clEnqueueReadBuffer ( commandQueue, DOTBuffer, CL_TRUE, 0, sizeof (ValueType ), &gpu_result, 0, NULL, NULL );
    }

    // cleanup
    opencl_memfreedev ( DOTBuffer );
    free ( host_buffer );

    return gpu_result;
}

template <typename ValueType>
ValueType openclnrm2 ( cl_mem vector1, int size, cl_context opencl_context, cl_kernel opencl_kernel_dot, cl_kernel opencl_kernel_reduction, cl_command_queue commandQueue, size_t globalThreads, const size_t localThreads, const size_t globalBlocks, int kerneltype, bool cpuFinalReduction, cl_int cpuFinalThreshold, cl_int maxThreads, const size_t max_work_group_size )
{
    ValueType result = opencldot<ValueType>( vector1, vector1, size, opencl_context, opencl_kernel_dot, opencl_kernel_reduction, commandQueue, globalThreads, localThreads, globalBlocks, kerneltype, cpuFinalReduction, cpuFinalThreshold, maxThreads, max_work_group_size );

    return sqrt ( result ); //TODO: sqrt on device, hence we need a special reduction kernel
}

void openclcsrvectormult ( cl_mem val, cl_mem col, cl_mem row, cl_mem vector1, cl_mem vector2, int size, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads )
{
    cl_event my_event;
    cl_int status;

    // Set arguements for kernel call
    status = clSetKernelArg ( opencl_kernel, 0, sizeof (cl_mem ), ( void * ) &val );
    checkOpenCLErr ( status, "Error: zeros (val)" );

    status = clSetKernelArg ( opencl_kernel, 1, sizeof (cl_mem ), ( void * ) &col );
    checkOpenCLErr ( status, "Error: zeros (col)" );

    status = clSetKernelArg ( opencl_kernel, 2, sizeof (cl_mem ), ( void * ) &row );
    checkOpenCLErr ( status, "Error: zeros (row)" );

    status = clSetKernelArg ( opencl_kernel, 3, sizeof (cl_mem ), ( void * ) &vector1 );
    checkOpenCLErr ( status, "Error: zeros (vector1)" );

    status = clSetKernelArg ( opencl_kernel, 4, sizeof (cl_mem ), ( void * ) &vector2 );
    checkOpenCLErr ( status, "Error: zeros (vector2)" );

    status = clSetKernelArg ( opencl_kernel, 5, sizeof (unsigned int ), ( void * ) &size );
    checkOpenCLErr ( status, "Error: zeros (size)" );

    // Enqueue kernel run
    status = clEnqueueNDRangeKernel ( commandQueue, opencl_kernel, 1, NULL, &globalThreads, &localThreads, 0, NULL, &my_event ); //&localThreads
    checkOpenCLErr ( status, "Error: Enqueueing kernel (csr vectormult) onto command queue (clEnqueueNDRangeKernel)" );

    /* Wait for the kernel call to finish execution */
    status = clWaitForEvents ( 1, &my_event );
    checkOpenCLErr ( status, "Error: Waiting for kernel run to finish (clWaitForEvents)" );
    status = clReleaseEvent ( my_event );
    checkOpenCLErr ( status, "Error: Enqueueing kernel (csr vectormult) onto command queue (clEnqueueNDRangeKernel)" );
}

template float opencldot<float>( cl_mem vector1, cl_mem vector2, int size, cl_context opencl_context, cl_kernel opencl_kernel_dot, cl_kernel opencl_kernel_reduction, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads, const size_t globalBlocks, int kerneltype, bool cpuFinalReduction, cl_int cpuFinalThreshold, cl_int maxThreads, const size_t max_work_group_size );
template float openclnrm2<float>( cl_mem vector1, int size, cl_context opencl_context, cl_kernel opencl_kernel_dot, cl_kernel opencl_kernel_reduction, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads, const size_t globalBlocks, int kerneltype, bool cpuFinalReduction, cl_int cpuFinalThreshold, cl_int maxThreads, const size_t max_work_group_size );
template void openclscale ( cl_mem vector, int size, float a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );
template void openclaxpy ( cl_mem vector1, cl_mem vector2, int size, float a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );
template void openclscaledadd ( cl_mem vector1, cl_mem vector2, int size, float a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

template double opencldot<double>( cl_mem vector1, cl_mem vector2, int size, cl_context opencl_context, cl_kernel opencl_kernel_dot, cl_kernel opencl_kernel_reduction, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads, const size_t globalBlocks, int kerneltype, bool cpuFinalReduction, cl_int cpuFinalThreshold, cl_int maxThreads, const size_t max_work_group_size );
template double openclnrm2<double>( cl_mem vector1, int size, cl_context opencl_context, cl_kernel opencl_kernel_dot, cl_kernel opencl_kernel_reduction, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads, const size_t globalBlocks, int kerneltype, bool cpuFinalReduction, cl_int cpuFinalThreshold, cl_int maxThreads, const size_t max_work_group_size );
template void openclscale ( cl_mem vector, int size, double a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );
template void openclaxpy ( cl_mem vector1, cl_mem vector2, int size, double a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );
template void openclscaledadd ( cl_mem vector1, cl_mem vector2, int size, double a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

#else

#    define ERROR LOG_ERROR("ERROR no OPENCL support");  exit(-1);

#endif
