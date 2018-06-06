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

#include <assert.h>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>
#include "opencl/mem_opencl.h"

#ifdef WITH_OPENCL

#    ifdef __APPLE__
#        include <cl.h>
#    else
#        include <CL/cl.h>
#    endif
#    include "opencl/opencl_utils.h"

template<typename ValueType>
void opencl_memcpy2dev ( cl_mem dest, const ValueType *src, const int size, cl_command_queue commandQueue )
{

    if ( size > 0 )
    {
        cl_int status;
        cl_event my_event;
        status = clEnqueueWriteBuffer ( commandQueue, dest, CL_TRUE, 0, sizeof (ValueType ) * size, src, 0, NULL, &my_event );
        checkOpenCLErr ( status, "Error: memcpy2dev (clEnqueueWriteBuffer)" );
        status = clWaitForEvents ( 1, &my_event );
        checkOpenCLErr ( status, "Error: memcpy2dev (clWaitForEvents)" );
        status = clReleaseEvent ( my_event );
        checkOpenCLErr ( status, "Error: memcpy2dev (clReleaseEvent)" );
    }

}

template<typename ValueType>
void opencl_memcpy2host ( ValueType *dest, cl_mem src, const int size, cl_command_queue commandQueue )
{

    if ( size > 0 )
    {
        cl_int status;
        cl_event my_event;
        status = clEnqueueReadBuffer ( commandQueue, src, CL_TRUE, 0, sizeof (ValueType ) * size, dest, 0, NULL, &my_event );
        checkOpenCLErr ( status, "Error: memcpy2host (clEnqueueReadBuffer)" );
        status = clWaitForEvents ( 1, &my_event );
        checkOpenCLErr ( status, "Error: memcpy2host (clWaitForEvents)" );
        status = clReleaseEvent ( my_event );
        checkOpenCLErr ( status, "Error: memcpy2host (clReleaseEvent)" );
    }

}

template<typename ValueType>
void opencl_memcpydev ( cl_mem dest, cl_mem src, const int size, cl_command_queue commandQueue )
{

    if ( size > 0 )
    {
        cl_int status;
        cl_event my_event;
        status = clEnqueueCopyBuffer ( commandQueue, src, dest, 0, 0, sizeof (ValueType ) * size, 0, NULL, &my_event );
        checkOpenCLErr ( status, "Error: memcpydev (clEnqueueCopyBuffer)" );
        status = clWaitForEvents ( 1, &my_event );
        checkOpenCLErr ( status, "Error: memcpydev (clWaitForEvents)" );
        status = clReleaseEvent ( my_event );
        checkOpenCLErr ( status, "Error: memcpydev (clReleaseEvent)" );
    }

}

void opencl_memfreedev ( cl_mem p )
{

    cl_int status;
    status = clReleaseMemObject ( p );
    checkOpenCLErr ( status, "Error: memfreedev (clReleaseMemObject)" );
}

template<typename ValueType>
void opencl_memallocdev ( cl_mem &p, const unsigned int size, cl_context opencl_context )
{
    if ( size > 0 )
    {
        cl_int status;
        p = clCreateBuffer ( opencl_context, CL_MEM_READ_WRITE, sizeof (ValueType ) * size, NULL, &status );
        checkOpenCLErr ( status, "Error: memallocdev (clCreateBuffer)" );
    }

}

void opencl_memsetdev ( void *p, const int value, const size_t size, const size_t size_type )
{
    LOG_ERROR ( "memsetdev not implemented yet" );
    exit ( -1 );
}

template void opencl_memcpy2dev<double>( cl_mem dest, const double *src, const int size, cl_command_queue commandQueue );
template void opencl_memcpy2host<double>( double *dest, cl_mem src, const int size, cl_command_queue commandQueue );
template void opencl_memcpydev<double>( cl_mem dest, cl_mem src, const int size, cl_command_queue commandQueue );
template void opencl_memallocdev<double>( cl_mem &p, const unsigned int size, cl_context opencl_context );

template void opencl_memcpy2dev<float>( cl_mem dest, const float *src, const int size, cl_command_queue commandQueue );
template void opencl_memcpy2host<float>( float *dest, cl_mem src, const int size, cl_command_queue commandQueue );
template void opencl_memcpydev<float>( cl_mem dest, cl_mem src, const int size, cl_command_queue commandQueue );
template void opencl_memallocdev<float>( cl_mem &p, const unsigned int size, cl_context opencl_context );

template void opencl_memcpy2dev<int>( cl_mem dest, const int *src, const int size, cl_command_queue commandQueue );
template void opencl_memcpy2host<int>( int *dest, cl_mem src, const int size, cl_command_queue commandQueue );
template void opencl_memcpydev<int>( cl_mem dest, cl_mem src, const int size, cl_command_queue commandQueue );
template void opencl_memallocdev<int>( cl_mem &p, const unsigned int size, cl_context opencl_context );

#else

#    define ERROR LOG_ERROR("ERROR no OPENCL support");  exit(-1);

#endif
