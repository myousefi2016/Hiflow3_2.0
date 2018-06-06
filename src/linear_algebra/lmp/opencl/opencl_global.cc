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
#include <stdio.h>
#include <sys/time.h>
#include <fstream>
#include <string.h>

#include "opencl/opencl_utils.h"
#include "opencl/opencl_global.h"

#ifdef WITH_OPENCL

#    ifdef __APPLE__
#        include <cl.h>
#    else
#        include <CL/cl.h>
#    endif

#else

#    define ERROR LOG_ERROR("ERROR no OPENCL support");  exit(-1);

#endif

opencl_manager::opencl_manager ( )
{
    this->initialized = false;
}

opencl_manager::~opencl_manager ( )
{
}

void opencl_manager::init ( bool is_double )
{
#ifdef WITH_OPENCL

    this->cpuFinalReduction = 1; // perform final reduction on CPU or on DEVICE
    this->cpuFinalThreshold = 1; // threshold for the final reduction
    this->maxThreads = 128; // number of max threads for device
    this->maxBlocks = 64; // number of max blocks for device

    cl_int status;
    cl_device_id *cl_devices;
    cl_platform_id *cl_platforms;
    size_t deviceListSize;

    cl_uint num_opencl_platforms;
    cl_uint num_opencl_devices;
    char info_text[1024];
    cl_uint input;

    //Get OpenCL platform count
    status = clGetPlatformIDs ( 0, NULL, &num_opencl_platforms );
    checkOpenCLErr ( status, "Error: clGetPlatformIDs" );
    if ( num_opencl_platforms == 0 )
    {
        LOG_ERROR ( "ERROR - no available opencl platform(s)" ); // is this actually required here?
        exit ( -1 );
    }
    cl_platforms = new cl_platform_id[num_opencl_platforms];

    //Choose OpenCL platform
    status = clGetPlatformIDs ( num_opencl_platforms, cl_platforms, NULL );
    checkOpenCLErr ( status, "Error: clGetPlatformIDs" );
    printf ( "List of available OpenCL Platforms\n" );
    for ( cl_uint i = 0; i < num_opencl_platforms; i++ )
    {
        status = clGetPlatformInfo ( cl_platforms[i], CL_PLATFORM_NAME, 1024, &info_text, NULL );
        checkOpenCLErr ( status, "Error: clGetPlatformInfo" );
        printf ( "\t(%d)\tPlatform %d:\t%s\n", i, i, info_text );
    }

    printf ( "Which Platform do you want to use?\t" );
    std::cin >>input;
    num_opencl_platforms = input;

    //Get OpenCL devices count
    status = clGetDeviceIDs ( cl_platforms[num_opencl_platforms], CL_DEVICE_TYPE_ALL, 0, NULL, &num_opencl_devices );
    checkOpenCLErr ( status, "Error: clGetDeviceIDs" );
    if ( num_opencl_devices == 0 )
    {
        LOG_ERROR ( "ERROR - no available opencl device(s) found on this platform" ); // is this actually required here?
        exit ( -1 );
    }
    cl_devices = new cl_device_id[num_opencl_devices];

    //Choose OpenCL devices
    status = clGetDeviceIDs ( cl_platforms[num_opencl_platforms], CL_DEVICE_TYPE_ALL, num_opencl_devices, cl_devices, NULL );
    checkOpenCLErr ( status, "Error: clGetDeviceIDs" );
    printf ( "List of available OpenCL Devices for Platform %d\n", num_opencl_platforms );
    for ( cl_uint i = 0; i < num_opencl_devices; i++ )
    {
        clGetDeviceInfo ( cl_devices[i], CL_DEVICE_NAME, sizeof (info_text ), &info_text, NULL );
        printf ( "\t(%d)\tDevice %d:\t%s\n", i, i, info_text );
        cl_device_type type;
        clGetDeviceInfo ( cl_devices[i], CL_DEVICE_TYPE, sizeof (type ), &type, NULL );
        if ( type & CL_DEVICE_TYPE_CPU ) printf ( "\t\tType:\t%s\n", "CPU" );
        if ( type & CL_DEVICE_TYPE_GPU ) printf ( "\t\tType:\t%s\n", "GPU" );
        if ( type & CL_DEVICE_TYPE_ACCELERATOR ) printf ( "\t\tType:\t%s\n", "ACCELERATOR" );
        if ( type & CL_DEVICE_TYPE_DEFAULT ) printf ( "\t\tType:\t%s\n", "DEFAULT" );
    }

    printf ( "Which Device do you want to use?\t" );
    std::cin >>input;
    std::cout << std::endl;
    num_opencl_devices = input;

    //Get MAX_WORK_GROUP_SIZE
    status = clGetDeviceInfo ( cl_devices[num_opencl_devices], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof (size_t ), &max_work_group_size, NULL );
    checkOpenCLErr ( status, "Error: Getting MAX_WORK_GROUP_SIZE" );

    opencl_context = clCreateContext ( 0, 1, &cl_devices[num_opencl_devices], NULL, NULL, &status );
    checkOpenCLErr ( status, "Error: Creating Context. (clCreateContext)" );

    // First, get the size of device list data
    status = clGetContextInfo ( opencl_context, CL_CONTEXT_DEVICES, 0, NULL, &deviceListSize );
    checkOpenCLErr ( status, "Error: Getting Context Info (device list size, clGetContextInfo)" );

    devices = new cl_device_id[deviceListSize];
    status = clGetContextInfo ( opencl_context, CL_CONTEXT_DEVICES, deviceListSize, devices, NULL );
    checkOpenCLErr ( status, "Error: clGetContextInfo" );

    // Create command queue for the given device
    commandQueue = clCreateCommandQueue ( opencl_context, devices[0], 0, &status );
    checkOpenCLErr ( status, "Error: clCreateCommandQueue" );

    char *source = NULL;
    cl_device_type type;
    char vendor[1024];

    clGetDeviceInfo ( cl_devices[num_opencl_devices], CL_DEVICE_TYPE, sizeof (type ), &type, NULL );
    clGetDeviceInfo ( cl_devices[num_opencl_devices], CL_DEVICE_VENDOR, sizeof (vendor ), &vendor, NULL );

    // compile corresponding kernel source for the given device
    if ( type & CL_DEVICE_TYPE_CPU )
    {
        if ( is_double )
        {
            source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_cpu_double.cl" );
            LOG_INFO ( "opencl_global", "Loading CPU Kernel for double precision..." );
        }
        else
        {
            source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_cpu.cl" );
            LOG_INFO ( "opencl_global", "Loading CPU Kernel for single precision..." );
        }
        kerneltype = 2;
    }
    if ( type & CL_DEVICE_TYPE_GPU )
    {
        if ( vendor[0] == 'N' )
        {
            if ( is_double )
            {
                source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_nvidia_double.cl" );
                LOG_INFO ( "opencl_global", "Loading NVIDIA GPU Kernel for double precision..." );
            }
            else
            {
                source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_nvidia.cl" );
                LOG_INFO ( "opencl_global", "Loading NVIDIA GPU Kernel for single precision..." );
            }
            kerneltype = 0;
        }
        if ( vendor[0] == 'A' )
        {
            if ( is_double )
            {
                source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_ati_double.cl" );
                LOG_INFO ( "opencl_global", "Loading ATI GPU Kernel for double precision..." );
            }
            else
            {
                source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_ati.cl" );
                LOG_INFO ( "opencl_global", "Loading ATI GPU Kernel for single precision..." );
            }
            kerneltype = 1;
        }
    }

    // Compile kernel source
    my_program = clCreateProgramWithSource ( opencl_context, 1, ( const char ** ) &source, NULL, &status );
    checkOpenCLErr ( status, "Error: Loading Binary into cl_program (clCreateProgramWithBinary)" );

    status = clBuildProgram ( my_program, 1, devices, "-w -Werror -cl-mad-enable -cl-strict-aliasing", NULL, NULL );
    if ( status != CL_SUCCESS )
    {
        size_t length;
        char build_log[2048];
        printf ( "Error: Failed to build program executable!\n" );
        clGetProgramBuildInfo ( my_program, *devices, CL_PROGRAM_BUILD_LOG, sizeof (build_log ), build_log, &length );
        printf ( "%s\n", build_log );
    }
    checkOpenCLErr ( status, "Error: Building Program (clBuildProgram)" );

    // build kernels
    my_kernels[0] = clCreateKernel ( my_program, "opencl_setvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_setvalues)" );
    my_kernels[1] = clCreateKernel ( my_program, "opencl_getvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_getvalues)" );
    my_kernels[2] = clCreateKernel ( my_program, "opencl_zeros", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (zeros)" );
    my_kernels[3] = clCreateKernel ( my_program, "opencl_scale", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (scale)" );
    my_kernels[4] = clCreateKernel ( my_program, "opencl_axpy", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (axpy)" );
    my_kernels[5] = clCreateKernel ( my_program, "opencl_scaled_add", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (scaled_add)" );
    my_kernels[6] = clCreateKernel ( my_program, "opencl_dot", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (dot)" );
    my_kernels[7] = clCreateKernel ( my_program, "opencl_nrm2", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (nrm2)" );
    my_kernels[8] = clCreateKernel ( my_program, "opencl_csr_vectormult", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_csr_vectormult)" );
    my_kernels[9] = clCreateKernel ( my_program, "opencl_csr_vectormultadd", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_csr_vectormultadd)" );
    my_kernels[10] = clCreateKernel ( my_program, "Reduction", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_csr_vectormultadd)" );
    my_kernels[11] = clCreateKernel ( my_program, "DotProduct", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_csr_vectormultadd)" );
    my_kernels[12] = clCreateKernel ( my_program, "opencl_setblockvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_setblockvalues)" );
    my_kernels[13] = clCreateKernel ( my_program, "opencl_getblockvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_getblockvalues)" );
    my_kernels[14] = clCreateKernel ( my_program, "opencl_multvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_multvalues)" );

    this->initialized = true;

#else
    ERROR;
#endif
}

void opencl_manager::init ( bool is_double, int opencl_plat, int opencl_dev )
{
#ifdef WITH_OPENCL

    this->cpuFinalReduction = 1;
    this->cpuFinalThreshold = 1;
    this->maxThreads = 128;
    this->maxBlocks = 64;

    cl_int status;
    cl_device_id *cl_devices;
    cl_platform_id *cl_platforms;
    size_t deviceListSize;

    cl_uint num_opencl_platforms;
    cl_uint num_opencl_devices;
    char info_text[1024];
    cl_uint input;

    //Get OpenCL platform count
    status = clGetPlatformIDs ( 0, NULL, &num_opencl_platforms );
    checkOpenCLErr ( status, "Error: clGetPlatformIDs" );
    if ( num_opencl_platforms == 0 )
    {
        LOG_ERROR ( "ERROR - no available opencl platform(s)" ); // is this actually required here?
        exit ( -1 );
    }
    cl_platforms = new cl_platform_id[num_opencl_platforms];

    //Choose OpenCL platform
    status = clGetPlatformIDs ( num_opencl_platforms, cl_platforms, NULL );
    checkOpenCLErr ( status, "Error: clGetPlatformIDs" );

    if ( opencl_plat > num_opencl_platforms )
    {
        LOG_ERROR ( "ERROR - specified OpenCL platform does not exist" );
        exit ( -1 );
    }

    num_opencl_platforms = opencl_plat;

    //Get OpenCL devices count
    status = clGetDeviceIDs ( cl_platforms[num_opencl_platforms], CL_DEVICE_TYPE_ALL, 0, NULL, &num_opencl_devices );
    checkOpenCLErr ( status, "Error: clGetDeviceIDs" );
    if ( num_opencl_devices == 0 )
    {
        LOG_ERROR ( "ERROR - no available opencl device(s) found on this platform" ); // is this actually required here?
        exit ( -1 );
    }
    cl_devices = new cl_device_id[num_opencl_devices];

    //Choose OpenCL devices
    status = clGetDeviceIDs ( cl_platforms[num_opencl_platforms], CL_DEVICE_TYPE_ALL, num_opencl_devices, cl_devices, NULL );
    checkOpenCLErr ( status, "Error: clGetDeviceIDs" );

    if ( num_opencl_devices < opencl_dev )
    {
        LOG_ERROR ( "ERROR - specified OpenCL device does not exist" );
        exit ( -1 );
    }

    num_opencl_devices = opencl_dev;

    //Get MAX_WORK_GROUP_SIZE
    status = clGetDeviceInfo ( cl_devices[num_opencl_devices], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof (size_t ), &max_work_group_size, NULL );
    checkOpenCLErr ( status, "Error: Getting MAX_WORK_GROUP_SIZE" );

    opencl_context = clCreateContext ( 0, 1, &cl_devices[num_opencl_devices], NULL, NULL, &status );
    checkOpenCLErr ( status, "Error: Creating Context. (clCreateContext)" );

    // First, get the size of device list data
    status = clGetContextInfo ( opencl_context, CL_CONTEXT_DEVICES, 0, NULL, &deviceListSize );
    checkOpenCLErr ( status, "Error: Getting Context Info (device list size, clGetContextInfo)" );

    devices = new cl_device_id[deviceListSize];
    status = clGetContextInfo ( opencl_context, CL_CONTEXT_DEVICES, deviceListSize, devices, NULL );
    checkOpenCLErr ( status, "Error: clGetContextInfo" );

    commandQueue = clCreateCommandQueue ( opencl_context, devices[0], 0, &status );
    checkOpenCLErr ( status, "Error: clCreateCommandQueue" );

    char *source = NULL;
    cl_device_type type;
    char vendor[1024];

    clGetDeviceInfo ( cl_devices[num_opencl_devices], CL_DEVICE_TYPE, sizeof (type ), &type, NULL );
    clGetDeviceInfo ( cl_devices[num_opencl_devices], CL_DEVICE_VENDOR, sizeof (vendor ), &vendor, NULL );
    if ( type & CL_DEVICE_TYPE_CPU )
    {
        if ( is_double )
        {
            source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_cpu_double.cl" );
            LOG_INFO ( "opencl_global", "Loading CPU Kernel for double precision..." );
        }
        else
        {
            source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_cpu.cl" );
            LOG_INFO ( "opencl_global", "Loading CPU Kernel for single precision..." );
        }
        kerneltype = 2;
    }
    if ( type & CL_DEVICE_TYPE_GPU )
    {
        if ( vendor[0] == 'N' )
        {
            if ( is_double )
            {
                source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_nvidia_double.cl" );
                LOG_INFO ( "opencl_global", "Loading NVIDIA GPU Kernel for double precision..." );
            }
            else
            {
                source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_nvidia.cl" );
                LOG_INFO ( "opencl_global", "Loading NVIDIA GPU Kernel for single precision..." );
            }
            kerneltype = 0;
        }
        if ( vendor[0] == 'A' )
        {
            if ( is_double )
            {
                source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_ati_double.cl" );
                LOG_INFO ( "opencl_global", "Loading ATI GPU Kernel for double precision..." );
            }
            else
            {
                source = load_program_source ( "/home/dlukarski/hiflow3/hiflow/src/linear_algebra/lmp/opencl/kernels/kernels_ati.cl" );

                LOG_INFO ( "opencl_global", "Loading ATI GPU Kernel for single precision..." );
            }
            kerneltype = 1;
        }
    }

    my_program = clCreateProgramWithSource ( opencl_context, 1, ( const char ** ) &source, NULL, &status );
    checkOpenCLErr ( status, "Error: Loading Binary into cl_program (clCreateProgramWithBinary)" );

    status = clBuildProgram ( my_program, 1, devices, "-w -Werror -cl-mad-enable -cl-strict-aliasing", NULL, NULL );
    if ( status != CL_SUCCESS )
    {
        size_t length;
        char build_log[2048];
        printf ( "Error: Failed to build program executable!\n" );
        clGetProgramBuildInfo ( my_program, *devices, CL_PROGRAM_BUILD_LOG, sizeof (build_log ), build_log, &length );
        printf ( "%s\n", build_log );
    }
    checkOpenCLErr ( status, "Error: Building Program (clBuildProgram)" );

    // build kernels
    my_kernels[0] = clCreateKernel ( my_program, "opencl_setvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_setvalues)" );
    my_kernels[1] = clCreateKernel ( my_program, "opencl_getvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_getvalues)" );
    my_kernels[2] = clCreateKernel ( my_program, "opencl_zeros", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (zeros)" );
    my_kernels[3] = clCreateKernel ( my_program, "opencl_scale", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (scale)" );
    my_kernels[4] = clCreateKernel ( my_program, "opencl_axpy", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (axpy)" );
    my_kernels[5] = clCreateKernel ( my_program, "opencl_scaled_add", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (scaled_add)" );
    my_kernels[6] = clCreateKernel ( my_program, "opencl_dot", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (dot)" );
    my_kernels[7] = clCreateKernel ( my_program, "opencl_nrm2", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (nrm2)" );
    my_kernels[8] = clCreateKernel ( my_program, "opencl_csr_vectormult", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_csr_vectormult)" );
    my_kernels[9] = clCreateKernel ( my_program, "opencl_csr_vectormultadd", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_csr_vectormultadd)" );
    my_kernels[10] = clCreateKernel ( my_program, "Reduction", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_csr_vectormultadd)" );
    my_kernels[11] = clCreateKernel ( my_program, "DotProduct", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_csr_vectormultadd)" );
    my_kernels[12] = clCreateKernel ( my_program, "opencl_setblockvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_setblockvalues)" );
    my_kernels[13] = clCreateKernel ( my_program, "opencl_getblockvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_getblockvalues)" );
    my_kernels[14] = clCreateKernel ( my_program, "opencl_multvalues", &status );
    checkOpenCLErr ( status, "Error: Creating Kernel from program (opencl_multvalues)" );

    this->initialized = true;

#else
    ERROR;
#endif
}

void opencl_manager::clear ( )
{
#ifdef WITH_OPENCL

    cl_int status;

    //release kernels
    if ( my_kernels[0] )
    {
        status = clReleaseKernel ( my_kernels[0] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[1] )
    {
        status = clReleaseKernel ( my_kernels[1] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[2] )
    {
        status = clReleaseKernel ( my_kernels[2] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[3] )
    {
        status = clReleaseKernel ( my_kernels[3] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[4] )
    {
        status = clReleaseKernel ( my_kernels[4] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[5] )
    {
        status = clReleaseKernel ( my_kernels[5] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[6] )
    {
        status = clReleaseKernel ( my_kernels[6] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[7] )
    {
        status = clReleaseKernel ( my_kernels[7] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[8] )
    {
        status = clReleaseKernel ( my_kernels[8] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[9] )
    {
        status = clReleaseKernel ( my_kernels[9] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[10] )
    {
        status = clReleaseKernel ( my_kernels[10] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[11] )
    {
        status = clReleaseKernel ( my_kernels[11] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[12] )
    {
        status = clReleaseKernel ( my_kernels[12] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[13] )
    {
        status = clReleaseKernel ( my_kernels[13] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }
    if ( my_kernels[14] )
    {
        status = clReleaseKernel ( my_kernels[14] );
        checkOpenCLErr ( status, "Error: In clReleaseKernel\n" );
    }

    if ( my_program )
    {
        status = clReleaseProgram ( my_program );
        checkOpenCLErr ( status, "Error: In clReleaseProgram" );
    }

    if ( commandQueue )
    {
        status = clReleaseCommandQueue ( commandQueue );
        checkOpenCLErr ( status, "Error: In clReleaseCommandQueue" );
    }

    if ( opencl_context )
    {
        status = clReleaseContext ( opencl_context );
        checkOpenCLErr ( status, "Error: In clReleaseContext" );
    }

    delete[] devices;

#else
    ERROR;
#endif
}
#ifdef WITH_OPENCL

cl_context opencl_manager::get_context ( )
{
    return opencl_context;
}

cl_command_queue opencl_manager::get_command_queue ( )
{
    return commandQueue;
}

cl_kernel opencl_manager::get_kernel ( int id )
{
    return my_kernels[id];
}

cl_device_id opencl_manager::get_device ( )
{
    return devices[0];
}
#endif
