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

#include "opencl/platform_management_opencl.h"

#include <assert.h>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>
#include <stdio.h>

#include "opencl/opencl_utils.h"
#include "opencl/opencl_global.h"
#include "../la_global.h"

#ifdef WITH_OPENCL

#    ifdef __APPLE__
#        include <cl.h>
#    else
#        include <CL/cl.h>
#    endif

#else

#    define ERROR LOG_ERROR("ERROR no OPENCL support");  exit(-1);

#endif

namespace hiflow
{
    namespace la
    {

        void print_platform_opencl_info ( void )
        {
#ifdef WITH_OPENCL

            cl_int status;
            cl_uint num_opencl_platforms;
            cl_platform_id *cl_platforms;
            char info_text[1024];

            // Get OpenCL platform count
            status = clGetPlatformIDs ( 0, NULL, &num_opencl_platforms );
            checkOpenCLErr ( status, "Error: clGetPlatformIDs" );

            if ( num_opencl_platforms == 0 )
            {
                LOG_ERROR ( "ERROR print_platform_opencl_info() - no available opencl platform(s)" );
                exit ( -1 );
            }

            // allocate mem for platforms
            cl_platforms = new cl_platform_id[num_opencl_platforms];

            // get platform info for each platform
            status = clGetPlatformIDs ( num_opencl_platforms, cl_platforms, NULL );
            checkOpenCLErr ( status, "Error: clGetPlatformIDs" );
            for ( cl_uint i = 0; i < num_opencl_platforms; i++ )
            {
                status = clGetPlatformInfo ( cl_platforms[i], CL_PLATFORM_NAME, 1024, &info_text, NULL );
                checkOpenCLErr ( status, "Error: clGetPlatformInfo" );
                printf ( "Platform %d:   %s\n", i, info_text );

                // check for OpenCL devices
                cl_uint num_opencl_devices;
                cl_device_id *cl_devices;
                status = clGetDeviceIDs ( cl_platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &num_opencl_devices );
                checkOpenCLErr ( status, "Error: clGetDeviceIDs" ); //TODO evtl ohne Error check wegen 0 devices
                if ( num_opencl_devices == 0 )
                {
                    printf ( "No OpenCL devices." );
                }
                else
                {
                    //allocate mem for devices
                    cl_devices = new cl_device_id[num_opencl_devices];
                    status = clGetDeviceIDs ( cl_platforms[i], CL_DEVICE_TYPE_ALL, num_opencl_devices, cl_devices, &num_opencl_devices );

                    for ( int device = 0; device < num_opencl_devices; device++ )
                    {
                        printf ( "Device %d of platform %d:\n", device, i );
                        // CL_DEVICE_NAME
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_NAME, sizeof (info_text ), &info_text, NULL );
                        printf ( "  CL_DEVICE_NAME:   %s\n", info_text );

                        // CL_DEVICE_VENDOR
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_VENDOR, sizeof (info_text ), &info_text, NULL );
                        printf ( "  CL_DEVICE_VENDOR:   %s\n", info_text );

                        // CL_DRIVER_VERSION
                        clGetDeviceInfo ( cl_devices[device], CL_DRIVER_VERSION, sizeof (info_text ), &info_text, NULL );
                        printf ( "  CL_DRIVER_VERSION:   %s\n", info_text );

                        // CL_DEVICE_INFO
                        cl_device_type type;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_TYPE, sizeof (type ), &type, NULL );
                        if ( type & CL_DEVICE_TYPE_CPU ) printf ( "  CL_DEVICE_TYPE:   %s\n", "CL_DEVICE_TYPE_CPU" );
                        if ( type & CL_DEVICE_TYPE_GPU ) printf ( "  CL_DEVICE_TYPE:   %s\n", "CL_DEVICE_TYPE_GPU" );
                        if ( type & CL_DEVICE_TYPE_ACCELERATOR ) printf ( "  CL_DEVICE_TYPE:   %s\n", "CL_DEVICE_TYPE_ACCELERATOR" );
                        if ( type & CL_DEVICE_TYPE_DEFAULT ) printf ( "  CL_DEVICE_TYPE:   %s\n", "CL_DEVICE_TYPE_DEFAULT" );

                        // CL_DEVICE_MAX_COMPUTE_UNITS
                        cl_uint compute_units;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof (compute_units ), &compute_units, NULL );
                        printf ( "  CL_DEVICE_MAX_COMPUTE_UNITS:\t\t%u\n", compute_units );

                        // CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
                        size_t workitem_dims;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof (workitem_dims ), &workitem_dims, NULL );
                        printf ( "  CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS:\t%u\n", ( unsigned int ) workitem_dims );

                        // CL_DEVICE_MAX_WORK_ITEM_SIZES
                        size_t workitem_size[3];
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof (workitem_size ), &workitem_size, NULL );
                        printf ( "  CL_DEVICE_MAX_WORK_ITEM_SIZES:\t%u / %u / %u \n", ( unsigned int ) workitem_size[0], ( unsigned int ) workitem_size[1], ( unsigned int ) workitem_size[2] );

                        // CL_DEVICE_MAX_WORK_GROUP_SIZE
                        unsigned int workgroup_size;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof (workgroup_size ), &workgroup_size, NULL );
                        printf ( "  CL_DEVICE_MAX_WORK_GROUP_SIZE:\t%u\n", workgroup_size );

                        // CL_DEVICE_MAX_CLOCK_FREQUENCY
                        cl_uint clock_frequency;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof (clock_frequency ), &clock_frequency, NULL );
                        printf ( "  CL_DEVICE_MAX_CLOCK_FREQUENCY:\t%u MHz\n", clock_frequency );

                        // CL_DEVICE_ADDRESS_BITS
                        cl_uint addr_bits;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_ADDRESS_BITS, sizeof (addr_bits ), &addr_bits, NULL );
                        printf ( "  CL_DEVICE_ADDRESS_BITS:\t\t%u\n", addr_bits );

                        // CL_DEVICE_MAX_MEM_ALLOC_SIZE
                        cl_ulong max_mem_alloc_size;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof (max_mem_alloc_size ), &max_mem_alloc_size, NULL );
                        printf ( "  CL_DEVICE_MAX_MEM_ALLOC_SIZE:\t\t%u MByte\n", ( unsigned int ) ( max_mem_alloc_size / ( 1024 * 1024 ) ) );

                        // CL_DEVICE_GLOBAL_MEM_SIZE
                        cl_ulong mem_size;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof (mem_size ), &mem_size, NULL );
                        printf ( "  CL_DEVICE_GLOBAL_MEM_SIZE:\t\t%u MByte\n", ( unsigned int ) ( mem_size / ( 1024 * 1024 ) ) );

                        // CL_DEVICE_ERROR_CORRECTION_SUPPORT
                        cl_bool error_correction_support;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof (error_correction_support ), &error_correction_support, NULL );
                        printf ( "  CL_DEVICE_ERROR_CORRECTION_SUPPORT:\t%s\n", error_correction_support == CL_TRUE ? "yes" : "no" );

                        // CL_DEVICE_LOCAL_MEM_TYPE
                        cl_device_local_mem_type local_mem_type;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_LOCAL_MEM_TYPE, sizeof (local_mem_type ), &local_mem_type, NULL );
                        printf ( "  CL_DEVICE_LOCAL_MEM_TYPE:\t\t%s\n", local_mem_type == 1 ? "local" : "global" );

                        // CL_DEVICE_LOCAL_MEM_SIZE
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_LOCAL_MEM_SIZE, sizeof (mem_size ), &mem_size, NULL );
                        printf ( "  CL_DEVICE_LOCAL_MEM_SIZE:\t\t%u KByte\n", ( unsigned int ) ( mem_size / 1024 ) );

                        // CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof (mem_size ), &mem_size, NULL );
                        printf ( "  CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE:\t%u KByte\n", ( unsigned int ) ( mem_size / 1024 ) );

                        // CL_DEVICE_QUEUE_PROPERTIES
                        cl_command_queue_properties queue_properties;
                        clGetDeviceInfo ( cl_devices[device], CL_DEVICE_QUEUE_PROPERTIES, sizeof (queue_properties ), &queue_properties, NULL );
                        if ( queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE )printf ( "  CL_DEVICE_QUEUE_PROPERTIES:\t\t%s\n", "CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE" );
                        if ( queue_properties & CL_QUEUE_PROFILING_ENABLE )printf ( "  CL_DEVICE_QUEUE_PROPERTIES:\t\t%s\n", "CL_QUEUE_PROFILING_ENABLE" );

                    }
                }
            }

            delete[] cl_platforms;

#else
            ERROR;
#endif

        }

        // init platform opencl

        void init_platform_opencl ( struct SYSTEM &my_system )
        {
#ifdef WITH_OPENCL

            my_system.opencl_initialized = false;
            my_system.my_manager = new opencl_manager;
            my_system.my_manager->init ( my_system.is_double );
            my_system.opencl_initialized = true;

#else
            ERROR;
#endif
        }

        // init platform opencl

        void init_platform_opencl ( struct SYSTEM &my_system, int opencl_plat, int opencl_dev )
        {
#ifdef WITH_OPENCL

            my_system.opencl_initialized = false;
            my_system.my_manager = new opencl_manager;
            my_system.my_manager->init ( my_system.is_double, opencl_plat, opencl_dev );
            my_system.opencl_initialized = true;

#else
            ERROR;
#endif
        }

        // stop platform opencl

        void stop_platform_opencl ( struct SYSTEM &my_system )
        {
#ifdef WITH_OPENCL

            my_system.my_manager->clear ( );
            my_system.opencl_initialized = false;

#else
            ERROR;
#endif
        }

    }
}
