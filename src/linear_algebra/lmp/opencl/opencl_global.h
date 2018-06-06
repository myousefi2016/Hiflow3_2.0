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

#ifndef __OPENCL_GLOBAL_H
#    define __OPENCL_GLOBAL_H

#    include "config.h"

#    ifdef WITH_OPENCL

#        ifdef __APPLE__
#            include <cl.h>
#        else
#            include <CL/cl.h>
#        endif

#    endif

class opencl_manager
{
  private:

#    ifdef WITH_OPENCL
    cl_context opencl_context; //Note: cl_... seems to be already a pointer (see CL/cl.h)
    cl_command_queue commandQueue;
    cl_program my_program;
    cl_kernel my_kernels[15];
    cl_device_id *devices;
#    endif
    bool initialized;
  public:
    opencl_manager ( );
    ~opencl_manager ( );

    void init ( bool is_double );
    void init ( bool is_double, int opencl_plat, int opencl_dev );
    bool is_initialized ( );
    void clear ( );

#    ifdef WITH_OPENCL
    cl_context get_context ( );
    cl_command_queue get_command_queue ( );
    cl_kernel get_kernel ( int id );
    cl_device_id get_device ( );
    cl_int kerneltype;
    cl_int cpuFinalReduction;
    cl_int cpuFinalThreshold;
    cl_int maxThreads;
    cl_int maxBlocks;
    size_t max_work_group_size;
#    endif

};

#endif
