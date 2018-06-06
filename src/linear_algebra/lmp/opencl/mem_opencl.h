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

#ifndef __MEM_OPENCL_H
#    define __MEM_OPENCL_H

#    include "config.h"

#    ifdef WITH_OPENCL

#        include "opencl_utils.h"

template<typename ValueType>
void opencl_memcpy2dev ( cl_mem dest, const ValueType *src, const int size, cl_command_queue commandQueue );

template<typename ValueType>
void opencl_memcpy2host ( ValueType *dest, cl_mem src, const int size, cl_command_queue commandQueue );

template<typename ValueType>
void opencl_memcpydev ( cl_mem dest, cl_mem src, const int size, cl_command_queue commandQueue );

void opencl_memfreedev ( cl_mem p );

template<typename ValueType>
void opencl_memallocdev ( cl_mem &p, const unsigned int size, cl_context opencl_context );

void opencl_memsetdev ( void *p, const int value, const size_t size, const size_t size_type );

#    endif

#endif
