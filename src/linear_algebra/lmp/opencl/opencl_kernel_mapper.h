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

#ifndef __OPENCL_KERNEL_MAPPER_H
#    define __OPENCL_KERNEL_MAPPER_H
#    include "config.h"
void openclsetvalues ( cl_mem vector, cl_mem index_buffer, cl_mem value_buffer, int size, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

void openclgetvalues ( cl_mem vector, cl_mem index_buffer, cl_mem value_buffer, int size, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

void openclsetblockvalues ( cl_mem vector, cl_mem value_buffer, int start_i, int start_sub_vec, int size, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

void openclgetblockvalues ( cl_mem vector, cl_mem value_buffer, int start_i, int end_i, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

void openclmultvalues ( cl_mem vector, cl_mem value_buffer, int size, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

void openclzeros ( cl_mem vector, int size, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

template <typename ValueType>
void openclscale ( cl_mem vector, int size, ValueType a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

template <typename ValueType>
void openclaxpy ( cl_mem vector1, cl_mem vector2, int size, ValueType a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

template <typename ValueType>
void openclscaledadd ( cl_mem vector1, cl_mem vector2, int size, ValueType a, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

template <typename ValueType>
ValueType opencldot_cpu ( cl_mem vector1, cl_mem vector2, int size, cl_context opencl_context, cl_kernel opencl_kernel_dot, cl_kernel opencl_kernel_reduction, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads, const size_t globalBlocks, int kerneltype, bool cpuFinalReduction, cl_int cpuFinalThreshold, cl_int maxThreads, const size_t max_work_group_size );

template <typename ValueType>
ValueType opencldot ( cl_mem vector1, cl_mem vector2, int size, cl_context opencl_context, cl_kernel opencl_kernel_dot, cl_kernel opencl_kernel_reduction, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads, const size_t globalBlocks, int kerneltype, bool cpuFinalReduction, cl_int cpuFinalThreshold, cl_int maxThreads, const size_t max_work_group_size );

template <typename ValueType>
ValueType openclnrm2 ( cl_mem vector1, int size, cl_context opencl_context, cl_kernel opencl_kernel_dot, cl_kernel opencl_kernel_reduction, cl_command_queue commandQueue, size_t globalThreads, const size_t localThreads, const size_t globalBlocks, int kerneltype, bool cpuFinalReduction, cl_int cpuFinalThreshold, cl_int maxThreads, const size_t max_work_group_size );

void openclcsrvectormult ( cl_mem val, cl_mem col, cl_mem row, cl_mem vector1, cl_mem vector2, int size, cl_kernel opencl_kernel, cl_command_queue commandQueue, const size_t globalThreads, const size_t localThreads );

void getNumBlocksAndThreads ( int n, cl_int &blocks, cl_int &threads, cl_int kerneltype, cl_int maxThreads );

void create_reduction_pass_counts ( int count, int max_group_size, int max_groups, int max_work_items, int *pass_count, size_t **group_counts, size_t **work_item_counts, int **operation_counts, int **entry_counts );

unsigned int nextPow2 ( unsigned int x );

bool isPow2 ( unsigned int x );

#endif
