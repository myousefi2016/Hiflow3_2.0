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

#include "../lmp_mem.h"
#include "../lmp_log.h"

#include <assert.h>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

#ifndef NOCUDA

#include <cuda.h>

#include <cuda_runtime_api.h>

#include "cuda/cuda_utils.h"

#include "cuda/mem_cuda.h"

#else

#define ERROR LOG_ERROR("no Cuda support");  exit(-1);

#endif

template<typename ValueType>
void memcpy2dev(ValueType *dest, const ValueType *src, const size_t size)
{
#ifndef NOCUDA

  cudaThreadSynchronize();

  cudaMemcpy(dest, src, sizeof(ValueType)*size, cudaMemcpyHostToDevice) ;
  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template<typename ValueType>
void memcpy2host(ValueType *dest, const ValueType *src, const size_t size)
{
#ifndef NOCUDA

  cudaThreadSynchronize();

  cudaMemcpy(dest , src, sizeof(ValueType)*size, cudaMemcpyDeviceToHost);
  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template<typename ValueType>
void memcpydev(ValueType *dest, const ValueType *src, const size_t size)
{
#ifndef NOCUDA

  cudaThreadSynchronize();

  cudaMemcpy(dest, src, sizeof(ValueType)*size, cudaMemcpyDeviceToDevice);
  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template<typename ValueType>
void memfreedev(ValueType *p)
{
#ifndef NOCUDA

  cudaThreadSynchronize();

  CUDA_SAFE_CALL(cudaFree(p)) ;
#else
  ERROR ;
#endif
}

void memallocdev(void **p, const size_t size, const size_t size_type)
{
#ifndef NOCUDA

  cudaThreadSynchronize();

  cudaMalloc((void**) p, size*size_type) ;
  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

void memsetdev(void *p, const int value, const size_t size, const size_t size_type)
{
#ifndef NOCUDA

  cudaThreadSynchronize();

  cudaMemset(p, value, size*size_type);
  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudasetvalues(const int *index, const size_t size, const ValueType *values, ValueType *buffer,
                   int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = size / thread_block_size ;
  if (thread_block_size*thread_block < (int)(size))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudasetvalues<<<dimGrid,dimBlock>>>(index, size, values, buffer) ;

  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudagetvalues(const int *index, const size_t size, ValueType *values, const ValueType *buffer,
                   int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = size / thread_block_size ;
  if (thread_block_size*thread_block < (int)(size) )
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudagetvalues<<<dimGrid,dimBlock>>>(index, size, values, buffer) ;

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudasetblockvalues(const int start_i, const int start_sub_vec, const int size, const ValueType *values, ValueType *buffer,
                   int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (size) / thread_block_size ;
  if (thread_block_size*thread_block < (size))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudasetblockvalues<<<dimGrid,dimBlock>>>(start_i, start_sub_vec, size, values, buffer) ;

  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudagetblockvalues(const int start_i, const int end_i, ValueType *values, const ValueType *buffer,
                   int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (end_i-start_i) / thread_block_size ;
  if (thread_block_size*thread_block < (end_i-start_i))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudagetblockvalues<<<dimGrid,dimBlock>>>(start_i, end_i, values, buffer) ;

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudaaddblockvalues(const int start_i, const int start_sub_vec, const int size, const ValueType *values, ValueType *buffer,
                        const ValueType weight, int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (size) / thread_block_size ;
  if (thread_block_size*thread_block < (size))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudaaddblockvalues<<<dimGrid,dimBlock>>>(start_i, start_sub_vec, size, values, buffer, weight) ;

  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

void cuda_sync_threads(void)
{
#ifndef NOCUDA

  cudaThreadSynchronize();

#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudamultvalues(const int size, const ValueType *values, ValueType *buffer, int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (size) / thread_block_size ;
  if (thread_block_size*thread_block < (size))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudamultvalues<<<dimGrid,dimBlock>>>(size, values, buffer) ;

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudacastfromdouble(const int size, ValueType *dest, const double *src, int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (size) / thread_block_size ;
  if (thread_block_size*thread_block < (size))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudacastfromdouble<<<dimGrid,dimBlock>>>(size, dest, src);

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudacastfromfloat(const int size, ValueType *dest, const float *src, int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (size) / thread_block_size ;
  if (thread_block_size*thread_block < (size))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudacastfromfloat<<<dimGrid,dimBlock>>>(size, dest, src);

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudacasttodouble(const int size, double *dest, const ValueType *src, int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (size) / thread_block_size ;
  if (thread_block_size*thread_block < (size))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudacasttodouble<<<dimGrid,dimBlock>>>(size, dest, src);

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudacasttofloat(const int size, float *dest, const ValueType *src, int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (size) / thread_block_size ;
  if (thread_block_size*thread_block < (size))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudacasttofloat<<<dimGrid,dimBlock>>>(size, dest, src);

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template<typename ValueType>
void cudaswapdiagelemtorowfront(ValueType *val, int *col, const int *row, const int n_rows, int thread_block_size)
{
#ifndef NOCUDA

  int thread_block = (n_rows) / thread_block_size ;
  if (thread_block_size * thread_block < (n_rows))
    thread_block++ ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  cudaThreadSynchronize();

  kernel_cudaswapdiagelemtorowfront<<<dimGrid,dimBlock>>>(val, col, row, n_rows);

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template void memcpy2dev<double   >(double    *dest, const double    *src, const size_t size);
template void memcpy2dev<float    >(float     *dest, const float     *src, const size_t size);
template void memcpy2dev<int      >(int       *dest, const int       *src, const size_t size);
template void memcpy2dev<short int>(short int *dest, const short int *src, const size_t size);

template void memcpy2host<double   >(double    *dest, const double    *src, const size_t size);
template void memcpy2host<float    >(float     *dest, const float     *src, const size_t size);
template void memcpy2host<int      >(int       *dest, const int       *src, const size_t size);
template void memcpy2host<short int>(short int *dest, const short int *src, const size_t size);

template void memcpydev<double   >(double    *dest, const double    *src, const size_t size);
template void memcpydev<float    >(float     *dest, const float     *src, const size_t size);
template void memcpydev<int      >(int       *dest, const int       *src, const size_t size);
template void memcpydev<short int>(short int *dest, const short int *src, const size_t size);

template void memfreedev<double   >(double    *dest);
template void memfreedev<float    >(float     *dest);
template void memfreedev<int      >(int       *dest);
template void memfreedev<short int>(short int *dest);

template void cudagetvalues(const int *index, const size_t size, double    *values, const double    *buffer,
                            int thread_block_size);
template void cudagetvalues(const int *index, const size_t size, float     *values, const float     *buffer,
                            int thread_block_size);
template void cudagetvalues(const int *index, const size_t size, int       *values, const int       *buffer,
                            int thread_block_size);
template void cudagetvalues(const int *index, const size_t size, short int *values, const short int *buffer,
                            int thread_block_size);

template void cudasetvalues(const int *index, const size_t size, const double    *values, double    *buffer,
                            int thread_block_size);
template void cudasetvalues(const int *index, const size_t size, const float     *values, float     *buffer,
                            int thread_block_size);
template void cudasetvalues(const int *index, const size_t size, const int       *values, int       *buffer,
                            int thread_block_size);
template void cudasetvalues(const int *index, const size_t size, const short int *values, short int *buffer,
                            int thread_block_size);

template void cudagetblockvalues(const int start_i, const int end_i,  double    *values, const double    *buffer,
                            int thread_block_size);
template void cudagetblockvalues(const int start_i, const int end_i,  float     *values, const float     *buffer,
                            int thread_block_size);
template void cudagetblockvalues(const int start_i, const int end_i,  int       *values, const int       *buffer,
                            int thread_block_size);
template void cudagetblockvalues(const int start_i, const int end_i,  short int *values, const short int *buffer,
                            int thread_block_size);

template void cudasetblockvalues(const int start_i, const int start_sub_vec, const int size,  const double    *values, double    *buffer,
                            int thread_block_size);
template void cudasetblockvalues(const int start_i, const int start_sub_vec, const int size,  const float     *values, float     *buffer,
                            int thread_block_size);
template void cudasetblockvalues(const int start_i, const int start_sub_vec, const int size,  const int       *values, int       *buffer,
                            int thread_block_size);
template void cudasetblockvalues(const int start_i, const int start_sub_vec, const int size,  const short int *values, short int *buffer,
                            int thread_block_size);

template void cudaaddblockvalues(const int start_i, const int start_sub_vec, const int size,  const double *values, double *buffer,
                                 const double weight, int thread_block_size);
template void cudaaddblockvalues(const int start_i, const int start_sub_vec, const int size,  const float *values, float *buffer,
                                 const float weight, int thread_block_size);
template void cudaaddblockvalues(const int start_i, const int start_sub_vec, const int size,  const int *values, int *buffer,
                                 const int weight, int thread_block_size);
template void cudaaddblockvalues(const int start_i, const int start_sub_vec, const int size,  const short int *values, short int *buffer,
                                 const short int weight, int thread_block_size);

template void cudamultvalues(const int size, const double *values, double *buffer, int thread_block_size);
template void cudamultvalues(const int size, const float *values, float *buffer, int thread_block_size);
template void cudamultvalues(const int size, const int *values, int *buffer, int thread_block_size);
template void cudamultvalues(const int size, const short int *values, short int *buffer, int thread_block_size);

template void cudacastfromdouble<double>(const int, double*, const double*, int);
template void cudacastfromfloat<double> (const int, double*, const float*,  int);
template void cudacasttodouble<double>  (const int, double*, const double*, int);
template void cudacasttofloat<double>   (const int, float*,  const double*, int);

template void cudacastfromdouble<float>(const int, float*,  const double*, int);
template void cudacastfromfloat<float> (const int, float*,  const float*,  int);
template void cudacasttodouble<float>  (const int, double*, const float*,  int);
template void cudacasttofloat<float>   (const int, float*,  const float*,  int);

template void cudaswapdiagelemtorowfront<double>(double*, int*, const int*, const int, int);
template void cudaswapdiagelemtorowfront<float> (float* , int*, const int*, const int, int);
