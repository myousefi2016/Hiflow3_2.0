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

#include "../la_global.h"
#include "cuda/platform_management_gpu.h"
#include "../lmp_log.h"

#include <assert.h>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

#ifndef NOCUDA

#include <cuda.h>

#include <cublas.h>

#include "cuda/cuda_utils.h"

#else

#define ERROR LOG_ERROR("no Cuda support");  exit(-1);

#endif

void print_platform_gpu_info(void)
{
#ifndef NOCUDA
  int device;
  int dCount;
  //      cudaError_t cuda_status ;
  cudaGetDeviceCount(&dCount);
  for (device = 0; device < dCount; ++device)
    {
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device);

      LOG_INFO("lmpLAtoolbox GPU platform","device:" << device);
      LOG_INFO("lmpLAtoolbox GPU platform","name:" << deviceProp.name);
      LOG_INFO("lmpLAtoolbox GPU platform","totalGlobalMem [bytes]: " << deviceProp.totalGlobalMem);
      LOG_INFO("lmpLAtoolbox GPU platform","sharedMemPerBlock [bytes]:" << deviceProp.sharedMemPerBlock);
      LOG_INFO("lmpLAtoolbox GPU platform","regsPerBlock:" << deviceProp.regsPerBlock);
      LOG_INFO("lmpLAtoolbox GPU platform","warpSize:" << deviceProp.warpSize);
      LOG_INFO("lmpLAtoolbox GPU platform","memPitch:" << deviceProp.memPitch);
      LOG_INFO("lmpLAtoolbox GPU platform","maxThreadsPerBlock:" << deviceProp.maxThreadsPerBlock);
      LOG_INFO("lmpLAtoolbox GPU platform","maxThreadsDim[0]:" << deviceProp.maxThreadsDim[0]);
      LOG_INFO("lmpLAtoolbox GPU platform","maxThreadsDim[1]:" << deviceProp.maxThreadsDim[1]);
      LOG_INFO("lmpLAtoolbox GPU platform","maxThreadsDim[2]:" << deviceProp.maxThreadsDim[2]);
      LOG_INFO("lmpLAtoolbox GPU platform","maxGridSize[0]:" << deviceProp.maxGridSize[0]);
      LOG_INFO("lmpLAtoolbox GPU platform","maxGridSize[1]:" << deviceProp.maxGridSize[1]);
      LOG_INFO("lmpLAtoolbox GPU platform","maxGridSize[2]:" << deviceProp.maxGridSize[2]);
      LOG_INFO("lmpLAtoolbox GPU platform","totalConstMem [bytes]: " << deviceProp.totalConstMem);
      LOG_INFO("lmpLAtoolbox GPU platform","compute capability:" << deviceProp.major << "." << deviceProp.minor);
      LOG_INFO("lmpLAtoolbox GPU platform","clockRate [kHz]: " << deviceProp.clockRate);
      LOG_INFO("lmpLAtoolbox GPU platform","textureAlignment" << deviceProp.textureAlignment);

    }
#else
  ERROR ;
#endif
}

// init platform GPU
void init_platform_gpu(bool &GPU_CUBLAS)
{
#ifndef NOCUDA
  // Select the next available GPU/CUDA device
  int my_dev = -1 ;
  int device;
  int dCount;
  cudaError_t cuda_status ;

  //  print_platform_gpu_info();

  LOG_INFO("init_platform"," checking for GPUs ...");

  cudaGetDeviceCount(&dCount);

  for (device = 0; device < dCount; ++device)
    {
      cuda_status = cudaSetDevice(device);

      if (cuda_status == cudaSuccess) {

        if ( GPU_CUBLAS )
            cublasInit(); // init CUDA BLAS

        if ( cublasGetError() == CUBLAS_STATUS_SUCCESS) {
          my_dev = device ;
          break ;
        }
      }
    }

  if (my_dev == -1) {
    LOG_INFO("init_platform" ,"no available GPU device(s)");
    exit(-1);
  }

  LOG_INFO("init_platform", " checking - done");

  // clear the last errors
  cudaGetLastError();

  LOG_INFO("init_platform", "Using GPU device=" << my_dev);

  cudaThreadExit();
  cudaSetDevice(my_dev);
  Cuda_Error_Check ;

  LOG_INFO("init_platform", "set device -  done");

  if ( GPU_CUBLAS )
    {

      LOG_INFO("init_platform", "init CUBLAS ..." );

      cublasInit(); // init CUDA BLAS
      Cuda_Blas_Error_Check ;

      LOG_INFO("init_platform", "init CUBLAS - done" );
    }

#else
  ERROR ;
#endif
}

// init platform GPU
void init_platform_gpu(bool &GPU_CUBLAS, int gpu_dev)
{
#ifndef NOCUDA

  //  print_platform_gpu_info();

  LOG_INFO("init_platform", "Using GPU device=" << gpu_dev );

  cudaThreadExit();
  cudaSetDevice(gpu_dev);
  Cuda_Error_Check ;

  LOG_INFO("init_platform", "set device -  done");

  if ( GPU_CUBLAS )
    {
      LOG_INFO("init_platform", "init CUBLAS ..." );

      cublasInit(); // init CUDA BLAS
      Cuda_Blas_Error_Check ;

      LOG_INFO("init_platform", "init CUBLAS - done" );
    }

#else
  ERROR ;
#endif
}

// stop platform GPU
void stop_platform_gpu(bool &GPU_CUBLAS)
{
#ifndef NOCUDA

  LOG_INFO("init_platform", "stop platform CUDA/CUBLAS ... ");

  if ( GPU_CUBLAS )
    {
      cublasShutdown(); // shutdown CUDA BLAS
      Cuda_Blas_Error_Check ;
    }

  cudaThreadExit();

  LOG_INFO("init_platform", "stop platform CUDA/CUBLAS ... done");

#else
  ERROR ;
#endif
}
