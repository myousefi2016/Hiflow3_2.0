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

#ifndef __CUDA_UTILS_H
#    define __CUDA_UTILS_H

#    include <iostream>
#    include "../lmp_log.h"

// Error handler for regular CUDA kernels
#    define Cuda_Error_Check {\
    cudaError_t error; \
    if ( (error = cudaGetLastError() ) != cudaSuccess) {\
      LOG_ERROR("CUDA ERROR:" << cudaGetErrorString(error) << " , file:" << __FILE__ << " , line:" <<  __LINE__);  \
      exit(-1);\
    }\
  }

// Error handler for CUDA BLAS routines
#    define Cuda_Blas_Error_Check {\
    cublasStatus error;\
    if ( (error = cublasGetError() ) != CUBLAS_STATUS_SUCCESS) {\
      LOG_ERROR("CUBLAS ERROR:" << error << " , file:" << __FILE__ << " , line:" <<  __LINE__);  \
      exit(-1);\
    }\
  }

// Macro checking for errors in CUDA calls
#    define CUDA_SAFE_CALL( cuda_function_call ) {\
    cudaError_t error = cuda_function_call; \
    if (error != cudaSuccess) {\
      LOG_ERROR("CUDA ERROR:" << cudaGetErrorString(error) << " , file: " << __FILE__ << " , line: " << __LINE__); \
      exit(-1);\
    }\
  }

#endif
