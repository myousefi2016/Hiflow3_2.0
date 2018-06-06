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

/// @author Dimitar Lukarski, Martin Wlotzka

#include "cuda/lmatrix_csr_gpu.h"
#include "../lmp_log.h"
#include "../lmp_mem.h"

#include <assert.h>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

#if (!defined(NOCUDA) && !defined(CUDA_NO_SM_13_DOUBLE_INTRINSICS))
//#ifndef NOCUDA

#include <cuda.h>

#include <cuda_runtime_api.h>

#include "cuda/cuda_utils.h"

#include "cuda/cuda_spmv_kernel_csr_scalar_tex.h"

#else

#define ERROR LOG_ERROR("no Cuda (CUDA_NO_SM_13_DOUBLE_INTRINSICS) support");  exit(-1);

#endif

using namespace hiflow::la;

template <typename ValueType>
GPUscalartex_CSR_lMatrix<ValueType>::GPUscalartex_CSR_lMatrix(int init_nnz, int init_num_row, int init_num_col, std::string init_name)
{
#ifndef NOCUDA
  this->Init(init_nnz, init_num_row, init_num_col, init_name);
  this->implementation_name_ = "Scalar + Tex Cache" ;
  this->implementation_id_ = SCALAR_TEX ;
  // thread_block_size, thread_block are init in the Init();
#else
  ERROR ;
#endif
}

template <typename ValueType>
GPUscalartex_CSR_lMatrix<ValueType>::GPUscalartex_CSR_lMatrix()
{
#ifndef NOCUDA
  this->implementation_name_ = "Scalar + Tex Cache" ;
  this->implementation_id_ = SCALAR_TEX ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
GPUscalartex_CSR_lMatrix<ValueType>::~GPUscalartex_CSR_lMatrix()
{
}

template <typename ValueType>
void GPUscalartex_CSR_lMatrix<ValueType>::VectorMult(const lVector<ValueType> &invec, lVector<ValueType> *outvec) const
{

  const GPU_lVector<ValueType> *casted_invec = dynamic_cast<const GPU_lVector<ValueType>*> (&invec) ;
  GPU_lVector<ValueType> *casted_outvec = dynamic_cast<GPU_lVector<ValueType>*> (outvec) ;

  if ((casted_invec == NULL) && (casted_outvec == NULL))  {
    LOG_ERROR("ERROR GPUscalartex_CSR_lMatrix<ValueType>::VectorMult unsupported in or out vector");
    this->print();
    invec.print();
    outvec->print();
    exit(-1);
  }

  this->VectorMult(*casted_invec, casted_outvec) ;

}

template <typename ValueType>
void GPUscalartex_CSR_lMatrix<ValueType>::VectorMult(const GPU_lVector<ValueType> &invec, GPU_lVector<ValueType> *outvec) const
{
#if (!defined(NOCUDA) && !defined(CUDA_NO_SM_13_DOUBLE_INTRINSICS))
  //#ifndef NOCUDA

  dim3 dimBlock(this->thread_block_size_) ;
  dim3 dimGrid(this->thread_block_);

  assert(invec .get_size() > 0);
  assert(outvec->get_size() > 0);
  assert(invec.get_size()   == this->get_num_col());
  assert(outvec->get_size() == this->get_num_row());

#if 0
  int indices[3] = {0,1,2};
  ValueType values[3];
  invec.GetValues(indices, 3, values);
  std::cout << "---CUDA DEBUG--- " << values[0] << " " << values[1] << " " << values[2] << std::endl;

  ValueType mat[9];
  memcpy2host(mat, this->matrix.val, 9);
  for (int i=0; i < 9; ++i) std::cout << "---CUDA DEBUG--- mat " << i << " " << mat[i] << std::endl;

  int col[9];
  memcpy2host(col, this->matrix.col, 9);
  for (int i=0; i < 9; ++i) std::cout << "---CUDA DEBUG--- col " << i << " " << col[i] << std::endl;

  int row[4];
  memcpy2host(row, this->matrix.row, 4);
  for (int i=0; i < 4; ++i) std::cout << "---CUDA DEBUG--- row " << i << " " << row[i] << std::endl;
#endif

  cudaThreadSynchronize();

  bind2tex(invec.buffer) ;

  kernel_spmv_csr_scalar_tex<<<dimGrid,dimBlock>>>((ValueType *) this->matrix.val, (int *) this->matrix.col , (int *) this->matrix.row,
                                                   (ValueType *) invec.buffer , (ValueType *) outvec->buffer , this->num_row_);

  unbind2tex(invec.buffer) ;

#if 0
  outvec->GetValues(indices, 3, values);
  std::cout << "---CUDA DEBUG--- " << values[0] << " " << values[1] << " " << values[2] << std::endl;
#endif

  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
void GPUscalartex_CSR_lMatrix<ValueType>::VectorMultAdd(const lVector<ValueType> &invec, lVector<ValueType> *outvec) const
{

  const GPU_lVector<ValueType> *casted_invec = dynamic_cast<const GPU_lVector<ValueType>*> (&invec) ;
  GPU_lVector<ValueType> *casted_outvec = dynamic_cast<GPU_lVector<ValueType>*> (outvec) ;

  if ((casted_invec == NULL) && (casted_outvec == NULL))  {
    LOG_ERROR("ERROR GPUscalartex_CSR_lMatrix<ValueType>::VectorMultAdd unsupported in or out vector");
    this->print();
    invec.print();
    outvec->print();
    exit(-1);
  }

  this->VectorMultAdd(*casted_invec, casted_outvec) ;

}

template <>
void GPUscalartex_CSR_lMatrix<float>::VectorMultAdd(const GPU_lVector<float> &invec, GPU_lVector<float> *outvec) const
{
#if (!defined(NOCUDA) && !defined(CUDA_NO_SM_13_DOUBLE_INTRINSICS))
  //#ifndef NOCUDA

  dim3 dimBlock(this->thread_block_size_) ;
  dim3 dimGrid(this->thread_block_);

  assert(invec  .get_size() >= 0);
  assert(outvec->get_size() >= 0);
  assert(invec  .get_size() == this->get_num_col());
  assert((outvec->get_size() == this->get_num_row()) ||
         ( invec .get_size() == 0) );

  if ( this->get_nnz() > 0) {
    bind2tex(invec.buffer) ;

    kernel_add_spmv_csr_scalar_tex<<<dimGrid,dimBlock>>>((float *) this->matrix.val, (int *) this->matrix.col , (int *) this->matrix.row,
                                                         (float *) invec.buffer , (float *) outvec->buffer , this->num_row_);

    unbind2tex(invec.buffer) ;
  }

  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template <>
void GPUscalartex_CSR_lMatrix<double>::VectorMultAdd(const GPU_lVector<double> &invec, GPU_lVector<double> *outvec) const
{
#if (!defined(NOCUDA) && !defined(CUDA_NO_SM_13_DOUBLE_INTRINSICS))
  //#ifndef NOCUDA

  dim3 dimBlock(this->thread_block_size_) ;
  dim3 dimGrid(this->thread_block_);

  assert(invec  .get_size() >= 0);
  assert(outvec->get_size() >= 0);
  assert(invec  .get_size() == this->get_num_col());
  assert((outvec->get_size() == this->get_num_row()) ||
         ( invec .get_size() == 0) );

  if ( this->get_nnz() > 0) {
    bind2tex(invec.buffer) ;

    kernel_add_spmv_csr_scalar_tex<<<dimGrid,dimBlock>>>((double *) this->matrix.val, (int *) this->matrix.col , (int *) this->matrix.row,
                                                         (double *) invec.buffer , (double *) outvec->buffer , this->num_row_);

    unbind2tex(invec.buffer) ;
  }

  Cuda_Error_Check ;
#else
  ERROR ;
#endif
}

template <>
void GPUscalartex_CSR_lMatrix<float>::VectorMultNoDiag(const GPU_lVector<float> &in,
                                                           GPU_lVector<float> *out) const
{
  assert(out != 0);
  assert(this->get_num_row() == out->get_size());
  assert(this->get_num_col() == in.get_size());

#if (!defined(NOCUDA) && !defined(CUDA_NO_SM_13_DOUBLE_INTRINSICS))

  dim3 dimBlock(this->thread_block_size_);
  dim3 dimGrid (this->thread_block_);

  cudaThreadSynchronize();

  bind2tex(in.buffer);

  kernel_spmvnd_csr_scalar_tex<<<dimGrid,dimBlock>>>(this->matrix.val, this->matrix.col, this->matrix.row,
                                                     in.buffer, out->buffer, this->get_num_row());

  unbind2tex(in.buffer);

#else
  ERROR ;
#endif
}

template <>
void GPUscalartex_CSR_lMatrix<double>::VectorMultNoDiag(const GPU_lVector<double> &in,
                                                           GPU_lVector<double> *out) const
{
  assert(out != 0);
  assert(this->get_num_row() == out->get_size());
  assert(this->get_num_col() == in.get_size());

#if (!defined(NOCUDA) && !defined(CUDA_NO_SM_13_DOUBLE_INTRINSICS))

  dim3 dimBlock(this->thread_block_size_);
  dim3 dimGrid (this->thread_block_);

  cudaThreadSynchronize();

  bind2tex(in.buffer);

  kernel_spmvnd_csr_scalar_tex<<<dimGrid,dimBlock>>>(this->matrix.val, this->matrix.col, this->matrix.row,
                                                     in.buffer, out->buffer, this->get_num_row());

  unbind2tex(in.buffer);

#else
  ERROR ;
#endif
}

template <>
void GPUscalartex_CSR_lMatrix<float>::VectorMultNoDiag(const hiflow::la::lVector<float> &in,
                                                           hiflow::la::lVector<float> *out) const
{
  assert(out != 0);
  assert(this->get_num_row() == out->get_size());
  assert(this->get_num_col() == in.get_size());

  const GPU_lVector<float> *in_gpu  = dynamic_cast<const GPU_lVector<float>*>(&in);
        GPU_lVector<float> *out_gpu = dynamic_cast<      GPU_lVector<float>*>(out);

  if ((in_gpu != 0) && (out_gpu != 0))
  {
    this->VectorMultNoDiag(*in_gpu, out_gpu);
  }
  else
  {
    LOG_ERROR("GPUscalartex_CSR_lMatrix<float>::VectorMultNoDiag called with non-GPU vector argument.");
    this->print();
    in.print();
    out->print();
  }
}

template <>
void GPUscalartex_CSR_lMatrix<double>::VectorMultNoDiag(const hiflow::la::lVector<double> &in,
                                                           hiflow::la::lVector<double> *out) const
{
  assert(out != 0);
  assert(this->get_num_row() == out->get_size());
  assert(this->get_num_col() == in.get_size());

  const GPU_lVector<double> *in_gpu  = dynamic_cast<const GPU_lVector<double>*>(&in);
        GPU_lVector<double> *out_gpu = dynamic_cast<      GPU_lVector<double>*>(out);

  if ((in_gpu != 0) && (out_gpu != 0))
  {
    this->VectorMultNoDiag(*in_gpu, out_gpu);
  }
  else
  {
    LOG_ERROR("GPUscalartex_CSR_lMatrix<double>::VectorMultNoDiag called with non-GPU vector argument.");
    this->print();
    in.print();
    out->print();
  }
}

template class GPUscalartex_CSR_lMatrix<float>;
template class GPUscalartex_CSR_lMatrix<double>;
