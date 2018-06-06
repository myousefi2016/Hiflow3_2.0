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

#include <assert.h>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

#ifndef NOCUDA

#include <cuda.h>

#include <cuda_runtime_api.h>

#include "cuda/cuda_utils.h"
#include "cuda/cuda_spmv_kernel_csr_scalar.h"

#else

#define ERROR LOG_ERROR("no Cuda support");  exit(-1);

#endif

using namespace hiflow::la;

template <typename ValueType>
GPUscalar_CSR_lMatrix<ValueType>::GPUscalar_CSR_lMatrix(int init_nnz, int init_num_row, int init_num_col, std::string init_name)
{
#ifndef NOCUDA
  this->Init(init_nnz, init_num_row, init_num_col, init_name);
  this->implementation_name_ = "Scalar" ;
  this->implementation_id_ = SCALAR ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
GPUscalar_CSR_lMatrix<ValueType>::GPUscalar_CSR_lMatrix()
{
#ifndef NOCUDA
  this->implementation_name_ = "Scalar" ;
  this->implementation_id_ = SCALAR ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
GPUscalar_CSR_lMatrix<ValueType>::~GPUscalar_CSR_lMatrix()
{
}

template <typename ValueType>
void GPUscalar_CSR_lMatrix<ValueType>::VectorMult(const lVector<ValueType> &invec, lVector<ValueType> *outvec) const
{

  const GPU_lVector<ValueType> *casted_invec = dynamic_cast<const GPU_lVector<ValueType>*> (&invec) ;
  GPU_lVector<ValueType> *casted_outvec = dynamic_cast<GPU_lVector<ValueType>*> (outvec) ;

  if ((casted_invec == NULL) && (casted_outvec == NULL)) {
    LOG_ERROR("ERROR GPUscalar_CSR_lMatrix<ValueType>::VectorMult unsupported in or out vector");
    this->print();
    invec.print();
    outvec->print();
    exit(-1);
  }

  this->VectorMult(*casted_invec, casted_outvec) ;
}

template <typename ValueType>
void GPUscalar_CSR_lMatrix<ValueType>::VectorMult(const GPU_lVector<ValueType> &invec, GPU_lVector<ValueType> *outvec) const
{
#ifndef NOCUDA

  dim3 dimBlock(this->thread_block_size_) ;
  dim3 dimGrid(this->thread_block_);

  assert(invec .get_size() > 0);
  assert(outvec->get_size() > 0);
  assert(invec. get_size()  == this->get_num_col());
  assert(outvec->get_size() == this->get_num_row());

  cudaThreadSynchronize();

  kernel_spmv_csr_scalar<<<dimGrid,dimBlock>>>((ValueType *) this->matrix.val, (int *) this->matrix.col , (int *) this->matrix.row,
                                               (ValueType *) invec.buffer , (ValueType *) outvec->buffer , this->num_row_);

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template <typename ValueType>
void GPUscalar_CSR_lMatrix<ValueType>::VectorMultAdd(const lVector<ValueType> &invec, lVector<ValueType> *outvec) const
{

  const GPU_lVector<ValueType> *casted_invec = dynamic_cast<const GPU_lVector<ValueType>*> (&invec) ;
  GPU_lVector<ValueType> *casted_outvec = dynamic_cast<GPU_lVector<ValueType>*> (outvec) ;

  if ((casted_invec == NULL) && (casted_outvec == NULL)) {
    LOG_ERROR("ERROR GPUscalar_CSR_lMatrix<ValueType>::VectorMultAdd unsupported in or out vector");
    this->print();
    invec.print();
    outvec->print();
    exit(-1);
  }

  this->VectorMultAdd(*casted_invec, casted_outvec) ;
}

template <>
void GPUscalar_CSR_lMatrix<float>::VectorMultAdd(const GPU_lVector<float> &invec, GPU_lVector<float> *outvec) const
{
#ifndef NOCUDA

  dim3 dimBlock(this->thread_block_size_) ;
  dim3 dimGrid(this->thread_block_);

  assert(invec  .get_size() >= 0);
  assert(outvec->get_size() >= 0);
  assert(invec  .get_size() == this->get_num_col());
  assert((outvec->get_size() == this->get_num_row()) ||
         ( invec .get_size() == 0) );

  if ( this->get_nnz() > 0) {
    kernel_add_spmv_csr_scalar<<<dimGrid,dimBlock>>>((float *) this->matrix.val, (int *) this->matrix.col , (int *) this->matrix.row,
                                                     (float *) invec.buffer , (float *) outvec->buffer , this->num_row_);
  }

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUscalar_CSR_lMatrix<double>::VectorMultAdd(const GPU_lVector<double> &invec, GPU_lVector<double> *outvec) const
{
#ifndef NOCUDA

  dim3 dimBlock(this->thread_block_size_) ;
  dim3 dimGrid(this->thread_block_);

  assert(invec  .get_size() >= 0);
  assert(outvec->get_size() >= 0);
  assert(invec  .get_size() == this->get_num_col());
  assert((outvec->get_size() == this->get_num_row()) ||
         ( invec .get_size() == 0) );

  if ( this->get_nnz() > 0) {
    kernel_add_spmv_csr_scalar<<<dimGrid,dimBlock>>>((double *) this->matrix.val, (int *) this->matrix.col , (int *) this->matrix.row,
                                                     (double *) invec.buffer , (double *) outvec->buffer , this->num_row_);
  }

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUscalar_CSR_lMatrix<float>::VectorMultNoDiag(const GPU_lVector<float> &in,
                                                        GPU_lVector<float> *out) const
{
  assert(out != 0);
  assert(this->get_num_row() == out->get_size());
  assert(this->get_num_col() == in.get_size());

#ifndef NOCUDA

  dim3 dimBlock(this->thread_block_size_);
  dim3 dimGrid (this->thread_block_);

  kernel_spmvnd_csr_scalar<<<dimGrid,dimBlock>>>(this->matrix.val, this->matrix.col, this->matrix.row,
                                                 in.buffer, out->buffer, this->get_num_row());
#else
  ERROR ;
#endif
}

template <>
void GPUscalar_CSR_lMatrix<double>::VectorMultNoDiag(const GPU_lVector<double> &in,
                                                        GPU_lVector<double> *out) const
{
  assert(out != 0);
  assert(this->get_num_row() == out->get_size());
  assert(this->get_num_col() == in.get_size());

#ifndef NOCUDA

  dim3 dimBlock(this->thread_block_size_);
  dim3 dimGrid (this->thread_block_);

  kernel_spmvnd_csr_scalar<<<dimGrid,dimBlock>>>(this->matrix.val, this->matrix.col, this->matrix.row,
                                                 in.buffer, out->buffer, this->get_num_row());
#else
  ERROR ;
#endif
}

template <>
void GPUscalar_CSR_lMatrix<float>::VectorMultNoDiag(const hiflow::la::lVector<float> &in,
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
    LOG_ERROR("GPUscalar_CSR_lMatrix<float>::VectorMultNoDiag called with non-GPU vector argument.");
    this->print();
    in.print();
    out->print();
  }
}

template <>
void GPUscalar_CSR_lMatrix<double>::VectorMultNoDiag(const hiflow::la::lVector<double> &in,
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
    LOG_ERROR("GPUscalar_CSR_lMatrix<double>::VectorMultNoDiag called with non-GPU vector argument.");
    this->print();
    in.print();
    out->print();
  }
}

template <>
void GPUscalar_CSR_lMatrix<float>::BlocksPsgauss_seidel(const lVector<float> &invec,
    lVector<float> *outvec, const int num_blocks) const
{
#ifndef NOCUDA

  const GPU_lVector<float> *casted_invec = dynamic_cast<const GPU_lVector<float>*> (&invec) ;
  GPU_lVector<float> *casted_outvec = dynamic_cast<GPU_lVector<float>*> (outvec) ;

  if ((casted_invec == NULL) && (casted_outvec == NULL)) {
    LOG_ERROR("ERROR GPUscalar_CSR_lMatrix<float>::VectorMultAdd unsupported in or out vector");
    this->print();
    invec.print();
    outvec->print();
    exit(-1);
  }

    // default values
  int thread_block_size = this->thread_block_size_ ;
  int thread_block = num_blocks/thread_block_size + 1 ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  assert(casted_invec ->get_size() >= 0);
  assert(casted_outvec->get_size() >= 0);
  assert(casted_invec ->get_size() == this->get_num_col());
  assert((casted_outvec->get_size() == this->get_num_row()) ||
         ( casted_invec->get_size() == 0) );

  if (casted_invec->get_size()  > 0) {

    int step_size = this->num_row_ / num_blocks;

    cudaThreadSynchronize();

    kernel_BlockPsgauss_seidel<<<dimGrid,dimBlock>>>((float *) this->matrix.val, (int *) this->matrix.col , (int *) this->matrix.row,
                                                     (float *) casted_invec->buffer , (float *) casted_outvec->buffer,
                                                     num_blocks, this->num_row_, step_size);
    cudaThreadSynchronize();

  }

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUscalar_CSR_lMatrix<double>::BlocksPsgauss_seidel(const lVector<double> &invec,
    lVector<double> *outvec, const int num_blocks) const
{
#ifndef NOCUDA

  const GPU_lVector<double> *casted_invec = dynamic_cast<const GPU_lVector<double>*> (&invec) ;
  GPU_lVector<double> *casted_outvec = dynamic_cast<GPU_lVector<double>*> (outvec) ;

  if ((casted_invec == NULL) && (casted_outvec == NULL)) {
    LOG_ERROR("ERROR GPUscalar_CSR_lMatrix<double>::VectorMultAdd unsupported in or out vector");
    this->print();
    invec.print();
    outvec->print();
    exit(-1);
  }

    // default values
  int thread_block_size = this->thread_block_size_ ;
  int thread_block = num_blocks/thread_block_size + 1 ;

  dim3 dimBlock(thread_block_size) ;
  dim3 dimGrid(thread_block);

  assert(casted_invec ->get_size() >= 0);
  assert(casted_outvec->get_size() >= 0);
  assert(casted_invec ->get_size() == this->get_num_col());
  assert((casted_outvec->get_size() == this->get_num_row()) ||
         ( casted_invec->get_size() == 0) );

  if (casted_invec->get_size()  > 0) {

    int step_size = this->num_row_ / num_blocks;

    cudaThreadSynchronize();

    kernel_BlockPsgauss_seidel<<<dimGrid,dimBlock>>>((double *) this->matrix.val, (int *) this->matrix.col , (int *) this->matrix.row,
                                                     (double *) casted_invec->buffer , (double *) casted_outvec->buffer,
                                                     num_blocks, this->num_row_, step_size);
    cudaThreadSynchronize();

  }

  Cuda_Error_Check ;

#else
  ERROR ;
#endif
}

template class GPUscalar_CSR_lMatrix<float>;
template class GPUscalar_CSR_lMatrix<double>;
