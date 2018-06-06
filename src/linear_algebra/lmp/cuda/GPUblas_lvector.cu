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

#include "../lvector.h"
#include "lvector_gpu.h"
#include "cuda_dot_kernel.h"

#include <assert.h>
#include <iostream>
#include <stdlib.h>

#include "../lmp_log.h"

#ifndef NOCUDA

#include <cuda.h>

#include <cublas.h>

#include "cuda/cuda_utils.h"

#else

#define ERROR LOG_ERROR("no Cuda support");  exit(-1);

#endif

using namespace hiflow::la;

template <typename ValueType>
GPUblas_lVector<ValueType>::GPUblas_lVector(const int size, const std::string name)
{
#ifndef NOCUDA
  this->Init(size, name) ;
  this->implementation_name_ = "CUDA/BLAS";
  this->implementation_id_ = BLAS ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
GPUblas_lVector<ValueType>::GPUblas_lVector()
{
#ifndef NOCUDA
  this->implementation_name_ = "CUDA/BLAS";
  this->implementation_id_ = BLAS ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
GPUblas_lVector<ValueType>::~GPUblas_lVector()
{
}

template <>
int GPUblas_lVector<double>::ArgMin(void) const
{
#ifndef NOCUDA
  int minimum = -1 ;

  assert(this->get_size() > 0);

  minimum = cublasIdamin (this->get_size(),  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

  return minimum ;
#else
  ERROR ;
#endif
}

template <>
int GPUblas_lVector<float>::ArgMin(void) const
{
#ifndef NOCUDA
  int minimum = -1 ;

  assert(this->get_size() > 0);

  minimum =  cublasIsamin (this->get_size(),  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

  return minimum ;
#else
  ERROR ;
#endif
}

template <>
int GPUblas_lVector<double>::ArgMax(void) const
{
#ifndef NOCUDA
   int maximum = -1 ;

   assert(this->get_size() > 0);

   maximum = cublasIdamax (this->get_size(),  this->buffer, 1) ;
   Cuda_Blas_Error_Check ;

   return maximum ;
#else
  ERROR ;
#endif
}

template <>
int GPUblas_lVector<float>::ArgMax(void) const
{
#ifndef NOCUDA
   int maximum = -1 ;

   assert(this->get_size() > 0);

   maximum =  cublasIsamax (this->get_size(),  this->buffer, 1) ;
   Cuda_Blas_Error_Check ;

   return maximum ;
#else
  ERROR ;
#endif
}

template <>
double GPUblas_lVector<double>::Norm2(void) const
{
#ifndef NOCUDA
  double norm2 = -1 ;

  assert(this->get_size() > 0);

  norm2 =  cublasDnrm2 (this->get_size(),  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

  return norm2 ;
#else
  ERROR ;
#endif
}

template <>
float GPUblas_lVector<float>::Norm2(void) const
{
#ifndef NOCUDA
  float norm2 = -1 ;

  assert(this->get_size() > 0);

  norm2 =  cublasSnrm2 (this->get_size(),  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

  return norm2 ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
ValueType GPUblas_lVector<ValueType>::NormMax(void) const
{
  LOG_ERROR("GPUblas_lVector<ValueType>::NormMax not implemented yet.");
  this->print();
  exit(-1);
  return 0.0;
}

template <typename ValueType>
ValueType GPUblas_lVector<ValueType>::Dot(const lVector<ValueType> &vec) const
{
  ValueType dot ;

  const GPU_lVector<ValueType> *casted_vec;

  if (casted_vec = dynamic_cast<const GPU_lVector<ValueType>*> (&vec)) {

    dot = this->Dot(*casted_vec) ;

  } else {
    LOG_ERROR("ERROR GPUblas_lVector::Dot unsupported vectors");
    this->print();
    vec.print();
    exit(-1);
 }

  return dot ;
}

template <>
double GPUblas_lVector<double>::Dot(const GPU_lVector<double> &vec) const
{
#ifndef NOCUDA
  double dot ;

  assert(this->get_size() > 0);
  assert(this->get_size() == vec.get_size()) ;

  dot =  cublasDdot (this->get_size(),  this->buffer, 1,
                     vec.buffer, 1) ;
  Cuda_Blas_Error_Check ;

  return dot ;
#else
  ERROR ;
#endif
}

#if 0
template <>
double GPUblas_lVector<double>::Dot(const GPU_lVector<double> &vec) const
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec.get_size());

  static const int block_size = 256;
  dim3 dimBlock(block_size);
  dim3 dimGrid(1);

  double *d_dot;

  cudaMalloc((void **) &d_dot, sizeof(double));
  Cuda_Error_Check;

  dot_kernel<block_size><<<dimGrid, dimBlock>>>(this->buffer, vec.buffer, d_dot, this->get_size());
  Cuda_Error_Check;

  double dot;
  cudaMemcpy(&dot, d_dot, sizeof(double), cudaMemcpyDeviceToHost);
  Cuda_Error_Check;

  cudaFree(d_dot);
  Cuda_Error_Check;

  return dot;
#else
  ERROR;
#endif
}
#endif

template <>
float GPUblas_lVector<float>::Dot(const GPU_lVector<float> &vec) const
{
#ifndef NOCUDA
  float dot ;

  assert(this->get_size() > 0);
  assert(this->get_size() == vec.get_size()) ;

  dot =  cublasSdot (this->get_size(),  this->buffer, 1,
                     vec.buffer, 1) ;
  Cuda_Blas_Error_Check ;

  return dot ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
void GPUblas_lVector<ValueType>::Axpy(const lVector<ValueType> &vec, const ValueType scalar)
{

  const GPU_lVector<ValueType> *casted_vec;

  if (casted_vec = dynamic_cast<const GPU_lVector<ValueType>*> (&vec)) {

    this->Axpy(*casted_vec, scalar) ;

  } else {
    LOG_ERROR("ERROR GPUblas_lVector::Axpy unsupported vectors");
    this->print();
    vec.print();
    exit(-1);
 }

}

template <>
void GPUblas_lVector<double>::Axpy(const GPU_lVector<double> &vec, const double scalar)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec.get_size());

  cublasDaxpy (this->get_size(),  scalar,  vec.buffer, 1,
               this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<float>::Axpy(const GPU_lVector<float> &vec, const float scalar)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec.get_size());

  cublasSaxpy (this->get_size(),  scalar,  vec.buffer, 1,
               this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <typename ValueType>
void GPUblas_lVector<ValueType>::ScaleAdd(const ValueType scalar, const lVector<ValueType> &vec)
{

  const GPU_lVector<ValueType> *casted_vec;

  if (casted_vec = dynamic_cast<const GPU_lVector<ValueType>*> (&vec)) {

    this->ScaleAdd(scalar, *casted_vec) ;

  } else {
    LOG_ERROR("ERROR GPUblas_lVector::ScaleAdd unsupported vectors");
    this->print();
    vec.print();
    exit(-1);
  }

}

template <>
void GPUblas_lVector<double>::ScaleAdd(const double scalar, const GPU_lVector<double> &vec)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec.get_size());

  // buffer = scalar*buffer ;
  cublasDscal (this->get_size(),  scalar,  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

  // buffer += vec.buffer ;
  cublasDaxpy (this->get_size(),  (double) (1.0),  vec.buffer, 1,
               this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<float>::ScaleAdd(const float scalar, const GPU_lVector<float> &vec)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec.get_size());

  // buffer = scalar*buffer ;
  cublasSscal (this->get_size(),  scalar,  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

  // buffer += vec.buffer ;
  cublasSaxpy (this->get_size(),  (float) (1.0),  vec.buffer, 1,
               this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<double>::Scale(const double scalar)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);

  cublasDscal (this->get_size(), scalar,  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<float>::Scale(const float scalar)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);

  cublasSscal (this->get_size(), scalar,  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
double GPUblas_lVector<double>::Norm1(void) const
{
#ifndef NOCUDA
  double sum = -1 ;

  assert(this->get_size() > 0);

  sum =  cublasDasum (this->get_size(),  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

  return sum ;
#else
  ERROR ;
#endif
}

template <>
float GPUblas_lVector<float>::Norm1(void) const
{
#ifndef NOCUDA
  float sum = -1 ;

  assert(this->get_size() > 0);

  sum =  cublasSasum (this->get_size(),  this->buffer, 1) ;
  Cuda_Blas_Error_Check ;

  return sum ;
#else
  ERROR ;
#endif
}

template <typename ValueType>
void GPUblas_lVector<ValueType>::Rot(lVector<ValueType> *vec, const ValueType &sc, const ValueType &ss)
{

  GPU_lVector<ValueType> *casted_vec;

  if (casted_vec = dynamic_cast<GPU_lVector<ValueType>*> (vec)) {

    this->Rot(casted_vec, sc, ss) ;

  } else {
    LOG_ERROR("ERROR GPUblas_lVector::Rot unsupported vectors");
    this->print();
    vec->print();
    exit(-1);
 }

}

template <>
void GPUblas_lVector<double>::Rot(GPU_lVector<double> *vec, const double &sc, const double &ss)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec->get_size());

  cublasDrot (this->get_size(),  vec->buffer, 1, this->buffer, 1, sc, ss);
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<float>::Rot(GPU_lVector<float> *vec, const float &sc, const float &ss)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec->get_size());

  cublasSrot (this->get_size(),  vec->buffer, 1, this->buffer, 1, sc, ss);
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<double>::Rotg(double *sa, double *sb , double *sc, double *ss) const
{
#ifndef NOCUDA

  cublasDrotg (sa, sb, sc, ss);
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<float>::Rotg(float *sa, float *sb , float *sc, float *ss) const
{
#ifndef NOCUDA

  cublasSrotg (sa, sb, sc, ss);
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <typename ValueType>
void GPUblas_lVector<ValueType>::Rotm(lVector<ValueType> *vec, const ValueType &sparam)
{

  GPU_lVector<ValueType> *casted_vec;

  if (casted_vec = dynamic_cast<GPU_lVector<ValueType>*> (vec)) {

    this->Rotm(casted_vec, sparam) ;

  } else {
    LOG_ERROR("ERROR GPUblas_lVector::Dotm unsupported vectors");
    this->print();
    vec->print();
    exit(-1);
 }

}

template <>
void GPUblas_lVector<double>::Rotm(GPU_lVector<double> *vec, const double &sparam)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec->get_size());

  cublasDrotm (this->get_size(), vec->buffer, 1, this->buffer, 1, &sparam) ;
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<float>::Rotm(GPU_lVector<float> *vec, const float &sparam)
{
#ifndef NOCUDA

  assert(this->get_size() > 0);
  assert(this->get_size() == vec->get_size());

  cublasSrotm (this->get_size(), vec->buffer, 1, this->buffer, 1, &sparam) ;
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<double>::Rotmg(double *sd1, double *sd2, double *x1, const double &x2, double *sparam) const
{
#ifndef NOCUDA

  cublasDrotmg (sd1, sd2, x1, &x2, sparam);
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

template <>
void GPUblas_lVector<float>::Rotmg(float *sd1, float *sd2, float *x1, const float &x2, float *sparam) const
{
#ifndef NOCUDA

  cublasSrotmg (sd1, sd2, x1, &x2, sparam);
  Cuda_Blas_Error_Check ;

#else
  ERROR ;
#endif
}

//
//
//
//

template class GPUblas_lVector<double>;
template class GPUblas_lVector<float>;
