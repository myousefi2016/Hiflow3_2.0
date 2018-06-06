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

#ifndef __MEM_H
#    define __MEM_H

#    include <cstring>

/// @brief Cross-Platform Memory managements utility
/// @author Dimitar Lukarski
///
/// Provides platform and cross-platform memory transfer mechanism
/// Including - CPU, GPU

template<typename ValueType>
void memcpy2dev ( ValueType *dest, const ValueType *src, const size_t size );

template<typename ValueType>
void memcpy2host ( ValueType *dest, const ValueType *src, const size_t size );

template<typename ValueType>
void memcpydev ( ValueType *dest, const ValueType *src, const size_t size );

template<typename ValueType>
void memcpyhost ( ValueType *dest, const ValueType *src, const size_t size );

template<typename ValueType>
void memfreedev ( ValueType *p );

template<typename ValueType>
void cudasetvalues ( const int *index, const size_t size, const ValueType *values, ValueType *buffer,
                     int thread_block_size );

template<typename ValueType>
void cudagetvalues ( const int *index, const size_t size, ValueType *values, const ValueType *buffer,
                     int thread_block_size );

template<typename ValueType>
void cudasetblockvalues ( const int start_i, const int start_sub_vec, const int size, const ValueType *values, ValueType *buffer,
                          int thread_block_size );

template<typename ValueType>
void cudagetblockvalues ( const int start_i, const int end_i, ValueType *values, const ValueType *buffer,
                          int thread_block_size );

template<typename ValueType>
void cudaaddblockvalues ( const int start_i, const int start_sub_vec, const int size, const ValueType *values, ValueType *buffer,
                          const ValueType weight, int thread_block_size );

template<typename ValueType>
void cudamultvalues ( const int size, const ValueType *values, ValueType *buffer, int thread_block_size );

template<typename ValueType>
void cudacastfromdouble ( const int, ValueType*, const double*, int );

template<typename ValueType>
void cudacastfromfloat ( const int, ValueType*, const float*, int );

template<typename ValueType>
void cudacasttodouble ( const int, double*, const ValueType*, int );

template<typename ValueType>
void cudacasttofloat ( const int, float*, const ValueType*, int );

template<typename ValueType>
void cudaswapdiagelemtorowfront ( ValueType*, int*, const int*, const int, int );

void cuda_sync_threads ( void );

// no template
// size_type = sizeof(ValueType)
void memallocdev ( void **p, const size_t size, const size_t size_type );

void memsetdev ( void *p, const int value, const size_t size, const size_t size_type );
void memsethost ( void *p, const int value, const size_t size, const size_t size_type );

#endif
