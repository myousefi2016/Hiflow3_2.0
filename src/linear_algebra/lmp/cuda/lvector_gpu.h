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

#ifndef __LVECTOR_GPU_H
#    define __LVECTOR_GPU_H

#    include <iostream>
#    include <cstring>

#    include "../lvector_cpu.h"

/// @brief Provides the base vector class for GPU/Cuda Platform
/// @author Dimitar Lukarski
///
/// GPU_lVector maintains the vector;
/// access functions and copy operators (to and from a cpu vector)

template <typename ValueType>
class GPU_lVector : public hiflow::la::lVector<ValueType>
{
  public:
    GPU_lVector ( );
    virtual ~GPU_lVector ( );

    virtual void Init ( const int size,
                        const std::string name );
    virtual void Clear ( void );

    virtual void Zeros ( void );

    virtual void Reorder ( const int *index );

    virtual hiflow::la::lVector<ValueType> &operator= ( const hiflow::la::lVector<ValueType> &vec2 );
    virtual void CopyFrom ( const hiflow::la::lVector<ValueType>& vec );
    virtual void CopyTo ( hiflow::la::lVector<ValueType>& vec ) const;

    virtual void CastFrom ( const hiflow::la::lVector<double>& vec );
    virtual void CastFrom ( const hiflow::la::lVector<float>& vec );

    virtual void CastFrom ( const GPU_lVector<double>& other );
    virtual void CastFrom ( const GPU_lVector<float>& other );

    virtual void CastFrom ( const CPU_lVector<double>& other );
    virtual void CastFrom ( const CPU_lVector<float>& other );

    virtual void CastTo ( hiflow::la::lVector<double>& vec ) const;
    virtual void CastTo ( hiflow::la::lVector<float>& vec ) const;

    virtual void CastTo ( GPU_lVector<double>& vec ) const;
    virtual void CastTo ( GPU_lVector<float>& vec ) const;

    virtual void CastTo ( CPU_lVector<double>& vec ) const;
    virtual void CastTo ( CPU_lVector<float>& vec ) const;

    virtual hiflow::la::lVector<ValueType> *extract_subvector ( const int start_i,
                                                                const int end_i ) const;

    virtual void partial_replace_subvector ( const int start_i,
                                             const int start_sub_vec,
                                             const int size,
                                             const hiflow::la::lVector<ValueType> &sub_vec );

    virtual void partial_add_subvector ( const int start_i,
                                         const int start_sub_vec,
                                         const int size,
                                         const ValueType weight,
                                         const hiflow::la::lVector<ValueType> &sub_vec );

    virtual void partial_replace_subvector ( const int start_i,
                                             const int start_sub_vec,
                                             const int size,
                                             const GPU_lVector<ValueType> &sub_vec );

    virtual void partial_add_subvector ( const int start_i,
                                         const int start_sub_vec,
                                         const int size,
                                         const ValueType weight,
                                         const GPU_lVector<ValueType> &sub_vec );

    virtual void Sync ( void ) const;

    virtual void add_value ( const int i, ValueType val );
    virtual void add_values ( const int* indices, int length, const ValueType* values );

    virtual void SetValues ( const int *index,
                             const int size,
                             const ValueType *values );
    virtual void GetValues ( const int *index,
                             const int size,
                             ValueType *values ) const;

    virtual void GetBlockValues ( const int start_i,
                                  const int end_i,
                                  ValueType *values ) const;

    virtual void SetBlockValues ( const int start_i,
                                  const int end_i,
                                  const ValueType *values );

    virtual void set_thread_blocks ( void );
    virtual void set_thread_blocks ( int thread_block, int thread_block_size );

    virtual void set_indexset ( const int *indexset,
                                const int size );

    virtual void get_indexset ( int *indexset ) const;

    virtual void GetIndexedValues ( ValueType *values ) const;
    virtual void SetIndexedValues ( const ValueType *values );
    virtual int get_indexset_size ( void ) const;

    ValueType *buffer;

    virtual ValueType* GetBuffer ( ) const
    {
        return this->buffer;
    }

    virtual void ElementWiseMult ( const hiflow::la::lVector<ValueType> &vec ) const;
    virtual void ElementWiseMult ( const GPU_lVector<ValueType> &vec ) const;

  protected:
    int *indexset_;
    int indexset_size_;

    int thread_block_size_;
    int thread_block_;

};

/// @brief Provides wrapper to GPU/CUBLAS (NVIDIA) implementation
/// of the Blas 1 routines
/// @author Dimitar Lukarski

template <typename ValueType>
class GPUblas_lVector : public GPU_lVector<ValueType>
{
  public:

    GPUblas_lVector ( );
    GPUblas_lVector ( const int size,
                      const std::string name );
    virtual ~GPUblas_lVector ( );

    virtual int ArgMin ( void ) const;
    virtual int ArgMax ( void ) const;
    virtual ValueType Norm1 ( void ) const;
    virtual ValueType Norm2 ( void ) const;
    virtual ValueType NormMax ( void ) const;
    virtual ValueType Dot ( const hiflow::la::lVector<ValueType> &vec ) const;
    virtual ValueType Dot ( const GPU_lVector<ValueType> &vec ) const;
    virtual void Axpy ( const hiflow::la::lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void Axpy ( const GPU_lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void ScaleAdd ( const ValueType scalar,
                            const hiflow::la::lVector<ValueType> &vec );
    virtual void ScaleAdd ( const ValueType scalar,
                            const GPU_lVector<ValueType> &vec );
    virtual void Scale ( const ValueType scalar );
    virtual void Rot ( hiflow::la::lVector<ValueType> *vec,
                       const ValueType &sc, const ValueType &ss );
    virtual void Rot ( GPU_lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rotg ( ValueType *sa,
                        ValueType *sb,
                        ValueType *sc,
                        ValueType *ss ) const;
    virtual void Rotm ( hiflow::la::lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotm ( GPU_lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotmg ( ValueType *sd1,
                         ValueType *sd2,
                         ValueType *x1,
                         const ValueType &x2,
                         ValueType *sparam ) const;

};

#endif
