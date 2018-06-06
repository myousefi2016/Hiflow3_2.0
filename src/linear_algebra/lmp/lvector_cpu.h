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

#ifndef __LVECTOR_CPU_H
#    define __LVECTOR_CPU_H

#    include <iostream>
#    include <cstring>

#    include "lvector.h"

/// @brief Provides the base vector class for CPU
/// @author Dimitar Lukarski
///
/// CPU_lVector maintains the buffer of the vector;
/// provides accessing functions and the
/// assignment operator (assign vector from cpu and gpu)

template <typename ValueType>
class CPU_lVector : public hiflow::la::lVector<ValueType>
{
  public:
    CPU_lVector ( );
    virtual ~CPU_lVector ( );

    virtual void Init ( const int size,
                        const std::string name );
    virtual void Clear ( void );
    virtual void Zeros ( void );

    virtual void ReadFile ( const char* filename );
    virtual void WriteFile ( const char* filename );

    virtual void Reorder ( const int *index );

    virtual hiflow::la::lVector<ValueType> &operator= ( const hiflow::la::lVector<ValueType> &vec2 );
    virtual void CopyFrom ( const hiflow::la::lVector<ValueType>& vec );

    virtual void CopyTo ( hiflow::la::lVector<ValueType>& vec ) const;

    virtual void CastFrom ( const hiflow::la::lVector<double>& vec );
    virtual void CastFrom ( const hiflow::la::lVector<float>& vec );

    virtual void CastFrom ( const CPU_lVector<double>& vec ) = 0;
    virtual void CastFrom ( const CPU_lVector<float>& vec ) = 0;

    virtual void CastTo ( hiflow::la::lVector<double>& vec ) const;
    virtual void CastTo ( hiflow::la::lVector<float>& vec ) const;

    virtual void CastTo ( CPU_lVector<double>& vec ) const = 0;
    virtual void CastTo ( CPU_lVector<float>& vec ) const = 0;

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
                                             const CPU_lVector<ValueType> &sub_vec );

    virtual void partial_add_subvector ( const int start_i,
                                         const int start_sub_vec,
                                         const int size,
                                         const ValueType weight,
                                         const CPU_lVector<ValueType> &sub_vec );

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

    virtual void set_indexset ( const int *indexset,
                                const int size );

    virtual void get_indexset ( int *indexset ) const;

    virtual void GetIndexedValues ( ValueType *values ) const;
    virtual void SetIndexedValues ( const ValueType *values );

    virtual int get_indexset_size ( void ) const;

    ValueType *buffer;

    virtual void ElementWiseMult ( const hiflow::la::lVector<ValueType> &vec ) const;
    virtual void ElementWiseMult ( const CPU_lVector<ValueType> &vec ) const;

    virtual ValueType* GetBuffer ( ) const
    {
        return this->buffer;
    }

  protected:
    int *indexset_;
    int indexset_size_;

};

/// @brief Provides CPU naive/simple sequential implementation
/// of the Blas 1 routines
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUsimple_lVector : public CPU_lVector<ValueType>
{
  public:

    CPUsimple_lVector ( );
    CPUsimple_lVector ( const int size,
                        const std::string name );
    virtual ~CPUsimple_lVector ( );

    virtual int ArgMin ( void ) const;
    virtual int ArgMax ( void ) const;
    virtual ValueType Norm1 ( void ) const;
    virtual ValueType Norm2 ( void ) const;
    virtual ValueType NormMax ( void ) const;
    virtual ValueType Dot ( const hiflow::la::lVector<ValueType> &vec ) const;
    virtual ValueType Dot ( const CPU_lVector<ValueType> &vec ) const;
    virtual void Axpy ( const hiflow::la::lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void Axpy ( const CPU_lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void ScaleAdd ( const ValueType scalar,
                            const hiflow::la::lVector<ValueType> &vec );
    virtual void ScaleAdd ( const ValueType scalar,
                            const CPU_lVector<ValueType> &vec );
    virtual void Scale ( const ValueType scalar );
    virtual void Rot ( hiflow::la::lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rot ( CPU_lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rotg ( ValueType *sa,
                        ValueType *sb,
                        ValueType *sc,
                        ValueType *ss ) const;
    virtual void Rotm ( hiflow::la::lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotm ( CPU_lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotmg ( ValueType *sd1,
                         ValueType *sd2,
                         ValueType *x1,
                         const ValueType &x2,
                         ValueType *sparam ) const;

    virtual void CastFrom ( const CPU_lVector<double>& vec );
    virtual void CastFrom ( const CPU_lVector<float>& vec );

    virtual void CastTo ( CPU_lVector<double>& vec ) const;
    virtual void CastTo ( CPU_lVector<float>& vec ) const;
};

/// @brief Provides CPU OpenMP parallel implementation
/// of the Blas 1 routines
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUopenmp_lVector : public CPU_lVector<ValueType>
{
  public:

    CPUopenmp_lVector ( );
    CPUopenmp_lVector ( const int size,
                        const std::string name );
    virtual ~CPUopenmp_lVector ( );

    virtual void CloneFrom ( const hiflow::la::lVector<ValueType>& vec );

    virtual void ElementWiseMult ( const CPU_lVector<ValueType> &vec ) const;

    virtual int ArgMin ( void ) const;
    virtual int ArgMax ( void ) const;
    virtual ValueType Norm1 ( void ) const;
    virtual ValueType Norm2 ( void ) const;
    virtual ValueType NormMax ( void ) const;
    virtual ValueType Dot ( const hiflow::la::lVector<ValueType> &vec ) const;
    virtual ValueType Dot ( const CPU_lVector<ValueType> &vec ) const;
    virtual void Axpy ( const hiflow::la::lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void Axpy ( const CPU_lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void ScaleAdd ( const ValueType scalar,
                            const hiflow::la::lVector<ValueType> &vec );
    virtual void ScaleAdd ( const ValueType scalar,
                            const CPU_lVector<ValueType> &vec );
    virtual void Scale ( const ValueType scalar );
    virtual void Rot ( hiflow::la::lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rot ( CPU_lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rotg ( ValueType *sa,
                        ValueType *sb,
                        ValueType *sc,
                        ValueType *ss ) const;
    virtual void Rotm ( hiflow::la::lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotm ( CPU_lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotmg ( ValueType *sd1,
                         ValueType *sd2,
                         ValueType *x1,
                         const ValueType &x2,
                         ValueType *sparam ) const;

    virtual void set_num_threads ( void );
    virtual void set_num_threads ( int num_thread );

    int num_threads ( void ) const
    {
        return this->num_threads_;
    }

    virtual void CastFrom ( const CPU_lVector<double>& other );
    virtual void CastFrom ( const CPU_lVector<float>& other );

    virtual void CastTo ( CPU_lVector<double>& vec ) const;
    virtual void CastTo ( CPU_lVector<float>& vec ) const;

  protected:
    int num_threads_;

};

/// @brief Provides wrapper to CPU Intel/MKL implementation
/// of the Blas 1 routines
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUmkl_lVector : public CPU_lVector<ValueType>
{
  public:

    CPUmkl_lVector ( );
    CPUmkl_lVector ( const int size,
                     const std::string name );
    virtual ~CPUmkl_lVector ( );

    virtual void CloneFrom ( const hiflow::la::lVector<ValueType>& vec );

    virtual int ArgMin ( void ) const;
    virtual int ArgMax ( void ) const;
    virtual ValueType Norm1 ( void ) const;
    virtual ValueType Norm2 ( void ) const;
    virtual ValueType NormMax ( void ) const;
    virtual ValueType Dot ( const hiflow::la::lVector<ValueType> &vec ) const;
    virtual ValueType Dot ( const CPU_lVector<ValueType> &vec ) const;
    virtual void Axpy ( const hiflow::la::lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void Axpy ( const CPU_lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void ScaleAdd ( const ValueType scalar,
                            const hiflow::la::lVector<ValueType> &vec );
    virtual void ScaleAdd ( const ValueType scalar,
                            const CPU_lVector<ValueType> &vec );
    virtual void Scale ( const ValueType scalar );
    virtual void Rot ( hiflow::la::lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rot ( CPU_lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rotg ( ValueType *sa,
                        ValueType *sb,
                        ValueType *sc,
                        ValueType *ss ) const;
    virtual void Rotm ( hiflow::la::lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotm ( CPU_lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotmg ( ValueType *sd1,
                         ValueType *sd2,
                         ValueType *x1,
                         const ValueType &x2,
                         ValueType *sparam ) const;

    virtual void set_num_threads ( void );
    virtual void set_num_threads ( int num_thread );

    int num_threads ( void ) const
    {
        return this->num_threads_;
    }

    virtual void CastFrom ( const CPU_lVector<double>& other );
    virtual void CastFrom ( const CPU_lVector<float>& other );

    virtual void CastTo ( CPU_lVector<double>& vec ) const;
    virtual void CastTo ( CPU_lVector<float>& vec ) const;

  protected:
    int num_threads_;
};

/// @brief Provides wrapper to CPU Cblas/Cblas implementation
/// of the Blas 1 routines
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUcblas_lVector : public CPU_lVector<ValueType>
{
  public:

    CPUcblas_lVector ( );
    CPUcblas_lVector ( const int size,
                       const std::string name );

    virtual ~CPUcblas_lVector ( );

    virtual int ArgMin ( void ) const;
    virtual int ArgMax ( void ) const;
    virtual ValueType Norm1 ( void ) const;
    virtual ValueType Norm2 ( void ) const;
    virtual ValueType NormMax ( void ) const;
    virtual ValueType Dot ( const hiflow::la::lVector<ValueType> &vec ) const;
    virtual ValueType Dot ( const CPU_lVector<ValueType> &vec ) const;
    virtual void Axpy ( const hiflow::la::lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void Axpy ( const CPU_lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void ScaleAdd ( const ValueType scalar,
                            const hiflow::la::lVector<ValueType> &vec );
    virtual void ScaleAdd ( const ValueType scalar,
                            const CPU_lVector<ValueType> &vec );
    virtual void Scale ( const ValueType scalar );
    virtual void Rot ( hiflow::la::lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rot ( CPU_lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );
    virtual void Rotg ( ValueType *sa,
                        ValueType *sb,
                        ValueType *sc,
                        ValueType *ss ) const;
    virtual void Rotm ( hiflow::la::lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotm ( CPU_lVector<ValueType> *vec,
                        const ValueType &sparam );
    virtual void Rotmg ( ValueType *sd1,
                         ValueType *sd2,
                         ValueType *x1,
                         const ValueType &x2,
                         ValueType *sparam ) const;

    virtual void CastFrom ( const CPU_lVector<double>& other );
    virtual void CastFrom ( const CPU_lVector<float>& other );

    virtual void CastTo ( CPU_lVector<double>& vec ) const;
    virtual void CastTo ( CPU_lVector<float>& vec ) const;

};

#endif
