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

/// @author Nico Trost, Benedikt Galler, Dimitar Lukarski, Martin Wlotzka

#ifndef __LVECTOR_OPENCL_H
#    define __LVECTOR_OPENCL_H

#    include "config.h"

#    include <iostream>
#    include <cstring>

#    include "lvector_cpu.h"

#    ifdef WITH_OPENCL

#        ifdef __APPLE__
#            include <cl.h>
#        else
#            include <CL/cl.h>
#        endif
#        include "opencl/opencl_global.h"

#    endif

using namespace hiflow::la;

/// @brief Provides the base vector class for OPENCL Platform
/// @author Nico Trost, Benedikt Galler, Dimitar Lukarski
///
/// OPENCL_lVector maintains the buffer of the vector;
/// assignment operator (assign vector from cpu and OPENCL)

//namespace hiflow {
//namespace la {

template <typename ValueType>
class OPENCL_lVector : public hiflow::la::lVector<ValueType>
{
  public:
    OPENCL_lVector ( );

#    ifdef WITH_OPENCL
    OPENCL_lVector ( opencl_manager *man );
#    endif
    virtual ~OPENCL_lVector ( );

    virtual void Init ( const int size,
                        const std::string name );
#    ifdef WITH_OPENCL
    virtual void set_opencl_manager ( opencl_manager *man );
#    endif

    virtual void Clear ( void );
    virtual void Zeros ( void );

    virtual void Reorder ( const int *index );

    virtual void CopyFrom ( const lVector<ValueType>& vec );
    virtual void CopyTo ( lVector<ValueType>& vec ) const;

    virtual lVector<ValueType> &operator= ( const lVector<ValueType> &vec );

    virtual lVector<ValueType> *CloneWithoutContent ( ) const;

    virtual void CastFrom ( const hiflow::la::lVector<double>& vec );
    virtual void CastFrom ( const hiflow::la::lVector<float>& vec );

    virtual void CastTo ( hiflow::la::lVector<double>& vec ) const;
    virtual void CastTo ( hiflow::la::lVector<float>& vec ) const;

    virtual lVector<ValueType> *extract_subvector ( const int start_i,
                                                    const int end_i ) const;

    virtual void partial_replace_subvector ( const int start_i,
                                             const int start_sub_vec,
                                             const int size,
                                             const lVector<ValueType> &sub_vec );

    virtual void partial_add_subvector ( const int start_i,
                                         const int start_sub_vec,
                                         const int size,
                                         const ValueType weight,
                                         const lVector<ValueType> &sub_vec );

    virtual void partial_replace_subvector ( const int start_i,
                                             const int start_sub_vec,
                                             const int size,
                                             const OPENCL_lVector<ValueType> &sub_vec );

    virtual void partial_add_subvector ( const int start_i,
                                         const int start_sub_vec,
                                         const int size,
                                         const ValueType weight,
                                         const OPENCL_lVector<ValueType> &sub_vec );

    virtual void Sync ( void ) const;

    virtual void add_value ( const int i, ValueType val );
    virtual void add_values ( const int* indices, int length, const ValueType* values );

    virtual void SetValues ( const int *index,
                             const int size,
                             const ValueType *values );

    virtual void GetValues ( const int *index,
                             const int size,
                             ValueType *values ) const;

#    ifdef WITH_OPENCL
    virtual void set_thread_blocks ( cl_device_id my_device, const int size );
#    endif

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

#    ifdef WITH_OPENCL
    cl_mem buffer;
    opencl_manager *my_manager;
    bool manager_initialized;
#    endif

    virtual void ElementWiseMult ( const lVector<ValueType> &vec ) const;
    virtual void ElementWiseMult ( const OPENCL_lVector<ValueType> &vec ) const;

    virtual ValueType Norm1 ( void ) const;

    virtual ValueType Norm2 ( void ) const;

    virtual ValueType NormMax ( void ) const
    {
        LOG_ERROR ( "OPENCL_lVector<ValueType>::NormMax() not implemented" );
        exit ( -1 );
    }

    virtual ValueType Dot ( const lVector<ValueType> &vec ) const;
    virtual ValueType Dot ( const OPENCL_lVector<ValueType> &vec ) const;

    virtual void Axpy ( const lVector<ValueType> &vec,
                        const ValueType scalar );
    virtual void Axpy ( const OPENCL_lVector<ValueType> &vec,
                        const ValueType scalar );

    virtual void ScaleAdd ( const ValueType scalar,
                            const lVector<ValueType> &vec );
    virtual void ScaleAdd ( const ValueType scalar,
                            const OPENCL_lVector<ValueType> &vec );

    virtual void Scale ( const ValueType scalar );

    virtual void Rot ( lVector<ValueType> *vec,
                       const ValueType &sc, const ValueType &ss );

    virtual void Rot ( OPENCL_lVector<ValueType> *vec,
                       const ValueType &sc,
                       const ValueType &ss );

    virtual void Rotg ( ValueType *sa,
                        ValueType *sb,
                        ValueType *sc,
                        ValueType *ss ) const;

    virtual void Rotm ( lVector<ValueType> *vec,
                        const ValueType &sparam );

    virtual void Rotm ( OPENCL_lVector<ValueType> *vec,
                        const ValueType &sparam );

    virtual void Rotmg ( ValueType *sd1,
                         ValueType *sd2,
                         ValueType *x1,
                         const ValueType &x2,
                         ValueType *sparam ) const;

    virtual int ArgMin ( void ) const;

    virtual int ArgMax ( void ) const;

    template <typename OtherValueType>
    void CastFrom ( const hiflow::la::lVector<OtherValueType>& other );

    template <typename OtherValueType>
    void CastFrom ( const OPENCL_lVector<OtherValueType>& other );

    template <typename OtherValueType>
    void CastFrom ( const CPU_lVector<OtherValueType>& other );

  protected:

    size_t local_threads;
    size_t global_threads;
    size_t global_blocks;

#    ifdef WITH_OPENCL
    cl_mem indexset_;
#    endif
    int indexset_size_;

};

#endif
