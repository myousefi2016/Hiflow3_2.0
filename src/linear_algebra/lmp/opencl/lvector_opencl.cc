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

#include "../init_vec_mat.h"
#include "../lmp_log.h"

#include "opencl/lvector_opencl.h"
#include "opencl/mem_opencl.h"
#include "opencl/opencl_utils.h"
#include "opencl/opencl_global.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifdef WITH_OPENCL

#    ifdef __APPLE__
#        include <cl.h>
#    else
#        include <CL/cl.h>
#    endif
#    include "opencl/mem_opencl.h"
#    include "opencl/opencl_kernel_mapper.h"
#else

#    define ERROR LOG_ERROR("ERROR no OPENCL support");  exit(-1);

#endif

using namespace hiflow::la;

// class OPENCL_lVector

template <typename ValueType>
OPENCL_lVector<ValueType>::OPENCL_lVector ( )
{

    this->platform_name_ = "OPENCL";
    this->platform_id_ = OPENCL;
    this->implementation_name_ = "OpenCL Stream";

#ifdef WITH_OPENCL
    this->manager_initialized = false;
#endif

    this->indexset_size_ = 0;

}

#ifdef WITH_OPENCL

template <typename ValueType>
OPENCL_lVector<ValueType>::OPENCL_lVector ( opencl_manager *man )
{

    this->platform_name_ = "OPENCL";
    this->platform_id_ = OPENCL;
    this->implementation_name_ = "OpenCL Stream";
    this->implementation_id_ = OPEN_CL;
    this->my_manager = new opencl_manager;
    my_manager = man;
    this->manager_initialized = true;

    this->indexset_size_ = 0;

}

template <typename ValueType>
void OPENCL_lVector<ValueType>::set_opencl_manager ( opencl_manager *man )
{

    if ( !manager_initialized ) this->my_manager = new opencl_manager;
    my_manager = man;
    this->manager_initialized = true;

}
#endif

template <typename ValueType>
OPENCL_lVector<ValueType>::~OPENCL_lVector ( )
{
    this->Clear ( );
}

template <typename ValueType>
lVector<ValueType> *OPENCL_lVector<ValueType>::CloneWithoutContent ( ) const
{
    std::string cloned_name;
    cloned_name = "clone from ";
    cloned_name.append ( this->name_ );

    lVector<ValueType> *new_vector;

#ifdef WITH_OPENCL
    new_vector = new OPENCL_lVector<ValueType>( this->my_manager );
#else
    ERROR;
#endif

    new_vector->Init ( this->get_size ( ), cloned_name );
    return new_vector;
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Clear ( )
{
#ifdef WITH_OPENCL

    if ( this->get_size ( ) > 0 )
        opencl_memfreedev ( this->buffer );

    this->size_ = 0;
    this->name_ = "";

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Zeros ( )
{
#ifdef WITH_OPENCL

    int size = this->get_size ( );

    if ( size > 0 )
        openclzeros ( this->buffer, size,
                      my_manager->get_kernel ( 2 ), my_manager->get_command_queue ( ),
                      this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Reorder ( const int *index )
{

    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::Reorder() not implemented yet" );
    this->print ( );
    exit ( -1 );

}

#ifdef WITH_OPENCL

template <typename ValueType>
void OPENCL_lVector<ValueType>::set_thread_blocks ( cl_device_id my_device, const int size )//TODO, so far only for GPU
{

    if ( size > 0 )
    {

        cl_int threads, blocks, globalThreads;
        int maxBlocks = my_manager->maxBlocks;
        int maxThreads = my_manager->maxThreads;
        int kerneltype = my_manager->kerneltype;

        // THREAD management by NVIDIA
        // This software contains source code provided by NVIDIA Corporation.

        if ( kerneltype == 0 )
        {
            // Device = NVIDIA GPU
            threads = ( size < maxThreads * 2 ) ? nextPow2 ( ( size + 1 ) / 2 ) : maxThreads;
            blocks = ( size + ( threads * 2 - 1 ) ) / ( threads * 2 );
            blocks = ( blocks < maxBlocks ) ? blocks : maxBlocks;
        }

        if ( kerneltype == 1 )
        {
            // Device = ATI GPU
            threads = ( size < maxThreads ) ? nextPow2 ( size ) : maxThreads;
            blocks = ( size + threads - 1 ) / threads;
        }

        if ( kerneltype == 2 )
        {
            // DEVICE = CPU
            threads = maxThreads;
            blocks = ( size + threads - 1 ) / threads;
        }

        globalThreads = ( size / threads ) * threads;
        if ( size % threads != 0 ) globalThreads += threads;

        this->local_threads = threads;
        this->global_threads = globalThreads;
        this->global_blocks = blocks;
    }
    else
    {
        LOG_ERROR ( "Invalid size" );
    }

}
#endif

template <typename ValueType>
void OPENCL_lVector<ValueType>::Init ( const int size, const std::string name )
{
#ifdef WITH_OPENCL

    assert ( manager_initialized );
    assert ( size >= 0 );

    if ( this->get_size ( ) > 0 )
        this->Clear ( );

    if ( size > 0 )
        opencl_memallocdev<ValueType>( this->buffer, size, my_manager->get_context ( ) );

    this->set_thread_blocks ( my_manager->get_device ( ), size );
    this->name_ = name;
    this->size_ = size;

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::add_value ( const int i, ValueType val )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::add_value not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::add_values ( const int* indices, int length, const ValueType* values )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::add_values not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::SetValues ( const int *index, const int size, const ValueType *values )
{
#ifdef WITH_OPENCL

    assert ( values != NULL );
    assert ( this->get_size ( ) > 0 );
    assert ( size != 0 );

    for ( int i = 0; i < size; ++i )
    {
        assert ( index[i] >= 0 );
        assert ( index[i] < this->get_size ( ) );
    }

    // copy the index set to the device
    cl_mem index_buffer;
    opencl_memallocdev<int>( index_buffer, size, my_manager->get_context ( ) );
    opencl_memcpy2dev<int>( index_buffer, index, size, my_manager->get_command_queue ( ) );

    // copy to the device the values
    cl_mem value_buffer;
    opencl_memallocdev<ValueType>( value_buffer, size, my_manager->get_context ( ) );
    opencl_memcpy2dev<ValueType>( value_buffer, values, size, my_manager->get_command_queue ( ) );

    openclsetvalues ( this->buffer, index_buffer, value_buffer, size,
                      my_manager->get_kernel ( 0 ), my_manager->get_command_queue ( ),
                      this->global_threads, this->local_threads );

    opencl_memfreedev ( index_buffer );
    opencl_memfreedev ( value_buffer );
#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::GetValues ( const int *index, const int size, ValueType *values ) const
{
#ifdef WITH_OPENCL

    assert ( values != NULL );
    assert ( this->get_size ( ) > 0 );
    assert ( size != 0 );

    for ( int i = 0; i < size; ++i )
    {
        assert ( index[i] >= 0 );
        assert ( index[i] < this->get_size ( ) );
    }

    // copy the index set to the device
    cl_mem index_buffer;
    opencl_memallocdev<int>( index_buffer, size, my_manager->get_context ( ) );
    opencl_memcpy2dev<int>( index_buffer, index, size, my_manager->get_command_queue ( ) );
    // create buffer for the values on the device
    cl_mem value_buffer;
    opencl_memallocdev<ValueType>( value_buffer, size, my_manager->get_context ( ) );

    openclgetvalues ( this->buffer, index_buffer, value_buffer, size,
                      my_manager->get_kernel ( 1 ), my_manager->get_command_queue ( ),
                      this->global_threads, this->local_threads );

    //copy value buffer to host
    opencl_memcpy2host<ValueType>( values, value_buffer, size, my_manager->get_command_queue ( ) );

    opencl_memfreedev ( index_buffer );
    opencl_memfreedev ( value_buffer );
#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::GetBlockValues ( const int start_i, const int end_i, ValueType *values ) const
{
#ifdef WITH_OPENCL

    assert ( start_i >= 0 );
    assert ( start_i <= end_i );
    assert ( this->get_size ( ) > 0 );
    assert ( end_i <= this->get_size ( ) );
    assert ( values != NULL );

    int size = end_i - start_i;

    if ( size > 0 )
    {
        cl_mem value_buffer;

        opencl_memallocdev<ValueType>( value_buffer, size, my_manager->get_context ( ) );

        openclgetblockvalues ( this->buffer, value_buffer, start_i, end_i,
                               my_manager->get_kernel ( 13 ), my_manager->get_command_queue ( ),
                               this->global_threads, this->local_threads );

        // copy the values to the host
        opencl_memcpy2host ( values, value_buffer, size,
                             my_manager->get_command_queue ( ) );

        opencl_memfreedev ( value_buffer );
    }
    else
    {
        values = NULL;
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::SetBlockValues ( const int start_i, const int end_i, const ValueType *values )
{
#ifdef WITH_OPENCL

    assert ( start_i >= 0 );
    assert ( start_i <= end_i );
    assert ( this->get_size ( ) >= 0 );
    assert ( end_i <= this->get_size ( ) );
    assert ( values != NULL );

    int size = end_i - start_i;

    if ( ( size > 0 ) && ( this->get_size ( ) > 0 ) )
    {
        cl_mem value_buffer;

        opencl_memallocdev<ValueType>( value_buffer, size, my_manager->get_context ( ) );
        opencl_memcpy2dev<ValueType>( value_buffer, values, size, my_manager->get_command_queue ( ) );

        openclsetblockvalues ( this->buffer, value_buffer, start_i, 0, end_i - start_i,
                               my_manager->get_kernel ( 12 ), my_manager->get_command_queue ( ),
                               this->global_threads, this->local_threads );

        opencl_memfreedev ( value_buffer );
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
lVector<ValueType> *OPENCL_lVector<ValueType>::extract_subvector ( const int start_i, const int end_i ) const
{
#ifdef WITH_OPENCL

    assert ( start_i >= 0 );
    assert ( start_i < end_i );
    assert ( this->get_size ( ) > 0 );
    assert ( end_i <= this->get_size ( ) );

    lVector<ValueType> *sub_vector;
    std::string sub_vec_name;
    sub_vec_name = "sub vecrix from ";
    sub_vec_name.append ( this->get_name ( ) );

    sub_vector = new OPENCL_lVector<ValueType>( this->my_manager );
    sub_vector->Init ( end_i - start_i, sub_vec_name );

    if ( const OPENCL_lVector<ValueType> *casted_vec =
         dynamic_cast < const OPENCL_lVector<ValueType>* > ( sub_vector ) )
    {

        openclgetblockvalues ( this->buffer, casted_vec->buffer, start_i, end_i,
                               my_manager->get_kernel ( 13 ), my_manager->get_command_queue ( ),
                               this->global_threads, this->local_threads );

    }
    else
    {

        LOG_ERROR ( "ERROR in OPENCL_lVector<ValueType>::extract_subvector(); unsupported vector type" );
        this->print ( );
        sub_vector->print ( );
        exit ( -1 );
    }

    return sub_vector;

#else
    ERROR;
#endif
}

template <typename ValueType>
lVector<ValueType> &OPENCL_lVector<ValueType>::operator= ( const lVector<ValueType> &vec )
{
    if ( this == &vec )
        return *this;

    this->CopyFrom ( vec );
    return *this;
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::CopyFrom ( const lVector<ValueType>& vec )
{
#ifdef WITH_OPENCL
    if ( this != &vec )
    {
        // OPENCL = OPENCL  (with same opencl manager)
        if ( const OPENCL_lVector<ValueType> *casted_vec =
             dynamic_cast < const OPENCL_lVector<ValueType>* > ( &vec ) )
        {

            assert ( this->get_size ( ) == casted_vec->get_size ( ) );

            if ( this->get_size ( ) > 0 )
            {
                opencl_memcpydev<ValueType>( this->buffer, casted_vec->buffer, this->get_size ( ), my_manager->get_command_queue ( ) );
            }
        }
        else
        {
            // OPENCL = CPU
            if ( const CPU_lVector<ValueType> *casted_vec =
                 dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
            {

                assert ( this->get_size ( ) == casted_vec->get_size ( ) );

                if ( this->get_size ( ) > 0 )
                {
                    opencl_memcpy2dev<ValueType>( this->buffer, casted_vec->buffer, this->get_size ( ), my_manager->get_command_queue ( ) );
                }
            }
            else
            {
                // unsupported type
                LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::operator= unsupported vector type" );
                exit ( -1 );
                this->print ( );
                vec.print ( );
                exit ( -1 );
            }
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::CastFrom ( const lVector<double>& vec )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::CastFrom<double> not implemented yet" )
            this->print ( );
    ;
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::CastFrom ( const lVector<float>& vec )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::CastFrom<float> not implemented yet" );
    this->print ( );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::CastTo ( lVector<double>& vec ) const
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::CastTo<double> not implemented yet" )
            this->print ( );
    ;
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::CastTo ( lVector<float>& vec ) const
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::CastTo<float> not implemented yet" );
    this->print ( );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::ElementWiseMult ( const lVector<ValueType> &vec ) const
{
#ifdef WITH_OPENCL

    if ( const OPENCL_lVector<ValueType> *casted_vec =
         dynamic_cast < const OPENCL_lVector<ValueType>* > ( &vec ) )
    {

        this->ElementWiseMult ( *casted_vec );

    }
    else
    {
        std::cerr << "ERROR in CPU_lVector<ValueType>::ElementWiseMult() unsupported vector type\n";
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::ElementWiseMult ( const OPENCL_lVector<ValueType> &vec ) const
{
#ifdef WITH_OPENCL

    assert ( this->get_size ( ) == vec.get_size ( ) );

    if ( this->get_size ( ) > 0 )
        openclmultvalues ( vec.buffer, this->buffer, this->get_size ( ),
                           my_manager->get_kernel ( 14 ), my_manager->get_command_queue ( ),
                           this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Scale ( const ValueType scalar )
{
#ifdef WITH_OPENCL

    if ( this->get_size ( ) > 0 )
        openclscale ( this->buffer, this->get_size ( ), scalar,
                      my_manager->get_kernel ( 3 ), my_manager->get_command_queue ( ),
                      this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Axpy ( const OPENCL_lVector<ValueType> &vec, const ValueType scalar )
{
#ifdef WITH_OPENCL

    assert ( this->get_size ( ) == vec.get_size ( ) );

    if ( this->get_size ( ) > 0 )
        openclaxpy ( this->buffer, vec.buffer, this->get_size ( ), scalar,
                     my_manager->get_kernel ( 4 ), my_manager->get_command_queue ( ),
                     this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Axpy ( const lVector<ValueType> &vec, const ValueType scalar )
{
#ifdef WITH_OPENCL

    assert ( this->get_size ( ) == vec.get_size ( ) );

    if ( const OPENCL_lVector<ValueType> *casted_vec =
         dynamic_cast < const OPENCL_lVector<ValueType>* > ( &vec ) )
    {

        this->Axpy ( *casted_vec, scalar );

    }
    else
    {
        LOG_ERROR ( "ERROR OPENCL_lVector::Axpy unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::ScaleAdd ( const ValueType scalar, const OPENCL_lVector<ValueType> &vec )
{
#ifdef WITH_OPENCL

    assert ( this->get_size ( ) == vec.get_size ( ) );

    if ( this->get_size ( ) > 0 )
        openclscaledadd ( this->buffer, vec.buffer, this->get_size ( ), scalar,
                          my_manager->get_kernel ( 5 ), my_manager->get_command_queue ( ),
                          this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::ScaleAdd ( const ValueType scalar, const lVector<ValueType> &vec )
{
#ifdef WITH_OPENCL

    assert ( this->get_size ( ) == vec.get_size ( ) );

    if ( const OPENCL_lVector<ValueType> *casted_vec =
         dynamic_cast < const OPENCL_lVector<ValueType>* > ( &vec ) )
    {

        this->ScaleAdd ( scalar, *casted_vec );

    }
    else
    {
        LOG_ERROR ( "ERROR OPENCL_lVector::ScaleAdd unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
ValueType OPENCL_lVector<ValueType>::Dot ( const OPENCL_lVector<ValueType> &vec ) const
{
    ValueType result = 0;

#ifdef WITH_OPENCL

    assert ( this->get_size ( ) == vec.get_size ( ) );

    if ( this->get_size ( ) > 0 )
        result = opencldot<ValueType>( this->buffer, vec.buffer, this->get_size ( ),
            my_manager->get_context ( ), my_manager->get_kernel ( 11 ),
            my_manager->get_kernel ( 10 ), my_manager->get_command_queue ( ),
            this->global_threads, this->local_threads,
            this->global_blocks, my_manager->kerneltype,
            my_manager->cpuFinalReduction, my_manager->cpuFinalThreshold,
            my_manager->maxThreads, my_manager->max_work_group_size );

#else
    ERROR;
#endif
    return result;
}

template <typename ValueType>
ValueType OPENCL_lVector<ValueType>::Dot ( const lVector<ValueType> &vec ) const
{
    ValueType result = 0;

#ifdef WITH_OPENCL
    assert ( this->get_size ( ) == vec.get_size ( ) );
    if ( const OPENCL_lVector<ValueType> *casted_vec =
         dynamic_cast < const OPENCL_lVector<ValueType>* > ( &vec ) )
    {

        result = this->Dot ( *casted_vec );

    }
    else
    {
        LOG_ERROR ( "ERROR OPENCL_lVector::Dot unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

#else
    ERROR;
#endif
    return result;
}

template <typename ValueType>
ValueType OPENCL_lVector<ValueType>::Norm2 ( void ) const
{
    ValueType result = 0;
#ifdef WITH_OPENCL

    if ( this->get_size ( ) > 0 )
        result = openclnrm2<ValueType>( this->buffer, this->get_size ( ),
            my_manager->get_context ( ), my_manager->get_kernel ( 11 ),
            my_manager->get_kernel ( 10 ), my_manager->get_command_queue ( ),
            this->global_threads, this->local_threads, this->global_blocks,
            my_manager->kerneltype, my_manager->cpuFinalReduction,
            my_manager->cpuFinalThreshold, my_manager->maxThreads,
            my_manager->max_work_group_size );

#else
    ERROR;
#endif
    return result;
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::CopyTo ( lVector<ValueType>& vec ) const
{
#ifdef WITH_OPENCL

    if ( const OPENCL_lVector<ValueType> *casted_vec = dynamic_cast < const OPENCL_lVector<ValueType>* > ( &vec ) )
    {

        assert ( this->get_size ( ) == casted_vec->get_size ( ) );
        opencl_memcpydev<ValueType>( casted_vec->buffer, this->buffer, casted_vec->get_size ( ),
                my_manager->get_command_queue ( ) );

    }
    else
    {
        if ( const CPU_lVector<ValueType> *casted_vec = dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
        {

            assert ( this->get_size ( ) == casted_vec->get_size ( ) );
            opencl_memcpy2host<ValueType>( casted_vec->buffer, this->buffer, casted_vec->get_size ( ),
                    my_manager->get_command_queue ( ) );

        }
        else
        {
            LOG_ERROR ( "ERROR OPENCL_lVector::CopyTo unsupported vectors" );
            this->print ( );
            vec.print ( );
            exit ( -1 );
        }
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::partial_replace_subvector ( const int start_i, const int start_sub_vec, const int size, const lVector<ValueType> &sub_vec )
{
#ifdef WITH_OPENCL

    if ( const OPENCL_lVector<ValueType> *casted_vec =
         dynamic_cast < const OPENCL_lVector<ValueType>* > ( &sub_vec ) )
    {

        this->partial_replace_subvector ( start_i,
                                          start_sub_vec,
                                          size,
                                          *casted_vec );

    }
    else
    {
        std::cerr << "ERROR in OPENCL_lVector<ValueType>::partial_replace_subvector() unsupported vector type\n";
        this->print ( );
        sub_vec.print ( );
        exit ( -1 );
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::partial_replace_subvector ( const int start_i, const int start_sub_vec, const int size, const OPENCL_lVector<ValueType> &sub_vec )
{
#ifdef WITH_OPENCL

    assert ( start_i >= 0 );
    assert ( start_sub_vec >= 0 );
    assert ( size > 0 );
    assert ( start_sub_vec + size <= sub_vec.get_size ( ) );
    assert ( start_i + size <= this-> get_size ( ) );

    if ( size > 0 )
        openclsetblockvalues ( this->buffer, sub_vec.buffer, start_i, start_sub_vec, size,
                               my_manager->get_kernel ( 12 ), my_manager->get_command_queue ( ),
                               this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::partial_add_subvector ( const int start_i, const int start_sub_vec, const int size, const ValueType weight, const lVector<ValueType> &sub_vec )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::partial_add_subvector not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::partial_add_subvector ( const int start_i, const int start_sub_vec, const int size, const ValueType weight, const OPENCL_lVector<ValueType> &sub_vec )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::partial_add_subvector not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Sync ( void ) const
{
    // do nothing
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::set_indexset ( const int *indexset, const int size )
{
#ifdef WITH_OPENCL

    assert ( size >= 0 );

    for ( int i = 0; i < size; ++i )
    {
        assert ( indexset[i] >= 0 );
        assert ( indexset[i] < this->get_size ( ) );
    }

    if ( this->indexset_size_ > 0 )
    {
        opencl_memfreedev ( this->indexset_ );
        this->indexset_size_ = 0;
    }

    this->indexset_size_ = size;

    if ( size > 0 )
    {
        opencl_memallocdev<int>( this->indexset_, size,
                this->my_manager->get_context ( ) );
        opencl_memcpy2dev<int>( this->indexset_, indexset, size,
                this->my_manager->get_command_queue ( ) );
    }

#else
    ERROR;
#endif

}

template <typename ValueType>
void OPENCL_lVector<ValueType>::get_indexset ( int *indexset ) const
{
#ifdef WITH_OPENCL

    if ( this->get_indexset_size ( ) > 0 )
    {

        assert ( indexset != NULL );

        opencl_memcpy2host<int>( indexset, this->indexset_, this->get_indexset_size ( ), this->my_manager->get_command_queue ( ) );

    }
    else
    {
        indexset = NULL;
    }

#else
    ERROR;
#endif

}

template <typename ValueType>
void OPENCL_lVector<ValueType>::GetIndexedValues ( ValueType *values ) const
{
#ifdef WITH_OPENCL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_indexset_size ( ) >= 0 );

    if ( this->get_indexset_size ( ) > 0 )
    {
        cl_mem dev_values;

        // allocate the values on the device
        opencl_memallocdev<ValueType>( dev_values, this->indexset_size_, my_manager->get_context ( ) );

        // get the values from the vector on the values buffer (on the device)
        openclgetvalues ( this->buffer, this->indexset_, dev_values, this->indexset_size_,
                          this->my_manager->get_kernel ( 1 ), this->my_manager->get_command_queue ( ),
                          this->global_threads, this->local_threads );

        // copy the values (on the device) to the values (on the host/CPU)
        opencl_memcpy2host<ValueType>( values, dev_values, this->indexset_size_,
                my_manager->get_command_queue ( ) );

        opencl_memfreedev ( dev_values );
    }
    else
    {
        values = NULL;
    }

#else
    ERROR;
#endif

}

template <typename ValueType>
void OPENCL_lVector<ValueType>::SetIndexedValues ( const ValueType *values )
{
#ifdef WITH_OPENCL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_indexset_size ( ) >= 0 );

    if ( this->get_indexset_size ( ) > 0 )
    {

        cl_mem dev_values;

        // allocate the values on the device
        opencl_memallocdev<ValueType>( dev_values, this->indexset_size_, my_manager->get_context ( ) );
        // copy the values on the device
        opencl_memcpy2dev<ValueType>( dev_values, values, this->indexset_size_, this->my_manager->get_command_queue ( ) );

        // set the values from the vector on the values buffer
        openclsetvalues ( this->buffer, this->indexset_, dev_values, this->indexset_size_,
                          this->my_manager->get_kernel ( 0 ), this->my_manager->get_command_queue ( ),
                          this->global_threads, this->local_threads );

        opencl_memfreedev ( dev_values );
    }

#else
    ERROR;
#endif

}

template <typename ValueType>
int OPENCL_lVector<ValueType>::get_indexset_size ( void ) const
{
    return this->indexset_size_;
}

template <typename ValueType>
ValueType OPENCL_lVector<ValueType>::Norm1 ( void ) const
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::Sum not implemented yet" );
    exit ( -1 );
    return -1.0;
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Rot ( lVector<ValueType> *vec, const ValueType &sc, const ValueType &ss )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::Rot not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Rot ( OPENCL_lVector<ValueType> *vec, const ValueType &sc, const ValueType &ss )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::Rot not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Rotg ( ValueType *sa, ValueType *sb, ValueType *sc, ValueType *ss ) const
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::Rotg not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Rotm ( lVector<ValueType> *vec, const ValueType &sparam )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::Rotm not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Rotm ( OPENCL_lVector<ValueType> *vec, const ValueType &sparam )
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::Rotm not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_lVector<ValueType>::Rotmg ( ValueType *sd1, ValueType *sd2, ValueType *x1, const ValueType &x2, ValueType *sparam ) const
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::Rotmg not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
int OPENCL_lVector<ValueType>::ArgMin ( void ) const
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::ArgMin() not implemented yet" );
    exit ( -1 );
    return -1;
}

template <typename ValueType>
int OPENCL_lVector<ValueType>::ArgMax ( void ) const
{
    LOG_ERROR ( "ERROR OPENCL_lVector<ValueType>::ArgMax() not implemented yet" );
    exit ( -1 );
    return -1;
}

template <typename ValueType>
template <typename OtherValueType>
void OPENCL_lVector<ValueType>::CastFrom ( const hiflow::la::lVector<OtherValueType>& other )
{
    if ( const OPENCL_lVector<OtherValueType>* other_opencl_vec
         = dynamic_cast < const OPENCL_lVector<OtherValueType>* > ( &other ) )
    {
        this->CastFrom ( *other_opencl_vec );
    }
    else if ( const CPU_lVector<OtherValueType>* other_cpu_vec
              = dynamic_cast < const CPU_lVector<OtherValueType>* > ( &other ) )
    {
        this->CastFrom ( *other_cpu_vec );
    }
    else
    {
        LOG_ERROR ( "OPENCL_lVector<ValueType>::CastFrom<OtherValueType> unsupported other vector" );
        this->print ( );
        other.print ( );
    }
}

template <typename ValueType>
template <typename OtherValueType>
void OPENCL_lVector<ValueType>::CastFrom ( const OPENCL_lVector<OtherValueType>& other )
{
    ValueType *vals = new ValueType[this->get_size ( )];
    OtherValueType *other_vals = new OtherValueType[other.get_size ( )];

    opencl_memcpy2host<OtherValueType>( other_vals, other.buffer, other.get_size ( ), other.my_manager->get_command_queue ( ) );

    for ( int i = 0; i < this->get_size ( ); ++i )
        vals[i] = static_cast < ValueType > ( other_vals[i] );

    opencl_memcpy2dev<ValueType>( this->buffer, vals, this->get_size ( ) );

    delete [] vals;
    delete [] other_vals;
}

template <typename ValueType>
template <typename OtherValueType>
void OPENCL_lVector<ValueType>::CastFrom ( const CPU_lVector<OtherValueType>& other )
{
    ValueType *vals = new ValueType[this->get_size ( )];

    for ( int i = 0; i < this->get_size ( ); ++i )
        vals[i] = static_cast < ValueType > ( other.buffer[i] );

    opencl_memcpy2dev<ValueType>( this->buffer, vals, this->get_size ( ) );

    delete [] vals;
}

template class OPENCL_lVector<float>;
template class OPENCL_lVector<double>;
