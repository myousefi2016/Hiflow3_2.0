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

#include "cuda/lvector_gpu.h"
#include "../lmp_mem.h"
#include "../init_vec_mat.h"
#include "../lmp_log.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

using namespace hiflow::la;

// class GPU_lVector

template <typename ValueType>
GPU_lVector<ValueType>::GPU_lVector ( )
{
    this->platform_name_ = "GPU (CUDA)";
    this->platform_id_ = GPU;
    this->indexset_size_ = 0;
    this->buffer = NULL;
    this->indexset_ = NULL;
}

template <typename ValueType>
GPU_lVector<ValueType>::~GPU_lVector ( )
{
    this->Clear ( );
}

template <typename ValueType>
void GPU_lVector<ValueType>::Clear ( )
{

    if ( this->indexset_ != NULL )
        memfreedev ( this->indexset_ );
    this->indexset_ = NULL;

    if ( this->buffer != NULL )
        memfreedev ( this->buffer );
    this->buffer = NULL;

    this->indexset_size_ = 0;
    this->size_ = 0;
    this->name_ = "";
}

template <typename ValueType>
void GPU_lVector<ValueType>::Zeros ( )
{
    memsetdev ( this->buffer, 0, this->get_size ( ), sizeof (ValueType ) );
}

template <typename ValueType>
void GPU_lVector<ValueType>::Reorder ( const int *index )
{

    LOG_ERROR ( "ERROR GPU_lVector<ValueType>::Reorder() not implemented yet" );
    this->print ( );
    exit ( -1 );

}

template <typename ValueType>
void GPU_lVector<ValueType>::set_thread_blocks ( int thread_block, int thread_block_size )
{
    this->thread_block_ = thread_block;
    this->thread_block_size_ = thread_block_size;
}

template <typename ValueType>
void GPU_lVector<ValueType>::set_thread_blocks ( void )
{
    // default values
    int thread_block_size = 256;
    int thread_block = this->get_size ( ) / thread_block_size;
    if ( thread_block_size * thread_block < this->get_size ( ) )
        thread_block++;

    this->set_thread_blocks ( thread_block, thread_block_size );
}

template <typename ValueType>
void GPU_lVector<ValueType>::Init ( const int size, const std::string name )
{

    assert ( size >= 0 );

    this->Clear ( );

    memallocdev ( ( void** ) & this->buffer, size, sizeof (ValueType ) );
    memsetdev ( this->buffer, 0, size, sizeof (ValueType ) );

    this->name_ = name;
    this->size_ = size;
    this->set_thread_blocks ( );

    this->Zeros ( );

}

template <typename ValueType>
void GPU_lVector<ValueType>::add_value ( const int i, ValueType val )
{
    assert ( this->get_size ( ) > 0 );

    LOG_ERROR ( "ERROR GPU_lVector::add_value is not supported for this platform" );
    this->print ( );
    exit ( -1 );
}

template <typename ValueType>
void GPU_lVector<ValueType>::add_values ( const int* indices, int length, const ValueType* values )
{
    assert ( this->get_size ( ) > 0 );

    LOG_ERROR ( "ERROR GPU_lVector::add_values is not supported for this platform" );
    this->print ( );
    exit ( -1 );
}

template <typename ValueType>
void GPU_lVector<ValueType>::SetValues ( const int *index, const int size, const ValueType *values )
{
    int thread_block_size = this->thread_block_size_;

    int *dev_index;
    ValueType *dev_values;

    assert ( values != NULL );
    assert ( this->get_size ( ) > 0 );
    assert ( size != 0 );

    for ( int i = 0; i < size; ++i )
    {
        assert ( index[i] >= 0 );
        assert ( index[i] < this->get_size ( ) );
    }

    // copy the index set to the device
    memallocdev ( ( void** ) &dev_index, size, sizeof (int ) );
    memcpy2dev ( dev_index, index, size );

    // copy to the device the values
    memallocdev ( ( void** ) &dev_values, size, sizeof (ValueType ) );
    memcpy2dev ( dev_values, values, size );

    cudasetvalues ( dev_index, size, dev_values, this->buffer, thread_block_size );

    memfreedev ( dev_index );
    memfreedev ( dev_values );

}

template <typename ValueType>
void GPU_lVector<ValueType>::GetValues ( const int *index, const int size, ValueType *values ) const
{
    int thread_block_size = this->thread_block_size_;

    int *dev_index;
    ValueType *dev_values;

    assert ( values != NULL );
    assert ( this->get_size ( ) > 0 );
    assert ( size > 0 );

    for ( int i = 0; i < size; ++i )
    {
        assert ( index[i] >= 0 );
        assert ( index[i] < this->get_size ( ) );
    }

    // copy to the device the index set
    memallocdev ( ( void** ) &dev_index, size, sizeof (int ) );
    memcpy2dev ( dev_index, index, size );

    memallocdev ( ( void** ) &dev_values, size, sizeof (ValueType ) );

    cudagetvalues ( dev_index, size, dev_values, this->buffer, thread_block_size );

    // copy the values to the host
    memcpy2host ( values, dev_values, size );

    memfreedev ( dev_index );
    memfreedev ( dev_values );

}

template <typename ValueType>
int GPU_lVector<ValueType>::get_indexset_size ( void ) const
{
    return this->indexset_size_;
}

template <typename ValueType>
void GPU_lVector<ValueType>::set_indexset ( const int *indexset,
                                            const int size )
{
    assert ( size >= 0 );

    for ( int i = 0; i < size; ++i )
    {
        assert ( indexset[i] >= 0 );
        assert ( indexset[i] < this->get_size ( ) );
    }

    if ( this->indexset_size_ > 0 )
    {
        memfreedev ( this->indexset_ );
        this->indexset_size_ = 0;
    }

    this->indexset_size_ = size;

    if ( size > 0 )
    {
        memallocdev ( ( void** ) & this->indexset_, size, sizeof (ValueType ) );
        memcpy2dev ( this->indexset_, indexset, size );
    }

}

template <typename ValueType>
void GPU_lVector<ValueType>::get_indexset ( int *indexset ) const
{

    if ( this->get_indexset_size ( ) > 0 )
    {

        assert ( indexset != NULL );

        memcpy2host ( indexset, this->indexset_, this->get_indexset_size ( ) );

    }
    else
    {
        indexset = NULL;
    }

}

template <typename ValueType>
void GPU_lVector<ValueType>::GetIndexedValues ( ValueType *values ) const
{

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_indexset_size ( ) >= 0 );

    if ( this->get_indexset_size ( ) > 0 )
    {
        int thread_block_size = this->thread_block_size_;

        ValueType *dev_values;

        // allocate the values on the device
        memallocdev ( ( void** ) &dev_values, this->indexset_size_, sizeof (ValueType ) );

        // get the values from the vector on the values buffer (on the device)
        cudagetvalues ( this->indexset_, this->indexset_size_, dev_values, this->buffer, thread_block_size );

        // copy the values (on the device) to the values (on the host/CPU)
        memcpy2host ( values, dev_values, this->indexset_size_ );

        memfreedev ( dev_values );
    }
    else
    {
        values = NULL;
    }

}

template <typename ValueType>
void GPU_lVector<ValueType>::SetIndexedValues ( const ValueType *values )
{

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_indexset_size ( ) >= 0 );

    if ( this->get_indexset_size ( ) > 0 )
    {
        int thread_block_size = this->thread_block_size_;

        ValueType *dev_values;

        // allocate the values on the device
        memallocdev ( ( void** ) &dev_values, this->indexset_size_, sizeof (ValueType ) );
        // copy the values on the device
        memcpy2dev ( dev_values, values, this->indexset_size_ );

        // set the values from the vector on the values buffer
        cudasetvalues ( this->indexset_, this->indexset_size_, dev_values, this->buffer, thread_block_size );

        memfreedev ( dev_values );
    }

}

template <typename ValueType>
void GPU_lVector<ValueType>::SetBlockValues ( const int start_i,
                                              const int end_i,
                                              const ValueType *values )
{
    assert ( start_i >= 0 );
    assert ( start_i <= end_i );
    assert ( this->get_size ( ) >= 0 );
    assert ( end_i <= this->get_size ( ) );
    assert ( values != NULL );

    int size = end_i - start_i;
    int thread_block_size = this->thread_block_size_;

    if ( ( size > 0 ) && ( this->get_size ( ) > 0 ) )
    {
        ValueType *dev_values;

        memallocdev ( ( void** ) &dev_values, size, sizeof (ValueType ) );
        memcpy2dev ( dev_values, values, size );

        cudasetblockvalues ( start_i, 0, end_i - start_i, dev_values, this->buffer, thread_block_size );

        memfreedev ( dev_values );
    }
}

template <typename ValueType>
void GPU_lVector<ValueType>::GetBlockValues ( const int start_i,
                                              const int end_i,
                                              ValueType *values ) const
{
    assert ( start_i >= 0 );
    assert ( start_i <= end_i );
    assert ( this->get_size ( ) > 0 );
    assert ( end_i <= this->get_size ( ) );
    assert ( values != NULL );

    int size = end_i - start_i;
    int thread_block_size = this->thread_block_size_;

    if ( size > 0 )
    {
        ValueType *dev_values;

        memallocdev ( ( void** ) &dev_values, size, sizeof (ValueType ) );

        cudagetblockvalues ( start_i, end_i, dev_values, this->buffer, thread_block_size );

        // copy the values to the host
        memcpy2host ( values, dev_values, size );

        memfreedev ( dev_values );
    }
    else
    {
        values = NULL;
    }

}

template <typename ValueType>
lVector<ValueType> &GPU_lVector<ValueType>::operator= ( const lVector<ValueType> &vec2 )
{

    if ( this == &vec2 )
        return *this;

    this->CopyFrom ( vec2 );
    return *this;

}

template <typename ValueType>
void GPU_lVector<ValueType>::CopyFrom ( const lVector<ValueType>& vec )
{

    if ( this != &vec )
    {

        // GPU = GPU
        if ( const GPU_lVector<ValueType> *casted_vec =
             dynamic_cast < const GPU_lVector<ValueType>* > ( &vec ) )
        {

            assert ( this->get_size ( ) == casted_vec->get_size ( ) );

            if ( this->get_size ( ) > 0 )
            {
                memcpydev ( this->buffer, casted_vec->buffer, this->get_size ( ) );
            }

        }
        else
        {

            // GPU = CPU
            if ( const CPU_lVector<ValueType> *casted_vec =
                 dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
            {

                assert ( this->get_size ( ) == casted_vec->get_size ( ) );

                if ( this->get_size ( ) > 0 )
                {
                    memcpy2dev ( this->buffer, casted_vec->buffer, this->get_size ( ) );
                }

            }
            else
            {

                // unsupported type
                LOG_ERROR ( "ERROR GPU_lVector<ValueType>::CopyFrom() unsupported vector type" );
                this->print ( );
                vec.print ( );
                exit ( -1 );
            }
        }
    }
}

template <typename ValueType>
void GPU_lVector<ValueType>::CopyTo ( lVector<ValueType>& vec ) const
{

    if ( this != &vec )
    {

        // GPU = GPU
        if ( const GPU_lVector<ValueType> *casted_vec =
             dynamic_cast < const GPU_lVector<ValueType>* > ( &vec ) )
        {

            assert ( this->get_size ( ) == casted_vec->get_size ( ) );

            if ( this->get_size ( ) > 0 )
            {
                memcpydev ( casted_vec->buffer, this->buffer, this->get_size ( ) );
            }

        }
        else
        {

            // CPU = GPU
            if ( const CPU_lVector<ValueType> *casted_vec =
                 dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
            {

                assert ( this->get_size ( ) == casted_vec->get_size ( ) );

                if ( this->get_size ( ) > 0 )
                {
                    memcpy2host ( casted_vec->buffer, this->buffer, this->get_size ( ) );
                }

            }
            else
            {

                // unsupported type
                LOG_ERROR ( "ERROR GPU_lVector<ValueType>::CopyTo() unsupported vector type" );
                this->print ( );
                vec.print ( );
                exit ( -1 );
            }
        }
    }
}

template <typename ValueType>
lVector<ValueType> *GPU_lVector<ValueType>::extract_subvector ( const int start_i,
                                                                const int end_i ) const
{
    assert ( start_i >= 0 );
    assert ( start_i < end_i );
    assert ( this->get_size ( ) > 0 );
    assert ( end_i <= this->get_size ( ) );

    lVector<ValueType> *sub_vector;
    std::string sub_vec_name;
    sub_vec_name = "sub vecrix from ";
    sub_vec_name.append ( this->get_name ( ) );

    sub_vector = init_vector<ValueType>( end_i - start_i,
            sub_vec_name,
            this->get_platform ( ),
            this->get_implementation ( ) );

    if ( const GPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const GPU_lVector<ValueType>* > ( sub_vector ) )
    {

        int thread_block_size = this->thread_block_size_;

        cudagetblockvalues ( start_i, end_i, casted_vec->buffer, this->buffer, thread_block_size );

        // actually this cannot happen
    }
    else
    {
        LOG_ERROR ( "ERROR in GPU_lVector<ValueType>::extract_subvector(); unsupported vector type" );
        this->print ( );
        sub_vector->print ( );
    }

    return sub_vector;
}

template <typename ValueType>
void GPU_lVector<ValueType>::partial_replace_subvector ( const int start_i,
                                                         const int start_sub_vec,
                                                         const int size,
                                                         const lVector<ValueType> &sub_vec )
{

    if ( const GPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const GPU_lVector<ValueType>* > ( &sub_vec ) )
    {

        this->partial_replace_subvector ( start_i,
                                          start_sub_vec,
                                          size,
                                          *casted_vec );

    }
    else
    {
        LOG_ERROR ( "ERROR in GPU_lVector<ValueType>::partial_replace_subvector() unsupported vector type" );
        this->print ( );
        sub_vec.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void GPU_lVector<ValueType>::partial_replace_subvector ( const int start_i,
                                                         const int start_sub_vec,
                                                         const int size,
                                                         const GPU_lVector<ValueType> &sub_vec )
{
    assert ( start_i >= 0 );
    assert ( start_sub_vec >= 0 );
    assert ( size > 0 );
    assert ( start_sub_vec + size <= sub_vec.get_size ( ) );
    assert ( start_i + size <= this-> get_size ( ) );

    int thread_block_size = this->thread_block_size_;

    if ( size > 0 )
        cudasetblockvalues ( start_i, start_sub_vec, size, sub_vec.buffer, this->buffer, thread_block_size );

}

template <typename ValueType>
void GPU_lVector<ValueType>::partial_add_subvector ( const int start_i,
                                                     const int start_sub_vec,
                                                     const int size,
                                                     const ValueType weight,
                                                     const lVector<ValueType> &sub_vec )
{
    if ( const GPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const GPU_lVector<ValueType>* > ( &sub_vec ) )
    {

        this->partial_add_subvector ( start_i,
                                      start_sub_vec,
                                      size,
                                      weight,
                                      *casted_vec );

    }
    else
    {
        LOG_ERROR ( "ERROR in GPU_lVector<ValueType>::partial_add_subvector() unsupported vector type" );
        this->print ( );
        sub_vec.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void GPU_lVector<ValueType>::partial_add_subvector ( const int start_i,
                                                     const int start_sub_vec,
                                                     const int size,
                                                     const ValueType weight,
                                                     const GPU_lVector<ValueType> &sub_vec )
{
    assert ( start_i >= 0 );
    assert ( start_sub_vec >= 0 );
    assert ( size > 0 );
    assert ( start_sub_vec + size <= sub_vec.get_size ( ) );
    assert ( start_i + size <= this-> get_size ( ) );

    int thread_block_size = this->thread_block_size_;

    if ( this->get_size ( ) > 0 )
        cudaaddblockvalues ( start_i, start_sub_vec, size, sub_vec.buffer, this->buffer, weight, thread_block_size );

}

template <typename ValueType>
void GPU_lVector<ValueType>::Sync ( void ) const
{
    cuda_sync_threads ( );
}

template <typename ValueType>
void GPU_lVector<ValueType>::ElementWiseMult ( const lVector<ValueType> &vec ) const
{

    if ( const GPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const GPU_lVector<ValueType>* > ( &vec ) )
    {

        this->ElementWiseMult ( *casted_vec );

    }
    else
    {
        LOG_ERROR ( "ERROR in CPU_lVector<ValueType>::ElementWiseMult() unsupported vector type" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void GPU_lVector<ValueType>::ElementWiseMult ( const GPU_lVector<ValueType> &vec ) const
{

    assert ( this->get_size ( ) == vec.get_size ( ) );

    int thread_block_size = this->thread_block_size_;

    if ( this->get_size ( ) > 0 )
        cudamultvalues ( this->get_size ( ), vec.buffer, this->buffer, thread_block_size );

}

template <typename ValueType>
void GPU_lVector<ValueType>::CastFrom ( const hiflow::la::lVector<double>& other )
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    if ( const GPU_lVector<double>* other_gpu_vec
         = dynamic_cast < const GPU_lVector<double>* > ( &other ) )
    {
        this->CastFrom ( *other_gpu_vec );
    }
    else if ( const CPU_lVector<double>* other_cpu_vec
              = dynamic_cast < const CPU_lVector<double>* > ( &other ) )
    {
        this->CastFrom ( *other_cpu_vec );
    }
    else
    {
        LOG_ERROR ( "GPU_lVector<ValueType>::CastFrom<double> called with unsupported vector argument." );
        this->print ( );
        other.print ( );
    }
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastFrom ( const GPU_lVector<double>& other )
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    cudacastfromdouble ( this->get_size ( ), this->buffer, other.buffer, this->thread_block_size_ );
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastFrom ( const CPU_lVector<double>& other )
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    ValueType *vals = new ValueType[this->get_size ( )];

    for ( int i = 0; i < this->get_size ( ); ++i )
    {
        vals[i] = static_cast < ValueType > ( other.buffer[i] );
    }

    memcpy2dev ( this->buffer, vals, this->get_size ( ) );

    delete [] vals;
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastFrom ( const hiflow::la::lVector<float>& other )
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    if ( const GPU_lVector<float>* other_gpu_vec
         = dynamic_cast < const GPU_lVector<float>* > ( &other ) )
    {
        this->CastFrom ( *other_gpu_vec );
    }
    else if ( const CPU_lVector<float>* other_cpu_vec
              = dynamic_cast < const CPU_lVector<float>* > ( &other ) )
    {
        this->CastFrom ( *other_cpu_vec );
    }
    else
    {
        LOG_ERROR ( "GPU_lVector<ValueType>::CastFrom<float> called with unsupported vector argument." );
        this->print ( );
        other.print ( );
    }
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastFrom ( const GPU_lVector<float>& other )
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    cudacastfromfloat ( this->get_size ( ), this->buffer, other.buffer, this->thread_block_size_ );
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastFrom ( const CPU_lVector<float>& other )
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    ValueType *vals = new ValueType[this->get_size ( )];

    for ( int i = 0; i < this->get_size ( ); ++i )
    {
        vals[i] = static_cast < ValueType > ( other.buffer[i] );
    }

    memcpy2dev ( this->buffer, vals, this->get_size ( ) );

    delete [] vals;
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastTo ( hiflow::la::lVector<double>& other ) const
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    if ( GPU_lVector<double>* other_gpu_vec
         = dynamic_cast < GPU_lVector<double>* > ( &other ) )
    {
        this->CastTo ( *other_gpu_vec );
    }
    else if ( CPU_lVector<double>* other_cpu_vec
              = dynamic_cast < CPU_lVector<double>* > ( &other ) )
    {
        this->CastTo ( *other_cpu_vec );
    }
    else
    {
        LOG_ERROR ( "GPU_lVector<ValueType>::CastTo<double> called with unsupported vector argument." );
        this->print ( );
        other.print ( );
    }
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastTo ( GPU_lVector<double>& other ) const
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    cudacasttodouble ( this->get_size ( ), other.buffer, this->buffer, this->thread_block_size_ );
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastTo ( CPU_lVector<double>& other ) const
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    ValueType *vals = new ValueType[this->get_size ( )];

    for ( int i = 0; i < this->get_size ( ); ++i )
    {
        vals[i] = static_cast < ValueType > ( other.buffer[i] );
    }

    memcpy2dev ( this->buffer, vals, this->get_size ( ) );

    delete [] vals;
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastTo ( hiflow::la::lVector<float>& other ) const
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    if ( GPU_lVector<float>* other_gpu_vec
         = dynamic_cast < GPU_lVector<float>* > ( &other ) )
    {
        this->CastTo ( *other_gpu_vec );
    }
    else if ( CPU_lVector<float>* other_cpu_vec
              = dynamic_cast < CPU_lVector<float>* > ( &other ) )
    {
        this->CastTo ( *other_cpu_vec );
    }
    else
    {
        LOG_ERROR ( "GPU_lVector<ValueType>::CastTo<float> called with unsupported vector argument." );
        this->print ( );
        other.print ( );
    }
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastTo ( GPU_lVector<float>& other ) const
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    cudacasttofloat ( this->get_size ( ), other.buffer, this->buffer, this->thread_block_size_ );
}

template <typename ValueType>
void GPU_lVector<ValueType>::CastTo ( CPU_lVector<float>& other ) const
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    ValueType *vals = new ValueType[this->get_size ( )];

    memcpy2host ( vals, this->buffer, this->get_size ( ) );

    for ( int i = 0; i < this->get_size ( ); ++i )
    {
        other.buffer[i] = static_cast < float > ( vals[i] );
    }

    delete [] vals;
}

template class GPU_lVector<double>;
template class GPU_lVector<float>;
