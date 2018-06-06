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

/// @author Dimitar Lukarski, Simon Gawlok, Martin Wlotzka

#include "lvector_cpu.h"
#include "lmp_mem.h"
#include "init_vec_mat.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <typeinfo>

#include "lmp_log.h"

using namespace hiflow::la;

// class CPU_lVector

template <typename ValueType>
CPU_lVector<ValueType>::CPU_lVector ( )
{
    this->platform_name_ = "CPU (x86)";
    this->platform_id_ = CPU;
    this->indexset_size_ = 0;
    this->buffer = NULL;
    this->indexset_ = NULL;
}

template <typename ValueType>
CPU_lVector<ValueType>::~CPU_lVector ( )
{
    this->Clear ( );
}

template <typename ValueType>
void CPU_lVector<ValueType>::Clear ( )
{

    if ( this->indexset_ != NULL )
        delete [] this->indexset_;
    this->indexset_ = NULL;

    if ( this->buffer != NULL )
        delete [] this->buffer;
    this->buffer = NULL;

    this->indexset_size_ = 0;
    this->size_ = 0;
    this->name_ = "";
}

template <typename ValueType>
void CPU_lVector<ValueType>::Zeros ( )
{
    memsethost ( this->buffer, 0, this->get_size ( ), sizeof (ValueType ) );
}

template <typename ValueType>
void CPU_lVector<ValueType>::Reorder ( const int *index )
{

    ValueType *old_buffer = new ValueType[this->get_size ( )];

    assert ( old_buffer != NULL );

    // copy the buffer
    for ( int i = 0, i_e = this->get_size ( ); i != i_e; ++i )
        old_buffer[i] = this->buffer[i];

    // permute
    for ( int i = 0, i_e = this->get_size ( ); i != i_e; ++i )
        this->buffer[ i ] = old_buffer[index[i]];
    //    this->buffer[ index[i] ] = old_buffer[i];

    delete [] old_buffer;

}

template <typename ValueType>
void CPU_lVector<ValueType>::Init ( const int size, const std::string name )
{

    assert ( size >= 0 );

    this->Clear ( );

    this->buffer = new ValueType[size];
    assert ( this->buffer != NULL );

    this->name_ = name;
    this->size_ = size;

    this->Zeros ( );

}

template <typename ValueType>
void CPU_lVector<ValueType>::add_value ( const int i, ValueType val )
{
    assert ( 0 <= i );
    assert ( i < this->get_size ( ) );

    this->buffer[i] += val;
}

template <typename ValueType>
void CPU_lVector<ValueType>::add_values ( const int* indices, int length, const ValueType* values )
{
    assert ( 0 <= length );
    assert ( !( ( length > 0 ) && ( indices == 0 ) ) );
    assert ( !( ( length > 0 ) && ( values == 0 ) ) );

    for ( int i ( 0 ); i != length; ++i )
    {
        assert ( indices[i] < this->get_size ( ) );
        this->buffer[indices[i]] += values[i];
    }
}

template <typename ValueType>
void CPU_lVector<ValueType>::SetValues ( const int *index, const int size, const ValueType *values )
{
    assert ( 0 <= size );
    assert ( !( ( size > 0 ) && ( index == 0 ) ) );
    assert ( !( ( size > 0 ) && ( values == 0 ) ) );

    for ( int i ( 0 ); i != size; ++i )
    {
        const int current = index[i];
        assert ( current >= 0 );
        assert ( current < this->get_size ( ) );
        this->buffer[ current ] = values[i];
    }
}

template <typename ValueType>
void CPU_lVector<ValueType>::GetValues ( const int *index, const int size, ValueType *values ) const
{
    assert ( 0 <= size );
    assert ( !( ( size > 0 ) && ( index == 0 ) ) );
    assert ( !( ( size > 0 ) && ( values == 0 ) ) );

    for ( int i ( 0 ); i != size; ++i )
    {
        const int current = index[i];
        assert ( current >= 0 );
        assert ( current < this->get_size ( ) );
        values[i] = this->buffer[ current ];
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::SetBlockValues ( const int start_i,
                                              const int end_i,
                                              const ValueType *values )
{
    assert ( start_i >= 0 );
    assert ( start_i <= end_i );
    assert ( this->get_size ( ) >= 0 );
    assert ( end_i <= this->get_size ( ) );

    if ( ( this->get_size ( ) > 0 ) &&
         ( start_i < end_i ) &&
         ( values != NULL ) )
    {
        for ( int i = 0, i_e = end_i - start_i; i != i_e; ++i )
        {
            this->buffer[start_i + i] = values[i];
        }
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::GetBlockValues ( const int start_i,
                                              const int end_i,
                                              ValueType *values ) const
{
    assert ( start_i >= 0 );
    assert ( start_i <= end_i );
    assert ( this->get_size ( ) >= 0 );
    assert ( end_i <= this->get_size ( ) );

    if ( ( this->get_size ( ) > 0 ) &&
         ( start_i < end_i ) )
    {
        for ( int i = 0, i_e = end_i - start_i; i != i_e; ++i )
        {
            values[i] = this->buffer[start_i + i];
        }
    }
    else
    {
        values = NULL;
    }

}

template <typename ValueType>
int CPU_lVector<ValueType>::get_indexset_size ( void ) const
{
    return this->indexset_size_;
}

template <typename ValueType>
void CPU_lVector<ValueType>::set_indexset ( const int *indexset,
                                            const int size )
{
    assert ( size >= 0 );

    if ( this->indexset_size_ > 0 )
    {
        delete [] this->indexset_;
        this->indexset_size_ = 0;
    }

    this->indexset_size_ = size;

    if ( size > 0 )
    {
        this->indexset_ = new int[size];
        assert ( this->indexset_ != NULL );

        memcpyhost ( this->indexset_, indexset, size );
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::get_indexset ( int *indexset ) const
{

    if ( this->get_indexset_size ( ) > 0 )
    {

        assert ( indexset != NULL );

        memcpyhost ( indexset, this->indexset_, this->get_indexset_size ( ) );

    }
    else
    {
        indexset = NULL;
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::GetIndexedValues ( ValueType *values ) const
{
    assert ( this->get_size ( ) >= 0 );
    assert ( this->get_indexset_size ( ) >= 0 );

    if ( ( this->get_indexset_size ( ) > 0 ) &&
         ( this->get_size ( ) > 0 ) )
    {
        this->GetValues ( this->indexset_, this->indexset_size_, values );
    }
    else
    {
        values = NULL;
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::SetIndexedValues ( const ValueType *values )
{
    assert ( this->get_size ( ) >= 0 );
    assert ( this->get_indexset_size ( ) >= 0 );

    if ( ( this->get_indexset_size ( ) > 0 ) &&
         ( this->get_size ( ) > 0 ) &&
         ( values != NULL ) )
    {
        this->SetValues ( this->indexset_, this->indexset_size_, values );
    }

}

template <typename ValueType>
lVector<ValueType> &CPU_lVector<ValueType>::operator= ( const lVector<ValueType> &vec2 )
{
    if ( this == &vec2 )
        return *this;

    this->CopyFrom ( vec2 );
    return *this;

}

template <typename ValueType>
void CPU_lVector<ValueType>::CopyFrom ( const lVector<ValueType>& vec )
{

    if ( this != &vec )
    {

        // CPU = CPU
        if ( const CPU_lVector<ValueType> *casted_vec =
             dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
        {

            assert ( this->get_size ( ) == casted_vec->get_size ( ) );

            if ( this->get_size ( ) > 0 )
            {
                memcpyhost ( this->buffer, casted_vec->buffer, this->get_size ( ) );
            }

        }
        else
        {

            // if vec is not a CPU vector
            vec.CopyTo ( *this );

        }
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::CopyTo ( lVector<ValueType>& vec ) const
{
    vec.CopyFrom ( *this );
}

template <typename ValueType>
void CPU_lVector<ValueType>::CastFrom ( const lVector<double>& other )
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    if ( const CPU_lVector<double> *other_cpu =
         dynamic_cast < const CPU_lVector<double>* > ( &other ) )
    {
        // CPU from CPU
        this->CastFrom ( *other_cpu );
    }
    else
    {
        // CPU from non-CPU via non-CPU to CPU
        other.CastTo ( *this );
    }
}

template <typename ValueType>
void CPU_lVector<ValueType>::CastFrom ( const lVector<float>& other )
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    if ( const CPU_lVector<float> *other_cpu =
         dynamic_cast < const CPU_lVector<float>* > ( &other ) )
    {
        // CPU from CPU
        this->CastFrom ( *other_cpu );
    }
    else
    {
        // CPU from non-CPU via non-CPU to CPU
        other.CastTo ( *this );
    }
}

template <typename ValueType>
void CPU_lVector<ValueType>::CastTo ( lVector<double>& other ) const
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    if ( CPU_lVector<double> *other_cpu =
         dynamic_cast < CPU_lVector<double>* > ( &other ) )
    {
        // CPU to CPU
        this->CastTo ( *other_cpu );
    }
    else
    {
        // CPU to non-CPU via non-CPU from CPU
        other.CastFrom ( *this );
    }
}

template <typename ValueType>
void CPU_lVector<ValueType>::CastTo ( lVector<float>& other ) const
{
    assert ( this->get_size ( ) == other.get_size ( ) );

    if ( CPU_lVector<float> *other_cpu =
         dynamic_cast < CPU_lVector<float>* > ( &other ) )
    {
        // CPU to CPU
        this->CastTo ( *other_cpu );
    }
    else
    {
        // CPU to non-CPU via non-CPU from CPU
        other.CastFrom ( *this );
    }
}

template <typename ValueType>
lVector<ValueType> *CPU_lVector<ValueType>::extract_subvector ( const int start_i,
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

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( sub_vector ) )
    {

        this->GetBlockValues ( start_i,
                               end_i,
                               casted_vec->buffer );

        // actually this cannot happen
    }
    else
    {
        LOG_ERROR ( "CPU_lVector<ValueType>::extract_subvector(); unsupported vector type" );
        this->print ( );
        sub_vector->print ( );
        exit ( -1 );
    }

    return sub_vector;
}

template <typename ValueType>
void CPU_lVector<ValueType>::partial_replace_subvector ( const int start_i,
                                                         const int start_sub_vec,
                                                         const int size,
                                                         const lVector<ValueType> &sub_vec )
{

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &sub_vec ) )
    {

        this->partial_replace_subvector ( start_i,
                                          start_sub_vec,
                                          size,
                                          *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPU_lVector<ValueType>::partial_replace_subvector() unsupported vector type" );
        this->print ( );
        sub_vec.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::partial_replace_subvector ( const int start_i,
                                                         const int start_sub_vec,
                                                         const int size,
                                                         const CPU_lVector<ValueType> &sub_vec )
{
    assert ( start_i >= 0 );
    assert ( start_sub_vec >= 0 );
    assert ( size > 0 );
    assert ( start_sub_vec + size <= sub_vec.get_size ( ) );
    assert ( start_i + size <= this-> get_size ( ) );

    for ( int i = 0; i != size; ++i )
        this->buffer[start_i + i] = sub_vec.buffer[start_sub_vec + i];

}

template <typename ValueType>
void CPU_lVector<ValueType>::partial_add_subvector ( const int start_i,
                                                     const int start_sub_vec,
                                                     const int size,
                                                     const ValueType weight,
                                                     const lVector<ValueType> &sub_vec )
{
    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &sub_vec ) )
    {

        this->partial_add_subvector ( start_i,
                                      start_sub_vec,
                                      size,
                                      weight,
                                      *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPU_lVector<ValueType>::partial_add_subvector() unsupported vector type" );
        this->print ( );
        sub_vec.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::partial_add_subvector ( const int start_i,
                                                     const int start_sub_vec,
                                                     const int size,
                                                     const ValueType weight,
                                                     const CPU_lVector<ValueType> &sub_vec )
{
    assert ( start_i >= 0 );
    assert ( start_sub_vec >= 0 );
    assert ( size > 0 );
    assert ( start_sub_vec + size <= sub_vec.get_size ( ) );
    assert ( start_i + size <= this-> get_size ( ) );

    for ( int i = 0; i != size; ++i )
        this->buffer[start_i + i] += weight * sub_vec.buffer[start_sub_vec + i];

}

template <typename ValueType>
void CPU_lVector<ValueType>::ElementWiseMult ( const lVector<ValueType> &vec ) const
{

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        this->ElementWiseMult ( *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPU_lVector<ValueType>::ElementWiseMult() unsupported vector type" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::ElementWiseMult ( const CPU_lVector<ValueType> &vec ) const
{

    assert ( this->get_size ( ) == vec.get_size ( ) );

    for ( int i = 0, i_e = this->get_size ( ); i != i_e; ++i )
        this->buffer[i] *= vec.buffer[i];

}

template <typename ValueType>
void CPU_lVector<ValueType>::ReadFile ( const char* filename )
{
    std::string line;
    ValueType val;
    std::ifstream myfile;

    myfile.open ( filename, std::ifstream::in );

    if ( myfile.is_open ( ) )
    {

        int row = 0;

        // find number of lines = number of row
        while ( myfile.good ( ) )
        {
            myfile >> val;
            row++;
        }
        myfile.close ( );

        this->Init ( row - 1, filename );
        this->Zeros ( );
    }

    myfile.open ( filename );

    if ( myfile.is_open ( ) )
    {

        int row = 0;

        while ( myfile.good ( ) )
        {
            myfile >> val;

            this->buffer[row] = val;

            row++;

            if ( row == this->get_size ( ) )
                break;

        }

    }
    else
    {

        LOG_ERROR ( "CPU_lVector<ValueType>::ReadFile() Unable to open file: " << filename );
        this->print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void CPU_lVector<ValueType>::WriteFile ( const char* filename )
{
    std::string line;
    std::ofstream myfile;

    myfile.open ( filename, std::ifstream::out );

    if ( myfile.is_open ( ) )
    {

        for ( int row = 0, row_e = this->get_size ( ); row != row_e; ++row )
            myfile << std::scientific << this->buffer[row] << std::endl;

        myfile.close ( );

    }
    else
    {

        LOG_ERROR ( "CPU_lVector<ValueType>::WriteFile() Unable to open file: " << filename );
        this->print ( );
        exit ( -1 );
    }

}

template class CPU_lVector<double>;
template class CPU_lVector<float>;
