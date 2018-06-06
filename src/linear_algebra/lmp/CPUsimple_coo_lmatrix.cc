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

#include "lvector.h"
#include "lvector_cpu.h"
#include "lmatrix.h"
#include "lmatrix_coo.h"
#include "lmatrix_coo_cpu.h"
#include "lmp_log.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>

using namespace hiflow::la;

template <typename ValueType>
CPUsimple_COO_lMatrix<ValueType>::CPUsimple_COO_lMatrix ( int init_nnz, int init_num_row, int init_num_col, std::string init_name )
{
    this->Init ( init_nnz, init_num_row, init_num_col, init_name );
    this->implementation_name_ = "sequential naive/simple";
    this->implementation_id_ = NAIVE;
}

template <typename ValueType>
CPUsimple_COO_lMatrix<ValueType>::CPUsimple_COO_lMatrix ( )
{
    this->implementation_name_ = "sequential naive/simple";
    this->implementation_id_ = NAIVE;
}

template <typename ValueType>
CPUsimple_COO_lMatrix<ValueType>::~CPUsimple_COO_lMatrix ( )
{
}

template <typename ValueType>
void CPUsimple_COO_lMatrix<ValueType>::VectorMult ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_COO_lMatrix<ValueType>::VectorMult unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMult ( *casted_invec, casted_outvec );

}

template <typename ValueType>
void CPUsimple_COO_lMatrix<ValueType>::VectorMult ( const CPU_lVector<ValueType> &invec, CPU_lVector<ValueType> *outvec ) const
{

    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( outvec->get_size ( ) == this->get_num_row ( ) );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        outvec->buffer[i] = 0;

    for ( int i = 0, e = this->get_nnz ( ); i != e; ++i )
        outvec->buffer[this->matrix.row[i]] += this->matrix.val[i] * invec.buffer[this->matrix.col[i]];

}

template <typename ValueType>
void CPUsimple_COO_lMatrix<ValueType>::CastFrom ( const CPU_COO_lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        this->matrix.val[i] = static_cast < ValueType > ( other.matrix.val[i] );
    }
}

template <typename ValueType>
void CPUsimple_COO_lMatrix<ValueType>::CastFrom ( const CPU_COO_lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        this->matrix.val[i] = static_cast < ValueType > ( other.matrix.val[i] );
    }
}

template <typename ValueType>
void CPUsimple_COO_lMatrix<ValueType>::CastTo ( CPU_COO_lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        other.matrix.val[i] = static_cast < double > ( this->matrix.val[i] );
    }
}

template <typename ValueType>
void CPUsimple_COO_lMatrix<ValueType>::CastTo ( CPU_COO_lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        other.matrix.val[i] = static_cast < float > ( this->matrix.val[i] );
    }
}

template class CPUsimple_COO_lMatrix<float>;
template class CPUsimple_COO_lMatrix<double>;
