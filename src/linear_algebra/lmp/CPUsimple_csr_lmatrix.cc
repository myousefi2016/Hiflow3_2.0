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

/// @author Dimitar Lukarski, Simon Gawlok

#include "la_global.h"
#include "lvector.h"
#include "lvector_cpu.h"
#include "lmatrix.h"
#include "lmatrix_csr.h"
#include "lmatrix_csr_cpu.h"
#include "lmp_log.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>

const int DEBUG_LEVEL = 2;

using namespace hiflow::la;

template <typename ValueType>
CPUsimple_CSR_lMatrix<ValueType>::CPUsimple_CSR_lMatrix ( int init_nnz, int init_num_row, int init_num_col, std::string init_name )
{
    this->Init ( init_nnz, init_num_row, init_num_col, init_name );
    this->implementation_name_ = "sequential naive/simple";
    this->implementation_id_ = NAIVE;
}

template <typename ValueType>
CPUsimple_CSR_lMatrix<ValueType>::CPUsimple_CSR_lMatrix ( )
{
    this->implementation_name_ = "sequential naive/simple";
    this->implementation_id_ = NAIVE;
}

template <typename ValueType>
CPUsimple_CSR_lMatrix<ValueType>::~CPUsimple_CSR_lMatrix ( )
{
    this->Clear ( );
}

template <typename ValueType>
void CPUsimple_CSR_lMatrix<ValueType>::VectorMult ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_CSR_lMatrix<ValueType>::VectorMult unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMult ( *casted_invec, casted_outvec );

}

template <typename ValueType>
void CPUsimple_CSR_lMatrix<ValueType>::VectorMult ( const CPU_lVector<ValueType> &invec, CPU_lVector<ValueType> *outvec ) const
{

    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( outvec->get_size ( ) == this->get_num_row ( ) );

    for ( size_t i = 0, e_i = this->get_num_row ( ); i != e_i; ++i )
    {
        //outvec->buffer[i] = 0 ;
        ValueType res = 0.;
        for ( size_t j = this->matrix.row[i], e_j = this->matrix.row[i + 1]; j != e_j; ++j )
        {
            /*
              LOG_DEBUG(2, " row " << j);
              LOG_DEBUG(2, " col " << this->matrix.col[j]);
              LOG_DEBUG(2, " matrix val " << this->matrix.val[j]);
              LOG_DEBUG(2, " vector size " << invec.get_size());
              LOG_DEBUG(2, " vector val " << invec.buffer[ this->matrix.col[j] ]);
             */
            if ( this->matrix.col[j] < 0 )
            {

                std::cout << j << " : " << this->matrix.col[j] << std::endl;
            }

            assert ( this->matrix.col[j] >= 0 );
            assert ( this->matrix.col[j] < invec.get_size ( ) );

            res += this->matrix.val[j] * invec.buffer[ this->matrix.col[j] ];

        }
        assert ( i < outvec->get_size ( ) );
        outvec->buffer[i] = res;
    }

}

template <typename ValueType>
void CPUsimple_CSR_lMatrix<ValueType>::CastFrom ( const CPU_CSR_lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    for ( int i = 0, e_i = this->get_nnz ( ); i != e_i; ++i )
    {
        this->matrix.val[i] = static_cast < ValueType > ( other.matrix.val[i] );
    }
}

template <typename ValueType>
void CPUsimple_CSR_lMatrix<ValueType>::CastFrom ( const CPU_CSR_lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    for ( int i = 0, e_i = this->get_nnz ( ); i != e_i; ++i )
    {
        this->matrix.val[i] = static_cast < ValueType > ( other.matrix.val[i] );
    }
}

template <typename ValueType>
void CPUsimple_CSR_lMatrix<ValueType>::CastTo ( CPU_CSR_lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    for ( int i = 0, e_i = this->get_nnz ( ); i != e_i; ++i )
    {
        other.matrix.val[i] = static_cast < double > ( this->matrix.val[i] );
    }
}

template <typename ValueType>
void CPUsimple_CSR_lMatrix<ValueType>::CastTo ( CPU_CSR_lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    for ( int i = 0, e_i = this->get_nnz ( ); i != e_i; ++i )
    {
        other.matrix.val[i] = static_cast < float > ( this->matrix.val[i] );
    }
}

template class CPUsimple_CSR_lMatrix<float>;
template class CPUsimple_CSR_lMatrix<double>;
