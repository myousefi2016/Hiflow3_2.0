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

#include "config.h"

#include "../lmatrix_coo_cpu.h"
#include "../lmp_log.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>

#ifdef WITH_OPENMP

#    include <omp.h>
#    include <sched.h>

#else

#    define ERROR LOG_ERROR("no OpenMP support");  exit(-1);

#endif

using namespace hiflow::la;

template <typename ValueType>
CPUopenmp_COO_lMatrix<ValueType>::CPUopenmp_COO_lMatrix ( int init_nnz, int init_num_row, int init_num_col, std::string init_name )
{
#ifdef WITH_OPENMP
    this->Init ( init_nnz, init_num_row, init_num_col, init_name );
    this->implementation_name_ = "parallel OpenMP";
    this->implementation_id_ = OPENMP;
    this->set_num_threads ( );
#else
    ERROR;
#endif
}

template <typename ValueType>
CPUopenmp_COO_lMatrix<ValueType>::CPUopenmp_COO_lMatrix ( )
{
#ifdef WITH_OPENMP
    this->implementation_name_ = "parallel OpenMP";
    this->implementation_id_ = OPENMP;
    this->set_num_threads ( );
#else
    ERROR;
#endif
}

template <typename ValueType>
CPUopenmp_COO_lMatrix<ValueType>::~CPUopenmp_COO_lMatrix ( )
{
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::CloneFrom ( const hiflow::la::lMatrix<ValueType> &other )
{
    if ( this != &other )
    {
        // if it is not empty init() will clean it
        this->Init ( other.get_nnz ( ), other.get_num_row ( ), other.get_num_col ( ), other.get_name ( ) );

        this->CopyStructureFrom ( other );

        this->CopyFrom ( other );
    }

    const CPUopenmp_COO_lMatrix<ValueType> * omp_other = dynamic_cast < const CPUopenmp_COO_lMatrix<ValueType> * > ( &other );
    if ( omp_other != 0 )
    {
        this->set_num_threads ( omp_other->num_threads ( ) );
    }
    else
    {
        const CPUmkl_COO_lMatrix<ValueType> * mkl_other = dynamic_cast < const CPUmkl_COO_lMatrix<ValueType> * > ( &other );
        if ( mkl_other != 0 )
        {
            this->set_num_threads ( mkl_other->num_threads ( ) );
        }
        // no else
    }
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::set_num_threads ( void )
{
#ifdef WITH_OPENMP
    // default value
    this->set_num_threads ( 8 );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::set_num_threads ( int num_thread )
{
#ifdef WITH_OPENMP
    this->num_threads_ = num_thread;
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::VectorMult ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
#ifdef WITH_OPENMP

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "ERROR CPUopenmp_COO_lMatrix<ValueType>::VectorMult unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMult ( *casted_invec, casted_outvec );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::VectorMult ( const CPU_lVector<ValueType> &invec, CPU_lVector<ValueType> *outvec ) const
{
    assert ( invec.get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec.get_size ( ) == this->get_num_col ( ) );
    assert ( outvec->get_size ( ) == this->get_num_row ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < this->get_num_row ( ); ++i )
        outvec->buffer[i] = 0;

#    pragma omp parallel for
    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        ValueType tmp = this->matrix.val[i] * invec.buffer[this->matrix.col[i]];

#    pragma omp atomic
        outvec->buffer[this->matrix.row[i]] += tmp;
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::VectorMultAdd ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
#ifdef WITH_OPENMP

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "ERROR CPUopenmp_COO_lMatrix<ValueType>::VectorMultAdd unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMultAdd ( *casted_invec, casted_outvec );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::VectorMultAdd ( const CPU_lVector<ValueType> &invec, CPU_lVector<ValueType> *outvec ) const
{
    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( ( outvec->get_size ( ) == this->get_num_row ( ) ) ||
             ( invec .get_size ( ) == 0 ) );

#ifdef WITH_OPENMP
    if ( this->get_nnz ( ) > 0 )
    {

        omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
        for ( int i = 0; i<this->get_nnz ( ); ++i )
        {
            ValueType tmp = this->matrix.val[i] * invec.buffer[this->matrix.col[i]];

#    pragma omp atomic
            outvec->buffer[this->matrix.row[i]] += tmp;
        }
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::CastFrom ( const CPU_COO_lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        this->matrix.val[i] = static_cast < ValueType > ( other.matrix.val[i] );
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::CastFrom ( const CPU_COO_lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        this->matrix.val[i] = static_cast < ValueType > ( other.matrix.val[i] );
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::CastTo ( CPU_COO_lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        other.matrix.val[i] = static_cast < double > ( this->matrix.val[i] );
    }
#else
    ERROR;
#endif

}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::CastTo ( CPU_COO_lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        other.matrix.val[i] = static_cast < float > ( this->matrix.val[i] );
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_COO_lMatrix<ValueType>::VectorMultNoDiag ( const CPU_lVector<ValueType>& in,
                                                          CPU_lVector<ValueType>* out ) const
{
    assert ( in.get_size ( ) == this->get_num_col ( ) );
    assert ( out->get_size ( ) == this->get_num_row ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < this->get_num_row ( ); ++i )
    {
        out->buffer[i] = 0.0;
    }

#    pragma omp parallel for
    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        if ( this->matrix.row[i] != this->matrix.col[i] )
        {
            out->buffer[this->matrix.row[i]] += this->matrix.val[i] * in.buffer[this->matrix.col[i]];
        }
    }
#else
    ERROR;
#endif
}

template class CPUopenmp_COO_lMatrix<float>;
template class CPUopenmp_COO_lMatrix<double>;
