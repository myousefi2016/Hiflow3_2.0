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
#include "../lmatrix_csr_cpu.h"
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
CPUopenmp_CSR_lMatrix<ValueType>::CPUopenmp_CSR_lMatrix ( int init_nnz, int init_num_row, int init_num_col, std::string init_name )
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
CPUopenmp_CSR_lMatrix<ValueType>::CPUopenmp_CSR_lMatrix ( )
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
CPUopenmp_CSR_lMatrix<ValueType>::~CPUopenmp_CSR_lMatrix ( )
{
    this->Clear ( );
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::CloneFrom ( const hiflow::la::lMatrix<ValueType> &other )
{
    if ( this != &other )
    {
        // if it is not empty init() will clean it 
        this->Init ( other.get_nnz ( ), other.get_num_row ( ), other.get_num_col ( ), other.get_name ( ) );

        this->CopyStructureFrom ( other );

        this->CopyFrom ( other );
    }

    const CPUopenmp_CSR_lMatrix<ValueType> * omp_other = dynamic_cast < const CPUopenmp_CSR_lMatrix<ValueType> * > ( &other );
    if ( omp_other != 0 )
    {
        this->set_num_threads ( omp_other->num_threads ( ) );
    }
    else
    {
        const CPUmkl_CSR_lMatrix<ValueType> * mkl_other = dynamic_cast < const CPUmkl_CSR_lMatrix<ValueType> * > ( &other );
        if ( mkl_other != 0 )
        {
            this->set_num_threads ( mkl_other->num_threads ( ) );
        }
        // no else
    }
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::set_num_threads ( void )
{
#ifdef WITH_OPENMP
    // default value
    this->set_num_threads ( 8 );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::set_num_threads ( int num_thread )
{
#ifdef WITH_OPENMP
    this->num_threads_ = num_thread;
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::VectorMult ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
#ifdef WITH_OPENMP

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUopenmp_CSR_lMatrix<ValueType>::VectorMult unsupported in or out vector" );
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
void CPUopenmp_CSR_lMatrix<ValueType>::VectorMult ( const CPU_lVector<ValueType> &invec, CPU_lVector<ValueType> *outvec ) const
{
    assert ( invec.get_size ( ) > 0 );
    assert ( outvec->get_size ( ) > 0 );
    assert ( invec.get_size ( ) == this->get_num_col ( ) );
    assert ( outvec->get_size ( ) == this->get_num_row ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < this->get_num_row ( ); ++i )
    {
        outvec->buffer[i] = 0;
        for ( int j = this->matrix.row[i]; j != this->matrix.row[i + 1]; ++j )
            outvec->buffer[i] += this->matrix.val[j] * invec.buffer[ this->matrix.col[j] ];
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::VectorMultAdd ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
#ifdef WITH_OPENMP

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) || ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUopenmp_CSR_lMatrix<ValueType>::VectorMultAdd unsupported in or out vector" );
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
void CPUopenmp_CSR_lMatrix<ValueType>::VectorMultAdd ( const CPU_lVector<ValueType> &invec, CPU_lVector<ValueType> *outvec ) const
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
        for ( int i = 0; i<this->get_num_row ( ); ++i )
            for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
                outvec->buffer[i] += this->matrix.val[j] * invec.buffer[ this->matrix.col[j] ];

    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::VectorMultNoDiag ( const CPU_lVector<ValueType> &in, CPU_lVector<ValueType> *out ) const
{
    assert ( in.get_size ( ) == this->get_num_col ( ) );
    assert ( out->get_size ( ) == this->get_num_row ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < this->get_num_row ( ); ++i )
    {
        ValueType res = 0.0;
        for ( int j = this->matrix.row[i]; j < this->matrix.row[i + 1]; ++j )
        {
            if ( i != this->matrix.col[j] )
            {
                res += this->matrix.val[j] * in.buffer[this->matrix.col[j]];
            }
        }
        out->buffer[i] = res;
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::BlocksPsgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec,
                                                              const int num_blocks ) const
{
#ifdef WITH_OPENMP

    int start_i = 0;
    int end_i = 0;
    int step_i = invec.get_size ( ) / num_blocks;

    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for private(start_i,end_i)
    for ( int np = 0; np < num_blocks; ++np )
    {
        start_i = np*step_i;
        end_i = ( np + 1 ) * step_i;

        if ( np == num_blocks - 1 ) end_i = invec.get_size ( );

        this->BlockPsgauss_seidel ( invec, outvec, start_i, end_i );

    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::CastFrom ( const CPU_CSR_lMatrix<double>& other )
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
void CPUopenmp_CSR_lMatrix<ValueType>::get_add_values ( const int* rows,
                                                        int num_rows,
                                                        const int* cols,
                                                        int num_cols,
                                                        const int* cols_target,
                                                        int num_cols_target,
                                                        ValueType* values ) const
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < num_rows; ++i )
    {
        int k = 0;
        const int row_offset = i*num_cols_target;
        for ( int j = this->matrix.row[rows[i]]; k < num_cols && j<this->matrix.row[rows[i] + 1]; ++j )
        {
            if ( cols[k] == this->matrix.col[j] )
            {
                values[row_offset + cols_target[k]] += this->matrix.val[j];
                ++k;
            }
            else if ( cols[k] < this->matrix.col[j] )
            {
                while ( cols[k] < this->matrix.col[j] )
                {
                    ++k;
                }
                if ( cols[k] == this->matrix.col[j] )
                {
                    values[row_offset + cols_target[k]] += this->matrix.val[j];
                    ++k;
                }
            }
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::VectorMultAdd_submatrix ( const int* rows,
                                                                 int num_rows,
                                                                 const int* cols,
                                                                 int num_cols,
                                                                 const int* cols_input,
                                                                 const ValueType* in_values,
                                                                 ValueType* out_values ) const
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < num_rows; ++i )
    {
        int k = 0;
        for ( int j = this->matrix.row[rows[i]]; k < num_cols && j<this->matrix.row[rows[i] + 1]; ++j )
        {
            if ( cols[k] == this->matrix.col[j] )
            {
                out_values[i] += this->matrix.val[j] * in_values[cols_input[k]];
                ++k;
            }
            else if ( cols[k] < this->matrix.col[j] )
            {
                while ( cols[k] < this->matrix.col[j] )
                {
                    ++k;
                }
                if ( cols[k] == this->matrix.col[j] )
                {
                    out_values[i] += this->matrix.val[j] * in_values[cols_input[k]];
                    ++k;
                }
            }
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::VectorMultAdd_submatrix_vanka ( const int* rows,
                                                                       int num_rows,
                                                                       const CPU_lVector< ValueType > &invec,
                                                                       ValueType* out_values ) const
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < num_rows; ++i )
    {
        for ( int j = this->matrix.row[rows[i]]; j<this->matrix.row[rows[i] + 1]; ++j )
        {
            out_values[i] += this->matrix.val[j] * invec.buffer[this->matrix.col[j]];
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::VectorMultAdd_submatrix_vanka ( const int* rows,
                                                                       int num_rows,
                                                                       const hiflow::la::lVector< ValueType > &invec,
                                                                       ValueType* out_values ) const
{
#ifdef WITH_OPENMP
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );

    if ( casted_invec == NULL )
    {
        LOG_ERROR ( "CPUopenmp_CSR_lMatrix<ValueType>::VectorMultAdd_submatrix_vanka unsupported in vector" );
        this->print ( );
        invec.print ( );
        exit ( -1 );
    }

    this->VectorMultAdd_submatrix_vanka ( rows, num_rows, *casted_invec, out_values );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::add_values ( const int* rows, int num_rows, const int* cols, int num_cols, const ValueType* values )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for
    for ( int i = 0; i < num_rows; ++i )
    {
        int k = 0;
        const int row_offset = i*num_cols;
        for ( int j = this->matrix.row[rows[i]]; k < num_cols && j<this->matrix.row[rows[i] + 1]; ++j )
        {
            if ( cols[k] == this->matrix.col[j] )
            {
                this->matrix.val[j] += values[row_offset + k];
                ++k;
            }
            else if ( cols[k] < this->matrix.col[j] )
            {
                while ( cols[k] < this->matrix.col[j] && k < num_cols )
                {
                    ++k;
                }
                if ( k >= num_cols )
                {
                    break;
                }
                if ( cols[k] == this->matrix.col[j] )
                {
                    this->matrix.val[j] += values[row_offset + k];
                    ++k;
                }
            }
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_CSR_lMatrix<ValueType>::CastFrom ( const CPU_CSR_lMatrix<float>& other )
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
void CPUopenmp_CSR_lMatrix<ValueType>::CastTo ( CPU_CSR_lMatrix<double>& other ) const
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
void CPUopenmp_CSR_lMatrix<ValueType>::CastTo ( CPU_CSR_lMatrix<float>& other ) const
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

template class CPUopenmp_CSR_lMatrix<float>;
template class CPUopenmp_CSR_lMatrix<double>;
