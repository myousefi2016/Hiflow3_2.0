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

#include "config.h"

#include "../lvector_cpu.h"
#include "../lmp_log.h"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <typeinfo>
#include <cmath>

#ifdef WITH_OPENMP

#    include <omp.h>
#    include <sched.h>

#else

#    define ERROR LOG_ERROR("no OpenMP support");  exit(-1);

#endif

using namespace hiflow::la;

template <typename ValueType>
CPUopenmp_lVector<ValueType>::CPUopenmp_lVector ( const int size, const std::string name )
{
#ifdef WITH_OPENMP
    this->Init ( size, name );
    this->implementation_name_ = "parallel OpenMP";
    this->implementation_id_ = OPENMP;
    this->set_num_threads ( );
#else
    ERROR;
#endif
}

template <typename ValueType>
CPUopenmp_lVector<ValueType>::CPUopenmp_lVector ( )
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
void CPUopenmp_lVector<ValueType>::set_num_threads ( void )
{
#ifdef WITH_OPENMP
    // default value
    this->set_num_threads ( 8 );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::set_num_threads ( int num_thread )
{
#ifdef WITH_OPENMP
    this->num_threads_ = num_thread;
#else
    ERROR;
#endif
}

template <typename ValueType>
CPUopenmp_lVector<ValueType>::~CPUopenmp_lVector ( )
{
    this->Clear ( );
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::CloneFrom ( const hiflow::la::lVector<ValueType>& other )
{
    if ( this != &other )
    {
        this->Clear ( );
        this->CopyStructureFrom ( other );
        this->CopyFrom ( other );
    }

    const CPUopenmp_lVector<ValueType> * omp_other = dynamic_cast < const CPUopenmp_lVector<ValueType> * > ( &other );
    if ( omp_other != 0 )
    {
        this->set_num_threads ( omp_other->num_threads ( ) );
    }
    else
    {
        const CPUmkl_lVector<ValueType> * mkl_other = dynamic_cast < const CPUmkl_lVector<ValueType> * > ( &other );
        if ( mkl_other != 0 )
        {
            this->set_num_threads ( mkl_other->num_threads ( ) );
        }
        // no else
    }
}

template <typename ValueType>
int CPUopenmp_lVector<ValueType>::ArgMin ( void ) const
{
#ifdef WITH_OPENMP

    int minimum = 0;
    int n = this->get_size ( );
    int num_threads = this->num_threads_;

    assert ( this->get_size ( ) > 0 );

    omp_set_num_threads ( this->num_threads_ );

    std::vector<int> min ( num_threads );
    std::vector<ValueType> vals ( num_threads );

#    pragma omp parallel firstprivate(n,num_threads)
    {
        int start, end, id = omp_get_thread_num ( ), loc_min;
        start = n / num_threads;
        ValueType val;

        if ( id + 1 == num_threads )
        {
            end = n;
        }
        else
        {
            end = ( id + 1 ) * start;
        }

        start = id*start;

        loc_min = start;
        val = std::abs ( this->buffer[start] );

        for ( int i = start + 1; i < end; i++ )
            if ( std::abs ( this->buffer[i] ) < val )
            {
                loc_min = i;
                val = std::abs ( this->buffer[i] );
            }

        min[id] = loc_min;
        vals[id] = val;
    }

    minimum = min[0];
    ValueType val = vals[0];

    for ( int i = 1; i < num_threads; i++ )
        if ( vals[i] < val )
        {
            minimum = min[i];
            val = vals[i];
        }

    return minimum;

#else
    ERROR;
    return -1;
#endif
}

template <typename ValueType>
int CPUopenmp_lVector<ValueType>::ArgMax ( void ) const
{
#ifdef WITH_OPENMP

    int maximum = 0;
    int n = this->get_size ( );
    int num_threads = this->num_threads_;

    assert ( this->get_size ( ) > 0 );

    omp_set_num_threads ( this->num_threads_ );

    std::vector<int> max ( num_threads );
    std::vector<ValueType> vals ( num_threads );

#    pragma omp parallel firstprivate(n,num_threads)
    {
        int start, end, id = omp_get_thread_num ( ), loc_max;
        start = n / num_threads;
        ValueType val;

        if ( id + 1 == num_threads )
        {
            end = n;
        }
        else
        {
            end = ( id + 1 ) * start;
        }

        start = id*start;

        loc_max = start;
        val = std::abs ( this->buffer[start] );

        for ( int i = start + 1; i < end; i++ )
            if ( fabs ( this->buffer[i] ) > val )
            {
                loc_max = i;
                val = std::abs ( this->buffer[i] );
            }

        max[id] = loc_max;
        vals[id] = val;
    }

    maximum = max[0];
    ValueType val = vals[0];

    for ( int i = 1; i < num_threads; i++ )
        if ( vals[i] > val )
        {
            maximum = max[i];
            val = vals[i];
        }

    return maximum;

#else
    ERROR;
    return -1;
#endif
}

template <typename ValueType>
ValueType CPUopenmp_lVector<ValueType>::Norm1 ( void ) const
{
#ifdef WITH_OPENMP
    ValueType nrm1 = 0.0;

    const int N = this->get_size ( );

    assert ( N > 0 );

    /* Version that uses loop unrolling by an unroll-factor of 5*/
    ValueType ntemp = 0.0;

    // compute overhead to unroll factor
    const int M = N % 5;

    omp_set_num_threads ( this->num_threads_ );

    // if N is a multiple of 5
    if ( M == 0 )
    {
#    pragma omp simd reduction(+:ntemp)
        for ( int i = 0; i < N; i += 5 )
        {
            ntemp += std::abs ( this->buffer[i] ) + std::abs ( this->buffer[i + 1] ) + std::abs ( this->buffer[i + 2] ) + std::abs ( this->buffer[i + 3] ) + std::abs ( this->buffer[i + 4] );
        }
    }
    else
    {
        // result for overhead to unroll factor
        for ( int i = 0; i < M; ++i )
        {
            ntemp += std::abs ( this->buffer[i] );
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#    pragma omp simd reduction(+:ntemp)
            for ( int i = M; i < N; i += 5 )
            {
                ntemp += std::abs ( this->buffer[i] ) + std::abs ( this->buffer[i + 1] ) + std::abs ( this->buffer[i + 2] ) + std::abs ( this->buffer[i + 3] ) + std::abs ( this->buffer[i + 4] );
            }
        }
    }
    nrm1 = ntemp;
    return nrm1;
#else
    ERROR;
    return -1.0;
#endif
}

template <typename ValueType>
ValueType CPUopenmp_lVector<ValueType>::Norm2 ( void ) const
{
#ifdef WITH_OPENMP
    ValueType norm = 0.0;

    assert ( this->get_size ( ) > 0 );

    norm = this->Dot ( *this );

    return sqrt ( norm );

#else
    ERROR;
    return -1.0;
#endif
}

template <typename ValueType>
ValueType CPUopenmp_lVector<ValueType>::NormMax ( void ) const
{
#ifdef WITH_OPENMP
    return std::abs ( this->buffer[this->ArgMax ( )] );
#else
    ERROR;
    return -1.0;
#endif
}

template <typename ValueType>
ValueType CPUopenmp_lVector<ValueType>::Dot ( const CPU_lVector<ValueType> &vec ) const
{
#ifdef WITH_OPENMP
    ValueType dot = 0.0;

    const int N = this->get_size ( );

    assert ( N >= 0 );
    assert ( N == vec.get_size ( ) );

    /*for (int i=0, e = this->get_size(); i != e ; ++i)
    dot += *this->buffer[i]*(vec.buffer)[i] ;*/

    /* Version that uses loop unrolling by an unroll-factor of 5*/
    ValueType dtemp = 0.0;

    // compute overhead to unroll factor
    const int M = N % 5;

    omp_set_num_threads ( this->num_threads_ );

    // if N is a multiple of 5
    if ( M == 0 )
    {
#    pragma omp simd reduction(+:dtemp)
        for ( int i = 0; i < N; i += 5 )
        {
            dtemp += this->buffer[i ] * ( vec.buffer )[i ]
                    + this->buffer[i + 1] * ( vec.buffer )[i + 1]
                    + this->buffer[i + 2] * ( vec.buffer )[i + 2]
                    + this->buffer[i + 3] * ( vec.buffer )[i + 3]
                    + this->buffer[i + 4] * ( vec.buffer )[i + 4];
        }
    }
    else
    {
        // result for overhead to unroll factor
        for ( int i = 0; i != M; ++i )
        {
            dtemp += this->buffer[i] * ( vec.buffer )[i];
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#    pragma omp simd reduction(+:dtemp)
            for ( int i = M; i < N; i += 5 )
            {
                dtemp += this->buffer[i ] * ( vec.buffer )[i ]
                        + this->buffer[i + 1] * ( vec.buffer )[i + 1]
                        + this->buffer[i + 2] * ( vec.buffer )[i + 2]
                        + this->buffer[i + 3] * ( vec.buffer )[i + 3]
                        + this->buffer[i + 4] * ( vec.buffer )[i + 4];
            }
        }
    }

    dot = dtemp;
    return dot;
#else
    ERROR;
    return -1.0;
#endif
}

template <typename ValueType>
ValueType CPUopenmp_lVector<ValueType>::Dot ( const lVector<ValueType> &vec ) const
{
#ifdef WITH_OPENMP
    ValueType dot = 0.0;

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        dot = this->Dot ( *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPUopenmp_lVector::Dot unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

    return dot;
#else
    ERROR;
    return -1.0;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Axpy ( const CPU_lVector<ValueType> &vec, const ValueType scalar )
{
#ifdef WITH_OPENMP

    const int N = this->get_size ( );
    assert ( N > 0 );
    assert ( N == vec.get_size ( ) );

    /*for (int i=0, e = this->get_size(); i != e ; ++i)
    t*his->buffer[i] += scalar*(vec.buffer)[i] ;*/

    /* Version that uses loop unrolling by an unroll-factor of 5*/

    // compute overhead to unroll factor
    const int M = N % 5;

    omp_set_num_threads ( this->num_threads_ );

    // if N is a multiple of 5
    if ( M == 0 )
    {
#    pragma omp simd
        for ( int i = 0; i < N; i += 5 )
        {
            this->buffer[i] += scalar * ( vec.buffer )[i];
            this->buffer[i + 1] += scalar * ( vec.buffer )[i + 1];
            this->buffer[i + 2] += scalar * ( vec.buffer )[i + 2];
            this->buffer[i + 3] += scalar * ( vec.buffer )[i + 3];
            this->buffer[i + 4] += scalar * ( vec.buffer )[i + 4];
        }
    }
    else
    {
        // result for overhead to unroll factor
        for ( int i = 0; i < M; ++i )
        {
            this->buffer[i] += scalar * ( vec.buffer )[i];
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#    pragma omp simd
            for ( int i = M; i < N; i += 5 )
            {
                this->buffer[i] += scalar * ( vec.buffer )[i];
                this->buffer[i + 1] += scalar * ( vec.buffer )[i + 1];
                this->buffer[i + 2] += scalar * ( vec.buffer )[i + 2];
                this->buffer[i + 3] += scalar * ( vec.buffer )[i + 3];
                this->buffer[i + 4] += scalar * ( vec.buffer )[i + 4];
            }
        }
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Axpy ( const lVector<ValueType> &vec, const ValueType scalar )
{
#ifdef WITH_OPENMP
    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        this->Axpy ( *casted_vec, scalar );

    }
    else
    {
        LOG_ERROR ( "CPUopenmp_lVector::Axpy unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::ScaleAdd ( const ValueType scalar, const CPU_lVector<ValueType> &vec )
{
#ifdef WITH_OPENMP

    const int N = this->get_size ( );
    assert ( N > 0 );
    assert ( N == vec.get_size ( ) );

    /* Version that uses loop unrolling by an unroll-factor of 5*/

    // compute overhead to unroll factor
    const int M = N % 5;

    omp_set_num_threads ( this->num_threads_ );

    // if N is a multiple of 5
    if ( M == 0 )
    {
#    pragma omp simd
        for ( int i = 0; i < N; i += 5 )
        {
            this->buffer[i] = scalar * ( this->buffer )[i] + ( vec.buffer )[i];
            this->buffer[i + 1] = scalar * ( this->buffer )[i + 1] + ( vec.buffer )[i + 1];
            this->buffer[i + 2] = scalar * ( this->buffer )[i + 2] + ( vec.buffer )[i + 2];
            this->buffer[i + 3] = scalar * ( this->buffer )[i + 3] + ( vec.buffer )[i + 3];
            this->buffer[i + 4] = scalar * ( this->buffer )[i + 4] + ( vec.buffer )[i + 4];
        }
    }
    else
    {
        // result for overhead to unroll factor
        for ( int i = 0; i < M; ++i )
        {
            this->buffer[i] = scalar * ( this->buffer )[i] + ( vec.buffer )[i];
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#    pragma omp simd
            for ( int i = M; i < N; i += 5 )
            {
                this->buffer[i] = scalar * ( this->buffer )[i] + ( vec.buffer )[i];
                this->buffer[i + 1] = scalar * ( this->buffer )[i + 1] + ( vec.buffer )[i + 1];
                this->buffer[i + 2] = scalar * ( this->buffer )[i + 2] + ( vec.buffer )[i + 2];
                this->buffer[i + 3] = scalar * ( this->buffer )[i + 3] + ( vec.buffer )[i + 3];
                this->buffer[i + 4] = scalar * ( this->buffer )[i + 4] + ( vec.buffer )[i + 4];
            }
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::ScaleAdd ( const ValueType scalar, const lVector<ValueType> &vec )
{
#ifdef WITH_OPENMP

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        this->ScaleAdd ( scalar, *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPUopenmp_lVector::ScaleAdd unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Scale ( const ValueType scalar )
{
#ifdef WITH_OPENMP

    const int N = this->get_size ( );
    assert ( N > 0 );

    /* Version that uses loop unrolling by an unroll-factor of 5*/

    // compute overhead to unroll factor
    const int M = N % 5;

    omp_set_num_threads ( this->num_threads_ );

    // if N is a multiple of 5
    if ( M == 0 )
    {
#    pragma omp simd
        for ( int i = 0; i < N; i += 5 )
        {
            this->buffer[i] *= scalar;
            this->buffer[i + 1] *= scalar;
            this->buffer[i + 2] *= scalar;
            this->buffer[i + 3] *= scalar;
            this->buffer[i + 4] *= scalar;
        }
    }
    else
    {
        // result for overhead to unroll factor
        for ( int i = 0; i < M; ++i )
        {
            this->buffer[i] *= scalar;
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#    pragma omp simd
            for ( int i = M; i < N; i += 5 )
            {
                this->buffer[i] *= scalar;
                this->buffer[i + 1] *= scalar;
                this->buffer[i + 2] *= scalar;
                this->buffer[i + 3] *= scalar;
                this->buffer[i + 4] *= scalar;
            }
        }
    }

#else
    ERROR;
#endif
}

// Thanks to Benedikt Galler

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Rot ( CPU_lVector<ValueType> *vec, const ValueType &sc, const ValueType &ss )
{
#ifdef WITH_OPENMP

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec->get_size ( ) );
    assert ( ss != 0 );
    assert ( sc != 1 );

    omp_set_num_threads ( this->num_threads_ );

    ValueType tmp;
    ValueType sin = ss;
    ValueType cos = sc;

    omp_set_num_threads ( this->num_threads_ );

#    pragma omp parallel for private(tmp) firstprivate(sin,cos)
    for ( int i = 0; i < this->get_size ( ); i++ )
    {
        tmp = cos * ( this->buffer )[i] + sin * ( vec->buffer )[i];
        ( vec->buffer )[i] = cos * ( vec->buffer )[i] - sin * ( this->buffer )[i];
        ( this->buffer )[i] = tmp;
    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Rot ( lVector<ValueType> *vec, const ValueType &sc, const ValueType &ss )
{
#ifdef WITH_OPENMP

    if ( CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < CPU_lVector<ValueType>* > ( vec ) )
    {

        this->Rot ( casted_vec, sc, ss );

    }
    else
    {
        LOG_ERROR ( "CPUopenmp_lVector::Rot unsupported vectors" );
        this->print ( );
        vec->print ( );
        exit ( -1 );
    }
#else
    ERROR;
#endif
}

// Thanks to Benedikt Galler

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Rotg ( ValueType *sa, ValueType *sb, ValueType *sc, ValueType *ss ) const
{
#ifdef WITH_OPENMP

    omp_set_num_threads ( this->num_threads_ );

    ValueType roe, scal, r, z, tmp_a, tmp_b, t0, t1;

    tmp_a = fabs ( *sa );
    tmp_b = fabs ( *sb );

    if ( tmp_a > tmp_b )
    {
        roe = *sa;
    }
    else
    {
        roe = *sb;
    }

    scal = tmp_a + tmp_b;

    if ( scal == 0 )
    {
        *sc = 1;
        *ss = *sa = *sb = 0;
    }
    else
    {
        t0 = tmp_a / scal;
        t1 = tmp_b / scal;
        r = scal * sqrt ( t0 * t0 + t1 * t1 );

        if ( roe < 0 ) r = -r;
        *sc = *sa / r;
        *ss = *sb / r;

        if ( tmp_a > tmp_b ) z = *ss;
        else if ( *sc != 0 ) z = 1 / *sc;
        else z = 1;

        *sa = r;
        *sb = z;
    }

#else
    ERROR;
#endif
}

// Thanks to Benedikt Galler

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Rotm ( CPU_lVector<ValueType> *vec, const ValueType &sparam )
{
#ifdef WITH_OPENMP

    omp_set_num_threads ( this->num_threads_ );

    ValueType flag = ( &sparam )[0];
    ValueType h11, h21, h12, h22, tmp_1, tmp_2;

    assert ( this->get_size ( ) > 0 );

    if ( flag != -2 )
    {

        if ( flag == -1 )
        {
            h11 = ( &sparam )[1];
            h21 = ( &sparam )[2];
            h12 = ( &sparam )[3];
            h22 = ( &sparam )[4];

#    pragma omp parallel for private(tmp_1,tmp_2) firstprivate(h11,h21,h12,h22)
            for ( int i = 0; i< this->get_size ( ); i++ )
            {
                tmp_1 = ( vec->buffer )[i];
                tmp_2 = ( this->buffer )[i];
                ( vec->buffer )[i] = tmp_1 * h11 + tmp_2 * h12;
                ( this->buffer )[i] = tmp_1 * h21 + tmp_2 * h22;
            }
        }

        if ( flag == 0 )
        {
            h21 = ( &sparam )[2];
            h12 = ( &sparam )[3];

#    pragma omp parallel for private(tmp_1,tmp_2) firstprivate(h21,h12)
            for ( int i = 0; i< this->get_size ( ); i++ )
            {
                tmp_1 = ( vec->buffer )[i];
                tmp_2 = ( this->buffer )[i];
                ( vec->buffer )[i] = tmp_1 + tmp_2 * h12;
                ( this->buffer )[i] = tmp_1 * h21 + tmp_2;
            }

        }

        if ( flag == 1 )
        {
            h11 = ( &sparam )[1];
            h22 = ( &sparam )[4];

#    pragma omp parallel for private(tmp_1,tmp_2) firstprivate(h11,h22)
            for ( int i = 0; i< this->get_size ( ); i++ )
            {
                tmp_1 = ( vec->buffer )[i];
                tmp_2 = ( this->buffer )[i];
                ( vec->buffer )[i] = tmp_1 * h11 + tmp_2;
                ( this->buffer )[i] = tmp_2 * h22 - tmp_1;
            }

        }

    }

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Rotm ( lVector<ValueType> *vec, const ValueType &sparam )
{
#ifdef WITH_OPENMP

    if ( CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < CPU_lVector<ValueType>* > ( vec ) )
    {

        this->Rotm ( casted_vec, sparam );

    }
    else
    {
        LOG_ERROR ( "CPUopenmp_lVector::Rotm unsupported vectors" );
        this->print ( );
        vec->print ( );
        exit ( -1 );

    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::Rotmg ( ValueType *sd1, ValueType *sd2, ValueType *x1, const ValueType &x2, ValueType *sparam ) const
{

    LOG_ERROR ( "ERROR CPUopenmp_lVector::Rotmg is not yet implemented" );
    this->print ( );
    exit ( -1 );

}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::ElementWiseMult ( const CPU_lVector<ValueType> &vec ) const
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp simd
    for ( int i = 0; i<this->get_size ( ); ++i )
        this->buffer[i] *= vec.buffer[i];
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::CastFrom ( const CPU_lVector<double>& vec )
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp simd
    for ( int i = 0; i < this->get_size ( ); ++i )
        ( this->buffer )[i] = static_cast < ValueType > ( vec.buffer[i] );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::CastFrom ( const CPU_lVector<float>& vec )
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp simd
    for ( int i = 0; i < this->get_size ( ); ++i )
        ( this->buffer )[i] = static_cast < ValueType > ( vec.buffer[i] );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::CastTo ( CPU_lVector<double>& vec ) const
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp simd
    for ( int i = 0; i < this->get_size ( ); ++i )
        vec.buffer[i] = static_cast < double > ( ( this->buffer )[i] );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUopenmp_lVector<ValueType>::CastTo ( CPU_lVector<float>& vec ) const
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#ifdef WITH_OPENMP
    omp_set_num_threads ( this->num_threads_ );

#    pragma omp simd
    for ( int i = 0; i < this->get_size ( ); ++i )
        vec.buffer[i] = static_cast < float > ( ( this->buffer )[i] );
#else
    ERROR;
#endif
}

//
//
//
//

template class CPUopenmp_lVector<double>;
template class CPUopenmp_lVector<float>;
