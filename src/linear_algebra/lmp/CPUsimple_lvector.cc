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
#include "lmp_log.h"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <typeinfo>
#include <cmath>

using namespace hiflow::la;

template <typename ValueType>
CPUsimple_lVector<ValueType>::CPUsimple_lVector ( const int size, const std::string name )
{
    this->Init ( size, name );
    this->implementation_name_ = "sequential naive/simple";
    this->implementation_id_ = NAIVE;
}

template <typename ValueType>
CPUsimple_lVector<ValueType>::CPUsimple_lVector ( )
{
    this->implementation_name_ = "sequential naive/simple";
    this->implementation_id_ = NAIVE;
}

template <typename ValueType>
CPUsimple_lVector<ValueType>::~CPUsimple_lVector ( )
{
    this->Clear ( );
}

template <typename ValueType>
int CPUsimple_lVector<ValueType>::ArgMin ( void ) const
{
    int minimum = 0;

    assert ( this->get_size ( ) > 0 );

    ValueType val = std::abs ( this->buffer[0] );

    for ( int i = 1; i<this->get_size ( ); i++ )
    {
        const ValueType current = std::abs ( this->buffer[i] );
        if ( current < val )
        {
            minimum = i;
            val = current;
        }
    }

    return minimum;
}

template <typename ValueType>
int CPUsimple_lVector<ValueType>::ArgMax ( void ) const
{
    int maximum = 0;

    assert ( this->get_size ( ) > 0 );

    ValueType val = std::abs ( this->buffer[0] );

    for ( int i = 1; i<this->get_size ( ); i++ )
    {
        const ValueType current = std::abs ( this->buffer[i] );
        if ( current > val )
        {
            maximum = i;
            val = current;
        }
    }

    return maximum;
}

template <typename ValueType>
ValueType CPUsimple_lVector<ValueType>::Norm1 ( void ) const
{
    ValueType nrm1 = 0.0;

    const int N = this->get_size ( );

    assert ( N > 0 );

    /* Version that uses loop unrolling by an unroll-factor of 5*/
    ValueType ntemp = 0.0;

    // compute overhead to unroll factor
    const int M = N % 5;

    // if N is a multiple of 5
    if ( M == 0 )
    {
#pragma clang loop vectorize(enable)
        for ( int i = 0; i < N; i += 5 )
        {
            ntemp += std::abs ( this->buffer[i] ) + std::abs ( this->buffer[i + 1] ) + std::abs ( this->buffer[i + 2] ) + std::abs ( this->buffer[i + 3] ) + std::abs ( this->buffer[i + 4] );
        }
    }
    else
    {
        // result for overhead to unroll factor
#pragma clang loop vectorize(enable)
        for ( int i = 0; i != M; ++i )
        {
            ntemp += std::abs ( this->buffer[i] );
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#pragma clang loop vectorize(enable)
            for ( int i = M; i < N; i += 5 )
            {
                ntemp += std::abs ( this->buffer[i] ) + std::abs ( this->buffer[i + 1] ) + std::abs ( this->buffer[i + 2] ) + std::abs ( this->buffer[i + 3] ) + std::abs ( this->buffer[i + 4] );
            }
        }
    }
    nrm1 = ntemp;
    return nrm1;
}

template <typename ValueType>
ValueType CPUsimple_lVector<ValueType>::Norm2 ( void ) const
{
    ValueType norm = 0.0;

    assert ( this->get_size ( ) > 0 );

    norm = this->Dot ( *this );

    return sqrt ( norm );
}

template <typename ValueType>
ValueType CPUsimple_lVector<ValueType>::NormMax ( void ) const
{
    return std::abs ( this->buffer[this->ArgMax ( )] );
}

template <typename ValueType>
ValueType CPUsimple_lVector<ValueType>::Dot ( const CPU_lVector<ValueType> &vec ) const
{
    ValueType dot = 0.0;

    const int N = this->get_size ( );

    assert ( N > 0 );
    assert ( N == vec.get_size ( ) );

    /*for (int i=0, e = this->get_size(); i != e ; ++i)
        dot += this->buffer[i]*(vec.buffer)[i] ;*/

    /* Version that uses loop unrolling by an unroll-factor of 5*/
    ValueType dtemp = 0.0;

    // compute overhead to unroll factor
    const int M = N % 5;

    // if N is a multiple of 5
    if ( M == 0 )
    {
#pragma clang loop vectorize(enable)
        for ( int i = 0; i < N; i += 5 )
        {
            dtemp += this->buffer[i] * ( vec.buffer )[i] + this->buffer[i + 1] * ( vec.buffer )[i + 1] + this->buffer[i + 2] * ( vec.buffer )[i + 2] + this->buffer[i + 3] * ( vec.buffer )[i + 3] + this->buffer[i + 4] * ( vec.buffer )[i + 4];
        }
    }
    else
    {
        // result for overhead to unroll factor
#pragma clang loop vectorize(enable)
        for ( int i = 0; i != M; ++i )
        {
            dtemp += this->buffer[i] * ( vec.buffer )[i];
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#pragma clang loop vectorize(enable)
            for ( int i = M; i < N; i += 5 )
            {
                dtemp += this->buffer[i] * ( vec.buffer )[i] + this->buffer[i + 1] * ( vec.buffer )[i + 1] + this->buffer[i + 2] * ( vec.buffer )[i + 2] + this->buffer[i + 3] * ( vec.buffer )[i + 3] + this->buffer[i + 4] * ( vec.buffer )[i + 4];
            }
        }
    }
    dot = dtemp;
    return dot;
}

template <typename ValueType>
ValueType CPUsimple_lVector<ValueType>::Dot ( const lVector<ValueType> &vec ) const
{
    ValueType dot = 0.0;

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        dot = this->Dot ( *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPUsimple_lVector::Dot unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

    return dot;
}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Axpy ( const CPU_lVector<ValueType> &vec, const ValueType scalar )
{

    const int N = this->get_size ( );
    assert ( N > 0 );
    assert ( N == vec.get_size ( ) );

    /*for (int i=0, e = this->get_size(); i != e ; ++i)
      this->buffer[i] += scalar*(vec.buffer)[i] ;*/

    /* Version that uses loop unrolling by an unroll-factor of 5*/

    // compute overhead to unroll factor
    const int M = N % 5;

    // if N is a multiple of 5
    if ( M == 0 )
    {
#pragma clang loop vectorize(enable)
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
#pragma clang loop vectorize(enable)
        for ( int i = 0; i != M; ++i )
        {
            this->buffer[i] += scalar * ( vec.buffer )[i];
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#pragma clang loop vectorize(enable)
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

}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Axpy ( const lVector<ValueType> &vec, const ValueType scalar )
{

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        this->Axpy ( *casted_vec, scalar );

    }
    else
    {
        LOG_ERROR ( "CPUsimple_lVector::Axpy unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::ScaleAdd ( const ValueType scalar, const CPU_lVector<ValueType> &vec )
{

    /*assert(this->get_size() > 0);
    assert(this->get_size() == vec.get_size());

      for (int i=0, e = this->get_size(); i != e ; ++i)
        this->buffer[i] = scalar*(this->buffer)[i] + (vec.buffer)[i] ;*/

    const int N = this->get_size ( );
    assert ( N > 0 );
    assert ( N == vec.get_size ( ) );

    /* Version that uses loop unrolling by an unroll-factor of 5*/

    // compute overhead to unroll factor
    const int M = N % 5;

    // if N is a multiple of 5
    if ( M == 0 )
    {
#pragma clang loop vectorize(enable)
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
#pragma clang loop vectorize(enable)
        for ( int i = 0; i != M; ++i )
        {
            this->buffer[i] = scalar * ( this->buffer )[i] + ( vec.buffer )[i];
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#pragma clang loop vectorize(enable)
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
}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::ScaleAdd ( const ValueType scalar, const lVector<ValueType> &vec )
{

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        this->ScaleAdd ( scalar, *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPUsimple_lVector::ScaleAdd unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Scale ( const ValueType scalar )
{

    /*assert(this->get_size() > 0);

    for (int i=0, e = this->get_size(); i != e ; ++i)
      this->buffer[i] *= scalar ; */

    const int N = this->get_size ( );
    assert ( N > 0 );

    /* Version that uses loop unrolling by an unroll-factor of 5*/

    // compute overhead to unroll factor
    const int M = N % 5;

    // if N is a multiple of 5
    if ( M == 0 )
    {
#pragma clang loop vectorize(enable)
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
#pragma clang loop vectorize(enable)
        for ( int i = 0; i != M; ++i )
        {
            this->buffer[i] *= scalar;
        }

        // result for rest of vectors if length is greater than the unroll factor
        if ( N > 5 )
        {
#pragma clang loop vectorize(enable)
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

}

// Thanks to Benedikt Galler
//
// no complex values

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Rot ( CPU_lVector<ValueType> *vec, const ValueType &sc, const ValueType &ss )
{

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec->get_size ( ) );
    assert ( ss != 0 );
    assert ( sc != 1 );

    ValueType tmp;

    for ( size_t i = 0, e = this->get_size ( ); i != e; ++i )
    {
        tmp = sc * ( this->buffer )[i] + ss * ( vec->buffer )[i];
        ( vec->buffer )[i] = sc * ( vec->buffer )[i] - ss * ( this->buffer )[i];
        ( this->buffer )[i] = tmp;
    }

}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Rot ( lVector<ValueType> *vec, const ValueType &sc, const ValueType &ss )
{
    if ( CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < CPU_lVector<ValueType>* > ( vec ) )
    {

        this->Rot ( casted_vec, sc, ss );

    }
    else
    {
        LOG_ERROR ( "ERROR CPUsimple_lVector::Rot unsupported vectors" );
        this->print ( );
        vec->print ( );
        exit ( -1 );
    }

}

// Thanks to Benedikt Galler
//
// no complex values

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Rotg ( ValueType *sa, ValueType *sb, ValueType *sc, ValueType *ss ) const
{

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
}

// Thanks to Benedikt Galler
//
// zero-based sparam

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Rotm ( CPU_lVector<ValueType> *vec, const ValueType &sparam )
{

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

            for ( size_t i = 0, i_e = this->get_size ( ); i != i_e; i++ )
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

            for ( int i = 0, i_e = this->get_size ( ); i != i_e; i++ )
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

            for ( int i = 0, i_e = this->get_size ( ); i != i_e; i++ )
            {
                tmp_1 = ( vec->buffer )[i];
                tmp_2 = ( this->buffer )[i];
                ( vec->buffer )[i] = tmp_1 * h11 + tmp_2;
                ( this->buffer )[i] = tmp_2 * h22 - tmp_1;
            }

        }

    }

}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Rotm ( lVector<ValueType> *vec, const ValueType &sparam )
{
    if ( CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < CPU_lVector<ValueType>* > ( vec ) )
    {

        this->Rotm ( casted_vec, sparam );

    }
    else
    {
        LOG_ERROR ( "ERROR CPUsimple_lVector::Rotm unsupported vectors" );
        this->print ( );
        vec->print ( );
        exit ( -1 );

    }

}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::Rotmg ( ValueType *sd1, ValueType *sd2, ValueType *x1, const ValueType &x2, ValueType *sparam ) const
{

    LOG_ERROR ( "ERROR CPUsimple_lVector::Rotmg is not yet implemented" );
    this->print ( );
    exit ( -1 );

}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::CastFrom ( const CPU_lVector<double>& vec )
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#pragma clang loop vectorize(enable)
    for ( int i = 0; i < this->get_size ( ); ++i )
        ( this->buffer )[i] = static_cast < ValueType > ( vec.buffer[i] );
}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::CastFrom ( const CPU_lVector<float>& vec )
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#pragma clang loop vectorize(enable)
    for ( int i = 0; i < this->get_size ( ); ++i )
        ( this->buffer )[i] = static_cast < ValueType > ( vec.buffer[i] );
}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::CastTo ( CPU_lVector<double>& vec ) const
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#pragma clang loop vectorize(enable)
    for ( int i = 0; i < this->get_size ( ); ++i )
        vec.buffer[i] = static_cast < double > ( ( this->buffer )[i] );
}

template <typename ValueType>
void CPUsimple_lVector<ValueType>::CastTo ( CPU_lVector<float>& vec ) const
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

#pragma clang loop vectorize(enable)
    for ( int i = 0; i < this->get_size ( ); ++i )
        vec.buffer[i] = static_cast < float > ( ( this->buffer )[i] );
}

//
//
//
//

template class CPUsimple_lVector<double>;
template class CPUsimple_lVector<float>;
