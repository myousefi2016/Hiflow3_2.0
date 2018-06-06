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

#include "../lvector_cpu.h"
#include "../lmp_log.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>

#ifdef WITH_MKL

#    include <mkl.h>
#    include <mkl_cblas.h>

#else

#    define ERROR LOG_ERROR("no Intel MKL support");  exit(-1);

#endif

using namespace hiflow::la;

template <typename ValueType>
CPUmkl_lVector<ValueType>::CPUmkl_lVector ( const int size, const std::string name )
{
#ifdef WITH_MKL
    this->Init ( size, name );
    this->implementation_name_ = "Intel/MKL";
    this->implementation_id_ = MKL;
    this->set_num_threads ( );
#else
    ERROR;
#endif
}

template <typename ValueType>
CPUmkl_lVector<ValueType>::CPUmkl_lVector ( )
{
#ifdef WITH_MKL
    this->implementation_name_ = "Intel/MKL";
    this->implementation_id_ = MKL;
    this->set_num_threads ( );
#else
    ERROR;
#endif
}

template <typename ValueType>
CPUmkl_lVector<ValueType>::~CPUmkl_lVector ( )
{
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::CloneFrom ( const hiflow::la::lVector<ValueType>& other )
{
    if ( this != &other )
    {
        this->Clear ( );
        this->CopyStructureFrom ( other );
        this->CopyFrom ( other );
    }

    const CPUmkl_lVector<ValueType> * mkl_other = dynamic_cast < const CPUmkl_lVector<ValueType> * > ( &other );
    if ( mkl_other != 0 )
    {
        this->set_num_threads ( mkl_other->num_threads ( ) );
    }
    else
    {
        const CPUopenmp_lVector<ValueType> * omp_other = dynamic_cast < const CPUopenmp_lVector<ValueType> * > ( &other );
        if ( omp_other != 0 )
        {
            this->set_num_threads ( omp_other->num_threads ( ) );
        }
        // no else
    }
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::set_num_threads ( void )
{
#ifdef WITH_MKL
    // default value
    this->set_num_threads ( 1 );
#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::set_num_threads ( int num_thread )
{
#ifdef WITH_MKL
    this->num_threads_ = num_thread;
#else
    ERROR;
#endif
}

template <>
int CPUmkl_lVector<double>::ArgMin ( void ) const
{
#ifdef WITH_MKL
    int minimum = -1;

    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    minimum = cblas_idamin ( this->get_size ( ), this->buffer, 1 );

    return minimum;
#else
    ERROR;
    return -1;
#endif
}

template <>
int CPUmkl_lVector<float>::ArgMin ( void ) const
{
#ifdef WITH_MKL
    int minimum = -1;

    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    minimum = cblas_isamin ( this->get_size ( ), this->buffer, 1 );

    return minimum;
#else
    ERROR;
    return -1;
#endif
}

template <>
int CPUmkl_lVector<double>::ArgMax ( void ) const
{
#ifdef WITH_MKL
    int maximum = -1;

    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    maximum = cblas_idamax ( this->get_size ( ), this->buffer, 1 );

    return maximum;
#else
    ERROR;
    return -1;
#endif
}

template <>
int CPUmkl_lVector<float>::ArgMax ( void ) const
{
#ifdef WITH_MKL
    int maximum = -1;

    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    maximum = cblas_isamax ( this->get_size ( ), this->buffer, 1 );

    return maximum;
#else
    ERROR;
    return -1;
#endif
}

template <>
double CPUmkl_lVector<double>::Norm1 ( void ) const
{
#ifdef WITH_MKL
    double nrm1 = -1;

    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    nrm1 = cblas_dasum ( this->get_size ( ), this->buffer, 1 );

    return nrm1;
#else
    ERROR;
    return -1.0;
#endif
}

template <>
float CPUmkl_lVector<float>::Norm1 ( void ) const
{
#ifdef WITH_MKL
    float nrm1 = -1;

    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    nrm1 = cblas_sasum ( this->get_size ( ), this->buffer, 1 );

    return nrm1;
#else
    ERROR;
    return -1.0;
#endif
}

template <>
double CPUmkl_lVector<double>::Norm2 ( void ) const
{
#ifdef WITH_MKL
    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    double dot = cblas_ddot ( this->get_size ( ), this->buffer, 1, this->buffer, 1 );

    return sqrt ( dot );
#else
    ERROR;
    return -1.0;
#endif
}

template <>
float CPUmkl_lVector<float>::Norm2 ( void ) const
{
#ifdef WITH_MKL
    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    float dot = cblas_sdot ( this->get_size ( ), this->buffer, 1, this->buffer, 1 );

    return sqrt ( dot );
#else
    ERROR;
    return -1.0;
#endif
}

template <typename ValueType>
ValueType CPUmkl_lVector<ValueType>::NormMax ( void ) const
{
#ifdef WITH_MKL
    return std::abs ( this->buffer[this->ArgMax ( )] );
#else
    ERROR;
    return -1.0;
#endif
}

template <typename ValueType>
ValueType CPUmkl_lVector<ValueType>::Dot ( const lVector<ValueType> &vec ) const
{
    ValueType dot = 0.0;

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        dot = this->Dot ( *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPUmkl_lVector::Dot unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

    return dot;
}

template <>
double CPUmkl_lVector<double>::Dot ( const CPU_lVector<double> &vec ) const
{
#ifdef WITH_MKL
    double dot = 0.0;

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec.get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    dot = cblas_ddot ( this->get_size ( ), this->buffer, 1, vec.buffer, 1 );

    return dot;
#else
    ERROR;
    return -1.0;
#endif
}

template <>
float CPUmkl_lVector<float>::Dot ( const CPU_lVector<float> &vec ) const
{
#ifdef WITH_MKL
    float dot = 0.0;

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec.get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    dot = cblas_sdot ( this->get_size ( ), this->buffer, 1, vec.buffer, 1 );

    return dot;
#else
    ERROR;
    return -1.0;
#endif
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::Axpy ( const lVector<ValueType> &vec, const ValueType scalar )
{
    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        this->Axpy ( *casted_vec, scalar );

    }
    else
    {
        LOG_ERROR ( "CPUmkl_lVector::Axpy unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

}

template <>
void CPUmkl_lVector<double>::Axpy ( const CPU_lVector<double> &vec, const double scalar )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec.get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_daxpy ( this->get_size ( ), scalar, vec.buffer, 1, this->buffer, 1 );

#else
    ERROR;
#endif
}

template <>
void CPUmkl_lVector<float>::Axpy ( const CPU_lVector<float> &vec, const float scalar )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec.get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_saxpy ( this->get_size ( ), scalar, vec.buffer, 1, this->buffer, 1 );

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::ScaleAdd ( const ValueType scalar, const lVector<ValueType> &vec )
{

    if ( const CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < const CPU_lVector<ValueType>* > ( &vec ) )
    {

        this->ScaleAdd ( scalar, *casted_vec );

    }
    else
    {
        LOG_ERROR ( "CPUmkl_lVector::ScaleAdd unsupported vectors" );
        this->print ( );
        vec.print ( );
        exit ( -1 );
    }

}

template <>
void CPUmkl_lVector<double>::ScaleAdd ( const double scalar, const CPU_lVector<double> &vec )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec.get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_dscal ( this->get_size ( ), scalar, this->buffer, 1 );
    cblas_daxpy ( this->get_size ( ), ( double ) ( 1.0 ), vec.buffer, 1, this->buffer, 1 );

#else
    ERROR;
#endif
}

template <>
void CPUmkl_lVector<float>::ScaleAdd ( const float scalar, const CPU_lVector<float> &vec )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec.get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_sscal ( this->get_size ( ), scalar, this->buffer, 1 );
    cblas_saxpy ( this->get_size ( ), ( float ) ( 1.0 ), vec.buffer, 1, this->buffer, 1 );

#else
    ERROR;
#endif
}

template <>
void CPUmkl_lVector<double>::Scale ( const double scalar )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_dscal ( this->get_size ( ), scalar, this->buffer, 1 );

#else
    ERROR;
#endif
}

template <>
void CPUmkl_lVector<float>::Scale ( const float scalar )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_sscal ( this->get_size ( ), scalar, this->buffer, 1 );

#else
    ERROR;
#endif
}

// CPUmkl_lVector::Rot()
// Rot

template <typename ValueType>
void CPUmkl_lVector<ValueType>::Rot ( lVector<ValueType> *vec, const ValueType &sc, const ValueType &ss )
{

    if ( CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < CPU_lVector<ValueType>* > ( vec ) )
    {

        this->Rot ( casted_vec, sc, ss );

    }
    else
    {
        LOG_ERROR ( "CPUmkl_lVector::Rot unsupported vectors" );
        this->print ( );
        vec->print ( );
        exit ( -1 );
    }

}

// CPUmkl_lVector::Rot()
// Rot

template <>
void CPUmkl_lVector<double>::Rot ( CPU_lVector<double> *vec, const double &sc, const double &ss )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec->get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_drot ( this->get_size ( ), vec->buffer, 1, this->buffer, 1, sc, ss );

#else
    ERROR;
#endif
}

// CPUmkl_lVector::Rot()
// Rot

template <>
void CPUmkl_lVector<float>::Rot ( CPU_lVector<float> *vec, const float &sc, const float &ss )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec->get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_srot ( this->get_size ( ), vec->buffer, 1, this->buffer, 1, sc, ss );

#else
    ERROR;
#endif
}

// CPUmkl_lVector::Rotg()
// Rotg

template <>
void CPUmkl_lVector<double>::Rotg ( double *sa, double *sb, double *sc, double *ss ) const
{
#ifdef WITH_MKL

    mkl_set_num_threads ( this->num_threads_ );
    cblas_drotg ( sa, sb, sc, ss );

#else
    ERROR;
#endif
}

// CPUmkl_lVector::Rotg()
// Rotg

template <>
void CPUmkl_lVector<float>::Rotg ( float *sa, float *sb, float *sc, float *ss ) const
{
#ifdef WITH_MKL

    mkl_set_num_threads ( this->num_threads_ );
    cblas_srotg ( sa, sb, sc, ss );

#else
    ERROR;
#endif
}

// CPUmkl_lVector::Rotm()
// Rotm

template <typename ValueType>
void CPUmkl_lVector<ValueType>::Rotm ( lVector<ValueType> *vec, const ValueType &sparam )
{

    if ( CPU_lVector<ValueType> *casted_vec =
         dynamic_cast < CPU_lVector<ValueType>* > ( vec ) )
    {

        this->Rotm ( casted_vec, sparam );

    }
    else
    {
        LOG_ERROR ( "CPUmkl_lVector::Rotm unsupported vectors" );
        this->print ( );
        vec->print ( );
        exit ( -1 );
    }

}

// CPUmkl_lVector::Rotm()
// Rotm

template <>
void CPUmkl_lVector<double>::Rotm ( CPU_lVector<double> *vec, const double &sparam )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec->get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_drotm ( this->get_size ( ), vec->buffer, 1, this->buffer, 1, &sparam );

#else
    ERROR;
#endif
}

// CPUmkl_lVector::Rotm()
// Rotm

template <>
void CPUmkl_lVector<float>::Rotm ( CPU_lVector<float> *vec, const float &sparam )
{
#ifdef WITH_MKL

    assert ( this->get_size ( ) > 0 );
    assert ( this->get_size ( ) == vec->get_size ( ) );

    mkl_set_num_threads ( this->num_threads_ );
    cblas_srotm ( this->get_size ( ), vec->buffer, 1, this->buffer, 1, &sparam );

#else
    ERROR;
#endif
}

// CPUmkl_lVector::Rotmg()
// Rotmg

template <>
void CPUmkl_lVector<double>::Rotmg ( double *sd1, double *sd2, double *x1, const double &x2, double *sparam ) const
{
#ifdef WITH_MKL

    mkl_set_num_threads ( this->num_threads_ );
    cblas_drotmg ( sd1, sd2, x1, x2, sparam );

#else
    ERROR;
#endif
}

// CPUmkl_lVector::Rotmg()
// Rotmg

template <>
void CPUmkl_lVector<float>::Rotmg ( float *sd1, float *sd2, float *x1, const float &x2, float *sparam ) const
{
#ifdef WITH_MKL

    mkl_set_num_threads ( this->num_threads_ );
    cblas_srotmg ( sd1, sd2, x1, x2, sparam );

#else
    ERROR;
#endif
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::CastFrom ( const CPU_lVector<double>& vec )
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

    for ( int i = 0; i < this->get_size ( ); ++i )
        ( this->buffer )[i] = static_cast < ValueType > ( vec.buffer[i] );
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::CastFrom ( const CPU_lVector<float>& vec )
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

    for ( int i = 0; i < this->get_size ( ); ++i )
        ( this->buffer )[i] = static_cast < ValueType > ( vec.buffer[i] );
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::CastTo ( CPU_lVector<double>& vec ) const
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

    for ( int i = 0; i < this->get_size ( ); ++i )
        vec.buffer[i] = static_cast < double > ( ( this->buffer )[i] );
}

template <typename ValueType>
void CPUmkl_lVector<ValueType>::CastTo ( CPU_lVector<float>& vec ) const
{
    assert ( this->get_size ( ) == vec.get_size ( ) );

    for ( int i = 0; i < this->get_size ( ); ++i )
        vec.buffer[i] = static_cast < float > ( ( this->buffer )[i] );
}

//
//
//
//

template class CPUmkl_lVector<double>;
template class CPUmkl_lVector<float>;
