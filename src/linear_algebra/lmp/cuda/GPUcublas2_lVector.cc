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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-08-03

#include "GPUcublas2_lVector.h"
#include "../lmp_log.h"
#include <assert.h>
#include <iostream>
#include <stdlib.h>

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        GPUcublas2_lVector<ValueType>::GPUcublas2_lVector ( )
        {
#ifdef WITH_CUDA
            this->implementation_name_ = "cublas2";
            this->implementation_id_ = CUBLAS2;
            cublasCreate ( &this->handle_ );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        GPUcublas2_lVector<ValueType>::GPUcublas2_lVector ( int size,
                                                            std::string const &name )
        {
#ifdef WITH_CUDA
            this->Init ( size, name );
            this->implementation_name_ = "cublas2";
            this->implementation_id_ = CUBLAS2;
            cublasCreate ( &this->handle_ );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        GPUcublas2_lVector<ValueType>::~GPUcublas2_lVector ( )
        {
#ifdef WITH_CUDA
            cublasDestroy ( this->handle_ );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        int GPUcublas2_lVector<ValueType>::ArgMin ( ) const
        {
#ifdef WITH_CUDA
            int result = 0;
            ArgMin ( &result );
            return result;
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
            return -1;
#endif
        }

        template <>
        void GPUcublas2_lVector<double>::ArgMin ( int *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            cublasIdamin ( this->handle_, this->get_size ( ), this->buffer, 1, result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::ArgMin ( int *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            cublasIsamin ( this->handle_, this->get_size ( ), this->buffer, 1, result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        int GPUcublas2_lVector<ValueType>::ArgMax ( ) const
        {
#ifdef WITH_CUDA
            int result = 0;
            ArgMax ( &result );
            return result;
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
            return -1;
#endif
        }

        template <>
        void GPUcublas2_lVector<double>::ArgMax ( int *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            cublasIdamax ( this->handle_, this->get_size ( ), this->buffer, 1, result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::ArgMax ( int *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            cublasIsamax ( this->handle_, this->get_size ( ), this->buffer, 1, result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        ValueType GPUcublas2_lVector<ValueType>::Norm1 ( ) const
        {
#ifdef WITH_CUDA
            ValueType result = 0.0;
            Norm1 ( &result );
            return result;
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
            return -1.0;
#endif
        }

        template <>
        void GPUcublas2_lVector<double>::Norm1 ( double *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            cublasDasum ( this->handle_, this->get_size ( ), this->buffer, 1, result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Norm1 ( float *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            cublasSasum ( this->handle_, this->get_size ( ), this->buffer, 1, result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        ValueType GPUcublas2_lVector<ValueType>::Norm2 ( ) const
        {
#ifdef WITH_CUDA
            ValueType result = 0.0;
            Norm2 ( &result );
            return result;
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
            return -1.0;
#endif
        }

        template <>
        void GPUcublas2_lVector<double>::Norm2 ( double *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            cublasDnrm2 ( this->handle_, this->get_size ( ), this->buffer, 1, result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Norm2 ( float *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            cublasSnrm2 ( this->handle_, this->get_size ( ), this->buffer, 1, result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        ValueType GPUcublas2_lVector<ValueType>::NormMax ( ) const
        {
#ifdef WITH_CUDA
            ValueType result = 0.0;
            NormMax ( &result );
            return result;
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
            return -1.0;
#endif
        }

        template <typename ValueType>
        void GPUcublas2_lVector<ValueType>::NormMax ( ValueType *result ) const
        {
#ifdef WITH_CUDA
            *result = std::abs ( this->buffer[this->ArgMax ( )] );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        ValueType GPUcublas2_lVector<ValueType>::Dot (
                                                       const lVector<ValueType> &vec ) const
        {
            ValueType result;
            Dot ( vec, &result );
            return result;
        }

        template <typename ValueType>
        void GPUcublas2_lVector<ValueType>::Dot ( const lVector<ValueType> &vec,
                                                  ValueType *result ) const
        {
            const GPU_lVector<ValueType> *casted_vec;

            if ( casted_vec = dynamic_cast < const GPU_lVector<ValueType> * > ( &vec ) )
            {
                *result = this->Dot ( *casted_vec );
            }
            else
            {
                LOG_ERROR ( "ERROR GPUcublas2_lVector::Dot unsupported vectors" );
                this->print ( );
                vec.print ( );
                exit ( -1 );
            }
        }

        template <typename ValueType>
        ValueType GPUcublas2_lVector<ValueType>::Dot (
                                                       const GPU_lVector<ValueType> &vec ) const
        {
#ifdef WITH_CUDA
            ValueType result;
            Dot ( vec, &result );
            return result;
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
            return -1.0;
#endif
        }

        template <>
        void GPUcublas2_lVector<double>::Dot ( const GPU_lVector<double> &vec,
                                               double *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec.get_size ( ) );
            cublasDdot ( this->handle_, this->get_size ( ), this->buffer, 1, vec.buffer, 1,
                         result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Dot ( const GPU_lVector<float> &vec,
                                              float *result ) const
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec.get_size ( ) );
            cublasSdot ( this->handle_, this->get_size ( ), this->buffer, 1, vec.buffer, 1,
                         result );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        void GPUcublas2_lVector<ValueType>::Axpy ( const lVector<ValueType> &vec,
                                                   const ValueType scalar )
        {
            const GPU_lVector<ValueType> *casted_vec;

            if ( casted_vec = dynamic_cast < const GPU_lVector<ValueType> * > ( &vec ) )
            {
                this->Axpy ( *casted_vec, scalar );
            }
            else
            {
                LOG_ERROR ( "ERROR GPUcublas2_lVector::Axpy unsupported vectors" );
                this->print ( );
                vec.print ( );
                exit ( -1 );
            }
        }

        template <>
        void GPUcublas2_lVector<double>::Axpy ( const GPU_lVector<double> &vec,
                                                const double scalar )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec.get_size ( ) );

            cublasDaxpy ( this->handle_, this->get_size ( ), &scalar, vec.buffer, 1,
                          this->buffer, 1 );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Axpy ( const GPU_lVector<float> &vec,
                                               const float scalar )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec.get_size ( ) );

            cublasSaxpy ( this->handle_, this->get_size ( ), &scalar, vec.buffer, 1,
                          this->buffer, 1 );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        void GPUcublas2_lVector<ValueType>::ScaleAdd ( const ValueType scalar,
                                                       const lVector<ValueType> &vec )
        {
            const GPU_lVector<ValueType> *casted_vec;

            if ( casted_vec = dynamic_cast < const GPU_lVector<ValueType> * > ( &vec ) )
            {
                this->ScaleAdd ( scalar, *casted_vec );
            }
            else
            {
                LOG_ERROR ( "ERROR GPUcublas2_lVector::ScaleAdd unsupported vectors" );
                this->print ( );
                vec.print ( );
                exit ( -1 );
            }
        }

        template <>
        void GPUcublas2_lVector<double>::ScaleAdd ( const double scalar,
                                                    const GPU_lVector<double> &vec )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec.get_size ( ) );

            cublasDscal ( this->handle_, this->get_size ( ), &scalar, this->buffer, 1 );

            double one = 1.0;
            cublasDaxpy ( this->handle_, this->get_size ( ), &one, vec.buffer, 1,
                          this->buffer, 1 );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::ScaleAdd ( const float scalar,
                                                   const GPU_lVector<float> &vec )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec.get_size ( ) );

            cublasSscal ( this->handle_, this->get_size ( ), &scalar, this->buffer, 1 );

            float one = 1.0;
            cublasSaxpy ( this->handle_, this->get_size ( ), &one, vec.buffer, 1,
                          this->buffer, 1 );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<double>::Scale ( const double scalar )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );

            cublasDscal ( this->handle_, this->get_size ( ), &scalar, this->buffer, 1 );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Scale ( const float scalar )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );

            cublasSscal ( this->handle_, this->get_size ( ), &scalar, this->buffer, 1 );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        void GPUcublas2_lVector<ValueType>::Rot ( lVector<ValueType> *vec,
                                                  const ValueType &sc,
                                                  const ValueType &ss )
        {
            GPU_lVector<ValueType> *casted_vec;

            if ( casted_vec = dynamic_cast < GPU_lVector<ValueType> * > ( vec ) )
            {
                this->Rot ( casted_vec, sc, ss );
            }
            else
            {
                LOG_ERROR ( "ERROR GPUcublas2_lVector::Rot unsupported vectors" );
                this->print ( );
                vec->print ( );
                exit ( -1 );
            }
        }

        template <>
        void GPUcublas2_lVector<double>::Rot ( GPU_lVector<double> *vec, const double &sc,
                                               const double &ss )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec->get_size ( ) );

            cublasDrot ( this->handle_, this->get_size ( ), vec->buffer, 1, this->buffer, 1,
                         &sc, &ss );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Rot ( GPU_lVector<float> *vec, const float &sc,
                                              const float &ss )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec->get_size ( ) );

            cublasSrot ( this->handle_, this->get_size ( ), vec->buffer, 1, this->buffer, 1,
                         &sc, &ss );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<double>::Rotg ( double *sa, double *sb, double *sc,
                                                double *ss ) const
        {
#ifdef WITH_CUDA
            cublasDrotg ( this->handle_, sa, sb, sc, ss );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Rotg ( float *sa, float *sb, float *sc,
                                               float *ss ) const
        {
#ifdef WITH_CUDA
            cublasSrotg ( this->handle_, sa, sb, sc, ss );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        void GPUcublas2_lVector<ValueType>::Rotm ( lVector<ValueType> *vec,
                                                   const ValueType &sparam )
        {
            GPU_lVector<ValueType> *casted_vec;

            if ( casted_vec = dynamic_cast < GPU_lVector<ValueType> * > ( vec ) )
            {
                this->Rotm ( casted_vec, sparam );
            }
            else
            {
                LOG_ERROR ( "ERROR GPUcublas2_lVector::Dotm unsupported vectors" );
                this->print ( );
                vec->print ( );
                exit ( -1 );
            }
        }

        template <>
        void GPUcublas2_lVector<double>::Rotm ( GPU_lVector<double> *vec,
                                                const double &sparam )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec->get_size ( ) );

            cublasDrotm ( this->handle_, this->get_size ( ), vec->buffer, 1, this->buffer, 1,
                          &sparam );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Rotm ( GPU_lVector<float> *vec,
                                               const float &sparam )
        {
#ifdef WITH_CUDA
            assert ( this->get_size ( ) > 0 );
            assert ( this->get_size ( ) == vec->get_size ( ) );

            cublasSrotm ( this->handle_, this->get_size ( ), vec->buffer, 1, this->buffer, 1,
                          &sparam );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<double>::Rotmg ( double *sd1, double *sd2, double *x1,
                                                 const double &x2, double *sparam ) const
        {
#ifdef WITH_CUDA
            cublasDrotmg ( this->handle_, sd1, sd2, x1, &x2, sparam );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_lVector<float>::Rotmg ( float *sd1, float *sd2, float *x1,
                                                const float &x2, float *sparam ) const
        {
#ifdef WITH_CUDA
            cublasSrotmg ( this->handle_, sd1, sd2, x1, &x2, sparam );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template class GPUcublas2_lVector<double>;
        template class GPUcublas2_lVector<float>;

    } // namespace la
} // namespace hiflow
