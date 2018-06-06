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

#ifndef SRC_LINEAR_ALGEBRA_LMP_CUDA_GPUCUBLAS2_LVECTOR_H_
#    define SRC_LINEAR_ALGEBRA_LMP_CUDA_GPUCUBLAS2_LVECTOR_H_

#    include "config.h"
#    include "lvector_gpu.h"

#    ifdef WITH_CUDA
#        include <cuda.h>
#        include <cublas_v2.h>
#        include "cuda_utils.h"
#    endif

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        class GPUcublas2_lVector : public GPU_lVector<ValueType>
        {
          public:
            GPUcublas2_lVector ( );
            GPUcublas2_lVector ( int size, std::string const &name );
            virtual ~GPUcublas2_lVector ( );

            virtual int ArgMin ( ) const;
            virtual void ArgMin ( int *result ) const;

            virtual int ArgMax ( ) const;
            virtual void ArgMax ( int *result ) const;

            virtual ValueType Norm1 ( ) const;
            virtual void Norm1 ( ValueType *result ) const;

            virtual ValueType Norm2 ( ) const;
            virtual void Norm2 ( ValueType *result ) const;

            virtual ValueType NormMax ( ) const;
            virtual void NormMax ( ValueType *result ) const;

            virtual ValueType Dot ( const hiflow::la::lVector<ValueType> &vec ) const;
            virtual void Dot ( const hiflow::la::lVector<ValueType> &vec,
                               ValueType *result ) const;

            virtual ValueType Dot ( const GPU_lVector<ValueType> &vec ) const;
            virtual void Dot ( const GPU_lVector<ValueType> &vec, ValueType *result ) const;

            virtual void Axpy ( const hiflow::la::lVector<ValueType> &vec,
                                const ValueType scalar );
            virtual void Axpy ( const GPU_lVector<ValueType> &vec, const ValueType scalar );
            virtual void ScaleAdd ( const ValueType scalar,
                                    const hiflow::la::lVector<ValueType> &vec );
            virtual void ScaleAdd ( const ValueType scalar,
                                    const GPU_lVector<ValueType> &vec );
            virtual void Scale ( const ValueType scalar );
            virtual void Rot ( hiflow::la::lVector<ValueType> *vec, const ValueType &sc,
                               const ValueType &ss );
            virtual void Rot ( GPU_lVector<ValueType> *vec, const ValueType &sc,
                               const ValueType &ss );
            virtual void Rotg ( ValueType *sa, ValueType *sb, ValueType *sc,
                                ValueType *ss ) const;
            virtual void Rotm ( hiflow::la::lVector<ValueType> *vec,
                                const ValueType &sparam );
            virtual void Rotm ( GPU_lVector<ValueType> *vec, const ValueType &sparam );
            virtual void Rotmg ( ValueType *sd1, ValueType *sd2, ValueType *x1,
                                 const ValueType &x2, ValueType *sparam ) const;

          protected:
#    ifdef WITH_CUDA
            cublasHandle_t handle_;
#    endif
        };

    } // namespace la
} // namespace hiflow

#endif /* SRC_LINEAR_ALGEBRA_LMP_CUDA_GPUCUBLAS2_LVECTOR_H_ */
