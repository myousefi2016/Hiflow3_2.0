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
/// @date 2015-07-31

#ifndef SRC_LINEAR_ALGEBRA_LMP_CUDA_GPUCUBLAS2_CSR_LMATRIX_H_
#    define SRC_LINEAR_ALGEBRA_LMP_CUDA_GPUCUBLAS2_CSR_LMATRIX_H_

#    include "config.h"
#    include "lmatrix_csr_gpu.h"

#    ifdef WITH_CUDA
#        include <cuda.h>
#        include <cusparse_v2.h>
#    endif

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        class GPUcublas2_CSR_lMatrix : public GPU_CSR_lMatrix<ValueType>
        {
          public:
            GPUcublas2_CSR_lMatrix ( );

            GPUcublas2_CSR_lMatrix ( int init_nnz, int init_num_row, int init_num_col,
                                     std::string const &init_name );

            virtual ~GPUcublas2_CSR_lMatrix ( );

            virtual void VectorMultAdd ( const hiflow::la::lVector<ValueType> &invec,
                                         hiflow::la::lVector<ValueType> *outvec ) const;

            virtual void VectorMultAdd ( const GPU_lVector<ValueType> &invec,
                                         GPU_lVector<ValueType> *outvec ) const;

            virtual void VectorMult ( const hiflow::la::lVector<ValueType> &invec,
                                      hiflow::la::lVector<ValueType> *outvec ) const;

            virtual void VectorMult ( const GPU_lVector<ValueType> &invec,
                                      GPU_lVector<ValueType> *outvec ) const;

            virtual void VectorMultNoDiag ( const lVector<ValueType> &in,
                                            lVector<ValueType> *out ) const;

          protected:
#    ifdef WITH_CUDA
            cusparseHandle_t handle_;
#    endif
        };

    } // namespace la
} // namespace hiflow

#endif /* SRC_LINEAR_ALGEBRA_LMP_CUDA_GPUCUBLAS2_CSR_LMATRIX_H_ */
