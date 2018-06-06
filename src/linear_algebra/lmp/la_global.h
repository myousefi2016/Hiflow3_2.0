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

#ifndef __LA_GLOBAL_H
#    define __LA_GLOBAL_H

/// @brief Global Platform, Implementation, Matrix Format IDs
/// @author Dimitar Lukarski
///
/// Global ID for Platform, Implementation, Matrix Format, Preconditioners
/// System Parameter Descriptor

#    ifdef WITH_OPENCL
#        ifdef __APPLE__
#            include <cl.h>
#        else
#            include <CL/cl.h>
#        endif
#        include "opencl/opencl_global.h"
#    endif

namespace hiflow
{
    namespace la
    {

        enum PLATFORM
        {
            CPU = 0, GPU = 1, OPENCL = 2
        };

        enum IMPLEMENTATION
        {
            NAIVE, BLAS, OPENMP, MKL, SCALAR, SCALAR_TEX, OPEN_CL, CUBLAS2
        };

        enum MATRIX_FORMAT
        {
            DENSE, CSR, COO, ELL
        };

        enum MATRIX_FREE_PRECOND
        {
            NOPRECOND, JACOBI, GAUSS_SEIDEL, SGAUSS_SEIDEL, SOR, SSOR, ILU, ILU2
        };

        /// @brief System Parameter Descriptor

        struct SYSTEM
        {
            enum PLATFORM Platform;

            bool GPU_CUBLAS;

            bool Initialized;

            int rank;

#    ifdef WITH_OPENCL
            opencl_manager *my_manager;
            bool opencl_initialized;
#    endif
            bool is_double;

        };

    } // namespace la
} // namespace hiflow

#endif
