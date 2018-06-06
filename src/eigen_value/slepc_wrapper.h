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

/// @author Philipp Gerstner

#ifndef HIFLOW_EIGENVALUE_SLEPC_WRAPPER_H_
#    define HIFLOW_EIGENVALUE_SLEPC_WRAPPER_H_

#    include "slepc.h"

namespace hiflow
{
    namespace la
    {

        /// Definition of PETSc matrix wrapper class
        namespace slepc
        {

            struct EPS_wrapper
            {
                ::EPS eps_;
            };

            struct ST_wrapper
            {
                ::ST st_;
            };

            struct SVD_wrapper
            {
                ::SVD svd_;
            };

        } // namespace petsc
    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_PETSC_MATRIX_INTERFACE_H_
