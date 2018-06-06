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

#ifndef HIFLOW_LINEARALGEBRA_PETSC_WRAPPER_H_
#    define HIFLOW_LINEARALGEBRA_PETSC_WRAPPER_H_

#    include "petsc.h"

namespace hiflow
{
    namespace la
    {

        /// Definition of PETSc matrix wrapper class
        namespace petsc
        {

            struct Mat_wrapper
            {
                ::Mat mat_;
            };

            struct Vec_wrapper
            {
                ::Vec vec_;
            };

            struct KSP_wrapper
            {
                ::KSP ksp_;
            };

            struct PC_wrapper
            {
                ::PC pc_;
            };

            struct Val_wrapper
            {
                ::PetscScalar val_;
                ::PetscReal real_;
                ::PetscReal imag_;
            };

        } // namespace petsc
    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_PETSC_MATRIX_INTERFACE_H_
