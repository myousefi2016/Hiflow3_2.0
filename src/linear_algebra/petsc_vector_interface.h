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
/// @date 2015-11-27

#ifndef HIFLOW_LINEARALGEBRA_PETSC_VECTOR_INTERFACE_H_
#    define HIFLOW_LINEARALGEBRA_PETSC_VECTOR_INTERFACE_H_

#    include "petsc.h"

namespace hiflow
{
    namespace la
    {

        /// Definition of PETSc vector wrapper class
        namespace petsc
        {

            struct Vec_wrapper
            {
                ::Vec vec_;
            };

        } // namespace petsc
    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_PETSC_VECTOR_INTERFACE_H_
