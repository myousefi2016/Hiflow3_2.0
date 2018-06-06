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
/// @date 2015-12-02

#ifndef HIFLOW_LINEARSOLVER_PETSC_PRECONDITIONER_H_
#    define HIFLOW_LINEARSOLVER_PETSC_PRECONDITIONER_H_

#    include <iostream>
#    include <cmath>
#    include "config.h"
#    include "common/log.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Base class for all PETSc preconditioners in HiFlow.

        template <class LAD>
        class PETScPreconditioner
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            PETScPreconditioner ( )
            {
            }

            virtual ~PETScPreconditioner ( )
            {
            }

            /// Sets up the operator, e.g. the system matrix.
            /// If implemented, it must be invoked before
            /// @c ApplyPreconditioner is called.
            /// For instance it could compute an ILU decomposition of a GlobalMatrix.
            /// @param op system matrix
            virtual void SetupOperator ( OperatorType& op ) = 0;

            /// Erase possibly allocated data.

            virtual void Clear ( )
            {
            }

            /// Print possibly info

            virtual void Print ( std::ostream& out = std::cout ) const
            {
            };
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PETSC_PRECONDITIONER_H_
