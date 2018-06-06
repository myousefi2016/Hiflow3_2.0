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

#ifndef HIFLOW_LINEARSOLVER_PETSC_LINEAR_SOLVER_H_
#    define HIFLOW_LINEARSOLVER_PETSC_LINEAR_SOLVER_H_

#    include <cstdlib>
#    include "common/log.h"
#    include "linear_algebra/petsc_environment.h"
#    include "linear_algebra/petsc_la_descriptor.h"
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/petsc_preconditioner.h"

// TODO: Fix order dependency of next inclusion
#    include "common/smart_pointers.h"

namespace hiflow
{
    namespace la
    {

        /// Forwarding PETSc KSP class
        namespace petsc
        {
            class KSP_wrapper;
        }

        /// @brief Base class for wrappers to solvers in PETSc library

        template <class LAD>
        class PETScLinearSolver : public PETScPreconditioner<LAD>, public LinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            PETScLinearSolver ( petsc::KSP_wrapper* ptr_ksp_wrapper )
            : ptr_ksp_wrapper_ ( ptr_ksp_wrapper )
            {
                PETScEnvironment::initialize ( );
            }

            PETScLinearSolver ( )
            {
            }

            virtual ~PETScLinearSolver ( )
            {
            };
            /// Initialize convergence control
            /// @param[in] maxits Maximum number of iterations
            /// @param[in] abstol Absolute tolerance for residual
            /// @param[in] reltol Relative tolerance for residual
            virtual void InitControl ( int maxits, DataType abstol, DataType reltol ) = 0;

            virtual void SetupOperator ( OperatorType& op ) = 0;

            /// Solves a linear system.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if solver succeeded
            virtual LinearSolverState Solve ( const VectorType& b, VectorType* x ) = 0;

            /// Clear allocated data

            virtual void Clear ( )
            {
            }

            /// @return residual

            DataType res ( ) const
            {
                return this->res_;
            }

            /// @return Number of iterations for last solve.

            int iter ( ) const
            {
                return this->iter_;
            }

          protected:
            /// PETSc solver object Krylov subspace method (KSP)
            hiflow::scoped_ptr<petsc::KSP_wrapper> ptr_ksp_wrapper_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PETSC_LINEAR_SOLVER_H_
