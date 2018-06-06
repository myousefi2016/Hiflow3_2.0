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
/// @date 2015-12-21

#ifndef HIFLOW_LINEARSOLVER_PETSC_GENERAL_KSP_H_
#    define HIFLOW_LINEARSOLVER_PETSC_GENERAL_KSP_H_

#    include <mpi.h>
#    include <cstdlib>
#    include "common/log.h"
#    include "linear_solver/petsc_linear_solver.h"

namespace hiflow
{
    namespace la
    {
        namespace petsc
        {

            enum KSPType
            {
                CG,
                GMRES,
                FGMRES,
                PREONLY
            };

        } // namespace petsc

        /// @brief Wrapper class for general KSP implementations of PETSc

        template <class LAD>
        class PETScGeneralKSP : public PETScLinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            PETScGeneralKSP ( const MPI_Comm& comm, const OperatorType& op, petsc::KSPType ksp_type );

            ~PETScGeneralKSP ( );

            /// Initialize convergence control
            /// @param[in] maxits Maximum number of iterations
            /// @param[in] abstol Absolute tolerance for residual
            /// @param[in] reltol Relative tolerance for residual
            /// @param[in] divtol Divergence tolerance for residual
            void InitControl ( int maxits, double abstol, double reltol, double divtol );

            /// Initialize convergence control
            /// @param[in] maxits Maximum number of iterations
            /// @param[in] abstol Absolute tolerance for residual
            /// @param[in] reltol Relative tolerance for residual
            void InitControl ( int maxits, DataType abstol, DataType reltol );

            /// Set relative tolerance. Function is needed in nonlinear solver when
            /// a forcing strategy is applied
            void SetRelativeTolerance ( DataType reltol );

            void SetupOperator ( OperatorType& op );

            /// Solves a linear system.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if solver succeeded
            LinearSolverState Solve ( const VectorType& b, VectorType* x );

            /// Clear allocated data
            void Clear ( );

          private:
            /// Flag if solver is initialized
            bool initialized_;

            /// MPI communicator
            MPI_Comm comm_;

            /// System matrix
            const OperatorType& op_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PETSC_GENERAL_KSP_H_
