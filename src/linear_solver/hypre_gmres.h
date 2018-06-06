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

#ifndef HIFLOW_LINEARSOLVER_HYPRE_GMRES_H_
#    define HIFLOW_LINEARSOLVER_HYPRE_GMRES_H_

#    include <mpi.h>
#    include <cstdlib>

#    include "common/log.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_solver/hypre_linear_solver.h"

namespace hiflow
{
    namespace la
    {
        /// @author Simon Gawlok

        /// @brief Wrapper class for GMRES implementation of Hypre
        /// A linear solver is in particular a preconditioner.

        template<class LAD>
        class HypreGMRES : public HypreLinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            HypreGMRES ( );

            ~HypreGMRES ( );

            /// Set maximum basis size of Krylov subspace

            void InitParameter ( const int size_basis, std::string method = "RightPreconditioning" )
            {
                this->basis_size_ = size_basis;
                this->precond_method_ = method;
            }

            /// Solves a linear system.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if solver succeeded
            LinearSolverState Solve ( const VectorType& b, VectorType* x );

            /// Clear allocated data
            void Clear ( );

            /// Initialize solver object 
            void Init ( );

            /// Build solver + preconditioner 
            void Build ( const VectorType& b, VectorType* x );

            /// Destroy solver object 

            void DestroySolver ( );

#    ifdef WITH_HYPRE
            /// Get pointer to solve function of preconditioner

            HYPRE_PtrToSolverFcn get_solve_function ( )
            {
                return ( HYPRE_PtrToSolverFcn ) HYPRE_ParCSRGMRESSolve;
            }

            /// Get pointer to setup function of preconditioner

            HYPRE_PtrToSolverFcn get_setup_function ( )
            {
                return ( HYPRE_PtrToSolverFcn ) HYPRE_ParCSRGMRESSetup;
            }

            /// Get hypre preconditioner object

            HYPRE_Solver& get_solver ( )
            {
                return this->solver_;
            }

#    endif

          private:
            /// Maximum size of Krylov subspace basis
            int basis_size_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_HYPRE_GMRES_H_
