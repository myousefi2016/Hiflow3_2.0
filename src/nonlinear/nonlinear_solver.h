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

#ifndef HIFLOW_NONLINEAR_NONLINEAR_SOLVER_H_
#    define HIFLOW_NONLINEAR_NONLINEAR_SOLVER_H_

#    include <cstdlib>
#    include "linear_solver/linear_solver.h"
#    include "common/iterate_control.h"
#    include "nonlinear/nonlinear_problem.h"
#    include "common/log.h"

namespace hiflow
{

    /// Enumerator @em NonlinearSolverState as return value for the nonlinear solvers.

    enum NonlinearSolverState
    {
        kNonlinearSolverSuccess = 0,
        kNonlinearSolverExceeded,
        kNonlinearSolverInitError,
        kNonlinearSolverError
    };

    /// @brief Base class for all nonlinear solvers in HiFlow.
    /// @author Tobias Hahn
    ///
    /// Solves for x in F(x)=y with nonlinear F

    template<class LAD>
    class NonlinearSolver
    {
      public:
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;
        typedef typename LAD::DataType DataType;

        /// Enumerator @em NonlinearSolverParameter as place-holder for custom parameters.

        enum NonlinearSolverParameter
        {
        };

        NonlinearSolver ( );
        virtual ~NonlinearSolver ( );

        void InitControl ( int maxits, double abstol, double reltol, double divtol );

        /// Sets up the problem for the nonlinear solver to solve.
        /// @param op Nonlinear problem

        void SetOperator ( NonlinearProblem<LAD>& op )
        {
            this->op_ = &op;
        }

        /// Returns the problem for the nonlinear solver to solve.

        NonlinearProblem<LAD>* GetOperator ( )
        {
            return this->op_;
        }

        /// Sets up the linear solver.
        /// @param linsolve LinearSolver

        void SetLinearSolver ( la::LinearSolver<LAD>& linsolve )
        {
            this->linsolve_ = &linsolve;
        }

        /// Sets up paramaters.

        virtual NonlinearSolverState InitParameter ( NonlinearSolverParameter param )
        {
            return kNonlinearSolverSuccess;
        }

        virtual void InitParameter ( VectorType* residual, MatrixType* matrix )
        {
        }

        /// Solves a nonlinear system.
        /// @param b right hand side vector
        /// @param x solution vector
        /// @return status if solver succeeded
        virtual NonlinearSolverState Solve ( const VectorType& y, VectorType* x ) = 0;

        /// Returns the residual

        DataType GetResidual ( ) const
        {
            return this->residual_;
        }

        /// @return linear solver

        const la::LinearSolver<LAD>& linsolve ( )
        {
            return *( this->linsolve_ );
        }
        /// @return iterate control

        const IterateControl& control ( ) const
        {
            return this->control_;
        }

        IterateControl& control ( )
        {
            return this->control_;
        }

        /// Returns current solver step count
        /// @return number of current iteration step

        int iter ( )
        {
            return iter_;
        }

        /// Erase possibly allocated data.

        virtual void Clear ( )
        {
        }

      protected:
        /// Underlying linear solver
        la::LinearSolver<LAD>* linsolve_;
        /// Nonlinear problem to be solved
        NonlinearProblem<LAD>* op_;
        /// Iteration control
        IterateControl control_;
        /// Residual value
        DataType residual_;
        /// Iteration counter
        int iter_;
    };

    /// standard constructor

    template<class DataType>
    NonlinearSolver<DataType>::NonlinearSolver ( )
    {
        this->linsolve_ = NULL;
        this->op_ = NULL;
    }

    /// destructor

    template<class DataType>
    NonlinearSolver<DataType>::~NonlinearSolver ( )
    {
    }

    /// Sets up linear control.
    /// @param maxits maximum number of iteration steps
    /// @param abstol absolute tolerance of residual to converge
    /// @param reltol relative tolerance of residual to converge
    /// @param divtol relative tolerance of residual to diverge

    template<class DataType>
    void NonlinearSolver<DataType>::InitControl ( int maxits, double abstol,
                                                  double reltol, double divtol )
    {
        this->control_.Init ( maxits, abstol, reltol, divtol );
        LOG_INFO ( "Maximum iterations", maxits );
        LOG_INFO ( "Absolute tolerance [convergence]", abstol );
        LOG_INFO ( "Relative tolerance [convergence]", reltol );
        LOG_INFO ( "Relative tolerance [divergence]", divtol );
    }

} // namespace hiflow

#endif  // HIFLOW_NONLINEARSOLVER_NONLINEAR_SOLVER_H_
