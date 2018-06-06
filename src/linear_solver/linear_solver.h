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

/// @author Chandramowli Subramanian, Philipp Gerstner

#ifndef HIFLOW_LINEARSOLVER_LINEAR_SOLVER_H_
#    define HIFLOW_LINEARSOLVER_LINEAR_SOLVER_H_

#    include <cstdlib>

#    include <mpi.h>
#    include "common/iterate_control.h"
#    include "common/log.h"
#    include "linear_solver/preconditioner.h"
#    include "nonlinear/nonlinear_problem.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Base class for all linear solvers in HiFlow.
        ///
        /// A linear solver is in particular a preconditioner.

        template<class LAD>
        class LinearSolver : virtual public Preconditioner<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            LinearSolver ( )
            : res_ ( 0.0 ),
            iter_ ( 0 ),
            maxits_ ( 1000 ),
            reltol_ ( 1e-8 ),
            abstol_ ( 1e-15 ),
            divtol_ ( 1e6 ),
            filter_solution_ ( false ),
            precond_ ( NULL ),
            non_lin_op_ ( NULL )
            {
                this->SetMethod ( "NoPreconditioning" );
                this->state_ = false;
                this->print_level_ = 0;
                this->modified_op_ = false;
                this->reuse_ = true;
                this->op_ = NULL;
                this->use_solver_op_ = false;
                this->is_critical_hypre_solver_ = false;
                this->print_level_ = 0;
            };

            virtual ~LinearSolver ( )
            {
            };

            virtual void InitControl ( int maxits, double abstol, double reltol, double divtol );

            /// Set relative tolerance. Function is needed in nonlinear solver when
            /// a forcing strategy is applied

            virtual void SetRelativeTolerance ( double reltol )
            {
                this->reltol_ = reltol;
            }

            /// Sets up paramaters.

            virtual void InitParameter ( )
            {
            }

            /// Sets up paramaters.

            virtual void InitParameter ( std::string const& type )
            {
            }

            /// Sets up a preconditioner for the linear solver.
            /// @param precond preconditioner

            virtual void SetupPreconditioner ( Preconditioner<LAD>& precond )
            {
                this->precond_ = &precond;
                this->SetMethod ( "RightPreconditioning" );
            }

            /// Sets up the operator, e.g. the system matrix.
            /// If implemented, it must be invoked before
            /// @c ApplyPreconditioner is called.
            /// For instance it could compute an ILU decomposition of a GlobalMatrix.
            /// Sets up the operator for the linear solver.
            /// @param op linear operator

            virtual void SetupOperator ( OperatorType& op )
            {
                this->op_ = &op;
                this->SetModifiedOperator ( true );

                if ( this->precond_ != NULL )
                {
                    if ( this->precond_->GetUseSolverOperator ( ) )
                    {
                        this->precond_->SetupOperator ( op );
                        this->precond_->SetUseSolverOperator ( true );
                    }
                }
            }

            /// Left, Right, or no preconditining
            /// \param[in] method "NoPreconditioning" or "Preconditioning" or "LeftPreconditioning" or "RightPreconditioning"

            virtual void SetMethod ( const std::string& method )
            {
                this->precond_method_ = method;
            }

            /// @return Type of preconditioning

            virtual const std::string& Method ( ) const
            {
                return this->precond_method_;
            }

            /// Solves a linear system.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if solver succeeded
            virtual LinearSolverState Solve ( const VectorType& b, VectorType* x ) = 0;

            /// Prepare Linear solver as preconditioner

            virtual void Build ( )
            {
                assert ( this->op_ != NULL );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Build Solver", 1 );
                }

                if ( this->precond_ != NULL )
                {
                    if ( !this->precond_->GetState ( ) )
                    {
                        this->BuildPreconditioner ( );
                    }
                }
                this->SetState ( true );
                this->SetModifiedOperator ( false );
            }

            /// Build the preconditioner, such that it is ready to be used inside of a solver.
            /// The Build function of the underlying preconditioner is always called, even if precond_setup == true

            virtual void BuildPreconditioner ( )
            {
                assert ( this->op_ != NULL );
                assert ( this->precond_ != NULL );

                if ( this->precond_->GetOperator ( ) == NULL )
                {
                    this->precond_->SetupOperator ( *this->op_ );
                    this->precond_->SetUseSolverOperator ( true );
                }

                this->precond_->Build ( );
                this->precond_->SetState ( true );
                this->precond_->SetModifiedOperator ( false );
            }

            inline LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );

            /// Clear allocated data

            virtual void Clear ( )
            {
                this->res_ = 0.;
                this->iter_ = 0;
                this->op_ = NULL;
                this->precond_ = NULL;
                this->non_lin_op_ = NULL;
                this->filter_solution_ = false;
                this->print_level_ = 0;
                this->modified_op_ = false;
                this->state_ = false;
                this->reuse_ = true;
                this->use_solver_op_ = false;
                this->is_critical_hypre_solver_ = false;
                this->maxits_ = 1000;
                this->abstol_ = 1e-15;
                this->reltol_ = 1e-8;
                this->divtol_ = 1e6;
            }

            /// Set NonlinearProblem operator which provides ApplyFilter function
            /// \param[in] nlp Pointer to nonlinear Problem

            virtual void SetupNonLinProblem ( NonlinearProblem<LAD> * nlp )
            {
                this->non_lin_op_ = nlp;
                this->filter_solution_ = true;
            }

            /// @return pointer to preconditioner

            virtual Preconditioner<LAD>* GetPreconditioner ( )
            {
                return this->precond_;
            }

            /// @return iterate control

            virtual IterateControl& control ( )
            {
                return this->control_;
            }

            /// @return residual

            virtual DataType res ( ) const
            {
                return this->res_;
            }

            /// @return Number of iterations for last solve.

            virtual int iter ( ) const
            {
                return this->iter_;
            }

          protected:

            /// Pointer to preconditioner
            Preconditioner<LAD>* precond_;

            /// Pointer to nonlinear problem in order to apply solution filtering
            NonlinearProblem<LAD> *non_lin_op_;

            /// Flag whether solution should be filtered
            bool filter_solution_;

            /// Convergence control
            IterateControl control_;

            /// Residual norm
            DataType res_;

            /// Number of iterations
            int iter_;

            /// Maximum number of iterations
            int maxits_;

            /// Relative convergence tolerance
            DataType reltol_;

            /// Absolute convergence tolerance
            DataType abstol_;

            /// Absolute divergence tolerance
            DataType divtol_;

            /// Left, right or no preconditioning
            std::string precond_method_;

        };

        /// Sets up linear control.
        /// @param maxits maximum number of iteration steps
        /// @param abstol absolute tolerance of residual to converge
        /// @param reltol relative tolerance of residual to converge
        /// @param divtol relative tolerance of residual to diverge

        template<class DataType>
        void LinearSolver<DataType>::InitControl ( int maxits, double abstol, double reltol, double divtol )
        {
            this->control_.Init ( maxits, abstol, reltol, divtol );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Maximum iterations", maxits );
                LOG_INFO ( "Absolute tolerance [convergence]", abstol );
                LOG_INFO ( "Relative tolerance [convergence]", reltol );
                LOG_INFO ( "Relative tolerance [divergence]", divtol );
            }

            this->maxits_ = maxits;
            this->abstol_ = abstol;
            this->reltol_ = reltol;
            this->divtol_ = divtol;
        }

        /// Applying a linear solver as a preconditioner corresponds to solving.
        /// @param b right hand side vector
        /// @param x solution vector

        template<class DataType>
        inline LinearSolverState LinearSolver<DataType>::ApplyPreconditioner ( const VectorType& b, VectorType* x )
        {
            return this->Solve ( b, x );
        }

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_LINEAR_SOLVER_H_
