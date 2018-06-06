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

#ifndef HIFLOW_LINEARSOLVER_PRECONDITIONER_H_
#    define HIFLOW_LINEARSOLVER_PRECONDITIONER_H_

#    include <iostream>
#    include <cassert>
#    include "common/log.h"

namespace hiflow
{
    namespace la
    {

        /// Enumerator @em LinearSolverState as return value for the preconditioners
        /// and linear solvers.

        enum LinearSolverState
        {
            kSolverSuccess = 0,
            kSolverExceeded,
            kSolverError
        };

        /// @brief Base class for all preconditioners and linear solvers in HiFlow.

        template<class LAD>
        class Preconditioner
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// Constructor

            Preconditioner ( )
            : state_ ( false ), reuse_ ( true ), modified_op_ ( false ), use_solver_op_ ( false ), print_level_ ( 0 ), is_critical_hypre_solver_ ( false )
            {
            }

            /// Destructor

            virtual ~Preconditioner ( )
            {
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
                this->SetUseSolverOperator ( false );
            }

            virtual void SetupVector ( const VectorType& vec )
            {
            };

            /// Sets up paramaters.

            virtual void InitParameter ( )
            {
            }

            /// Build the preconditioner such that is ready to be used inside of the solver

            virtual void Build ( void )
            {
                assert ( this->op_ != NULL );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Build Solver", 1 );
                }

                this->SetState ( true );
                this->SetModifiedOperator ( false );
            }

            /// Applies the preconditioner which is possibly an inexact linear solver.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if preconditioning succeeded
            virtual LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType * x ) = 0;

            /// Erase possibly allocated data.

            virtual void Clear ( )
            {
                this->op_ = NULL;
                this->modified_op_ = false;
                this->state_ = false;
                this->reuse_ = true;
                this->print_level_ = 0;
                this->use_solver_op_ = false;
                this->is_critical_hypre_solver_ = false;
            }

            /// Destroy solver object, if external libraries like Hypre are used

            virtual void DestroySolver ( )
            {
            };

            /// Set level of information in log file
            /// 0: no output, 1: initial+final residual + iteration count,
            /// 2: CPU time, 3: additional information and parameters (if available)
            /// \param[in] print_level Level

            virtual void SetPrintLevel ( int print_level )
            {
                this->print_level_ = print_level;
            }

            /// Print possibly info

            virtual void Print ( std::ostream &out = std::cout ) const
            {
            };

            /// Set state of the preconditioner, i.e. whether it is ready to use or not.
            /// In case of reuse_ == false, this function always sets the state to false
            /// @param[in] bool state

            virtual void SetState ( bool state )
            {
                this->state_ = state;
            }

            /// Get State of the preconditioner
            /// @return state

            virtual bool GetState ( )
            {
                return this->state_;
            }

            /// Set flag whether preconditioner should be resued by further calls to the outer solving routine.
            /// Usually, this option should always be set to true. However, there might be MPI communicator problems when reusing
            /// too many BoomerAMG preconditioners at the same time .
            /// @param[in] bool flag

            virtual void SetReuse ( bool flag )
            {
                this->reuse_ = flag;
            }

            /// Get reuse flag
            /// @return  flag

            virtual bool GetReuse ( )
            {
                return this->reuse_;
            }

            /// Set flag whether solver operator should be used
            /// @param[in] bool flag

            virtual void SetUseSolverOperator ( bool flag )
            {
                this->use_solver_op_ = flag;
            }

            /// Return flag if solver operator should be used
            /// @return flag

            virtual bool GetUseSolverOperator ( )
            {
                return this->use_solver_op_;
            }
            /// Set status of operator
            /// @param[in] bool flag

            virtual void SetModifiedOperator ( bool flag )
            {
                this->modified_op_ = flag;
                if ( flag )
                    this->SetState ( false );
            }

            /// Get status of operator
            /// @param[in] bool flag

            virtual bool GetModifiedOperator ( )
            {
                return this->modified_op_;
            }

            /// Return pointer to operator
            /// @return pointer

            virtual OperatorType * GetOperator ( )
            {
                return this->op_;
            }

            /// Set critical flag
            /// @param[in] flag

            virtual void SetCritical ( bool flag )
            {
                this->is_critical_hypre_solver_ = flag;
            }

            /// Return critical flag
            /// @return flag

            virtual bool IsCritical ( )
            {
                return this->is_critical_hypre_solver_;
            }

          protected:
            /// Pointer to operator
            OperatorType* op_;

            /// Flag if operator has changed
            bool modified_op_;

            /// Flag if preconditioner is set up, i.e. ready to use inside of the solver. This flag is set to false, if either the operator or some parameters have changed
            bool state_;

            /// Flag if preconditioner should be reused by several calls to outer solve routine. Only valid, if no changes for the oprator have been made. Default is true
            bool reuse_;

            /// Flag if the solver's operator is used
            bool use_solver_op_;

            /// Print Level
            int print_level_;

            /// Indicating whether preconditioner gives rise to MPI Problems, if it is reuse. E.g. Hypre CG
            bool is_critical_hypre_solver_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PRECONDITIONER_H_
