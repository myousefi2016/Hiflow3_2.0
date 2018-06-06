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

#ifndef HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_PARASAILS_H_
#    define HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_PARASAILS_H_

#    include <mpi.h>
#    include <cstdlib>

#    include "common/log.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_solver/hypre_linear_solver.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        class HyprePreconditionerParaSails : public HyprePreconditioner<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            HyprePreconditionerParaSails ( );

            HyprePreconditionerParaSails ( MPI_Comm &comm );

            ~HyprePreconditionerParaSails ( );

            void Init ( MPI_Comm &comm );

            //          void SetupOperator ( OperatorType& op );

            // Clear allocated data
            void Clear ( );

            /// Solves a linear system.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if solver succeeded
            LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );

#    ifdef WITH_HYPRE
            // Get pointer to solve function of preconditioner

            HYPRE_PtrToSolverFcn get_solve_function ( )
            {
                return ( HYPRE_PtrToSolverFcn ) HYPRE_ParCSRParaSailsSolve;
            }

            // Get pointer to setup function of preconditioner

            HYPRE_PtrToSolverFcn get_setup_function ( )
            {
                return ( HYPRE_PtrToSolverFcn ) HYPRE_ParCSRParaSailsSetup;
            }

            // Get hypre preconditioner object

            HYPRE_Solver& get_solver ( )
            {
                return this->solver_;
            }

            /// Setup Preconditioner 

            void Build ( )
            {
                assert ( this->op_ != NULL );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Build Solver", 1 );
                }

                // Create dummy vector, only needed for interface;	
                HYPRE_ParVector tmp;
                HYPRE_ParCSRParaSailsSetup ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), tmp, tmp );

                this->SetInitialized ( true );
                this->SetModifiedParam ( false );
                this->SetModifiedOperator ( false );
                this->SetState ( true );
            }

            /// Destroy solver object 

            void DestroySolver ( )
            {
                HYPRE_ParaSailsDestroy ( this->solver_ );
                this->SetInitialized ( false );
            }

            // Set the threshold and levels parameter

            int SetParams ( DataType thresh, int nlevels )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Threshold", thresh );
                    LOG_INFO ( "Number of levels", nlevels );
                }
                return HYPRE_ParaSailsSetParams ( this->solver_, thresh, nlevels );
            }

            // Set the filter parameter

            int SetFilter ( DataType filter )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Filter", filter );
                }
                return HYPRE_ParaSailsSetFilter ( this->solver_, filter );
            }

            // Set the symmetry parameter

            int SetSym ( int sym )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Symmetry", sym );
                }
                return HYPRE_ParaSailsSetSym ( this->solver_, sym );
            }

            // Set the load balance parameter

            int SetLoadbal ( DataType loadbal )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Load balancing parameter", loadbal );
                }
                return HYPRE_ParaSailsSetLoadbal ( this->solver_, loadbal );
            }

            // Set the pattern reuse

            int SetReuse ( int reuse )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Reuse", reuse );
                }
                return HYPRE_ParaSailsSetReuse ( this->solver_, reuse );
            }

            // Set the logging parameter

            int SetLogging ( int logging )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Logging parameter", logging );
                }
                return HYPRE_ParaSailsSetLogging ( this->solver_, logging );
            }

#    endif

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_PARASAILS_H_
