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

#ifndef HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_EUCLID_H_
#    define HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_EUCLID_H_

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

        /// @brief Wrapper class for BoomerAMG implementation of Hypre
        /// A linear solver is in particular a preconditioner.

        template<class LAD>
        class HyprePreconditionerEuclid : public HyprePreconditioner<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            HyprePreconditionerEuclid ( MPI_Comm &comm );

            ~HyprePreconditionerEuclid ( );

            //          void SetupOperator ( OperatorType& op );

            /// Clear allocated data
            void Clear ( );

            /// Solves a linear system.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if solver succeeded
            LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );

#    ifdef WITH_HYPRE
            /// Get pointer to solve function of preconditioner

            HYPRE_PtrToSolverFcn get_solve_function ( )
            {
                return ( HYPRE_PtrToSolverFcn ) HYPRE_EuclidSolve;
            }

            /// Get pointer to setup function of preconditioner

            HYPRE_PtrToSolverFcn get_setup_function ( )
            {
                return ( HYPRE_PtrToSolverFcn ) HYPRE_EuclidSetup;
            }

            /// Get pointer to solver

            HYPRE_Solver& get_solver ( )
            {
                return this->solver_;
            }

            /// Setup Preconditioner 

            void Build ( )
            {
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Build Solver", 1 );
                }
                assert ( this->op_ != NULL );
                // Create dummy vector, only needed for interface;	
                HYPRE_ParVector tmp;
                HYPRE_EuclidSetup ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), tmp, tmp );

                this->SetInitialized ( true );
                this->SetModifiedParam ( false );
                this->SetModifiedOperator ( false );
                this->SetState ( true );
            }

            /// Destroy solver object 

            void DestroySolver ( )
            {
                HYPRE_EuclidDestroy ( this->solver_ );
                this->SetInitialized ( false );
            }

            /// Set level k for ILU(k) factorization, default: 1

            int SetLevel ( int level = 1 )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Level for ILU(k) factorization", level );
                }
                return HYPRE_EuclidSetLevel ( this->solver_, level );
            }

            /// Use block Jacobi ILU instead of PILU

            int SetBJ ( int bj = 0 )
            {
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Use block Jacobi", bj );
                }
                this->SetModifiedParam ( true );
                return HYPRE_EuclidSetBJ ( this->solver_, bj );
            }

            /// Define a drop tolerance for ILU(k)

            int SetSparseA ( DataType sparse_A )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Drop tolerance for ILU(k)", sparse_A );
                }
                return HYPRE_EuclidSetSparseA ( this->solver_, sparse_A );
            }

            /// If row scale not equal 0, values are scaled prior to factorization so that 
            /// largest value in any row is +1 or -1.

            int SetRowScale ( int row_scale )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Row scaling", row_scale );
                }
                return HYPRE_EuclidSetRowScale ( this->solver_, row_scale );
            }

            /// uses ILUT and defines a drop tolerance relative to the largest absolute value
            /// of any entry in the row being factored. 

            int SetILUT ( DataType drop_tol )
            {
                this->SetModifiedParam ( true );
                if ( this->print_level_ > 2 )
                {
                    LOG_INFO ( "Relative drop tolerance for ILUT", drop_tol );
                }
                return HYPRE_EuclidSetILUT ( this->solver_, drop_tol );
            }

#    endif

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_EUCLID_H_
