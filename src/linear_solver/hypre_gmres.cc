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

/// @author Simon Gawlok

#include "hypre_gmres.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        HypreGMRES<LAD>::HypreGMRES ( )
        : HypreLinearSolver<LAD>( ),
        basis_size_ ( 100 )
        {
            this->is_critical_hypre_solver_ = true;
        }

        template<class LAD>
        HypreGMRES<LAD>::~HypreGMRES ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void HypreGMRES<LAD>::Init ( )
        {
            assert ( this->op_ != NULL );

#ifdef WITH_HYPRE   
            if ( this->initialized_ )
            {
                HYPRE_ParCSRGMRESDestroy ( this->solver_ );
            }
            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
            }

            MPI_Comm_dup ( this->op_->comm ( ), &this->comm_ );

            HYPRE_ParCSRGMRESCreate ( this->comm_, &( this->solver_ ) );
            this->SetInitialized ( true );
            this->SetState ( false );
#endif   
        }

        template<class LAD>
        void HypreGMRES<LAD>::Build ( const VectorType& b, VectorType* x )
        {
#ifdef WITH_HYPRE 
            assert ( this->op_ != NULL );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Build Solver", 1 );
            }
            if ( this->hypre_precond_ == NULL )
            {
                HYPRE_ParCSRGMRESSetup ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), *( b.GetParVector ( ) ), *( x->GetParVector ( ) ) );
            }
            else
            {
                HYPRE_GMRESSetPrecond ( this->solver_, this->hypre_precond_->get_solve_function ( ), this->hypre_precond_->get_setup_function ( ), this->hypre_precond_->get_solver ( ) );
                HYPRE_ParCSRGMRESSetup ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), *( b.GetParVector ( ) ), *( x->GetParVector ( ) ) );

                // TODO is the preconditioenr really already set up at this point?
                this->hypre_precond_->SetState ( true );
                this->hypre_precond_->SetModifiedOperator ( false );
            }
            this->SetModifiedOperator ( false );
            this->SetState ( true );
#endif   
        }

        template<class LAD>
        LinearSolverState HypreGMRES<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
#ifdef WITH_HYPRE
            if ( !this->GetReuse ( ) || !this->initialized_ )
            {
                this->Init ( );
            }

            HYPRE_ParCSRGMRESSetMaxIter ( this->solver_, this->maxits_ ); /* max iterations */
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Maximum iterations", this->maxits_ );
            }
            HYPRE_ParCSRGMRESSetTol ( this->solver_, this->reltol_ ); /* relative conv. tolerance */
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Relative tolerance [convergence]", this->reltol_ );
            }
            HYPRE_ParCSRGMRESSetAbsoluteTol ( this->solver_, this->abstol_ ); /* absolute conv. tolerance */
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Absolute tolerance [convergence]", this->abstol_ );
            }
            HYPRE_ParCSRGMRESSetKDim ( this->solver_, this->basis_size_ );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Maximum size of Krylov basis", this->basis_size_ );
            }

            //HYPRE_GMRESSetPrintLevel(this->solver_, 2); /* print solve info */
            HYPRE_ParCSRGMRESSetLogging ( this->solver_, 1 ); /* needed to get run info later */

            if ( this->hypre_precond_ != NULL )
            {
                if ( !this->GetReuse ( ) || !this->hypre_precond_->GetReuse ( ) || !this->GetState ( ) || !this->hypre_precond_->GetState ( ) )
                {
                    this->Build ( b, x );
                }
            }
            else
            {
                if ( !this->GetReuse ( ) || !this->GetState ( ) )
                {
                    this->Build ( b, x );
                }
            }
            HYPRE_ParCSRGMRESSolve ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), *( b.GetParVector ( ) ), *( x->GetParVector ( ) ) );
            HYPRE_ParCSRGMRESGetNumIterations ( this->solver_, &( this->iter_ ) );
            HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm ( this->solver_, &( this->res_ ) );

            return kSolverSuccess;
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class LAD>
        void HypreGMRES<LAD>::DestroySolver ( )
        {
#ifdef WITH_HYPRE   
            HYPRE_ParCSRGMRESDestroy ( this->solver_ );

            this->SetInitialized ( false );
            if ( this->hypre_precond_ != NULL )
            {
                if ( this->hypre_precond_->IsCritical ( ) )
                {
                    this->hypre_precond_->DestroySolver ( );
                }
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class LAD>
        void HypreGMRES<LAD>::Clear ( )
        {
#ifdef WITH_HYPRE
            if ( this->initialized_ )
            {
                HYPRE_ParCSRGMRESDestroy ( this->solver_ );
            }
            HypreLinearSolver<LAD>::Clear ( );
            this->is_critical_hypre_solver_ = true;
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template class HypreGMRES<LADescriptorHypreD>;
    } // namespace la
} // namespace hiflow
