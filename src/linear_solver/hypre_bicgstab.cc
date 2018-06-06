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

#include "hypre_bicgstab.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        HypreBiCGSTAB<LAD>::HypreBiCGSTAB ( )
        : HypreLinearSolver<LAD> ( )
        {
            this->is_critical_hypre_solver_ = true;
        }

        template<class LAD>
        HypreBiCGSTAB<LAD>::~HypreBiCGSTAB ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void HypreBiCGSTAB<LAD>::Init ( )
        {
            assert ( this->op_ != NULL );

#ifdef WITH_HYPRE   
            if ( this->initialized_ )
            {
                HYPRE_ParCSRBiCGSTABDestroy ( this->solver_ );
            }
            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
            }
            MPI_Comm_dup ( this->op_->comm ( ), &this->comm_ );
            HYPRE_ParCSRBiCGSTABCreate ( this->comm_, &( this->solver_ ) );
            this->SetInitialized ( true );
#endif   
        }

        template<class LAD>
        void HypreBiCGSTAB<LAD>::Build ( const VectorType& b, VectorType* x )
        {
#ifdef WITH_HYPRE 
            assert ( this->op_ != NULL );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Build Solver", 1 );
            }

            if ( this->hypre_precond_ == NULL )
            {
                HYPRE_ParCSRBiCGSTABSetup ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), *( b.GetParVector ( ) ), *( x->GetParVector ( ) ) );
            }
            else
            {
                HYPRE_BiCGSTABSetPrecond ( this->solver_, this->hypre_precond_->get_solve_function ( ), this->hypre_precond_->get_setup_function ( ), this->hypre_precond_->get_solver ( ) );
                HYPRE_ParCSRBiCGSTABSetup ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), *( b.GetParVector ( ) ), *( x->GetParVector ( ) ) );

                // TODO is the preconditioenr really already set up at this point?
                this->hypre_precond_->SetState ( true );
                this->hypre_precond_->SetModifiedOperator ( false );
            }
            this->SetModifiedOperator ( false );
            this->SetState ( true );
#endif   
        }

        template<class LAD>
        LinearSolverState HypreBiCGSTAB<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
#ifdef WITH_HYPRE
            if ( !this->GetReuse ( ) || !this->initialized_ )
            {
                this->Init ( );
            }

            HYPRE_BiCGSTABSetMaxIter ( this->solver_, this->maxits_ ); /* max iterations */
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Maximum iterations", this->maxits_ );
            }
            HYPRE_BiCGSTABSetTol ( this->solver_, this->reltol_ ); /* relative conv. tolerance */
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Relative tolerance [convergence]", this->reltol_ );
            }
            HYPRE_BiCGSTABSetAbsoluteTol ( this->solver_, this->abstol_ ); /* absolute conv. tolerance */
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Absolute tolerance [convergence]", this->abstol_ );
            }
            //HYPRE_BiCGSTABSetPrintLevel(solver_, 2); /* print solve info */
            HYPRE_BiCGSTABSetLogging ( this->solver_, 1 ); /* needed to get run info later */

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

            HYPRE_ParCSRBiCGSTABSolve ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), *( b.GetParVector ( ) ), *( x->GetParVector ( ) ) );
            HYPRE_BiCGSTABGetNumIterations ( this->solver_, &( this->iter_ ) );
            HYPRE_BiCGSTABGetFinalRelativeResidualNorm ( this->solver_, &( this->res_ ) );
            return kSolverSuccess;
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class LAD>
        void HypreBiCGSTAB<LAD>::DestroySolver ( )
        {
#ifdef WITH_HYPRE   
            HYPRE_ParCSRBiCGSTABDestroy ( this->solver_ );

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
        void HypreBiCGSTAB<LAD>::Clear ( )
        {
#ifdef WITH_HYPRE
            if ( this->initialized_ )
            {
                HYPRE_ParCSRBiCGSTABDestroy ( this->solver_ );
            }
            HypreLinearSolver<LAD>::Clear ( );
            this->is_critical_hypre_solver_ = true;
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template class HypreBiCGSTAB<LADescriptorHypreD>;
    } // namespace la
} // namespace hiflow
