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

#include "hypre_preconditioner_parasails.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        HyprePreconditionerParaSails<LAD>::HyprePreconditionerParaSails ( )
        : HyprePreconditioner<LAD>( )
        {
        }

        template<class LAD>
        HyprePreconditionerParaSails<LAD>::HyprePreconditionerParaSails ( MPI_Comm &comm )
        : HyprePreconditioner<LAD>( )
        {
#ifdef WITH_HYPRE
            MPI_Comm_dup ( comm, &this->comm_ );
            HYPRE_ParaSailsCreate ( this->comm_, &( this->solver_ ) );
            this->initialized_ = true;
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class LAD>
        LinearSolverState HyprePreconditionerParaSails<LAD>::ApplyPreconditioner ( const VectorType& b, VectorType* x )
        {
#ifdef WITH_HYPRE
            if ( !this->GetState ( ) )
            {
                HYPRE_ParCSRParaSailsSetup ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), *( b.GetParVector ( ) ), *( x->GetParVector ( ) ) );

                this->SetModifiedOperator ( false );
                this->SetState ( true );
            }
            HYPRE_ParCSRParaSailsSolve ( this->solver_, *( this->op_->GetParCSRMatrix ( ) ), *( b.GetParVector ( ) ), *( x->GetParVector ( ) ) );
            return kSolverSuccess;
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif   
        }

        template<class LAD>
        HyprePreconditionerParaSails<LAD>::~HyprePreconditionerParaSails ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void HyprePreconditionerParaSails<LAD>::Init ( MPI_Comm &comm )
        {
            this->comm_ = MPI_COMM_NULL;
            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
            }
            this->SetInitialized ( false );
#ifdef WITH_HYPRE
            MPI_Comm_dup ( comm, &this->comm_ );
            HYPRE_ParaSailsCreate ( this->comm_, &( this->solver_ ) );
            this->SetInitialized ( true );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        /*
        template<class LAD>
        void HyprePreconditionerParaSails<LAD>::SetupOperator( OperatorType& op )
        {
            LOG_ERROR( "Operator is set by solver!" );
            exit( -1 );
        }
         */
        template<class LAD>
        void HyprePreconditionerParaSails<LAD>::Clear ( )
        {
#ifdef WITH_HYPRE
            if ( this->initialized_ )
            {
                HYPRE_ParaSailsDestroy ( this->solver_ );
            }

            HyprePreconditioner<LAD>::Clear ( );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template class HyprePreconditionerParaSails<LADescriptorHypreD>;
    } // namespace la
} // namespace hiflow
