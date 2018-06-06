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

/// \author Simon Gawlok

#include "preconditioner_parallel_naive.h"

#include <cassert>

#include "common/log.h"
#include "linear_algebra/la_descriptor.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        PreconditionerParallelNaive<LAD>::PreconditionerParallelNaive ( )
        : PreconditionerBlockJacobi<LAD>( )
        {
        }

        template<class LAD>
        PreconditionerParallelNaive<LAD>::~PreconditionerParallelNaive ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void PreconditionerParallelNaive<LAD>::InitParameter ( int num_iter, DataType omega )
        {
            assert ( num_iter >= 1 );
            this->maxits_ = num_iter;
            omega_ = omega;
        }

        template<class LAD>
        void PreconditionerParallelNaive<LAD>::SetupOperator ( OperatorType& op )
        {
            this->op_ = &op;
            this->precond_diag_->SetupOperator ( op );

            this->SetModifiedOperator ( true );
            this->precond_diag_->SetModifiedOperator ( true );
        }

        template<class LAD>
        void PreconditionerParallelNaive<LAD>::SetPreconditioner ( PreconditionerBlockJacobi<LAD>& precond_diag )
        {
            this->precond_diag_ = &precond_diag;
        }

        template<class LAD>
        void PreconditionerParallelNaive<LAD>::Build ( )
        {
            assert ( this->precond_diag_ != NULL );
            assert ( this->op_ != NULL );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Build Solver", 1 );
            }

            this->precond_diag_->Build ( );

            this->SetModifiedOperator ( false );
            this->SetState ( true );
            this->precond_diag_->SetModifiedOperator ( false );
            this->precond_diag_->SetState ( true );
        }

        template<class LAD>
        LinearSolverState PreconditionerParallelNaive<LAD>::ApplyPreconditioner ( const VectorType& b,
                                                                                  VectorType* x )
        {

            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            // auxiliary vectors
            VectorType h1, h2, x_old;
            h1.CloneFromWithoutContent ( b );
            h2.CloneFromWithoutContent ( b );
            x_old.CloneFromWithoutContent ( b );

            x->Update ( );

            for ( int i = 0; i < this->maxits_; ++i )
            {
                x_old.CopyFrom ( *x );

                // Compute h1 = op.offdiag * x
                h1.Zeros ( );
                this->op_->VectorMultOffdiag ( *x, &h1 );

                // Compute h2 = b - h1
                h2.CopyFrom ( b );
                h2.Axpy ( h1, static_cast < DataType > ( -1.0 ) );

                // Apply Preconditioner of diagonal block
                this->precond_diag_->ApplyPreconditioner ( h2, x );

                x->Scale ( omega_ );
                x->Axpy ( x_old, static_cast < DataType > ( 1. - omega_ ) );
                x->Update ( );
            }

            return kSolverSuccess;

            h1.Clear ( );
            h2.Clear ( );
        }

        template<class LAD>
        void PreconditionerParallelNaive<LAD>::Clear ( )
        {
            this->precond_diag_->Clear ( );
            Preconditioner<LAD>::Clear ( );
        }

        /// template instantiation
        template class PreconditionerParallelNaive<LADescriptorCoupledD>;
        // template class PreconditionerParallelNaive<LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow
