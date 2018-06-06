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

/// \author Michael Schick

#include "polynomial_chaos/pc_meanbased_preconditioner.h"
#include <cassert>
#include "common/log.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class LAD>
        MeanbasedPreconditioner<LAD>::MeanbasedPreconditioner ( ) : la::Preconditioner<LAD>( )
        {
        }

        template<class LAD>
        MeanbasedPreconditioner<LAD>::~MeanbasedPreconditioner ( )
        {
        }

        template<class LAD>
        void MeanbasedPreconditioner<LAD>::Build ( )
        {
            assert ( this->op_ != NULL );

            // set the matrix to be used as the operator
            matrix_.CloneFrom ( *this->op_->GetModes ( )->at ( 0 ) );

            linear_solver_.SetupOperator ( matrix_ );
            linear_solver_.Build ( );

            this->SetModifiedOperator ( false );
            this->SetState ( true );
        }

        template<class LAD>
        la::LinearSolverState MeanbasedPreconditioner<LAD>::ApplyPreconditioner ( const VectorType& b,
                                                                                  VectorType* x )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            la::LinearSolverState state;
            for ( int mode = 0; mode < b.NModes ( ); ++mode )
            {
                state = this->linear_solver_.Solve ( *b.Mode ( mode ), x->Mode ( mode ) );
            }
            return state;
        }

        /// template instantiation
        template class MeanbasedPreconditioner<la::LADescriptorPolynomialChaosD>;

    } // namespace stochastic
} // namespace hiflow
