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

/// @author Dimitar Lukarski, Chandramowli Subramanian

#include "preconditioner_bjacobi.h"
#include "preconditioner_bjacobi_standard.h"
#include "linear_algebra/lmp/lpreconditioner.h"

#include <cassert>

#include "common/log.h"
#include "linear_algebra/la_descriptor.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        PreconditionerBlockJacobiStand<LAD>::PreconditionerBlockJacobiStand ( )
        : PreconditionerBlockJacobi<LAD>( )
        {
            this->localPrecond_ = NULL;
        }

        template<class LAD>
        PreconditionerBlockJacobiStand<LAD>::~PreconditionerBlockJacobiStand ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Init_Jacobi ( const VectorType& vec )
        {
            this->Clear ( );
            lPreconditioner_Jacobi<typename LAD::DataType> *lp = new lPreconditioner_Jacobi<typename LAD::DataType>;
            lp->Init ( vec.interior ( ) );
            this->localPrecond_ = lp;
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Init_GaussSeidel ( )
        {
            this->Clear ( );
            lPreconditioner_GaussSeidel<typename LAD::DataType> *lp = new lPreconditioner_GaussSeidel<typename LAD::DataType>;
            lp->Init ( );
            this->localPrecond_ = lp;
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Init_SymmetricGaussSeidel ( )
        {
            this->Clear ( );
            lPreconditioner_SymmetricGaussSeidel<typename LAD::DataType> *lp = new lPreconditioner_SymmetricGaussSeidel<typename LAD::DataType>;
            lp->Init ( );
            this->localPrecond_ = lp;
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Init_SOR ( const typename LAD::DataType omega )
        {
            this->Clear ( );
            lPreconditioner_SOR<typename LAD::DataType> *lp = new lPreconditioner_SOR<typename LAD::DataType>;
            lp->Init ( omega );
            this->localPrecond_ = lp;
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Init_SSOR ( const typename LAD::DataType omega )
        {
            this->Clear ( );
            lPreconditioner_SSOR<typename LAD::DataType> *lp = new lPreconditioner_SSOR<typename LAD::DataType>;
            lp->Init ( omega );
            this->localPrecond_ = lp;
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Init_ILUp ( const int ilu_p )
        {
            this->Clear ( );
            lPreconditioner_ILUp<typename LAD::DataType> *lp = new lPreconditioner_ILUp<typename LAD::DataType>;
            lp->Init ( ilu_p );
            this->localPrecond_ = lp;
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::SetupOperator ( OperatorType& op )
        {
            this->op_ = &op;
            this->localPrecond_->SetupOperator ( op.diagonal ( ) );

            this->SetModifiedOperator ( true );
            this->localPrecond_->SetModifiedOperator ( true );
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Build ( )
        {
            assert ( this->localPrecond_ != NULL );
            assert ( this->op_ != NULL );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Build Solver", 1 );
            }
            this->localPrecond_->Build ( );
            this->localPrecond_->SetModifiedOperator ( false );
            this->localPrecond_->SetState ( true );

            this->SetModifiedOperator ( false );
            this->SetState ( true );
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Clear ( )
        {
            if ( this->localPrecond_ != NULL )
            {
                this->localPrecond_->Clear ( );
                delete this->localPrecond_;
            }
            this->localPrecond_ = NULL;
            Preconditioner<LAD>::Clear ( );
        }

        template<class LAD>
        LinearSolverState PreconditionerBlockJacobiStand<LAD>::ApplyPreconditioner ( const VectorType& b, VectorType* x )
        {
            assert ( this->op_ != NULL );
            assert ( this->localPrecond_ != NULL );

            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }
            this->localPrecond_->ApplylPreconditioner ( b.interior ( ), &( x->interior ( ) ) );

            return kSolverSuccess;
        }

        template<class LAD>
        void PreconditionerBlockJacobiStand<LAD>::Print ( std::ostream &out ) const
        {
            this->localPrecond_->print ( out );
        }

        /// template instantiation
        template class PreconditionerBlockJacobiStand<LADescriptorCoupledD>;
        // template class PreconditionerBlockJacobiStand<LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow
