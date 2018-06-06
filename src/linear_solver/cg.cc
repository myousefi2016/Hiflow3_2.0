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

/// @author Chandramowli Subramanian

#include <cassert>
#include <cmath>

#include "cg.h"
#include "linear_algebra/la_descriptor.h"
#include "common/log.h"
#include <iomanip>

namespace hiflow
{
    namespace la
    {

        /// standard constructor

        template<class LAD>
        CG<LAD>::CG ( )
        : LinearSolver<LAD>( )
        {
            this->SetMethod ( "NoPreconditioning" );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Linear solver", "CG" );
                LOG_INFO ( "Preconditioning", this->Method ( ) );
            }
        }

        /// destructor

        template<class LAD>
        CG<LAD>::~CG ( )
        {
        }

        /// Sets parameters of the solution process
        /// @param method "NoPreconditioning" or "Preconditioning" -- whether to use preconditioning or not.

        template<class LAD>
        void CG<LAD>::InitParameter ( std::string method )
        {
            // chose method_
            this->precond_method_ = method;
            assert ( ( this->Method ( ) == "NoPreconditioning" ) ||
                     ( this->Method ( ) == "Preconditioning" ) );
        }

        /// Solves the linear system.
        /// @param [in] b right hand side vector
        /// @param [in,out] x start and solution vector

        template<class LAD>
        LinearSolverState CG<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
            if ( this->Method ( ) == "NoPreconditioning" )
                return this->SolveNoPrecond ( b, x );
            else if ( this->Method ( ) == "Preconditioning" )
                return this->SolvePrecond ( b, x );
            else if ( this->Method ( ) == "RightPreconditioning" )
                return this->SolvePrecond ( b, x );
            else if ( this->Method ( ) == "LeftPreconditioning" )
                return this->SolvePrecond ( b, x );
            else
                return kSolverError;
        }

        /// Solve without preconditioning.

        template<class LAD>
        LinearSolverState CG<LAD>::SolveNoPrecond ( const VectorType& b, VectorType* x )
        {
            assert ( this->Method ( ) == "NoPreconditioning" );
            assert ( this->op_ != NULL );

            IterateControl::State conv = IterateControl::kIterate;

            // needed vectors
            VectorType r, p, Ap;
            r.CloneFromWithoutContent ( b );
            p.CloneFromWithoutContent ( b );
            Ap.CloneFromWithoutContent ( b );

            // needed values
            DataType alpha, beta, ressquared;

            // initialization step
            this->iter_ = 0;

            this->op_->VectorMult ( *x, &r );

            r.ScaleAdd ( b, static_cast < DataType > ( -1. ) );

            ressquared = r.Dot ( r );
            this->res_ = sqrt ( ressquared );
            conv = this->control ( ).Check ( this->iter_, this->res_ );

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "CG", " without preconditioning" );
                LOG_INFO ( "CG", " starts with residual norm " << this->res_ );
            }

            p.CopyFrom ( r );
            beta = ressquared;

            // main loop
            while ( conv == IterateControl::kIterate )
            {
                ++( this->iter_ );
                this->op_->VectorMult ( p, &Ap );
                alpha = beta / ( Ap.Dot ( p ) );
                x->Axpy ( p, alpha );
                r.Axpy ( Ap, -alpha );

                ressquared = r.Dot ( r );
                this->res_ = sqrt ( ressquared );
                if ( this->print_level_ > 1 )
                {
                    LOG_INFO ( "CG residual (iteration " << this->iter_ << ")", this->res_ );
                }

                conv = this->control ( ).Check ( this->iter_, this->res_ );
                if ( conv != IterateControl::kIterate )
                {
                    break;
                }

                beta = ressquared / beta;
                p.ScaleAdd ( r, beta );
                beta = ressquared;
            }

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "CG", " without preconditioning ended after " << this->iter_ << " iterations " );
                LOG_INFO ( "CG", " with residual norm " << this->res_ );
            }

            if ( conv == IterateControl::kFailureDivergenceTol ||
                 conv == IterateControl::kFailureMaxitsExceeded )
                return kSolverExceeded;
            else
                return kSolverSuccess;
        }

        /// Solve with preconditioning.

        template<class LAD>
        LinearSolverState CG<LAD>::SolvePrecond ( const VectorType& b, VectorType* x )
        {
            assert ( this->Method ( ) != "NoPreconditioning" );
            assert ( this->op_ != NULL );
            assert ( this->precond_ != NULL );

            if ( !this->precond_->GetState ( ) || !this->precond_->GetReuse ( ) )
            {
                this->BuildPreconditioner ( );
            }

            IterateControl::State conv = IterateControl::kIterate;

            // needed vectors
            VectorType r, p, z, Ap;
            r.CloneFromWithoutContent ( b );
            p.CloneFromWithoutContent ( b );
            z.CloneFromWithoutContent ( b );
            Ap.CloneFromWithoutContent ( b );

            // needed values
            DataType alpha, beta, gamma, ressquared;

            // initialization step
            this->iter_ = 0;

            this->op_->VectorMult ( *x, &r );
            r.ScaleAdd ( b, static_cast < DataType > ( -1. ) );
            ressquared = r.Dot ( r );
            this->res_ = sqrt ( ressquared );
            conv = this->control ( ).Check ( this->iter_, this->res_ );
            z.Zeros ( );
            this->precond_->ApplyPreconditioner ( r, &z );
            p.CopyFrom ( z );

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "CG", " with preconditioning" );
                LOG_INFO ( "CG", " starts with residual norm " << this->res_ );
            }

            beta = r.Dot ( z );

            // main loop
            while ( conv == IterateControl::kIterate )
            {
                ++( this->iter_ );
                this->op_->VectorMult ( p, &Ap );
                alpha = beta / ( Ap.Dot ( p ) );
                x->Axpy ( p, alpha );
                r.Axpy ( Ap, -alpha );

                ressquared = r.Dot ( r );
                this->res_ = sqrt ( ressquared );
                if ( this->print_level_ > 1 )
                {
                    LOG_INFO ( "CG residual (iteration " << this->iter_ << ")", this->res_ );
                }

                conv = this->control ( ).Check ( this->iter_, this->res_ );
                if ( conv != IterateControl::kIterate )
                {
                    break;
                }

                z.Zeros ( );
                this->precond_->ApplyPreconditioner ( r, &z );
                gamma = r.Dot ( z );
                beta = gamma / beta;
                p.ScaleAdd ( z, beta );
                beta = gamma;
            }

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "CG", " with preconditioning ended after " << this->iter_ << " iterations " );
                LOG_INFO ( "CG", " with residual norm " << this->res_ );
            }

            if ( conv == IterateControl::kFailureDivergenceTol ||
                 conv == IterateControl::kFailureMaxitsExceeded )
                return kSolverExceeded;
            else
                return kSolverSuccess;
        }

        /// template instantiation
        template class CG<LADescriptorCoupledD>;
        template class CG<LADescriptorCoupledS>;
        template class CG<LADescriptorPolynomialChaosD>;
#ifdef WITH_HYPRE
        template class CG<LADescriptorHypreD>;
        template class CG<LADescriptorPolynomialChaosExpansionD>;
#endif

    } // namespace la
} // namespace hiflow
