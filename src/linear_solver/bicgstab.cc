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

#include <cassert>
#include <cmath>

#include "bicgstab.h"
#include "linear_algebra/la_descriptor.h"
#include "common/log.h"

namespace hiflow
{
    namespace la
    {

        /// standard constructor

        template<class LAD>
        BiCGSTAB<LAD>::BiCGSTAB ( )
        : LinearSolver<LAD>( )
        {
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Linear solver", "BiCGSTAB" );
            }
        }

        /// destructor

        template<class LAD>
        BiCGSTAB<LAD>::~BiCGSTAB ( )
        {
            this->op_ = NULL;
            this->precond_ = NULL;
        }

        /// Sets parameters of the solution process
        /// @param method "NoPreconditioning" or "Preconditioning" -- whether to use preconditioning or not.

        template<class LAD>
        void BiCGSTAB<LAD>::InitParameter ( std::string method )
        {
            // chose method_
            this->precond_method_ = method;
            assert ( ( this->Method ( ) == "NoPreconditioning" ) ||
                     ( this->Method ( ) == "RightPreconditioning" ) ||
                     ( this->Method ( ) == "LeftPreconditioning" ) );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Preconditioning", this->Method ( ) );
            }
        }

        /// Solves the linear system.
        /// @param [in] b right hand side vector
        /// @param [in,out] x start and solution vector

        template<class LAD>
        LinearSolverState BiCGSTAB<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
            if ( this->Method ( ) == "NoPreconditioning" )
                return this->SolveNoPrecond ( b, x );
            else if ( this->Method ( ) == "RightPreconditioning" )
                return this->SolvePrecondRight ( b, x );
            else if ( this->Method ( ) == "LeftPreconditioning" )
                return this->SolvePrecondLeft ( b, x );
            else
                return kSolverError;
        }

        /// Solve without preconditioning.

        template<class LAD>
        LinearSolverState BiCGSTAB<LAD>::SolveNoPrecond ( const VectorType& b, VectorType* x )
        {
            assert ( this->Method ( ) == "NoPreconditioning" );
            assert ( this->op_ != NULL );

            IterateControl::State conv = IterateControl::kIterate;

            // needed vectors
            VectorType r, r0, p, Ap, s, t;
            r.CloneFromWithoutContent ( b );
            r0.CloneFromWithoutContent ( b );
            p.CloneFromWithoutContent ( b );
            Ap.CloneFromWithoutContent ( b );
            s.CloneFromWithoutContent ( b );
            t.CloneFromWithoutContent ( b );

            // needed values
            DataType alpha, beta, ressquared, rho, rho_old, omega, sigma, tol1;

            // initialization step
            int iter = 0;
            bool restart = true;

            // restart option for stability purposes
            while ( restart )
            {
                restart = false;

                // clear data structures
                r.Zeros ( );
                r0.Zeros ( );
                p.Zeros ( );
                Ap.Zeros ( );
                s.Zeros ( );
                t.Zeros ( );

                // compute residual
                // r0 = b - Ax
                this->op_->VectorMult ( *x, &r0 );
                r0.ScaleAdd ( b, static_cast < DataType > ( -1. ) );

                // r = r0
                r.CloneFrom ( r0 );
                ressquared = r.Dot ( r );

                // rho = rho_old = (r0, r0)
                rho = rho_old = ressquared;
                this->res_ = sqrt ( ressquared );

                // check for trivial convergence
                conv = this->control ( ).Check ( iter, this->res ( ) );

                if ( this->print_level_ > 0 )
                {
                    LOG_INFO ( "BiCGSTAB", " without preconditioning" );
                    LOG_INFO ( "BiCGSTAB", " (re)starts with residual norm " << this->res_ );
                }

                // p0 = r0
                p.CopyFrom ( r );

                // main loop
                while ( conv == IterateControl::kIterate && !restart )
                {
                    // increment interation count
                    ++iter;

                    // compute Ap = A*p
                    this->op_->VectorMult ( p, &Ap );
                    // compute sigma = (Ap, r0)
                    sigma = Ap.Dot ( r0 );

                    // compute stability tolerance 1: tol1 = 1e-15*norm2(Ap)*norm2(r0)
                    tol1 = static_cast < DataType > ( 1.e-15 ) * Ap.Norm2 ( ) * r0.Norm2 ( );

                    //check for stability
                    if ( std::abs ( sigma ) > tol1 )
                    { // stable case
                        // alpha = rho / sigma
                        alpha = rho / sigma;

                        // s = r - alpha*Ap
                        s.CloneFrom ( r );
                        s.Axpy ( Ap, -alpha );

                        // Stability criterion 2: norm2(s) > 1.e-15
                        if ( s.Norm2 ( ) > static_cast < DataType > ( 1.e-15 ) )
                        { // stable case
                            // t = A*s
                            this->op_->VectorMult ( s, &t );

                            // omega = (t, s) / (t, t)
                            omega = t.Dot ( s ) / t.Dot ( t );

                            // x = x + alpha*p + omega*s
                            x->Axpy ( p, alpha );
                            x->Axpy ( s, omega );

                            // r = s - omega*t
                            r.CloneFrom ( s );
                            r.Axpy ( t, -omega );

                            //check for convergence
                            ressquared = r.Dot ( r );
                            this->res_ = sqrt ( ressquared );
                            if ( this->print_level_ > 1 )
                            {
                                LOG_INFO ( "BiCGSTAB residual (iteration " << iter << ")", this->res_ );
                            }

                            conv = this->control ( ).Check ( iter, this->res ( ) );
                            if ( conv != IterateControl::kIterate )
                            {
                                break;
                            }

                            // rho = (r, r0);
                            rho = r.Dot ( r0 );

                            // beta = (alpha * rho) / (omega * rho_old)
                            beta = ( alpha * rho ) / ( omega * rho_old );

                            // rho_old = rho
                            rho_old = rho;

                            //p = r + beta*(p - omega*Ap)
                            p.ScaleAdd ( r, beta );
                            p.Axpy ( Ap, -omega * beta );
                        }
                        else
                        { // do stable update
                            // x = x + alpha*p
                            x->Axpy ( p, alpha );

                            // r = s
                            r.CloneFrom ( s );

                            //check for convergence
                            ressquared = r.Dot ( r );
                            this->res_ = sqrt ( ressquared );
                            if ( this->print_level_ > 1 )
                            {
                                LOG_INFO ( "BiCGSTAB residual (iteration " << iter << ")", this->res_ );
                            }

                            conv = this->control ( ).Check ( iter, this->res ( ) );
                            if ( conv != IterateControl::kIterate )
                            {
                                break;
                            }
                        }
                    }
                    else
                    { // restart
                        restart = true;
                    }
                }
            }

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "BiCGSTAB", " without preconditioning ended after " << iter << " iterations " );
                LOG_INFO ( "BiCGSTAB", " with residual norm " << this->res_ );
            }

            if ( conv == IterateControl::kFailureDivergenceTol ||
                 conv == IterateControl::kFailureMaxitsExceeded )
                return kSolverExceeded;
            else
                return kSolverSuccess;
        }

        /// Solve with preconditioning.

        template<class LAD>
        LinearSolverState BiCGSTAB<LAD>::SolvePrecondRight ( const VectorType& b, VectorType* x )
        {
            assert ( this->Method ( ) == "RightPreconditioning" );
            assert ( this->op_ != NULL );
            assert ( this->precond_ != NULL );

            if ( !this->precond_->GetState ( ) || !this->precond_->GetReuse ( ) )
            {
                this->BuildPreconditioner ( );
            }

            IterateControl::State conv = IterateControl::kIterate;

            // needed vectors
            VectorType r, r0, p, Ap, s, t, h, h1;
            r.CloneFromWithoutContent ( b );
            r0.CloneFromWithoutContent ( b );
            p.CloneFromWithoutContent ( b );
            Ap.CloneFromWithoutContent ( b );
            s.CloneFromWithoutContent ( b );
            t.CloneFromWithoutContent ( b );
            h.CloneFromWithoutContent ( b );
            h1.CloneFromWithoutContent ( b );

            // needed values
            DataType alpha, beta, ressquared, rho, rho_old, omega, sigma, tol1;

            // initialization step
            int iter = 0;
            bool restart = true;

            // restart option for stability purposes
            while ( restart )
            {
                restart = false;

                // clear data structures
                r.Zeros ( );
                r0.Zeros ( );
                p.Zeros ( );
                Ap.Zeros ( );
                s.Zeros ( );
                t.Zeros ( );
                h.Zeros ( );
                h1.Zeros ( );

                // compute residual
                // r0 = b - Ax
                this->op_->VectorMult ( *x, &r0 );
                r0.ScaleAdd ( b, static_cast < DataType > ( -1. ) );

                // r = r0
                r.CloneFrom ( r0 );
                ressquared = r.Dot ( r );

                // rho = rho_old = (r0, r0)
                rho = rho_old = ressquared;
                this->res_ = sqrt ( ressquared );

                // check for trivial convergence
                conv = this->control ( ).Check ( iter, this->res ( ) );

                if ( this->print_level_ > 0 )
                {
                    LOG_INFO ( "BiCGSTAB", " with right preconditioning" );
                    LOG_INFO ( "BiCGSTAB", " (re)starts with residual norm " << this->res_ );
                }

                // p0 = r0
                p.CopyFrom ( r );

                // main loop
                while ( conv == IterateControl::kIterate && !restart )
                {
                    // increment interation count
                    ++iter;

                    // compute Ap = A*P*p
                    h.Zeros ( );
                    this->precond_->ApplyPreconditioner ( p, &h );
                    this->op_->VectorMult ( h, &Ap );
                    // compute sigma = (Ap, r0)
                    sigma = Ap.Dot ( r0 );

                    // compute stability tolerance 1: tol1 = 5e-16*norm2(Ap)*norm2(r0)
                    tol1 = static_cast < DataType > ( 5.e-16 ) * Ap.Norm2 ( ) * r0.Norm2 ( );

                    //check for stability
                    if ( std::abs ( sigma ) > tol1 )
                    { // stable case
                        // alpha = rho / sigma
                        alpha = rho / sigma;

                        // s = r - alpha*Ap
                        s.CloneFrom ( r );
                        s.Axpy ( Ap, -alpha );

                        // Stability criterion 2: norm2(s) > 5e-16
                        if ( s.Norm2 ( ) > static_cast < DataType > ( 5.e-16 ) )
                        { // stable case
                            // t = A*P*s
                            h.Zeros ( );
                            this->precond_->ApplyPreconditioner ( s, &h );
                            this->op_->VectorMult ( h, &t );

                            // omega = (t, s) / (t, t)
                            omega = t.Dot ( s ) / t.Dot ( t );

                            // x = x + P*(alpha*p + omega*s)
                            h.Zeros ( );
                            h1.Zeros ( );
                            h1.Axpy ( p, alpha );
                            h1.Axpy ( s, omega );
                            this->precond_->ApplyPreconditioner ( h1, &h );
                            x->Axpy ( h, static_cast < DataType > ( 1. ) );

                            // r = s - omega*t
                            r.CloneFrom ( s );
                            r.Axpy ( t, -omega );

                            //check for convergence
                            ressquared = r.Dot ( r );
                            this->res_ = sqrt ( ressquared );
                            if ( this->print_level_ > 1 )
                            {
                                LOG_INFO ( "BiCGSTAB residual (iteration " << iter << ")", this->res_ );
                            }

                            conv = this->control ( ).Check ( iter, this->res ( ) );
                            if ( conv != IterateControl::kIterate )
                            {
                                break;
                            }

                            // rho = (r, r0);
                            rho = r.Dot ( r0 );

                            // beta = (alpha * rho) / (omega * rho_old)
                            beta = ( alpha * rho ) / ( omega * rho_old );

                            // rho_old = rho
                            rho_old = rho;

                            //p = r + beta*(p - omega*Ap)
                            p.ScaleAdd ( r, beta );
                            p.Axpy ( Ap, -omega * beta );
                        }
                        else
                        { // do stable update
                            // x = x + P*(alpha*p)
                            h1.Zeros ( );
                            h1.Axpy ( p, alpha );
                            h.Zeros ( );
                            this->precond_->ApplyPreconditioner ( h1, &h );
                            x->Axpy ( h, static_cast < DataType > ( 1. ) );

                            // r = s
                            r.CloneFrom ( s );

                            //check for convergence
                            ressquared = r.Dot ( r );
                            this->res_ = sqrt ( ressquared );
                            if ( this->print_level_ > 1 )
                            {
                                LOG_INFO ( "BiCGSTAB residual (iteration " << iter << ")", this->res_ );
                            }

                            conv = this->control ( ).Check ( iter, this->res ( ) );
                            if ( conv != IterateControl::kIterate )
                            {
                                break;
                            }
                        }
                    }
                    else
                    { // restart
                        restart = true;
                    }
                }
            }

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "BiCGSTAB", " with right preconditioning ended after " << iter << " iterations " );
                LOG_INFO ( "BiCGSTAB", " with residual norm " << this->res_ );
            }

            if ( conv == IterateControl::kFailureDivergenceTol ||
                 conv == IterateControl::kFailureMaxitsExceeded )
                return kSolverExceeded;
            else
                return kSolverSuccess;
        }

        /// Solve with left preconditioning.

        template<class LAD>
        LinearSolverState BiCGSTAB<LAD>::SolvePrecondLeft ( const VectorType& b, VectorType* x )
        {
            assert ( this->Method ( ) == "LeftPreconditioning" );

            assert ( this->op_ != NULL );
            assert ( this->precond_ != NULL );

            LOG_ERROR ( "BiCGSTAB::SolvePrecondLeft: Not implemented yet." );
            LOG_ERROR ( "Returning solver error...\n" );
            return kSolverError;
        }

        /// template instantiation
        template class BiCGSTAB<LADescriptorCoupledD>;
        template class BiCGSTAB<LADescriptorCoupledS>;
#ifdef WITH_HYPRE
        template class BiCGSTAB<LADescriptorHypreD>;
#endif

    } // namespace la
} // namespace hiflow
