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

/// @author Michael Schick

#include "polynomial_chaos/pc_gmres.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "common/log.h"
#include "linear_algebra/la_descriptor.h"

using namespace hiflow::la;

namespace hiflow
{
    namespace polynomialchaos
    {

        /// standard constructor

        template<class LAD>
        PCGMRES<LAD>::PCGMRES ( )
        : GMRES<LAD>( )
        {
        }

        /// destructor

        template<class LAD>
        PCGMRES<LAD>::~PCGMRES ( )
        {
            this->op_ = NULL;
            this->precond_ = NULL;
        }

        template<class LAD>
        void PCGMRES<LAD>::SetupResidualAssembler ( ResidualAssembler<LAD>& assembler )
        {
            res_assembler_ = &assembler;
        }

        /// Applies Givens rotation.
        /// @param cs cos(phi)
        /// @param sn sin(phi)
        /// @param dx first coordinate
        /// @param dy second coordinate

        template<class LAD>
        void PCGMRES<LAD>::ApplyPlaneRotation ( const DataType& cs, const DataType& sn,
                                                DataType* dx, DataType* dy ) const
        {
            DataType temp = cs * ( *dx ) + sn * ( *dy );
            *dy = -sn * ( *dx ) + cs * ( *dy );
            *dx = temp;
        }

        /// Generates Givens rotation.
        /// @param dx first coordinate
        /// @param dy second coordinate
        /// @param cs cos(phi)
        /// @param sn sin(phi)

        template<class LAD>
        void PCGMRES<LAD>::GeneratePlaneRotation ( const DataType& dx, const DataType& dy,
                                                   DataType* cs, DataType* sn ) const
        {
            if ( dy == static_cast < DataType > ( 0. ) )
            {
                *cs = static_cast < DataType > ( 1. );
                *sn = static_cast < DataType > ( 0. );
            }
            else if ( std::abs ( dy ) > std::abs ( dx ) )
            {
                assert ( dy != static_cast < DataType > ( 0. ) );
                DataType temp = dx / dy;
                assert ( sqrt ( static_cast < DataType > ( 1. ) + temp * temp ) !=
                         static_cast < DataType > ( 0. ) );
                *sn = static_cast < DataType > ( 1. ) /
                        sqrt ( static_cast < DataType > ( 1. ) + temp * temp );
                *cs = temp * ( *sn );
            }
            else
            {
                assert ( dx != static_cast < DataType > ( 0. ) );
                DataType temp = dy / dx;
                assert ( sqrt ( static_cast < DataType > ( 1. ) + temp * temp ) !=
                         static_cast < DataType > ( 0. ) );
                *cs = static_cast < DataType > ( 1. ) /
                        sqrt ( static_cast < DataType > ( 1. ) + temp * temp );
                *sn = temp * ( *cs );
            }
        }

        /// Solve without preconditioning.

        template<class LAD>
        LinearSolverState PCGMRES<LAD>::SolveNoPrecond ( const VectorType& b,
                                                         VectorType* x )
        {
            assert ( this->Method ( ) == "NoPreconditioning" );

            assert ( this->op_ != NULL );
            assert ( this->size_basis ( ) > 0 );

            IterateControl::State conv = IterateControl::kIterate;

            // Hessenberg matrix
            la::SeqDenseMatrix<DataType> H;
            H.Resize ( ( this->size_basis ( ) ), this->size_basis ( ) + 1 );

            // Allocate array of pointer for Krylov subspace basis
            VectorType** V = new VectorType*[this->size_basis ( ) + 1];
            for ( int i = 0; i < this->size_basis ( ) + 1; ++i )
            {
                V[i] = new VectorType;
                V[i]->CloneFromWithoutContent ( b );
            }

            VectorType w; // temporary vector
            w.CloneFromWithoutContent ( b );

            std::vector<DataType> g ( this->size_basis ( ) + 1 ); // rhs of least squares problem
            std::vector<DataType> cs ( this->size_basis ( ) + 1 ); // Givens rotations
            std::vector<DataType> sn ( this->size_basis ( ) + 1 ); // Givens rotations

            int iter = 0;

            // compute residual V[0] = b - Ax
            if ( solve_type_ == "Free" )
            {
                this->res_assembler_->assemble_residual ( b, x, *V[0] );
            }
            else
            {
                this->op_->VectorMult ( *x, V[0] );
                V[0]->ScaleAdd ( b, static_cast < DataType > ( -1. ) );
            }

            this->res_ = V[0]->Norm2 ( );
            conv = this->control ( ).Check ( iter, this->res ( ) );

            VectorType zero_rhs;
            zero_rhs.CloneFrom ( *x );
            zero_rhs.Zeros ( );

            LOG_INFO ( "PCGMRES", " without preconditioning" );
            LOG_INFO ( "PCGMRES", " starts with residual norm " << this->res_ );

            // main loop
            while ( conv == IterateControl::kIterate )
            {
                g.assign ( g.size ( ), static_cast < DataType > ( 0. ) ); // g = 0
                H.Zeros ( );

                assert ( this->res_ != static_cast < DataType > ( 0. ) );
                V[0]->Scale ( static_cast < DataType > ( 1. ) / this->res_ ); // norm residual
                g[0] = this->res_;

                for ( int j = 0; j < this->size_basis ( ); ++j )
                {
                    ++iter;

                    // w = Av_j
                    if ( solve_type_ == "Free" )
                    {
                        this->res_assembler_->assemble_residual ( zero_rhs, V[j], w );
                        w.Scale ( -1.0 );
                    }
                    else
                    {
                        this->op_->VectorMult ( *V[j], &w );
                    }

                    // -- start building Hessenberg matrix H --
                    // vectors in V are ONB of Krylov subspace K_i(A,V[0])
                    for ( int i = 0; i <= j; ++i )
                    {
                        H ( j, i ) = w.Dot ( *V[i] );
                        w.Axpy ( *V[i], static_cast < DataType > ( -1. ) * H ( j, i ) );
                    }

                    H ( j, j + 1 ) = w.Norm2 ( );
                    assert ( H ( j, j + 1 ) != static_cast < DataType > ( 0. ) );

                    w.Scale ( static_cast < DataType > ( 1. ) / H ( j, j + 1 ) );
                    V[j + 1]->CopyFrom ( w );
                    // -- end building Hessenberg matrix H --

                    // apply old Givens rotation on old H entries
                    for ( int k = 0; k < j; ++k )
                        this->ApplyPlaneRotation ( cs[k], sn[k], &H ( j, k ), &H ( j, k + 1 ) );

                    // determine new Givens rotation for actual iteration i
                    this->GeneratePlaneRotation ( H ( j, j ), H ( j, j + 1 ), &cs[j], &sn[j] );

                    // apply Givens rotation on new H element
                    this->ApplyPlaneRotation ( cs[j], sn[j], &H ( j, j ), &H ( j, j + 1 ) );

                    // update g for next dimension -> g[j+1] is norm of actual residual
                    this->ApplyPlaneRotation ( cs[j], sn[j], &g[j], &g[j + 1] );

                    this->res_ = std::abs ( g[j + 1] );
                    conv = this->control ( ).Check ( iter, this->res_ );

                    if ( conv != IterateControl::kIterate )
                    {
                        this->UpdateSolution ( V, H, g, j, x );
                        break;
                    }
                } // for (int j = 0; j < this->size_basis(); ++j)

                // setup for restart
                if ( conv == IterateControl::kIterate )
                {
                    // -> update solution
                    this->UpdateSolution ( V, H, g, this->size_basis_ - 1, x ); // x = x + Vy
                    // -> compute residual Ax-b
                    if ( solve_type_ == "Free" )
                    {
                        this->res_assembler_->assemble_residual ( b, x, *V[0] );
                    }
                    else
                    {
                        this->op_->VectorMult ( *x, V[0] );
                        V[0]->ScaleAdd ( b, static_cast < DataType > ( -1. ) );
                    }
                }
            } // while (conv == IterateControl::kIterate)

            // deallocate Krylov subspace basis V
            for ( int j = 0; j < this->size_basis_ + 1; ++j )
            {
                delete V[j];
            }
            delete[] V;

            this->iter_ = iter;

            LOG_INFO ( "PCGMRES", " without preconditioning ended after " << iter << " iterations " );
            LOG_INFO ( "PCGMRES", " with residual norm " << this->res_ );

            if ( conv == IterateControl::kFailureDivergenceTol ||
                 conv == IterateControl::kFailureMaxitsExceeded )
                return kSolverExceeded;
            else
                return kSolverSuccess;
        }

        /// Solve with right preconditioning.

        template<class LAD>
        LinearSolverState PCGMRES<LAD>::SolveRight ( const VectorType& b, VectorType* x )
        {
            assert ( this->Method ( ) == "RightPreconditioning" );

            assert ( this->op_ != NULL );
            assert ( this->precond_ != NULL );
            assert ( this->size_basis ( ) > 0 );

            if ( !this->precond_->GetState ( ) || !this->precond_->GetReuse ( ) )
            {
                this->BuildPreconditioner ( );
            }

            IterateControl::State conv = IterateControl::kIterate;

            // Hessenberg matrix
            SeqDenseMatrix<DataType> H;
            H.Resize ( ( this->size_basis ( ) ), this->size_basis ( ) + 1 );

            // Allocate array of pointer for Krylov subspace basis
            VectorType** V = new VectorType*[this->size_basis ( ) + 1];
            for ( int i = 0; i < this->size_basis ( ) + 1; ++i )
            {
                V[i] = new VectorType;
                V[i]->CloneFromWithoutContent ( b );
            }

            VectorType w, z;
            w.CloneFromWithoutContent ( b );
            z.CloneFromWithoutContent ( b );

            std::vector<DataType> g ( this->size_basis ( ) + 1 ); // rhs of least squares problem
            std::vector<DataType> cs ( this->size_basis ( ) + 1 ); // Givens rotations
            std::vector<DataType> sn ( this->size_basis ( ) + 1 ); // Givens rotations

            int iter = 0;

            // compute residual V[0] = b - Ax
            if ( solve_type_ == "Free" )
            {
                this->res_assembler_->assemble_residual ( b, x, *V[0] );
            }
            else
            {
                this->op_->VectorMult ( *x, V[0] );
                V[0]->ScaleAdd ( b, static_cast < DataType > ( -1. ) );
            }

            this->res_ = V[0]->Norm2 ( );
            conv = this->control ( ).Check ( iter, this->res ( ) );

            LOG_INFO ( "PCGMRES", " with right preconditioning" );
            LOG_INFO ( "PCGMRES", " starts with residual norm " << this->res_ );

            VectorType zero_rhs;
            zero_rhs.CloneFrom ( *x );
            zero_rhs.Zeros ( );

            // main loop
            while ( conv == IterateControl::kIterate )
            {
                g.assign ( g.size ( ), static_cast < DataType > ( 0. ) ); // g = 0
                H.Zeros ( );

                assert ( this->res_ != static_cast < DataType > ( 0. ) );
                V[0]->Scale ( static_cast < DataType > ( 1. ) / this->res_ ); // norm residual
                g[0] = this->res_;

                for ( int j = 0; j < this->size_basis ( ); ++j )
                {
                    ++iter;

                    this->precond_->ApplyPreconditioner ( *V[j], &z ); // z ~= M^-1 v_j
                    // w = A z
                    if ( solve_type_ == "Free" )
                    {
                        this->res_assembler_->assemble_residual ( zero_rhs, &z, w );
                        w.Scale ( -1.0 );
                    }
                    else
                    {
                        this->op_->VectorMult ( z, &w );
                    }

                    // -- start building Hessenberg matrix H --
                    // vectors in V are ONB of Krylov subspace K_i(A,V[0])
                    for ( int i = 0; i <= j; ++i )
                    {
                        H ( j, i ) = w.Dot ( *V[i] );
                        w.Axpy ( *V[i], static_cast < DataType > ( -1. ) * H ( j, i ) );
                    }

                    H ( j, j + 1 ) = w.Norm2 ( );
                    assert ( H ( j, j + 1 ) != static_cast < DataType > ( 0. ) );

                    w.Scale ( static_cast < DataType > ( 1. ) / H ( j, j + 1 ) );
                    V[j + 1]->CopyFrom ( w );
                    // -- end building Hessenberg matrix H --

                    // apply old Givens rotation on old H entries
                    for ( int k = 0; k < j; ++k )
                        this->ApplyPlaneRotation ( cs[k], sn[k], &H ( j, k ), &H ( j, k + 1 ) );

                    // determine new Givens rotation for actual iteration i
                    this->GeneratePlaneRotation ( H ( j, j ), H ( j, j + 1 ), &cs[j], &sn[j] );

                    // apply Givens rotation on new H element
                    this->ApplyPlaneRotation ( cs[j], sn[j], &H ( j, j ), &H ( j, j + 1 ) );

                    // update g for next dimension -> g[j+1] is norm of actual residual
                    this->ApplyPlaneRotation ( cs[j], sn[j], &g[j], &g[j + 1] );

                    this->res_ = std::abs ( g[j + 1] );
                    conv = this->control ( ).Check ( iter, this->res_ );

                    if ( conv != IterateControl::kIterate )
                    {
                        z.Zeros ( );
                        this->UpdateSolution ( V, H, g, j, &z );
                        this->precond_->ApplyPreconditioner ( z, &w );
                        x->Axpy ( w, static_cast < DataType > ( 1. ) );
                        break;
                    }
                } // for (int j = 0; j < this->size_basis(); ++j)

                // setup for restart
                if ( conv == IterateControl::kIterate )
                {
                    // -> update solution
                    z.Zeros ( );
                    this->UpdateSolution ( V, H, g, this->size_basis_ - 1, &z );
                    this->precond_->ApplyPreconditioner ( z, &w );
                    x->Axpy ( w, static_cast < DataType > ( 1. ) );
                    // -> compute residual Ax-b
                    if ( solve_type_ == "Free" )
                    {
                        this->res_assembler_->assemble_residual ( b, x, *V[0] );
                    }
                    else
                    {
                        this->op_->VectorMult ( *x, V[0] );
                        V[0]->ScaleAdd ( b, static_cast < DataType > ( -1. ) );
                    }
                }
            } // while (conv == IterateControl::kIterate)

            // deallocate Krylov subspace basis V
            for ( int j = 0; j < this->size_basis_ + 1; ++j )
            {
                delete V[j];
            }
            delete[] V;

            this->iter_ = iter;

            LOG_INFO ( "PCGMRES", " with right preconditioning ended after " << iter << " iterations " );
            LOG_INFO ( "PCGMRES", " with residual norm " << this->res_ );

            if ( conv == IterateControl::kFailureDivergenceTol ||
                 conv == IterateControl::kFailureMaxitsExceeded )
                return kSolverExceeded;
            else
                return kSolverSuccess;

        }

        /// Updates solution: x = x + Vy with y solution of least squares problem.
        /// @param V Krylov subspace basis
        /// @param H Hessenberg matrix
        /// @param g rhs of least squares problem
        /// @param k iteration step
        /// @param x solution vector

        template<class LAD>
        void PCGMRES<LAD>::UpdateSolution ( const VectorType * const * V,
                                            const la::SeqDenseMatrix<DataType>& H,
                                            const std::vector<DataType>& g,
                                            int k,
                                            VectorType* x ) const
        {
            std::vector<DataType> y ( g );

            // back substitution
            for ( int i = k; i >= 0; --i )
            {
                assert ( H ( i, i ) != static_cast < DataType > ( 0. ) );
                y[i] /= H ( i, i );

                const DataType temp = y[i];
                for ( int j = i - 1; j >= 0; --j )
                    y[j] -= H ( i, j ) * temp;
            }

            // compute solution
            for ( int j = 0; j <= k; ++j )
                x->Axpy ( *V[j], y[j] );
        }

        /// template instantiation
        template class PCGMRES<la::LADescriptorPolynomialChaosD>;

    } // namespace
} // namespace
