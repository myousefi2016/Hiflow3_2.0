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

/// @author Hendryk Bockelmann, Chandramowli Subramanian

#include "fgmres.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "common/log.h"
#include "linear_algebra/la_descriptor.h"
#include "linear_algebra/seq_dense_matrix.h"

namespace hiflow
{
    namespace la
    {

        /// standard constructor

        template<class LAD>
        FGMRES<LAD>::FGMRES ( )
        : GMRES<LAD>( )
        {
            this->size_basis_ = 0;
            this->SetMethod ( "RightPreconditioning" );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Linear solver", "FGMRES" );
            }
        }

        /// destructor

        template<class LAD>
        FGMRES<LAD>::~FGMRES ( )
        {
            this->op_ = NULL;
            this->precond_ = NULL;
        }

        /// Solve with left preconditioning.

        template<class LAD>
        LinearSolverState FGMRES<LAD>::SolveLeft ( const VectorType& b, VectorType* x )
        {
            assert ( this->Method ( ) == "LeftPreconditioning" );

            assert ( this->op_ != NULL );
            assert ( this->precond_ != NULL );
            assert ( this->size_basis ( ) > 0 );

            LOG_ERROR ( "FGMRES::SolveLeft: Not implemented yet." );
            LOG_ERROR ( "Returning solver error..." );
            return kSolverError;
        }

        /// Solve with right preconditioning.

        template<class LAD>
        LinearSolverState FGMRES<LAD>::SolveRight ( const VectorType& b, VectorType* x )
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

            // compute really used basis size as minimum of maximum iterations and
            // given basis size
            const int basis_size_actual = std::min ( this->size_basis ( ), this->control ( ).maxits ( ) );

            // Hessenberg matrix
            SeqDenseMatrix<DataType> H;
            H.Resize ( basis_size_actual, basis_size_actual + 1 );

            // Allocate array of pointer for Krylov subspace basis
            VectorType** V = new VectorType*[basis_size_actual + 1];
            for ( int i = 0; i != basis_size_actual + 1; ++i )
            {
                V[i] = new VectorType;
            }
            // preconditioned Krylov subspace basis
            VectorType** Z = new VectorType*[basis_size_actual + 1];
            for ( int i = 0; i != basis_size_actual + 1; ++i )
            {
                Z[i] = new VectorType;
            }

            VectorType w;
            w.CloneFromWithoutContent ( b );

            std::vector<DataType> g ( basis_size_actual + 1 ); // rhs of least squares problem
            std::vector<DataType> cs ( basis_size_actual + 1 ); // Givens rotations
            std::vector<DataType> sn ( basis_size_actual + 1 ); // Givens rotations

            int iter = 0;

            // compute residual V[0] = b - Ax
            V[0]->CloneFromWithoutContent ( b );
            this->op_->VectorMult ( *x, V[0] );
            V[0]->ScaleAdd ( b, static_cast < DataType > ( -1. ) );

            if ( this->filter_solution_ )
            {
                V[0]->Update ( );
                this->non_lin_op_->ApplyFilter ( *V[0] );
            }

            this->res_ = V[0]->Norm2 ( );
            conv = this->control ( ).Check ( iter, this->res_ );

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "FGMRES", " with right preconditioning" );
                LOG_INFO ( "FGMRES", " starts with residual norm " << this->res_ );
            }

            // main loop
            while ( conv == IterateControl::kIterate )
            {
                g.assign ( g.size ( ), static_cast < DataType > ( 0. ) ); // g = 0
                H.Zeros ( );

                assert ( this->res_ != static_cast < DataType > ( 0. ) );
                V[0]->Scale ( static_cast < DataType > ( 1. ) / this->res_ ); // norm residual
                g[0] = this->res_;

                for ( int j = 0; j != basis_size_actual; ++j )
                {
                    ++iter;
                    Z[j]->CloneFromWithoutContent ( b );
                    //Z[j]->Zeros();
                    this->precond_->ApplyPreconditioner ( *V[j], Z[j] ); // Z[j] ~= M^-1 v_j
                    if ( this->filter_solution_ )
                    {
                        Z[j]->Update ( );
                        this->non_lin_op_->ApplyFilter ( *Z[j] );
                    }
                    this->op_->VectorMult ( *Z[j], &w ); // w = A Z[j]
                    if ( this->filter_solution_ )
                    {
                        w.Update ( );
                        this->non_lin_op_->ApplyFilter ( w );
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

                    V[j + 1]->CloneFrom ( w );
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
                    if ( this->print_level_ > 1 )
                    {
                        LOG_INFO ( "FGMRES residual (iteration " << iter << ")", this->res_ );
                    }

                    conv = this->control ( ).Check ( iter, this->res_ );

                    if ( conv != IterateControl::kIterate )
                    {
                        this->UpdateSolution ( Z, H, g, j, x ); // x = x + Zy
                        if ( this->filter_solution_ )
                        {
                            x->Update ( );
                            this->non_lin_op_->ApplyFilter ( *x );
                        }
                        break;
                    }
                } // for (int j = 0; j < this->size_basis(); ++j)

                // setup for restart
                if ( conv == IterateControl::kIterate )
                {
                    this->UpdateSolution ( Z, H, g, basis_size_actual - 1, x ); // x = x + Zy
                    if ( this->filter_solution_ )
                    {
                        x->Update ( );
                        this->non_lin_op_->ApplyFilter ( *x );
                    }

                    this->op_->VectorMult ( *x, V[0] );
                    V[0]->Scale ( static_cast < DataType > ( -1. ) );
                    V[0]->Axpy ( b, static_cast < DataType > ( 1. ) );
                    if ( this->filter_solution_ )
                    {
                        V[0]->Update ( );
                        this->non_lin_op_->ApplyFilter ( *V[0] );
                    }
                }
            } // while (conv == IterateControl::kIterate)

            // deallocate Krylov subspace basis V and Z
            for ( int j = 0; j != basis_size_actual + 1; ++j )
            {
                delete V[j];
                delete Z[j];
            }
            delete[] V;
            delete[] Z;

            this->iter_ = iter;

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "FGMRES", " with right preconditioning ended after " << iter << " iterations " );
                LOG_INFO ( "FGMRES", " with residual norm " << this->res_ );
            }

            if ( conv == IterateControl::kFailureDivergenceTol ||
                 conv == IterateControl::kFailureMaxitsExceeded )
                return kSolverExceeded;
            else
                return kSolverSuccess;
        }

        /// template instantiation
        template class FGMRES<LADescriptorCoupledD>;
        template class FGMRES<LADescriptorCoupledS>;
#ifdef WITH_HYPRE
        template class FGMRES<LADescriptorHypreD>;
        template class FGMRES<LADescriptorPolynomialChaosExpansionD>;
#endif

    } // namespace la
} // namespace hiflow
