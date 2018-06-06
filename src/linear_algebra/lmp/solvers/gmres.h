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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-08-06

#ifndef HIFLOW_LINEAR_ALGEBRA_LMP_SOLVERS_GMRES_H_
#    define HIFLOW_LINEAR_ALGEBRA_LMP_SOLVERS_GMRES_H_

#    include "../lpreconditioner.h"
#    include "../lmp_mem.h"
#    include "ScopedTimer.h"

namespace hiflow
{
    namespace la
    {

        /// Generalized minimal residual method (GMRES)

        template <typename ValueType, bool Timer = false >
        struct GMRESFunctor
        {
            // Default Constructor

            GMRESFunctor ( int basis_size ) : basis_size ( basis_size )
            {
                if ( Timer )
                {
                    timeTotal = new ScopedTimer ( "conjugate gradient, Total", false );
                    timePreconditioner =
                            new ScopedTimer ( "conjugate gradient, Preconditioner", false );
                }
            }

            // Destructor

            ~GMRESFunctor ( )
            {
                if ( Timer )
                {
                    delete timeTotal;
                    delete timePreconditioner;
                }
            }

            int operator() ( lVector<ValueType> *x, const lVector<ValueType> *b,
                    const lMatrix<ValueType> *matrix, const ValueType rel_eps,
                    const ValueType abs_eps, int &max_iter, const int print,
                    lPreconditioner<ValueType> *lPrecond = NULL ) const;

          protected:
            /// Applies Givens rotation.
            /// @param cs cos(phi)
            /// @param sn sin(phi)
            /// @param dx first coordinate
            /// @param dy second coordinate
            void apply_plane_rotation ( ValueType cs, ValueType sn, ValueType *dx,
                                        ValueType *dy ) const;

            /// Generates Givens rotation.
            /// @param dx first coordinate
            /// @param dy second coordinate
            /// @param cs cos(phi)
            /// @param sn sin(phi)
            void generate_plane_rotation ( ValueType dx, ValueType dy, ValueType *cs,
                                           ValueType *sn ) const;

            /// Updates solution: x = x + Vy with y solution of least squares problem.
            /// @param V Krylov subspace basis
            /// @param H Hessenberg matrix
            /// @param g rhs of least squares problem
            /// @param k iteration step
            /// @param x solution vector
            void update_solution ( std::vector<lVector<ValueType> *> const &v,
                                   const SeqDenseMatrix<ValueType> &H,
                                   const std::vector<ValueType> &g, int k,
                                   lVector<ValueType> *x ) const;

            void reset_all ( ) const
            {
                reset<Timer>( timeTotal );
                reset<Timer>( timePreconditioner );
            }

            ScopedTimer *timeTotal;
            ScopedTimer *timePreconditioner;

            int basis_size;
        };

        template <typename ValueType, bool Timer>
        int GMRESFunctor<ValueType, Timer>::operator() (
                lVector<ValueType> *x, const lVector<ValueType> *b,
                const lMatrix<ValueType> *matrix, const ValueType rel_eps,
                const ValueType abs_eps, int &max_iter, const int print,
                lPreconditioner<ValueType> *lPrecond ) const
        {
            reset_all ( );
            start<Timer>( timeTotal );

            lVector<ValueType> *w = b->CloneWithoutContent ( );

            std::vector<ValueType> g ( basis_size + 1, 0. ); // rhs of least square problem
            std::vector<ValueType> cs ( basis_size + 1, 0. ); // Givens rotations
            std::vector<ValueType> sn ( basis_size + 1, 0. ); // Givens rotations

            // Hessenberg matrix
            SeqDenseMatrix<ValueType> H;
            H.Resize ( basis_size, basis_size + 1 );

            // Krylov subspace basis
            std::vector<lVector<ValueType> *> v ( basis_size + 1 );
            for ( int i = 0; i < basis_size + 1; ++i )
            {
                v[i] = b->CloneWithoutContent ( );
            }

            // Compute residual v[0] = b - Ax
            matrix->VectorMult ( *x, v[0] );
            v[0]->ScaleAdd ( -1.0, *b );

            ValueType norm = v[0]->Norm2 ( );

            if ( print > -1 )
            {
                LOG_INFO ( "GMRESFunctor", " with " << ( lPrecond ? "" : "no" )
                           << " preconditioner" );
                LOG_INFO ( "GMRESFunctor", " Init, norm = " << norm );
            }

            int iter = 0;
            ValueType conv_crit = std::max ( rel_eps * norm, abs_eps );
            for (; iter < max_iter; )
            {
                g.assign ( basis_size + 1, 0.0 );
                H.Zeros ( );

                v[0]->Scale ( 1.0 / norm ); // norm residual
                g[0] = norm;

                for ( int j = 0; j < basis_size; ++j )
                {
                    ++iter;

                    matrix->VectorMult ( *v[j], w );

                    // Start building Hessenberg matrix H
                    // vectors in V are ONB of Krylov subspace K_i(A,V[0])
                    for ( int i = 0; i <= j; ++i )
                    {
                        H ( j, i ) = w->Dot ( *v[i] );
                        w->Axpy ( *v[i], -1.0 * H ( j, i ) );
                    }

                    H ( j, j + 1 ) = w->Norm2 ( );

                    w->Scale ( 1.0 / H ( j, j + 1 ) );
                    v[j + 1]->CopyFrom ( *w );
                    // End building Hessenberg matrix H

                    // apply old Givens rotation on old H entries
                    for ( int k = 0; k < j; ++k )
                        this->apply_plane_rotation ( cs[k], sn[k], &H ( j, k ), &H ( j, k + 1 ) );

                    // determine new Givens rotation for actual iteration i
                    this->generate_plane_rotation ( H ( j, j ), H ( j, j + 1 ), &cs[j], &sn[j] );

                    // apply Givens rotation on new H element
                    this->apply_plane_rotation ( cs[j], sn[j], &H ( j, j ), &H ( j, j + 1 ) );

                    // update g for next dimension -> g[j+1] is norm of actual residual
                    this->apply_plane_rotation ( cs[j], sn[j], &g[j], &g[j + 1] );

                    norm = std::abs ( g[j + 1] );

                    if ( print > 0 )
                        LOG_INFO ( "gmres", "iteration=" << iter << ", norm = " << norm );

                    // Convergence reached?
                    if ( norm <= conv_crit )
                    {
                        this->update_solution ( v, H, g, j, x );
                        break;
                    }
                } // j

                // Convergence reached?
                if ( norm <= conv_crit ) break;

            } // iter

            // Free Krylov subspace basis
            for ( int j = 0; j < basis_size + 1; ++j ) delete v[j];

            if ( print > -1 )
                LOG_INFO ( "gmres", "End after: iterations=" << iter << ", norm = " << norm );
            max_iter = iter;

            stop<Timer>( timeTotal );
            return 0;
        }

        template <typename ValueType, bool Timer>
        void GMRESFunctor<ValueType, Timer>::apply_plane_rotation ( ValueType cs,
                                                                    ValueType sn,
                                                                    ValueType *dx,
                                                                    ValueType *dy ) const
        {
            ValueType temp = cs * ( *dx ) + sn * ( *dy );
            *dy = -sn * ( *dx ) + cs * ( *dy );
            *dx = temp;
        }

        template <typename ValueType, bool Timer>
        void GMRESFunctor<ValueType, Timer>::generate_plane_rotation (
                                                                       ValueType dx, ValueType dy, ValueType *cs, ValueType *sn ) const
        {
            if ( dy == 0.0 )
            {
                *cs = 1.0;
                *sn = 0.0;
            }
            else if ( std::abs ( dy ) > std::abs ( dx ) )
            {
                assert ( dy != 0.0 );
                ValueType temp = dx / dy;
                assert ( sqrt ( 1.0 + temp * temp ) != 0.0 );
                *sn = 1.0 / sqrt ( 1.0 + temp * temp );
                *cs = temp * ( *sn );
            }
            else
            {
                assert ( dx != 0.0 );
                ValueType temp = dy / dx;
                assert ( sqrt ( 1.0 + temp * temp ) != 0.0 );
                *cs = 1.0 / sqrt ( 1.0 + temp * temp );
                *sn = temp * ( *cs );
            }
        }

        template <typename ValueType, bool Timer>
        void GMRESFunctor<ValueType, Timer>::update_solution (
                                                               std::vector<lVector<ValueType> *> const &v,
                                                               const SeqDenseMatrix<ValueType> &H, const std::vector<ValueType> &g, int k,
                                                               lVector<ValueType> *x ) const
        {
            std::vector<ValueType> y ( g );

            // back substitution
            for ( int i = k; i >= 0; --i )
            {
                assert ( H ( i, i ) != 0.0 );
                y[i] /= H ( i, i );

                const ValueType temp = y[i];
                for ( int j = i - 1; j >= 0; --j ) y[j] -= H ( i, j ) * temp;
            }

            // compute solution
            for ( int j = 0; j <= k; ++j ) x->Axpy ( *v[j], y[j] );
        }

    } // namespace hiflow
} // namespace la

#endif /* HIFLOW_LINEAR_ALGEBRA_LMP_SOLVERS_GMRES_H_ */
