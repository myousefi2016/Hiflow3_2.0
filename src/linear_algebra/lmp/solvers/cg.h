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
/// @date 2015-07-23

#ifndef HIFLOW_LINEAR_ALGEBRA_LMP_SOLVERS_CG_H_
#    define HIFLOW_LINEAR_ALGEBRA_LMP_SOLVERS_CG_H_

#    include "../lpreconditioner.h"
#    include "../lmp_mem.h"
#    include "ScopedTimer.h"

namespace hiflow
{
    namespace la
    {

        /// Function wrapper calling CGFunctor for backward compatibility
        template <typename ValueType>
        int cg ( lVector<ValueType> *x, const lVector<ValueType> *b,
                 const lMatrix<ValueType> *matrix, const ValueType rel_eps,
                 const ValueType abs_eps, int &max_iter, const int print,
                 lPreconditioner<ValueType> *lPrecond = NULL );

        /// CGFunctor doing the real work

        template <typename ValueType, bool Timer>
        struct CGFunctor
        {
            // Default Constructor

            CGFunctor ( )
            {
                if ( Timer )
                {
                    timeNrm2 = new ScopedTimer ( "conjugate gradient, Nrm2", false );
                    timeDot = new ScopedTimer ( "conjugate gradient, Dot", false );
                    timeAxpy = new ScopedTimer ( "conjugate gradient, Axpy", false );
                    timeScaleAdd = new ScopedTimer ( "conjugate gradient, ScaleAdd", false );
                    timeVectorMult = new ScopedTimer ( "conjugate gradient, VectorMult", false );
                    timeTotal = new ScopedTimer ( "conjugate gradient, Total", false );
                    timePreconditioner =
                            new ScopedTimer ( "conjugate gradient, Preconditioner", false );
                }
            }

            // Destructor

            ~CGFunctor ( )
            {
                if ( Timer )
                {
                    delete timeNrm2;
                    delete timeDot;
                    delete timeAxpy;
                    delete timeScaleAdd;
                    delete timeVectorMult;
                    delete timeTotal;
                    delete timePreconditioner;
                }
            }

            int operator() ( lVector<ValueType> *x, const lVector<ValueType> *b,
                    const lMatrix<ValueType> *matrix, const ValueType rel_eps,
                    const ValueType abs_eps, int &max_iter, const int print,
                    lPreconditioner<ValueType> *lPrecond = NULL ) const;

          private:

            void resetAll ( ) const
            {
                reset<Timer>( timeNrm2 );
                reset<Timer>( timeDot );
                reset<Timer>( timeAxpy );
                reset<Timer>( timeScaleAdd );
                reset<Timer>( timeVectorMult );
                reset<Timer>( timeTotal );
                reset<Timer>( timePreconditioner );
            }

            ScopedTimer *timeNrm2;
            ScopedTimer *timeDot;
            ScopedTimer *timeAxpy;
            ScopedTimer *timeScaleAdd;
            ScopedTimer *timeVectorMult;
            ScopedTimer *timeTotal;
            ScopedTimer *timePreconditioner;
        };

        template <typename ValueType, bool Timer>
        int CGFunctor<ValueType, Timer>::operator() (
                lVector<ValueType> *x, const lVector<ValueType> *b,
                const lMatrix<ValueType> *matrix, const ValueType rel_eps,
                const ValueType abs_eps, int &max_iter, const int print,
                lPreconditioner<ValueType> *lPrecond ) const
        {
            resetAll ( );

            start<Timer>( timeTotal );

            if ( print > 0 )
            {
                LOG_INFO ( "pcg", "Platform input data:" );
                x->print ( );
                b->print ( );
                matrix->print ( );
            }

            lVector<ValueType> *p = x->CloneWithoutContent ( );
            lVector<ValueType> *q = x->CloneWithoutContent ( );
            lVector<ValueType> *r = x->CloneWithoutContent ( );

            start<Timer>( timeVectorMult );
            matrix->VectorMult ( *x, r );
            stop<Timer>( timeVectorMult );

            start<Timer>( timeScaleAdd );
            r->ScaleAdd ( -1.0, *b ); // r = b - Ax
            stop<Timer>( timeScaleAdd );

            lVector<ValueType> *z = NULL;
            if ( lPrecond )
            {
                z = x->CloneWithoutContent ( );

                start<Timer>( timePreconditioner );
                lPrecond->ApplylPreconditioner ( *r, z );
                stop<Timer>( timePreconditioner );

                lPrecond->print ( );
            }

            start<Timer>( timeDot );
            ValueType rho = r->Dot ( *r );
            stop<Timer>( timeDot );

            start<Timer>( timeNrm2 );
            ValueType norm = std::sqrt ( rho );
            stop<Timer>( timeNrm2 );

            if ( print > -1 )
            {
                LOG_INFO ( "cg", "CG with " << ( lPrecond ? "" : "No" ) << " Preconditioner" );
                LOG_INFO ( "cg", "Init, rho = " << rho << ", norm = " << norm );
            }

            int iter = 0;
            ValueType beta = 0.0;
            ValueType rho_old = 0.0;
            ValueType conv_crit = std::max ( rel_eps * norm, abs_eps );
            for (; iter < max_iter; ++iter )
            {
                start<Timer>( timePreconditioner );
                if ( lPrecond ) lPrecond->ApplylPreconditioner ( *r, z );
                stop<Timer>( timePreconditioner );

                start<Timer>( timeScaleAdd );
                p->ScaleAdd ( beta, *r ); // p = beta*p + r
                stop<Timer>( timeScaleAdd );

                start<Timer>( timeVectorMult );
                matrix->VectorMult ( *p, q );
                stop<Timer>( timeVectorMult );

                start<Timer>( timeDot );
                ValueType alpha = rho / ( p->Dot ( *q ) );
                stop<Timer>( timeDot );

                start<Timer>( timeAxpy );
                x->Axpy ( *p, alpha );
                stop<Timer>( timeAxpy );

                start<Timer>( timeAxpy );
                r->Axpy ( *q, -alpha );
                stop<Timer>( timeAxpy );

                rho_old = rho;

                start<Timer>( timeDot );
                rho = r->Dot ( *r );
                stop<Timer>( timeDot );

                beta = rho / rho_old;

                start<Timer>( timeNrm2 );
                norm = std::sqrt ( rho );
                stop<Timer>( timeNrm2 );

                if ( print > 0 )
                    LOG_INFO ( "cg", "iter = " << iter << ", rho = " << rho
                               << ", norm = " << norm );

                // Convergence reached?
                if ( norm <= conv_crit ) break;

            } // iter

            if ( print > -1 )
                LOG_INFO ( "cg", "End after: iterations=" << iter << ", rho = " << rho
                           << ", norm = " << norm );
            max_iter = iter;

            delete p;
            delete q;
            delete r;

            stop<Timer>( timeTotal );
            return 0;
        }

    } // namespace hiflow
} // namespace la

#endif /* HIFLOW_LINEAR_ALGEBRA_LMP_SOLVERS_CG_H_ */
