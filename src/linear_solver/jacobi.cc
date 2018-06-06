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

#include "linear_solver/jacobi.h"
#include "linear_algebra/la_descriptor.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        Jacobi<LAD>::Jacobi ( )
        : LinearSolver<LAD>( ),
        solve_mode_ ( 0 ),
        smooth_iter_ ( 1 ),
        inner_iter_ ( 1 ),
        w_ ( 1.0 ),
        async_ ( false ),
        inv_diag_ ( 0 ),
        csr_mat_ ( 0 ),
        y_ ( 0 )
        {
        }

        template<class LAD>
        Jacobi<LAD>::~Jacobi ( )
        {
            if ( inv_diag_ != 0 ) delete inv_diag_;
            if ( y_ != 0 ) delete y_;
        }

        template<class LAD>
        void Jacobi<LAD>::Prepare ( OperatorType& op, const VectorType& sample )
        {
            // try to cast to a CSR matrix
            csr_mat_ = dynamic_cast < CSR_lMatrix<DataType>* > ( &( op.diagonal ( ) ) );

            // for CSR matrices: shift diagonal elements to front in matrix values array
            if ( csr_mat_ )
            {
                csr_mat_->SwapDiagElementsToRowFront ( );
            }

            // delete old inv_diag if existent
            if ( inv_diag_ != 0 ) delete inv_diag_;
            if ( y_ != 0 ) delete y_;

            // create local vectors from sample
            inv_diag_ = sample.interior ( ).CloneWithoutContent ( );
            y_ = sample.interior ( ).CloneWithoutContent ( );

            op.diagonal ( ).extract_invdiagelements ( 0, op.nrows_local ( ), inv_diag_ );

            this->SetupOperator ( op );
        }

        template<class LAD>
        void Jacobi<LAD>::SetSolveMode ( const int mode )
        {
            assert ( mode >= 0 );
            assert ( mode <= 3 );

            // 0: solve normal
            // 1: solve damped
            // 2: smooth normal
            // 3 smooth damped
            solve_mode_ = mode;
        }

        template<class LAD>
        void Jacobi<LAD>::SetDampingParameter ( const DataType w )
        {
            assert ( w > 0.0 );

            w_ = w;
        }

        template<class LAD>
        void Jacobi<LAD>::SetNumIter ( const int iter )
        {
            assert ( iter > 0 );

            smooth_iter_ = iter;
        }

        template<class LAD>
        void Jacobi<LAD>::SetInnerIter ( const int iter )
        {
            assert ( iter > 0 );

            inner_iter_ = iter;
        }

        template<class LAD>
        LinearSolverState Jacobi<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
            assert ( x != 0 );

            switch ( solve_mode_ )
            {
                case 0:
                    return SolveNormal ( b, x );
                    break;

                case 1:
                    return SolveDamped ( b, x );
                    break;

                case 2:
                    return SmoothNormal ( b, x );
                    break;

                case 3:
                    return SmoothDamped ( b, x );
                    break;

                default:
                    LOG_ERROR ( "Jacobi: unknown solve mode" );
                    exit ( -1 );
                    break;
            }
            return kSolverExceeded;
        }

        template<class LAD>
        LinearSolverState Jacobi<LAD>::SolveNormal ( const VectorType& b, VectorType* x )
        {
            assert ( x != 0 );
            assert ( this->op_ != 0 );
            assert ( inv_diag_ != 0 );
            assert ( y_ != 0 );

            this->iter_ = 0;

            // create residual vector
            VectorType res;
            res.CloneFromWithoutContent ( b );

            // compute residual
            this->op_->VectorMult ( *x, &res );
            res.ScaleAdd ( b, static_cast < DataType > ( -1.0 ) );
            this->res_ = res.Norm2 ( );
            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "Jacobi", " ================================================" );
                LOG_INFO ( "Jacobi", " SolveNormal starts with residual norm " << this->res_ );
            }

            IterateControl::State iter_ctrl = this->control ( ).Check ( this->iter_, this->res_ );

            while ( iter_ctrl == IterateControl::kIterate )
            {
                // do iteration
                if ( async_ ) x->UpdateGhost ( );
                for ( int i = 0; i < inner_iter_; ++i )
                {
                    DoNormalIteration ( b, x );
                }

                // compute residual
                this->op_->VectorMult ( *x, &res );
                res.ScaleAdd ( b, static_cast < DataType > ( -1.0 ) );
                this->res_ = res.Norm2 ( );

                // check
                iter_ctrl = this->control ( ).Check ( this->iter_, this->res_ );
            }
            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "Jacobi", " stops after " << this->iter_ << " iterations" );
                LOG_INFO ( "Jacobi", " with residual norm " << this->res_ );
                LOG_INFO ( "Jacobi", " ================================================" );
            }

            if ( ( iter_ctrl == IterateControl::kFailureDivergenceTol ) ||
                 ( iter_ctrl == IterateControl::kFailureMaxitsExceeded ) )
                return kSolverExceeded;

            return kSolverSuccess;
        }

        template<class LAD>
        LinearSolverState Jacobi<LAD>::SolveDamped ( const VectorType& b, VectorType* x )
        {
            assert ( x != 0 );
            assert ( this->op_ != 0 );
            assert ( inv_diag_ != 0 );
            assert ( y_ != 0 );
            assert ( w_ > 0.0 );

            this->iter_ = 0;

            // create residual vector
            VectorType res;
            res.CloneFromWithoutContent ( b );

            // compute residual
            this->op_->VectorMult ( *x, &res );
            res.ScaleAdd ( b, static_cast < DataType > ( -1.0 ) );
            this->res_ = res.Norm2 ( );
            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "Jacobi", " ================================================" );
                LOG_INFO ( "Jacobi", " SolveDamped starts with residual norm " << this->res_ );
            }

            IterateControl::State iter_ctrl = this->control ( ).Check ( this->iter_, this->res_ );

            while ( iter_ctrl == IterateControl::kIterate )
            {
                // do iteration
                if ( async_ ) x->UpdateGhost ( );
                for ( int i = 0; i < inner_iter_; ++i )
                {
                    DoDampedIteration ( b, x );
                }

                // compute residual
                this->op_->VectorMult ( *x, &res );
                res.ScaleAdd ( b, static_cast < DataType > ( -1.0 ) );
                this->res_ = res.Norm2 ( );

                // check
                iter_ctrl = this->control ( ).Check ( this->iter_, this->res_ );
            }
            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "Jacobi", " stops after " << this->iter_ << " iterations" );
                LOG_INFO ( "Jacobi", " with residual norm " << this->res_ );
                LOG_INFO ( "Jacobi", " ================================================" );
            }

            if ( ( iter_ctrl == IterateControl::kFailureDivergenceTol ) ||
                 ( iter_ctrl == IterateControl::kFailureMaxitsExceeded ) )
                return kSolverExceeded;

            return kSolverSuccess;
        }

        template<class LAD>
        LinearSolverState Jacobi<LAD>::SmoothNormal ( const VectorType& b, VectorType* x )
        {
            assert ( smooth_iter_ > 0 );

            for ( int n = 0; n < smooth_iter_; ++n )
            {
                if ( async_ ) x->UpdateGhost ( );
                for ( int i = 0; i < inner_iter_; ++i )
                {
                    DoNormalIteration ( b, x );
                }
            }

            return kSolverSuccess;
        }

        template<class LAD>
        LinearSolverState Jacobi<LAD>::SmoothDamped ( const VectorType& b, VectorType* x )
        {
            assert ( smooth_iter_ > 0 );

            for ( int n = 0; n < smooth_iter_; ++n )
            {
                if ( async_ ) x->UpdateGhost ( );
                for ( int i = 0; i < inner_iter_; ++i )
                {
                    DoDampedIteration ( b, x );
                }
            }

            return kSolverSuccess;
        }

        template<class LAD>
        void Jacobi<LAD>::DoNormalIteration ( const VectorType& b, VectorType* x )
        {
            assert ( x != 0 );
            assert ( this->op_ != 0 );
            assert ( inv_diag_ != 0 );
            assert ( y_ != 0 );

            ++( this->iter_ );

            // start ghost update
            if ( !async_ )
            {
                x->SendBorder ( );
                x->ReceiveGhost ( );
            }

            // multiply diagonal block without diagonal entries 
            this->op_->diagonal ( ).VectorMultNoDiag ( x->interior ( ), y_ );

            // multiply offdiagonal block
            if ( !async_ ) x->WaitForRecv ( );
            this->op_->offdiagonal ( ).VectorMultAdd ( x->ghost ( ), y_ );

            y_->ScaleAdd ( -1.0, b.interior ( ) );

            // multiply with inverse diagonal
            y_->ElementWiseMult ( *inv_diag_ );

            x->interior ( ).CopyFrom ( *y_ );

            if ( !async_ ) x->WaitForSend ( );
        }

        template<class LAD>
        void Jacobi<LAD>::DoDampedIteration ( const VectorType& b, VectorType* x )
        {
            assert ( x != 0 );
            assert ( this->op_ != 0 );
            assert ( inv_diag_ != 0 );
            assert ( y_ != 0 );
            assert ( w_ > 0.0 );

            ++( this->iter_ );

            // start ghost update
            if ( !async_ )
            {
                x->SendBorder ( );
                x->ReceiveGhost ( );
            }

            // multiply diagonal block without diagonal entries 
            this->op_->diagonal ( ).VectorMultNoDiag ( x->interior ( ), y_ );

            // multiply offdiagonal block
            if ( !async_ ) x->WaitForRecv ( );
            this->op_->offdiagonal ( ).VectorMultAdd ( x->ghost ( ), y_ );

            y_->ScaleAdd ( -1.0, b.interior ( ) );

            // multiply with inverse diagonal
            y_->ElementWiseMult ( *inv_diag_ );

            // damping
            x->interior ( ).Scale ( 1.0 - w_ );

            // damping
            x->interior ( ).Axpy ( *y_, w_ );

            if ( !async_ ) x->WaitForSend ( );
        }

        /// template instantiation
        template class Jacobi<LADescriptorCoupledD>;
        template class Jacobi<LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow
