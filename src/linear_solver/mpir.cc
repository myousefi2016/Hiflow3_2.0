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

#include "mpir.h"
#include "linear_algebra/la_descriptor.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD_high, class LAD_low>
        MPIR<LAD_high, LAD_low>::MPIR ( )
        : LinearSolver<LAD_high>( ),
        max_ech_ ( 0 ), max_ecl_ ( 0 ), n_ech_ ( 0 ), n_ecl_ ( 0 ),
        max_iter_preh_ ( 0 ), max_iter_posth_ ( 0 ),
        max_iter_prel_ ( 0 ), max_iter_postl_ ( 0 ),
        max_iter_ech_ ( 0 ), max_iter_ecl_ ( 0 ),
        accu_iter_preh_ ( 0 ), accu_iter_posth_ ( 0 ),
        accu_iter_prel_ ( 0 ), accu_iter_postl_ ( 0 ),
        accu_iter_ech_ ( 0 ), accu_iter_ecl_ ( 0 ),
        abs_tol_mpir_ ( 0.0 ),
        abs_tol_preh_ ( 0.0 ), abs_tol_posth_ ( 0.0 ),
        abs_tol_prel_ ( 0.0 ), abs_tol_postl_ ( 0.0 ),
        abs_tol_ech_ ( 0.0 ), abs_tol_ecl_ ( 0.0 ),
        rel_tol_mpir_ ( 0.0 ),
        rel_tol_preh_ ( 0.0 ), rel_tol_posth_ ( 0.0 ),
        rel_tol_prel_ ( 0.0 ), rel_tol_postl_ ( 0.0 ),
        rel_tol_ech_ ( 0.0 ), rel_tol_ecl_ ( 0.0 ),
        div_tol_mpir_ ( 0.0 ),
        div_tol_preh_ ( 0.0 ), div_tol_posth_ ( 0.0 ),
        div_tol_prel_ ( 0.0 ), div_tol_postl_ ( 0.0 ),
        div_tol_ech_ ( 0.0 ), div_tol_ecl_ ( 0.0 ),
        pre_smoother_high_ ( 0 ), post_smoother_high_ ( 0 ),
        pre_smoother_low_ ( 0 ), post_smoother_low_ ( 0 ),
        error_correction_solver_high_ ( 0 ),
        error_correction_solver_low_ ( 0 ),
        Ah_ ( 0 ), Al_ ( 0 ), resh_ ( 0 ), corh_ ( 0 ), resl_ ( 0 ), corl_ ( 0 ),
        res0_ ( 0.0 ),
        use_pre_smoother_high_ ( false ),
        use_post_smoother_high_ ( false ),
        use_pre_smoother_low_ ( false ),
        use_post_smoother_low_ ( false )
        {
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Linear solver", "MPIR" );
            }
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::SetupPreconditioner ( Preconditioner<LAD_high>& precond )
        {

        }

        template<class LAD_high, class LAD_low>
        LinearSolverState MPIR<LAD_high, LAD_low>::Solve ( const VectorType_high& b, VectorType_high* x )
        {
            assert ( x != 0 );
            assert ( Ah_ != 0 );
            assert ( Al_ != 0 );
            assert ( resh_ != 0 );
            assert ( corh_ != 0 );
            assert ( resl_ != 0 );
            assert ( corl_ != 0 );

            // reset iteration counters
            this->iter_ = 0;
            accu_iter_preh_ = 0;
            accu_iter_posth_ = 0;
            accu_iter_prel_ = 0;
            accu_iter_postl_ = 0;
            accu_iter_ech_ = 0;
            accu_iter_ecl_ = 0;

            if ( use_pre_smoother_high_ )
            {
                assert ( pre_smoother_high_ != 0 );
                ResetHighPrecPreSmootherParameters ( );
                pre_smoother_high_->SetupOperator ( *Ah_ );
            }

            if ( use_pre_smoother_low_ )
            {
                assert ( pre_smoother_low_ != 0 );
                ResetLowPrecPreSmootherParameters ( );
                pre_smoother_low_->SetupOperator ( *Al_ );
            }

            assert ( error_correction_solver_low_ != 0 );
            error_correction_solver_low_->SetupOperator ( *Al_ );
            ResetLowPrecErrCorSolverParameters ( );

            if ( use_post_smoother_low_ )
            {
                assert ( post_smoother_low_ != 0 );
                ResetLowPrecPostSmootherParameters ( );
                post_smoother_low_->SetupOperator ( *Al_ );
            }

            assert ( error_correction_solver_high_ != 0 );
            ResetHighPrecErrCorSolverParameters ( );
            error_correction_solver_high_->SetupOperator ( *Ah_ );

            if ( use_post_smoother_high_ )
            {
                assert ( post_smoother_high_ != 0 );
                ResetHighPrecPostSmootherParameters ( );
                post_smoother_high_->SetupOperator ( *Ah_ );
            }

            // compute high precision initial residual
            Ah_->VectorMult ( *x, resh_ );
            resh_->ScaleAdd ( b, -1.0 );
            this->res_ = resh_->Norm2 ( );
            res0_ = this->res_;
            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "MPIR", " ===================================================================" );
                LOG_INFO ( "MPIR", " start with residual norm " << this->res_ );
                LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
            }

            while ( ContinueIteration ( ) )
            {
                // high precision pre-smoothing
                if ( use_pre_smoother_high_ )
                {
                    if ( this->print_level_ > 2 )
                    {
                        LOG_INFO ( "MPIR", " high precision pre-smoother" );
                    }
                    pre_smoother_high_->Solve ( b, x );
                    if ( this->print_level_ > 2 )
                    {
                        LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                    }
                    accu_iter_preh_ += pre_smoother_high_->iter ( );
                    this->iter_ += pre_smoother_high_->iter ( );

                    // recompute high precision residual
                    Ah_->VectorMult ( *x, resh_ );
                    resh_->ScaleAdd ( b, -1.0 );
                    this->res_ = resh_->Norm2 ( );
                }

                // error correction
                if ( ( n_ecl_ < max_ecl_ ) &&
                     ( this->res_ > abs_tol_ecl_ ) &&
                     ( this->res_ > res0_ * rel_tol_ecl_ ) &&
                     ( this->res_ < div_tol_ecl_ ) )
                    // do low precision error correction
                {
                    ++n_ecl_;

                    // transfer high precision residual to low precision
                    resl_->CastInteriorFrom ( *resh_ );
                    corl_->Zeros ( );

                    // low precision pre-smoothing
                    if ( use_pre_smoother_low_ )
                    {
                        if ( this->print_level_ > 2 )
                        {
                            LOG_INFO ( "MPIR", " low precision pre-smoother" );
                        }
                        pre_smoother_low_->Solve ( *resl_, corl_ );
                        if ( this->print_level_ > 2 )
                        {
                            LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                        }
                        accu_iter_prel_ += pre_smoother_low_->iter ( );
                        this->iter_ += pre_smoother_low_->iter ( );
                    }

                    // compute correction
                    if ( this->print_level_ > 2 )
                    {
                        LOG_INFO ( "MPIR", " low precision error correction" );
                    }
                    error_correction_solver_low_->Solve ( *resl_, corl_ );
                    if ( this->print_level_ > 2 )
                    {
                        LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                    }
                    accu_iter_ecl_ += error_correction_solver_low_->iter ( );
                    this->iter_ += error_correction_solver_low_->iter ( );

                    // low precision post-smoothing
                    if ( use_post_smoother_low_ )
                    {
                        if ( this->print_level_ > 2 )
                        {
                            LOG_INFO ( "MPIR", " low precision post-smoother" );
                        }
                        post_smoother_low_->Solve ( *resl_, corl_ );
                        if ( this->print_level_ > 2 )
                        {
                            LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                        }
                        accu_iter_postl_ += post_smoother_low_->iter ( );
                        this->iter_ += post_smoother_low_->iter ( );
                    }

                    // transfer low precision correction to high precision
                    DataType_low corl_norm = corl_->Norm2 ( );
                    if ( this->print_level_ > 1 )
                    {
                        LOG_INFO ( "MPIR", "low precision correction norm = " << corl_norm );
                    }
                    corh_->CastInteriorFrom ( *corl_ );
                }
                else if ( ( n_ech_ < max_ech_ ) &&
                          ( this->res_ > abs_tol_ech_ ) &&
                          ( this->res_ > res0_ * rel_tol_ech_ ) &&
                          ( this->res_ < div_tol_ech_ ) )
                    // do high precision error correction
                {
                    ++n_ech_;

                    corh_->Zeros ( );
                    if ( this->print_level_ > 2 )
                    {
                        LOG_INFO ( "MPIR", " high precision error correction" );
                    }
                    error_correction_solver_high_->Solve ( *resh_, corh_ );
                    if ( this->print_level_ > 2 )
                    {
                        LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                    }

                    accu_iter_ech_ += error_correction_solver_high_->iter ( );
                    this->iter_ += error_correction_solver_high_->iter ( );
                }

                DataType_high corh_norm = corh_->Norm2 ( );
                if ( this->print_level_ > 1 )
                {
                    LOG_INFO ( "MPIR", "high precision correction norm = " << corh_norm );
                }
                x->Axpy ( *corh_, 1.0 );

                // high precision post-smoothing
                if ( use_post_smoother_high_ )
                {
                    if ( this->print_level_ > 2 )
                    {
                        LOG_INFO ( "MPIR", " high precision post-smoother" );
                    }
                    post_smoother_high_->Solve ( b, x );
                    if ( this->print_level_ > 2 )
                    {
                        LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                    }
                    accu_iter_posth_ += post_smoother_high_->iter ( );
                    this->iter_ += post_smoother_high_->iter ( );
                }

                // compute high precision residual
                Ah_->VectorMult ( *x, resh_ );
                resh_->ScaleAdd ( b, -1.0 );
                this->res_ = resh_->Norm2 ( );
                res0_ = this->res_;
                if ( this->print_level_ > 1 )
                {
                    LOG_INFO ( "MPIR", " residual norm " << this->res_ << " after error correction" );
                    LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                }
            }

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "MPIR", " stops after " << this->iter_ << " total iterations" );
                LOG_INFO ( "MPIR", " with residual norm " << this->res_ );
                LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                LOG_INFO ( "MPIR", " number of low precision error correction loops: " << n_ecl_ );
                LOG_INFO ( "MPIR", " accumulated low precision pre-smoother iterations: " << accu_iter_prel_ );
                LOG_INFO ( "MPIR", " accumulated low precision error correction solver iterations: " << accu_iter_ecl_ );
                LOG_INFO ( "MPIR", " accumulated low precision post-smoother iterations: " << accu_iter_postl_ );
                LOG_INFO ( "MPIR", " -------------------------------------------------------------------" );
                LOG_INFO ( "MPIR", " number of high precision error correction loops: " << n_ech_ );
                LOG_INFO ( "MPIR", " accumulated high precision pre-smoother iterations: " << accu_iter_preh_ );
                LOG_INFO ( "MPIR", " accumulated high precision error correction solver iterations: " << accu_iter_ech_ );
                LOG_INFO ( "MPIR", " accumulated high precision post-smoother iterations: " << accu_iter_posth_ );
                LOG_INFO ( "MPIR", " ===================================================================" );
            }

            if ( ( this->res_ < abs_tol_mpir_ ) || ( this->res_ < res0_ * rel_tol_mpir_ ) )
                return kSolverSuccess;

            return kSolverExceeded;

        }

        template<class LAD_high, class LAD_low>
        bool MPIR<LAD_high, LAD_low>::ContinueIteration ( ) const
        {
            return (this->res_ > abs_tol_mpir_ ) &&
                    ( this->res_ > res0_ * rel_tol_mpir_ ) &&
                    ( this->res_ < div_tol_mpir_ ) &&
                    (
                    (
                    ( n_ech_ < max_ech_ ) &&
                    ( this->res_ > abs_tol_ech_ ) &&
                    ( this->res_ > res0_ * rel_tol_ech_ ) &&
                    ( this->res_ < div_tol_ech_ )
                    ) ||
                    (
                    ( n_ecl_ < max_ecl_ ) &&
                    ( this->res_ > abs_tol_ecl_ ) &&
                    ( this->res_ > res0_ * rel_tol_ecl_ ) &&
                    ( this->res_ < div_tol_ecl_ )
                    )
                    );
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::SetMpirParameters ( const int max_ech,
                                                          const int max_ecl,
                                                          const DataType_high abs_tol,
                                                          const DataType_high rel_tol,
                                                          const DataType_high div_tol )
        {
            max_ech_ = max_ech;
            max_ecl_ = max_ecl;
            abs_tol_mpir_ = abs_tol;
            rel_tol_mpir_ = rel_tol;
            div_tol_mpir_ = div_tol;
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::SetHighPrecPreSmootherParameters ( const int max_iter,
                                                                         const DataType_high abs_tol,
                                                                         const DataType_high rel_tol,
                                                                         const DataType_high div_tol )
        {
            max_iter_preh_ = max_iter;
            abs_tol_preh_ = abs_tol;
            rel_tol_preh_ = rel_tol;
            div_tol_preh_ = div_tol;
            ResetHighPrecPreSmootherParameters ( );
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::ResetHighPrecPreSmootherParameters ( )
        {
            if ( pre_smoother_high_ != 0 )
            {
                pre_smoother_high_->InitControl ( max_iter_preh_, abs_tol_preh_, rel_tol_preh_, div_tol_preh_ );
            }
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::SetHighPrecPostSmootherParameters ( const int max_iter,
                                                                          const DataType_high abs_tol,
                                                                          const DataType_high rel_tol,
                                                                          const DataType_high div_tol )
        {
            max_iter_posth_ = max_iter;
            abs_tol_posth_ = abs_tol;
            rel_tol_posth_ = rel_tol;
            div_tol_posth_ = div_tol;
            ResetHighPrecPostSmootherParameters ( );
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::ResetHighPrecPostSmootherParameters ( )
        {
            if ( post_smoother_high_ != 0 )
            {
                post_smoother_high_->InitControl ( max_iter_posth_, abs_tol_posth_, rel_tol_posth_, div_tol_posth_ );
            }
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::SetLowPrecPreSmootherParameters ( const int max_iter,
                                                                        const DataType_low abs_tol,
                                                                        const DataType_low rel_tol,
                                                                        const DataType_low div_tol )
        {
            max_iter_prel_ = max_iter;
            abs_tol_prel_ = abs_tol;
            rel_tol_prel_ = rel_tol;
            div_tol_prel_ = div_tol;
            ResetLowPrecPreSmootherParameters ( );
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::ResetLowPrecPreSmootherParameters ( )
        {
            if ( pre_smoother_low_ != 0 )
            {
                pre_smoother_low_->InitControl ( max_iter_prel_, abs_tol_prel_, rel_tol_prel_, div_tol_prel_ );
            }
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::SetLowPrecPostSmootherParameters ( const int max_iter,
                                                                         const DataType_low abs_tol,
                                                                         const DataType_low rel_tol,
                                                                         const DataType_low div_tol )
        {
            max_iter_postl_ = max_iter;
            abs_tol_postl_ = abs_tol;
            rel_tol_postl_ = rel_tol;
            div_tol_postl_ = div_tol;
            ResetLowPrecPostSmootherParameters ( );
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::ResetLowPrecPostSmootherParameters ( )
        {
            if ( post_smoother_low_ != 0 )
            {
                post_smoother_low_->InitControl ( max_iter_postl_, abs_tol_postl_, rel_tol_postl_, div_tol_postl_ );
            }
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::SetHighPrecErrCorSolverParameters ( const int max_iter,
                                                                          const DataType_high abs_tol,
                                                                          const DataType_high rel_tol,
                                                                          const DataType_high div_tol )
        {
            max_iter_ech_ = max_iter;
            abs_tol_ech_ = abs_tol;
            rel_tol_ech_ = rel_tol;
            div_tol_ech_ = div_tol;
            ResetHighPrecErrCorSolverParameters ( );
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::ResetHighPrecErrCorSolverParameters ( )
        {
            if ( error_correction_solver_high_ != 0 )
            {
                error_correction_solver_high_->InitControl ( max_iter_ech_, abs_tol_ech_, rel_tol_ech_, div_tol_ech_ );
            }
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::SetLowPrecErrCorSolverParameters ( const int max_iter,
                                                                         const DataType_low abs_tol,
                                                                         const DataType_low rel_tol,
                                                                         const DataType_low div_tol )
        {
            max_iter_ecl_ = max_iter;
            abs_tol_ecl_ = abs_tol;
            rel_tol_ecl_ = rel_tol;
            div_tol_ecl_ = div_tol;
            ResetLowPrecErrCorSolverParameters ( );
        }

        template<class LAD_high, class LAD_low>
        void MPIR<LAD_high, LAD_low>::ResetLowPrecErrCorSolverParameters ( )
        {
            if ( error_correction_solver_low_ != 0 )
            {
                error_correction_solver_low_->InitControl ( max_iter_ecl_, abs_tol_ecl_, rel_tol_ecl_, div_tol_ecl_ );
            }
        }

        template class MPIR<LADescriptorCoupledD, LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow