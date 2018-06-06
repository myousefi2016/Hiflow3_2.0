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

/// @author Martin Wlotzka

#ifndef HILFOW_LINEARSOLVER_MPIR_H_
#    define HILFOW_LINEARSOLVER_MPIR_H_

#    include "linear_solver/linear_solver.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Mixed precision iterative refinement linear solver.

        template<class LAD_high, class LAD_low>
        class MPIR : public LinearSolver<LAD_high>
        {
          public:

            // high precision matrix, vector and data type
            typedef typename LAD_high::MatrixType MatrixType_high;
            typedef typename LAD_high::VectorType VectorType_high;
            typedef typename LAD_high::DataType DataType_high;

            // low precision matrix, vector and data type
            typedef typename LAD_low::MatrixType MatrixType_low;
            typedef typename LAD_low::VectorType VectorType_low;
            typedef typename LAD_low::DataType DataType_low;

            MPIR ( );

            virtual ~MPIR ( )
            {
            }

            // inherited from LinearSolver base class

            virtual void SetupOperator ( MatrixType_high& op )
            {
                Ah_ = &op;
            }

            void SetHighPrecMatrix ( MatrixType_high& mh )
            {
                Ah_ = &mh;
            }

            void SetLowPrecMatrix ( MatrixType_low& ml )
            {
                Al_ = &ml;
            }

            void SetHighPrecResVector ( VectorType_high& vh )
            {
                resh_ = &vh;
            }

            void SetHighPrecCorVector ( VectorType_high& vh )
            {
                corh_ = &vh;
            }

            void SetLowPrecResVector ( VectorType_low& vl )
            {
                resl_ = &vl;
            }

            void SetLowPrecCorVector ( VectorType_low& vl )
            {
                corl_ = &vl;
            }

            virtual void SetupPreconditioner ( Preconditioner<LAD_high>& precond );

            virtual LinearSolverState Solve ( const VectorType_high& b, VectorType_high* x );

            /// Sets the relative tolerance.
            /// Needed by Inexact Newton Methods
            /// @param reltol relative tolerance of residual to converge

            void SetRelativeTolerance ( double reltol )
            {
                int maxits = this->control_.maxits ( );
                double atol = this->control_.absolute_tol ( );
                double dtol = this->control_.divergence_tol ( );
                this->control_.Init ( maxits, atol, reltol, dtol );
            }

            void SetHighPrecPreSmoother ( LinearSolver<LAD_high>& s )
            {
                pre_smoother_high_ = &s;
            }

            void SetHighPrecPostSmoother ( LinearSolver<LAD_high>& s )
            {
                post_smoother_high_ = &s;
            }

            void SetLowPrecPreSmoother ( LinearSolver<LAD_low>& s )
            {
                pre_smoother_low_ = &s;
            }

            void SetLowPrecPostSmoother ( LinearSolver<LAD_low>& s )
            {
                post_smoother_low_ = &s;
            }

            void SetHighPrecErrCorSolver ( LinearSolver<LAD_high>& ecs )
            {
                error_correction_solver_high_ = &ecs;
            }

            void SetLowPrecErrCorSolver ( LinearSolver<LAD_low>& ecs )
            {
                error_correction_solver_low_ = &ecs;
            }

            void SetMpirParameters ( const int max_ech,
                                     const int max_ecl,
                                     const DataType_high abs_tol,
                                     const DataType_high rel_tol,
                                     const DataType_high div_tol );

            void SetHighPrecPreSmootherParameters ( const int max_iter,
                                                    const DataType_high abs_tol,
                                                    const DataType_high rel_tol,
                                                    const DataType_high div_tol );
            void ResetHighPrecPreSmootherParameters ( );

            void SetHighPrecPostSmootherParameters ( const int max_iter,
                                                     const DataType_high abs_tol,
                                                     const DataType_high rel_tol,
                                                     const DataType_high div_tol );
            void ResetHighPrecPostSmootherParameters ( );

            void SetLowPrecPreSmootherParameters ( const int max_iter,
                                                   const DataType_low abs_tol,
                                                   const DataType_low rel_tol,
                                                   const DataType_low div_tol );
            void ResetLowPrecPreSmootherParameters ( );

            void SetLowPrecPostSmootherParameters ( const int max_iter,
                                                    const DataType_low abs_tol,
                                                    const DataType_low rel_tol,
                                                    const DataType_low div_tol );
            void ResetLowPrecPostSmootherParameters ( );

            void SetHighPrecErrCorSolverParameters ( const int max_iter,
                                                     const DataType_high abs_tol,
                                                     const DataType_high rel_tol,
                                                     const DataType_high div_tol );
            void ResetHighPrecErrCorSolverParameters ( );

            void SetLowPrecErrCorSolverParameters ( const int max_iter,
                                                    const DataType_low abs_tol,
                                                    const DataType_low rel_tol,
                                                    const DataType_low div_tol );
            void ResetLowPrecErrCorSolverParameters ( );

            void UseHighPrecPreSmoother ( const bool use )
            {
                use_pre_smoother_high_ = use;
            }

            void UseHighPrecPostSmoother ( const bool use )
            {
                use_post_smoother_high_ = use;
            }

            void UseLowPrecPreSmoother ( const bool use )
            {
                use_pre_smoother_low_ = use;
            }

            void UseLowPrecPostSmoother ( const bool use )
            {
                use_post_smoother_low_ = use;
            }

          protected:

            bool ContinueIteration ( ) const;

            int max_ech_; // max no. of error corrections in high precision
            int max_ecl_; // max no. of error corrections in low precision

            int n_ech_; // no. of error corrections done in high precision
            int n_ecl_; // no. of error corrections done in low precision

            int max_iter_preh_; // max iter per error correction loop for pre_smoother_high_
            int max_iter_posth_; // max iter per error correction loop for post_smoother_high_
            int max_iter_prel_; // max iter per error correction loop for pre_smoother_low_
            int max_iter_postl_; // max iter per error correction loop for post_smoother_low_
            int max_iter_ech_; // max iter per error correction loop for error_correction_solver_high_
            int max_iter_ecl_; // max iter per error correction loop for error_correction_solver_low_

            int accu_iter_preh_; // accumulated no. of iterations of pre_smoother_high_
            int accu_iter_posth_; // accumulated no. of iterations of post_smoother_high_
            int accu_iter_prel_; // accumulated no. of iterations of pre_smoother_low_
            int accu_iter_postl_; // accumulated no. of iterations of post_smoother_low_
            int accu_iter_ech_; // accumulated no. of iterations of error_correction_solver_high_
            int accu_iter_ecl_; // accumulated no. of iterations of error_correction_solver_low_

            DataType_high abs_tol_mpir_; // overall absolute tolerance
            DataType_high abs_tol_preh_; // absolute tolerance pre_smoother_high_
            DataType_high abs_tol_posth_; // absolute tolerance post_smoother_high_
            DataType_low abs_tol_prel_; // absolute tolerance pre_smoother_low_
            DataType_low abs_tol_postl_; // absolute tolerance post_smoother_low_
            DataType_high abs_tol_ech_; // absolute tolerance error_correction_solver_high_
            DataType_low abs_tol_ecl_; // absolute tolerance error_correction_solver_low_

            DataType_high rel_tol_mpir_; // overall relative tolerance
            DataType_high rel_tol_preh_; // relative tolerance pre_smoother_high_
            DataType_high rel_tol_posth_; // relative tolerance post_smoother_high_
            DataType_low rel_tol_prel_; // relative tolerance pre_smoother_low_
            DataType_low rel_tol_postl_; // relative tolerance post_smoother_low_
            DataType_high rel_tol_ech_; // relative tolerance error_correction_solver_high_
            DataType_low rel_tol_ecl_; // relative tolerance error_correction_solver_low_

            DataType_high div_tol_mpir_; // overall divergence tolerance
            DataType_high div_tol_preh_; // divergence tolerance pre_smoother_high_
            DataType_high div_tol_posth_; // divergence tolerance post_smoother_high_
            DataType_low div_tol_prel_; // divergence tolerance pre_smoother_low_
            DataType_low div_tol_postl_; // divergence tolerance post_smoother_low_
            DataType_high div_tol_ech_; // divergence tolerance error_correction_solver_high_
            DataType_low div_tol_ecl_; // divergence tolerance error_correction_solver_low_

            LinearSolver<LAD_high>* pre_smoother_high_;
            LinearSolver<LAD_high>* post_smoother_high_;

            LinearSolver<LAD_low>* pre_smoother_low_;
            LinearSolver<LAD_low>* post_smoother_low_;

            LinearSolver<LAD_high>* error_correction_solver_high_;
            LinearSolver<LAD_low>* error_correction_solver_low_;

            MatrixType_high* Ah_; // high precision matrix
            MatrixType_low* Al_; // low precision matrix

            VectorType_high *resh_, *corh_; // high precision residual and correction vector
            VectorType_low *resl_, *corl_; // low precision residual and correction vector

            DataType_high res0_; // initial residual norm, only used in high precision

            bool use_pre_smoother_high_;
            bool use_post_smoother_high_;
            bool use_pre_smoother_low_;
            bool use_post_smoother_low_;
        };

    } // namespace la
} // namespace hiflow

#endif
