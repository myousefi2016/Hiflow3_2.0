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

#ifndef HIFLOW_LINEARSOLVER_JACOBI_
#    define HIFLOW_LINEARSOLVER_JACOBI_

#    include "linear_solver/linear_solver.h"
#    include "linear_algebra/lmp/lvector.h"
#    include "linear_algebra/lmp/lmatrix_csr.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        class Jacobi : public LinearSolver<LAD>
        {
          public:

            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            Jacobi ( );

            virtual ~Jacobi ( );

            void Prepare ( OperatorType& op, const VectorType& sample );

            void SetSolveMode ( const int mode );
            void SetNumIter ( const int iter );
            void SetInnerIter ( const int iter );
            void SetDampingParameter ( const DataType w );

            void SetAsyncFlag ( const bool async )
            {
                async_ = async;
            }

            virtual LinearSolverState Solve ( const VectorType& b, VectorType* x );

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

            LinearSolverState SolveNormal ( const VectorType& b, VectorType* x );
            LinearSolverState SolveDamped ( const VectorType& b, VectorType* x );
            LinearSolverState SmoothNormal ( const VectorType& b, VectorType* x );
            LinearSolverState SmoothDamped ( const VectorType& b, VectorType* x );

            void DoNormalIteration ( const VectorType& b, VectorType* x );
            void DoDampedIteration ( const VectorType& b, VectorType* x );

          protected:

            int solve_mode_;
            int smooth_iter_;
            int inner_iter_;

            DataType w_;

            bool async_;

            lVector<DataType> *inv_diag_;
            CSR_lMatrix<DataType> *csr_mat_;
            lVector<DataType> *y_;

        };

    } // namespace la
} // namespace hiflow

#endif
