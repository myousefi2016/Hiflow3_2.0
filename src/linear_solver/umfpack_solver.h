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

#ifndef HIFLOW_LINEARSOLVER_UMFPACK_SOLVER_H_
#    define HIFLOW_LINEARSOLVER_UMFPACK_SOLVER_H_

#    include "linear_algebra/la_descriptor.h"
#    include "linear_solver/linear_solver.h"
#    include "umfpack.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Wrapper for UMFPACK solver.
        ///
        /// Detailed description.

        template<class LAD>
        class UmfpackSolver : public LinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            UmfpackSolver ( );
            virtual ~UmfpackSolver ( );

            // no use of iterate control

            void InitControl ( int maxits, double abstol, double reltol, double divtol )
            {
            }

            /// Sets the relative tolerance.
            /// Needed by Inexact Newton Methods
            /// @param reltol relative tolerance of residual to converge

            void SetRelativeTolerance ( double reltol )
            {
            }

            void InitParameter ( );
            void SetupOperator ( OperatorType& op );

            void FactorizeSymbolic ( );
            void FactorizeNumeric ( );
            LinearSolverState SolveFactorized ( const VectorType& b, VectorType* x );
            LinearSolverState SolveFactorized ( const std::vector<DataType>& b, std::vector<DataType>* x );

            LinearSolverState Solve ( const VectorType& b, VectorType* x );
            LinearSolverState Solve ( const std::vector<DataType>& b, std::vector<DataType>* x );

            virtual void Build ( );

            void Clear ( );

          private:

            // local matrix in coordinate format
            std::vector<int> row_; // local row indices
            std::vector<int> col_; // local column indices
            std::vector<DataType> val_; // local values

            /// size of the matrix
            int umf_n_;
            /// number of nonzeros
            int umf_nnz_;
            /// CRS storage
            int* Ap_, *Ai_;
            bool Ap_b_, Ai_b_;
            DataType* Ax_;
            bool Ax_b_;

            DataType* Control;
            DataType* Info;

            bool is_symb_factorized_;
            bool is_num_factorized_;

            void* symbolic_, *numeric_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_UMFPACK_SOLVER_H_
