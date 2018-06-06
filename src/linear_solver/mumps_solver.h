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

/// @author Chandramowli Subramanian, Martin Wlotzka

#ifndef HIFLOW_LINEARSOLVER_MUMPS_SOLVER_H_
#    define HIFLOW_LINEARSOLVER_MUMPS_SOLVER_H_

#    include <cassert>
#    include <vector>

#    include <mpi.h>

#    include "linear_solver/linear_solver.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Wrapper for MUMPS solver.
        ///
        /// Detailed description.

        template<class LAD, class MUMPS_STRUCTURE>
        class MumpsSolver : public LinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            MumpsSolver ( );
            virtual ~MumpsSolver ( );

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

            void InitParameter ( const MPI_Comm& comm /* more params */ );
            void SetupOperator ( OperatorType& op );

            void Build ( );

            void FactorizeSymbolic ( );
            void FactorizeNumeric ( );
            LinearSolverState SolveFactorized ( const VectorType& b, VectorType* x );
            LinearSolverState Solve ( const VectorType& b, VectorType* x );

            void Clear ( );

            bool transposed ( ) const
            {
                return this->transposed_;
            }

            bool& transposed ( )
            {
                return this->transposed_;
            }

            /// @return MPI communicator

            const MPI_Comm& comm ( ) const
            {
                return this->comm_;
            }

          private:
            void CreateLocalMatrix ( );

            MPI_Comm comm_;

            MUMPS_STRUCTURE mumps_struc_; // MUMPS structure for single or double precision
            bool mumps_struc_initialized_;
            bool transposed_; // indicates if transposed system is to be solved

            // local matrix in coordinate format
            std::vector<int> row_; // local row indices
            std::vector<int> col_; // local column indices
            std::vector<DataType> val_; // local values
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_MUMPS_SOLVER_H_
