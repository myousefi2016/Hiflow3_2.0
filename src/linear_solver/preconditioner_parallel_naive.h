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

/// \author Simon Gawlok

#ifndef HIFLOW_LINEARSOLVER_PRECONDITIONER_PARALLEL_NAIVE_H_
#    define HIFLOW_LINEARSOLVER_PRECONDITIONER_PARALLEL_NAIVE_H_

#    include <vector>

#    include "config.h"
#    include "linear_solver/preconditioner_bjacobi.h"

namespace hiflow
{
    namespace la
    {

        /// \author Simon Gawlok
        /// \brief Parallel naive preconditioner interface. Uses splitting of 
        /// matrix A = D + OD with D diagonal and OD offdiagonal part of A. 
        /// Preconditioning is done via the Jacobi iteration x = D^{-1}(b - OD*x) 
        /// where D^{-1} is approximated by the application of any block-preconditioner.
        ///

        template<class LAD>
        class PreconditionerParallelNaive : public PreconditionerBlockJacobi<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// standard constructor
            PreconditionerParallelNaive ( );
            /// destructor
            virtual ~PreconditionerParallelNaive ( );

            /// Inits parameters for parallel naive preconditioner.
            /// \param num_iter NUmber of Jacobi iterations
            /// \param omega Relaxation parameter for Jacobi iteration.
            void InitParameter ( int num_iter, DataType omega = 1. );

            /// Inits the operator i.e. sets up the local matrix.
            /// @param op linear operator to be preconditioned
            void SetupOperator ( OperatorType& op );

            /// Set preconditioner for diagonal part of matrix
            void SetPreconditioner ( PreconditionerBlockJacobi<LAD>& precond_diag );

            /// Build the precond 
            void Build ( );

            /// Applies the parallel naive preconditioner.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if preconditioning succeeded
            LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );

            /// Clears allocated data.
            void Clear ( );

          protected:

            PreconditionerBlockJacobi<LAD>* precond_diag_; // preconditioner for diagonal block

            DataType omega_; // relaxation parameter

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PRECONDITIONER_PARALLEL_NAIVE_H_
