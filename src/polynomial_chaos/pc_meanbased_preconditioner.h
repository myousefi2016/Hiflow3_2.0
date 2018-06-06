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

#ifndef HIFLOW_POLYNOMIALCHAOS_MEANBASED_PRECONDITIONER_H_
#    define HIFLOW_POLYNOMIALCHAOS_MEANBASED_PRECONDITIONER_H_

/// \file pc_meanbased_preconditioner.h
/// \brief Preconditioner (mean based) for Polynomial Chaos.
/// Currently, it only support LU factorization of mean matrix
/// provided by the external library UMFPACK
/// \author Michael Schick

#    include <vector>

#    include "config.h"
#    include "linear_solver/preconditioner.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_algebra/coupled_matrix.h"
#    include "linear_solver/umfpack_solver.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class LAD>
        class MeanbasedPreconditioner : public la::Preconditioner<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// Default constructor
            MeanbasedPreconditioner ( );
            /// Destructor
            virtual ~MeanbasedPreconditioner ( );
            /// Build the preconditioner
            virtual void Build ( );

            /// Applies the mean based preconditioner.
            virtual la::LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );

          private:
            /// Blocksolver Umfpack
            la::UmfpackSolver<la::LADescriptorCoupledD> linear_solver_;
            /// Copy of operator matrix
            la::CoupledMatrix<DataType> matrix_;

        };

    } // namespace stochastic
} // namespace hiflow

#endif  // HIFLOW_POLYNOMIALCHAOS_MEANBASED_PRECONDITIONER_H_
