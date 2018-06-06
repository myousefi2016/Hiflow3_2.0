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

#ifndef HIFLOW_POLYNOMIALCHAOS_PCGMRES_H_
#    define HIFLOW_POLYNOMIALCHAOS_PCGMRES_H_

/// @brief GMRES solver for Polynomial Chaos Galerkin projected system Ax = b
///
/// Currently support no preconditioning and right preconditioning. Can be used
/// with assembled matrix A or as a matrix free approach

/// @author Michael Schick

#    include <string>
#    include <vector>
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/gmres.h"
#    include "polynomial_chaos/pc_residual_assembler.h"
#    include "linear_algebra/seq_dense_matrix.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class LAD>
        class PCGMRES : public la::GMRES<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// Default constructor
            PCGMRES ( );
            /// Default destructor
            virtual ~PCGMRES ( );

            /// Set up residual assembler for matrix free approach
            virtual void SetupResidualAssembler ( ResidualAssembler<LAD>& assembler );

          protected:
            /// Applies the Givens Rotation
            virtual void ApplyPlaneRotation ( const DataType& cs, const DataType& sn,
                                              DataType* dx, DataType* dy ) const;
            /// Computes the Givens Rotation
            virtual void GeneratePlaneRotation ( const DataType& dx, const DataType& dy,
                                                 DataType* cs, DataType* sn ) const;
            /// Solve Ax=b without preconditioning
            virtual la::LinearSolverState SolveNoPrecond ( const VectorType& b, VectorType* x );
            /// Solve with right preconditioning
            virtual la::LinearSolverState SolveRight ( const VectorType& b, VectorType* x );

            /// Update the iterate with Krylov data
            virtual void UpdateSolution ( const VectorType * const* V,
                                          const la::SeqDenseMatrix<DataType>& H,
                                          const std::vector<DataType>& g,
                                          int k,
                                          VectorType* x ) const;

            /// Matrix free mode or with assembled matrix mode
            std::string solve_type_;
            /// Residual assembler for matrix free approach
            ResidualAssembler<LAD>* res_assembler_;
        };

    }
}

#endif  // HIFLOW_POLYNOMIALCHAOS_PCGMRES_H_
