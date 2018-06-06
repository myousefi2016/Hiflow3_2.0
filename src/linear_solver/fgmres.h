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

/// @author Hendryk Bockelmann, Chandramowli Subramanian

#ifndef HIFLOW_LINEARSOLVER_FGMRES_H_
#    define HIFLOW_LINEARSOLVER_FGMRES_H_

#    include <string>
#    include <vector>
#    include "linear_solver/gmres.h"
#    include "linear_solver/linear_solver.h"

namespace hiflow
{
    namespace la
    {

        template<class DataType> class SeqDenseMatrix;

        /// @brief Flexible GMRES solver
        ///
        /// Flexible GMRES solver for linear systems Ax=b.

        template<class LAD>
        class FGMRES : public GMRES<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            FGMRES ( );
            virtual ~FGMRES ( );

          private:
            LinearSolverState SolveLeft ( const VectorType& b, VectorType* x );
            LinearSolverState SolveRight ( const VectorType& b, VectorType* x );
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_FGMRES_H_
