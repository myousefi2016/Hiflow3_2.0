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

/// @author Tobias Hahn

#ifndef HIFLOW_LINEARSOLVER_LINEAR_SOLVER_CREATOR_H_
#    define HIFLOW_LINEARSOLVER_LINEAR_SOLVER_CREATOR_H_

#    include "linear_solver/linear_solver.h"
#    include "common/property_tree.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Creator base class for linear solvers in HiFlow.

        template<class LAD>
        class LinearSolverCreator
        {
          public:
            virtual LinearSolver<LAD>* params ( const PropertyTree& c ) = 0;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_LINEAR_SOLVER_CREATOR_H_
