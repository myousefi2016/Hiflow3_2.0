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

#ifndef HIFLOW_TOOLS_POINT_EVALUATOR
#    define HIFLOW_TOOLS_POINT_EVALUATOR

/// \author Jonas Kratzke, Jonathan Schwegler, Simon Gawlok, Philipp Gerstner
#    include "space/vector_space.h"
#    include "mesh/mesh.h"
#    include "space/visualization_tools.h"
#    include "fem/cell_transformation.h"
#    include "mesh/geometric_search.h"

namespace hiflow
{

    /// \brief Evaluates a Function in _one_ point with help of
    /// the GridGeometricSearch.

    template<class DataType>
    class PointEvaluator
    {
        // Type of function for evaluation.
        typedef boost::function3<void,
        const mesh::Entity&, // cell
        const std::vector<Coordinate>&, // reference coordinates
        std::vector<DataType>& // values of function at the points
        > EvalFunction;

      public:
        PointEvaluator ( const VectorSpace<DataType>& space );

        /// Evaluates a Function. Returns false if the point was not found.
        /// Furthermore, the cells containing the requested point are returned.
        bool evaluate_fun ( const EvalFunction& fun,
                            const std::vector<Coordinate>& point,
                            DataType& value,
                            std::vector<int>& cells ) const;

        bool evaluate_fun ( const EvalFunction& fun,
                            const std::vector<Coordinate>& point,
                            DataType& value ) const;

        /// Evaluates a Function and communicates the value to all processes.
        /// Returns false if point was not found (and value = 0).
        bool evaluate_fun_global ( const EvalFunction& fun,
                                   const std::vector<Coordinate>& point,
                                   DataType& value,
                                   std::vector<int>& cells,
                                   const MPI_Comm& comm ) const;

        bool evaluate_fun_global ( const EvalFunction& fun,
                                   const std::vector<Coordinate>& point,
                                   DataType& value,
                                   const MPI_Comm& comm ) const;

        /// Set trial cells which are searched for the given point before the complete grid is searched.
        /// Caution: Using trial cells might give a different return value if all of the following conditions are satisfied:
        /// (i) the 'point' lies on the interface of several cells
        /// (ii) not all of these cells are contained in the list of trial cells
        /// (iii) the function 'fun' is discontinuous in the given 'point'
        void set_trial_cells ( const std::vector<int>& trial_cells );

      private:
        const VectorSpace<DataType>& space_;
        std::vector<int> trial_cells_;
    };
}

#endif
