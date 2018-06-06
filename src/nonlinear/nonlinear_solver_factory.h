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

#ifndef HIFLOW_NONLINEAR_NONLINEAR_SOLVER_FACTORY_H_
#    define HIFLOW_NONLINEAR_NONLINEAR_SOLVER_FACTORY_H_

#    include <cstdlib>

#    include "common/log.h"
#    include "nonlinear/nonlinear_solver.h"
#    include "nonlinear/nonlinear_solver_creator.h"
#    include "nonlinear/newton.h"

namespace hiflow
{

    /// @brief Factory for nonlinear solvers in HiFlow.
    /// @author Tobias Hahn

    template<class LAD>
    class NonlinearSolverFactory
    {
        typedef std::map< std::string, NonlinearSolverCreator<LAD>* > products_t;
        products_t products;
      public:
        /// Register built-in nonlinear solvers on construction.

        NonlinearSolverFactory ( )
        {
            this->Register ( "Newton", new Newtoncreator<LAD>( ) );
        }

        /// Register new product in factory.

        bool Register ( const std::string id, NonlinearSolverCreator<LAD>* cr )
        {
            return products.insert ( typename products_t::value_type ( id, cr ) ).second;
        }

        /// Get new NonlinearSolverCreator object of given type.

        NonlinearSolverCreator<LAD>* Get ( const std::string & id ) const
        {
            typename products_t::const_iterator it = products.find ( id );
            if ( it != products.end ( ) )
                return it->second;
            else
            {
                LOG_ERROR ( "NonlinearSolverFactory::Get: No solver of this name registered." );
                return NULL;
            }
        }
    };

} // namespace hiflow

#endif  // HIFLOW_NONLINEAR_NONLINEAR_SOLVER_FACTORY_H_
