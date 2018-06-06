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

/// @author Tobias Hahn, Simon Gawlok

#ifndef HIFLOW_LINEARSOLVER_LINEAR_SOLVER_FACTORY_H_
#    define HIFLOW_LINEARSOLVER_LINEAR_SOLVER_FACTORY_H_

#    include <cstdlib>
#    include "common/log.h"
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/linear_solver_creator.h"
#    include "linear_solver/cg.h"
#    include "linear_solver/gmres.h"
#    include "linear_solver/bicgstab.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Factory for linear solvers in HiFlow.

        template<class LAD>
        class LinearSolverFactory
        {
            typedef std::map< std::string, LinearSolverCreator<LAD>* > products_t;
            products_t products;
          public:
            /// Register built-in linear solvers on construction.

            LinearSolverFactory ( )
            {
                this->Register ( "CG", new CGcreator<LAD>( ) );
                this->Register ( "GMRES", new GMREScreator<LAD>( ) );
                this->Register ( "BiCGSTAB", new BiCGSTABcreator<LAD>( ) );
                //this->Register( "MUMPS", new MUMPScreator<LAD>() );
            }

            /// Register new product in factory.

            bool Register ( const std::string id, LinearSolverCreator<LAD>* cr )
            {
                return products.insert ( typename products_t::value_type ( id, cr ) ).second;
            }

            /// Get new LinearSolverCreator object of given type.

            LinearSolverCreator<LAD>* Get ( const std::string & id ) const
            {
                typename products_t::const_iterator it = products.find ( id );
                if ( it != products.end ( ) )
                    return it->second;
                else
                {
                    LOG_ERROR ( "LinearSolverFactory::Get: No solver of this name registered." );
                    return NULL;
                }
            }
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_LINEAR_SOLVER_FACTORY_H_
