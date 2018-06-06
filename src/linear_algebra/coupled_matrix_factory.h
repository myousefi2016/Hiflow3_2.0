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

/// \author Carmen Straub

#ifndef HIFLOW_LINEARALGEBRA_COUPLED_MATRIX_FACTORY_H_
#    define HIFLOW_LINEARALGEBRA_COUPLED_MATRIX_FACTORY_H_

#    include <cstdlib>
#    include "common/log.h"
#    include "linear_algebra/coupled_matrix.h"
#    include "common/property_tree.h"

namespace hiflow
{
    namespace la
    {

        template <class DataType> class CoupledMatrixCreator;

        /// @brief Factory for coupled matrices in HiFlow.

        template<class DataType>
        class CoupledMatrixFactory
        {
            typedef std::map< std::string, CoupledMatrixCreator<DataType>* > products_t;
            products_t products;
          public:
            /// Register built-in coupled matrices on construction.

            CoupledMatrixFactory ( )
            {
                this->Register ( "CoupledMatrix", new CoupledMatrixCreator<DataType>( ) );
            }

            /// Register new product in factory.

            bool Register ( const std::string id, CoupledMatrixCreator<DataType>* cr )
            {
                return products.insert ( typename products_t::value_type ( id, cr ) ).second;
            }

            /// Get new CoupledMatrixCreator object of given type.

            CoupledMatrixCreator<DataType>* Get ( const std::string & id ) const
            {
                typename products_t::const_iterator it = products.find ( id );
                if ( it != products.end ( ) )
                    return it->second;
                else
                {
                    LOG_ERROR ( "CoupledMatrixFactory::Get: No format of this name registered." );
                    return NULL;
                }
            }
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_COUPLED_MATRIX_FACTORY_H_
