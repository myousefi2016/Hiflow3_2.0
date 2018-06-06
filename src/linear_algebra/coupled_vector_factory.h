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

//
/// \author Carmen Straub

#ifndef HIFLOW_LINEARALGEBRA_COUPLED_VECTOR_FACTORY_H_
#    define HIFLOW_LINEARALGEBRA_COUPLED_VECTOR_FACTORY_H_

#    include <cstdlib>
#    include "common/log.h"
#    include "linear_algebra/coupled_vector.h"
#    include "common/property_tree.h"

namespace hiflow
{
    namespace la
    {

        template <class DataType> class CoupledVectorCreator;

        /// @brief Factory for coupled vectors in HiFlow.

        template<class DataType>
        class CoupledVectorFactory
        {
            typedef std::map< std::string, CoupledVectorCreator<DataType>* > products_t;
            products_t products;
          public:
            /// Register built-in coupled vectors on construction.

            CoupledVectorFactory ( )
            {
                this->Register ( "CoupledVector", new CoupledVectorCreator<DataType>( ) );
            }

            /// Register new product in factory.

            bool Register ( const std::string id, CoupledVectorCreator<DataType>* cr )
            {
                return products.insert ( typename products_t::value_type ( id, cr ) ).second;
            }

            /// Get new CoupledVectorCreator object of given type.

            CoupledVectorCreator<DataType>* Get ( const std::string & id ) const
            {
                typename products_t::const_iterator it = products.find ( id );
                if ( it != products.end ( ) )
                    return it->second;
                else
                {
                    LOG_ERROR ( "CoupledVectorFactory::Get: No format of this name registered." );
                    return NULL;
                }
            }
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_COUPLED_VECTOR_FACTORY_H_
