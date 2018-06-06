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

#ifndef HIFLOW_QUADRATURE_QUADRATURE_FACTORY_H_
#    define HIFLOW_QUADRATURE_QUADRATURE_FACTORY_H_

#    include <cstdlib>

#    include "common/log.h"
#    include "quadrature/quadrature.h"

namespace hiflow
{

    /// @brief Creator base class for quadratures in HiFlow.
    /// @author Tobias Hahn

    template <class DataType>
    class QuadratureCreator
    {
      public:
        virtual Quadrature<DataType>* params ( const PropertyTree& c ) = 0;
    };

    /// @brief Concrete creator for quadratures in HiFlow.
    /// @author Tobias Hahn

    template <class DataType>
    class DefaultQuadratureCreator : public QuadratureCreator<DataType>
    {
      public:

        Quadrature<DataType>* params ( const PropertyTree& c )
        {
            Quadrature<DataType>* newQuad = new Quadrature<DataType>( );

            if ( c.contains ( "Size" ) )
                newQuad->set_quadrature ( c["Name"].template get<std::string>( ).c_str ( ), c["Size"].template get<int>( ) );
            else if ( c.contains ( "Order" ) )
                newQuad->set_quadrature_by_order ( c["Name"].template get<std::string>( ).c_str ( ), c["Order"].template get<int>( ) );

            const char* cellType = c["CellType"].template get<std::string>( ).c_str ( );

            if ( !strcmp ( cellType, "Hexahedron" ) )
                newQuad->set_cell_type ( mesh::CellType::HEXAHEDRON );

            else if ( !strcmp ( cellType, "Quadrilateral" ) )
                newQuad->set_cell_type ( mesh::CellType::QUADRILATERAL );

            else if ( !strcmp ( cellType, "Triangle" ) )
                newQuad->set_cell_type ( mesh::CellType::TRIANGLE );

            else if ( !strcmp ( cellType, "Pyramid" ) )
                newQuad->set_cell_type ( mesh::CellType::PYRAMID );

            return newQuad;
        }
    };

    /// @brief Factory for quadratures in HiFlow.
    /// @author Tobias Hahn

    template <class DataType>
    class QuadratureFactory
    {
        typedef std::map< std::string, QuadratureCreator<DataType>* > products_t;
        products_t products;
      public:
        /// Register built-in quadratures on construction.

        QuadratureFactory ( )
        {
            this->Register ( "Default", new DefaultQuadratureCreator<DataType>( ) );
        }

        /// Register new product in factory.

        bool Register ( const std::string id, QuadratureCreator<DataType>* cr )
        {
            return products.insert ( typename products_t::value_type ( id, cr ) ).second;
        }

        /// Get new QuadratureCreator object of given type.

        QuadratureCreator<DataType>* Get ( const std::string & id ) const
        {
            typename products_t::const_iterator it = products.find ( id );
            if ( it != products.end ( ) )
                return it->second;
            else
            {
                LOG_ERROR ( "QuadratureFactory::Get: No quadrature of this name registered." );
                return NULL;
            }
        }
    };

} // namespace hiflow

#endif  // HIFLOW_QUADRATURE_QUADRATURE_FACTORY_H_
