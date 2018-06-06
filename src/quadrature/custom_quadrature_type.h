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

#ifndef HIFLOW_QUADRATURE_CUSTOM_QUADRATURETYPE_H
#    define HIFLOW_QUADRATURE_CUSTOM_QUADRATURETYPE_H

#    include <vector>
#    include <map>
#    include <iostream>

#    include <boost/function.hpp>

#    include "quadraturetype.h"

namespace hiflow
{

    ///
    /// \class CustomQuadratureType
    /// \brief QuadratureType where user can set the points and weights arbitrarily.
    /// \author Staffan Ronnas
    ///

    template<class DataType>
    class CustomQuadratureType : public QuadratureType<DataType>
    {
      public:

        /// Constructor

        CustomQuadratureType ( const std::vector<DataType>& xc, const std::vector<DataType>& yc,
                               const std::vector<DataType>& zc, const std::vector<DataType>& wgt )
        {
            this->x_.push_back ( xc );
            this->y_.push_back ( yc );
            this->z_.push_back ( zc );
            this->w_.push_back ( wgt );
            this->index_field_.insert ( std::make_pair ( wgt.size ( ), 0 ) );

            // order_size_map left empty
        }

        virtual ~CustomQuadratureType ( )
        {
        }

        /// Polymorphic cloning function, returning a deep copy of itself.

        virtual QuadratureType<DataType>* clone ( ) const
        {
            return new CustomQuadratureType ( *this );
        }
    };

} // namespace hiflow

#endif
