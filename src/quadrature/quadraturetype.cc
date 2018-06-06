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

#include "quadraturetype.h"
#include <cassert>

namespace hiflow
{

    template<class DataType>
    QuadratureType<DataType>::QuadratureType ( )
    {
    }

    template<class DataType>
    const std::vector<DataType>* QuadratureType<DataType>::x ( int size ) const
    {
        return &x_[find_quadrature ( size )];
    }

    template<class DataType>
    const std::vector<DataType>* QuadratureType<DataType>::y ( int size ) const
    {
        return &y_[find_quadrature ( size )];
    }

    template<class DataType>
    const std::vector<DataType>* QuadratureType<DataType>::z ( int size ) const
    {
        return &z_[find_quadrature ( size )];
    }

    template<class DataType>
    const std::vector<DataType>* QuadratureType<DataType>::w ( int size ) const
    {
        return &w_[find_quadrature ( size )];
    }

    template<class DataType>
    int QuadratureType<DataType>::size_for_order ( int order ) const
    {
        if ( order < 0 )
        {
            throw "There are no quadrature rules with negative orders!\n";
        }
        else if ( order >= order_size_map_.size ( ) )
        {
            throw "Requested order of quadrature rule not available.\n";
        }
        return order_size_map_[order];
    }

    template<class DataType>
    int QuadratureType<DataType>::find_quadrature ( int size ) const
    {
        std::map<int, int>::const_iterator it = index_field_.find ( size );
        if ( it != index_field_.end ( ) )
            return it->second;
        else
            assert ( 0 );

        return -1;
    }

    // template instanciation
    template class QuadratureType<double>;
    template class QuadratureType<float>;

} // namespace hiflow
