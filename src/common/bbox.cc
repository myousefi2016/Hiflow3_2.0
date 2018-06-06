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

#include "bbox.h"

#include <iostream>
#include <limits>
#include <cmath>
#include <assert.h>

namespace hiflow
{

    template<class DataType>
    BBox<DataType>::BBox ( int dim )
    {
        extents_.resize ( 2 * dim );
        for ( int i = 0; i < dim; ++i )
        {
            extents_[i * 2] = std::numeric_limits<DataType>::max ( );
            extents_[i * 2 + 1] = -std::numeric_limits<DataType>::max ( );
        }
    }

    template<class DataType>
    BBox<DataType>::BBox ( DataType* extents, int dim )
    {
        extents_.assign ( &extents[0], &extents[2 * dim] );
        assert ( extents_.size ( ) == 2 * dim );
    }

    template<class DataType>
    BBox<DataType>::BBox ( const std::vector<DataType>& extents )
    : extents_ ( extents )
    {
        assert ( !extents.empty ( ) );
        assert ( !( extents.size ( ) % 2 ) );
    }

    template<class DataType>
    DataType BBox<DataType>::min ( int dir ) const
    {
        assert ( extents_.size ( ) > 2 * dir + 1 );
        return extents_[2 * dir];
    }

    template<class DataType>
    DataType BBox<DataType>::max ( int dir ) const
    {
        assert ( extents_.size ( ) > 2 * dir + 1 );
        return extents_[2 * dir + 1];
    }

    template<class DataType>
    void BBox<DataType>::add_point ( const std::vector<DataType>& pt )
    {
        assert ( extents_.size ( ) == 2 * pt.size ( ) );
        for ( int i = 0; i < extents_.size ( ) / 2; ++i )
        {
            extents_[2 * i] = std::min ( extents_[2 * i], pt[i] );
            extents_[2 * i + 1] = std::max ( extents_[2 * i + 1], pt[i] );
        }
    }

    template<class DataType>
    void BBox<DataType>::add_points ( const std::vector<DataType>& pts )
    {
        int dim = extents_.size ( ) / 2;
        assert ( !( pts.size ( ) % dim ) );
        int num_pts = pts.size ( ) / dim;
        for ( int n = 0; n < num_pts; ++n )
        {
            for ( int i = 0; i < dim; ++i )
            {
                extents_[2 * i] = std::min ( extents_[2 * i], pts[n * dim + i] );
                extents_[2 * i + 1] = std::max ( extents_[2 * i + 1], pts[n * dim + i] );
            }
        }
    }

    template<class DataType>
    void BBox<DataType>::uniform_extension ( DataType extension )
    {
        for ( int i = 0; i < extents_.size ( ) / 2; ++i )
        {
            extents_[2 * i] -= extension;
            extents_[2 * i + 1] += extension;
        }
    }

    template<class DataType>
    bool BBox<DataType>::intersects ( const BBox<DataType>& other ) const
    {
        assert ( extents_.size ( ) == 2 * other.get_dim ( ) );
        bool intersection = true;
        for ( int i = 0; i < extents_.size ( ) / 2; ++i )
        {
            if ( ( extents_[2 * i] > other.max ( i ) ) || extents_[2 * i + 1] < other.min ( i ) )
            {
                return false;
            }
        }
        return intersection;
    }

    template<class DataType>
    void BBox<DataType>::print ( std::ostream& os ) const
    {
        os << "[";
        for ( int i = 0; i < extents_.size ( ) / 2; ++i )
        {
            os << extents_[2 * i]
                    << " .. "
                    << extents_[2 * i + 1];
            if ( i < extents_.size ( ) / 2 - 1 )
            {
                os << " x ";
            }
        }
        os << "]\n";
    }

    template<class DataType>
    int BBox<DataType>::get_dim ( ) const
    {
        return extents_.size ( ) / 2;
    }

    template<class DataType>
    std::vector<DataType> BBox<DataType>::get_extents ( ) const
    {
        return extents_;
    }

    template<class DataType>
    std::vector<DataType> BBox<DataType>::get_vertices ( ) const
    {
        int dim = extents_.size ( ) / 2;
        int num_vertices = 1;
        for ( int d = 0; d < dim; ++d )
        {
            num_vertices *= 2;
        }
        std::vector<DataType> vertices ( num_vertices * dim );
        for ( int n = 0; n < num_vertices; ++n )
        {
            int fac = 1;
            for ( int d = 0; d < dim; ++d )
            {
                vertices[n * dim + d] = extents_[2 * d + ( n / fac ) % 2];
                fac *= 2;
            }
        }
        return vertices;
    }

    template<class DataType>
    DataType BBox<DataType>::compute_diagonal ( ) const
    {
        DataType diagonal = 0;
        for ( int i = 0; i < extents_.size ( ) / 2; ++i )
        {
            diagonal += ( extents_[2 * i + 1] - extents_[2 * i] )*( extents_[2 * i + 1] - extents_[2 * i] );
        }
        diagonal = std::sqrt ( diagonal );

        return diagonal;
    }

    template class BBox<double>;
    template class BBox<float>;
}

