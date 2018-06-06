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

#include "bsphere.h"

#include <iostream>
#include <limits>
#include <cmath>
#include <assert.h>

namespace hiflow
{

    template<class DataType>
    BSphere<DataType>::BSphere ( const std::vector<DataType>& origin, DataType radius )
    : origin_ ( origin ), radius_ ( radius )
    {
        assert ( !origin.empty ( ) );
    }

    template<class DataType>
    BSphere<DataType>::BSphere ( int dim, const std::vector<DataType>& pts )
    {
        assert ( !( pts.size ( ) % dim ) );
        int num_pts = pts.size ( ) / dim;
        // Jack Ritter's algorithm [1990] to approximate the smallest bounding sphere

        // pick some point
        std::vector<DataType> point_a ( pts.begin ( ), pts.begin ( ) + dim );
        // determine the point with biggest distance
        std::vector<DataType> point_b;
        DataType dist = 0;
        for ( int n = 1; n < num_pts; ++n )
        {
            DataType dist_temp = 0;
            for ( int d = 0; d < dim; ++d )
            {
                dist_temp += ( point_a[d] - pts[n * dim + d] )*( point_a[d] - pts[n * dim + d] );
            }
            dist_temp = std::sqrt ( dist_temp );
            if ( dist_temp > dist )
            {
                point_b.assign ( pts.begin ( ) + n*dim, pts.begin ( ) + ( n + 1 ) * dim );
                dist = dist_temp;
            }
        }
        // determine the point with biggest distance to that point
        for ( int n = 0; n < num_pts; ++n )
        {
            DataType dist_temp = 0;
            for ( int d = 0; d < dim; ++d )
            {
                dist_temp += ( point_b[d] - pts[n * dim + d] )*( point_b[d] - pts[n * dim + d] );
            }
            dist_temp = std::sqrt ( dist_temp );
            if ( dist_temp > dist )
            {
                point_a.assign ( pts.begin ( ) + n*dim, pts.begin ( ) + ( n + 1 ) * dim );
                dist = dist_temp;
            }
        }
        // now we have determined the two points of the set with biggest distance
        // initialize the bounding sphere with the midpoint as origin
        origin_ = point_a;
        for ( int d = 0; d < dim; ++d )
        {
            origin_[d] = 0.5 * ( origin_[d] + point_b[d] );
        }
        // and the radius is half the value of the distance
        radius_ = 0.5 * dist;
        // all other points can just be added with non-fixed origin
        add_points ( pts, false );
    }

    template<class DataType>
    void BSphere<DataType>::add_points ( const std::vector<DataType>& pts, bool fixed_origin )
    {
        int dim = origin_.size ( );
        assert ( !( pts.size ( ) % dim ) );
        int num_pts = pts.size ( ) / dim;

        for ( int n = 0; n < num_pts; ++n )
        {
            DataType dist = 0;
            for ( int d = 0; d < dim; ++d )
            {
                dist += ( origin_[d] - pts[n * dim + d] )*( origin_[d] - pts[n * dim + d] );
            }
            if ( fixed_origin )
            {
                radius_ = std::max ( radius_, std::sqrt ( dist ) );
            }
            else
            {
                if ( dist > radius_ )
                {
                    radius_ = 0.5 * ( dist + radius_ );
                    DataType factor = radius_ / dist;
                    for ( int d = 0; d < dim; ++d )
                    {
                        origin_[d] = factor * origin_[d] + ( 1 - factor ) * pts[n * dim + d];
                    }
                }
            }
        }
    }

    template<class DataType>
    void BSphere<DataType>::radial_extension ( DataType extension )
    {
        radius_ += extension;
        assert ( radius_ >= 0. );
    }

    template<class DataType>
    bool BSphere<DataType>::intersects ( const BSphere<DataType>& other ) const
    {
        assert ( origin_.size ( ) == other.get_dim ( ) );

        std::vector<DataType> other_origin = other.get_origin ( );
        DataType dist = 0;
        for ( int d = 0; d < origin_.size ( ); ++d )
        {
            dist += ( origin_[d] - other_origin[d] )*( origin_[d] - other_origin[d] );
        }
        dist = std::sqrt ( dist );

        return (dist < radius_ + other.get_radius ( ) );
    }

    template<class DataType>
    bool BSphere<DataType>::intersects ( const BBox<DataType>& box ) const
    {
        assert ( origin_.size ( ) == box.get_dim ( ) );
        int dim = origin_.size ( );

        std::vector<DataType> nearest_box_point ( dim );
        // do a projection on the box
        for ( int d = 0; d < dim; d++ )
        {
            if ( origin_[d] <= box.min ( d ) )
            {
                nearest_box_point[d] = box.min ( d );
            }
            else if ( origin_[d] >= box.max ( d ) )
            {
                nearest_box_point[d] = box.max ( d );
            }
            else
            {
                nearest_box_point[d] = origin_[d];
            }
        }
        DataType dist = 0;
        for ( int d = 0; d < dim; ++d )
        {
            dist += ( origin_[d] - nearest_box_point[d] )*( origin_[d] - nearest_box_point[d] );
        }
        dist = std::sqrt ( dist );

        return (dist <= radius_ );
    }

    template<class DataType>
    bool BSphere<DataType>::contains ( const BBox<DataType>& box ) const
    {
        assert ( origin_.size ( ) == box.get_dim ( ) );

        bool contains = true;
        std::vector<DataType> box_vertices = box.get_vertices ( );
        int dim = origin_.size ( );
        assert ( !( box_vertices.size ( ) % dim ) );
        for ( int n = 0; n < box_vertices.size ( ) / dim; ++n )
        {
            DataType dist = 0;
            for ( int d = 0; d < dim; ++d )
            {
                dist += ( origin_[d] - box_vertices[n * dim + d] )*( origin_[d] - box_vertices[n * dim + d] );
            }
            if ( std::sqrt ( dist ) > radius_ )
            {
                return false;
            }
        }
        return contains;
    }

    template<class DataType>
    BBox<DataType> BSphere<DataType>::bbox ( ) const
    {
        BBox<DataType> bbox ( origin_.size ( ) );
        bbox.add_point ( origin_ );
        bbox.uniform_extension ( radius_ );
        return bbox;
    }

    template<class DataType>
    int BSphere<DataType>::get_dim ( ) const
    {
        return origin_.size ( );
    }

    template<class DataType>
    std::vector<DataType> BSphere<DataType>::get_origin ( ) const
    {
        return origin_;
    }

    template<class DataType>
    DataType BSphere<DataType>::get_radius ( ) const
    {
        return radius_;
    }

    template<class DataType>
    DataType BSphere<DataType>::get_diameter ( ) const
    {
        return radius_ * 2;
    }

    template<class DataType>
    void BSphere<DataType>::print ( std::ostream& os ) const
    {
        os << "[(";
        for ( int i = 0; i < origin_.size ( ); ++i )
        {
            os << origin_[i];
            if ( i < origin_.size ( ) - 1 )
            {
                os << ", ";
            }
        }
        os << "), "
                << radius_
                << "]\n";
    }

    template class BSphere<double>;
    template class BSphere<float>;
}

