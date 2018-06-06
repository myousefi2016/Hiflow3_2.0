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

#ifndef HIFLOW_COMMON_BSPHERE_H
#    define HIFLOW_COMMON_BSPHERE_H

/// \author Jonas Kratzke

#    include <iosfwd>
#    include <vector>
#    include "bbox.h"

namespace hiflow
{
    ///
    /// \brief   A simple sphere in any dimension.
    ///
    /// \details This class constructs a bounding sphere around a set
    /// of points. Points can be added by the member function
    /// add_point(s). One can choose, whether the origin of the sphere
    /// is to be kept or a simple algorithm by Jack Ritter[1990] should
    /// be applied to approximately get the smallest bounding sphere.
    /// It can be checked, whether two bounding spheres intersect.
    ///

    template<class DataType>
    class BSphere
    {
      public:

        /// \brief   Construction by an origin and radius
        BSphere ( const std::vector<DataType>& origin, DataType radius = DataType ( 0 ) );

        /// \brief   Construction by a set of points
        BSphere ( int dim, const std::vector<DataType>& pts );

        /// \brief   Extension to additional points, given sequentially
        void add_points ( const std::vector<DataType>& pts, bool fixed_origin = true );

        /// \brief   Radial extension
        void radial_extension ( DataType extension );

        /// \brief   Check if two spheres intersect
        bool intersects ( const BSphere<DataType>& other ) const;

        /// \brief   Check if a box intersects the sphere
        bool intersects ( const BBox<DataType>& box ) const;

        /// \brief   Check if a box is completely contained in the sphere
        bool contains ( const BBox<DataType>& box ) const;

        /// \brief   Compute bounding box
        BBox<DataType> bbox ( ) const;

        int get_dim ( ) const;
        std::vector<DataType> get_origin ( ) const;
        DataType get_radius ( ) const;
        DataType get_diameter ( ) const;

        void print ( std::ostream& os ) const;

      private:
        std::vector<DataType> origin_;
        DataType radius_;
    };
}

#endif
