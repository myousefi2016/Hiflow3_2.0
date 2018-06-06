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

#ifndef HIFLOW_COMMON_BBOX_H
#    define HIFLOW_COMMON_BBOX_H

/// \author Staffan Ronnas, Jonas Kratzke

#    include <iosfwd>
#    include <vector>

namespace hiflow
{
    ///
    /// \brief   A simple rectangular box in any dimension.
    ///
    /// \details This class constructs a bounding box around a set
    /// of points. Points can be added by the member function
    /// add_point(s). It can be checked, whether two bounding
    /// boxes intersect.
    ///

    template<class DataType>
    class BBox
    {
      public:

        /// \brief   Constructs an empty box
        BBox ( int dim );

        /// \brief   Construction by an array of initial extents
        BBox ( DataType* extents, int dim );

        /// \brief   Construction by a vector of initial extents
        BBox ( const std::vector<DataType>& extents );

        DataType min ( int dir ) const;
        DataType max ( int dir ) const;

        /// \brief   Extension to an additional point
        void add_point ( const std::vector<DataType>& pt );

        /// \brief   Extension to additional points, given sequentially
        void add_points ( const std::vector<DataType>& pts );

        /// \brief   Extension with a constant value in every direction
        void uniform_extension ( DataType extension );

        /// \brief   Check if two boxes intersect
        bool intersects ( const BBox<DataType>& other ) const;

        /// \brief   Compute the diagonal of the box
        DataType compute_diagonal ( ) const;

        void print ( std::ostream& os ) const;
        int get_dim ( ) const;
        std::vector<DataType> get_extents ( ) const;
        std::vector<DataType> get_vertices ( ) const;

      private:
        std::vector<DataType> extents_;
    };
}

#endif
