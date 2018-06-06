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

#ifndef HIFLOW_COMMON_GRID_H
#    define HIFLOW_COMMON_GRID_H

/// \author Jonas Kratzke, Jonathan Schwegler

#    include "common/bbox.h"
#    include "common/bsphere.h"
#    include "mesh/cell_type.h"

#    include <iosfwd>
#    include <vector>

namespace hiflow
{

    /// \brief Description of a regular 2d or 3d grid.
    ///
    /// \details Given a bounding box, this class constructs a
    /// regular grid in its extensions. In general, the grid
    /// can be of the form of all available cell_tags.
    ///
    /// Currently full functionality is only given for rectangular
    /// 2d and 3d grids.
    ///

    template<class DataType>
    class Grid
    {
      public:
        /// \brief   Constructs a rectangular grid with a BBox and a number of intervals on each axis.
        Grid ( const std::vector<int> num_intervals, const BBox<DataType> bbox );

        /// \brief   Constructs a regular grid with a BBox and a number of intervals on each axis.
        Grid ( const mesh::CellType::Tag cell_tag, const std::vector<int> num_intervals, const BBox<DataType> bbox );

        int get_gdim ( ) const;
        std::vector<int> get_num_intervals ( ) const;
        int get_num_points ( ) const;
        int get_num_cells ( ) const;
        BBox<DataType> get_bbox ( ) const;

        /// \brief   Spacing on the respective axis.
        DataType delta ( int i ) const;

        /// \brief   Get the ids of the vertices of cell i
        const std::vector<int>& vertices_of_cell ( int i );

        /// \brief   Get the coordinates of the grid
        const std::vector<DataType>& coords ( );

        /// \brief   Get the coordinates of a vertex
        std::vector<DataType> coords ( std::vector<int>& indices ) const;

        /// \brief   Get a rectangular grid cell as a box
        BBox<DataType> box ( std::vector<int>& indices ) const;

        /// \brief   Find cells that intersect a given box
        void intersect ( const BBox<DataType>& bbox, std::vector<int>& cells ) const;

        /// \brief   Find cells that intersect a given sphere
        void intersect ( const BSphere<DataType>& sphere, std::vector<int>& cells ) const;

        /// \brief   Find the cell containing the point pt
        int cell_with_point ( const std::vector<DataType> pt ) const;

      private:
        void init ( );
        void compute_coords_cells ( );
        void compute1d_Line ( );
        void compute2d_Quad ( );
        void compute3d_Hex ( );
        void compute2d_Tri ( );
        void compute3d_Tet ( );
        void compute3d_Pyr ( );
        int vertex_index_Line ( int i ) const;
        int vertex_index_Quad ( int i, int j ) const;
        int vertex_index_Tri ( int i, int j ) const;
        int vertex_index_Hex ( int i, int j, int k ) const;
        int vertex_index_Tet ( int i, int j, int k ) const;
        int vertex_index_Pyr ( int i, int j, int k ) const;
        int cell_index_Line ( int i ) const;
        int cell_index_Quad ( int i, int j ) const;
        int cell_index_Tri ( int i, int j ) const;
        int cell_index_Hex ( int i, int j, int k ) const;
        int cell_index_Tet ( int i, int j, int k ) const;
        int cell_index_Pyr ( int i, int j, int k ) const;
        const int gdim_;
        const std::vector<int> num_intervals_;
        int num_points_;
        int num_cells_;
        const BBox<DataType> bbox_;
        std::vector< std::vector<int> > cells_;
        std::vector<DataType> coords_;
        const mesh::CellType::Tag tag_;
    };
}

#endif
