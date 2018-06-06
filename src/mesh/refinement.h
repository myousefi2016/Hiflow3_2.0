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

/// \author Staffan Ronnas

#ifndef HIFLOW_MESH_REFINEMENT_H
#    define HIFLOW_MESH_REFINEMENT_H

#    include <map>
#    include <tr1/unordered_map>
#    include <utility>
#    include <vector>

#    include <boost/function.hpp>

#    include "mesh/types.h"

namespace hiflow
{
    namespace mesh
    {

        class RefinementTree;
        class CellType;
        class Entity;

        /// \brief Describes a refined cell in terms of its parent cell.

        class RefinedCell
        {
          public:
            const CellType& cell_type ( ) const;

            bool is_refined ( ) const;
            bool is_root ( ) const;
            int parent ( ) const;

            int child_index ( int i ) const;
            int num_children ( ) const;
            int sub_cell_number ( ) const;
            int quadrant_number ( ) const;

            const std::vector<int>& vertex_numbers ( ) const;

            void set_coordinates ( const std::vector<Coordinate>& coords );
            void get_coordinates ( std::vector<Coordinate>& coords ) const;

          private:
            friend class RefinementTree;

            /// ctor for root cell
            explicit RefinedCell ( const CellType& cell_type );

            /// ctor for non-root cell
            RefinedCell ( const CellType& cell_type,
                          const std::vector<int>& vertex_numbers,
                          int parent, int sub_cell_number, int quadrant_number = -1 );

            void replace_child ( int old_child, int new_child );

            // Cell type for this cell (always defined)
            const CellType* cell_type_;

            // Children of this cell
            std::vector<int> children_;

            // Parent of this cell
            // = -1 for root cell
            int parent_;

            // Local number of cell in parent CellType.
            int sub_cell_number_;

            // Local vertex numbers in parent CellType.
            std::vector<int> vertex_numbers_;

            // additional index which is needed when building a refinement tree that matches an existing p4est tree
            int quadrant_number_;

            // Coordinates of vertices that form this cell
            std::vector<Coordinate> vertices_coords_;

        };

        /// \brief Describes an arbitrarily deep adaptive refinement of a cell.

        class RefinementTree
        {
          public:
            /// \brief constructor for Refinement tree containing root cell only
            RefinementTree ( TDim tdim, const CellType& root_cell_type );
            /// \brief constructor for Refinement tree containing root cell and given number of uniform refinement levels with
            /// each cell having num_children childrens
            RefinementTree ( TDim tdim, const CellType& root_cell_type, int levels, int num_children );

            void add_subtree ( int parent, const std::vector<int>& sub_cell_numbers );
            void add_subtree ( int parent, const std::vector<int>& sub_cell_numbers, const std::vector<int>& quadrant_numbers );
            //      void remove_subtree(int parent);

            const RefinedCell& cell ( int i ) const;
            int num_cells ( ) const;
            int num_children ( int parent ) const;

            bool is_root ( int cell ) const;
            bool is_refined ( int cell ) const;

          private:
            RefinedCell& cell_ref ( int i );

            TDim tdim_;
            std::vector<RefinedCell> refined_cells_;

            friend void compute_refined_cells_recursive ( const Entity& root_cell,
                                                          RefinementTree& tree,
                                                          std::vector<Coordinate>& refined_cell_coordinates,
                                                          std::vector<int>& refined_cell_sizes,
                                                          std::vector<int>& sub_cell_numbers,
                                                          std::vector<int>& tree_node_numbers,
                                                          bool leaf_only,
                                                          RecursiveRefinedGeometryFunction coord_fun );
        };

        void compute_refined_vertices ( const CellType& cell_type, GDim gdim,
                                        const std::vector<Coordinate>& cell_coords,
                                        std::vector<Coordinate>& refined_coords );

        void standard_refined_geometry_fun ( const Entity& cell, std::vector<Coordinate>& refined_coords );

        void periodic_refined_geometry_fun ( const Entity& cell, std::vector<Coordinate>& refined_coords );

        void recursive_refined_geometry_fun ( const Entity& root_cell, const RefinedCell& cell, std::vector<Coordinate>& refined_coords );

        void recursive_periodic_refined_geometry_fun ( const Entity& root_cell, const RefinedCell& cell, std::vector<Coordinate>& refined_coords );

        void compute_refined_cells ( const Entity& cell,
                                     const RefinementTree& tree,
                                     std::vector<Coordinate>& refined_cell_coordinates,
                                     std::vector<int>& refined_cell_sizes,
                                     std::vector<int>& sub_cell_numbers,
                                     RefinedGeometryFunction coord_fun = &standard_refined_geometry_fun );

        /// \brief Computes the vertex coordinates and sub-cells to be created in a refinement tree.
        /// Here, all new cells are returned, i.e. also non-leaf cells. The vector tree_node_numbers
        /// stores the index of the respective cell according to its entry in the array refined_cells of tree.
        void compute_refined_cells_recursive ( const Entity& root_cell,
                                               RefinementTree& tree,
                                               std::vector<Coordinate>& refined_cell_coordinates,
                                               std::vector<int>& refined_cell_sizes,
                                               std::vector<int>& sub_cell_numbers,
                                               std::vector<int>& tree_node_numbers,
                                               bool leaf_only,
                                               RecursiveRefinedGeometryFunction coord_fun = &recursive_refined_geometry_fun );
    }
} // namespace hiflow

#endif /* _REFINEMENT_H_ */
