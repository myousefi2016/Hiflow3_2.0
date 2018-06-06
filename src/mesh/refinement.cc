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

#include "refinement.h"

#include <algorithm>
#include <cassert>
#include <stack>
#include <utility>

#include "common/log.h"

#include "cell_type.h"
#include "entity.h"
#include "mesh_database.h"
#include "types.h"
#include "mesh/periodicity_tools.h"

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace mesh
    {

        //////////////// RefinedCell ////////////////

        const CellType& RefinedCell::cell_type ( ) const
        {
            return *cell_type_;
        }

        bool RefinedCell::is_refined ( ) const
        {
            return num_children ( ) != 0;
        }

        bool RefinedCell::is_root ( ) const
        {
            return parent_ == -1;
        }

        int RefinedCell::parent ( ) const
        {
            return parent_;
        }

        int RefinedCell::child_index ( int i ) const
        {
            range_check ( children_, i );
            return children_[i];
        }

        int RefinedCell::num_children ( ) const
        {
            return children_.size ( );
        }

        int RefinedCell::sub_cell_number ( ) const
        {
            return sub_cell_number_;
        }

        int RefinedCell::quadrant_number ( ) const
        {
            return quadrant_number_;
        }

        const std::vector<int>& RefinedCell::vertex_numbers ( ) const
        {
            assert ( !is_root ( ) );
            return vertex_numbers_;
        }

        void RefinedCell::set_coordinates ( const std::vector<Coordinate>& coords )
        {
            vertices_coords_ = coords;
        }

        void RefinedCell::get_coordinates ( std::vector<Coordinate>& coords ) const
        {
            coords = vertices_coords_;
        }

        RefinedCell::RefinedCell ( const CellType& cell_type )
        : cell_type_ ( &cell_type ), parent_ ( -1 ), sub_cell_number_ ( 0 )
        {
        }

        RefinedCell::RefinedCell ( const CellType& cell_type,
                                   const std::vector<int>& vertex_numbers,
                                   int parent, int sub_cell_number, int quadrant_number )
        : cell_type_ ( &cell_type ), parent_ ( parent ),
        sub_cell_number_ ( sub_cell_number ), vertex_numbers_ ( vertex_numbers ), quadrant_number_ ( quadrant_number )
        {
        }

        void RefinedCell::replace_child ( int old_child, int new_child )
        {
            std::vector<int>::iterator it = std::find ( children_.begin ( ), children_.end ( ), old_child );
            assert ( it != children_.end ( ) );
            *it = new_child;
        }

        //////////////// RefinementTree ////////////////

        RefinementTree::RefinementTree ( TDim tdim, const CellType& root_cell_type )
        : tdim_ ( tdim )
        {
            refined_cells_.push_back ( RefinedCell ( root_cell_type ) );
        }

        RefinementTree::RefinementTree ( TDim tdim, const CellType& root_cell_type, int levels, int num_children )
        : tdim_ ( tdim )
        {
            refined_cells_.push_back ( RefinedCell ( root_cell_type ) );

            std::vector<int> cur_leafs ( 1, 0 );
            std::vector<int> sub_cell_numbers ( num_children );
            for ( int l = 0; l < num_children; ++l )
            {
                sub_cell_numbers[l] = l + 1;
            }

            // Loop over requested uniform refinement steps
            for ( int jl = 0; jl < levels; ++jl )
            {
                std::vector<int> tmp_leafs;
                // Loop over current leafs
                for ( int jn = 0; jn < cur_leafs.size ( ); ++jn )
                {
                    this->add_subtree ( cur_leafs[jn], sub_cell_numbers );

                    for ( int jc = 1; jc <= num_children; ++jc )
                    {
                        tmp_leafs.push_back ( this->num_cells ( ) - jc );
                    }
                }
                cur_leafs.clear ( );
                cur_leafs.insert ( cur_leafs.end ( ), tmp_leafs.begin ( ), tmp_leafs.end ( ) );
            }
        }

        void RefinementTree::add_subtree ( int parent, const std::vector<int>& sub_cell_numbers )
        {
            std::vector<int> quadrant_numbers ( sub_cell_numbers.size ( ), -1 );
            add_subtree ( parent, sub_cell_numbers, quadrant_numbers );
        }

        void RefinementTree::add_subtree ( int parent, const std::vector<int>& sub_cell_numbers, const std::vector<int>& quadrant_numbers )
        {
            // get reference to parent cell (stupid name, since we use
            // another reference in loop further down)
            const RefinedCell& parent_cell_1 = cell_ref ( parent );
            assert ( !parent_cell_1.is_refined ( ) );

            const CellType& parent_type = parent_cell_1.cell_type ( );

            const int num_sub_cells = sub_cell_numbers.size ( );
            assert ( num_sub_cells == quadrant_numbers.size ( ) );

            for ( int i = 0; i < num_sub_cells; ++i )
            {
                // need to get parent_cell on every iteration since
                // its position may change during push_back
                RefinedCell& parent_cell = cell_ref ( parent );
                const int child_index = refined_cells_.size ( );
                parent_cell.children_.push_back ( child_index );

                const int cell_number = sub_cell_numbers[i];
                const int quadrant_number = quadrant_numbers[i];

                // find cell type and vertex numbers in parent cell type.
                const std::vector<int> vertex_numbers =
                        parent_type.local_vertices_of_entity ( tdim_, cell_number );
                const CellType& cell_type = CellType::get_instance ( tdim_, vertex_numbers.size ( ) );

                // NB: we do not fill the parent_facets entry of the
                // RefinedCell structure, since it does not seem to be
                // very useful. If it turns out not to be needed, it
                // should be removed as soon as possible.

                refined_cells_.push_back ( RefinedCell ( cell_type, vertex_numbers, parent, cell_number, quadrant_number ) );
            }
        }

#if 0 // does not compile due to swap() conflict with const CellType& member of RefinedCell.
        // TODO Find work-around for this

        void RefinementTree::remove_subtree ( int parent )
        {
            assert ( refined_cells_[parent].is_refined ( ) );
            typedef std::vector<int>::iterator IndexIterator;

            RefinedCell& parent_cell = cell_ref ( parent );
            for ( int i = 0; i < parent_cell.children_.size ( ); ++i )
            {
                int child_index = parent_cell.children_[i];

                if ( is_refined ( child_index ) )
                {
                    // recursively remove subtree under me
                    remove_subtree ( child_index );
                }

                // renew child_index since it might have changed during recursive call
                child_index = parent_cell.children_[i];
                const RefinedCell& child = cell_ref ( child_index );

                const int last_index = refined_cells_.size ( ) - 1;
                RefinedCell& last_cell = refined_cells_.back ( );

                refined_cells_[last_cell.parent_].replace_child ( last_index, child_index );

                std::swap ( refined_cells_[child_index], last_cell );

                refined_cells_.pop_back ( );
            }
            parent_cell.children_.clear ( );
        }
#endif

        const RefinedCell& RefinementTree::cell ( int i ) const
        {
            assert ( i >= 0 );
            assert ( i < static_cast < int > ( refined_cells_.size ( ) ) );
            return refined_cells_[i];
        }

        RefinedCell& RefinementTree::cell_ref ( int i )
        {
            assert ( i >= 0 );
            assert ( i < static_cast < int > ( refined_cells_.size ( ) ) );
            return refined_cells_[i];
        }

        int RefinementTree::num_cells ( ) const
        {
            return refined_cells_.size ( );
        }

        int RefinementTree::num_children ( int parent ) const
        {
            assert ( parent >= 0 );
            assert ( parent < static_cast < int > ( refined_cells_.size ( ) ) );
            return refined_cells_[parent].children_.size ( );
        }

        bool RefinementTree::is_root ( int cell ) const
        {
            assert ( cell >= 0 );
            assert ( cell < static_cast < int > ( refined_cells_.size ( ) ) );
            return cell == 0;
        }

        bool RefinementTree::is_refined ( int cell ) const
        {
            assert ( cell >= 0 );
            assert ( cell < static_cast < int > ( refined_cells_.size ( ) ) );
            return refined_cells_[cell].is_refined ( );
        }

        void compute_refined_vertices ( const CellType& cell_type, GDim gdim,
                                        const std::vector<Coordinate>& cell_coords,
                                        std::vector<Coordinate>& refined_coords )
        {

            const int num_refined_vertices = cell_type.num_vertices ( );
            refined_coords.resize ( num_refined_vertices * gdim, 0.0 );

            for ( int v = 0; v < num_refined_vertices; ++v )
            {
                const std::vector<int>& refined_vertex = cell_type.refined_vertex ( v );
                const int num_super_vertices = refined_vertex.size ( );

                // compute v:th refined vertex as barycenter of super vertices
                for ( int i = 0; i < num_super_vertices; ++i )
                {
                    const int super_v = refined_vertex[i];
                    for ( int c = 0; c < gdim; ++c )
                    {
                        refined_coords[v * gdim + c] += cell_coords[gdim * super_v + c];
                    }
                }

                // normalize barycenter coordinates
                for ( int c = 0; c < gdim; ++c )
                {
                    refined_coords[v * gdim + c] /= static_cast < Coordinate > ( num_super_vertices );
                }
            }
        }

        void standard_refined_geometry_fun ( const Entity& cell, std::vector<Coordinate>& refined_coords )
        {
            assert ( !cell.mesh ( )->is_periodic ( ) );

            std::vector<Coordinate> coords;
            cell.get_coordinates ( coords );
            compute_refined_vertices ( cell.cell_type ( ), cell.gdim ( ), coords, refined_coords );
        }

        void recursive_refined_geometry_fun ( const Entity& root_cell, const RefinedCell& cell, std::vector<Coordinate>& refined_coords )
        {
            assert ( !root_cell.mesh ( )->is_periodic ( ) );

            std::vector<Coordinate> coords;
            cell.get_coordinates ( coords );
            compute_refined_vertices ( cell.cell_type ( ), root_cell.gdim ( ), coords, refined_coords );
        }

        void periodic_refined_geometry_fun ( const Entity& cell, std::vector<Coordinate>& refined_coords )
        {
            assert ( cell.mesh ( )->is_periodic ( ) );

            std::vector<Coordinate> coords;
            cell.get_coordinates ( coords );
            const std::vector<MasterSlave>& period = cell.mesh ( )->get_period ( );

            coords = unperiodify ( coords, cell.gdim ( ), period );
            compute_refined_vertices ( cell.cell_type ( ), cell.gdim ( ), coords, refined_coords );
            refined_coords = periodify ( refined_coords, cell.gdim ( ), period );
        }

        void recursive_periodic_refined_geometry_fun ( const Entity& root_cell, const RefinedCell& cell, std::vector<Coordinate>& refined_coords )
        {
            assert ( root_cell.mesh ( )->is_periodic ( ) );

            std::vector<Coordinate> coords;
            cell.get_coordinates ( coords );
            const std::vector<MasterSlave>& period = root_cell.mesh ( )->get_period ( );

            coords = unperiodify ( coords, root_cell.gdim ( ), period );
            compute_refined_vertices ( cell.cell_type ( ), root_cell.gdim ( ), coords, refined_coords );
            refined_coords = periodify ( refined_coords, root_cell.gdim ( ), period );
        }

        void compute_refined_cells ( const Entity& cell,
                                     const RefinementTree& tree,
                                     std::vector<Coordinate>& refined_cell_coordinates,
                                     std::vector<int>& refined_cell_sizes,
                                     std::vector<int>& sub_cell_numbers,
                                     RefinedGeometryFunction coord_fun )
        {
            // Non-recursive version -> computes only direct sub-cells.

            assert ( tree.is_refined ( 0 ) ); // root of tree should be refined

            // Use top-level RefinedCell from RefinementTree.
            const RefinedCell& parent_cell = tree.cell ( 0 );
            const GDim gdim = cell.gdim ( );

            // Compute coords of refined vertices.
            std::vector<Coordinate> refined_coords;
            coord_fun ( cell, refined_coords );

            for ( int i = 0; i < parent_cell.num_children ( ); ++i )
            {
                // get child cell
                const int child_index = parent_cell.child_index ( i );
                const RefinedCell& child_cell = tree.cell ( child_index );

                // get refined coordinates for this child cell
                const std::vector<int>& vertex_numbers = child_cell.vertex_numbers ( );
                std::vector<Coordinate> child_coords ( vertex_numbers.size ( ) * gdim );

                // copy the refined vertex coordinates of this cell into separate vector
                for ( int i = 0; i < static_cast < int > ( vertex_numbers.size ( ) ); ++i )
                {
                    for ( int c = 0; c < gdim; ++c )
                    {
                        child_coords[i * gdim + c] = refined_coords[vertex_numbers[i] * gdim + c];
                    }
                }

                assert ( !child_cell.is_refined ( ) );

                // add to output variables.
                refined_cell_sizes.push_back ( vertex_numbers.size ( ) );
                refined_cell_coordinates.insert ( refined_cell_coordinates.end ( ),
                                                  child_coords.begin ( ), child_coords.end ( ) );
                sub_cell_numbers.push_back ( child_cell.sub_cell_number ( ) );
            }
        }

        void compute_refined_cells_recursive ( const Entity& root_cell,
                                               RefinementTree& tree,
                                               std::vector<Coordinate>& refined_cell_coordinates,
                                               std::vector<int>& refined_cell_sizes,
                                               std::vector<int>& sub_cell_numbers,
                                               std::vector<int>& tree_node_numbers,
                                               bool leaf_only,
                                               RecursiveRefinedGeometryFunction coord_fun )
        {
            assert ( tree.is_refined ( 0 ) ); // root of tree should be refined
            sub_cell_numbers.clear ( );
            tree_node_numbers.clear ( );

            const GDim gdim = root_cell.gdim ( );

            std::stack<int> parent_cell_index;
            parent_cell_index.push ( 0 ); // root cell index

            //std::stack< std::vector<Coordinate> > parent_cell_coords;
            std::vector<Coordinate> coords;
            root_cell.get_coordinates ( coords );
            //parent_cell_coords.push(coords);

            RefinedCell& root_ref_cell = tree.cell_ref ( 0 );
            root_ref_cell.set_coordinates ( coords );

            // traverse tree depth-first
            while ( !parent_cell_index.empty ( ) )
            {
                const int parent_index = parent_cell_index.top ( );
                const RefinedCell& parent_cell = tree.cell ( parent_index );

                parent_cell_index.pop ( );

                //const std::vector<Coordinate>& parent_coords = parent_cell_coords.top();
                //const CellType& parent_type = parent_cell.cell_type();

                // compute refined vertices
                std::vector<Coordinate> refined_coords;
                coord_fun ( root_cell, parent_cell, refined_coords );
                //compute_refined_vertices(parent_type, gdim, parent_coords, refined_coords);

                LOG_DEBUG ( 2, "refined_coords = " << string_from_range ( refined_coords.begin ( ), refined_coords.end ( ) ) );

                //parent_cell_coords.pop();

                for ( int i = 0; i < parent_cell.num_children ( ); ++i )
                {
                    // get child cell
                    const int child_index = parent_cell.child_index ( i );
                    RefinedCell& child_cell = tree.cell_ref ( child_index );

                    // get refined coordinates for this child cell
                    const std::vector<int>& vertex_numbers = child_cell.vertex_numbers ( );
                    std::vector<Coordinate> child_coords ( vertex_numbers.size ( ) * gdim );

                    // copy the refined vertex coordinates of this cell into separate vector
                    for ( int i = 0; i < static_cast < int > ( vertex_numbers.size ( ) ); ++i )
                    {
                        for ( int c = 0; c < gdim; ++c )
                        {
                            child_coords[i * gdim + c] = refined_coords[vertex_numbers[i] * gdim + c];
                        }
                    }

                    // store data of child cell
                    if ( !leaf_only || !child_cell.is_refined ( ) )
                    {
                        refined_cell_sizes.push_back ( vertex_numbers.size ( ) );
                        refined_cell_coordinates.insert ( refined_cell_coordinates.end ( ), child_coords.begin ( ), child_coords.end ( ) );
                        sub_cell_numbers.push_back ( child_cell.sub_cell_number ( ) );
                        tree_node_numbers.push_back ( child_index );
                    }

                    if ( child_cell.is_refined ( ) )
                    {
                        // child is not a leaf cell
                        // add child cell to stack
                        parent_cell_index.push ( child_index );
                        child_cell.set_coordinates ( child_coords );

                        // add child coords to stack (swap for performance)
                        //parent_cell_coords.push(std::vector<Coordinate>());
                        //parent_cell_coords.top().swap(child_coords);
                    }
                }
            }

            assert ( parent_cell_index.empty ( ) );
        }
    }
} // namespace hiflow
