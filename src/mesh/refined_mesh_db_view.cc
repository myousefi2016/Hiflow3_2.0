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

#include "refined_mesh_db_view.h"

#include <algorithm>

#include "common/log.h"

#include "cell_type.h"
#include "entity.h"
#include "iterator.h"
#include "refinement.h"

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace mesh
    {

        RefinedMeshDbView::RefinedMeshDbView ( TDim tdim, GDim gdim, MeshDatabasePtr db,
                                               const std::vector<Id>& cells,
                                               const std::vector<ConstMeshPtr>& history,
                                               const std::vector<int>& parent_history_indices,
                                               const std::vector<EntityNumber>& parent_cell_indices,
                                               const std::vector<int>& sub_cell_numbers,
                                               std::vector<MasterSlave> period )
        : MeshDbView ( tdim, gdim, db, cells, period ), history_ ( history )
        {
            const int num_cells = cells.size ( );

            assert ( static_cast < int > ( parent_history_indices.size ( ) ) == num_cells );
            assert ( static_cast < int > ( parent_cell_indices.size ( ) ) == num_cells );
            assert ( static_cast < int > ( sub_cell_numbers.size ( ) ) == num_cells );

            // Cells will contain the id:s of the cells in the order given
            // to the MeshDbView constructor, but the entities_[tdim]
            // vector will have been sorted. In order to map the data on
            // the cells to the new ordering, we compute the permutation
            // that will sort the cells vector, (this is in fact the id ->
            // index map), and then map the data vectors with this permutation.
            std::vector<int> permutation;
            compute_sorting_permutation ( cells, permutation );

            // setup attributes
            std::vector<int> sorted_history_indices ( num_cells );
            std::vector<EntityNumber> sorted_parent_cell_indices ( num_cells );
            std::vector<int> sorted_sub_cell_numbers ( num_cells );

            permute_vector ( permutation, parent_history_indices, sorted_history_indices );
            permute_vector ( permutation, parent_cell_indices, sorted_parent_cell_indices );
            permute_vector ( permutation, sub_cell_numbers, sorted_sub_cell_numbers );

            // register attributes with mesh
            parent_history_index_attr_ = AttributePtr ( new IntAttribute ( sorted_history_indices ) );
            add_attribute ( std::string ( "__parent_history_index__" ), tdim, parent_history_index_attr_ );
            parent_cell_index_attr_ = AttributePtr ( new IntAttribute ( sorted_parent_cell_indices ) );
            add_attribute ( std::string ( "__parent_cell_index__" ), tdim, parent_cell_index_attr_ );
            sub_cell_number_attr_ = AttributePtr ( new IntAttribute ( sorted_sub_cell_numbers ) );
            add_attribute ( std::string ( "__sub_cell_number__" ), tdim, sub_cell_number_attr_ );
        }

        RefinedMeshDbView::RefinedMeshDbView ( TDim tdim, GDim gdim, std::vector<MasterSlave> period )
        : MeshDbView ( tdim, gdim, period ), history_ ( std::vector<ConstMeshPtr>( 0 ) )
        {
        }

        /// \return Parent cell id of cell with index cell_index, or -1 if the cell has no parent

        Id RefinedMeshDbView::get_parent_cell_id ( EntityNumber cell_index ) const
        {
            const int parent_history_index = parent_history_index_attr_->get_int_value ( cell_index );
            if ( parent_history_index >= 0 )
            {
                const EntityNumber parent_cell_index = parent_cell_index_attr_->get_int_value ( cell_index );
                const ConstMeshPtr parent_mesh = get_ancestor_mesh ( parent_history_index );
                return parent_mesh->get_id ( tdim ( ), parent_cell_index );
            }
            else
            {
                return -1;
            }
        }

        Entity RefinedMeshDbView::get_parent_cell ( EntityNumber cell_index ) const
        {
            const int parent_history_index = parent_history_index_attr_->get_int_value ( cell_index );
            if ( parent_history_index >= 0 )
            {
                const EntityNumber parent_cell_index = parent_cell_index_attr_->get_int_value ( cell_index );
                const ConstMeshPtr parent_mesh = get_ancestor_mesh ( parent_history_index );
                return parent_mesh->get_entity ( tdim ( ), parent_cell_index );
            }
            else
            {
                assert ( false );
            }

            return Entity ( );
        }

        bool RefinedMeshDbView::cell_has_parent ( EntityNumber cell_index ) const
        {
            const int parent_history_index = parent_history_index_attr_->get_int_value ( cell_index );
            return parent_history_index >= 0;
        }

        ConstMeshPtr RefinedMeshDbView::get_ancestor_mesh ( int history_index ) const
        {
            assert ( history_index >= 0 );
            assert ( history_index < static_cast < int > ( history_.size ( ) ) );
            return history_[history_index];
        }

        std::vector<ConstMeshPtr> RefinedMeshDbView::get_all_ancestors ( ) const
        {
            return history_;
        }

        void RefinedMeshDbView::set_history ( std::vector<ConstMeshPtr> history )
        {
            history_ = history;
        }

        std::pair<ConstMeshPtr, EntityNumber> RefinedMeshDbView::get_cell_parent ( EntityNumber cell_index ) const
        {
            const int parent_history_index = parent_history_index_attr_->get_int_value ( cell_index );
            if ( parent_history_index >= 0 )
            {
                const ConstMeshPtr parent_mesh = get_ancestor_mesh ( parent_history_index );
                const EntityNumber parent_cell_index = parent_cell_index_attr_->get_int_value ( cell_index );
                return std::make_pair ( parent_mesh, parent_cell_index );
            }
            else
            {
                return std::make_pair ( MeshPtr ( ), -1 );
            }
        }

        MeshPtr RefinedMeshDbView::refine ( std::vector<int>& refinements ) const
        {
            assert ( static_cast < int > ( refinements.size ( ) ) == num_entities ( tdim ( ) ) );

            const TDim td = tdim ( );
            const GDim gd = gdim ( );

            // History index of current mesh
            const EntityCount num_cells = num_entities ( td );

            // Robust coarsening: only coarsen a cell if all children of
            // the cell are marked to be coarsened. This requires some
            // preprocessing of the refinements vector. First all
            // ancestors of cells that are not marked to be coarsened are
            // traversed and put into the array forbidden_ancestors. These
            // ancestors cannot be coarsened to. Then, a loop over all
            // cells that are marked to be coarsened is performed, in
            // which the function find_coarsenable_parent is used to
            // find out whether the parent of the cell can be coarsened to.
            SortedArray<Id> forbidden_ancestors;
            for ( int c = 0; c < num_cells; ++c )
            {
                if ( refinements[c] >= 0 )
                {
                    mark_forbidden_ancestor_cells ( c, forbidden_ancestors );
                }
            }

            // Correct refinement mark for cells that are marked to be coarsened.
            std::tr1::unordered_map< EntityNumber, std::pair<int, EntityNumber> > coarsening_data;
            for ( int c = 0; c < num_cells; ++c )
            {
                if ( refinements[c] < 0 )
                {
                    int parent_history_index;
                    EntityNumber parent_cell_index;
                    if ( find_coarsenable_parent ( c, forbidden_ancestors,
                                                   &parent_history_index, &parent_cell_index ) )
                    {
                        coarsening_data[c] = std::make_pair ( parent_history_index, parent_cell_index );
                    }
                    else
                    {
                        // do not refine cell
                        refinements[c] = 0;
                    }
                }
            }

            // History index of current mesh
            const int current_history_index = history_.size ( );
            LOG_DEBUG ( 1, "Current history index = " << current_history_index );

            // Cells to be added to refined mesh
            std::vector<Id> new_cells;

            // Vectors for parent data
            std::vector<int> parent_history_indices;
            std::vector<EntityNumber> parent_cell_indices;
            std::vector<int> sub_cell_numbers;

            SortedArray<Id> visited_coarsened_parents;

            // Loop over cells and refine, coarsen, or keep according to refinements marker
            for ( int c = 0; c < num_cells; ++c )
            {
                const Entity current_cell = get_entity ( td, c );
                const int r = refinements[c];

                if ( r > 0 )
                { // refinement (subdivision)
                    LOG_DEBUG ( 2, "Refining cell " << c << " with refinement " << r );

                    // compute children cells
                    const RefinementTree* tree = current_cell.cell_type ( ).refinement_tree ( r - 1 );
                    std::vector<int> local_sub_cell_numbers;
                    const std::vector<Id> children_cells = compute_refined_cells ( current_cell, tree, local_sub_cell_numbers );
                    new_cells.insert ( new_cells.end ( ), children_cells.begin ( ), children_cells.end ( ) );
                    sub_cell_numbers.insert ( sub_cell_numbers.end ( ), local_sub_cell_numbers.begin ( ), local_sub_cell_numbers.end ( ) );

                    // set parent data
                    for ( int rc = 0; rc < static_cast < int > ( children_cells.size ( ) ); ++rc )
                    {
                        parent_history_indices.push_back ( current_history_index );
                        parent_cell_indices.push_back ( c );
                    }

                }
                else if ( r < 0 )
                { // coarsening
                    const Id parent_id = get_parent_cell_id ( c );
                    if ( !visited_coarsened_parents.find_insert ( parent_id ) )
                    {
                        // each parent cell should only be added once
                        new_cells.push_back ( parent_id );

                        // set parent data
                        assert ( coarsening_data.count ( c ) != 0 );
                        std::pair<int, EntityNumber> parent_data = coarsening_data[c];

                        if ( parent_data.first > 0 )
                        {
                            // the parent cell is not in root mesh -> get its sub_cell_number

                            // get pointer to ancestor mesh (it is a RefinedMeshDbView* since it is not the root mesh).
                            const RefinedMeshDbView* ancestor_mesh =
                                    static_cast < const RefinedMeshDbView* > ( get_ancestor_mesh ( parent_data.first ).get ( ) );

                            // get parent_history_index from attribute
                            const int parent_history_index =
                                    ancestor_mesh->parent_history_index_attr_->get_int_value ( parent_data.second );
                            parent_history_indices.push_back ( parent_history_index );

                            // get parent_cell_index from attribute
                            const int parent_cell_index =
                                    ancestor_mesh->parent_cell_index_attr_->get_int_value ( parent_data.second );
                            parent_cell_indices.push_back ( parent_cell_index );

                            // get sub_cell_number from attribute
                            const int sub_cell_number =
                                    ancestor_mesh->sub_cell_number_attr_->get_int_value ( parent_data.second );
                            sub_cell_numbers.push_back ( sub_cell_number );
                        }
                        else
                        {
                            // ancestor cell in root mesh -> has not been refined.
                            parent_history_indices.push_back ( -1 );
                            parent_cell_indices.push_back ( -1 );
                            sub_cell_numbers.push_back ( -1 );
                        }
                    }
                }
                else
                { // keep cell
                    const Id cell_id = current_cell.id ( );
                    new_cells.push_back ( cell_id );

                    // inherit parent data from current cell
                    int parent_history_index;
                    current_cell.get ( "__parent_history_index__", &parent_history_index );
                    parent_history_indices.push_back ( parent_history_index );

                    EntityNumber parent_cell_index;
                    current_cell.get ( "__parent_cell_index__", &parent_cell_index );
                    parent_cell_indices.push_back ( parent_cell_index );

                    int sub_cell_number;
                    current_cell.get ( "__sub_cell_number__", &sub_cell_number );
                    sub_cell_numbers.push_back ( sub_cell_number );
                }
            }

            // allocate and return new mesh
            std::vector<ConstMeshPtr> refined_history ( history_ );
            refined_history.push_back ( ConstMeshPtr ( this ) );

            MeshPtr refined_mesh = MeshPtr ( new RefinedMeshDbView ( td, gd, get_db ( ), new_cells,
                                                                     refined_history,
                                                                     parent_history_indices,
                                                                     parent_cell_indices,
                                                                     sub_cell_numbers,
                                                                     period_ ) );

            return refined_mesh;
        }

        void RefinedMeshDbView::mark_forbidden_ancestor_cells ( EntityNumber cell_index,
                                                                SortedArray<Id>& ancestors ) const
        {
            const TDim td = tdim ( );
            const RefinedMeshDbView* mesh_ptr = this;

            // Starting from current mesh, loop over and mark ancestors of
            // given cell. Break out of loop when the root is reached, or
            // when an already visited parent is found.
            while ( true )
            {
                int parent_history_index = mesh_ptr->parent_history_index_attr_->get_int_value ( cell_index );
                if ( parent_history_index >= 0 )
                {
                    ConstMeshPtr parent_mesh = get_ancestor_mesh ( parent_history_index );
                    const int parent_cell_index = mesh_ptr->parent_cell_index_attr_->get_int_value ( cell_index );
                    const Id parent_id = parent_mesh->get_id ( td, parent_cell_index );

                    if ( ancestors.find_insert ( parent_id ) )
                    {
                        // Parent already marked, and therefore all its
                        // ancestors too. Hence, we are done.
                        break;
                    }

                    if ( parent_history_index > 0 )
                    {
                        // Update variables for next iteration. Static
                        // cast is OK since we know that pointers in
                        // history_ are RefinedMeshDbView* for parent_history_index > 0.
                        mesh_ptr = static_cast < const RefinedMeshDbView* > ( parent_mesh.get ( ) );
                        cell_index = parent_cell_index;
                    }
                    else
                    {
                        // Parent is in first mesh in the history, so we stop
                        // after having visited it.
                        break;
                    }
                }
                else
                {
                    // Cell is root cell, so we don't have any ancestors.
                    break;
                }
            }
        }

        bool RefinedMeshDbView::find_coarsenable_parent ( EntityNumber cell_index,
                                                          const SortedArray<Id>& forbidden_ancestors,
                                                          int* parent_history_index,
                                                          EntityNumber* parent_cell_index ) const
        {

            const TDim td = tdim ( );
            const RefinedMeshDbView* mesh_ptr = this;
            bool can_be_coarsened = false;

            // Check if parent of current cell (if it exists) is an ancestor that can
            // be the target of coarsening. If so, return true, and store attributes
            // of parent in parent_history_index and parent_cell_index.
            int p_history_index = mesh_ptr->parent_history_index_attr_->get_int_value ( cell_index );
            if ( p_history_index >= 0 )
            {
                ConstMeshPtr parent_mesh = get_ancestor_mesh ( p_history_index );
                const int p_cell_index = mesh_ptr->parent_cell_index_attr_->get_int_value ( cell_index );
                const Id parent_id = parent_mesh->get_id ( td, p_cell_index );

                if ( !forbidden_ancestors.find ( parent_id, 0 ) )
                {
                    can_be_coarsened = true;

                    if ( parent_history_index != 0 )
                    {
                        *parent_history_index = p_history_index;
                    }

                    if ( parent_cell_index != 0 )
                    {
                        *parent_cell_index = p_cell_index;
                    }
                }
            }

            return can_be_coarsened;
        }

        bool RefinedMeshDbView::is_boundary_facet ( EntityNumber facet_index ) const
        {
            LOG_DEBUG ( 3, "Is facet (index) " << facet_index << " a boundary facet?" );
            const TDim facet_tdim = tdim ( ) - 1;

            assert ( facet_index >= 0 );
            assert ( facet_index < num_entities ( facet_tdim ) );

            // get facet entity
            const Entity facet_entity = get_entity ( facet_tdim, facet_index );

            // get incident cell(s) of facet
            const std::vector<Entity> incident_cells ( begin_incident ( facet_entity, tdim ( ) ),
                                                       end_incident ( facet_entity, tdim ( ) ) );

            bool is_boundary = false;

            if ( incident_cells.size ( ) == 1 )
            {
                Entity cell_entity = incident_cells.front ( );
                std::pair<ConstMeshPtr, EntityNumber> parent_mesh_index = get_cell_parent ( cell_entity.index ( ) );

                ConstMeshPtr parent_mesh = parent_mesh_index.first;
                const EntityNumber parent_cell_index = parent_mesh_index.second;

                if ( parent_cell_index < 0 )
                {
                    // Cell is a root cell -- hence the facet could be on the boundary,
                    // or be a "big" inner facet.
                    is_boundary = true;

                    // Quick fix to exclude "big" inner facets: material number must be
                    // non-negative.
                    if ( facet_entity.get_material_number ( ) < 0 )
                    {
                        is_boundary = false;
                    }
                }
                else
                {
                    // cell has parent -- check if parent facet is on the boundary.

                    // get cell type of parent cell
                    const Entity parent_cell_entity = parent_mesh->get_entity ( tdim ( ), parent_cell_index );
                    const CellType& parent_cell_type = parent_cell_entity.cell_type ( );

                    // get local vertices of cell
                    const std::vector<int> cell_vertices ( cell_entity.begin_vertex_ids ( ), cell_entity.end_vertex_ids ( ) );
                    LOG_DEBUG ( 3, "Cell vertices: " << string_from_range ( cell_vertices.begin ( ), cell_vertices.end ( ) ) );

                    const CellType& cell_type = cell_entity.cell_type ( );

                    // get facet vertices
                    const std::vector<int> facet_vertices ( facet_entity.begin_vertex_ids ( ), facet_entity.end_vertex_ids ( ) );
                    LOG_DEBUG ( 3, "Facet vertices: " << string_from_range ( facet_vertices.begin ( ), facet_vertices.end ( ) ) );

                    // get local sub cell number of cell in parent cell
                    EntityNumber sub_cell_number;
                    cell_entity.get ( "__sub_cell_number__", &sub_cell_number );
                    LOG_DEBUG ( 3, "Sub cell number: " << sub_cell_number );

                    // TODO: can we exit early if sub_cell_number == 0 here?

                    // Get local facet number in child cell type by linear search, comparing the vertex sets.
                    EntityNumber facet_number_in_child = -1;
                    for ( int i = 0; i < cell_type.num_regular_entities ( facet_tdim ); ++i )
                    {
                        const std::vector<int> trial_facet_vertices
                                = cell_type.vertices_of_entity ( facet_tdim, i, cell_vertices );

                        if ( compare_as_sets ( trial_facet_vertices, facet_vertices ) )
                        {
                            facet_number_in_child = i;
                            break;
                        }
                    }

                    assert ( facet_number_in_child != -1 );

                    LOG_DEBUG ( 3, "Facet number in child cell type: " << facet_number_in_child );

                    // get facet number in parent cell type
                    const int facet_number_in_parent = parent_cell_type.sub_entities_of_cell ( facet_tdim, sub_cell_number )[facet_number_in_child];
                    LOG_DEBUG ( 3, "Facet number in parent cell type: " << facet_number_in_parent );
                    assert ( facet_number_in_parent >= 0 );

                    // If child facet is regular, then we have already found the facet number to use.
                    EntityNumber parent_facet_number = facet_number_in_parent;
                    if ( parent_facet_number >= parent_cell_type.num_regular_entities ( facet_tdim ) )
                    {
                        // Else, if facet is an irregular facet in the parent, we have to get its regular parent
                        parent_facet_number = parent_cell_type.regular_parent ( facet_tdim, facet_number_in_parent );
                        LOG_DEBUG ( 3, "Parent facet number: " << parent_facet_number );
                        assert ( parent_facet_number == -1 || parent_facet_number < parent_cell_type.num_regular_entities ( facet_tdim ) );
                    }

                    if ( parent_facet_number != -1 )
                    {
                        IncidentEntityIterator parent_facet_it = parent_cell_entity.begin_incident ( facet_tdim );
                        for ( int k = 0; k < parent_facet_number; ++k )
                        {
                            ++parent_facet_it;
                        }
                        assert ( parent_facet_it != parent_cell_entity.end_incident ( facet_tdim ) );
                        LOG_DEBUG ( 3, "Recursive call of is_boundary_facet!" );
                        is_boundary = parent_mesh->is_boundary_facet ( parent_facet_it->index ( ) );
                    }
                    else
                    {
                        // In this case, the facet is interior in the parent, and hence not on the boundary.
                        is_boundary = false;
                    }
                }
            }

            if ( is_boundary )
            {
                LOG_DEBUG ( 3, "Facet with index " << facet_index << " is a boundary facet." );
            }
            else
            {
                LOG_DEBUG ( 3, "Facet with index " << facet_index << " is NOT a boundary facet." );
            }
            return is_boundary;
        }
    }
} // namespace hiflow
