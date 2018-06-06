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

#include "mesh_db_view.h"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

#include "common/log.h"

#include "cell_type.h"
#include "communication.h"
#include "connectivity.h"
#include "entity.h"
#include "iterator.h"
#include "mesh_builder.h"
#include "mesh_database.h"
#include "mpi_communication.h"
#include "refined_mesh_db_view.h"
#include "refinement.h"
#include "common/hdf5_tools.h"
#include "mesh/periodicity_tools.h"

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace mesh
    {

        MeshDbView::MeshDbView ( TDim tdim, GDim gdim, MeshDatabasePtr db, const std::vector<Id>& cells, std::vector<MasterSlave> period )
        : Mesh ( tdim, gdim, period ),
        db_ ( db ),
        ref_geom_fun_ ( &standard_refined_geometry_fun ),
        rec_ref_geom_fun_ ( &recursive_refined_geometry_fun ),
        entities_ ( new SortedArray<Id>[tdim + 1] ),
        incidence_connections_ ( new Connectivity[( tdim + 1 )*( tdim + 1 )] )
        {
            LOG_DEBUG ( 1, "Creating mesh [tdim = " << tdim << ", gdim = " << gdim << "]" );
            // Assignment to entities_[tdim] performs sorting on a copy of
            // the cells vector (in SortedArray constructor).
            entities_[tdim] = SortedArray<Id>( cells );

            // initialize D->0 connectivity directly
            compute_entity_vertex_connectivity ( tdim );

            LOG_DEBUG ( 3, "Cells:" << string_from_range ( entities_[tdim].begin ( ),
                                                           entities_[tdim].end ( ) ) );
            if ( is_periodic ( ) )
            {
                ref_geom_fun_ = &periodic_refined_geometry_fun;
                rec_ref_geom_fun_ = &recursive_periodic_refined_geometry_fun;
            }
        }

        MeshDbView::MeshDbView ( TDim tdim, GDim gdim, std::vector<MasterSlave> period )
        : Mesh ( tdim, gdim, period ),
        ref_geom_fun_ ( &standard_refined_geometry_fun ),
        rec_ref_geom_fun_ ( &recursive_refined_geometry_fun ),
        entities_ ( new SortedArray<Id>[tdim + 1] ),
        incidence_connections_ ( new Connectivity[( tdim + 1 )*( tdim + 1 )] )
        {
            if ( is_periodic ( ) )
            {
                ref_geom_fun_ = &periodic_refined_geometry_fun;
                rec_ref_geom_fun_ = &recursive_periodic_refined_geometry_fun;
            }
        }

        MeshDbView::~MeshDbView ( )
        {
            delete[] incidence_connections_;
        }

        EntityIterator MeshDbView::begin ( TDim entity_dim ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim ( ) );

            const Entity entity = Entity ( ConstMeshPtr ( this ), entity_dim, 0 );
            return EntityIterator ( entity, 0 );
        }

        EntityIterator MeshDbView::end ( TDim entity_dim ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim ( ) );
            const EntityCount end_entity = num_entities ( entity_dim );
            const Entity entity = Entity ( ConstMeshPtr ( this ), entity_dim, end_entity );
            return EntityIterator ( entity, end_entity );
        }

        IncidentEntityIterator MeshDbView::begin_incident ( const Entity& entity, TDim entity_dim ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim ( ) );

            const TDim d1 = entity.tdim ( );
            assert ( d1 >= 0 );
            assert ( d1 <= tdim ( ) );

            if ( d1 == 0 )
            {
                assert ( entity_dim != 0 );
            }

            const Connectivity& connectivity = get_incidence_connectivity ( d1, entity_dim );
            const EntityNumber index = entity.index ( );

            EntityNumber new_entity_index;

            // TODO(Staffan) When can it happen that an incidence is empty?
            if ( connectivity.begin ( index ) != connectivity.end ( index ) )
            {
                new_entity_index = *( connectivity.begin ( index ) );
            }
            else
            {
                new_entity_index = -1;
            }

            const Entity new_entity = Entity ( ConstMeshPtr ( this ), entity_dim, new_entity_index );
            return IncidentEntityIterator ( new_entity, connectivity.begin ( index ) );

        }

        IncidentEntityIterator MeshDbView::end_incident ( const Entity& entity, TDim entity_dim ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim ( ) );

            const TDim d1 = entity.tdim ( );
            assert ( d1 >= 0 );
            assert ( d1 <= tdim ( ) );

            if ( d1 == 0 )
            {
                assert ( entity_dim != 0 );
            }
            const Connectivity& connectivity = get_incidence_connectivity ( d1, entity_dim );
            const EntityNumber index = entity.index ( );

            const Entity new_entity = Entity ( ConstMeshPtr ( this ), entity_dim, -1 );
            return IncidentEntityIterator ( new_entity, connectivity.end ( index ) );
        }

        Id MeshDbView::get_id ( TDim entity_dim, EntityNumber index ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim ( ) );
            assert ( index >= 0 );
#ifndef NDEBUG
            if ( index >= num_entities ( entity_dim ) )
            {
                LOG_DEBUG ( 3, "dim = " << entity_dim << ", index = "
                            << index << " num_entities = " << num_entities ( entity_dim ) );
            }
#endif
            assert ( index < num_entities ( entity_dim ) );
            return get_entities ( entity_dim )[index];
        }

        std::vector<Id> MeshDbView::get_vertex_ids ( TDim entity_dim, EntityNumber index ) const
        {
            const Id entity_id = get_id ( entity_dim, index );
            return std::vector<Id>( db_->begin_vertices ( entity_dim, entity_id ),
                    db_->end_vertices ( entity_dim, entity_id ) );
        }

        MaterialNumber MeshDbView::get_material_number ( TDim entity_dim, EntityNumber index ) const
        {
            return db_->get_material_number ( entity_dim, get_id ( entity_dim, index ) );
        }

        void MeshDbView::set_material_number ( TDim entity_dim, EntityNumber index, MaterialNumber material )
        {
            return db_->set_material_number ( entity_dim, get_id ( entity_dim, index ), material );
        }

        std::vector<Coordinate> MeshDbView::get_coordinates ( TDim entity_dim, EntityNumber index ) const
        {
            // only of size gdim if entity_dim == 0
            std::vector<Coordinate> coords ( gdim ( ) );
            if ( entity_dim == 0 )
            {
                coords = db_->get_coordinates ( get_id ( 0, index ) );
            }
            else
            {
                coords = db_->get_coordinates ( get_vertex_ids ( entity_dim, index ) );
            }

            bool displacement_x = has_attribute ( std::string ( "__displacement_vertices_x__" ), 0 );
            if ( displacement_x )
            {
                if ( entity_dim != 0 )
                {
                    VertexIdIterator vii = db_->begin_vertices ( entity_dim, get_id ( entity_dim, index ) );
                    int k = 0;
                    for (; vii != db_->end_vertices ( entity_dim, get_id ( entity_dim, index ) ); ++vii )
                    {
                        int vertex_index;
                        bool found = find_entity ( 0, *vii, &vertex_index );
                        assert ( found );

                        AttributePtr attr_ptr = get_attribute ( std::string ( "__displacement_vertices_x__" ), 0 );
                        coords[k] += attr_ptr->get_double_value ( vertex_index );
                        if ( has_attribute ( std::string ( "__displacement_vertices_y__" ), 0 ) )
                        {
                            assert ( gdim ( ) >= 2 );
                            AttributePtr attr_ptr = get_attribute ( std::string ( "__displacement_vertices_y__" ), 0 );
                            coords[k + 1] += attr_ptr->get_double_value ( vertex_index );
                            if ( has_attribute ( std::string ( "__displacement_vertices_z__" ), 0 ) )
                            {
                                assert ( gdim ( ) == 3 );
                                AttributePtr attr_ptr = get_attribute ( std::string ( "__displacement_vertices_z__" ), 0 );
                                coords[k + 2] += attr_ptr->get_double_value ( vertex_index );
                            }
                        }
                        k += gdim ( );
                    }
                }
                else
                { // entity_dim == 0
                    AttributePtr attr_ptr = get_attribute ( std::string ( "__displacement_vertices_x__" ), 0 );
                    coords[0] += attr_ptr->get_double_value ( index );
                    if ( has_attribute ( std::string ( "__displacement_vertices_y__" ), 0 ) )
                    {
                        assert ( gdim ( ) >= 2 );
                        AttributePtr attr_ptr = get_attribute ( std::string ( "__displacement_vertices_y__" ), 0 );
                        coords[1] += attr_ptr->get_double_value ( index );
                        if ( has_attribute ( std::string ( "__displacement_vertices_z__" ), 0 ) )
                        {
                            assert ( gdim ( ) == 3 );
                            AttributePtr attr_ptr = get_attribute ( std::string ( "__displacement_vertices_z__" ), 0 );
                            coords[2] += attr_ptr->get_double_value ( index );
                        }
                    }
                }
            }
            return coords;
        }

        Id MeshDbView::get_parent_cell_id ( EntityNumber cell_index ) const
        {
            return -1;
        }

        Entity MeshDbView::get_parent_cell ( EntityNumber cell_index ) const
        {
            assert ( false );
            return Entity ( );
        }

        bool MeshDbView::cell_has_parent ( EntityNumber cell_index ) const
        {
            return false;
        }

        std::vector<Id> MeshDbView::get_children_cell_ids ( EntityNumber cell_index ) const
        {
            NOT_YET_IMPLEMENTED;
            std::vector<Id> tmp;
            return tmp;
        }

        std::vector<EntityNumber> MeshDbView::get_sibling_cell_indices ( EntityNumber cell_index ) const
        {
            std::vector<EntityNumber> siblings ( 0 );
            return siblings;
        }

        Entity MeshDbView::get_entity ( TDim entity_dim, EntityNumber entity_number ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim ( ) );

            assert ( entity_number >= 0 );
            assert ( entity_number < num_entities ( entity_dim ) );

            return Entity ( ConstMeshPtr ( this ), entity_dim, entity_number );
        }

        EntityCount MeshDbView::num_entities ( TDim entity_dim ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim ( ) );
            const std::vector<Id>& entities = get_entities ( entity_dim );
            return entities.size ( );
        }

        EntityCount MeshDbView::num_incident_entities ( const Entity& entity, TDim entity_dim ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim ( ) );

            const TDim d1 = entity.tdim ( );
            assert ( d1 >= 0 );
            assert ( d1 <= tdim ( ) );

            if ( d1 == 0 )
            {
                assert ( entity_dim != 0 );
            }

            const Connectivity& connectivity = get_incidence_connectivity ( d1, entity_dim );
            const EntityNumber index = entity.index ( );
            return connectivity.num_connections ( index );
        }

        Connectivity& MeshDbView::get_incidence_connectivity ( TDim d1, TDim d2 ) const
        {
            assert ( d1 >= 0 );
            assert ( d1 <= tdim ( ) );
            assert ( d2 >= 0 );
            assert ( d2 <= tdim ( ) );
            assert ( d1 > 0 || d2 > 0 ); // 0->0 incidence not allowed

            const int index = incidence_connection_index ( d1, d2 );
            Connectivity& connectivity = incidence_connections_[index];

            if ( !has_computed_connectivity ( d1, d2 ) )
            {
                // d->0 connectivities are computed in get_entity() and should not be required here.
                assert ( d2 != 0 );
                LOG_DEBUG ( 2, "Start computing connectivity " << d1 << ", " << d2 );

                // should computed incidence be converted to local indices?
                bool should_convert_to_local = false;

                // TODO(Staffan) Can this be cleaned up?
                if ( d1 == 0 )
                {
                    const SortedArray<Id>& vertices = get_entities ( 0 );
                    const SortedArray<Id>& d2_entities = get_entities ( d2 );
                    db_->vertex_entity_incidence ( d2, vertices, d2_entities, connectivity );
                    should_convert_to_local = true;
                }
                else if ( d1 == tdim ( ) && d2 < tdim ( ) )
                {
                    // Compute d2-entities if necessary
                    ( void ) get_entities ( d2 );
                    assert ( has_computed_connectivity ( d1, d2 ) );
                    should_convert_to_local = false;
                    // the connectivity reference now points to the object that was computed in get_entities()
                }
                else if ( d1 > d2 )
                {
                    const SortedArray<Id>& entities_d1 = get_entities ( d1 );
                    const SortedArray<Id>& entities_d2 = get_entities ( d2 );
                    db_->compute_incident_entities ( d1, d2, entities_d1, entities_d2, connectivity );
                    should_convert_to_local = true;
                }
                else if ( d1 == d2 )
                {
                    const SortedArray<Id>& entities_d = get_entities ( d1 );
                    db_->compute_incident_entities ( d1, entities_d, connectivity );
                    should_convert_to_local = true;
                }
                else if ( d1 < d2 )
                {
                    const Connectivity& d2_d1_connectivity = get_incidence_connectivity ( d2, d1 );

                    const SortedArray<Id>& entities_d1 = get_entities ( d1 );
                    const EntityCount num_d1_entities = entities_d1.size ( );

                    d2_d1_connectivity.transpose ( num_d1_entities, connectivity );

                    should_convert_to_local = false;
                }
                else
                {
                    assert ( 0 );
                }

                if ( should_convert_to_local )
                {
                    // convert newly computed connectivity to local indices
                    assert ( has_computed_connectivity ( d2, 0 ) );
                    convert_connectivity_to_indices ( d2, connectivity );
                }

                LOG_DEBUG ( 2, "End computing connectivity " << d1 << ", " << d2
                            << " -> index = " << index );
            }
            return connectivity;
        }

        SortedArray<Id>& MeshDbView::get_entities ( TDim d ) const
        {
            assert ( d >= 0 );
            assert ( d <= tdim ( ) );
            assert ( !entities_[tdim ( )].empty ( ) );

            if ( entities_[d].empty ( ) )
            {
                if ( d == 0 )
                {
                    initialize_vertices ( );
                }
                else
                {
                    std::pair<TDim, TDim> dims ( tdim ( ), d );
                    Connectivity& D_d_connectivity = incidence_connections_[incidence_connection_index ( tdim ( ), d )];
                    db_->build ( d, entities_[tdim ( )], entities_[d], D_d_connectivity );
                    compute_entity_vertex_connectivity ( d );
                    convert_connectivity_to_indices ( d, D_d_connectivity );
                    LOG_DEBUG ( 2, "End computing connectivity " << tdim ( ) << ", " << d
                                << " -> index = " << incidence_connection_index ( tdim ( ), d ) );
                }
            }

            assert ( !entities_[d].empty ( ) );
            return entities_[d];
        }

        int MeshDbView::incidence_connection_index ( TDim d1, TDim d2 ) const
        {
            return d1 * ( tdim ( ) + 1 ) + d2;
        }

        bool MeshDbView::has_computed_connectivity ( TDim d1, TDim d2 ) const
        {
            return incidence_connections_[incidence_connection_index ( d1, d2 )].num_entities ( ) != 0;
        }

        void MeshDbView::initialize_vertices ( ) const
        {
            // loop cells
            assert ( entities_[0].empty ( ) );
            assert ( !entities_[tdim ( )].empty ( ) );

            const TDim mesh_tdim = tdim ( );

            for ( size_t c = 0; c != entities_[mesh_tdim].size ( ); ++c )
            {
                // loop vertices of cell
                const Id cell_id = get_id ( mesh_tdim, c );
                VertexIdIterator it = db_->begin_vertices ( mesh_tdim, cell_id );
                const VertexIdIterator end = db_->end_vertices ( mesh_tdim, cell_id );
                for (; it != end; ++it )
                {
                    // Insert v into entities_[0] in sorted order. Since
                    // entities_[0] is sorted, we can perform binary
                    // search to determine if and where to insert v
                    const Id vertex_id = *it;
                    if ( !entities_[0].find ( vertex_id, 0 ) )
                    {
                        entities_[0].insert ( vertex_id );
                    }
                }
            }
        }

        void MeshDbView::compute_entity_vertex_connectivity ( TDim d ) const
        {
            LOG_DEBUG ( 2, "Computing entity-vertex connectivity for dimension " << d );
            const SortedArray<Id>& vertices = get_entities ( 0 );
            const SortedArray<Id>& entities = entities_[d];
            assert ( !entities.empty ( ) );

            const int index = incidence_connection_index ( d, 0 );
            Connectivity& entity_vertex_connectivity = incidence_connections_[index];

            // compute connectivity cell -> vertex id:s
            for ( SortedArray<Id>::const_iterator it = entities.begin ( );
                  it != entities.end ( ); ++it )
            {
                const Id entity_id = *it;
                std::vector<Id> entity_vertices ( db_->begin_vertices ( d, entity_id ),
                                                  db_->end_vertices ( d, entity_id ) );
                // loop over vertices and convert from id to index
                for ( std::vector<Id>::iterator v_it = entity_vertices.begin ( ); v_it != entity_vertices.end ( ); ++v_it )
                {
                    EntityNumber index;
                    const bool found = vertices.find ( *v_it, &index );
                    assert ( found );
                    *v_it = index;
                }
                entity_vertex_connectivity.add_connections ( entity_vertices );
            }

        }

        bool MeshDbView::find_entity ( TDim entity_dim, Id id, EntityNumber* entity_number ) const
        {
            const SortedArray<Id>& entities = get_entities ( entity_dim );
            return entities.find ( id, entity_number );
        }

        bool MeshDbView::find_vertex ( const std::vector<Coordinate>& coords, EntityNumber* v_index ) const
        {
            assert ( static_cast < int > ( coords.size ( ) ) == gdim ( ) );

            // check if it exists in database
            Id v_id;
            const bool is_in_db = db_->has_vertex ( vec2ptr ( coords ), v_id );

            // check if it exists in this Mesh
            if ( is_in_db )
            {
                return find_entity ( 0, v_id, v_index );
            }

            return false;
        }

        MeshPtr MeshDbView::refine ( ) const
        {
            std::vector<int> all_cells ( num_entities ( tdim ( ) ), 1 );
            return refine ( all_cells );
        }

        MeshPtr MeshDbView::refine ( const std::string& attribute_name ) const
        {
            assert ( has_attribute ( attribute_name, tdim ( ) ) );
            AttributePtr attr = get_attribute ( attribute_name, tdim ( ) );
            std::vector<EntityNumber> refinements;
            refinements.reserve ( num_entities ( tdim ( ) ) );

            for ( EntityIterator it = begin ( tdim ( ) ); it != end ( tdim ( ) ); ++it )
            {
                const EntityNumber index = it->index ( );
                const int refinement_type = attr->get_int_value ( index );
                refinements.push_back ( refinement_type );
            }
            return refine ( refinements );
        }

        MeshPtr MeshDbView::refine ( std::vector<int>& refinements ) const
        {
            assert ( static_cast < int > ( refinements.size ( ) ) == num_entities ( tdim ( ) ) );

            const TDim td = tdim ( );
            const GDim gd = gdim ( );
            const EntityCount num_cells = num_entities ( td );

            // History index of current mesh (root)
            const int current_history_index = 0;

            // Cells to be added to refined mesh
            std::vector<Id> new_cells;

            // Vectors for parent data
            std::vector<int> parent_history_indices;
            std::vector<EntityNumber> parent_cell_indices;
            std::vector<int> sub_cell_numbers;

            // Loop over cells and refine or keep according to the refinement marker
            for ( int c = 0; c < num_cells; ++c )
            {
                const Entity current_cell = get_entity ( td, c );
                const int r = refinements[c];

                if ( r > 0 )
                { // refinement (subdivision)
                    LOG_DEBUG ( 2, "Refining cell " << c << " with refinement " << r );

                    const RefinementTree* tree = current_cell.cell_type ( ).refinement_tree ( r - 1 );

                    std::vector<int> local_sub_cell_numbers;
                    const std::vector<Id> children_cells = compute_refined_cells ( current_cell, tree, local_sub_cell_numbers );
                    new_cells.insert ( new_cells.end ( ), children_cells.begin ( ), children_cells.end ( ) );
                    sub_cell_numbers.insert ( sub_cell_numbers.end ( ), local_sub_cell_numbers.begin ( ), local_sub_cell_numbers.end ( ) );

                    // set parent data
                    for ( size_t rc = 0; rc != children_cells.size ( ); ++rc )
                    {
                        parent_history_indices.push_back ( current_history_index );
                        parent_cell_indices.push_back ( c );
                    }

                }
                else
                { // keep cell
                    const Id cell_id = current_cell.id ( );
                    new_cells.push_back ( cell_id );

                    // no parent
                    parent_history_indices.push_back ( -1 );
                    parent_cell_indices.push_back ( -1 );
                    sub_cell_numbers.push_back ( 0 );
                }
            }

            // create history with only current mesh
            std::vector<ConstMeshPtr> refined_history ( 1, ConstMeshPtr ( this ) );

            // create refined mesh
            MeshPtr refined_mesh = MeshPtr ( new RefinedMeshDbView ( td,
                                                                     gd,
                                                                     db_,
                                                                     new_cells,
                                                                     refined_history,
                                                                     parent_history_indices,
                                                                     parent_cell_indices,
                                                                     sub_cell_numbers,
                                                                     period_ ) );
            return refined_mesh;
        }

        MeshPtr MeshDbView::refine_uniform_seq ( int num_ref_steps ) const
        {
            assert ( num_ref_steps > 0 );
            MeshPtr tmp_mesh = this->refine ( );
            for ( int r = 1; r < num_ref_steps; ++r )
            {
                tmp_mesh = tmp_mesh->refine ( );
            }
            return tmp_mesh;
        }

        MeshPtr MeshDbView::refine_uniform_seq ( ) const
        {
            return this->refine_uniform_seq ( 1 );
        }

        bool MeshDbView::is_boundary_facet ( EntityNumber facet_index ) const
        {
            LOG_DEBUG ( 3, "Is facet (index) " << facet_index << " a boundary facet?" );
            const TDim facet_tdim = tdim ( ) - 1;
            assert ( facet_index >= 0 );
            assert ( facet_index < num_entities ( facet_tdim ) );

            // get facet entity
            const Entity facet_entity = get_entity ( facet_tdim, facet_index );

            // get incident cells of facet
            const std::vector<Entity> incident_cells ( begin_incident ( facet_entity, tdim ( ) ), end_incident ( facet_entity, tdim ( ) ) );
            LOG_DEBUG ( 3, "Incident cells: " << string_from_range ( incident_cells.begin ( ), incident_cells.end ( ) ) );

            // if only one cell is connected to the facet, it is boundary,
            // else it is none boundary
            if ( incident_cells.size ( ) == 1 )
            {
                LOG_DEBUG ( 3, "Facet with index " << facet_index << " is a boundary facet." );
            }
            else
            {
                LOG_DEBUG ( 3, "Facet with index " << facet_index << " is NOT a boundary facet." );
            }
            return (incident_cells.size ( ) == 1 );
        }

        //        namespace

        std::vector<MaterialNumber> get_cell_facet_materials ( const MeshDatabasePtr& db, const Entity& cell )
        {
            assert ( cell.tdim ( ) == db->tdim ( ) );
            const TDim facet_tdim = db->tdim ( ) - 1;
            const CellType& cell_type = cell.cell_type ( );
            const int num_facets = cell_type.num_regular_entities ( facet_tdim );

            const std::vector<Id> cell_vertices ( cell.begin_vertex_ids ( ), cell.end_vertex_ids ( ) );

            std::vector<MaterialNumber> material_numbers;

            for ( int f = 0; f < num_facets; ++f )
            {

                // get vertices of facet f
                const std::vector<Id> facet_vertices =
                        cell_type.vertices_of_entity ( facet_tdim, f, cell_vertices );

                Id facet_id;
                if ( db->has_entity ( facet_tdim, facet_vertices, facet_id ) )
                { // does facet exist in database?
                    const MaterialNumber material = db->get_material_number ( facet_tdim, facet_id );

                    if ( material >= 0 )
                    { // does it have non-negative material id?
                        if ( material_numbers.empty ( ) )
                        { // create vector only if a valid material id was found
                            material_numbers.resize ( num_facets, -1 );
                        }
                        material_numbers[f] = material;
                    }
                }
            }

            return material_numbers;
        }

        void inherit_facet_material_numbers ( const MeshDatabasePtr& db,
                                              const CellType& parent_cell_type,
                                              const std::vector<MaterialNumber>& parent_facet_materials,
                                              int sub_cell_number,
                                              const std::vector<Id>& sub_cell_vertices,
                                              const RefinementTree* tree )
        {
            assert ( tree != 0 );
            const TDim cell_tdim = db->tdim ( );
            const TDim facet_tdim = db->tdim ( ) - 1;
            assert ( facet_tdim >= 0 );

            // cell type needed to extract vertices of facets
            const CellType& cell_type = CellType::get_instance ( cell_tdim, sub_cell_vertices.size ( ) );
            const int num_facets = cell_type.num_regular_entities ( facet_tdim );

            // local facet numbers of the sub-cell in the parent cell type
            const std::vector<int> facets_in_parent = parent_cell_type.sub_entities_of_cell ( facet_tdim, sub_cell_number );
            assert ( static_cast < int > ( facets_in_parent.size ( ) ) == num_facets );

            for ( int f = 0; f < num_facets; ++f )
            {

                // get number of f:th parent facet from the refinement tree
                const int facet_number = facets_in_parent[f];
                const int parent_facet_number =
                        parent_cell_type.regular_parent ( facet_tdim, facet_number );

                if ( parent_facet_number == -1 )
                {
                    // interior facet
                    continue;
                }

                const MaterialNumber parent_material_number = parent_facet_materials[parent_facet_number];

                // Check whether the parent facet material number has been set, and should be inherited
                if ( parent_material_number >= 0 )
                {

                    // get vertices of sub-facet
                    std::vector<Id> sub_facet_vertices = cell_type.vertices_of_entity ( facet_tdim, f, sub_cell_vertices );

                    // add sub-facet to database
                    const Id sub_facet_id = db->add_entity ( facet_tdim, sub_facet_vertices );

                    // set material number for sub-facet
                    db->set_material_number ( facet_tdim, sub_facet_id, parent_material_number );
                }
            }
        }
        //        }

        std::vector<Id> MeshDbView::compute_refined_cells ( const Entity& cell,
                                                            const RefinementTree* tree,
                                                            std::vector<int>& sub_cell_numbers ) const
        {
            const GDim gd = gdim ( );
            const TDim td = tdim ( );

            // call mesh::compute_refined_cells to get coordinates and cell sizes
            std::vector<Coordinate> refined_cell_coordinates;
            std::vector<int> refined_cell_sizes;
            ::hiflow::mesh::compute_refined_cells ( cell,
                                                    *tree,
                                                    refined_cell_coordinates,
                                                    refined_cell_sizes,
                                                    sub_cell_numbers, ref_geom_fun_ );

            // add refined vertices to database
            const EntityCount num_refined_vertices = refined_cell_coordinates.size ( ) / gd;
            std::vector<Id> refined_vertex_ids ( num_refined_vertices );
            for ( int v = 0; v < num_refined_vertices; ++v )
            {
                refined_vertex_ids[v] = db_->add_vertex ( vec2ptr ( refined_cell_coordinates ) + gd * v );
            }

            // add refined cells to database
            const EntityCount num_children_cells = refined_cell_sizes.size ( );
            assert ( num_children_cells > 1 );

            // get material numbers on cell facets
            std::vector<MaterialNumber> parent_facet_materials;

            if ( td > 1 )
            {
                parent_facet_materials = get_cell_facet_materials ( db_, cell );
            }

            // vector for new cell id:s
            std::vector<Id> new_cells ( num_children_cells );

            int offset = 0;
            for ( int rc = 0; rc < num_children_cells; ++rc )
            {

                // get vertex id:s
                const int cell_size = refined_cell_sizes[rc];
                const std::vector<Id> refined_cell_vertices ( refined_vertex_ids.begin ( ) + offset,
                                                              refined_vertex_ids.begin ( ) + offset + cell_size );
                offset += cell_size;

                LOG_DEBUG ( 2, "Refined cell " << rc << " = "
                            << string_from_range ( refined_cell_vertices.begin ( ), refined_cell_vertices.end ( ) ) );

                // add to database and return vector
                const Id cell_id = db_->add_entity ( td, refined_cell_vertices );
                new_cells[rc] = cell_id;

                // inherit cell material number from parent
                db_->set_material_number ( td, cell_id, get_material_number ( td, cell.index ( ) ) );

                // inherit facet material numbers from parent facets (see helper function above)
                if ( !parent_facet_materials.empty ( ) )
                {
                    inherit_facet_material_numbers ( db_,
                                                     cell.cell_type ( ),
                                                     parent_facet_materials,
                                                     sub_cell_numbers[rc],
                                                     refined_cell_vertices,
                                                     tree );
                }
            }
            return new_cells;
        }

        void MeshDbView::replace_vertex ( const std::vector<Coordinate>& destination, Id vertex_index )
        {
            assert ( !destination.empty ( ) );
            assert ( static_cast < int > ( destination.size ( ) ) == gdim ( ) );
            // the local vertex_index has to be mapped to the global vertex_id with the help of entities
            db_->replace_vertex ( vec2ptr ( destination ), entities_[0][vertex_index] );
        }

        void MeshDbView::move_vertex ( const std::vector<Coordinate>& displacement, Id vertex_index )
        {
            assert ( !displacement.empty ( ) );
            assert ( static_cast < int > ( displacement.size ( ) ) == gdim ( ) );
            // the local vertex_index has to be mapped to the global vertex_id with the help of entities
            db_->move_vertex ( vec2ptr ( displacement ), entities_[0][vertex_index] );
        }

        void MeshDbView::move_vertices ( const std::vector<Coordinate>& displacements )
        {
            assert ( !displacements.empty ( ) );
            const EntityCount num_verts = displacements.size ( ) / gdim ( );
            assert ( static_cast < int > ( displacements.size ( ) ) == num_verts * gdim ( ) );

            for ( int v = 0; v < num_verts; ++v )
            {
                const std::vector<Coordinate> pt = std::vector<Coordinate>( displacements.begin ( ) + gdim ( ) * v,
                        displacements.begin ( ) + gdim ( ) * ( v + 1 ) );
                // the local vertex_index has to be mapped to the global vertex_id with the help of entities
                db_->move_vertex ( vec2ptr ( pt ), entities_[0][v] );
            }
        }

        MeshPtr MeshDbView::extract_boundary_mesh ( ) const
        {
            const TDim mesh_tdim = tdim ( );

            // does not work for 1D meshes at the moment
            assert ( mesh_tdim >= 2 );

            const EntityIterator facet_begin = begin ( mesh_tdim - 1 );
            const EntityIterator facet_end = end ( mesh_tdim - 1 );

            LOG_DEBUG ( 3, "Facets to check if it is boundary:\n " << string_from_range ( facet_begin, facet_end ) );

            // collect boundary facets
            std::vector<Id> boundary_facets; // id:s of boundary facets (used to create boundary mesh)
            std::vector<EntityNumber> mesh_facet_indices; // indices of boundary facets in current mesh (used to map from bdy mesh to current mesh)

            for ( EntityIterator it = facet_begin; it != facet_end; ++it )
            {
                if ( is_boundary_facet ( it->index ( ) ) )
                {
                    boundary_facets.push_back ( it->id ( ) );
                    mesh_facet_indices.push_back ( it->index ( ) );
                }
            }

            LOG_DEBUG ( 3, "Facets that are boundary:\n " << string_from_range ( boundary_facets.begin ( ), boundary_facets.end ( ) ) );
            if ( !boundary_facets.empty ( ) )
            {
                // compute permutation that sorts boundary_facets (as will be done in the constructor)
                std::vector<int> boundary_facet_permutation;
                compute_sorting_permutation ( boundary_facets, boundary_facet_permutation );

                // permute mesh_facet_indices according to this permutation
                std::vector<EntityNumber> sorted_mesh_facet_indices ( boundary_facets.size ( ) );
                permute_vector ( boundary_facet_permutation, mesh_facet_indices, sorted_mesh_facet_indices );

                // create attribute to map from boundary mesh to current mesh
                AttributePtr mesh_facets ( new IntAttribute ( sorted_mesh_facet_indices ) );

                MeshPtr bdy_mesh ( new MeshDbView ( mesh_tdim - 1, gdim ( ), db_, boundary_facets ) );
                bdy_mesh->add_attribute ( "_mesh_facet_index_", mesh_tdim - 1, mesh_facets );

                return bdy_mesh;
            }
            else
            {
                std::cerr << "Boundary mesh is empty!!!\n";
                return MeshPtr ( 0 );
            }
        }

        GeometricSearchPtr MeshDbView::get_geometric_search ( )
        {
            if ( !geometric_search_ )
            {
                geometric_search_.reset ( new GridGeometricSearch ( this ) );
            }
            return geometric_search_;
        }

        void MeshDbView::reset_geometric_search ( )
        {
            geometric_search_.reset ( new GridGeometricSearch ( this ) );
        }

        void MeshDbView::set_refined_geometry_function ( RefinedGeometryFunction f )
        {
            if ( f )
            {
                ref_geom_fun_ = f;
            }
            else
            {
                ref_geom_fun_ = &standard_refined_geometry_fun;
            }
        }

        void MeshDbView::set_recursive_refined_geometry_function ( RecursiveRefinedGeometryFunction f )
        {
            if ( f )
            {
                rec_ref_geom_fun_ = f;
            }
            else
            {
                rec_ref_geom_fun_ = &recursive_refined_geometry_fun;
            }
        }

        void MeshDbView::convert_connectivity_to_indices ( TDim d, Connectivity& connectivity ) const
        {
            typedef Connectivity::ConnectionIterator CIterator;
            for ( int i = 0; i < connectivity.num_entities ( ); ++i )
            {
                for ( CIterator it = connectivity.begin ( i ); it != connectivity.end ( i ); ++it )
                {
                    EntityNumber local_index;
                    const bool entity_found = find_entity ( d, *it, &local_index );
                    assert ( entity_found );

                    // replace Id with local index
                    *it = local_index;
                }
            }
        }

        int MeshDbView::num_global_cells ( const MPI_Comm& comm ) const
        {
            int num_cells = 0;
            int num_local_cells = this->num_local_cells ( );
            if ( this->tdim ( ) == 2 )
            {
                MPI_Allreduce ( &num_local_cells, &num_cells, 1, MPI_INT, MPI_SUM, comm );
            }
            else if ( this->tdim ( ) == 3 )
            {
                MPI_Allreduce ( &num_local_cells, &num_cells, 1, MPI_INT, MPI_SUM, comm );
            }
            return num_cells;
        }

        int MeshDbView::num_local_cells ( ) const
        {
            return this->num_entities ( this->tdim ( ) ) - this->num_ghost_cells ( );
        }

        int MeshDbView::num_ghost_cells ( ) const
        {
            int num_ghost = 0;
            for ( EntityNumber jc = 0; jc < this->num_entities ( this->tdim ( ) ); ++jc )
            {
                if ( !this->cell_is_local ( jc ) )
                {
                    num_ghost++;
                }
            }
            return num_ghost;
        }

        bool MeshDbView::cell_is_local ( EntityNumber index ) const
        {
            assert ( index >= 0 );
            assert ( index < this->num_entities ( this->tdim ( ) ) );

            if ( !this->has_attribute ( "_remote_index_", this->tdim ( ) ) )
            {
                return true;
            }

            int remote_index = -1;
            this->get_attribute_value ( "_remote_index_", this->tdim ( ), index, &remote_index );
            if ( remote_index >= 0 )
            {
                return false;
            }
            return true;
        }

        void MeshDbView::save ( std::string filename, const MPI_Comm& comm ) const
        {
#ifdef WITH_HDF5
            Mesh::save ( filename, comm );

            H5FilePtr file_ptr ( new H5File ( filename, "w", comm ) );

            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "MeshDbView";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "w" ) );

            for ( int dim = 0; dim < tdim ( ) + 1; ++dim )
            {
                std::stringstream entity_name;
                entity_name << "entities" << dim;
                write_array_parallel<Id>( group_ptr, entity_name.str ( ), entities_[dim].data ( ), comm );
            }

            for ( int i = 0; i < ( tdim ( ) + 1 ) * ( tdim ( ) + 1 ); ++i )
            {
                Connectivity& conny = incidence_connections_[i];
                //convert connectivity to vector
                std::vector<int> conn_as_vec ( 0 );
                std::vector<int> offsets ( conny.num_entities ( ) );
                for ( int ent = 0; ent < conny.num_entities ( ); ++ent )
                {
                    int num_conns = 0;
                    for ( ConnectionIterator inc = conny.begin ( ent ); inc != conny.end ( ent ); ++inc )
                    {
                        conn_as_vec.push_back ( *inc );
                        ++num_conns;
                    }
                    offsets[ent] = num_conns;
                }

                std::stringstream conn_data_name;
                conn_data_name << "incidence_connections" << i;
                write_array_parallel<int>( group_ptr, conn_data_name.str ( ), conn_as_vec, comm );
                std::stringstream offset_data_name;
                offset_data_name << "incidence_connections_offsets" << i;
                write_array_parallel<int>( group_ptr, offset_data_name.str ( ), offsets, comm );
            }
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        void MeshDbView::load ( std::string filename, const MPI_Comm& comm )
        {
#ifdef WITH_HDF5
            Mesh::load ( filename, comm );

            H5FilePtr file_ptr ( new H5File ( filename, "r", comm ) );

            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "MeshDbView";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "r" ) );

            for ( int dim = 0; dim < tdim ( ) + 1; ++dim )
            {
                std::stringstream entity_name;
                entity_name << "entities" << dim;
                read_array_parallel<Id>( group_ptr, entity_name.str ( ), entities_[dim].data ( ), comm );
            }

            for ( int i = 0; i < ( tdim ( ) + 1 ) * ( tdim ( ) + 1 ); ++i )
            {
                std::vector<int> conn_as_vec ( 0 );
                std::vector<int> offsets ( 0 );
                std::stringstream conn_data_name;
                conn_data_name << "incidence_connections" << i;
                read_array_parallel<int>( group_ptr, conn_data_name.str ( ), conn_as_vec, comm );
                std::stringstream offset_data_name;
                offset_data_name << "incidence_connections_offsets" << i;
                read_array_parallel<int>( group_ptr, offset_data_name.str ( ), offsets, comm );
                int curr_pos = 0;
                int next_pos = 0;
                //assert connectivities are empty
                incidence_connections_[i].clear ( );
                for ( int j = 0; j < offsets.size ( ); ++j )
                {
                    next_pos += offsets[j];
                    incidence_connections_[i].add_connections ( std::vector<int>( conn_as_vec.begin ( ) + curr_pos, conn_as_vec.begin ( ) + next_pos ) );
                    curr_pos = next_pos;
                }
            }
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        void MeshDbView::copy_from ( const MeshPtr mesh )
        {
            boost::intrusive_ptr<MeshDbView> mesh_dbview = boost::static_pointer_cast<MeshDbView> ( mesh );

            const int tdim = mesh->tdim ( );
            this->entities_.reset ( new SortedArray<Id>[tdim + 1] );
            this->incidence_connections_ = new Connectivity[( tdim + 1 )*( tdim + 1 )];

            // copy some pointers
            this->geometric_search_ = mesh_dbview->geometric_search_;
            this->ref_geom_fun_ = mesh_dbview->ref_geom_fun_;
            this->rec_ref_geom_fun_ = mesh_dbview->rec_ref_geom_fun_;
            this->db_ = mesh_dbview->db_;

            // copy entities
            for ( int dim = 0; dim < tdim + 1; ++dim )
            {
                this->entities_[dim].data ( ) = mesh_dbview->entities_[dim].data ( );
            }

            // copy incidence connections
            for ( int i = 0; i < ( tdim + 1 ) * ( tdim + 1 ); ++i )
            {
                Connectivity& conny = mesh_dbview->incidence_connections_[i];
                //convert connectivity to vector
                std::vector<int> conn_as_vec ( 0 );
                std::vector<int> offsets ( conny.num_entities ( ) );
                for ( int ent = 0; ent < conny.num_entities ( ); ++ent )
                {
                    int num_conns = 0;
                    for ( ConnectionIterator inc = conny.begin ( ent ); inc != conny.end ( ent ); ++inc )
                    {
                        conn_as_vec.push_back ( *inc );
                        ++num_conns;
                    }
                    offsets[ent] = num_conns;
                }

                int curr_pos = 0;
                int next_pos = 0;
                for ( int j = 0; j < offsets.size ( ); ++j )
                {
                    next_pos += offsets[j];
                    this->incidence_connections_[i].add_connections ( std::vector<int>( conn_as_vec.begin ( ) + curr_pos, conn_as_vec.begin ( ) + next_pos ) );
                    curr_pos = next_pos;
                }
            }

            // copy member variables defined in Mesh class
            Mesh::copy_from ( mesh );
        }

        void MeshDbView::deep_copy_from ( const MeshPtr mesh )
        {
            // standard copy
            this->copy_from ( mesh );

            // copy database
            MeshDatabasePtr tmp_db ( new MeshDatabase ( this->tdim ( ), this->gdim ( ) ) );
            boost::intrusive_ptr<MeshDbView> mesh_dbview = boost::static_pointer_cast<MeshDbView> ( mesh );
            tmp_db->deep_copy_from ( mesh_dbview->get_db ( ) );

            // update database pointer
            this->db_ = tmp_db;
        }

        //////////////// MeshDbViewBuilder ////////////////
        typedef MeshDbViewBuilder::VertexHandle VertexHandle;
        typedef MeshDbViewBuilder::EntityHandle EntityHandle;

        MeshDbViewBuilder::MeshDbViewBuilder ( MeshDatabasePtr db, std::vector<MasterSlave> period )
        : MeshBuilder ( db->tdim ( ), db->gdim ( ), period ), db_ ( db )
        {
            assert ( db_ != 0 );
        }

        MeshDbViewBuilder::MeshDbViewBuilder ( TDim tdim, GDim gdim, std::vector<MasterSlave> period )
        : MeshBuilder ( tdim, gdim, period ), db_ ( new MeshDatabase ( tdim, gdim ) )
        {
            assert ( db_ != 0 );
        }

        MeshDbViewBuilder::MeshDbViewBuilder ( TDim tdim, GDim gdim, int flag, std::vector<MasterSlave> period )
        : MeshBuilder ( tdim, gdim, period )
        {
        }

        MeshDbViewBuilder::MeshDbViewBuilder ( const MeshDbView& mesh )
        : MeshBuilder ( mesh.tdim ( ), mesh.gdim ( ), mesh.get_period ( ) ), db_ ( mesh.db_ )
        {
            // add cells from mesh to builder
            cells_ = mesh.get_entities ( tdim ( ) );
        }

#if 0

        MeshDbViewBuilder::MeshDbViewBuilder ( TDim tdim, GDim gdim, DatabaseConnectionPtr connection )
        : MeshBuilder ( tdim, gdim ), db_ ( new MeshDatabase ( tdim, gdim, connection ) )
        {
            assert ( db_ != 0 );
        }
#endif

        VertexHandle MeshDbViewBuilder::add_vertex ( const std::vector<Coordinate>& coordinates )
        {
            assert ( !coordinates.empty ( ) );
            assert ( static_cast < int > ( coordinates.size ( ) ) == gdim ( ) );
            if ( period_.size ( ) > 0 )
                return db_->add_vertex ( vec2ptr ( periodify ( coordinates, gdim ( ), period_ ) ) );
            else
                return db_->add_vertex ( vec2ptr ( coordinates ) );
        }

        std::vector<VertexHandle> MeshDbViewBuilder::add_vertices ( const std::vector<Coordinate>& coordinates )
        {
            const EntityCount num_verts = coordinates.size ( ) / gdim ( );
            assert ( static_cast < int > ( coordinates.size ( ) ) == num_verts * gdim ( ) );

            std::vector<VertexHandle> vertices ( num_verts );

            for ( int v = 0; v < num_verts; ++v )
            {
                const std::vector<Coordinate> pt = std::vector<Coordinate>( coordinates.begin ( ) + gdim ( ) * v,
                        coordinates.begin ( ) + gdim ( ) * ( v + 1 ) );
                if ( period_.size ( ) > 0 )
                    vertices[v] = db_->add_vertex ( vec2ptr ( periodify ( pt, gdim ( ), period_ ) ) );
                else
                    vertices[v] = db_->add_vertex ( vec2ptr ( pt ) );
            }
            return vertices;
        }

        EntityHandle MeshDbViewBuilder::add_entity ( TDim tdim, const std::vector<VertexHandle>& vertices )
        {
            assert ( !vertices.empty ( ) );
            const EntityHandle entity_id = db_->add_entity ( tdim, vertices );
            if ( tdim == this->tdim ( ) )
            {
                cells_.push_back ( entity_id );
            }
            return entity_id;
        }

        std::vector<EntityHandle> MeshDbViewBuilder::add_entities ( TDim tdim,
                                                                    const std::vector<VertexHandle>& vertices,
                                                                    const std::vector<int>& sizes )
        {
            const EntityCount num_entities = sizes.size ( );
            std::vector<EntityHandle> entities ( num_entities );

            int offset = 0;
            for ( int e = 0; e < num_entities; ++e )
            {
                const int entity_size = sizes[e];
                entities[e] = db_->add_entity ( tdim, std::vector<VertexHandle>( vertices.begin ( ) + offset,
                                                vertices.begin ( ) + offset + entity_size ) );
                offset += entity_size;
            }

            if ( tdim == this->tdim ( ) )
            {
                cells_.insert ( cells_.end ( ), entities.begin ( ), entities.end ( ) );
            }

            return entities;
        }

        void MeshDbViewBuilder::set_material_number ( TDim tdim, EntityHandle entity, MaterialNumber material )
        {
            db_->set_material_number ( tdim, entity, material );
        }

        void MeshDbViewBuilder::clear ( )
        {
            cells_.clear ( );
        }

        MeshPtr MeshDbViewBuilder::build ( )
        {
            if ( !cells_.empty ( ) )
            {
                return MeshPtr ( new MeshDbView ( db_->tdim ( ), db_->gdim ( ), db_, cells_, period_ ) );
            }
            else
            {
                return MeshPtr ( 0 );
            }
        }
    }
} // namespace hiflow
