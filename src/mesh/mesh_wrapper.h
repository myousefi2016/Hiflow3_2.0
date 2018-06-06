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

#ifndef HIFLOW_MESH_WRAPPER_H
#    define HIFLOW_MESH_WRAPPER_H

#    include <map>
#    include <vector>
#    include <string>

#    include "common/sorted_array.h"
#    include "mesh.h"
#    include "mesh_builder.h"
#    include "iterator.h"
#    include "types.h"

/// \file mesh_wrapper.h
/// \brief Wrapper class to faciliate individual implementation of abstract base class mesh.
///
/// \author Teresa Beck
///
/// Wrapper class used to faciliate individualized implementations of the abstract base class mesh.
/// A MeshWrapper wraps a MeshPtr by forewarding the respective functions of the mesh behind
/// to its respective mesh implementation. This enables the derivation of new mesh implementations
/// without completely redefineing every single pure virtual mesh function of the base class.
/// Instead, the functions are forewarded to the mesh implementation of the mesh
/// object behind the mesh pointer. Thus, additional functionalities can easily be amended or added
/// in the individual implementation.

namespace hiflow
{
    namespace mesh
    {
        class Connectivity;
        class MeshDbViewBuilder;

        typedef SharedPtr<MeshDatabase>::Type MeshDatabasePtr;

        class MeshWrapper : public Mesh
        {
          public:

            MeshWrapper ( MeshPtr mesh, TDim tdim, GDim gdim ) :
            Mesh ( tdim, gdim ),
            mymesh_ ( mesh )
            {
            }

            ~MeshWrapper ( )
            {
                ;
            }

            virtual EntityIterator begin ( TDim entity_dim ) const
            {
                assert ( entity_dim >= 0 );
                assert ( entity_dim <= tdim ( ) );

                const Entity entity = Entity ( ConstMeshPtr ( this ), entity_dim, 0 );
                return EntityIterator ( entity, 0 );
            }

            virtual EntityIterator end ( TDim entity_dim ) const
            {
                assert ( entity_dim >= 0 );
                assert ( entity_dim <= tdim ( ) );
                const EntityCount end_entity = num_entities ( entity_dim );
                const Entity entity = Entity ( ConstMeshPtr ( this ), entity_dim, end_entity );
                return EntityIterator ( entity, end_entity );
            }

            virtual IncidentEntityIterator begin_incident ( const Entity& entity, TDim entity_dim ) const
            {
                IncidentEntityIterator it = mymesh_->begin_incident ( entity, entity_dim );
                Entity wrap_incident_entity ( ConstMeshPtr ( this ), entity_dim, it->index ( ) );
                return IncidentEntityIterator ( wrap_incident_entity, it.get_incidence_iterator ( ) );
            }

            virtual IncidentEntityIterator end_incident ( const Entity& entity, TDim entity_dim ) const
            {

                IncidentEntityIterator it = mymesh_->end_incident ( entity, entity_dim );
                Entity wrap_incident_entity ( ConstMeshPtr ( this ), entity_dim, -1 );
                return IncidentEntityIterator ( wrap_incident_entity, it.get_incidence_iterator ( ) );
            }

            virtual Id get_id ( TDim entity_dim, EntityNumber index ) const
            {
                return mymesh_->get_id ( entity_dim, index );
            }

            virtual std::vector<Id> get_vertex_ids ( TDim entity_dim, EntityNumber index ) const
            {
                return mymesh_->get_vertex_ids ( entity_dim, index );
            }

            virtual MaterialNumber get_material_number ( TDim entity_dim, EntityNumber index ) const
            {
                return mymesh_->get_material_number ( entity_dim, index );
            }

            virtual void set_material_number ( TDim entity_dim, EntityNumber index, MaterialNumber material )
            {
                mymesh_->set_material_number ( entity_dim, index, material );
            }

            virtual std::vector<Coordinate> get_coordinates ( TDim entity_dim, EntityNumber index ) const
            {
                return mymesh_->get_coordinates ( entity_dim, index );
            }

            virtual Id get_parent_cell_id ( EntityNumber cell_index ) const
            {
                return mymesh_->get_parent_cell_id ( cell_index );
            }

            virtual Entity get_parent_cell ( EntityNumber cell_index ) const
            {
                return mymesh_->get_parent_cell ( cell_index );
            }

            virtual bool cell_has_parent ( EntityNumber cell_index ) const
            {
                return mymesh_->cell_has_parent ( cell_index );
            }

            virtual std::vector<Id> get_children_cell_ids ( EntityNumber cell_index ) const
            {
                return mymesh_->get_children_cell_ids ( cell_index );
            }

            virtual std::vector<EntityNumber> get_sibling_cell_indices ( EntityNumber cell_index ) const
            {
                return mymesh_->get_sibling_cell_indices ( cell_index );
            }

            virtual EntityCount num_entities ( TDim entity_dim ) const
            {
                return mymesh_->num_entities ( entity_dim );
            }

            virtual EntityCount num_incident_entities ( const Entity& entity, TDim entity_dim ) const
            {
                return mymesh_->num_incident_entities ( entity, entity_dim );
            }

            virtual Entity get_entity ( TDim entity_dim, EntityNumber entity_number ) const
            {
                assert ( entity_dim >= 0 );
                assert ( entity_dim <= tdim ( ) );
                assert ( entity_number >= 0 );
                assert ( entity_number < num_entities ( entity_dim ) );
                return Entity ( ConstMeshPtr ( this ), entity_dim, entity_number );
            }

            virtual bool find_entity ( TDim entity_dim, Id id, EntityNumber* entity_number ) const
            {
                return mymesh_->find_entity ( entity_dim, id, entity_number );
            }

            virtual bool find_vertex ( const std::vector<Coordinate>& coords, EntityNumber* v_index ) const
            {
                return mymesh_->find_vertex ( coords, v_index );
            }

            virtual MeshPtr refine ( ) const
            {
                return mymesh_->refine ( );
            }

            virtual MeshPtr refine ( const std::string& attribute_name ) const
            {
                return mymesh_->refine ( attribute_name );
            }

            virtual MeshPtr refine ( std::vector<int>& refinements ) const
            {
                return mymesh_->refine ( refinements );
            }

            virtual MeshPtr refine_uniform_seq ( int num_ref_steps ) const
            {
                return mymesh_->refine_uniform_seq ( num_ref_steps );
            }

            virtual MeshPtr refine_uniform_seq ( ) const
            {
                return mymesh_->refine_uniform_seq ( );
            }

            virtual void replace_vertex ( const std::vector<Coordinate>& destination, Id vertex_index )
            {
                mymesh_->replace_vertex ( destination, vertex_index );
            }

            virtual void move_vertex ( const std::vector<Coordinate>& displacement, Id vertex_index )
            {
                mymesh_->move_vertex ( displacement, vertex_index );
            }

            virtual void move_vertices ( const std::vector<Coordinate>& displacements )
            {
                mymesh_->move_vertices ( displacements );
            }

            virtual MeshPtr extract_boundary_mesh ( ) const
            {
                return mymesh_->extract_boundary_mesh ( );
            }

            virtual GeometricSearchPtr get_geometric_search ( )
            {
                return mymesh_->get_geometric_search ( );
            }

            virtual void reset_geometric_search ( )
            {
                mymesh_->reset_geometric_search ( );
            }

            virtual void set_refined_geometry_function ( RefinedGeometryFunction f )
            {
                mymesh_->set_refined_geometry_function ( f );
            }

            virtual int num_global_cells ( const MPI_Comm& comm ) const
            {
                return mymesh_->num_global_cells ( comm );
            }

            virtual int num_local_cells ( ) const
            {
                return mymesh_->num_local_cells ( );
            }

            virtual int num_ghost_cells ( ) const
            {
                return mymesh_->num_ghost_cells ( );
            }

            virtual bool cell_is_local ( EntityNumber index ) const
            {
                return mymesh_->cell_is_local ( index );
            }

            virtual void copy_from ( const MeshPtr mesh ) const
            {
                return mymesh_->copy_from ( mesh );
            }

            virtual void deep_copy_from ( const MeshPtr mesh ) const
            {
                return mymesh_->deep_copy_from ( mesh );
            }

          protected:

            virtual bool is_boundary_facet ( const EntityNumber facet_index ) const
            {
                return mymesh_->is_boundary_facet ( facet_index );
            }
            RefinedGeometryFunction ref_geom_fun_;

            //   private:
            friend class MeshDbViewBuilder;
            friend class PeriodicMeshDbView;
            MeshPtr mymesh_;

        };

    }
}

#endif
