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

/// \author Thomas Gengenbach and Staffan Ronnas

#ifndef HIFLOW_MESH_DB_VIEW_H
#    define HIFLOW_MESH_DB_VIEW_H

#    include <map>
#    include <vector>
#    include <string>

#    include "common/sorted_array.h"

#    include "mesh/mesh.h"
#    include "mesh/mesh_builder.h"
#    include "mesh/types.h"
#    include "mesh/geometric_search.h"

namespace hiflow
{
    namespace mesh
    {
        class Connectivity;
        class MeshDbViewBuilder;

        typedef SharedPtr<MeshDatabase>::Type MeshDatabasePtr;

        // helper functions for compute_refined_cells -> should not be used elsewhere
        /// \brief Get material numbers on facets of a cell in a vector. If no facets have any material numbers, the return vector is empty.
        std::vector<MaterialNumber> get_cell_facet_materials ( const MeshDatabasePtr& db, const Entity& cell );

        /// \brief Transfer the material numbers from parent facets to  exterior child facets, if they are non-negative.
        void inherit_facet_material_numbers ( const MeshDatabasePtr& db,
                                              const CellType& parent_cell_type,
                                              const std::vector<MaterialNumber>& parent_facet_materials,
                                              int sub_cell_number,
                                              const std::vector<Id>& sub_cell_vertices,
                                              const RefinementTree* tree );
        ///
        /// \brief A view of some of the entities of a MeshDatabase.
        /// This class implements the Mesh interface.
        /// \see Mesh.

        class MeshDbView : public Mesh
        {
          public:
            /// \brief Constructs MeshDbView with topological and
            /// geometrical dimensions tdim and gdim respectively, which
            /// is connected to a MeshDatabase object and contains the given cells.
            // TODO add argument cell_is_local, see MeshPXest
            MeshDbView ( TDim tdim, GDim gdim, MeshDatabasePtr db, const std::vector<Id>& cells, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            MeshDbView ( TDim tdim, GDim gdim, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            ~MeshDbView ( );

            virtual EntityIterator begin ( TDim entity_dim ) const;
            virtual EntityIterator end ( TDim entity_dim ) const;

            virtual IncidentEntityIterator begin_incident ( const Entity& entity, TDim entity_dim ) const;
            virtual IncidentEntityIterator end_incident ( const Entity& entity, TDim entity_dim ) const;

            virtual Id get_id ( TDim entity_dim, EntityNumber index ) const;
            virtual std::vector<Id> get_vertex_ids ( TDim entity_dim, EntityNumber index ) const;
            virtual MaterialNumber get_material_number ( TDim entity_dim, EntityNumber index ) const;
            virtual void set_material_number ( TDim entity_dim, EntityNumber index, MaterialNumber material );
            virtual std::vector<Coordinate> get_coordinates ( TDim entity_dim, EntityNumber index ) const;

            virtual Id get_parent_cell_id ( EntityNumber cell_index ) const;
            virtual Entity get_parent_cell ( EntityNumber cell_index ) const;
            virtual bool cell_has_parent ( EntityNumber cell_index ) const;
            virtual std::vector<Id> get_children_cell_ids ( EntityNumber cell_index ) const;
            virtual std::vector<EntityNumber> get_sibling_cell_indices ( EntityNumber cell_index ) const;

            virtual EntityCount num_entities ( TDim entity_dim ) const;
            virtual EntityCount num_incident_entities ( const Entity& entity, TDim entity_dim ) const;

            virtual Entity get_entity ( TDim entity_dim, EntityNumber entity_number ) const;

            virtual bool find_entity ( TDim entity_dim, Id id, EntityNumber* entity_number ) const;
            virtual bool find_vertex ( const std::vector<Coordinate>& coords, EntityNumber* v_index ) const;

            virtual MeshPtr refine ( ) const;
            virtual MeshPtr refine ( const std::string& attribute_name ) const;
            virtual MeshPtr refine ( std::vector<int>& refinements ) const;
            virtual MeshPtr refine_uniform_seq ( int num_ref_steps ) const;
            virtual MeshPtr refine_uniform_seq ( ) const;

            virtual void replace_vertex ( const std::vector<Coordinate>& destination, Id vertex_index );
            virtual void move_vertex ( const std::vector<Coordinate>& displacement, Id vertex_index );
            virtual void move_vertices ( const std::vector<Coordinate>& displacements );

            virtual int num_global_cells ( const MPI_Comm& comm ) const;
            virtual int num_local_cells ( ) const;
            virtual int num_ghost_cells ( ) const;

            virtual bool cell_is_local ( EntityNumber index ) const;

            // TODO: this class should have private Pointer to the boundary mesh
            // if the boundary already has been extracted, it should not be extracted again.
            virtual MeshPtr extract_boundary_mesh ( ) const;

            virtual GeometricSearchPtr get_geometric_search ( );

            virtual void reset_geometric_search ( );

            virtual void set_refined_geometry_function ( RefinedGeometryFunction f );

            virtual void set_recursive_refined_geometry_function ( RecursiveRefinedGeometryFunction f );

            MeshDatabasePtr get_db ( ) const
            {
                return db_;
            }

            void set_db ( MeshDatabasePtr db )
            {
                db_ = db;
            }

            virtual void save ( std::string filename, const MPI_Comm& comm ) const;

            virtual void load ( std::string filename, const MPI_Comm& comm );

            virtual void copy_from ( const MeshPtr mesh );

            virtual void deep_copy_from ( const MeshPtr mesh );

          protected:
            /// \brief Determines if facet is on the boundary
            virtual bool is_boundary_facet ( const EntityNumber facet_index ) const;

            std::vector<Id> compute_refined_cells ( const Entity& cell,
                                                    const RefinementTree* tree,
                                                    std::vector<int>& sub_cell_numbers ) const;

            void compute_entity_vertex_connectivity ( TDim d ) const;

            SortedArray<Id>& get_entities ( TDim d ) const;

            // Map from local number to database id of the contained
            // entities of dimension dim, with 0 <= dim <= tdim()
            mutable ScopedArray< SortedArray<Id> >::Type entities_;

            // Map from local cell number to database id of the incident
            // entities of dimension dim, with 0 < dim < tdim().
            mutable Connectivity* incidence_connections_;

            RefinedGeometryFunction ref_geom_fun_;
            RecursiveRefinedGeometryFunction rec_ref_geom_fun_;

            MeshDatabasePtr db_;

          private:
            friend class MeshDbViewBuilder;

            void initialize_vertices ( ) const;

            Connectivity& get_incidence_connectivity ( TDim d1, TDim d2 ) const;

            int incidence_connection_index ( TDim d1, TDim d2 ) const;

            bool has_computed_connectivity ( TDim d1, TDim d2 ) const;

            void convert_connectivity_to_indices ( TDim d, Connectivity& connectivity ) const;

            GeometricSearchPtr geometric_search_;
        };

        /// \brief A MeshBuilder that creates a MeshDbView
        /// \see MeshBuilder
        /// \see MeshDbView

        class MeshDbViewBuilder : public MeshBuilder
        {
          public:
            typedef MeshBuilder::VertexHandle VertexHandle;
            typedef MeshBuilder::EntityHandle EntityHandle;

            /// \brief Constructor to use when building a Mesh with an existing MeshDatabase
            explicit MeshDbViewBuilder ( MeshDatabasePtr db, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Constructor to initialize builder with existing MeshDbView mesh
            explicit MeshDbViewBuilder ( const MeshDbView& mesh );

            /// \brief Constructor to use when building a Mesh without an existing MeshDatabase
            MeshDbViewBuilder ( TDim tdim, GDim gdim, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Dummy constructor, only use in combination with p4est mesh
            MeshDbViewBuilder ( TDim tdim, GDim gdim, int flag, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            virtual VertexHandle add_vertex ( const std::vector<Coordinate>& coordinates );
            virtual std::vector<VertexHandle> add_vertices ( const std::vector<Coordinate>& coordinates );

            virtual EntityHandle add_entity ( TDim tdim, const std::vector<VertexHandle>& vertices );
            virtual std::vector<EntityHandle> add_entities ( TDim tdim,
                                                             const std::vector<VertexHandle>& vertices,
                                                             const std::vector<int>& sizes );

            virtual void set_material_number ( TDim tdim, EntityHandle entity, MaterialNumber material );

            virtual void clear ( );

            virtual MeshPtr build ( );

            MeshDatabasePtr get_db ( ) const
            {
                return db_;
            }

          protected:
            MeshDatabasePtr db_;
            std::vector<EntityHandle> cells_;
        };
    }
} // namespace hiflow

#endif
