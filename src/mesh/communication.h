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

#ifndef HIFLOW_MESH_COMMUNICATION_H
#    define HIFLOW_MESH_COMMUNICATION_H

#    include <map>
#    include <vector>
#    include <tr1/unordered_map>

#    include <mpi.h>

#    include "common/sorted_array.h"
#    include "mesh/entity.h"
#    include "mesh/partitioning.h"
#    include "mesh/types.h"

namespace hiflow
{
    namespace mesh
    {
        class Mesh;
        class MeshBuilder;

        void convert_sizes_to_offsets ( const std::vector<int>& sizes,
                                        std::vector<int>& offsets );

        // \brief Information about vertex on another sub-domain

        struct SharedVertex
        {
            SubDomainId sub_domain;
            Id remote_vertex_id;

            bool operator== ( const SharedVertex& v ) const
            {
                return sub_domain == v.sub_domain
                        && remote_vertex_id == v.remote_vertex_id;
            }
        };

        /// \brief Information about shared vertices
        /// Mapping based on vertex id:s.

        class SharedVertexTable
        {

            struct SharedVertexCmp
            {

                bool operator() ( const SharedVertex& v1, const SharedVertex& v2 ) const
                {
                    return v1.sub_domain < v2.sub_domain;
                }
            };

          public:

            typedef SortedArray<SharedVertex, SharedVertexCmp> SharedVertexArray;
            typedef std::tr1::unordered_map< Id, SharedVertexArray > Table;
            typedef std::tr1::unordered_map< Id, SharedVertexArray >::const_iterator const_iterator;
            typedef SortedArray<SubDomainId>::const_iterator NeighborIterator;

            void add_shared_vertex ( Id local_vertex, SharedVertex shared_vertex );
            int num_shared_vertices ( Id local_vertex ) const;
            const SharedVertex& shared_vertex ( Id local_vertex, int i ) const;
            const_iterator begin ( ) const;
            const_iterator end ( ) const;

            NeighborIterator begin_neighbor_domains ( ) const;
            NeighborIterator end_neighbor_domains ( ) const;
            int num_neighbor_domains ( ) const;

            // check if vertex exists in table
            bool has_vertex ( Id local_vertex ) const;

          private:
            Table shared_verts_;
            SortedArray<SubDomainId> neighbor_domains_;
        };

        std::ostream& operator<< ( std::ostream& os, const SharedVertexTable& shared_vertex_table );

        /// \brief Structure for exchange of entity data

        struct EntityPackage
        {
            TDim tdim;
            GDim gdim;
            std::vector<Coordinate> coords;

            // size = num_entities + 1
            std::vector<int> offsets;

            // size = \sum entity_size
            std::vector<int> connections;

            // size = num_entities
            std::vector<MaterialNumber> material_numbers;
        };

        /// \brief Packs all entities of a given dimension in a mesh.
        template<typename EntityIteratorType>
        void pack_entities ( const EntityIteratorType& begin, const EntityIteratorType& end,
                             EntityPackage* entities );

        /// \brief Helper function for packing entire mesh
        void pack_entities ( const Mesh& mesh, TDim tdim, EntityPackage* entities );

        /// \brief Helper function for packing distributed entities
        void pack_distributed_entities ( const Mesh& mesh,
                                         TDim tdim,
                                         const CellDistribution& distribution,
                                         std::vector<EntityPackage>& entities );

        /// \brief Unpacks the entities in an EntityPackage into a MeshBuilder.
        std::vector<int> unpack_entities ( const EntityPackage& entities, MeshBuilder& builder );

        /// \brief Identifies which entities in the local mesh are ghost entities in neighboring sub-domains
        typedef std::map< SubDomainId, SortedArray<EntityNumber> > GhostEntityMap;
        void find_ghost_cells_entities_to_send ( const Mesh& mesh,
                                                 const SharedVertexTable& shared_bdy_vertices,
                                                 GhostEntityMap& entities_to_send );

        void find_ghost_facets_entities_to_send ( const Mesh& mesh,
                                                  const SharedVertexTable& shared_bdy_vertices,
                                                  GhostEntityMap& entities_to_send );

        /// \brief Update a shared vertex table for all vertices in given
        /// mesh through MPI communication.
        void update_shared_vertex_table ( const Mesh& mesh,
                                          const MPI_Comm& comm,
                                          SharedVertexTable& table );

        /// \brief Communicates ghost cells.
        void communicate_ghost_cells ( const Mesh& mesh,
                                       const MPI_Comm& comm,
                                       const SharedVertexTable& shared_bdy_vertices,
                                       MeshBuilder& builder,
                                       std::vector<SubDomainId>& sub_domains,
                                       std::vector<EntityNumber>& remote_indices );

        /// \brief Communicates facets of ghost cells.
        void communicate_ghost_cell_facets ( const Mesh& mesh,
                                             const MPI_Comm& comm,
                                             const SharedVertexTable& shared_bdy_vertices,
                                             MeshBuilder& builder );

        /// \brief Broadcast a mesh.
        ///
        /// \param mesh  pointer to mesh to broadcast. Must be non-zero on process root only.
        void broadcast_mesh ( const Mesh* mesh, const MPI_Comm& comm, MeshBuilder& builder, int root );

        /// \brief Abstract base class for communicators

        class Communicator
        {
          public:

            virtual ~Communicator ( )
            {
            }

            virtual void communicate ( /*const*/ EntityPackage* sent_entities,
                                       EntityPackage* received_entities ) const = 0;
        };

        /// \brief Extended interface for non-blocking communicators

        class NonBlockingCommunicator : public Communicator
        {
          public:
            virtual void wait ( ) const = 0;
        };

        //////////////// pack_entities() template ////////////////

        template<typename EntityIteratorType>
        void pack_entities ( const EntityIteratorType& begin, const EntityIteratorType& end,
                             EntityPackage* entities )
        {
            const TDim tdim = begin->tdim ( );
            const GDim gdim = begin->gdim ( );
            assert ( tdim >= 0 );
            assert ( entities != 0 );

            entities->tdim = tdim;
            entities->gdim = gdim;

            std::tr1::unordered_map<Id, int> vertex_map;

            for ( EntityIteratorType it = begin; it != end; ++it )
            {
                assert ( it->tdim ( ) == tdim );
                assert ( it->gdim ( ) == gdim );

                entities->offsets.push_back ( entities->connections.size ( ) );
                for ( int v = 0; v < it->num_vertices ( ); ++v )
                {
                    const Id v_id = it->vertex_id ( v );
                    if ( vertex_map.find ( v_id ) == vertex_map.end ( ) )
                    {
                        Entity ent_temp = *it;
                        std::vector<Coordinate> v_coords;
                        ent_temp.get_coordinates ( v_coords, v );
                        const int v_index = entities->coords.size ( ) / gdim;
                        entities->coords.insert ( entities->coords.end ( ), v_coords.begin ( ), v_coords.end ( ) );
                        vertex_map[v_id] = v_index;
                    }
                    entities->connections.push_back ( vertex_map[v_id] );
                }
                entities->material_numbers.push_back ( it->get_material_number ( ) );
            }
            // last offset
            entities->offsets.push_back ( entities->connections.size ( ) );
        }

        template<typename EntityIteratorType>
        void pack_entities ( const EntityIteratorType& begin, int num_ent,
                             EntityPackage* entities )
        {
            const TDim tdim = begin->tdim ( );
            const GDim gdim = begin->gdim ( );
            assert ( tdim >= 0 );
            assert ( entities != 0 );

            entities->tdim = tdim;
            entities->gdim = gdim;

            std::tr1::unordered_map<Id, int> vertex_map;
            int count = 0;
            for ( EntityIteratorType it = begin; count != num_ent; ++it )
            {
                assert ( it->tdim ( ) == tdim );
                assert ( it->gdim ( ) == gdim );

                entities->offsets.push_back ( entities->connections.size ( ) );
                for ( int v = 0; v < it->num_vertices ( ); ++v )
                {
                    const Id v_id = it->vertex_id ( v );
                    if ( vertex_map.find ( v_id ) == vertex_map.end ( ) )
                    {
                        Entity ent_temp = *it;
                        std::vector<Coordinate> v_coords;
                        ent_temp.get_coordinates ( v_coords, v );
                        const int v_index = entities->coords.size ( ) / gdim;
                        entities->coords.insert ( entities->coords.end ( ), v_coords.begin ( ), v_coords.end ( ) );
                        vertex_map[v_id] = v_index;
                    }
                    entities->connections.push_back ( vertex_map[v_id] );
                }
                entities->material_numbers.push_back ( it->get_material_number ( ) );
                count++;
            }
            // last offset
            entities->offsets.push_back ( entities->connections.size ( ) );
        }

    }
} // namespace hiflow

#endif /* _COMMUNICATION_H_ */
