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

#include "communication.h"

#include <algorithm>
#include <boost/iterator/filter_iterator.hpp>
#include <numeric>
#include <functional>
#include <tr1/unordered_map>
#include <vector>

#include "common/log.h"
#include "common/sorted_array.h"

#include "iterator.h"
#include "mesh.h"
#include "mesh_builder.h"
#include "mpi_communication.h" // needed for ghost cells

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace mesh
    {

        // Helper function

        void convert_sizes_to_offsets ( const std::vector<int>& sizes,
                                        std::vector<int>& offsets )
        {
            offsets.clear ( );
            offsets.resize ( sizes.size ( ) );
            offsets[0] = 0;
            std::partial_sum ( sizes.begin ( ), sizes.end ( ) - 1, offsets.begin ( ) + 1 );
        }

        //////////////// SharedVertexTable implementation ////////////////

        void SharedVertexTable::add_shared_vertex ( Id local_vertex, SharedVertex shared_vertex )
        {
            Table::iterator it = shared_verts_.find ( local_vertex );

            if ( it != shared_verts_.end ( ) )
            {
                // local index exists already
                if ( !it->second.find ( shared_vertex, 0 ) )
                {
                    // shared vertex exists already
                    // (NB: search only based on sub_domain)
                    it->second.insert ( shared_vertex );

                }
                // TODO if shared_vertex has same sub_domain but different remote_vertex_id,
                // it is not added. Check that this is desired behavior.
            }
            else
            {
                // first time we are adding shared_vertex for local_vertex
                shared_verts_.insert ( std::make_pair ( local_vertex,
                                                        SharedVertexArray ( 1, shared_vertex ) ) );
            }

            // add sub_domain to neighbor_domains_ if it does not already exist
            neighbor_domains_.find_insert ( shared_vertex.sub_domain );
        }

        SharedVertexTable::const_iterator SharedVertexTable::begin ( ) const
        {
            return shared_verts_.begin ( );
        }

        SharedVertexTable::const_iterator SharedVertexTable::end ( ) const
        {
            return shared_verts_.end ( );
        }

        int SharedVertexTable::num_shared_vertices ( Id local_vertex ) const
        {
            Table::const_iterator it = shared_verts_.find ( local_vertex );

            if ( it != shared_verts_.end ( ) )
            {
                return it->second.size ( );
            }

            return 0;
        }

        const SharedVertex& SharedVertexTable::shared_vertex ( Id local_vertex, int i ) const
        {
            assert ( i >= 0 );
            assert ( i < num_shared_vertices ( local_vertex ) );
            Table::const_iterator it = shared_verts_.find ( local_vertex );
            return it->second[i];
        }

        bool SharedVertexTable::has_vertex ( Id local_vertex ) const
        {
            Table::const_iterator it = shared_verts_.find ( local_vertex );
            return it != shared_verts_.end ( );
        }

        SharedVertexTable::NeighborIterator SharedVertexTable::begin_neighbor_domains ( ) const
        {
            return neighbor_domains_.begin ( );
        }

        SharedVertexTable::NeighborIterator SharedVertexTable::end_neighbor_domains ( ) const
        {
            return neighbor_domains_.end ( );
        }

        int SharedVertexTable::num_neighbor_domains ( ) const
        {
            return neighbor_domains_.size ( );
        }

        std::ostream& operator<< ( std::ostream& os, const SharedVertexTable& shared_vertex_table )
        {
            for ( SharedVertexTable::const_iterator it = shared_vertex_table.begin ( );
                  it != shared_vertex_table.end ( ); ++it )
            {
                os << "Local vertex id " << it->first << " is shared with \n";
                for ( int sv = 0; sv < shared_vertex_table.num_shared_vertices ( it->first ); ++sv )
                {
                    SharedVertex v = shared_vertex_table.shared_vertex ( it->first, sv );
                    os << "\t (Subdomain " << v.sub_domain
                            << ", Remote vertex id " << v.remote_vertex_id << ")\n";
                }
            }
            return os;
        }

        void pack_entities ( const Mesh& mesh, TDim tdim, EntityPackage* entities )
        {
            assert ( tdim <= mesh.tdim ( ) );
            pack_entities ( mesh.begin ( tdim ), mesh.end ( tdim ), entities );
        }

        void pack_distributed_entities ( const Mesh& mesh,
                                         TDim tdim,
                                         const CellDistribution& distribution,
                                         std::vector<EntityPackage>& entities )
        {
            const GDim gdim = mesh.gdim ( );
            assert ( tdim >= 0 );
            assert ( tdim <= mesh.tdim ( ) );
            assert ( tdim <= gdim );

            const int num_procs = distribution.num_cells.size ( );

            std::vector<int> proc_offsets ( num_procs + 1, 0 );
            for ( int i = 0; i < num_procs; ++i )
            {
                proc_offsets[i + 1] = proc_offsets[i] + distribution.num_cells[i];
            }

            entities.clear ( );
            entities.resize ( num_procs );

            for ( int p = 0; p < num_procs; ++p )
            {
                entities[p].tdim = tdim;
                entities[p].gdim = gdim;
            }

            std::vector<std::tr1::unordered_map<Id, int> > vertex_map ( num_procs );

            for ( int p = 0; p < num_procs; ++p )
            {
                for ( int c = proc_offsets[p]; c < proc_offsets[p + 1]; ++c )
                {
                    entities[p].offsets.push_back ( entities[p].connections.size ( ) );

                    const Entity ent_temp = mesh.get_entity ( tdim, distribution.cell_indices[c] );

                    // get vertex indices, adding vertices if they do not already exist
                    for ( int v = 0; v < ent_temp.num_vertices ( ); ++v )
                    {
                        const Id v_id = ent_temp.vertex_id ( v );
                        if ( vertex_map[p].find ( v_id ) == vertex_map[p].end ( ) )
                        {
                            std::vector<Coordinate> v_coords;
                            ent_temp.get_coordinates ( v_coords, v );
                            const int v_index = entities[p].coords.size ( ) / gdim;
                            entities[p].coords.insert ( entities[p].coords.end ( ), v_coords.begin ( ), v_coords.end ( ) );
                            vertex_map[p][v_id] = v_index;
                        }
                        entities[p].connections.push_back ( vertex_map[p][v_id] );
                    }
                    entities[p].material_numbers.push_back ( ent_temp.get_material_number ( ) );
                }
                // last offset
                entities[p].offsets.push_back ( entities[p].connections.size ( ) );
            }
        }

        std::vector<int> unpack_entities ( const EntityPackage& entities, MeshBuilder& builder )
        {
            assert ( entities.gdim == builder.gdim ( ) );

            // add all vertices to builder and receive corresponding VertexHandle:s
            const std::vector<MeshBuilder::VertexHandle> vertex_ids
                    = builder.add_vertices ( entities.coords );

            const TDim tdim = entities.tdim;
            assert ( tdim <= builder.tdim ( ) );

            // add all entities to builder, after translating according to received VertexHandle:s
            std::vector<int>::const_iterator it = entities.offsets.begin ( );
            std::vector<int>::const_iterator next = it;
            ++next;
            std::vector<int>::const_iterator end = entities.offsets.end ( );

            std::vector<MeshBuilder::VertexHandle> vertex_connections;
            std::vector<int> entity_sizes;
            entity_sizes.reserve ( entities.offsets.size ( ) );

            while ( next != end )
            {
                entity_sizes.push_back ( 0 );
                const int entity_begin = *it;
                const int entity_end = ( next == end ) ? entities.offsets.size ( ) - 1 : *next;
                for ( int k = entity_begin; k < entity_end; ++k )
                {
                    const MeshBuilder::VertexHandle v = entities.connections[k];
                    vertex_connections.push_back ( vertex_ids[v] );
                    ++entity_sizes.back ( );
                }
                ++it;
                ++next;
            }

            const std::vector<MeshBuilder::EntityHandle> entity_handles =
                    builder.add_entities ( tdim, vertex_connections, entity_sizes );

            // add material numbers
            if ( !entities.material_numbers.empty ( ) )
            {
                assert ( entities.material_numbers.size ( ) == entity_handles.size ( ) );
                for ( size_t i = 0; i != entity_handles.size ( ); ++i )
                {
                    LOG_DEBUG ( 3, "Setting material number of entity (" << tdim << ", " << entity_handles[i] << ") to " << entities.material_numbers[i] );
                    builder.set_material_number ( tdim, entity_handles[i], entities.material_numbers[i] );
                }
            }
            return entity_handles;
        }

        void find_ghost_cells_entities_to_send ( const Mesh& mesh,
                                                 const SharedVertexTable& shared_bdy_vertices,
                                                 GhostEntityMap& entities_to_send )
        {

            for ( EntityIterator ent_it = mesh.begin ( mesh.tdim ( ) ); ent_it != mesh.end ( mesh.tdim ( ) ); ++ent_it )
            {
                const Id ent_index = ent_it->index ( );
                for ( VertexIdIterator v_it = ent_it->begin_vertex_ids ( ); v_it != ent_it->end_vertex_ids ( ); ++v_it )
                {
                    const Id v_id = *v_it;
                    for ( int s = 0; s < shared_bdy_vertices.num_shared_vertices ( v_id ); ++s )
                    {
                        const SharedVertex& sv = shared_bdy_vertices.shared_vertex ( v_id, s );
                        const SubDomainId sub_domain = sv.sub_domain;
                        entities_to_send[sub_domain].find_insert ( ent_index );
                    }
                }
            }
        }

        void find_ghost_facets_entities_to_send ( const Mesh& mesh,
                                                  const SharedVertexTable& shared_bdy_vertices,
                                                  GhostEntityMap& entities_to_send )
        {

            for ( EntityIterator cell_it = mesh.begin ( mesh.tdim ( ) ), end_cell_it = mesh.end ( mesh.tdim ( ) );
                  cell_it != end_cell_it;
                  ++cell_it )
            {

                for ( VertexIdIterator v_it = cell_it->begin_vertex_ids ( ), end_v_it = cell_it->end_vertex_ids ( );
                      v_it != end_v_it;
                      ++v_it )
                {

                    const Id v_id = *v_it;
                    for ( int s = 0; s < shared_bdy_vertices.num_shared_vertices ( v_id ); ++s )
                    {

                        const SharedVertex& sv = shared_bdy_vertices.shared_vertex ( v_id, s );
                        const SubDomainId sub_domain = sv.sub_domain;

                        for ( IncidentEntityIterator facet_it = cell_it->begin_incident ( mesh.tdim ( ) - 1 ), end_facet_it = cell_it->end_incident ( mesh.tdim ( ) - 1 );
                              facet_it != end_facet_it;
                              ++facet_it )
                        {

                            const Id ent_index = facet_it->index ( );
                            entities_to_send[sub_domain].find_insert ( ent_index );

                        }

                    }

                }

            }

        }

        void update_shared_vertex_table ( const Mesh& mesh, const MPI_Comm& comm,
                                          SharedVertexTable& table )
        {
            // Precondition: gdim of meshes on all procs is equal

            const GDim gdim = mesh.gdim ( );

            // assume SubDomainId == rank in comm
            SubDomainId local_domain = -1;
            MPI_Comm_rank ( comm, &local_domain );
            int num_domains = -1;
            MPI_Comm_size ( comm, &num_domains );

            std::vector<Coordinate> local_coords;
            std::vector<Id> local_ids;

            // reserve memory
            local_coords.reserve ( mesh.num_entities ( 0 ) );
            local_ids.reserve ( mesh.num_entities ( 0 ) );

            // 1. find vertices to be sent -- all vertices in Mesh not already in table (local)
            for ( EntityIterator v_it = mesh.begin ( 0 ); v_it != mesh.end ( 0 ); ++v_it )
            {
                const Id v_id = v_it->id ( );
                if ( !table.has_vertex ( v_id ) )
                {
                    std::vector<Coordinate> v_coords;
                    v_it->get_coordinates ( v_coords );
                    local_coords.insert ( local_coords.end ( ), v_coords.begin ( ), v_coords.end ( ) );
                    local_ids.push_back ( v_id );
                }
            }

            assert ( local_coords.size ( ) == local_ids.size ( ) * gdim );

            // 2. exchange number of vertices to be sent on all procs (Allgather)
            std::vector<int> num_vertices ( num_domains );
            int num_local_vertices = local_ids.size ( );
            MPI_Allgather ( &num_local_vertices, 1, MPI_INT,
                            vec2ptr ( num_vertices ), 1, MPI_INT, comm );

            // 3. exchange coordinates (Allgatherv)
            std::vector<int> coord_sizes ( num_domains );
            // == num_vertices * gdim
            std::transform ( num_vertices.begin ( ), num_vertices.end ( ), coord_sizes.begin ( ),
                             std::bind1st ( std::multiplies<Coordinate>( ), gdim ) );

            std::vector<int> coord_offsets;
            convert_sizes_to_offsets ( coord_sizes, coord_offsets );
            const int sum_coord_sizes = coord_offsets.back ( ) + coord_sizes.back ( );
            std::vector<Coordinate> coords ( sum_coord_sizes );

            MPI_Allgatherv ( vec2ptr ( local_coords ), coord_sizes[local_domain], MPI_DOUBLE,
                             vec2ptr ( coords ), vec2ptr ( coord_sizes ), vec2ptr ( coord_offsets ), MPI_DOUBLE, comm );

            // 4. lookup all received vertices (local)
            std::vector<int> vertex_offsets;
            convert_sizes_to_offsets ( num_vertices, vertex_offsets );
            const int total_num_vertices = vertex_offsets.back ( ) + num_vertices.back ( );
            assert ( total_num_vertices * gdim == sum_coord_sizes );

            std::vector<Id> local_vertex_ids ( total_num_vertices, -1 );

            int v = 0;
            for ( int p = 0; p < num_domains; ++p )
            {
                // skip local vertices
                if ( p != local_domain )
                {
                    for (; v < vertex_offsets[p] + num_vertices[p]; ++v )
                    {
                        std::vector<Coordinate> pt ( &coords[gdim * v], &coords[gdim * ( v + 1 )] );
                        int v_index;
                        const bool found = mesh.find_vertex ( pt, &v_index );
                        if ( found )
                        {
                            local_vertex_ids[v] = mesh.get_id ( 0, v_index );
                        }
                    }
                }
                else
                {
                    // advance v, skipping local vertices
                    v += num_vertices[p];
                }

            }
            assert ( v == total_num_vertices );

            // 5. send/recv id:s of matched vertices to/from all other procs, -1 if not found (Alltoallv)
            std::vector<Id> recv_vertex_ids ( num_vertices[local_domain] * num_domains );
            std::vector<int> num_recv_vertex_ids ( num_domains, num_vertices[local_domain] );
            std::vector<int> recv_vertex_ids_offsets;
            convert_sizes_to_offsets ( num_recv_vertex_ids, recv_vertex_ids_offsets );

            MPI_Alltoallv ( vec2ptr ( local_vertex_ids ), vec2ptr ( num_vertices ),
                            vec2ptr ( vertex_offsets ), MPI_INT,
                            vec2ptr ( recv_vertex_ids ), vec2ptr ( num_recv_vertex_ids ),
                            vec2ptr ( recv_vertex_ids_offsets ), MPI_INT, comm );

            // 6. update table with received id:s (local)
            for ( int p = 0; p < num_domains; ++p )
            {
                // skip local vertices
                if ( p != local_domain )
                {
                    const int offset = recv_vertex_ids_offsets[p];
                    for ( int v = 0; v < num_vertices[local_domain]; ++v )
                    {
                        Id remote_vertex_id = recv_vertex_ids[offset + v];
                        if ( remote_vertex_id >= 0 )
                        {
                            const SharedVertex sv = { p, remote_vertex_id };
                            table.add_shared_vertex ( local_ids[v], sv );
                        }
                    }
                }
            }
        }

        /// \brief Iterate only over the entities with the provided
        /// indices.

        class SelectedIndicesFilter
        {
          public:

            SelectedIndicesFilter ( const SortedArray<EntityNumber>& indices )
            : indices_ ( indices )
            {
            }

            bool operator() ( const Entity& entity ) const
            {
                return indices_.find ( entity.index ( ), 0 );
            }

          private:
            const SortedArray<EntityNumber>& indices_;
        };

        void communicate_ghost_cells ( const Mesh& mesh,
                                       const MPI_Comm& comm,
                                       const SharedVertexTable& shared_bdy_vertices,
                                       MeshBuilder& builder,
                                       std::vector<SubDomainId>& sub_domains,
                                       std::vector<EntityNumber>& remote_indices )
        {
            typedef boost::filter_iterator<SelectedIndicesFilter, EntityIterator> SelectedEntitiesIterator;
            const int CELL_INDICES_TAG = 10;
            const TDim tdim = mesh.tdim ( );

            // compute ghost cells to send
            GhostEntityMap ghost_cells_map;
            find_ghost_cells_entities_to_send ( mesh, shared_bdy_vertices, ghost_cells_map );

            const int num_neighbor_domains = shared_bdy_vertices.num_neighbor_domains ( );
            int local_domain = -1;
            MPI_Comm_rank ( comm, &local_domain );

            ScopedArray<EntityPackage>::Type sent_entity_packages ( new EntityPackage[num_neighbor_domains] );
            ScopedArray<EntityPackage>::Type recv_entity_packages ( new EntityPackage[num_neighbor_domains] );
            ScopedArray<MpiNonBlockingPointToPoint*>::Type send_communicators (
                                                                                new MpiNonBlockingPointToPoint*[num_neighbor_domains] );
            ScopedArray<MpiNonBlockingPointToPoint*>::Type recv_communicators (
                                                                                new MpiNonBlockingPointToPoint*[num_neighbor_domains] );

            int k = 0;

            // pack entities

            // local indices of sent ghost cells
            std::vector< std::vector<EntityNumber> > sent_ghost_cell_indices ( num_neighbor_domains );

            for ( SharedVertexTable::NeighborIterator n_it = shared_bdy_vertices.begin_neighbor_domains ( );
                  n_it != shared_bdy_vertices.end_neighbor_domains ( ); ++n_it )
            {
                const SortedArray<EntityNumber>& ghost_cells_on_k = ghost_cells_map[*n_it];
                const SelectedEntitiesIterator begin ( SelectedIndicesFilter ( ghost_cells_on_k ),
                                                       mesh.begin ( tdim ), mesh.end ( tdim ) );
                const SelectedEntitiesIterator end ( SelectedIndicesFilter ( ghost_cells_on_k ),
                                                     mesh.end ( tdim ), mesh.end ( tdim ) );
                pack_entities ( begin, end, &sent_entity_packages[k] );

                for ( SelectedEntitiesIterator it = begin; it != end; ++it )
                {
                    sent_ghost_cell_indices[k].push_back ( it->index ( ) );
                }

                ++k;
            }

            std::vector<MPI_Request> cell_indices_send_requests ( num_neighbor_domains );
            k = 0;
            for ( SharedVertexTable::NeighborIterator n_it = shared_bdy_vertices.begin_neighbor_domains ( );
                  n_it != shared_bdy_vertices.end_neighbor_domains ( ); ++n_it )
            {
                // start all sends: must be done before receives
                send_communicators[k] = new MpiNonBlockingPointToPoint ( comm, local_domain, *n_it );
                send_communicators[k]->communicate ( &sent_entity_packages[k], 0 );

                // send also local indices of ghost cells
                MPI_Isend ( vec2ptr ( sent_ghost_cell_indices[k] ),
                            sent_ghost_cell_indices[k].size ( ),
                            MPI_INT, *n_it,
                            CELL_INDICES_TAG, comm, &cell_indices_send_requests[k] );

                ++k;
            }

            k = 0;
            for ( SharedVertexTable::NeighborIterator n_it = shared_bdy_vertices.begin_neighbor_domains ( );
                  n_it != shared_bdy_vertices.end_neighbor_domains ( ); ++n_it )
            {
                // start all receives
                recv_communicators[k] = new MpiNonBlockingPointToPoint ( comm, *n_it, local_domain );
                recv_communicators[k]->communicate ( 0, &recv_entity_packages[k] );

                ++k;
            }

            // treat received packages one after the other
            // TODO: can we avoid the ordering here?
            // TODO: what is gained by using nonblocking communication here?
            std::vector<MPI_Request> cell_indices_recv_requests ( num_neighbor_domains );
            k = 0;
            for ( SharedVertexTable::NeighborIterator n_it = shared_bdy_vertices.begin_neighbor_domains ( );
                  n_it != shared_bdy_vertices.end_neighbor_domains ( ); ++n_it )
            {
                recv_communicators[k]->wait ( );
                unpack_entities ( recv_entity_packages[k], builder );

                // append current sub-domain id to sub_domains vector
                const int num_cells = recv_entity_packages[k].offsets.size ( ) - 1;
                sub_domains.resize ( sub_domains.size ( ) + num_cells, *n_it );

                // receive ghost cell indices
                std::vector<EntityNumber> cell_indices ( num_cells );
                MPI_Irecv ( vec2ptr ( cell_indices ),
                            num_cells, MPI_INT, *n_it, CELL_INDICES_TAG,
                            comm, &cell_indices_recv_requests[k] );
                MPI_Status recv_status;
                MPI_Wait ( &cell_indices_recv_requests[k], &recv_status );
                remote_indices.insert ( remote_indices.end ( ), cell_indices.begin ( ), cell_indices.end ( ) );

                ++k;
            }

            // release send communicators
            for ( k = 0; k < num_neighbor_domains; ++k )
            {
                MPI_Status send_status;
                MPI_Wait ( &cell_indices_send_requests[k], &send_status );

                send_communicators[k]->wait ( );
                delete send_communicators[k];
                delete recv_communicators[k];
            }
            LOG_DEBUG ( 2, "Exit communicate_ghost_cells() on " << local_domain );
        }

        void communicate_ghost_cell_facets ( const Mesh& mesh,
                                             const MPI_Comm& comm,
                                             const SharedVertexTable& shared_bdy_vertices,
                                             MeshBuilder& builder )
        {
            typedef boost::filter_iterator<SelectedIndicesFilter, EntityIterator> SelectedEntitiesIterator;
            const int FACET_INDICES_TAG = 11;
            const TDim facet_tdim = mesh.tdim ( ) - 1;

            // compute ghost cells to send
            GhostEntityMap ghost_facets_map;
            find_ghost_facets_entities_to_send ( mesh, shared_bdy_vertices, ghost_facets_map );

            const int num_neighbor_domains = shared_bdy_vertices.num_neighbor_domains ( );
            int local_domain = -1;
            MPI_Comm_rank ( comm, &local_domain );

            ScopedArray<EntityPackage>::Type sent_entity_packages ( new EntityPackage[num_neighbor_domains] );
            ScopedArray<EntityPackage>::Type recv_entity_packages ( new EntityPackage[num_neighbor_domains] );
            ScopedArray<MpiNonBlockingPointToPoint*>::Type send_communicators (
                                                                                new MpiNonBlockingPointToPoint*[num_neighbor_domains] );
            ScopedArray<MpiNonBlockingPointToPoint*>::Type recv_communicators (
                                                                                new MpiNonBlockingPointToPoint*[num_neighbor_domains] );

            int k = 0;

            // pack entities

            for ( SharedVertexTable::NeighborIterator n_it = shared_bdy_vertices.begin_neighbor_domains ( );
                  n_it != shared_bdy_vertices.end_neighbor_domains ( ); ++n_it )
            {
                const SortedArray<EntityNumber>& ghost_facets_on_k = ghost_facets_map[*n_it];
                const SelectedEntitiesIterator begin ( SelectedIndicesFilter ( ghost_facets_on_k ),
                                                       mesh.begin ( facet_tdim ), mesh.end ( facet_tdim ) );
                const SelectedEntitiesIterator end ( SelectedIndicesFilter ( ghost_facets_on_k ),
                                                     mesh.end ( facet_tdim ), mesh.end ( facet_tdim ) );
                pack_entities ( begin, end, &sent_entity_packages[k] );

                ++k;
            }

            // Send entities

            k = 0;
            for ( SharedVertexTable::NeighborIterator n_it = shared_bdy_vertices.begin_neighbor_domains ( );
                  n_it != shared_bdy_vertices.end_neighbor_domains ( ); ++n_it )
            {
                // start all sends: must be done before receives
                send_communicators[k] = new MpiNonBlockingPointToPoint ( comm, local_domain, *n_it );
                send_communicators[k]->communicate ( &sent_entity_packages[k], 0 );

                ++k;
            }

            // Receive entities

            k = 0;
            for ( SharedVertexTable::NeighborIterator n_it = shared_bdy_vertices.begin_neighbor_domains ( );
                  n_it != shared_bdy_vertices.end_neighbor_domains ( ); ++n_it )
            {
                // start all receives
                recv_communicators[k] = new MpiNonBlockingPointToPoint ( comm, *n_it, local_domain );
                recv_communicators[k]->communicate ( 0, &recv_entity_packages[k] );

                ++k;
            }

            // Unpack entities
            // treat received packages one after the other
            // TODO: can we avoid the ordering here?
            // TODO: what is gained by using nonblocking communication here?
            k = 0;
            for ( SharedVertexTable::NeighborIterator n_it = shared_bdy_vertices.begin_neighbor_domains ( );
                  n_it != shared_bdy_vertices.end_neighbor_domains ( ); ++n_it )
            {
                recv_communicators[k]->wait ( );
                unpack_entities ( recv_entity_packages[k], builder );

                ++k;
            }

            // release send communicators
            for ( k = 0; k < num_neighbor_domains; ++k )
            {
                send_communicators[k]->wait ( );
                delete send_communicators[k];
                delete recv_communicators[k];
            }
            LOG_DEBUG ( 2, "Exit communicate_ghost_cell_facets() on " << local_domain );
        }

        void broadcast_mesh ( const Mesh* mesh, const MPI_Comm& comm, MeshBuilder& builder, int root )
        {
            int rank = -1;
            MPI_Comm_rank ( comm, &rank );
            EntityPackage sent_entities, recv_entities;
            const bool is_root = ( rank == root );

            if ( is_root )
            {
                assert ( mesh != 0 );
                pack_entities ( *mesh, mesh->tdim ( ), &sent_entities );
            }

            MpiBroadcast bcast ( comm, root );
            bcast.communicate ( &sent_entities, &recv_entities );

            unpack_entities ( recv_entities, builder );
        }

    }
} // namespace hiflow
