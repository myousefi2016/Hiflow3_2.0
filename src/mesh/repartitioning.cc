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

/// @author Simon Gawlok

#include <mpi.h>
#include <map>

// public headers
#include "common/macros.h"
#include "common/log.h"
#include "attributes.h"
#include "entity.h"
#include "iterator.h"
#include "mesh_db_view.h"
#include "partitioning.h"
#include "reader.h"
#include "common/timer.h"
#include "mesh/geometric_tools.h"

// private headers
#include "communication.h"
#include "mpi_communication.h"

#include "repartitioning.h"

#ifdef WITH_PARMETIS
#    include <parmetis.h>
#endif

const int DEBUG_LEVEL = 5;

namespace hiflow
{
    namespace mesh
    {
#ifdef WITH_PARMETIS

        //////////////// MpiAllToAll ////////////////

        MpiAllToAll::MpiAllToAll ( const MPI_Comm& mpi_communicator )
        : mpi_comm_ ( mpi_communicator )
        {
        }

        void MpiAllToAll::communicate ( std::vector<EntityPackage>& sent_entities,
                                        std::vector<EntityPackage>& received_entities ) const
        {
            LOG_DEBUG ( 2, "Communication starting on process " << rank ( ) );
            int num_procs = -1;
            MPI_Comm_size ( mpi_comm_, &num_procs );

            received_entities.clear ( );
            received_entities.resize ( num_procs );

            // set dimensions
            for ( int i = 0; i < num_procs; ++i )
            {
                received_entities[i].tdim = sent_entities[i].tdim;
                received_entities[i].gdim = sent_entities[i].gdim;
            }

            // broadcast coordinates vector size
            std::vector<int> num_coordinates_sent ( num_procs, 0 );

            for ( int i = 0; i < num_procs; ++i )
            {
                num_coordinates_sent[i] = sent_entities[i].coords.size ( );
            }

            std::vector<int> num_coordinates_recv ( num_procs, -1 );
            LOG_DEBUG ( 2, "AllToAll num_coordinates_sent" );
            MPI_Alltoall ( vec2ptr ( num_coordinates_sent ),
                           1, MPI_INT,
                           vec2ptr ( num_coordinates_recv ),
                           1, MPI_INT, mpi_comm_ );

            for ( int i = 0; i < num_procs; ++i )
            {
                received_entities[i].coords.resize ( num_coordinates_recv[i] );
                LOG_DEBUG ( 3, "On proc " << rank ( ) << ", "
                            << num_coordinates_recv[i]
                            << " coordinates will be received from process " << i );
            }

            // broadcast coordinates vector
            LOG_DEBUG ( 2, "AllToAll coordinates" );

            // Number of total coordinates to be receives
            int num_total_elements = 0;
            // Displacements for receiving coordinates
            std::vector<int> rdispl ( num_procs + 1, 0 );

            for ( int i = 0; i < num_procs; ++i )
            {
                // count number of total elements that will be received
                num_total_elements += num_coordinates_recv[i];
                rdispl[i + 1] = rdispl[i] + num_coordinates_recv[i];
            }

            // Create receive buffer
            std::vector<Coordinate> coordinates_recv ( num_total_elements, 0. );

            // Create data structures for sending
            std::vector<Coordinate> coordinates_sent;
            std::vector<int> sdispl ( num_procs + 1, 0 );
            for ( int i = 0; i < num_procs; ++i )
            {
                coordinates_sent.insert ( coordinates_sent.end ( ),
                                          sent_entities[i].coords.begin ( ),
                                          sent_entities[i].coords.end ( ) );
                sdispl[i + 1] = sdispl[i] + sent_entities[i].coords.size ( );
            }

            MPI_Alltoallv ( vec2ptr ( coordinates_sent ),
                            vec2ptr ( num_coordinates_sent ),
                            vec2ptr ( sdispl ),
                            MPI_DOUBLE,
                            vec2ptr ( coordinates_recv ),
                            vec2ptr ( num_coordinates_recv ),
                            vec2ptr ( rdispl ),
                            MPI_DOUBLE,
                            mpi_comm_ );

            // Unpack coordinates to received entities
            for ( int i = 0; i < num_procs; ++i )
            {
                for ( int j = 0; j < num_coordinates_recv[i]; ++j )
                {
                    received_entities[i].coords[j] = coordinates_recv[rdispl[i] + j];
                }
            }

            num_coordinates_sent.clear ( );
            num_coordinates_recv.clear ( );
            rdispl.clear ( );
            coordinates_recv.clear ( );
            coordinates_sent.clear ( );
            sdispl.clear ( );

            // scatter number of received entities for each proc
            LOG_DEBUG ( 2, "AllToAll offset sizes" );
            std::vector<int> offset_sizes_sent ( num_procs, 0 );
            std::vector<int> offset_sizes_recv ( num_procs, 0 );

            for ( int i = 0; i < num_procs; ++i )
            {
                offset_sizes_sent[i] = sent_entities[i].offsets.size ( );
            }

            MPI_Alltoall ( vec2ptr ( offset_sizes_sent ), 1, MPI_INT,
                           vec2ptr ( offset_sizes_recv ), 1, MPI_INT,
                           mpi_comm_ );

            for ( int i = 0; i < num_procs; ++i )
            {
                received_entities[i].offsets.resize ( offset_sizes_recv[i] );
                LOG_DEBUG ( 3, "On proc " << rank ( ) << ", "
                            << offset_sizes_recv[i]
                            << " offsets will be received from process " << i );
            }

            LOG_DEBUG ( 2, "AllToAll offset vector" );

            num_total_elements = 0;
            rdispl.resize ( num_procs + 1, 0 );
            for ( int i = 0; i < num_procs; ++i )
            {
                // count number of total elements that will be received
                num_total_elements += offset_sizes_recv[i];
                rdispl[i + 1] = rdispl[i] + offset_sizes_recv[i];
            }

            std::vector<int> received_offsets ( num_total_elements );

            sdispl.resize ( num_procs + 1, 0 );
            std::vector<int> sent_offsets;
            for ( int i = 0; i < num_procs; ++i )
            {
                sent_offsets.insert ( sent_offsets.end ( ),
                                      sent_entities[i].offsets.begin ( ),
                                      sent_entities[i].offsets.end ( ) );
                sdispl[i + 1] = sdispl[i] + sent_entities[i].offsets.size ( );

            }

            MPI_Alltoallv ( vec2ptr ( sent_offsets ),
                            vec2ptr ( offset_sizes_sent ),
                            vec2ptr ( sdispl ), MPI_INT,
                            vec2ptr ( received_offsets ),
                            vec2ptr ( offset_sizes_recv ),
                            vec2ptr ( rdispl ),
                            MPI_INT,
                            mpi_comm_ );

            // Unpack received offsets
            for ( int i = 0; i < num_procs; ++i )
            {
                for ( int j = 0; j < offset_sizes_recv[i]; ++j )
                {
                    received_entities[i].offsets[j] = received_offsets[rdispl[i] + j];
                }
            }

            offset_sizes_sent.clear ( );
            offset_sizes_recv.clear ( );
            received_offsets.clear ( );
            sent_offsets.clear ( );
            rdispl.clear ( );
            sdispl.clear ( );

            // scatter number of received entities for each proc
            LOG_DEBUG ( 2, "AllToAll material_number sizes" );
            std::vector<int> mat_sizes_sent ( num_procs, 0 );
            std::vector<int> mat_sizes_recv ( num_procs, 0 );

            for ( int i = 0; i < num_procs; ++i )
            {
                mat_sizes_sent[i] = sent_entities[i].material_numbers.size ( );
            }

            MPI_Alltoall ( vec2ptr ( mat_sizes_sent ), 1, MPI_INT,
                           vec2ptr ( mat_sizes_recv ), 1, MPI_INT,
                           mpi_comm_ );

            for ( int i = 0; i < num_procs; ++i )
            {
                received_entities[i].material_numbers.resize ( mat_sizes_recv[i] );
                LOG_DEBUG ( 3, "On proc " << rank ( ) << ", "
                            << mat_sizes_recv[i]
                            << " material numbers will be received from process " << i );
            }

            LOG_DEBUG ( 2, "AllToAll material_number vector" );

            num_total_elements = 0;
            rdispl.clear ( );
            rdispl.resize ( num_procs + 1, 0 );
            for ( int i = 0; i < num_procs; ++i )
            {
                // count number of total elements that will be received
                num_total_elements += mat_sizes_recv[i];
                rdispl[i + 1] = rdispl[i] + mat_sizes_recv[i];
            }

            std::vector<int> received_mat ( num_total_elements );

            sdispl.clear ( );
            sdispl.resize ( num_procs + 1, 0 );
            std::vector<int> sent_mat;
            for ( int i = 0; i < num_procs; ++i )
            {
                sent_mat.insert ( sent_mat.end ( ),
                                  sent_entities[i].material_numbers.begin ( ),
                                  sent_entities[i].material_numbers.end ( ) );
                sdispl[i + 1] = sdispl[i] + sent_entities[i].material_numbers.size ( );

            }

            MPI_Alltoallv ( vec2ptr ( sent_mat ),
                            vec2ptr ( mat_sizes_sent ),
                            vec2ptr ( sdispl ), MPI_INT,
                            vec2ptr ( received_mat ),
                            vec2ptr ( mat_sizes_recv ),
                            vec2ptr ( rdispl ),
                            MPI_INT,
                            mpi_comm_ );

            // Unpack received material numbers
            for ( int i = 0; i < num_procs; ++i )
            {
                for ( int j = 0; j < mat_sizes_recv[i]; ++j )
                {
                    received_entities[i].material_numbers[j] = received_mat[rdispl[i] + j];
                }
            }

            mat_sizes_sent.clear ( );
            mat_sizes_recv.clear ( );
            sent_mat.clear ( );
            received_mat.clear ( );
            rdispl.clear ( );
            sdispl.clear ( );

            // scatter number of received connections for each proc
            LOG_DEBUG ( 2, "AllToAll connections sizes" );
            std::vector<int> conn_sizes_sent ( num_procs, 0 );
            std::vector<int> conn_sizes_recv ( num_procs, 0 );

            for ( int i = 0; i < num_procs; ++i )
            {
                conn_sizes_sent[i] = sent_entities[i].connections.size ( );
            }

            MPI_Alltoall ( vec2ptr ( conn_sizes_sent ), 1, MPI_INT,
                           vec2ptr ( conn_sizes_recv ), 1, MPI_INT,
                           mpi_comm_ );

            for ( int i = 0; i < num_procs; ++i )
            {
                received_entities[i].connections.resize ( conn_sizes_recv[i] );
                LOG_DEBUG ( 3, "On proc " << rank ( ) << ", "
                            << conn_sizes_recv[i]
                            << " connections will be received from process " << i );
            }

            LOG_DEBUG ( 2, "AllToAll connections vector" );

            num_total_elements = 0;
            rdispl.clear ( );
            rdispl.resize ( num_procs + 1, 0 );
            for ( int i = 0; i < num_procs; ++i )
            {
                // count number of total elements that will be received
                num_total_elements += conn_sizes_recv[i];
                rdispl[i + 1] = rdispl[i] + conn_sizes_recv[i];
            }

            std::vector<int> received_conn ( num_total_elements );

            sdispl.clear ( );
            sdispl.resize ( num_procs + 1, 0 );
            std::vector<int> sent_conn;
            for ( int i = 0; i < num_procs; ++i )
            {
                sent_conn.insert ( sent_conn.end ( ),
                                   sent_entities[i].connections.begin ( ),
                                   sent_entities[i].connections.end ( ) );
                sdispl[i + 1] = sdispl[i] + sent_entities[i].connections.size ( );

            }

            MPI_Alltoallv ( vec2ptr ( sent_conn ),
                            vec2ptr ( conn_sizes_sent ),
                            vec2ptr ( sdispl ), MPI_INT,
                            vec2ptr ( received_conn ),
                            vec2ptr ( conn_sizes_recv ),
                            vec2ptr ( rdispl ),
                            MPI_INT,
                            mpi_comm_ );

            // Unpack received connections
            for ( int i = 0; i < num_procs; ++i )
            {
                for ( int j = 0; j < conn_sizes_recv[i]; ++j )
                {
                    received_entities[i].connections[j] = received_conn[rdispl[i] + j];
                }
            }

            conn_sizes_sent.clear ( );
            conn_sizes_recv.clear ( );
            sent_conn.clear ( );
            received_conn.clear ( );
            rdispl.clear ( );
            sdispl.clear ( );

            LOG_DEBUG ( 2, "Communication done on process " << rank ( ) );
        }

        int MpiAllToAll::rank ( ) const
        {
            int rank = -1;
            MPI_Comm_rank ( mpi_comm_, &rank );
            return rank;
        }

        ///////////////////////////////////////////

        /// \details Computes the dual graph of a distributed mesh, with two cells
        /// being considered as connected if they share a vertex.
        /// \param[in] mesh  The mesh
        /// \param[in] comm MPI communicator
        /// \param[out] graph  The computed graph (allocated by caller)
        /// \param[in] vtxdist Vertex, i.e., cell-distribution among the processors

        void compute_dual_graph_distributed ( const MeshPtr mesh,
                                              MPI_Comm& comm,
                                              Graph* graph,
                                              const std::vector<int>& vtxdist )
        {
            int my_rank = -1, num_ranks = -2;
            MPI_Comm_rank ( comm, &my_rank );
            MPI_Comm_size ( comm, &num_ranks );
            assert ( 0 <= my_rank );
            assert ( my_rank < num_ranks );

            assert ( graph != 0 );

            Timer dual_graph_timer;
            dual_graph_timer.start ( );

            // variables used in loop
            const TDim tdim = mesh->tdim ( );
            const EntityIterator end = mesh->end ( tdim );
            int num_graph_edges = 0;
            int c = 0;

            std::map<int, std::map<int, float> > adjacency_map;

            // Fill adjacency map with adjacency and weight information
            for ( EntityIterator it = mesh->begin ( tdim ); it != end; ++it )
            {
                int remote_index;
                mesh->get_attribute_value ( "_remote_index_", tdim, it->index ( ), &remote_index );
                if ( remote_index == -1 )
                {
                    adjacency_map[it->index ( )];
                    for ( int i = tdim; i--; )
                    {
                        for ( IncidentEntityIterator iit = mesh->begin_incident ( *it, i ),
                              e_it = mesh->end_incident ( *it, i );
                              iit != e_it; ++iit )
                        {
                            for ( IncidentEntityIterator iiit = mesh->begin_incident ( *iit, tdim ),
                                  e_iiit = mesh->end_incident ( *iit, tdim );
                                  iiit != e_iiit; ++iiit )
                            {
                                if ( it->index ( ) != iiit->index ( ) )
                                {
                                    int remote_index_2;
                                    mesh->get_attribute_value ( "_remote_index_", tdim, iiit->index ( ), &remote_index_2 );
                                    if ( remote_index_2 == -1 )
                                    {
                                        std::map<int, float>::iterator f_it = adjacency_map[it->index ( )].find ( iiit->index ( ) + vtxdist[my_rank] );
                                        if ( f_it == adjacency_map[it->index ( )].end ( ) )
                                        {
                                            adjacency_map[it->index ( )][iiit->index ( ) + vtxdist[my_rank]] = std::pow ( static_cast < float > ( tdim + 1 ), i );
                                        }
                                    }
                                    else
                                    {
                                        int subdomain;
                                        mesh->get_attribute_value ( "_sub_domain_", tdim, iiit->index ( ), &subdomain );
                                        std::map<int, float>::iterator f_it = adjacency_map[it->index ( )].find ( remote_index_2 + vtxdist[subdomain] );
                                        if ( f_it == adjacency_map[it->index ( )].end ( ) )
                                        {
                                            adjacency_map[it->index ( )][remote_index_2 + vtxdist[subdomain]] = std::pow ( static_cast < float > ( tdim + 1 ), i );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Translate adjacency_map to graph
            for ( std::map<int, std::map<int, float> >::iterator it = adjacency_map.begin ( ),
                  e_it = adjacency_map.end ( );
                  it != e_it;
                  ++it )
            {
                graph->node_offsets[c++] = num_graph_edges;
                for ( std::map<int, float>::iterator iit = it->second.begin ( ),
                      e_iit = it->second.end ( );
                      iit != e_iit;
                      ++iit )
                {
                    graph->adjacent_nodes.push_back ( iit->first );
                    graph->adjacent_weights.push_back ( static_cast < int > ( iit->second ) );
                    ++num_graph_edges;
                }
            }

            adjacency_map.clear ( );

            // set the last offset to the length of the adjacent_nodes array
            assert ( num_graph_edges == static_cast < int > ( graph->adjacent_nodes.size ( ) ) );
            graph->node_offsets[c] = num_graph_edges;

            dual_graph_timer.stop ( );
            LOG_INFO ( "Time for computing distributed dual graph", dual_graph_timer.get_duration ( ) );
        }

        /// \details Repartition already distributed mesh in parallel.
        /// The function returns the part of the mesh belonging to
        /// the local process.
        ///
        /// \param mmesh  Mesh to be repartitioned
        /// \param comm         MPI communicator used for communication.
        /// \param gpartitioner The graph partitioner to be employed, i.e., a ParMetisGraphPartitioner.
        /// \return  Shared pointer to local part of the repartitioned mesh.

        MeshPtr repartition_mesh ( const MeshPtr mesh,
                                   MPI_Comm& comm,
                                   const ParMetisGraphPartitioner* gpartitioner )
        {
            assert ( gpartitioner != 0 );

            int my_rank = -1, num_ranks = -2;
            MPI_Comm_rank ( comm, &my_rank );
            MPI_Comm_size ( comm, &num_ranks );
            assert ( 0 <= my_rank );
            assert ( my_rank < num_ranks );

            // check that mesh is non-null on all processes
            assert ( mesh != 0 );
            std::vector<MasterSlave> period = mesh->get_period ( );

            Graph mesh_graph;
            mesh_graph.clear ( );

            std::vector<int> vtxdist;

            mesh_graph.node_offsets.clear ( );
            mesh_graph.adjacent_nodes.clear ( );

            const TDim tdim = mesh->tdim ( );
            const GDim gdim = mesh->gdim ( );
            const EntityCount num_cells = mesh->num_entities ( tdim );
            int num_local_cells = num_cells;
            mesh_graph.node_offsets.resize ( num_cells + 1 );

            vtxdist.resize ( num_ranks + 1 );
            std::vector<int> local_cell_nums ( num_ranks );

            // Compute cell offsets of processes
            MPI_Allgather ( &num_local_cells,
                            1,
                            MPI_INT,
                            vec2ptr ( local_cell_nums ),
                            1,
                            MPI_INT,
                            comm );

            for ( int i = 0; i < local_cell_nums.size ( ); ++i )
            {
                vtxdist[i + 1] = vtxdist[i] + local_cell_nums[i];
            }

            Timer ghost_cell_computation;
            ghost_cell_computation.start ( );

            // Compute mesh with ghost cells in order to get cell numbers of neighbouring
            // processes
            SharedVertexTable shared_verts;
            MeshPtr mesh_with_ghost;

            mesh_with_ghost = compute_ghost_cells ( *mesh, comm, shared_verts );

            ghost_cell_computation.stop ( );
            LOG_INFO ( "Time to compute ghost cells", ghost_cell_computation.get_duration ( ) );

            // only master proc computes the graph
            compute_dual_graph_distributed ( mesh_with_ghost, comm, &mesh_graph, vtxdist );

            Timer partition_computation;
            partition_computation.start ( );

            // Edges are weighted
            int wgtflag = 1;

            // use C-style numbering
            int numflag = 0;

            // variable to store number of edge cuts
            int num_edge_cuts;

            int ncon = 1; // number of weights associated with each vertex
            int nparts = num_ranks;
            int* vsize = 0; // array for communication cost of vertices (not used)
            float* tpwgts = new float[ncon * nparts];
            for ( unsigned i = 0; i < ncon * nparts; ++i )
                tpwgts[i] = 1.0 / ( float ) ( nparts );
            float ubvec = 1.05; // max load imbalance (not used)
            MPI_Comm comm_copy = comm;
            int options[3] = { 0, 0, 0 };

            std::vector<int> partitioning ( vtxdist[my_rank + 1] - vtxdist[my_rank], 0 );
            int metis_return =
                    ParMETIS_V3_PartKway ( vec2ptr ( vtxdist ),
                                           vec2ptr ( mesh_graph.node_offsets ),
                                           vec2ptr ( mesh_graph.adjacent_nodes ),
                                           NULL,
                                           vec2ptr ( mesh_graph.adjacent_weights ),
                                           &wgtflag,
                                           &numflag, &ncon, &nparts, tpwgts, &ubvec, options,
                                           &num_edge_cuts, vec2ptr ( partitioning ),
                                           &comm_copy );
            assert ( metis_return == METIS_OK );

            delete [] tpwgts;

            partition_computation.stop ( );
            LOG_INFO ( "Time for computing repartitioning", partition_computation.get_duration ( ) );

            LOG_INFO ( "Number of edge cuts", num_edge_cuts );

            // Builder for the resulting mesh
            MeshBuilder* builder;

            builder = new MeshDbViewBuilder ( tdim, gdim, period );

            mesh_graph.node_offsets.clear ( );
            mesh_graph.adjacent_nodes.clear ( );
            mesh_graph.adjacent_weights.clear ( );
            vtxdist.clear ( );

            // prepare entities to send on all processors
            Timer fill_mesh_builder;
            fill_mesh_builder.start ( );

            std::vector<EntityPackage> sent_entities;
            std::vector<EntityPackage> recv_entities;
            std::vector<EntityPackage> sent_facet_entities;
            std::vector<EntityPackage> recv_facet_entities;
            std::vector<int> num_entities_on_proc;
            std::vector<int> num_facet_entities_on_proc;

            // compute distribution of cells from the partitioning
            CellDistribution distribution;
            compute_cell_distribution ( num_ranks, partitioning, &distribution );

            assert ( static_cast < int > ( distribution.num_cells.size ( ) ) == num_ranks );

            // compute distribution of facets from cell distribution
            CellDistribution facet_distribution;
            compute_facet_from_cell_distribution ( *mesh_with_ghost, distribution, &facet_distribution );

            // prepare sent entities
            pack_distributed_entities ( *mesh_with_ghost, tdim, distribution, sent_entities );
            pack_distributed_entities ( *mesh_with_ghost, tdim - 1, facet_distribution, sent_facet_entities );

            // Communicate cells
            Timer partition_communication;
            partition_communication.start ( );

            MpiAllToAll scatter ( comm );
            scatter.communicate ( sent_entities, recv_entities );

            partition_communication.stop ( );
            LOG_INFO ( "Time for cell communication (repartitioning)", partition_communication.get_duration ( ) );

            // unpack the received part of the mesh into a builder
            for ( int rank = 0; rank < num_ranks; ++rank )
            {
                unpack_entities ( recv_entities[rank], *builder );
            }

            // Communicate facets
            partition_communication.start ( );

            MpiAllToAll facet_scatter ( comm );
            facet_scatter.communicate ( sent_facet_entities, recv_facet_entities );

            partition_communication.stop ( );
            LOG_INFO ( "Time for facet communication (repartitioning)", partition_communication.get_duration ( ) );

            // unpack the received facets into a builder
            for ( int rank = 0; rank < num_ranks; ++rank )
            {
                unpack_entities ( recv_facet_entities[rank], *builder );
            }

            fill_mesh_builder.stop ( );
            LOG_INFO ( "Overall time to fill mesh builder object", fill_mesh_builder.get_duration ( ) );

            // build the local mesh
            MeshPtr recv_mesh ( builder->build ( ) );

            return recv_mesh;
        }
#endif
    }
}
