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

#include "partitioning.h"

#include <numeric>
#include <map>
#include <cmath>

#include "common/macros.h"
#include "common/log.h"
#include "config.h"
#include "mesh.h"
#include "iterator.h"
#include "types.h"

#ifndef DEBUG_LEVEL
#    define DEBUG_LEVEL 10
#endif

#ifdef WITH_METIS

#    if METIS_VERSION == 5
extern "C"
{
#        include <metis.h>
}

#    elif METIS_VERSION == 4
// workaround for problematic header files
typedef int idxtype;
extern "C"
{
    void METIS_PartGraphKway ( int *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype * );
}
#    else
#        error "Unknown METIS version " METIS_VERSION
#    endif
#endif

#ifdef WITH_PARMETIS
#    include <parmetis.h>
#endif

namespace hiflow
{
    namespace mesh
    {

        /// \details Computes the dual graph of a mesh, with two cells
        /// being considered as connected if they share a vertex.
        /// \param[in] mesh  The mesh
        /// \param[out] graph  The computed graph (allocated by caller)

        void compute_dual_graph ( const Mesh& mesh, Graph* graph )
        {
            assert ( graph != 0 );
            graph->node_offsets.clear ( );
            graph->adjacent_nodes.clear ( );

            const TDim tdim = mesh.tdim ( );
            const EntityCount num_cells = mesh.num_entities ( tdim );
            graph->node_offsets.resize ( num_cells + 1 );

            // variables used in loop
            const EntityIterator end = mesh.end ( tdim );
            int num_graph_edges = 0;
            int c = 0;

            std::map<int, std::map<int, float> > adjacency_map;

            // Fill adjacency map with adjacency and weight information
            for ( EntityIterator it = mesh.begin ( tdim ); it != end; ++it )
            {
                adjacency_map[it->index ( )];
                for ( int i = tdim; i--; )
                {
                    for ( IncidentEntityIterator iit = mesh.begin_incident ( *it, i ),
                          e_it = mesh.end_incident ( *it, i );
                          iit != e_it; ++iit )
                    {
                        for ( IncidentEntityIterator iiit = mesh.begin_incident ( *iit, tdim ),
                              e_iiit = mesh.end_incident ( *iit, tdim );
                              iiit != e_iiit; ++iiit )
                        {
                            if ( it->index ( ) != iiit->index ( ) )
                            {
                                std::map<int, float>::iterator f_it = adjacency_map[it->index ( )].find ( iiit->index ( ) );
                                if ( f_it == adjacency_map[it->index ( )].end ( ) )
                                {
                                    adjacency_map[it->index ( )][iiit->index ( )] = std::pow ( static_cast < float > ( tdim + 1 ), i );
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

            /*for ( EntityIterator it = mesh.begin( tdim ); it != end; ++it )
            {
                graph->node_offsets[c++] = num_graph_edges;
                for ( IncidentEntityIterator iit = mesh.begin_incident( *it, tdim );
                      iit != mesh.end_incident( *it, tdim ); ++iit )
                {
                    // insert entry in graph
                    graph->adjacent_nodes.push_back( iit->index( ) );
                    ++num_graph_edges;
                }
            }*/

            assert ( c == num_cells );

            // set the last offset to the length of the adjacent_nodes array
            assert ( num_graph_edges == static_cast < int > ( graph->adjacent_nodes.size ( ) ) );
            graph->node_offsets[c] = num_graph_edges;
        }

        void compute_dual_graph_facet ( const Mesh& mesh, Graph* graph )
        {
            assert ( graph != 0 );
            graph->node_offsets.clear ( );
            graph->adjacent_nodes.clear ( );

            const TDim tdim = mesh.tdim ( );
            const EntityCount num_cells = mesh.num_entities ( tdim );
            graph->node_offsets.resize ( num_cells + 1 );

            // variables used in loop
            const EntityIterator end = mesh.end ( tdim );
            int num_graph_edges = 0;
            int c = 0;

            for ( EntityIterator it = mesh.begin ( tdim ); it != end; ++it )
            {
                graph->node_offsets[c++] = num_graph_edges;
                const IncidentEntityIterator facet_end = mesh.end_incident ( *it, tdim - 1 );
                for ( IncidentEntityIterator facet_it = mesh.begin_incident ( *it, tdim - 1 );
                      facet_it != facet_end; ++facet_it )
                {
                    for ( IncidentEntityIterator neighbor_cell_it = mesh.begin_incident ( *facet_it, tdim );
                          neighbor_cell_it != mesh.end_incident ( *facet_it, tdim ); ++neighbor_cell_it )
                    {
                        if ( it->index ( ) != neighbor_cell_it->index ( ) )
                        {
                            // insert entry in graph
                            graph->adjacent_nodes.push_back ( neighbor_cell_it->index ( ) );
                            ++num_graph_edges;
                        }
                    }
                }
            }

            assert ( c == num_cells );

            // set the last offset to the length of the adjacent_nodes array
            assert ( num_graph_edges == static_cast < int > ( graph->adjacent_nodes.size ( ) ) );
            graph->node_offsets[c] = num_graph_edges;
        }

        /// \details Distributes a graph which is given on the master process
        /// to all processes in the MPI communicator. The nodes are simply
        /// distributed in contiguous blocks of size (no. nodes / no. procs),
        /// the last proc taking the rest.
        ///
        /// \param in The graph to be distributed, significant only on master process.
        /// \param master_rank MPI rank of the master process.
        /// \param comm MPI communicator of all processes taking a part of the graph.
        /// \param out Graph representing the local part which has been distributed to this process.
        /// The offsets are adjusted to local meaning, the adjacencies keep global meaning.
        /// The output graph can be used as local argument for ParMetisGraphPartitioner.

        void distribute_graph_naive ( const Graph& in,
                                      const int master_rank,
                                      const MPI_Comm& comm,
                                      Graph& out )
        {
            int rk = -1, sz = -2;
            MPI_Comm_rank ( comm, &rk );
            MPI_Comm_size ( comm, &sz );
            assert ( 0 <= rk );
            assert ( 0 <= master_rank );
            assert ( sz > rk );
            assert ( sz > master_rank );

            if ( sz == 1 )
            {
                // only one proc => copy input to output
                out = in;
                return;
            }
            out.clear ( );

            // master proc knows sizes
            int n_off = in.node_offsets.size ( );
            int n_adj = in.adjacent_nodes.size ( );
            int n_adj_wgts = in.adjacent_weights.size ( );
            MPI_Bcast ( &n_off, 1, MPI_INT, master_rank, comm );
            MPI_Bcast ( &n_adj, 1, MPI_INT, master_rank, comm );
            MPI_Bcast ( &n_adj_wgts, 1, MPI_INT, master_rank, comm );

            // helper variables
            std::vector<int> off_vec, adj_vec, adj_wgts;
            if ( rk == master_rank )
            {
                // master proc copies in to helper variables
                off_vec = in.node_offsets;
                adj_vec = in.adjacent_nodes;
                adj_wgts = in.adjacent_weights;
            }
            else
            {
                // other procs allocate helper variables
                off_vec.resize ( n_off, 0 );
                adj_vec.resize ( n_adj, 0 );
                adj_wgts.resize ( n_adj_wgts, 0 );
            }

            // bcast whole graph
            MPI_Bcast ( vec2ptr ( off_vec ), n_off, MPI_INT, master_rank, comm );
            MPI_Bcast ( vec2ptr ( adj_vec ), n_adj, MPI_INT, master_rank, comm );
            MPI_Bcast ( vec2ptr ( adj_wgts ), n_adj_wgts, MPI_INT, master_rank, comm );

            // nb. of nodes for naive distribution
            const int nodes_per_proc = ( n_off - 1 ) / sz; // integer division
            const int rest = n_off - 1 - ( nodes_per_proc * ( sz - 1 ) );
            assert ( nodes_per_proc > 0 );
            assert ( rest >= nodes_per_proc );

            int my_nb_nodes;
            if ( rk < sz - 1 ) my_nb_nodes = nodes_per_proc;
            else my_nb_nodes = rest;
            const int my_off = rk * nodes_per_proc;
            const int next_off = my_off + my_nb_nodes;

            // allocate output offset and adjacency vector
            out.node_offsets.resize ( my_nb_nodes + 1, 0 );
            out.adjacent_nodes.resize ( off_vec[next_off] - off_vec[my_off], 0 );
            out.adjacent_weights.resize ( off_vec[next_off] - off_vec[my_off], 0 );

            // copy my part of the offsets
            std::copy ( off_vec.begin ( ) + my_off,
                        off_vec.begin ( ) + next_off + 1,
                        out.node_offsets.begin ( ) );

            // correct offsets
            const int correction = out.node_offsets[0];
            for ( int i = 0; i <= my_nb_nodes; ++i )
                out.node_offsets[i] -= correction;

            // copy my part of the adjacency list
            std::copy ( adj_vec.begin ( ) + off_vec[my_off],
                        adj_vec.begin ( ) + off_vec[next_off],
                        out.adjacent_nodes.begin ( ) );

            // copy my part of the adjacency weights
            std::copy ( adj_wgts.begin ( ) + off_vec[my_off],
                        adj_wgts.begin ( ) + off_vec[next_off],
                        out.adjacent_weights.begin ( ) );
        }

        void compute_cell_distribution ( const int num_partitions,
                                         const std::vector<int>& partitioning,
                                         CellDistribution* distribution )
        {
            // The partitioning gives the map
            //
            // cell_index -> sub_domain_number
            //
            // which should be converted to a vector sorted according to
            // the sub_domain_number, containing the indices:s of the
            // cells.

            assert ( distribution != 0 );
            distribution->num_cells.clear ( );
            distribution->num_cells.resize ( num_partitions );

            // compute number of cells in each sub-domain
            for ( std::vector<int>::const_iterator it = partitioning.begin ( );
                  it != partitioning.end ( ); ++it )
            {
                ++distribution->num_cells[*it];
            }

            // compute offsets for each sub-domain in cell_ids vector
            std::vector<int> cell_offsets ( num_partitions );
            cell_offsets[0] = 0;
            std::partial_sum ( distribution->num_cells.begin ( ),
                               distribution->num_cells.end ( ) - 1, cell_offsets.begin ( ) + 1 );

            // compute cell_indices vector
            distribution->cell_indices.clear ( );
            distribution->cell_indices.resize ( partitioning.size ( ) );
            for ( int i = 0; i < static_cast < int > ( partitioning.size ( ) ); ++i )
            {
                const int sub_domain = partitioning[i];
                distribution->cell_indices[cell_offsets[sub_domain]] = i;
                ++cell_offsets[sub_domain];
            }
        }

        void compute_facet_from_cell_distribution ( const Mesh& mesh,
                                                    const CellDistribution& cell_distribution,
                                                    CellDistribution* facet_distribution )
        {
            int tdim = mesh.tdim ( );
            int num_partitions = cell_distribution.num_cells.size ( );
            facet_distribution->num_cells.clear ( );
            facet_distribution->num_cells.resize ( num_partitions, 0 );
            facet_distribution->cell_indices.clear ( );

            int offset = 0;
            // iterate subdomains
            for ( int p = 0; p < num_partitions; ++p )
            {
                // temporal container to collect facet distribution
                std::vector<int> temp_facet_list;
                // iterate cells of the current subdomain
                assert ( cell_distribution.cell_indices.size ( ) >= offset + cell_distribution.num_cells[p] );
                for ( int c = offset; c < offset + cell_distribution.num_cells[p]; ++c )
                {
                    // get current cell
                    Entity cell = mesh.get_entity ( tdim, cell_distribution.cell_indices[c] );
                    // add facet indices of the current cell to the temporal facet distribution
                    for ( IncidentEntityIterator it = cell.begin_incident ( tdim - 1 ); it != cell.end_incident ( tdim - 1 ); ++it )
                    {
                        temp_facet_list.push_back ( it->index ( ) );
                    }
                }
                offset += cell_distribution.num_cells[p];
                // prepare facet indices
                std::sort ( temp_facet_list.begin ( ), temp_facet_list.end ( ) );
                std::vector<int>::iterator last_unique = std::unique ( temp_facet_list.begin ( ), temp_facet_list.end ( ) );
                temp_facet_list.resize ( last_unique - temp_facet_list.begin ( ) );
                facet_distribution->num_cells[p] = temp_facet_list.size ( );
                facet_distribution->cell_indices.insert ( facet_distribution->cell_indices.end ( ), temp_facet_list.begin ( ), temp_facet_list.end ( ) );
            }
            assert ( facet_distribution->cell_indices.size ( ) == std::accumulate ( facet_distribution->num_cells.begin ( ), facet_distribution->num_cells.end ( ), 0 ) );
        }

        //////// NaiveGraphPartitioner ////////////////

        void NaiveGraphPartitioner::partition_graph ( const Graph& graph,
                                                      const int master_rank,
                                                      const std::vector<int>* nodal_weights,
                                                      const std::vector<int>* edge_weights,
                                                      const MPI_Comm& comm,
                                                      std::vector<int>* partitioning ) const
        {
            assert ( partitioning != 0 );

            int rk = -1, sz = -2;
            MPI_Comm_rank ( comm, &rk );
            MPI_Comm_size ( comm, &sz );
            assert ( 0 <= rk );
            assert ( 0 <= master_rank );
            assert ( sz > rk );
            assert ( sz > master_rank );

            partitioning->resize ( graph.num_nodes ( ) );
            const int nodes_per_partition = graph.num_nodes ( ) / sz;

            // give same number of nodes to all but the last partition
            int node_index = 0;
            for ( int p = 0; p < sz - 1; ++p )
            {
                for ( int i = 0; i < nodes_per_partition; ++i )
                {
                    ( *partitioning )[node_index] = p;
                    ++node_index;
                }
            }

            // assign rest of nodes to last partition
            for (; node_index < graph.num_nodes ( ); ++node_index )
            {
                ( *partitioning )[node_index] = sz - 1;
            }

            assert ( node_index == graph.num_nodes ( ) );

            LOG_INFO ( "NaiveGraphPartitioner", "partitioning done" );
        }

        //////// MetisGraphPartitioner ////////////////
#ifdef WITH_METIS

        void MetisGraphPartitioner::partition_graph ( const Graph& graph,
                                                      const int master_rank,
                                                      const std::vector<int>* nodal_weights,
                                                      const std::vector<int>* edge_weights,
                                                      const MPI_Comm& comm,
                                                      std::vector<int>* partitioning ) const
        {
            // NB const casts in this functions are necessary to use the C
            // interface. This wrapper only works with METIS configured to use
            // 32-bit integers for the indices (idx_t) and 64-bit doubles for the
            // floating point values. This corresponds to activating the compile option
            // METIS_USE_DOUBLEPRECISION, and deactivating METIS_USE_LONGINDEX .

            int rk = -1, sz = -2;
            MPI_Comm_rank ( comm, &rk );
            MPI_Comm_size ( comm, &sz );
            assert ( 0 <= rk );
            assert ( 0 <= master_rank );
            assert ( sz > rk );
            assert ( sz > master_rank );

            // METIS does not handle corner case sz == 1 in a
            // good way.
            if ( sz == 1 )
            {
                partitioning->clear ( );
                partitioning->resize ( graph.num_nodes ( ), 0 );
                return;
            }

            // setup arguments
            int num_nodes = graph.num_nodes ( );
            assert ( num_nodes > sz );

            int* xadj = const_cast < int* > ( vec2ptr ( graph.node_offsets ) );
            int* adjncy = const_cast < int* > ( vec2ptr ( graph.adjacent_nodes ) );

            // Should be 0 for no weights, 1 for edge weights, 2 for nodal
            // weights, and 3 for both nodal and edge weights. Bit
            // operations are used to set this up.
            int wgtflag = 0;

            int* vwgt = 0;
            if ( edge_weights )
            {
                vwgt = const_cast < int* > ( vec2ptr ( *edge_weights ) );
                wgtflag = wgtflag | 0x01;
            }

            int* adjwgt = 0;
            if ( nodal_weights )
            {
                vwgt = const_cast < int* > ( vec2ptr ( *nodal_weights ) );
                wgtflag = wgtflag | 0x02;
            }

#    if METIS_VERSION >= 5
            int options[METIS_NOPTIONS];
            METIS_SetDefaultOptions ( options );
            options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // minimize communication volume
            options[METIS_OPTION_NUMBERING] = 0; // C style numbering starting at 0
            options[METIS_OPTION_MINCONN] = 1; // minimize maximum connectivity of partitions
            options[METIS_OPTION_CONTIG] = 1; // force contiguous partitions
#    else
            // use default options
            int options[5] = { 0, 0, 0, 0, 0 };
#    endif

            // use C-style numbering
            int numflag = 0;

            // set size of partitioning vector
            partitioning->resize ( num_nodes );

            // variable to store number of edge cuts
            int num_edge_cuts;

            // call METIS
            // use PartGraphKWay partitioning function
            int nparts = sz;
#    if METIS_VERSION >= 5
            int num_constraints = 1; // number of weights associated with each vertex
            int* vsize = 0; // array for communication cost of vertices (not used)
            float* tpwgts = 0; // target partitioned weight (not used)
            float* ubvec = 0; // max load imbalance (not used)

            METIS_PartGraphKway ( &num_nodes, &num_constraints, xadj, adjncy,
                                  vwgt, vsize, adjwgt,
                                  &nparts, tpwgts, ubvec, options,
                                  &num_edge_cuts, vec2ptr ( *partitioning ) );
#    else

            METIS_PartGraphKway ( &num_nodes, xadj, adjncy, vwgt, adjwgt,
                                  &wgtflag, &numflag, &nparts, options,
                                  &num_edge_cuts, vec2ptr ( *partitioning ) );
#    endif

            LOG_INFO ( "MetisGraphPartitioner", "partitioning done" );

        }
#endif  // WITH_METIS

#ifdef WITH_PARMETIS

        void ParMetisGraphPartitioner::partition_graph ( const Graph& graph,
                                                         const int master_rank,
                                                         const std::vector<int>* nodal_weights,
                                                         const std::vector<int>* edge_weights,
                                                         const MPI_Comm& comm,
                                                         std::vector<int>* partitioning ) const
        {
            // NB const casts in this functions are necessary to use the C
            // interface. This wrapper only works with METIS configured to use
            // 32-bit integers for the indices (idx_t) and 64-bit doubles for the
            // floating point values. This corresponds to activating the compile option
            // METIS_USE_DOUBLEPRECISION, and deactivating METIS_USE_LONGINDEX .

            int rk = -1, sz = -2;
            MPI_Comm_rank ( comm, &rk );
            MPI_Comm_size ( comm, &sz );
            assert ( 0 <= rk );
            assert ( 0 <= master_rank );
            assert ( sz > rk );
            assert ( sz > master_rank );

            // resize partitioning on all procs; master proc knows num_nodes
            int num_nodes = graph.num_nodes ( );
            MPI_Bcast ( &num_nodes, 1, MPI_INT, master_rank, comm );

            // return immediately if no partitioning is needed
            if ( sz == 1 )
            {
                // make sure all nodes are assigned to proc 0
                // (we can't be sure when just resizing)
                partitioning->resize ( num_nodes );
                for ( unsigned i = 0; i < num_nodes; ++i ) ( *partitioning )[i] = 0;
                return;
            }

            // check that there are enough nodes
            assert ( num_nodes > sz );

            partitioning->clear ( );

            Graph distributed_graph;
            distribute_graph_naive ( graph, master_rank, comm, distributed_graph );
            assert ( distributed_graph.num_nodes ( ) > 0 );

            // naive node distribution
            const int nodes_per_proc = num_nodes / sz; // integer division
            std::vector<int> vtxdist ( sz + 1, 0 );
            for ( unsigned i = 0; i < sz; ++i )
                vtxdist[i] = i * nodes_per_proc;
            vtxdist[sz] = num_nodes;

            //int* xadj   = const_cast<int*>(vec2ptr(distributed_graph.node_offsets));
            //int* adjncy = const_cast<int*>(vec2ptr(distributed_graph.adjacent_nodes));

            // Should be 0 for no weights, 1 for edge weights, 2 for nodal
            // weights, and 3 for both nodal and edge weights. Bit
            // operations are used to set this up.
            int wgtflag = 0;

            int* vwgt = 0;
            if ( edge_weights )
            {
                vwgt = const_cast < int* > ( vec2ptr ( *edge_weights ) );
                wgtflag = wgtflag | 0x01;
            }

            int* adjwgt = 0;
            if ( nodal_weights )
            {
                vwgt = const_cast < int* > ( vec2ptr ( *nodal_weights ) );
                wgtflag = wgtflag | 0x02;
            }

            // use C-style numbering
            int numflag = 0;

            // variable to store number of edge cuts
            int num_edge_cuts;

            int ncon = 1; // number of weights associated with each vertex
            int nparts = sz;
            int* vsize = 0; // array for communication cost of vertices (not used)
            float* tpwgts = new float[ncon * nparts];
            for ( unsigned i = 0; i < ncon * nparts; ++i )
                tpwgts[i] = 1.0 / ( float ) ( nparts );
            float ubvec = 1.05; // max load imbalance (not used)

            int* my_part = new int[vtxdist[rk + 1] - vtxdist[rk]];
            MPI_Comm comm_copy = comm;
            int options[3] = { 0, 0, 0 };

            int metis_return =
                    ParMETIS_V3_PartKway ( vec2ptr ( vtxdist ),
                                           vec2ptr ( distributed_graph.node_offsets ),
                                           vec2ptr ( distributed_graph.adjacent_nodes ),
                                           vwgt, adjwgt, &wgtflag,
                                           &numflag, &ncon, &nparts, tpwgts, &ubvec, options,
                                           &num_edge_cuts, my_part,
                                           &comm_copy );
            assert ( metis_return == METIS_OK );

            delete [] tpwgts;

            // gather partitioning on master proc
            std::vector<int> recv_cts;
            if ( rk == master_rank )
            {
                recv_cts.resize ( sz );
                for ( unsigned i = 0; i < sz; ++i )
                    recv_cts[i] = vtxdist[i + 1] - vtxdist[i];
                partitioning->resize ( num_nodes );
            }

            MPI_Gatherv ( my_part, vtxdist[rk + 1] - vtxdist[rk], MPI_INT,
                          vec2ptr ( *partitioning ),
                          vec2ptr ( recv_cts ),
                          vec2ptr ( vtxdist ),
                          MPI_INT, master_rank, comm );

            delete [] my_part;

            LOG_INFO ( "ParMetisGraphPartitioner", "partitioning done" );
        }

        void ParMetisMeshPartitioner::partition_mesh ( const mesh::MeshPtr master_mesh,
                                                       const int master_rank,
                                                       const std::vector<int>* nodal_weights,
                                                       const std::vector<int>* edge_weights,
                                                       const MPI_Comm& comm,
                                                       std::vector<int>* partitioning ) const
        {
            LOG_DEBUG ( 1, "Not yet implemented." );
            quit_program ( );

            LOG_INFO ( "ParMetisMeshPartitioner", "partitioning done" );
        }
#endif // WITH_PARMETIS
    }
} // namespace hiflow
