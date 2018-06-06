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

#ifndef HIFLOW_MESH_PARTITIONING_H
#    define HIFLOW_MESH_PARTITIONING_H

/// \brief Classes and functions related to mesh partitioning
/// \author Staffan Ronnas

#    include <vector>

#    include "mesh/types.h"

#    include "mpi.h"

namespace hiflow
{
    namespace mesh
    {

        class Mesh;

        /// \brief Adjacency data structure describing a sparse graph
        ///
        /// \details Stores a sparse graph G = (V,E) of n nodes, where the nodes are implicitly
        /// numbered V = {0,1,...,n-1}. The nodes connected to node i (0 <= i < n) are stored in the
        /// subarray adjacent_nodes[node_offsets[i]],...,adjacent_nodes[node_offsets[i+1]-1]

        struct Graph
        {

            int num_nodes ( ) const
            {
                return node_offsets.size ( ) - 1;
            }

            void clear ( )
            {
                node_offsets.clear ( );
                adjacent_nodes.clear ( );
            }
            std::vector<int> node_offsets;
            std::vector<int> adjacent_nodes;
            std::vector<int> adjacent_weights;
        };

        /// \brief Compute the dual graph of a mesh
        void compute_dual_graph ( const Mesh& mesh, Graph* graph );
        void compute_dual_graph_facet ( const Mesh& mesh, Graph* graph );

        /// \brief Naively distributes a graph from master proc to all procs in the given MPI communicator.
        void distribute_graph_naive ( const Graph& in, const int master_rank, const MPI_Comm& comm, Graph& out );

        struct CellDistribution
        {
            // number of cells in each sub-domain
            std::vector<EntityCount> num_cells;

            // cell id:s sorted by sub-domain
            std::vector<Id> cell_indices;
        };

        /// \brief Compute the cell distribution corresponding to a partition
        void compute_cell_distribution ( const int num_partitions,
                                         const std::vector<int>& partitioning,
                                         CellDistribution* distribution );

        /// \brief Compute the facet distribution corresponding to a cell distribution
        void compute_facet_from_cell_distribution ( const Mesh& mesh,
                                                    const CellDistribution& cell_distribution,
                                                    CellDistribution* facet_distribution );

        /// \brief Abstract base class for GraphPartitioner:s.

        class GraphPartitioner
        {
          public:

            virtual ~GraphPartitioner ( )
            {
            }
            virtual void partition_graph ( const Graph& graph,
                                           const int master_rank,
                                           const std::vector<int>* nodal_weights,
                                           const std::vector<int>* edge_weights,
                                           const MPI_Comm& comm,
                                           std::vector<int>* partitioning ) const = 0;

            virtual bool is_collective ( ) const = 0;
        };

        /// \brief Partitioner which gives approximately that same number
        /// of cells to each partition, in the same order as the cell numbering.
        /// The connections between cells are not taken into account.

        class NaiveGraphPartitioner : public GraphPartitioner
        {
          public:
            virtual void partition_graph ( const Graph& graph,
                                           const int master_rank,
                                           const std::vector<int>* nodal_weights,
                                           const std::vector<int>* edge_weights,
                                           const MPI_Comm& comm,
                                           std::vector<int>* partitioning ) const;

            virtual bool is_collective ( ) const
            {
                return false;
            }
        };

        /// \brief Sequential partitioning using the METIS library

        class MetisGraphPartitioner : public GraphPartitioner
        {
          public:
            virtual void partition_graph ( const Graph& graph,
                                           const int master_rank,
                                           const std::vector<int>* nodal_weights,
                                           const std::vector<int>* edge_weights,
                                           const MPI_Comm& comm,
                                           std::vector<int>* partitioning ) const;

            virtual bool is_collective ( ) const
            {
                return false;
            }
        };

        /// \brief Parallel partitioning using the ParMETIS library. This is a collective routine w.r.t.
        /// the given MPI communicator.

        class ParMetisGraphPartitioner : public GraphPartitioner
        {
          public:
            virtual void partition_graph ( const Graph& graph,
                                           const int master_rank,
                                           const std::vector<int>* nodal_weights,
                                           const std::vector<int>* edge_weights,
                                           const MPI_Comm& comm,
                                           std::vector<int>* partitioning ) const;

            virtual bool is_collective ( ) const
            {
                return true;
            }
        };

        class MeshPartitioner
        {
          public:

            virtual ~MeshPartitioner ( )
            {
            }
            virtual void partition_mesh ( const mesh::MeshPtr master_mesh,
                                          const int master_rank,
                                          const std::vector<int>* nodal_weights,
                                          const std::vector<int>* edge_weights,
                                          const MPI_Comm& comm,
                                          std::vector<int>* partitioning ) const = 0;

            virtual bool is_collective ( ) const = 0;
        };

        class ParMetisMeshPartitioner : public MeshPartitioner
        {
          public:
            virtual void partition_mesh ( const mesh::MeshPtr master_mesh,
                                          const int master_rank,
                                          const std::vector<int>* nodal_weights,
                                          const std::vector<int>* edge_weights,
                                          const MPI_Comm& comm,
                                          std::vector<int>* partitioning ) const;

            virtual bool is_collective ( ) const
            {
                return true;
            }
        };
    }
} // namespace hiflow

#endif
