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

#include <boost/iterator/filter_iterator.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <vector>
#include <mpi.h>

#include "hiflow.h"

/// \author Staffan Ronnas

/// \brief This program demonstrates the partitioning of a mesh with Metis. The
/// mesh is read in on all process, and a partitioning is then computed using
/// the Metis graph partitioner on the master process, and communicated to all
/// processes. Note that the mesh itself is not distributed in this example.
/// NB: This program must run with at least two processes, otherwise Metis
/// crashes.

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESHES_DATADIR;

const int master_rank = 0;
const int DEBUG_LEVEL = 1;

struct PartitionFilter
{

    bool operator() ( const Entity& entity ) const
    {
        return partitioning[entity.index ( )] == rank;
    }

    std::vector<int> partitioning;
    int rank;
};

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );
    int num_partitions = -1;
    int rank = -1;

    MPI_Comm_size ( MPI_COMM_WORLD, &num_partitions );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    ostringstream rank_str;
    rank_str << rank;

    std::string debug_filename = string ( "dbg_output." ) + rank_str.str ( ) + string ( ".log" );
    std::ofstream debug_file ( debug_filename.c_str ( ) );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );

    LOG_DEBUG ( 1, "Starting global_partitioning_demo on process " << rank << " of " << num_partitions );

    std::string filename_in = std::string ( datadir ) + std::string ( "dfg_bench3d_cyl.inp" );

    MeshDbViewBuilder builder ( tdim, gdim );

    LOG_DEBUG ( 1, "Reading mesh from " << filename_in << " on process " << rank );
    //        VtkReader reader(&builder);
    UcdReader reader ( &builder );

    MeshPtr mesh;
    reader.read ( filename_in.c_str ( ), mesh );
    assert ( mesh != 0 );

    LOG_DEBUG ( 1, "Read mesh with " << mesh->num_entities ( tdim ) << " cells on process " << rank );

    PartitionFilter filter;
    filter.rank = rank;
    if ( rank == master_rank )
    {
        LOG_DEBUG ( 1, "Partitioning mesh on process " << rank );
        Graph graph;
        compute_dual_graph ( *mesh, &graph );
        MetisGraphPartitioner metis;
        metis.partition_graph ( graph, master_rank, 0, 0, MPI_COMM_WORLD, &filter.partitioning );
        assert ( static_cast < int > ( filter.partitioning.size ( ) ) == mesh->num_entities ( tdim ) );
    }
    else
    {
        filter.partitioning.resize ( mesh->num_entities ( tdim ) );
    }
    MPI_Bcast ( vec2ptr ( filter.partitioning ), mesh->num_entities ( tdim ), MPI_INT, master_rank, MPI_COMM_WORLD );

    typedef boost::filter_iterator<PartitionFilter, EntityIterator> PartitionCellIterator;
    PartitionCellIterator it ( filter, mesh->begin ( tdim ), mesh->end ( tdim ) );
    PartitionCellIterator end ( filter, mesh->end ( tdim ), mesh->end ( tdim ) );

    for (; it != end; ++it )
    {
        std::cout << "cell " << it->index ( ) << " belongs to process " << rank << "\n";
    }

    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );

    MPI_Finalize ( );

    return 0;
}
