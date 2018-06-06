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

/// \brief This program demonstrates how to partition and distribute a
/// mesh using the functionality in the "mesh toolbox". The function
/// partition_and_distribute is used to compute a partition of a mesh
/// read in on one process, and then distribute the parts to the other
/// processes.

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <vector>
#include <mpi.h>

#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESHES_DATADIR;

const int master_rank = 0;

const int DEBUG_LEVEL = 3;

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

    LOG_DEBUG ( 1, "Starting mesh_tools_partitioning_demo on process " << rank << " of " << num_partitions );

    MeshPtr master_mesh;

    if ( rank == master_rank )
    {
        std::string filename = std::string ( datadir ) + std::string ( "dfg_bench3d_cyl.inp" );

        LOG_DEBUG ( 1, "Reading mesh from " << filename << " on process " << rank );
        master_mesh = read_mesh_from_file ( filename, tdim, gdim, 0 );
        LOG_DEBUG ( 1, "Mesh has " << master_mesh->num_entities ( tdim ) << " cells" );
    }

    MeshPtr mesh_part = partition_and_distribute ( master_mesh, master_rank, MPI_COMM_WORLD );

    LOG_DEBUG ( 1, "Writing mesh part on process " << rank );
    VtkWriter writer;
    std::string filename_out = std::string ( "mesh_part_" ) + rank_str.str ( ) + std::string ( ".vtu" );
    writer.write ( filename_out.c_str ( ), *mesh_part );

    LOG_DEBUG ( 1, "Writing boundary of mesh part on process " << rank );
    std::string bdy_filename_out = std::string ( "bdy_mesh_part_" ) + rank_str.str ( ) + std::string ( ".vtu" );
    MeshPtr bdy_mesh_part = mesh_part->extract_boundary_mesh ( );
    writer.write ( bdy_filename_out.c_str ( ), *bdy_mesh_part );

    LOG_DEBUG ( 1, "Done on " << rank );

    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );

    MPI_Finalize ( );

    return 0;
}
