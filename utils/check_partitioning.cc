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

/// \author Thomas Gengenbach
#include <iomanip>
#include <fstream>
#include <string>

#include <mpi.h>

#include "mesh.h"
#include "dof.h"
#include "fem.h"
#include "common.h"
#include "space.h"
#include "config.h"
#include "tools.h"

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

// NB: if 2d mesh is checked, this has to be adapted!
const int tdim = 3;
const int gdim = 3;

bool check_partitioning ( MeshPtr mesh, int option, int num_partitions );

int main ( int argc, char *argv[] )
{
    MPI_Init ( &argc, &argv );
    int read_args = 1;

    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    // usage help
    if ( argc < 3 )
    {
        if ( rank == 0 )
        {
            std::cerr << "Usage : mpirun -n <num_partitions> ./check_partitioning <filename> <option>\n";
            std::cerr << "option = 0 -> NAIVE partitioner\n";
            std::cerr << "option = 1 -> METIS partitioner\n";
        }
        return 1;
    }
    // read in args
    const std::string in_filename ( argv[read_args++] );
    const int option ( std::atoi ( argv[read_args++] ) );

    // set num_partitions from mpirun argument
    int num_partitions = -1;
    MPI_Comm_size ( MPI_COMM_WORLD, &num_partitions );

    // set info log output
    std::string info_file = "info_check_partitioning.log";
    std::ofstream info_log ( info_file.c_str ( ) );
    LogKeeper::get_log ( "info" ).set_target ( &info_log );

    // read in mesh for which to check the partition
    MeshPtr mesh;
    if ( rank == 0 )
    {
        mesh = read_mesh_from_file ( in_filename, tdim, gdim, 0 );
        assert ( mesh != 0 );
    }

    // check partitioning
    const bool partitioning_ok = check_partitioning ( mesh, option, num_partitions );

    // output
    if ( partitioning_ok )
        std::cout << "Partitioning on process " << rank << " is defintely correct. (No isolated cells).\n";
    else
        std::cout << "Partitioning on process " << rank << " is probably NOT correct.\n";

    LogKeeper::get_log ( "info" ).set_target ( 0 );
    MPI_Finalize ( );
    return partitioning_ok ? 0 : -1;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

bool check_partitioning ( MeshPtr mesh, int option, int num_partitions )
{
    bool partitioning_ok = true;
    GraphPartitioner* p;

    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

#ifdef WITH_METIS
    if ( option == 0 )
    {
        if ( rank == 0 )
            std::cout << "Partitioning NAIVE...\n";
        p = new NaiveGraphPartitioner ( );
    }
    else
    {
        if ( rank == 0 )
            std::cout << "Partitioning with METIS...\n";
        p = new MetisGraphPartitioner ( );
    }
#else
    if ( option == 1 )
    {
        if ( rank == 0 )
            std::cerr << "METIS partitioner NOT available. Need to recompile Hiflow linking against METIS.\n";
    }
    if ( rank == 0 )
        std::cout << "Partitioning NAIVE...\n";
    p = new NaiveGraphPartitioner ( );
#endif

    // create local meshes
    MeshPtr local_mesh = partition_and_distribute ( mesh, 0, MPI_COMM_WORLD, p );
    assert ( local_mesh != 0 );

    // release memory
    delete p;

    // write local meshes to files
    std::string out_filename;
    std::stringstream ss;
    ss << "mesh_" << rank << ".vtu";
    VtkWriter writer;
    out_filename = ss.str ( );
    writer.write ( out_filename.c_str ( ), *local_mesh );

    // check for neighbor cells (one or less neighbors are a necessary
    // but not sufficient condition)
    //
    // TODO: Write a check, that also checks number of dofs, hence
    // delivers a sufficient condition as well.

    std::stringstream stream;
    for ( EntityIterator it = local_mesh->begin ( tdim ); it != local_mesh->end ( tdim ); ++it )
    {
        // Count boundary facets for each cell, if the number of
        // boundary facets is equal or 1 less than the number of total
        // facets of this cell, than it can be isolated.
        int boundary_facet_counter = 0;
        int number_of_facets = it->num_incident_entities ( tdim - 1 );

        for ( IncidentEntityIterator facet_it = it->begin_incident ( tdim - 1 );
              facet_it != it->end_incident ( tdim - 1 ); ++facet_it )
        {
            if ( local_mesh->is_boundary_facet ( facet_it->index ( ) ) )
            {
                boundary_facet_counter++;
            }
        }

        if ( boundary_facet_counter >= number_of_facets - 1 )
        {
            stream << it->index ( ) << " ";
            partitioning_ok = false;
        }
    }

    int num_proc = -1;
    MPI_Comm_size ( MPI_COMM_WORLD, &num_proc );

    for ( int i = 0; i < num_proc; ++i )
    {
        MPI_Barrier ( MPI_COMM_WORLD );
        if ( i == rank )
        {
            std::cerr << "=== PROCESS " << i << " ===\n";
            std::cerr << "Cells with indeces: " << stream.str ( ) << "could be isolated.\n\n";

        }
        MPI_Barrier ( MPI_COMM_WORLD );
    }
    MPI_Barrier ( MPI_COMM_WORLD );
    return partitioning_ok;
}
