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

/// Program that reads in a mesh on the master process, partitions it into as
/// many subdomains as there are processes, and writes the result as a .pvtu
/// file. The boundary is written as a separate .vtu file, which can be read in
/// on each process.

/// \author Staffan Ronnas

#include <iostream>
#include <string>

#include "hiflow.h"

using namespace hiflow;
using namespace hiflow::mesh;

const int MASTER_RANK = 0;

#define MASTER_WRITE(x) if (rank == MASTER_RANK) { std::cout << x; }

int main ( int argc, char *argv[] )
{
    MPI_Init ( &argc, &argv );

    int rank, num_proc;

    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &num_proc );

    int read_args = 1;

    if ( argc < 2 )
    {
        MASTER_WRITE ( "Usage : mpirun -np P " << argv[0] << " [dimension] filename\n" );
        return 1;
    }

    if ( argc == 2 )
    {
        MASTER_WRITE ( "Dimension 3 assumed!\n" );
    }

    TDim tdim = 3;
    if ( argc == 3 )
    {
        tdim = atoi ( argv[read_args++] );
    }

    MeshPtr in_mesh;

    std::string in_filename ( argv[read_args] );

    // Read mesh from file on master rank.
    if ( rank == MASTER_RANK )
    {
        // gdim == tdim is assumed
        MASTER_WRITE ( "Reading file " << in_filename << " on process " << MASTER_RANK << "\n" );
        in_mesh = read_mesh_from_file ( in_filename, tdim, tdim, 0 );

        // Write boundary mesh.
        int point_pos = in_filename.find_last_of ( '.' );
        std::string boundary_out_filename = in_filename.substr ( 0, point_pos ) + std::string ( "_bdy.vtu" );
        MASTER_WRITE ( "Writing boundary file " << boundary_out_filename << "\n" );
        MeshPtr bdy = in_mesh->extract_boundary_mesh ( );
        VtkWriter bdy_writer;
        bdy_writer.write ( boundary_out_filename.c_str ( ), *bdy );
    }

    // Partition and distribute mesh over processes.
#ifdef WITH_METIS
    MetisGraphPartitioner partitioner;
    MASTER_WRITE ( "Partitioning with METIS partitioner\n" );
#else
    NaiveGraphPartitioner partitioner;
    MASTER_WRITE ( "Partitioning with naive partitioner\n" );
#endif
    MeshPtr local_mesh = partition_and_distribute ( in_mesh, MASTER_RANK, MPI_COMM_WORLD, &partitioner );
    assert ( local_mesh != 0 );

    // Write distributed mesh in Parallel Vtk format.
    PVtkWriter writer ( MPI_COMM_WORLD );
    int point_pos = in_filename.find_last_of ( '.' );
    std::string out_filename = in_filename.substr ( 0, point_pos ) + std::string ( ".pvtu" );
    MASTER_WRITE ( "Writing parallel VTK file " << out_filename << "\n" );
    writer.write ( out_filename.c_str ( ), *local_mesh );

    MPI_Finalize ( );
    return 0;
}
