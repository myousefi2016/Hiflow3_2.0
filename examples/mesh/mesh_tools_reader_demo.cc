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

/// \author Staffan Ronnas, Thomas Gengenbach

/// \brief This program demonstrates how to read a mesh using the
/// functionality in the "mesh toolbox". The function
/// read_mesh_from_file() is used to read in the mesh given only its
/// dimension and filename. The type of file is determined from the
/// filename suffix.

#include <fstream>
#include <string>
#include <mpi.h>

#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

// set topological and geometrical dimension of meshes
const mesh::TDim tdim = 3;
const mesh::GDim gdim = 3;

static const char* datadir = MESHES_DATADIR;

int main ( int argc, char** argv )
{
    // setup mpi for parallel vtk reader
    MPI_Init ( &argc, &argv );

    // debug output
    std::ofstream debug_file ( "dbg_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_file );

    // read vtk file
    std::string filename_vtk = std::string ( datadir ) + std::string ( "unitcube_refinementlevel_3.vtu" );
    MeshPtr mesh_vtk = read_mesh_from_file ( filename_vtk, tdim, gdim, 0 );
    assert ( mesh_vtk != 0 );

    // read ucd (inp) file
    std::string filename_ucd = std::string ( datadir ) + std::string ( "unit_cube.inp" );
    MeshPtr mesh_ucd = read_mesh_from_file ( filename_ucd, tdim, gdim, 0 );
    assert ( mesh_ucd != 0 );

    // read parallel vtk file
    std::string filename_pvtk = std::string ( datadir ) + std::string ( "mat_test_unitsquare_8_cells.pvtu" );
    // to test the pvtk reader, run with:
    // mpirun -np 4 mesh_tools_reader_demo
    MPI_Comm comm = MPI_COMM_WORLD;
    MeshPtr mesh_pvtk = read_mesh_from_file ( filename_pvtk, tdim, gdim, &comm );
    assert ( mesh_pvtk != 0 );

    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );
    MPI_Finalize ( );
    return 0;
}
