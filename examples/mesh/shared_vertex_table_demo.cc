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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "hiflow.h"

/// \author Staffan Ronnas

/// \brief This program demonstrates the computation of a SharedVertexTable object, that can be used for communication between processes.

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESHES_DATADIR;

int main ( int argc, char** argv )
{
    // TODO: make this into a real test
    MPI_Init ( &argc, &argv );

    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    std::ofstream debug_file ( "dbg_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_file );
    std::string filename = std::string ( datadir ) + std::string ( "unitsquare_8_cells.pvtu" );
    //    std::string filename = std::string(datadir) + std::string("lunge_fein_3.pvtu");

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );

    ScopedPtr<Reader>::Type reader ( new PVtkReader ( mb, MPI_COMM_WORLD ) );

    MPI_Barrier ( MPI_COMM_WORLD );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    std::cout << "Read mesh with " << mesh->num_entities ( 0 )
            << " vertices on proc " << rank << "\n";

    MeshPtr bdy_mesh = mesh->extract_boundary_mesh ( );

    SharedVertexTable shared_verts;
    update_shared_vertex_table ( *mesh, MPI_COMM_WORLD, shared_verts );

    std::cout << "shared_verts computed with full mesh...\n";

    // we should get the same results, by only testing boundary vertices
    SharedVertexTable shared_verts_bdy;
    update_shared_vertex_table ( *bdy_mesh, MPI_COMM_WORLD, shared_verts_bdy );
    std::cout << "shared_verts computed with bdy mesh...\n";

    // check that we do get the same results
    for ( SharedVertexTable::const_iterator it = shared_verts.begin ( ); it != shared_verts.end ( ); ++it )
    {
        assert ( shared_verts_bdy.has_vertex ( it->first ) );
        assert ( shared_verts_bdy.num_shared_vertices ( it->first ) ==
                 shared_verts.num_shared_vertices ( it->first ) );
    }

    std::cout << "same results!\n";

    std::cout << "On proc " << rank
            << "\n" << shared_verts << "\n";

    delete mb;

    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );

    MPI_Finalize ( );
    return 0;
}
