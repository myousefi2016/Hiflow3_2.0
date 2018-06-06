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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESH_DATADIR;

int main ( int, char** )
{
    std::ofstream debug_file ( "dbg_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_file );

    std::string filename_vtk = std::string ( datadir ) + std::string ( "unitcube_refinementlevel_3.vtu" );
    std::string filename_ucd = std::string ( datadir ) + std::string ( "unit_cube.inp" );

    MeshBuilder * mb_vtk ( new MeshDbViewBuilder ( tdim, gdim ) );
    MeshBuilder * mb_ucd ( new MeshDbViewBuilder ( tdim, gdim ) );

    ScopedPtr<Reader>::Type reader_vtk ( new VtkReader ( mb_vtk ) );
    ScopedPtr<Reader>::Type reader_ucd ( new UcdReader ( mb_ucd ) );

    MeshPtr mesh_vtk;
    MeshPtr mesh_ucd;

    reader_vtk->read ( filename_vtk.c_str ( ), mesh_vtk );
    reader_ucd->read ( filename_ucd.c_str ( ), mesh_ucd );

    std::cout << "iteration over cells of mesh_ucd:\n";
    for ( EntityIterator it = mesh_ucd->begin ( tdim ); it != mesh_ucd->end ( tdim ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    //    mesh_ucd->compute_entities(2);
    std::cout << "iteration over faces of mesh_ucd:\n";
    for ( EntityIterator it = mesh_ucd->begin ( tdim - 1 ); it != mesh_ucd->end ( tdim - 1 ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    //    mesh_ucd->compute_entities(1);
    std::cout << "iteration over edges of mesh_ucd:\n";
    for ( EntityIterator it = mesh_ucd->begin ( tdim - 2 ); it != mesh_ucd->end ( tdim - 2 ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    std::cout << "iteration over cells of mesh_vtk:\n";
    for ( EntityIterator it = mesh_vtk->begin ( tdim ); it != mesh_vtk->end ( tdim ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    //    mesh_vtk->compute_entities(2);
    std::cout << "iteration over faces of mesh_vtk:\n";
    for ( EntityIterator it = mesh_vtk->begin ( tdim - 1 ); it != mesh_vtk->end ( tdim - 1 ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    //    mesh_vtk->compute_entities(1);
    std::cout << "iteration over edges of mesh_vtk:\n";
    for ( EntityIterator it = mesh_vtk->begin ( tdim - 2 ); it != mesh_vtk->end ( tdim - 2 ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    delete mb_vtk;
    delete mb_ucd;

    std::cout << "Test passed\n";
    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );
    return 0;
}
