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

    std::string filename = std::string ( datadir ) + std::string ( "unitsquare_8_cells_0.vtu" );

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    ScopedPtr<Reader>::Type reader ( new VtkReader ( mb ) );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    std::cout << "iteration over cells of mesh:\n";
    for ( EntityIterator it = mesh->begin ( tdim ); it != mesh->end ( tdim ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    std::cout << "iteration over faces of mesh:\n";
    for ( EntityIterator it = mesh->begin ( tdim - 1 ); it != mesh->end ( tdim - 1 ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    std::cout << "iteration over edges of mesh:\n";
    for ( EntityIterator it = mesh->begin ( tdim - 2 ); it != mesh->end ( tdim - 2 ); ++it )
    {
        std::cout << *it << std::endl;
    }
    std::cout << "\n\n";

    delete mb;

    // flush log here to avoid problems
    LogKeeper::get_log ( "debug" ).flush ( );
    return 0;
}
