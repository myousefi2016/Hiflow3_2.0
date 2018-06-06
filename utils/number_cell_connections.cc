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

#include "hiflow.h"

#include <cstdlib>
#include <tr1/array>

using namespace hiflow;
using namespace hiflow::mesh;

// Computes the size of the connectivity of each cell.

int main ( int argc, char *argv[] )
{
    int read_args = 1;

    if ( argc < 3 || argc > 4 )
    {
        std::cerr << "Usage : check_mesh_input filename tdim [gdim]\n";
        return 1;
    }

    const std::string in_filename ( argv[read_args++] );

    int tdim = atoi ( argv[read_args++] );

    assert ( tdim >= 1 );
    assert ( tdim <= 3 );

    int gdim = tdim;
    if ( argc == 4 )
    {
        gdim = atoi ( argv[read_args++] );
    }

    std::cout << "Reading mesh " << in_filename << "...\n";

    MeshPtr mesh = read_mesh_from_file ( in_filename, tdim, gdim, 0 );

    if ( !mesh )
    {
        std::cerr << "Failed to read mesh. Exiting.\n";
        return 2;
    }

    assert ( mesh != 0 );

    int min_conn = mesh->num_entities ( tdim ), max_conn = 0, total = 0;
    for ( EntityIterator it = mesh->begin ( tdim ), end_it = mesh->end ( tdim ); it != end_it; ++it )
    {
        const int num_conn = it->num_incident_entities ( tdim );
        //      std::cout << "Cell " << it->index() << " is connected to " << num_conn << " cells.\n";
        min_conn = std::min ( min_conn, num_conn );
        max_conn = std::max ( max_conn, num_conn );
        total += num_conn;
    }

    std::cout << "\nmin conn = " << min_conn << ", max conn = " << max_conn << ", mean conn = " << double(total ) / mesh->num_entities ( tdim ) << "\n";

    return 0;
}
