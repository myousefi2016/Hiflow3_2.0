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

#include <tr1/array>

using namespace hiflow;
using namespace hiflow::mesh;

const int tdim = 2;
const int gdim = 2;

bool check_mesh ( const Mesh& mesh );

// Checks the orientation of cells in a mesh file.
// Currently supports only triangles and tetrahedra.

int main ( int argc, char *argv[] )
{
    int read_args = 1;

    if ( argc < 2 )
    {
        std::cerr << "Usage : check_mesh_input filename\n";
        return 1;
    }

    const std::string in_filename ( argv[read_args++] );

    std::cout << "Reading mesh " << in_filename << "...\n";

    MeshPtr mesh = read_mesh_from_file ( in_filename, tdim, gdim, 0 );

    if ( !mesh )
    {
        std::cerr << "Failed to read mesh. Exiting.\n";
        return 2;
    }

    assert ( mesh != 0 );
    const bool mesh_ok = check_mesh ( *mesh );

    if ( mesh_ok )
        std::cout << "Mesh " << in_filename << " is correctly oriented.\n";
    else
        std::cout << "Mesh " << in_filename << " is NOT correctly oriented.\n";

    return mesh_ok ? 0 : -1;
}

bool check_mesh ( const Mesh& mesh )
{
    bool mesh_ok = true;
    for ( EntityIterator it = mesh.begin ( tdim ); it != mesh.end ( tdim ); ++it )
    {
        const CellType& cell_type = it->cell_type ( );
        std::vector<double> coords;
        it->get_coordinates ( coords );
        bool ok = cell_type.check_cell_geometry ( coords, it->gdim ( ) );

        if ( !ok )
        {
            std::cout << "Cell " << it->index ( ) << " has incorrect orientation\n";
            std::cout << "\t" << *it << "\n";
            mesh_ok = false;
        }
    }
    return mesh_ok;

}
