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
#include <string>

#include "mesh.h"

using namespace hiflow::mesh;

// TODO(Thomas): Implement conversion of lines too.

int main ( int argc, char *argv[] )
{
    int read_args = 1;

    if ( argc < 2 )
    {
        std::cerr << "Usage : inp2vtu [dimension] filename\n";
        return 1;
    }

    if ( argc == 2 )
        std::cout << "Dimension 3 assumed!\n";

    TDim tdim = 3;
    if ( argc == 3 )
    {
        tdim = atoi ( argv[read_args++] );
    }

    // gdim == tdim is assumed
    MeshDbViewBuilder mb ( tdim, tdim );

    UcdReader reader ( &mb );
    MeshPtr mesh;
    std::string in_filename ( argv[read_args] );
    std::clog << "Reading file " << in_filename << "\n";
    reader.read ( in_filename.c_str ( ), mesh );

    VtkWriter writer;
    int point_pos = in_filename.find_last_of ( '.' );
    std::string out_filename = in_filename.substr ( 0, point_pos ) + std::string ( ".vtu" );
    std::clog << "Writing file " << out_filename << "\n";
    writer.write ( out_filename.c_str ( ), *mesh );

    std::string boundary_out_filename = in_filename.substr ( 0, point_pos ) + std::string ( "_bdy.vtu" );
    std::clog << "Writing boundary file " << boundary_out_filename << "\n";
    MeshPtr bdy = mesh->extract_boundary_mesh ( );
    VtkWriter bdy_writer;
    bdy_writer.write ( boundary_out_filename.c_str ( ), *bdy );

    return 0;
}
