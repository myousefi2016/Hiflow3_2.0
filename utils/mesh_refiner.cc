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

#include <fstream>
#include <string>
#include <sstream>

using namespace hiflow;
using namespace hiflow::mesh;

int tdim = 3;
const int gdim = 3;

// Reads a mesh and refines it to the desired level.
// Output in same format as the input.

int main ( int argc, char *argv[] )
{
    int read_args = 1;

    if ( argc < 3 )
    {
        std::cerr << "Usage : mesh_refiner tdim filename level\n";
        return 1;
    }

    tdim = atoi ( argv[read_args++] );
    const std::string in_filename ( argv[read_args++] );
    const int level = atoi ( argv[read_args++] );

    MeshPtr mesh = read_mesh_from_file ( in_filename, tdim, gdim, 0 );
    assert ( mesh != 0 );

    for ( int r = 0; r < level; ++r )
    {
        mesh = mesh->refine ( );
        assert ( mesh != 0 );
    }

    assert ( mesh != 0 );

    bool write_vtk = false;
    std::size_t pos = in_filename.find_last_of ( "." );

    if ( pos != std::string::npos )
    {
        const std::string ending = in_filename.substr ( pos );

        if ( ending == ".vtu" )
        {
            write_vtk = true;
        }
    }

    std::string outfile_base = in_filename.substr ( 0, pos );
    std::ostringstream outfile_sstr;
    outfile_sstr << outfile_base << "_level_" << level;
    if ( write_vtk )
    {
        outfile_sstr << ".vtu";
        VtkWriter writer;
        writer.write ( outfile_sstr.str ( ).c_str ( ), *mesh );
    }
    else
    {
        outfile_sstr << ".inp";
        UcdWriter writer;
        writer.write ( outfile_sstr.str ( ).c_str ( ), *mesh );
    }

    return 0;
}
