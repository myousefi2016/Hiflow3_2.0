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

/// \brief This program measures the time it takes to perform uniform
/// refinement of a hexahedral mesh.

/// \author Staffan Ronnas

#include "hiflow.h"

using namespace hiflow::mesh;
using hiflow::Timer;

static const char* datadir = MESHES_DATADIR;

const TDim tdim = 3;
const GDim gdim = 3;

const int max_level = 6;

int main ( int argc, char *argv[] )
{
    // create Mesh Database
    MeshDbViewBuilder mb ( tdim, gdim );

    // create reader
    UcdReader reader ( &mb );

    // read mesh
    std::string filename = std::string ( datadir ) + std::string ( "unit_cube.inp" );
    MeshPtr mesh;
    reader.read ( filename.c_str ( ), mesh );

    // output
    VtkWriter writer;
    writer.write ( "original_cube.vtu", *mesh );

    MeshPtr refined_mesh;

    for ( int level = 0; level < max_level; ++level )
    {
        const EntityCount num_cells = mesh->num_entities ( tdim );
        Timer timer;
        refined_mesh = mesh->refine ( );
        timer.stop ( );

        std::cout << "Refinement to level " << level << " took "
                << timer.get_duration ( ) << " s.\t num refined cells = " << num_cells << std::endl;

        mesh = refined_mesh;
    }

    return 0;
}
