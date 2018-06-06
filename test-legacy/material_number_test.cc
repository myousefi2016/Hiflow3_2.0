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
#include <mpi.h>

#include "test.h"
#include "hiflow.h"

using namespace hiflow;
using namespace hiflow::mesh;

static const char* datadir = MESH_DATADIR;

const TDim tdim = 3;
const GDim gdim = 3;

int main ( int argc, char* argv[] )
{

    MPI_Init ( &argc, &argv );
    std::string filename = std::string ( datadir ) + std::string ( "mat_test_unitsquare_8_cells.pvtu" );

    MeshDbViewBuilder mb ( tdim, gdim );

    MeshPtr mesh;

    ScopedPtr<Reader>::Type reader ( new PVtkReader ( &mb, MPI_COMM_WORLD ) );
    reader->read ( filename.c_str ( ), mesh );

    // iterate over all entities of tdim = 3 and all boundary facets
    // getting the material numbers
    for ( EntityIterator it = mesh->begin ( tdim ); it != mesh->end ( tdim ); ++it )
    {
        MaterialNumber mat_num = it->get_material_number ( );
        //        std::cout << "Material Number of entity id " << it->id()
        //                  << " is " << mat_num << ".\n";
        TEST_EQUAL ( mat_num, 1234 );
    }

    MPI_Finalize ( );

    return 0;
}
