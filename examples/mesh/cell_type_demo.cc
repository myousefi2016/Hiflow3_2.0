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

/// \brief This program reads in a mesh containing two tetrahedra and
/// one hexahedra, refines all three cells, and writes out the resulting mesh.
/// \author Thomas Gengenbach

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "hiflow.h"

using namespace std;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;

const int DEBUG_LEVEL = 1;

static const char* datadir = MESHES_DATADIR;

int main ( int, char** )
{
    std::string filename = std::string ( datadir ) + std::string ( "cell_type_demo.inp" );
    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    Reader* reader = new UcdReader ( mb );
    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    MeshPtr refined_mesh = mesh->refine ( );

    std::string filename_out = std::string ( "cell_type_demo_refined.vtu" );
    Writer* writer = new VtkWriter ( );
    writer->write ( filename_out.c_str ( ), *refined_mesh );

    delete reader;
    delete writer;
    delete mb;
}
