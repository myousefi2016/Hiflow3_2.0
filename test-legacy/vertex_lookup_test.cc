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

#include <vector>
#include <iostream>

#include "mesh.h"
#include "common/log.h"
#include "../src/mesh/mesh_database.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;
const EntityCount num_verts = 8;

// Unit cube definition
Coordinate coords[num_verts][gdim] = {
    { 0., 0., 1.e-14 },
    { 0., 0., 0. },
    { 1., 0., -1. },
    { 1., 0., -1. },
    { 0., 1., 0.00001 },
    { -1., 1., 0. },
    { 0., 1., 0. },
    { 0., 0., 0. }
};

int main ( int, char ** )
{
    MeshDatabase db ( tdim, gdim );

    for ( int i = 0; i < num_verts; ++i )
    {
        const Id id = db.add_vertex ( &coords[i][0] );
        std::clog << "Added vertex " << id << ": ";
        std::clog << string_from_pointer_range ( &coords[i][0], &coords[i][gdim] ) << "\n";
    }

    return 0;
}
