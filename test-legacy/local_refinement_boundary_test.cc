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

// Test related to Bug #121

#include "test.h"
#include "hiflow.h"

using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 2;
const GDim gdim = 2;

static const char* datadir = MESH_DATADIR;

int main ( int argc, char *argv[] )
{
    const std::string filename = std::string ( datadir ) + "L_domain.inp";
    MeshPtr initial_mesh = read_mesh_from_file ( filename, tdim, gdim );

    // In Bug 121, this refinement causes segmentation fault.
    std::vector<int> ref ( 3, -1. );
    ref[0] = 0;
    ref[2] = 0;

    MeshPtr refined_mesh = initial_mesh->refine ( ref );

    MeshPtr ref_bdy_mesh = refined_mesh->extract_boundary_mesh ( );

    VtkWriter writer;
    writer.write ( "local_ref_bdy.vtu", *ref_bdy_mesh );

    return 0;
}
