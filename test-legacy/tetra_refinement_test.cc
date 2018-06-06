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

/// author Thomas Gengenbach

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "test.h"
#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;
const int NUM_REFINED_CELLS_TYPE0 = 1;
const int NUM_REFINED_CELLS_TYPE1 = 1;
const int NUM_REFINED_CELLS = NUM_REFINED_CELLS_TYPE0 + NUM_REFINED_CELLS_TYPE1;
const int DEBUG_LEVEL = 1;

static const char* datadir = MESH_DATADIR;

int main ( int, char** )
{
    std::ofstream debug_file ( "tetra_refinement_test_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );

    std::string filename = std::string ( datadir ) + std::string ( "2tetras.bdf.inp" );

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );

    ScopedPtr<Reader>::Type reader ( new UcdReader ( mb ) );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    const EntityCount num_cells = mesh->num_entities ( mesh->tdim ( ) );
    assert ( NUM_REFINED_CELLS <= num_cells );
    std::vector<int> refinement_types ( num_cells, -2 );
    for ( int i = 0; i < NUM_REFINED_CELLS_TYPE0; ++i )
    {
        refinement_types[i] = 1;
    }
    int offset = NUM_REFINED_CELLS_TYPE0;
    for ( int i = 0; i < NUM_REFINED_CELLS_TYPE1; ++i )
    {
        refinement_types[i + offset] = 2;
    }

    MeshPtr refined_mesh = mesh->refine ( refinement_types );

    LogKeeper::get_log ( "debug" ).flush ( );
    assert ( refined_mesh != 0 );

    ScopedPtr<Writer>::Type writer ( new VtkWriter ( ) );
    std::string filename_before = std::string ( "tetra_mesh_before_refinement.vtu" );
    std::string filename_after = std::string ( "tetra_mesh_after_refinement.vtu" );
    writer->write ( filename_before.c_str ( ), *mesh );
    writer->write ( filename_after.c_str ( ), *refined_mesh );

    const int expected_num_new_cells = 8 * NUM_REFINED_CELLS_TYPE0 + 4 * NUM_REFINED_CELLS_TYPE1;
    const int expected_num_refined_cells = ( num_cells - NUM_REFINED_CELLS ) + expected_num_new_cells;
    TEST_EQUAL ( refined_mesh->num_entities ( refined_mesh->tdim ( ) ),
                 expected_num_refined_cells );

    delete mb;

    return 0;
}
