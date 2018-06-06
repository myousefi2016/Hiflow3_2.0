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

#include "hiflow.h"

#include "mesh/refinement.h"

#include "test.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;
const int NUM_REFINEMENTS = 4;

const int DEBUG_LEVEL = 1;

static const char* datadir = MESH_DATADIR;

int main ( int, char** )
{
    std::ofstream debug_file ( "refinement_test_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );

    std::string filename = std::string ( datadir ) + std::string ( "unit_cube.inp" );

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    ScopedPtr<Reader>::Type reader ( new UcdReader ( mb ) );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    LOG_DEBUG ( 1, " === Creating refined mesh 1" );

    // first uniform refinement
    MeshPtr refined_mesh = mesh;
    for ( int i = 0; i < NUM_REFINEMENTS; ++i )
    {
        refined_mesh = refined_mesh->refine ( );
    }

    const int expected_num_new_cells = static_cast < int > ( std::pow ( static_cast < double > ( 8 ), NUM_REFINEMENTS ) );
    TEST_EQUAL ( refined_mesh->num_entities ( refined_mesh->tdim ( ) ),
                 expected_num_new_cells );
    LogKeeper::get_log ( "debug" ).flush ( );
    assert ( refined_mesh != 0 );

    ScopedPtr<Writer>::Type writer ( new VtkWriter ( ) );
    std::string filename_before = std::string ( "mesh_before_refinement.vtu" );
    std::string filename_after = std::string ( "mesh_after_refinement.vtu" );
    writer->write ( filename_before.c_str ( ), *mesh );
    writer->write ( filename_after.c_str ( ), *refined_mesh );

    LOG_DEBUG ( 1, " === Creating refined mesh 2" );
    LogKeeper::get_log ( "debug" ).flush ( );

    int new_num_cells = refined_mesh->num_entities ( refined_mesh->tdim ( ) );

    // then refine half of the cells (from the refined mesh)
    std::vector<int> refinement_types ( new_num_cells, -2 );
    int num_ref_cells = new_num_cells / 2;
    for ( int i = 0; i < num_ref_cells; ++i )
    {
        refinement_types[i] = 0;
    }

    refined_mesh = refined_mesh->refine ( refinement_types );

    const int expected_num_new_cells2 = ( new_num_cells - num_ref_cells ) + num_ref_cells / 8;
    TEST_EQUAL ( refined_mesh->num_entities ( refined_mesh->tdim ( ) ),
                 expected_num_new_cells2 );

    std::string filename_after2 = std::string ( "mesh_after_refinement2.vtu" );
    writer->write ( filename_after2.c_str ( ), *refined_mesh );
    delete mb;
    return 0;
}
