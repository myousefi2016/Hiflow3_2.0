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

/// \author Thomas Gengenbach, Staffan Ronnas

#include <iostream>
#include <fstream>

#include "test.h"
#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const int DEBUG_LEVEL = 1;

static const char* datadir = MESH_DATADIR;

int main ( int, char** )
{
    const TDim tdim = 3;
    const GDim gdim = 3;

    std::ofstream debug_file ( "boundary_extraction_test_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );

    LOG_DEBUG ( 1, "Start unit cube test" );
    std::string filename = std::string ( datadir ) + std::string ( "unit_cube.inp" );

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    ScopedPtr<Reader>::Type reader ( new UcdReader ( mb ) );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );
    assert ( mesh != 0 );

    MeshPtr bdy_mesh = mesh->extract_boundary_mesh ( );

    ScopedPtr<Writer>::Type writer ( new VtkWriter ( ) );
    std::string filename_out = std::string ( "boundary_mesh.vtu" );

    writer->write ( filename_out.c_str ( ), *bdy_mesh );

    LogKeeper::get_log ( "debug" ).flush ( );

#if 0 // use this for big mesh
    // refine
    std::vector<Id> refined_cells;
    for ( int i = 0; i < 50; ++i )
    {
        refined_cells.push_back ( i );
    }
    MeshPtr refined_mesh = mesh->refine ( refined_cells );
#endif

    LOG_DEBUG ( 1, "Refining unit cube" );
    MeshPtr refined_mesh = mesh->refine ( );

    TEST_EQUAL ( refined_mesh->num_entities ( tdim ), 8 );
    TEST_EQUAL ( refined_mesh->num_entities ( tdim - 1 ), 36 );

    LOG_DEBUG ( 1, "Extracting boundary mesh of refined unit cube" );
    MeshPtr refined_bdy_mesh = refined_mesh->extract_boundary_mesh ( );

    TEST_EQUAL ( refined_bdy_mesh->num_entities ( tdim - 1 ), 24 );

    assert ( refined_bdy_mesh != 0 );
    writer->write ( "refined_mesh.vtu", *refined_mesh );
    writer->write ( "refined_boundary_mesh.vtu", *refined_bdy_mesh );

    delete mb;
    LOG_DEBUG ( 1, "End unit cube test" );

    //////////////////////////////////////////////////////////////////
    //
    // Elaborate test with hanging nodes and weird refinements copied
    // from software hands on tutorial from Eva Ketelaer.
    //
    const TDim tdim_elaborate = 2;
    const GDim gdim_elaborate = 2;

    LOG_DEBUG ( 1, "Start elaborate test" );

    std::string filename_elaborate = std::string ( datadir ) + std::string ( "L_domain2.inp" );
    MeshBuilder * mb_elaborate ( new MeshDbViewBuilder ( tdim_elaborate, gdim_elaborate ) );
    ScopedPtr<Reader>::Type reader_elaborate ( new UcdReader ( mb_elaborate ) );

    MeshPtr mesh_elaborate;
    reader_elaborate->read ( filename_elaborate.c_str ( ), mesh_elaborate );
    assert ( mesh_elaborate.get ( ) != 0 );
    //step 1
    const EntityCount num_cells = mesh_elaborate->num_entities ( tdim_elaborate );
    AttributePtr ref_attr ( new IntAttribute ( std::vector<int>( num_cells, -1 ) ) );
    mesh_elaborate->add_attribute ( "refinement_type", tdim_elaborate, ref_attr );

    //refine all cells by refinement 1
    for ( EntityIterator cell_it = mesh_elaborate->begin ( tdim_elaborate ); cell_it != mesh_elaborate->end ( tdim_elaborate ); ++cell_it )
    {
        cell_it->set ( "refinement_type", 1 );
    }

    mesh_elaborate = mesh_elaborate->refine ( "refinement_type" );

    //step 2
    //search all cells with coordinate (0,0)
    const EntityCount num_cells2 = mesh_elaborate->num_entities ( tdim_elaborate );
    AttributePtr ref_attr2 ( new IntAttribute ( std::vector<int>( num_cells2, -1 ) ) );
    mesh_elaborate->add_attribute ( "refinement_type", tdim_elaborate, ref_attr2 );

    for ( EntityIterator cell_it = mesh_elaborate->begin ( tdim_elaborate ); cell_it != mesh_elaborate->end ( tdim_elaborate ); ++cell_it )
    {
        for ( IncidentEntityIterator vertex_it = cell_it->begin_incident ( 0 );
              vertex_it != cell_it->end_incident ( 0 ); ++vertex_it )
        {
            std::vector<Coordinate> coord_vertex;
            vertex_it->get_coordinates ( coord_vertex, 0 );
            //std::cout << "Coordinates of vertex with index: " << vertex_it->id() << std::endl;
            if ( coord_vertex[0] == 0. && coord_vertex[1] == 0. )
            {
                cell_it->set ( "refinement_type", 1 );
            }

        }
    }

    mesh_elaborate = mesh_elaborate->refine ( "refinement_type" );

    //step 3
    //search all cells with coordinate (0,0)
    const EntityCount num_cells3 = mesh_elaborate->num_entities ( tdim_elaborate );
    AttributePtr ref_attr3 ( new IntAttribute ( std::vector<int>( num_cells3, -1 ) ) );
    mesh_elaborate->add_attribute ( "refinement_type", tdim_elaborate, ref_attr3 );

    for ( EntityIterator cell_it = mesh_elaborate->begin ( tdim_elaborate ); cell_it != mesh_elaborate->end ( tdim_elaborate ); ++cell_it )
    {
        for ( IncidentEntityIterator vertex_it = cell_it->begin_incident ( 0 );
              vertex_it != cell_it->end_incident ( 0 ); ++vertex_it )
        {
            std::vector<Coordinate> coord_vertex;
            vertex_it->get_coordinates ( coord_vertex, 0 );
            //std::cout << "Coordinates of vertex with index: " << vertex_it->id() << std::endl;
            if ( coord_vertex[0] == 0. && coord_vertex[1] == 0. )
            {
                cell_it->set ( "refinement_type", 1 );
            }

        }
    }

    //refine all cells by refinement 4
    for ( EntityIterator cell_it = mesh_elaborate->begin ( tdim_elaborate ); cell_it != mesh_elaborate->end ( tdim_elaborate ); ++cell_it )
    {
        cell_it->set ( "refinement_type", 4 );
    }
    MeshPtr refined_mesh_elaborate;
    refined_mesh_elaborate = mesh_elaborate->refine ( "refinement_type" );

    assert ( refined_mesh_elaborate != 0 );

    MeshPtr boundary_mesh_elaborate = refined_mesh_elaborate->extract_boundary_mesh ( );

    assert ( boundary_mesh_elaborate != 0 );

    ScopedPtr<Writer>::Type writer_elaborate ( new VtkWriter ( ) );
    writer_elaborate->write ( "Ldomain2_refined_by_attribute_step3.vtu", *refined_mesh_elaborate );
    writer_elaborate->write ( "boundary_mesh_refined_by_attribute.vtu", *boundary_mesh_elaborate );

    TEST_EQUAL ( refined_mesh_elaborate->num_entities ( tdim_elaborate ), 42 );
    // Hanging nodes and so on...
    TEST_EQUAL ( refined_mesh_elaborate->num_entities ( tdim_elaborate - 1 ), 81 );

    TEST_EQUAL ( boundary_mesh_elaborate->num_entities ( tdim_elaborate - 1 ), 18 );

    LOG_DEBUG ( 1, "End elaborate test" );

    delete mb_elaborate;
}
