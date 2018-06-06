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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "hiflow.h"

/// author Staffan Ronnas

/// \brief This program demonstrates refining a mesh using an
/// attribute to specify the type of refinement for each cell.

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 2;
const GDim gdim = 2;

const int DEBUG_LEVEL = 1;
const int NUM_REFINED_CELLS = 1;

static const char* datadir = MESHES_DATADIR;

int find_edge_in_cell ( const Entity& edge, const Entity& cell );

int main ( int, char** )
{
    //refinement of quads
    std::string filename = std::string ( datadir ) + std::string ( "three_quads_L_2d.vtu" );

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    ScopedPtr<Reader>::Type reader ( new VtkReader ( mb ) );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    const EntityCount num_cells = mesh->num_entities ( tdim );
    AttributePtr ref_attr ( new IntAttribute ( std::vector<int>( num_cells, -1 ) ) );
    mesh->add_attribute ( "refinement_type", tdim, ref_attr );

    // refine cell 0 isotropically
    mesh->begin ( tdim )->set ( "refinement_type", 0 );

    // correct other cells
    for ( EntityIterator cell_it = mesh->begin ( tdim ); cell_it != mesh->end ( tdim ); ++cell_it )
    {
        int ref_type;
        cell_it->get ( "refinement_type", &ref_type );

        if ( ref_type == 0 )
        {
            for ( IncidentEntityIterator face_it = cell_it->begin_incident ( tdim - 1 ); face_it != cell_it->end_incident ( tdim - 1 ); ++face_it )
            {
                for ( IncidentEntityIterator neighbor_it = face_it->begin_incident ( tdim ); neighbor_it != face_it->end_incident ( tdim );
                      ++neighbor_it )
                {
                    int neighbor_ref_type;
                    neighbor_it->get ( "refinement_type", &neighbor_ref_type );
                    if ( neighbor_ref_type == -1 )
                    {
                        const int edge_number = find_edge_in_cell ( *face_it, *neighbor_it );
                        if ( edge_number % 2 == 1 )
                        {
                            neighbor_it->set ( "refinement_type", 1 );
                        }
                        else
                        {
                            neighbor_it->set ( "refinement_type", 2 );
                        }
                    }
                }
            }
        }
    }

    MeshPtr refined_mesh = mesh->refine ( "refinement_type" );

    ScopedPtr<Writer>::Type writer ( new VtkWriter ( ) );
    writer->write ( "three_quads_L_2d_refined_by_attribute.vtu", *refined_mesh );

    return 0;
}

int find_edge_in_cell ( const Entity& edge, const Entity& cell )
{
    SortedArray<Id> edge_verts = SortedArray<Id>( std::vector<Id>( edge.begin_vertex_ids ( ), edge.end_vertex_ids ( ) ) );
    int k = 0;
    for ( IncidentEntityIterator edge_it = cell.begin_incident ( tdim - 1 ); edge_it != cell.end_incident ( tdim - 1 ); ++edge_it )
    {
        bool found = true;
        for ( VertexIdIterator v_it = edge_it->begin_vertex_ids ( ); v_it != edge_it->end_vertex_ids ( ); ++v_it )
        {
            if ( !edge_verts.find ( *v_it, 0 ) )
            {
                found = false;
            }
        }
        if ( found )
        {
            return k;
        }
        ++k;
    }
    return false;
}
