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

/// \brief demo how two refine a mesh by attributes
/// \author Eva Ketelaer
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
    // read start mesh
    std::string filename = std::string ( datadir ) + std::string ( "L_domain2.inp" );

    MeshDbViewBuilder mb ( tdim, gdim );
    ScopedPtr<Reader>::Type reader ( new UcdReader ( &mb ) );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    // step 1 - refine globally
    int refinements_global = 2;

    // mesh = mesh->refine("refinement_type");
    for ( int i = 0; i < refinements_global; ++i )
    {
        mesh = mesh->refine ( );
    }

    // step 2 - find local cell which will be refined and refine them
    int refinements_local = 3;
    for ( int i = 0; i < refinements_local; ++i )
    {
        const EntityCount num_cells = mesh->num_entities ( tdim );
        // ATTENTION: refinement_type -1 stands for coarsening
        // Therefore put any negative number except -1 for example -10
        AttributePtr ref_attr ( new IntAttribute ( std::vector<int>( num_cells, -10 ) ) );
        mesh->add_attribute ( "refinement_type", tdim, ref_attr );

        // step 2a) find all cells with coordinate (0,0)
        // and define the refinement
        // loop over all cell in mesh
        for ( EntityIterator cell_it = mesh->begin ( tdim );
              cell_it != mesh->end ( tdim ); ++cell_it )
        {
            // loop over the vertecies of each cell
            for ( IncidentEntityIterator vertex_it = cell_it->begin_incident ( 0 );
                  vertex_it != cell_it->end_incident ( 0 );
                  ++vertex_it )
            {
                std::vector<Coordinate> coord_vertex;
                vertex_it->get_coordinates ( coord_vertex, 0 );
                // if one cell has a vertex with the coordinates (0,0)
                // and define refinement_type
                if ( coord_vertex[0] == 0. && coord_vertex[1] == 0. )
                {
                    cell_it->set ( "refinement_type", 0 );
                }
            }
        }

        // step 2b) find all the neighbors of the cell found in step 2a
        // to avoid hanging nodes and define the refinement
        // loop over all cell
        for ( EntityIterator cell_it = mesh->begin ( tdim );
              cell_it != mesh->end ( tdim ); ++cell_it )
        {
            // get the refinement_type of the cell
            int ref_type;
            cell_it->get ( "refinement_type", &ref_type );
            if ( ref_type == 0 )
            {
                // loop over faces of the cell
                for ( IncidentEntityIterator face_it = cell_it->begin_incident ( tdim - 1 );
                      face_it != cell_it->end_incident ( tdim - 1 );
                      ++face_it )
                {
                    // loop over neighbor cell of the faces
                    for ( IncidentEntityIterator neighbor_it =
                          face_it->begin_incident ( tdim );
                          neighbor_it !=
                          face_it->end_incident ( tdim );
                          ++neighbor_it )
                    {
                        // get the refinmenttype of the neighbor
                        int neighbor_ref_type;
                        neighbor_it->get ( "refinement_type", &neighbor_ref_type );
                        // if the refinement type differs from 0 find the local
                        // number of the egde
                        if ( neighbor_ref_type != 0 )
                        {
                            const int edge_number = find_edge_in_cell ( *face_it, *neighbor_it );
                            // depending on local edge number
                            // define refinement type for neighbor cell
                            if ( edge_number == 0 )
                            {
                                neighbor_it->set ( "refinement_type", 7 );
                            }
                            else if ( edge_number == 1 )
                            {
                                neighbor_it->set ( "refinement_type", 8 );
                            }
                            else if ( edge_number == 2 )
                            {
                                neighbor_it->set ( "refinement_type", 5 );
                            }
                            else if ( edge_number == 3 )
                            {
                                neighbor_it->set ( "refinement_type", 6 );
                            }
                        }
                    }
                }
            }
        }

        // step 2c) refine mesh depending on refinement type
        mesh = mesh->refine ( "refinement_type" );
    }

    // step 3 - find all quadrilateral cells and depending on quadrant
    // refine them to two triangles
    // this step was necesary because in the beginning it was not
    // possible to use mixed meshes
    const EntityCount num_cells = mesh->num_entities ( tdim );
    // ATTENTION: refinement_type -1 stands for coarsening
    // Therefore put any negative number except -1 for example -10
    AttributePtr ref_attr ( new IntAttribute ( std::vector<int>( num_cells, -10 ) ) );
    mesh->add_attribute ( "refinement_type", tdim, ref_attr );

    // step 3a) find all quadrilateral cells and define refinement typ
    // loop over all cells
    int num_quad = 0;
    for ( EntityIterator cell_it = mesh->begin ( tdim );
          cell_it != mesh->end ( tdim ); ++cell_it )
    {
        // check if type of cell is QUADRILATERAL
        if ( cell_it->cell_type ( ).tag ( ) == CellType::QUADRILATERAL )
        {
            num_quad++;
            // look if cell is in first quadrant
            // (for at least one vertex of cell is x-coord > 0 and y-coord > 0)
            bool first_quadrant = false;
            // loop over all vertices of cell
            for ( IncidentEntityIterator vertex_it = cell_it->begin_incident ( 0 );
                  vertex_it != cell_it->end_incident ( 0 );
                  ++vertex_it )
            {
                std::vector<Coordinate> coord_vertex;
                vertex_it->get_coordinates ( coord_vertex, 0 );
                if ( coord_vertex[0] > 0. && coord_vertex[1] > 0. )
                {
                    first_quadrant = true;
                }
            }
            // depending if in first quadrant or not define refinement type
            if ( first_quadrant )
            {
                cell_it->set ( "refinement_type", 3 );
            }
            else if ( !first_quadrant )
            {
                cell_it->set ( "refinement_type", 4 );
            }
        }
    }

    // step 3b) refine mesh depending on refinement type
    mesh = mesh->refine ( "refinement_type" );

    // extract boundary mesh
    MeshPtr boundary_mesh = mesh->extract_boundary_mesh ( );

    // output of refined mesh and the corresponding boundary mesh
    ScopedPtr<Writer>::Type writer ( new VtkWriter ( ) );
    writer->write ( "L_domain2_refined_by_attribute.vtu", *mesh );
    writer->write ( "L_domain2_boundary_refined_by_attribute.vtu", *boundary_mesh );

    return 0;
}

// function to find the local edge number of an edge in a cell

int find_edge_in_cell ( const Entity& edge, const Entity& cell )
{
    // sorted array of the ids of the vertices of the input edge
    SortedArray<Id> edge_verts =
            SortedArray<Id>( std::vector<Id>( edge.begin_vertex_ids ( ),
            edge.end_vertex_ids ( ) ) );
    int k = 0;
    // loop over all egdes of the cell
    for ( IncidentEntityIterator edge_it = cell.begin_incident ( tdim - 1 );
          edge_it != cell.end_incident ( tdim - 1 );
          ++edge_it )
    {
        bool found = true;
        // loop over all ids of the vertices
        for ( VertexIdIterator v_it = edge_it->begin_vertex_ids ( );
              v_it != edge_it->end_vertex_ids ( ); ++v_it )
        {
            // check if the vertex belongs to the input egde
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
    // if the edge is not egde of cell return false
    return false;
}
