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
#include <iomanip>

#include "hiflow.h"

/// \author Staffan Ronnas

/// \brief This program demonstrates simple iteration over entities in
/// a mesh. Both global iteration, using the Mesh::begin() and
/// Mesh::end() functions, as well as incident iteration, using
/// Mesh::begin_incident() / Mesh::end_incident(), are demonstrated.

using namespace hiflow;
using namespace hiflow::mesh;

std::vector<Coordinate> compute_barycenter ( const Entity& entity );

static const char* datadir = MESHES_DATADIR;

const int OUTPUT_COUNT = 10;

const TDim tdim = 3;
const GDim gdim = 3;

int main ( int argc, char *argv[] )
{

    std::cout << std::setprecision ( 3 );

    // create MeshBuilder
    MeshDbViewBuilder mb ( tdim, gdim );

    // create reader
    UcdReader reader ( &mb );

    // read mesh
    std::string filename = std::string ( datadir ) + std::string ( "dfg_bench3d_cyl.inp" );
    MeshPtr mesh;
    reader.read ( filename.c_str ( ), mesh );

    // output number of entities of each dimension
    for ( int d = 0; d <= tdim; ++d )
    {
        std::cout << "The mesh has " << mesh->num_entities ( d )
                << " entities of dimension " << d << "\n";
    }

    // iteration over cells
    std::cout << "\n\n === Iteration over cells === \n\n";
    for ( EntityIterator it = mesh->begin ( tdim ); it != mesh->end ( tdim ); ++it )
    {

        // compute barycenter of cell
        const std::vector<Coordinate> barycenter = compute_barycenter ( *it );

        if ( it->id ( ) < OUTPUT_COUNT )
        {
            std::cout << "Cell " << it->id ( ) << " has barycenter "
                    << string_from_range ( barycenter.begin ( ), barycenter.end ( ) ) << "\n";
        }
    }

    // iteration over faces of cells
    std::cout << "\n\n === Iteration over faces of cells === \n\n";
    for ( EntityIterator it = mesh->begin ( tdim ); it != mesh->end ( tdim ); ++it )
    {
        if ( it->id ( ) < OUTPUT_COUNT )
        {
            std::cout << "Cell " << it->id ( ) << " has faces: ";
            for ( IncidentEntityIterator iit = mesh->begin_incident ( *it, tdim - 1 );
                  iit != mesh->end_incident ( *it, tdim - 1 ); ++iit )
            {
                std::cout << iit->id ( ) << " ";
            }
            std::cout << "\n";
        }
    }

    // iteration over vertices of cells
    std::cout << "\n\n === Iteration over vertices of cells === \n\n";
    for ( EntityIterator it = mesh->begin ( tdim ); it != mesh->end ( tdim ); ++it )
    {
        if ( it->id ( ) < OUTPUT_COUNT )
        {
            std::cout << "Cell " << it->id ( ) << " has vertices: \n(\t";
            for ( IncidentEntityIterator iit = mesh->begin_incident ( *it, 0 );
                  iit != mesh->end_incident ( *it, 0 ); ++iit )
            {
                std::vector<Coordinate> coords;
                iit->get_coordinates ( coords );
                std::cout << iit->id ( ) << ": ";
                std::copy ( coords.begin ( ), coords.end ( ), std::ostream_iterator<Coordinate>( std::cout, " " ) );
                std::cout << "\t";
            }
            std::cout << ")\n\n";
        }
    }

    return 0;
}

std::vector<Coordinate> compute_barycenter ( const Entity& entity )
{
    std::vector<Coordinate> barycenter ( gdim, 0.0 );
    const int num_vertices = entity.num_vertices ( );

    // add coordinates
    for ( int v = 0; v < num_vertices; ++v )
    {
        std::vector<Coordinate> pt;
        entity.get_coordinates ( pt, v );
        for ( int c = 0; c < gdim; ++c )
        {
            barycenter[c] += pt[c];
        }
    }

    // normalize
    for ( int c = 0; c < gdim; ++c )
    {
        barycenter[c] /= static_cast < Coordinate > ( num_vertices );
    }

    return barycenter;
}
