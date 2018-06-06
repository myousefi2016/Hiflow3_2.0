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

/// \brief
/// This program demonstrates the use of boost::filter_iterator to
/// create a customized iterator, which loops over the vertices that
/// lie in the first octant. This is accomplished by creating a filter,
/// VerticesInFirstOctantFilter, which returns true for those vertices
/// that should be visited, and false for the others.

#include <iostream>
#include <iomanip>

#include <boost/iterator/filter_iterator.hpp>

#include "hiflow.h"

using namespace hiflow::mesh;

std::vector<Coordinate> compute_barycenter ( const Entity& entity );

static const char* datadir = MESHES_DATADIR;

const int OUTPUT_COUNT = 10;

const TDim tdim = 3;
const GDim gdim = 3;

struct VerticesInFirstOctantFilter
{

    bool operator() ( const Entity& entity ) const
    {
        std::vector<Coordinate> coords;
        entity.get_coordinates ( coords );
        for ( int i = 0; i < entity.gdim ( ); ++i )
        {
            if ( coords[i] < 0.0 )
            {
                return false;
            }
        }
        return true;
    }
};

int main ( int argc, char *argv[] )
{

    std::cout << std::setprecision ( 3 );

    // create Mesh Database
    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    // create reader
    Reader* reader = new UcdReader ( mb );

    // read mesh
    std::string filename = std::string ( datadir ) + std::string ( "dfg_bench3d_cyl.inp" );
    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    // output number of entities of each dimension
    for ( int d = 0; d <= tdim; ++d )
    {
        std::cout << "The mesh has " << mesh->num_entities ( d )
                << " entities of dimension " << d << "\n";
    }

    typedef boost::filter_iterator<VerticesInFirstOctantFilter, EntityIterator> FirstOctantIterator;

    FirstOctantIterator it ( VerticesInFirstOctantFilter ( ), mesh->begin ( 0 ), mesh->end ( 0 ) );
    FirstOctantIterator end ( VerticesInFirstOctantFilter ( ), mesh->end ( 0 ), mesh->end ( 0 ) );

    std::cout << "\n\nFiltered iterator example\n";
    for (; it != end; ++it )
    {
        std::cout << "vertex " << it->id ( ) << " lies in first octant\n";
    }

    delete reader;
    delete mb;

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
