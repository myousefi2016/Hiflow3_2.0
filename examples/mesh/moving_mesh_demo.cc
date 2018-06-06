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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "hiflow.h"

using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;
const int max_steps = 5;
static const char* datadir = MESHES_DATADIR;

void move_points ( MeshPtr mesh, const std::vector<double>& displacement );

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );
    MPI_Comm comm = MPI_COMM_WORLD;
    std::string filename = std::string ( datadir ) + std::string ( "unit_cube.inp" );
    MeshPtr mesh = read_mesh_from_file ( filename, tdim, gdim, &comm );

    for ( EntityIterator it = mesh->begin ( 0 ); it != mesh->end ( 0 ); ++it )
    {
        std::cout << *it << "\n";
    }

    std::vector<double> displacement_of_points ( mesh->gdim ( ) * mesh->num_entities ( 0 ) );

    for ( int t = 0; t < max_steps; ++t )
    {
        if ( t % 2 == 0 )
        {
            for ( std::vector<double>::iterator it = displacement_of_points.begin ( );
                  it != displacement_of_points.end ( ); it += mesh->gdim ( ) )
            {
                if ( it == displacement_of_points.begin ( ) + ( 3 * mesh->gdim ( ) ) )
                {
                    ( *it ) = 1.0;
                    ( *( it + 1 ) ) = 0.0;
                    ( *( it + 2 ) ) = 0.0;
                }
                else
                {
                    ( *it ) = 0.0;
                    ( *( it + 1 ) ) = 0.0;
                    ( *( it + 2 ) ) = 0.0;
                }
            }
        }
        else
        {
            for ( std::vector<double>::iterator it = displacement_of_points.begin ( );
                  it != displacement_of_points.end ( ); it += mesh->gdim ( ) )
            {
                ( *it ) = 0.0;
                ( *( it + 1 ) ) = 0.0;
                ( *( it + 2 ) ) = 0.0;
            }
        }
        move_points ( mesh, displacement_of_points );

        std::ostringstream time_str;
        time_str << t;
        std::string filename_out = std::string ( "moving_mesh_" ) + time_str.str ( ) + ".vtu";
        ScopedPtr<Writer>::Type writer ( new VtkWriter ( ) );
        writer->write ( filename_out.c_str ( ), *mesh );

        for ( EntityIterator it = mesh->begin ( 0 ); it != mesh->end ( 0 ); ++it )
        {
            std::cout << *it << "\n";
        }
    }

    MPI_Finalize ( );
}

void move_points ( MeshPtr mesh,
                   const std::vector<double>& displacement )
{
    assert ( static_cast < int > ( displacement.size ( ) ) == mesh->gdim ( ) * mesh->num_entities ( 0 ) );
    int gdim = mesh->gdim ( );
    // for now
    assert ( gdim == 3 );

    int num_points = mesh->num_entities ( 0 );
    std::vector<double> delta_x, delta_y, delta_z;
    delta_x.reserve ( num_points );
    delta_y.reserve ( num_points );
    delta_z.reserve ( num_points );

    for ( std::vector<double>::const_iterator it = displacement.begin ( );
          it != displacement.end ( ); it += gdim )
    {
        delta_x.push_back ( *it );
        delta_y.push_back ( *( it + 1 ) );
        delta_z.push_back ( *( it + 2 ) );
    }

    // TODO(Thomas): Move this functionality to mesh
    AttributePtr displacement_vertices_x = AttributePtr ( new DoubleAttribute ( delta_x ) );
    AttributePtr displacement_vertices_y = AttributePtr ( new DoubleAttribute ( delta_y ) );
    AttributePtr displacement_vertices_z = AttributePtr ( new DoubleAttribute ( delta_z ) );
    mesh->add_attribute ( std::string ( "__displacement_vertices_x__" ), 0, displacement_vertices_x );
    mesh->add_attribute ( std::string ( "__displacement_vertices_y__" ), 0, displacement_vertices_y );
    mesh->add_attribute ( std::string ( "__displacement_vertices_z__" ), 0, displacement_vertices_z );
}
