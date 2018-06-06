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

/// \author Jonathan Schwegler

/// \brief This program demonstrates the usage of
/// adapt_boundary_to_function(). For this it contains three implementations
/// of the class BoundaryDomainDescriptor.
/// \details Note that every BoundaryDomainDescriptor does not work with
/// every mesh. For example EllipsoidWithHole will always need
/// a starting mesh with a hole and with certain MaterialNumbers
/// on the boundaries (see implementation).
/// It's always good if the boundary of the initial mesh does lie
/// on the described boundary but not necessary.
#include <iostream>
#include <fstream>

#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const int DEBUG_LEVEL = 1;
const int GDIM = 2;

static const char* datadir = MESHES_DATADIR;

class Ellipsoid : public BoundaryDomainDescriptor
{
  public:

    //Constructor for a circle or a sphere

    Ellipsoid ( const Coordinate radius ) : a_ ( radius ), b_ ( radius ), c_ ( radius )
    {
    };

    //Constructor for a 2D ellipse

    Ellipsoid ( const Coordinate a, const Coordinate b ) : a_ ( a ), b_ ( b ), c_ ( 0. )
    {
    };

    //Constructor for a 3D ellipsoid

    Ellipsoid ( const Coordinate a, const Coordinate b, const Coordinate c ) :
    a_ ( a ), b_ ( b ), c_ ( c )
    {
    };

    //Implementation of the ellipsoid formula:
    // ((x-0.5)/a)^2 + ((y-0.5)/b)^2 + ((z-0.5)/c)^2 = 1.

    Coordinate eval_func ( const std::vector<Coordinate> &x, MaterialNumber mat_num ) const
    {
        Coordinate third_comp = 0.;
        if ( GDIM == 3 ) third_comp = ( x[2] - 0.5 )*( x[2] - 0.5 ) / ( c_ * c_ );
        return 1. - ( x[0] - 0.5 )*( x[0] - 0.5 ) / ( a_ * a_ )-( x[1] - 0.5 )*( x[1] - 0.5 ) / ( b_ * b_ ) - third_comp;
    }

    // Computation of the gradiant of the ellipsoid formula

    std::vector<Coordinate> eval_grad ( const std::vector<Coordinate> &x, MaterialNumber mat_num ) const
    {
        std::vector<Coordinate> grad ( GDIM );
        grad[0] = -2. * ( x[0] - 0.5 ) / ( a_ * a_ );
        grad[1] = -2. * ( x[1] - 0.5 ) / ( b_ * b_ );
        if ( GDIM == 3 )
            grad[2] = -2. * ( x[2] - 0.5 ) / ( c_ * c_ );
        return grad;
    }

    // Parameter describing the ellipsoid
    const Coordinate a_;
    const Coordinate b_;
    const Coordinate c_;
};

// Class for creating ellipsoids with holes.

class EllipsoidWithHole : public BoundaryDomainDescriptor
{
  public:

    //Constructor for a circle (or sphere) with a circle (or sphere) hole

    EllipsoidWithHole ( const Coordinate outer_radius,
                        const Coordinate inner_radius ) :
    outer_a_ ( outer_radius ), inner_a_ ( inner_radius ), outer_b_ ( outer_radius ),
    inner_b_ ( inner_radius ), outer_c_ ( outer_radius ), inner_c_ ( inner_radius )
    {
    };

    //Constructor for a 2D ellipse with a ellipse hole

    EllipsoidWithHole ( const Coordinate outer_a, const Coordinate outer_b,
                        const Coordinate inner_a, const Coordinate inner_b ) :
    outer_a_ ( outer_a ), inner_a_ ( inner_a ), outer_b_ ( outer_b ), inner_b_ ( inner_b ),
    outer_c_ ( 0. ), inner_c_ ( 0. )
    {
    };

    //Constructor for a 3D ellipsoid with a ellipsoid hole

    EllipsoidWithHole ( const Coordinate outer_a, const Coordinate outer_b,
                        const Coordinate outer_c, const Coordinate inner_a,
                        const Coordinate inner_b, const Coordinate inner_c ) :
    outer_a_ ( outer_a ), inner_a_ ( inner_a ), outer_b_ ( outer_b ), inner_b_ ( inner_b ),
    outer_c_ ( outer_c ), inner_c_ ( inner_c )
    {
    };

    //Implementation of the ellipsoid formulas.
    // Everything with MaterialNumber 11 will be mapped to the inner ellipsoid
    // Everything with MaterialNumber 12 to the outer one.

    Coordinate eval_func ( const std::vector<Coordinate> &x, MaterialNumber mat_num ) const
    {
        Coordinate a;
        Coordinate b;
        Coordinate c;
        if ( mat_num == 11 )
        {
            a = inner_a_;
            b = inner_b_;
            c = inner_c_;
        }
        else if ( mat_num == 12 )
        {
            a = outer_a_;
            b = outer_b_;
            c = outer_c_;
        }
        else
        {
            return 0.;
        }
        Coordinate third_comp = 0.;
        if ( GDIM == 3 ) third_comp = ( x[2] - 0.5 )*( x[2] - 0.5 ) / ( c * c );
        return 1. - ( x[0] - 0.5 )*( x[0] - 0.5 ) / ( a * a )-( x[1] - 0.5 )*( x[1] - 0.5 ) / ( b * b ) - third_comp;
    }

    // Computation of the corresponding gradients.

    std::vector<Coordinate> eval_grad ( const std::vector<Coordinate> &x, MaterialNumber mat_num ) const
    {
        Coordinate a;
        Coordinate b;
        Coordinate c;
        if ( mat_num == 11 )
        {
            a = inner_a_;
            b = inner_b_;
            c = inner_c_;
        }
        else if ( mat_num == 12 )
        {
            a = outer_a_;
            b = outer_b_;
            c = outer_c_;
        }
        else
        {
            return std::vector<Coordinate>( GDIM, 0. );
        }
        std::vector<Coordinate> grad ( GDIM );
        grad[0] = -2. * ( x[0] - 0.5 ) / ( a * a );
        grad[1] = -2. * ( x[1] - 0.5 ) / ( b * b );
        if ( GDIM == 3 )
            grad[2] = -2. * ( x[2] - 0.5 ) / ( c * c );
        return grad;
    }

    //Parameters to describe the to ellipsoids
    const Coordinate outer_a_;
    const Coordinate outer_b_;
    const Coordinate outer_c_;

    const Coordinate inner_a_;
    const Coordinate inner_b_;
    const Coordinate inner_c_;
};

// Descriptor for a heart <3 (in 2D)
// This class may behave weird for some starting meshs due to
// discontinuity of the zero set.
// Works fine with unit_square.inp

class Heart : public BoundaryDomainDescriptor
{
  public:

    Heart ( )
    {
    };

    // Function taken from http://www.wolframalpha.com/input/?i=first%20heart%20curve
    // (slightly modified)

    Coordinate eval_func ( const std::vector< Coordinate >& p, MaterialNumber mat_num ) const
    {
        assert ( GDIM == 2 );
        const Coordinate x = p[0] - 0.5;
        const Coordinate y = p[1] - 0.5;
        return pow ( ( x * x + y * y - 0.15 ), 3. ) - x * x * y * y*y;
    }

    std::vector<Coordinate> eval_grad ( const std::vector< Coordinate >& p, MaterialNumber mat_num ) const
    {
        const Coordinate x = p[0] - 0.5;
        const Coordinate y = p[1] - 0.5;
        std::vector<Coordinate> grad ( GDIM, 0. );
        grad[0] = 2. * x * 3. * pow ( ( x * x + y * y - 0.15 ), 2. ) - 2. * x * y * y*y;
        grad[1] = 2. * y * 3. * pow ( ( x * x + y * y - 0.15 ), 2. ) - 3. * x * x * y*y;
        return grad;
    }
};

int main ( int argc, char** argv )
{

    int refinement_level = 4;

    std::cout << "Using refinement level: " << refinement_level << std::endl;
    const TDim tdim = GDIM;
    const GDim gdim = GDIM;

    std::ofstream debug_file ( "bdd_demo_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_file );

    //Example for a circle
    std::string filename_circle = std::string ( datadir ) + std::string ( "unit_square_inner_square.inp" );
    MeshPtr circle_mesh = read_mesh_from_file ( filename_circle, tdim, gdim, 0 );

    Coordinate radius = 1.;
    Ellipsoid circle ( radius, radius );

    //initial fitting (if the boundary of the starting mesh does not coincide
    // with the coundary described by the Ellipsoid
    adapt_boundary_to_function ( circle_mesh, circle );
    for ( int r = 0; r < refinement_level; ++r )
    {
        circle_mesh = circle_mesh->refine ( );
        adapt_boundary_to_function ( circle_mesh, circle );
    }

    VtkWriter writer;
    std::string output_file_circle = std::string ( "circle_mesh.vtu" );
    writer.write ( output_file_circle.c_str ( ), *circle_mesh );

    //Example of a circle with a hole
    std::string filename_hole = std::string ( datadir ) + "unit_square_inner_hole.inp";
    MeshPtr hole_mesh = read_mesh_from_file ( filename_hole, tdim, gdim, 0 );

    Coordinate outer_radius = 1.;
    Coordinate inner_radius = 0.1;
    EllipsoidWithHole circle_with_hole ( outer_radius, inner_radius );

    adapt_boundary_to_function ( hole_mesh, circle_with_hole );
    for ( int r = 0; r < refinement_level; ++r )
    {
        hole_mesh = hole_mesh->refine ( );
        adapt_boundary_to_function ( hole_mesh, circle_with_hole );
    }

    std::string output_file_hole = std::string ( "circle_with_hole_mesh.vtu" );
    writer.write ( output_file_hole.c_str ( ), *hole_mesh );

    //Example of a heart shaped mesh
    std::string filename_heart = std::string ( datadir ) + "unit_square.inp";
    MeshPtr heart_mesh = read_mesh_from_file ( filename_heart, tdim, gdim, 0 );

    Heart heart;

    adapt_boundary_to_function ( heart_mesh, heart );
    for ( int r = 0; r < refinement_level; ++r )
    {
        heart_mesh = heart_mesh->refine ( );
        adapt_boundary_to_function ( heart_mesh, heart );
    }

    std::string output_file_heart = std::string ( "heart_mesh.vtu" );
    writer.write ( output_file_heart.c_str ( ), *heart_mesh );

    LogKeeper::get_log ( "debug" ).flush ( );
    return 0;
}
