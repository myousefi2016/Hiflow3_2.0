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

/// \brief This program demonstrates the possibility of defining
/// custom geometry during refinement. It transforms a unit square
/// into its circumscribing circle during refinement.

///  \author Staffan Ronnas

#include "hiflow.h"

#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 2;
const GDim gdim = 2;
static const char* datadir = MESHES_DATADIR;
const int REF_LEVEL = 5;

class RefineTowardsCircle
{
  public:

    RefineTowardsCircle ( double x, double y, double radius )
    : x_ ( x ), y_ ( y ), r_ ( radius )
    {
    }

    void add_boundary_vertex ( int v_id )
    {
        boundary_verts_.insert ( v_id );
    }

    void operator() ( const Entity& cell, std::vector<Coordinate>& coords )
    {
        // get necessary data
        std::vector<int> cell_vertex_ids ( cell.begin_vertex_ids ( ), cell.end_vertex_ids ( ) );
        std::vector<Coordinate> cell_coords;
        cell.get_coordinates ( cell_coords );
        const CellType& cell_type = cell.cell_type ( );
        const int num_refined_vertices = cell_type.num_vertices ( );
        coords.resize ( num_refined_vertices * gdim, 0.0 );

        // Compute refined vertices as barycenters of parent cell vertices. (Standard geometry).
        for ( int v = 0; v < num_refined_vertices; ++v )
        {
            bool should_be_mapped = true; // check if we should map the vertex onto the circle.

            const std::vector<int>& refined_vertex = cell_type.refined_vertex ( v );
            const int num_super_vertices = refined_vertex.size ( );

            // Only map mid-edge vertices.
            if ( refined_vertex.size ( ) != 2 )
            {
                should_be_mapped = false;
            }

            // compute v:th refined vertex as barycenter of super vertices
            for ( int i = 0; i < num_super_vertices; ++i )
            {
                const int super_v = refined_vertex[i];
                for ( int c = 0; c < gdim; ++c )
                {
                    coords[v * gdim + c] += cell_coords[gdim * super_v + c];
                }

                // only map the vertex if all its parent vertices are on the boundary
                should_be_mapped = should_be_mapped &&
                        boundary_verts_.find ( cell_vertex_ids.at ( super_v ) ) != boundary_verts_.end ( );
            }

            // normalize barycenter coordinates
            for ( int c = 0; c < gdim; ++c )
            {
                coords[v * gdim + c] /= static_cast < Coordinate > ( num_super_vertices );
            }

            // map the vertex onto the circle.
            if ( should_be_mapped )
            {
                map_boundary_vertex ( refined_vertex, &coords[v * gdim] );
            }
        }

    }

    void map_boundary_vertex ( const std::vector<int>& refined_vertex, double* vertex ) const
    {
        // Map input vertex in direction normal to face. The new
        // position is computed using quadratic formula to solve the
        // line-circle intersection problem.
        Vec<2, double> n = get_normal ( refined_vertex ); // get normal of face defined by refined_vertex
        Vec<2, double> center;
        center[0] = x_;
        center[1] = y_;
        Vec<2, double> pt ( vertex );
        Vec<2, double> delta = pt - center;
        if ( n[0] == 0. && n[1] == 0. )
        {
            return;
        }
        const double a = dot ( n, n );
        const double b = 2. * dot ( delta, n );
        const double c = dot ( delta, delta ) - r_*r_;
        const double t = -b / ( 2. * a ) + std::sqrt ( b * b - 4 * a * c ) / ( 2. * a );

        pt += t*n;

        vertex[0] = pt[0];
        vertex[1] = pt[1];
    }

    Vec<2, double> get_normal ( const std::vector<int>& refined_vertex ) const
    {
        Vec<2, double> n;
        n[0] = 0.;
        n[1] = 0.;

        assert ( refined_vertex[1] == ( ( refined_vertex[0] + 1 ) % 4 ) );

        switch ( refined_vertex[0] )
        {
            case 0:
                n[0] = 0.;
                n[1] = -1.; // bottom face
                return n;
            case 1:
                n[0] = 1.;
                n[1] = 0.; // right face
                return n;
            case 2:
                n[0] = 0.;
                n[1] = 1.; // top face
                return n;
            case 3:
                n[0] = -1.;
                n[1] = 0.; // left face
                return n;
            default:
                assert ( false );
                return n;
        }
    }

  private:
    double x_, y_, r_;
    std::set<int> boundary_verts_;
};

void update_bdy ( RefineTowardsCircle& ref_circle, MeshPtr mesh )
{
    // Update ref_circle with all vertices that lie on the boundary
    MeshPtr bdy_mesh = mesh->extract_boundary_mesh ( );

    for ( EntityIterator it = bdy_mesh->begin ( 0 ), end_it = bdy_mesh->end ( 0 ); it != end_it; ++it )
    {
        ref_circle.add_boundary_vertex ( it->id ( ) );
    }
}

void write_mesh ( MeshPtr mesh, int level )
{
    VtkWriter writer;
    std::ostringstream sstr;
    sstr << "unitsquare_circle.0" << level << ".vtu";
    writer.write ( sstr.str ( ).c_str ( ), *mesh );
}

int main ( int argc, char** argv )
{

    MPI_Init ( &argc, &argv );

    std::string filename = std::string ( datadir ) + std::string ( "unit_square.inp" );
    MeshPtr mesh = read_mesh_from_file ( filename, tdim, gdim, 0 );

    RefineTowardsCircle ref_circle ( 0.5, 0.5, 0.5 * std::sqrt ( 2. ) );

    write_mesh ( mesh, 0 );

    for ( int r = 0; r < REF_LEVEL; ++r )
    {
        update_bdy ( ref_circle, mesh );
        mesh->set_refined_geometry_function ( ref_circle );
        mesh = mesh->refine ( );

        write_mesh ( mesh, r + 1 );
    }

    MPI_Finalize ( );
    return 0;
}
