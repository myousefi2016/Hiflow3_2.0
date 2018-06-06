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

#include "hiflow.h"
#include <vector>
#include <iostream>
#include <fstream>

#include "test.h"

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

const int DIMENSION = 3;
const int MASTER_RANK = 0;

const int DEBUG_LEVEL = 2;

std::vector< std::vector<double> > TEST_POINTS ( 0 );
std::vector< std::vector<double> > PROJ_POINTS ( 0 );

// get test points for test with test_number. return number of test points

int get_test_points ( int test_number )
{
    if ( test_number == 0 )
    {
        int num_pts = 5;
        TEST_POINTS.resize ( num_pts, std::vector<double>( DIMENSION ) );

        PROJ_POINTS.resize ( num_pts, std::vector<double>( DIMENSION ) );

        TEST_POINTS[0][0] = 0.4899;
        TEST_POINTS[0][1] = 0.19;
        TEST_POINTS[0][2] = 0.5;

        PROJ_POINTS[0][0] = 0.45;
        PROJ_POINTS[0][1] = 0.19;
        PROJ_POINTS[0][2] = 0.41;

        TEST_POINTS[1][0] = 0.51;
        TEST_POINTS[1][1] = 0.1899;
        TEST_POINTS[1][2] = 0.2;

        PROJ_POINTS[1][0] = 0.51;
        PROJ_POINTS[1][1] = 0.15;
        PROJ_POINTS[1][2] = 0.2;

        TEST_POINTS[2][0] = 0.61;
        TEST_POINTS[2][1] = 0.3;
        TEST_POINTS[2][2] = 0.2;

        PROJ_POINTS[2][0] = 0.55;
        PROJ_POINTS[2][1] = 0.25;
        PROJ_POINTS[2][2] = 0.2;

        TEST_POINTS[3][0] = 99;
        TEST_POINTS[3][1] = 100;
        TEST_POINTS[3][2] = 101;

        PROJ_POINTS[3][0] = 2.5;
        PROJ_POINTS[3][1] = 0.41;
        PROJ_POINTS[3][2] = 0.41;

        TEST_POINTS[4][0] = 2;
        TEST_POINTS[4][1] = 0.3;
        TEST_POINTS[4][2] = 0.41;

        PROJ_POINTS[4][0] = 2;
        PROJ_POINTS[4][1] = 0.3;
        PROJ_POINTS[4][2] = 0.41;

        return num_pts;
    }
    else if ( test_number == 1 )
    {

        int num_pts = 6;
        TEST_POINTS.resize ( num_pts, std::vector<double>( DIMENSION ) );

        PROJ_POINTS.resize ( num_pts, std::vector<double>( DIMENSION ) );

        TEST_POINTS[0][0] = 0.52;
        TEST_POINTS[0][1] = 0.22;
        TEST_POINTS[0][2] = 0.1;

        PROJ_POINTS[0][0] = 0.537647794495919;
        PROJ_POINTS[0][1] = 0.231752237846649;
        PROJ_POINTS[0][2] = 0.100208280810009;

        TEST_POINTS[1][0] = 0.53;
        TEST_POINTS[1][1] = 0.1999;
        TEST_POINTS[1][2] = 0.6;

        PROJ_POINTS[1][0] = 0.548635152327107;
        PROJ_POINTS[1][1] = 0.194906524582428;
        PROJ_POINTS[1][2] = 0.41;

        TEST_POINTS[2][0] = 0.6;
        TEST_POINTS[2][1] = 0.3;
        TEST_POINTS[2][2] = 0.2;

        PROJ_POINTS[2][0] = 0.534117290657391;
        PROJ_POINTS[2][1] = 0.236548039516744;
        PROJ_POINTS[2][2] = 0.199702020849231;

        TEST_POINTS[3][0] = -0.1;
        TEST_POINTS[3][1] = 0.51;
        TEST_POINTS[3][2] = 0.51;

        PROJ_POINTS[3][0] = 0.0;
        PROJ_POINTS[3][1] = 0.41;
        PROJ_POINTS[3][2] = 0.41;

        TEST_POINTS[4][0] = 99;
        TEST_POINTS[4][1] = 100;
        TEST_POINTS[4][2] = 101;

        PROJ_POINTS[4][0] = 2.5;
        PROJ_POINTS[4][1] = 0.41;
        PROJ_POINTS[4][2] = 0.41;

        TEST_POINTS[5][0] = 0.54;
        TEST_POINTS[5][1] = 0.24;
        TEST_POINTS[5][2] = 0.2;

        PROJ_POINTS[5][0] = 0.53462356568747;
        PROJ_POINTS[5][1] = 0.235883167936807;
        PROJ_POINTS[5][2] = 0.199934602880438;
        return num_pts;
    }
    else
    {
        LOG_DEBUG ( 1, "ERROR: no test points for test number " << test_number );
        return 0;
    }
}

MeshPtr read_and_distribute_mesh ( std::string mesh_filename )
{
    LOG_DEBUG ( 2, "Loading mesh from file " << mesh_filename );
    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MeshPtr master_mesh;

    if ( rank == MASTER_RANK )
    {
        master_mesh = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0 );
        //default material number, so the geom search works correctly...
        set_default_material_number_on_bdy ( master_mesh, 0 );
    }
    MeshPtr local_mesh = partition_and_distribute ( master_mesh, MASTER_RANK, MPI_COMM_WORLD );
    SharedVertexTable shared_verts;
    return compute_ghost_cells ( *local_mesh, MPI_COMM_WORLD, shared_verts );
}

// for a given test and mesh compute the projected points

void project_points ( int test, MeshPtr mesh )
{
    GeometricSearchPtr grid_search = mesh->get_geometric_search ( );

    int num_test_points = get_test_points ( test );
    std::vector<double> proj_pt ( DIMENSION );
    for ( int j = 0; j < num_test_points; ++j )
    {
        int id; //local, unused variable
        proj_pt = grid_search->find_closest_point_parallel ( TEST_POINTS[j], id, MPI_COMM_WORLD );
        for ( int d = 0; d < DIMENSION; ++d )
        {
            LOG_DEBUG ( 3, "Diff expected and computed for point " << j << " and component " << d << ":" << proj_pt[d] - PROJ_POINTS[j][d] );
            TEST_EQUAL_EPS ( proj_pt[d], PROJ_POINTS[j][d], 10e-6 );
        }
    }
    LOG_DEBUG ( 1, "Successfully tested " << num_test_points << " points in test " << test );
}

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );
    //set LOG_DEBUG
    std::ofstream debug_file ( "closest_point__test.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );

    // first test
    std::stringstream mesh_name;
    mesh_name << MESH_DATADIR << "dfg_bench3d_rect.inp";
    MeshPtr mesh = read_and_distribute_mesh ( mesh_name.str ( ) );
    project_points ( 0, mesh );

    // second test
    mesh_name.str ( std::string ( ) ); //clear stringstream
    mesh_name << MESH_DATADIR << "dfg_bench3d_cyl.vtu";
    mesh = read_and_distribute_mesh ( mesh_name.str ( ) );
    project_points ( 1, mesh );

    LogKeeper::get_log ( "debug" ).flush ( );
    MPI_Finalize ( );
    return 0;
}
