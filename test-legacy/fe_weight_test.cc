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

/// \author Michael Schick

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include "common/macros.h"

#include "../src/fem/felagrange_hex.h"
#include "../src/fem/felagrange_tet.h"
#include "../src/fem/felagrange_tri.h"
#include "../src/fem/felagrange_quad.h"
#include "../src/fem/felagrange_line.h"
#include "test.h"

using namespace hiflow::doffem;

int main ( int argc, char** argv )
{
    // Test if sum of weights is equal to one
    // and if sum of derivatives is equal to zero

    int nb_3d_pts_line = 11;
    int nb_2d_pts_line = 11;
    int nb_1d_pts_line = 11;

    double h_3d = 1.0 / ( nb_3d_pts_line - 1. );
    double h_2d = 1.0 / ( nb_2d_pts_line - 1. );
    //double h_1d = 1.0/(nb_1d_pts_line - 1.);

    for ( int deg = 0; deg < 8; ++deg )
    {
        std::cout << "Degree: " << deg << std::endl;

        FELagrangeHex<double> fe_hex;
        FELagrangeTet<double> fe_tet;
        FELagrangeTri<double> fe_tri;
        FELagrangeQuad<double> fe_quad;
        FELagrangeLine<double> fe_line;

        fe_hex.set_fe_deg ( deg );
        fe_tet.set_fe_deg ( deg );
        fe_tri.set_fe_deg ( deg );
        fe_quad.set_fe_deg ( deg );
        fe_line.set_fe_deg ( deg );

        fe_hex.init ( );
        fe_tet.init ( );
        fe_tri.init ( );
        fe_quad.init ( );
        fe_line.init ( );

        // Test HEXAHEDRON

        for ( int k = 0; k < nb_3d_pts_line; ++k )
            for ( int j = 0; j < nb_3d_pts_line; ++j )
                for ( int i = 0; i < nb_3d_pts_line; ++i )
                {
                    std::vector<double> coord ( 3 );

                    coord[0] = h_3d*i;
                    coord[1] = h_3d*j;
                    coord[2] = h_3d*k;

                    assert ( coord[0] >= 0.0 && coord[0] <= 1.0 );
                    assert ( coord[1] >= 0.0 && coord[1] <= 1.0 );
                    assert ( coord[2] >= 0.0 && coord[2] <= 1.0 );

                    std::vector<double> weight ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_x ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_y ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_z ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_xx ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_xy ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_xz ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_yy ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_yz ( fe_hex.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_zz ( fe_hex.get_nb_dof_on_cell ( ) );

                    fe_hex.N ( coord, weight );
                    fe_hex.N_x ( coord, weight_x );
                    fe_hex.N_y ( coord, weight_y );
                    fe_hex.N_z ( coord, weight_z );
                    fe_hex.N_xx ( coord, weight_xx );
                    fe_hex.N_xy ( coord, weight_xy );
                    fe_hex.N_xz ( coord, weight_xz );
                    fe_hex.N_yy ( coord, weight_yy );
                    fe_hex.N_yz ( coord, weight_yz );
                    fe_hex.N_zz ( coord, weight_zz );

                    // Check
                    double sum = 0.0;
                    double sum_x = 0.0;
                    double sum_y = 0.0;
                    double sum_z = 0.0;
                    double sum_xx = 0.0;
                    double sum_xy = 0.0;
                    double sum_xz = 0.0;
                    double sum_yy = 0.0;
                    double sum_yz = 0.0;
                    double sum_zz = 0.0;

                    for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                    {
                        sum += weight[w];
                        sum_x += weight_x[w];
                        sum_y += weight_y[w];
                        sum_z += weight_z[w];
                        sum_xx += weight_xx[w];
                        sum_xy += weight_xy[w];
                        sum_xz += weight_xz[w];
                        sum_yy += weight_yy[w];
                        sum_yz += weight_yz[w];
                        sum_zz += weight_zz[w];
                    }

                    TEST_EQUAL_EPS ( sum, 1.0, 1.e-10 );
                    TEST_EQUAL_EPS ( sum_x, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_y, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_z, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_xx, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_xy, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_xz, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_yy, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_yz, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_zz, 0.0, 1.e-9 );
                }

        // Test TETRAHEDRON

        for ( int k = 0; k < nb_3d_pts_line; ++k )
            for ( int j = 0; j < nb_3d_pts_line - k; ++j )
                for ( int i = 0; i < nb_3d_pts_line - k - j; ++i )
                {
                    std::vector<double> coord ( 3 );

                    coord[2] = h_3d*k;
                    coord[1] = ( 1.0 - coord[2] ) / ( nb_3d_pts_line - 1. ) * j;
                    coord[0] = ( 1.0 - coord[2] - coord[1] ) / ( nb_3d_pts_line - 1. ) * i;

                    assert ( coord[0] >= 0.0 && coord[0] <= ( 1.0 - coord[1] - coord[2] ) );
                    assert ( coord[1] >= 0.0 && coord[1] <= ( 1.0 - coord[2] ) );
                    assert ( coord[2] >= 0.0 && coord[2] <= 1.0 );

                    std::vector<double> weight ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_x ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_y ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_z ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_xx ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_xy ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_xz ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_yy ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_yz ( fe_tet.get_nb_dof_on_cell ( ) );
                    std::vector<double> weight_zz ( fe_tet.get_nb_dof_on_cell ( ) );

                    fe_tet.N ( coord, weight );
                    fe_tet.N_x ( coord, weight_x );
                    fe_tet.N_y ( coord, weight_y );
                    fe_tet.N_z ( coord, weight_z );
                    fe_tet.N_xx ( coord, weight_xx );
                    fe_tet.N_xy ( coord, weight_xy );
                    fe_tet.N_xz ( coord, weight_xz );
                    fe_tet.N_yy ( coord, weight_yy );
                    fe_tet.N_yz ( coord, weight_yz );
                    fe_tet.N_zz ( coord, weight_zz );

                    // Check
                    double sum = 0.0;
                    double sum_x = 0.0;
                    double sum_y = 0.0;
                    double sum_z = 0.0;
                    double sum_xx = 0.0;
                    double sum_xy = 0.0;
                    double sum_xz = 0.0;
                    double sum_yy = 0.0;
                    double sum_yz = 0.0;
                    double sum_zz = 0.0;

                    for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                    {
                        sum += weight[w];
                        sum_x += weight_x[w];
                        sum_y += weight_y[w];
                        sum_z += weight_z[w];
                        sum_xx += weight_xx[w];
                        sum_xy += weight_xy[w];
                        sum_xz += weight_xz[w];
                        sum_yy += weight_yy[w];
                        sum_yz += weight_yz[w];
                        sum_zz += weight_zz[w];
                    }

                    if ( sum_z > 0.1 )
                    {
                        std::cout << "Degree: " << deg << std::endl;
                        std::cout << coord[0] << " " << coord[1] << "  " << coord[2] << std::endl;
                        for ( int w = 0; w<static_cast < int > ( weight_z.size ( ) ); ++w )
                            std::cout << weight_z[w] << "  ";

                        std::cout << std::endl;
                    }

                    TEST_EQUAL_EPS ( sum, 1.0, 1.e-10 );
                    TEST_EQUAL_EPS ( sum_x, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_y, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_z, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_xx, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_xy, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_xz, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_yy, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_yz, 0.0, 1.e-9 );
                    TEST_EQUAL_EPS ( sum_zz, 0.0, 1.e-9 );
                }

        // Test TRIANGLE

        for ( int j = 0; j < nb_2d_pts_line; ++j )
            for ( int i = 0; i < nb_2d_pts_line - j; ++i )
            {
                std::vector<double> coord ( 2 );

                coord[1] = h_2d*j;
                coord[0] = ( 1.0 - coord[1] ) / ( nb_2d_pts_line - 1. ) * i;

                assert ( coord[0] >= 0.0 && coord[0] <= ( 1.0 - coord[1] ) );
                assert ( coord[1] >= 0.0 && coord[1] <= 1.0 );

                std::vector<double> weight ( fe_tri.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_x ( fe_tri.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_y ( fe_tri.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_xx ( fe_tri.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_xy ( fe_tri.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_yy ( fe_tri.get_nb_dof_on_cell ( ) );

                fe_tri.N ( coord, weight );
                fe_tri.N_x ( coord, weight_x );
                fe_tri.N_y ( coord, weight_y );
                fe_tri.N_xx ( coord, weight_xx );
                fe_tri.N_xy ( coord, weight_xy );
                fe_tri.N_yy ( coord, weight_yy );

                // Check
                double sum = 0.0;
                double sum_x = 0.0;
                double sum_y = 0.0;
                double sum_xx = 0.0;
                double sum_xy = 0.0;
                double sum_yy = 0.0;

                for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                {
                    sum += weight[w];
                    sum_x += weight_x[w];
                    sum_y += weight_y[w];
                    sum_xx += weight_xx[w];
                    sum_xy += weight_xy[w];
                    sum_yy += weight_yy[w];
                }

                TEST_EQUAL_EPS ( sum, 1.0, 1.e-10 );
                TEST_EQUAL_EPS ( sum_x, 0.0, 1.e-9 );
                TEST_EQUAL_EPS ( sum_y, 0.0, 1.e-9 );
                TEST_EQUAL_EPS ( sum_xx, 0.0, 1.e-9 );
                TEST_EQUAL_EPS ( sum_xy, 0.0, 1.e-9 );
                TEST_EQUAL_EPS ( sum_yy, 0.0, 1.e-9 );
            }

        // Test QUADRILATERAL

        for ( int j = 0; j < nb_2d_pts_line; ++j )
            for ( int i = 0; i < nb_2d_pts_line; ++i )
            {
                std::vector<double> coord ( 2 );

                coord[1] = h_2d*j;
                coord[0] = h_2d*i;

                assert ( coord[0] >= 0.0 && coord[0] <= 1.0 );
                assert ( coord[1] >= 0.0 && coord[1] <= 1.0 );

                std::vector<double> weight ( fe_quad.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_x ( fe_quad.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_y ( fe_quad.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_xx ( fe_quad.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_xy ( fe_quad.get_nb_dof_on_cell ( ) );
                std::vector<double> weight_yy ( fe_quad.get_nb_dof_on_cell ( ) );

                fe_quad.N ( coord, weight );
                fe_quad.N_x ( coord, weight_x );
                fe_quad.N_y ( coord, weight_y );
                fe_quad.N_xx ( coord, weight_xx );
                fe_quad.N_xy ( coord, weight_xy );
                fe_quad.N_yy ( coord, weight_yy );

                // Check
                double sum = 0.0;
                double sum_x = 0.0;
                double sum_y = 0.0;
                double sum_xx = 0.0;
                double sum_xy = 0.0;
                double sum_yy = 0.0;

                for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                {
                    sum += weight[w];
                    sum_x += weight_x[w];
                    sum_y += weight_y[w];
                    sum_xx += weight_xx[w];
                    sum_xy += weight_xy[w];
                    sum_yy += weight_yy[w];
                }

                TEST_EQUAL_EPS ( sum, 1.0, 1.e-10 );
                TEST_EQUAL_EPS ( sum_x, 0.0, 1.e-9 );
                TEST_EQUAL_EPS ( sum_y, 0.0, 1.e-9 );
                TEST_EQUAL_EPS ( sum_xx, 0.0, 1.e-9 );
                TEST_EQUAL_EPS ( sum_xy, 0.0, 1.e-9 );
                TEST_EQUAL_EPS ( sum_yy, 0.0, 1.e-9 );
            }

        // Test LINE

        for ( int i = 0; i < nb_1d_pts_line; ++i )
        {
            std::vector<double> coord ( 1 );

            coord[0] = h_2d*i;

            assert ( coord[0] >= 0.0 && coord[0] <= 1.0 );

            std::vector<double> weight ( fe_line.get_nb_dof_on_cell ( ) );
            std::vector<double> weight_x ( fe_line.get_nb_dof_on_cell ( ) );
            std::vector<double> weight_xx ( fe_line.get_nb_dof_on_cell ( ) );

            fe_line.N ( coord, weight );
            fe_line.N_x ( coord, weight_x );
            fe_line.N_xx ( coord, weight_xx );

            // Check
            double sum = 0.0;
            double sum_x = 0.0;
            double sum_xx = 0.0;

            for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
            {
                sum += weight[w];
                sum_x += weight_x[w];
                sum_xx += weight_xx[w];

            }

            TEST_EQUAL_EPS ( sum, 1.0, 1.e-10 );
            TEST_EQUAL_EPS ( sum_x, 0.0, 1.e-9 );
            TEST_EQUAL_EPS ( sum_xx, 0.0, 1.e-9 );
        }

    }

    return 0;
}
