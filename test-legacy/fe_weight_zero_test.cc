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

#include "../src/fem/felagrange_line.h"
#include "../src/fem/felagrange_hex.h"
#include "../src/fem/felagrange_tet.h"
#include "../src/fem/felagrange_tri.h"
#include "../src/fem/felagrange_quad.h"
#include "test.h"

using namespace hiflow::doffem;

int main ( int argc, char** argv )
{
    // Test if shapefunctions are zero on dofs, which correspond to other shapefunctions

    for ( int deg = 1; deg < 8; ++deg )
    {
        std::cout << "Degree: " << deg << std::endl;

        FELagrangeLine<double> fe_line;
        FELagrangeHex<double> fe_hex;
        FELagrangeTet<double> fe_tet;
        FELagrangeTri<double> fe_tri;
        FELagrangeQuad<double> fe_quad;

        fe_line.set_fe_deg ( deg );
        fe_hex.set_fe_deg ( deg );
        fe_tet.set_fe_deg ( deg );
        fe_tri.set_fe_deg ( deg );
        fe_quad.set_fe_deg ( deg );

        fe_line.init ( );
        fe_hex.init ( );
        fe_tet.init ( );
        fe_tri.init ( );
        fe_quad.init ( );

        // Test HEXAHEDRON

        int nb_dof_line = deg + 1;

        for ( int k = 0; k < nb_dof_line; ++k )
            for ( int j = 0; j < nb_dof_line; ++j )
                for ( int i = 0; i < nb_dof_line; ++i )
                {
                    int shape_fct_index = ( i + j * nb_dof_line
                            + k * nb_dof_line * nb_dof_line );

                    std::vector<double> coord ( 3 );

                    double h = 1.0 / deg;

                    coord[0] = i*h;
                    coord[1] = j*h;
                    coord[2] = k*h;

                    std::vector<double> weight ( fe_hex.get_nb_dof_on_cell ( ) );

                    fe_hex.N ( coord, weight );

                    for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                        if ( w != shape_fct_index ) TEST_EQUAL_EPS ( weight[w], 0.0, 1.e-10 );
                }

        // Test TETRAHEDRON
        nb_dof_line = deg + 1;
        int offset = 0;

        for ( int k = 0; k < nb_dof_line; ++k )
            for ( int j = 0; j < nb_dof_line - k; ++j )
            {
                for ( int i = 0; i < nb_dof_line - k - j; ++i )
                {
                    int shape_fct_index = i + offset;

                    std::vector<double> coord ( 3 );

                    double h = 1.0 / deg;

                    coord[2] = h*k;
                    coord[1] = h*j;
                    coord[0] = h*i;

                    std::vector<double> weight ( fe_tet.get_nb_dof_on_cell ( ) );

                    fe_tet.N ( coord, weight );

                    for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                        if ( w != shape_fct_index ) TEST_EQUAL_EPS ( weight[w], 0.0, 1.e-9 );

                }
                offset += nb_dof_line - k - j;
            }

        // Test TRIANGLE

        nb_dof_line = deg + 1;
        offset = 0;

        for ( int j = 0; j < nb_dof_line; ++j )
        {
            for ( int i = 0; i < nb_dof_line - j; ++i )
            {
                int shape_fct_index = i + offset;

                double h = 1.0 / deg;

                std::vector<double> coord ( 2 );

                coord[1] = h*j;
                coord[0] = h*i;

                std::vector<double> weight ( fe_tri.get_nb_dof_on_cell ( ) );

                fe_tri.N ( coord, weight );

                for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                    if ( w != shape_fct_index ) TEST_EQUAL_EPS ( weight[w], 0.0, 1.e-9 );

            }
            offset += nb_dof_line - j;
        }

        // Test QUADRILATERAL

        nb_dof_line = deg + 1;

        for ( int j = 0; j < nb_dof_line; ++j )
            for ( int i = 0; i < nb_dof_line; ++i )
            {
                int shape_fct_index = i + j*nb_dof_line;

                std::vector<double> coord ( 2 );

                double h = 1.0 / deg;

                coord[1] = h*j;
                coord[0] = h*i;

                std::vector<double> weight ( fe_quad.get_nb_dof_on_cell ( ) );

                fe_quad.N ( coord, weight );

                for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                    if ( w != shape_fct_index ) TEST_EQUAL_EPS ( weight[w], 0.0, 1.e-9 );
            }

        // Test LINE

        nb_dof_line = deg + 1;

        for ( int j = 0; j < nb_dof_line; ++j )
        {
            int shape_fct_index = j;

            std::vector<double> coord ( 1 );

            double h = 1.0 / deg;

            coord[0] = h*j;

            std::vector<double> weight ( fe_line.get_nb_dof_on_cell ( ) );

            fe_line.N ( coord, weight );

            for ( int w = 0; w < static_cast < int > ( weight.size ( ) ); ++w )
                if ( w != shape_fct_index ) TEST_EQUAL_EPS ( weight[w], 0.0, 1.e-9 );
        }
    }

    return 0;
}
