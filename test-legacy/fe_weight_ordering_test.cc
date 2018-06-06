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

/// \author Michael Schick, Martin Baumann

#include <iostream>
#include <vector>
#include <cmath>

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
    // Test to check if the numbering in the fem module is consistent
    // i.e. check if weights are equal to 1.0 on the desired dof

    for ( int deg = 0; deg < 20; ++deg )
    {
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

        for ( int ind_dof = 0; ind_dof < fe_hex.get_nb_dof_on_cell ( ); ++ind_dof )
        {
            std::vector<double> weight ( fe_hex.get_nb_dof_on_cell ( ) );
            fe_hex.N ( fe_hex.get_coord_on_cell ( )[ind_dof], weight );

            for ( int i = 0; i<static_cast < int > ( weight.size ( ) ); ++i )
                if ( std::abs ( weight[i] - 1.0 ) < 1.0e-15 )
                {
                    if ( ind_dof != i )
                        std::cout << "========== TEST FAILED FOR HEXAHEDRON! with Degree " << deg << " ==========" << std::endl;

                    TEST_EQUAL ( ind_dof, i );
                }
        }

        for ( int ind_dof = 0; ind_dof < fe_tet.get_nb_dof_on_cell ( ); ++ind_dof )
        {
            std::vector<double> weight ( fe_tet.get_nb_dof_on_cell ( ) );

            fe_tet.N ( fe_tet.get_coord_on_cell ( )[ind_dof], weight );

            for ( int i = 0; i<static_cast < int > ( weight.size ( ) ); ++i )
                if ( std::abs ( weight[i] - 1.0 ) < 1.0e-15 )
                {
                    if ( ind_dof != i )
                        std::cout << "========== TEST FAILED FOR TETRAHEDRON! with Degree " << deg << " ==========" << std::endl;

                    TEST_EQUAL ( ind_dof, i );
                }
        }

        for ( int ind_dof = 0; ind_dof < fe_tri.get_nb_dof_on_cell ( ); ++ind_dof )
        {
            std::vector<double> weight ( fe_tri.get_nb_dof_on_cell ( ) );
            fe_tri.N ( fe_tri.get_coord_on_cell ( )[ind_dof], weight );

            for ( int i = 0; i<static_cast < int > ( weight.size ( ) ); ++i )
                if ( std::abs ( weight[i] - 1.0 ) < 1.0e-15 )
                {
                    if ( ind_dof != i )
                        std::cout << "========== TEST FAILED FOR TRIANGLE! with Degree " << deg << " ==========" << std::endl;

                    TEST_EQUAL ( ind_dof, i );
                }
        }

        for ( int ind_dof = 0; ind_dof < fe_quad.get_nb_dof_on_cell ( ); ++ind_dof )
        {
            std::vector<double> weight ( fe_quad.get_nb_dof_on_cell ( ) );
            fe_quad.N ( fe_quad.get_coord_on_cell ( )[ind_dof], weight );

            for ( int i = 0; i<static_cast < int > ( weight.size ( ) ); ++i )
                if ( std::abs ( weight[i] - 1.0 ) < 1.0e-15 )
                {
                    if ( ind_dof != i )
                        std::cout << "========== TEST FAILED FOR QUADRILATERAL! with Degree " << deg << " ==========" << std::endl;

                    TEST_EQUAL ( ind_dof, i );
                }
        }

        for ( int ind_dof = 0; ind_dof < fe_line.get_nb_dof_on_cell ( ); ++ind_dof )
        {
            std::vector<double> weight ( fe_line.get_nb_dof_on_cell ( ) );
            fe_line.N ( fe_line.get_coord_on_cell ( )[ind_dof], weight );

            for ( int i = 0; i<static_cast < int > ( weight.size ( ) ); ++i )
                if ( std::abs ( weight[i] - 1.0 ) < 1.0e-15 )
                {
                    if ( ind_dof != i )
                        std::cout << "========== TEST FAILED FOR LINE! with Degree " << deg << " ==========" << std::endl;

                    TEST_EQUAL ( ind_dof, i );
                }
        }
    }

    return 0;
}
