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

#include <fstream>

#include "test.h"
#include "hiflow.h"

using namespace hiflow;
using namespace hiflow::mesh;

void test_triangle ( );
void test_quad ( );
void test_tet ( );
void test_hex ( );

const int DEBUG_LEVEL = 1;

int main ( int argc, char *argv[] )
{
    std::ofstream debug_file ( "cell_type_test_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_file );

    test_triangle ( );
    test_quad ( );
    test_tet ( );
    test_hex ( );

    return 0;
}

namespace
{
    const int num_triangle_sub_cells = 11;
    const int triangle_sub_cell_parent_edges[num_triangle_sub_cells][3] = {
        { -1, -1, -1 },
        { 0, -1, 2 },
        { 1, -1, 0 },
        { 2, -1, 1 },
        { -1, -1, -1 },
        { -1, 1, -1 },
        { -1, 1, -1 },
        { -1, -1, 2 },
        { -1, -1, 2 },
        { 0, -1, -1 },
        { 0, -1, -1 }
    };
}

void test_triangle ( )
{
    const CellType& cell_type = CellType::get_instance ( CellType::TRIANGLE );
    LOG_DEBUG ( 1, "=== Triangle ===\n" << cell_type );

    const TDim facet_tdim = 1;
    for ( int i = 0; i < num_triangle_sub_cells; ++i )
    {
        const std::vector<int> edges_of_sub_cell = cell_type.sub_entities_of_cell ( facet_tdim, i );
        for ( int e = 0; e < 3; ++e )
        {
            const int edge = edges_of_sub_cell[e];
            TEST_EQUAL ( cell_type.regular_parent ( facet_tdim, edge ), triangle_sub_cell_parent_edges[i][e] );
        }
    }
}

void test_quad ( )
{
    // TODO: write tests (see beginning of table in Staffan's notes)
    const CellType& cell_type = CellType::get_instance ( CellType::QUADRILATERAL );
    LOG_DEBUG ( 1, "=== Quadrilateral ===\n" << cell_type );
}

void test_tet ( )
{
    // TODO: write tests (transfer from old cell_types)
    const CellType& cell_type = CellType::get_instance ( CellType::TETRAHEDRON );
    LOG_DEBUG ( 1, "=== Tetrahedron ===\n" << cell_type );
}

namespace
{
    const int num_hex_sub_cells = 9;
    const int hexahedron_sub_cell_parent_facets[num_hex_sub_cells][6] = {
        { -1, -1, -1, -1, -1, -1 }, // root cell
        { 0, 1, 2, -1, -1, -1 }, // bottom front left
        { 0, 1, -1, 3, -1, -1 }, // bottom front right
        { 0, -1, -1, 3, 4, -1 }, // bottom back right
        { 0, -1, 2, -1, 4, -1 }, // bottom back left
        {-1, 1, 2, -1, -1, 5 }, // top front left
        {-1, 1, -1, 3, -1, 5 }, // top front right
        {-1, -1, -1, 3, 4, 5 }, // top back right
        {-1, -1, 2, -1, 4, 5 } // top back left
    };
}

void test_hex ( )
{
    const CellType& cell_type = CellType::get_instance ( CellType::HEXAHEDRON );
    LOG_DEBUG ( 1, "=== Hexahedron ===\n" << cell_type );

    const TDim facet_tdim = 2;
    for ( int i = 0; i < num_hex_sub_cells; ++i )
    {
        const std::vector<int> facets_of_sub_cell = cell_type.sub_entities_of_cell ( facet_tdim, i );
        for ( int f = 0; f < 6; ++f )
        {
            const int facet = facets_of_sub_cell[f];
            TEST_EQUAL ( cell_type.regular_parent ( facet_tdim, facet ), hexahedron_sub_cell_parent_facets[i][f] );
        }
    }
}
