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

#include <cassert>
#include <fstream>
#include <cstdarg>
#include <map>
#include <utility>

#include "test.h"
#include "common/log.h"
#include "mesh/mesh_db_view.h"
#include "mesh/interface.h"
#include "mesh/writer.h"

/// author Staffan Ronnas

const int DEBUG_LEVEL = 2;

using namespace hiflow;
using namespace hiflow::mesh;

MeshPtr build_unit_square ( );
void test_unit_square ( );
void test_regular_refined_unit_square ( );
void test_irregular_refined_unit_square ( );

int main ( int argc, char *argv[] )
{
    std::ofstream debug_file ( "mesh_interface_test.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_file );

    // test of simple (non-refined) unit square: only boundary facets.
    test_unit_square ( );

    // test of level 1 regularly refined unit square: only boundary and regular facets.
    test_regular_refined_unit_square ( );

    // test of irregularly refined unit square: boundary, regular and irregular facets.

    // TODO: this test needs to be updated to new slave facet convention.
    //    test_irregular_refined_unit_square();

    return 0;
}

MeshPtr build_unit_square ( )
{
    MeshDbViewBuilder builder ( 2, 2 );

    int vertices[4];
    const double vertex_coords[4][2] = {
        {-1., -1. },
        { 1., -1. },
        { 1., 1. },
        {-1., 1. }
    };

    for ( int k = 0; k < 4; ++k )
    {
        vertices[k] = builder.add_vertex ( std::vector<double>( &vertex_coords[k][0],
                                           &vertex_coords[k][2] ) );
    }

    builder.add_entity ( 2, std::vector<int>( &vertices[0], &vertices[4] ) );
    return builder.build ( );
}

void test_unit_square ( )
{
    MeshPtr mesh = build_unit_square ( );
    InterfaceList interface_list = InterfaceList::create ( mesh );

    LOG_DEBUG ( 1, "=== Begin Unit square (no refinement) ====" );

    TEST_EQUAL ( interface_list.size ( ), 4 );

    for ( InterfaceList::const_iterator it = interface_list.begin ( );
          it != interface_list.end ( ); ++it )
    {
        LOG_DEBUG ( 2, *it );

        // all interfaces should be boundary
        TEST_EQUAL ( it->num_slaves ( ), 0 );

        InterfacePattern pattern = compute_interface_pattern ( mesh, *it );
        LOG_DEBUG ( 2, pattern );
        TEST_EQUAL ( pattern.master_facet_number ( ), it->master_facet_number ( ) );
        TEST_EQUAL ( pattern.num_slaves ( ), 0 );
        TEST_EQUAL ( pattern.orientation ( ), -1 );
    }
    LOG_DEBUG ( 1, "=== End Unit square (no refinement) ====" );
}

void test_regular_refined_unit_square ( )
{
    MeshPtr mesh = build_unit_square ( );
    mesh = mesh->refine ( );
    InterfaceList interface_list = InterfaceList::create ( mesh );

    LOG_DEBUG ( 1, "=== Begin Unit square (1 level of regular refinement) ====" );
    //    LOG_DEBUG(1, interface_list);

    TEST_EQUAL ( interface_list.size ( ), 12 );

    // Zero or one slave cells / interface
    // -1  -> exists, but no slave cells
    // -10 -> does not exist
    const int correct_slave_cells[4][4] = {
        { -1, 1, -10, -1 }, // cell 0
        { -1, -1, -10, -10 }, // cell 1
        { 1, -1, -1, -10 }, // cell 2
        { 0, 2, -1, -1 }
    }; // cell 3

    const int correct_orientations[4][4] = {
        { -1, 0, -10, -1 }, // cell 0
        { -1, -1, -10, -10 }, // cell 1
        { 3, -1, -1, -10 }, // cell 2
        { 3, 0, -1, -1 }
    }; // cell 3

    const int correct_slave_facets[4][4] = {
        { -1, 3, -10, -1 }, // cell 0
        { -1, -1, -10, -10 }, // cell 1
        { 2, -1, -1, -10 }, // cell 2
        { 2, 3, -1, -1 }
    }; // cell 3

    for ( InterfaceList::const_iterator it = interface_list.begin ( );
          it != interface_list.end ( ); ++it )
    {
        // check interface object
        LOG_DEBUG ( 2, *it );
        const int master_cell = it->master_index ( );
        const int master_facet_number = it->master_facet_number ( );
        TEST_NOT_EQUAL ( correct_slave_cells[master_cell][master_facet_number], -10 );
        if ( correct_slave_cells[master_cell][master_facet_number] == -1 )
        {
            TEST_EQUAL ( it->num_slaves ( ), 0 );
        }
        else
        {
            TEST_EQUAL ( it->num_slaves ( ), 1 );
            TEST_EQUAL ( it->slave_index ( 0 ), correct_slave_cells[master_cell][master_facet_number] );
        }

        // check pattern
        InterfacePattern pattern = compute_interface_pattern ( mesh, *it );
        LOG_DEBUG ( 2, pattern );
        TEST_EQUAL ( pattern.master_facet_number ( ), master_facet_number ); // always true

        if ( correct_slave_cells[master_cell][master_facet_number] == -1 )
        {
            TEST_EQUAL ( pattern.num_slaves ( ), 0 );
        }
        else
        {
            TEST_EQUAL ( pattern.num_slaves ( ), 1 );
            TEST_EQUAL ( pattern.slave_facet_number ( 0 ),
                         correct_slave_facets[master_cell][master_facet_number] );
        }
        TEST_EQUAL ( pattern.orientation ( ), correct_orientations[master_cell][master_facet_number] );
    }

    LOG_DEBUG ( 1, "=== End Unit square (1 level of regular refinement) ====" );
}

typedef std::map<std::pair<int, int>, Interface> InterfaceMap;
typedef std::map<std::pair<int, int>, InterfacePattern> PatternMap;

/// \brief Add an interface object to an InterfaceMap.
/// The variable arguments should list the cell indices of the slaves, one after the other.
/// num_slaves extra arguments should hence be provided (num_slaves can be 0).

void add_interface ( InterfaceMap& interfaces, int master_cell, int master_facet, int num_slaves, ... )
{
    using std::make_pair;
    interfaces.insert ( make_pair ( make_pair ( master_cell, master_facet ),
                                    Interface ( master_cell, master_facet ) ) );

    va_list vl;
    va_start ( vl, num_slaves );

    for ( int i = 0; i < num_slaves; ++i )
    {
        int slave_cell = va_arg ( vl, int );
        interfaces.find ( make_pair ( master_cell, master_facet ) )->second.add_slave ( slave_cell, i ); //TODO is i correct here??
    }

    va_end ( vl );
}

/// \brief Add an InterfacePattern object to a PatternMap.
/// The variable arguments should list the facet numbers of the slave facets, one after the other.
/// num_slaves extra arguments should hence be provided (num_slaves can be 0).

void add_pattern ( PatternMap& patterns, int master_cell, int master_facet, int orientation, int num_slaves, ... )
{
    using std::make_pair;
    patterns.insert ( make_pair ( make_pair ( master_cell, master_facet ),
                                  InterfacePattern ( ) ) );
    patterns[make_pair ( master_cell, master_facet )].set_master_facet_number ( master_facet );
    patterns[make_pair ( master_cell, master_facet )].set_orientation ( orientation );

    va_list vl;
    va_start ( vl, num_slaves );

    for ( int i = 0; i < num_slaves; ++i )
    {
        int slave_facet = va_arg ( vl, int );

        if ( num_slaves == 1 )
        {
            // regular interface
            patterns[make_pair ( master_cell, master_facet )].add_regular_slave_facet ( slave_facet );
        }
        else
        {
            // irregular inteface
            // TODO: find correct arguments for call below.
            patterns[make_pair ( master_cell, master_facet )].add_irregular_slave_facet ( -1, -2 );
        }
    }

    va_end ( vl );
}

/// List correct interfaces and patterns for the one-irregular square test case.

void build_correct_interfaces_and_patterns_irregular_square ( InterfaceMap& interfaces,
                                                              PatternMap& patterns )
{
    add_interface ( interfaces, 0, 0, 0 );
    add_interface ( interfaces, 0, 1, 1, 1 );
    add_interface ( interfaces, 0, 3, 0 );
    add_interface ( interfaces, 1, 0, 0 );
    add_interface ( interfaces, 1, 1, 0 );
    add_interface ( interfaces, 1, 2, 2, 3, 4 );
    add_interface ( interfaces, 2, 0, 1, 0 );
    add_interface ( interfaces, 2, 1, 2, 3, 6 );
    add_interface ( interfaces, 2, 2, 0 );
    add_interface ( interfaces, 2, 3, 0 );
    add_interface ( interfaces, 3, 1, 1, 4 );
    add_interface ( interfaces, 4, 1, 0 );
    add_interface ( interfaces, 5, 0, 1, 4 );
    add_interface ( interfaces, 5, 1, 0 );
    add_interface ( interfaces, 5, 2, 0 );
    add_interface ( interfaces, 6, 0, 1, 3 );
    add_interface ( interfaces, 6, 1, 1, 5 );
    add_interface ( interfaces, 6, 2, 0 );

    add_pattern ( patterns, 0, 0, -1, 0 );
    add_pattern ( patterns, 0, 1, 0, 1, 3 );
    add_pattern ( patterns, 0, 3, -1, 0 );
    add_pattern ( patterns, 1, 0, -1, 0 );
    add_pattern ( patterns, 1, 1, -1, 0 );
    add_pattern ( patterns, 1, 2, 1, 2, 4, 8 );
    add_pattern ( patterns, 2, 0, 3, 1, 2 );
    add_pattern ( patterns, 2, 1, 0, 2, 7, 15 );
    add_pattern ( patterns, 2, 2, -1, 0 );
    add_pattern ( patterns, 2, 3, -1, 0 );
    add_pattern ( patterns, 3, 1, 0, 1, 3 );
    add_pattern ( patterns, 4, 1, -1, 0 );
    add_pattern ( patterns, 5, 0, 3, 1, 2 );
    add_pattern ( patterns, 5, 1, -1, 0 );
    add_pattern ( patterns, 5, 2, -1, 0 );
    add_pattern ( patterns, 6, 0, 3, 1, 2 );
    add_pattern ( patterns, 6, 1, 0, 1, 3 );
    add_pattern ( patterns, 6, 2, -1, 0 );

}

/// \brief Checks an Interface object against the one in the provided map, without taking
/// into account the order of the slave cells.

std::vector<int> check_interface ( const InterfaceMap& correct_interfaces, const Interface& test_interface )
{
    using std::make_pair;
    const int master_cell = test_interface.master_index ( );
    const int master_facet = test_interface.master_facet_number ( );

    // check if interface exists and should not
    InterfaceMap::const_iterator it = correct_interfaces.find ( make_pair ( master_cell, master_facet ) );
    const bool interface_missing = ( it == correct_interfaces.end ( ) );
    TEST_NOT_EQUAL ( interface_missing, true );

    std::vector<int> order;

    // check that slave cells correspond (order not important)
    TEST_EQUAL ( test_interface.num_slaves ( ), it->second.num_slaves ( ) );
    for ( Interface::const_iterator s_it = test_interface.begin ( );
          s_it != test_interface.end ( ); ++s_it )
    {
        LOG_DEBUG ( 3, "Sought slave cell = " << *s_it );
        Interface::const_iterator find_it = find ( it->second.begin ( ), it->second.end ( ), *s_it );
        const bool slave_found = ( find_it != it->second.end ( ) );
        TEST_EQUAL ( slave_found, true );
        order.push_back ( find_it - it->second.begin ( ) );
    }

    return order;
}

/// \brief Checks one InterfacePattern object against the one in the
/// provided map, using the given order to permute the slave facets of the object in the map.

void check_pattern ( int master_cell, const PatternMap& correct_patterns, const InterfacePattern& test_pattern, const std::vector<int>& order )
{
    using std::make_pair;
    const int master_facet = test_pattern.master_facet_number ( );

    PatternMap::const_iterator it = correct_patterns.find ( make_pair ( master_cell, master_facet ) );
    assert ( it != correct_patterns.end ( ) );
    TEST_EQUAL ( test_pattern.master_facet_number ( ), it->second.master_facet_number ( ) );
    TEST_EQUAL ( test_pattern.num_slaves ( ), it->second.num_slaves ( ) );

    // check that slave facets correspond
    // NB: no check is performed that relative ordering corresponds between interface and pattern
    for ( int k = 0; k < test_pattern.num_slaves ( ); ++k )
    {
        TEST_EQUAL ( test_pattern.slave_facet_number ( k ), it->second.slave_facet_number ( order.at ( k ) ) );
    }

    TEST_EQUAL ( test_pattern.orientation ( ), it->second.orientation ( ) );
}

void test_irregular_refined_unit_square ( )
{
    MeshPtr mesh = build_unit_square ( );
    mesh = mesh->refine ( );
    int ref_list[4] = { -10, -10, 0, -10 };
    std::vector<int> refinements ( &ref_list[0], &ref_list[4] );
    mesh = mesh->refine ( refinements );

    VtkWriter writer;
    writer.write ( "mesh_interface_test_irregular_square.vtu", *mesh );

    LOG_DEBUG ( 1, "=== Begin Unit square (1-irregular refinement of cell 2) ====" );

    InterfaceList interface_list = InterfaceList::create ( mesh );
    TEST_EQUAL ( interface_list.size ( ), 18 );

    InterfaceMap correct_interfaces;
    PatternMap correct_patterns;
    build_correct_interfaces_and_patterns_irregular_square ( correct_interfaces, correct_patterns );
    assert ( correct_interfaces.size ( ) == 18 );
    assert ( correct_patterns.size ( ) == 18 );

    for ( InterfaceList::const_iterator it = interface_list.begin ( );
          it != interface_list.end ( ); ++it )
    {
        // check interface object
        LOG_DEBUG ( 2, *it );
        const std::vector<int> order = check_interface ( correct_interfaces, *it );

        // check pattern
        InterfacePattern pattern = compute_interface_pattern ( mesh, *it );
        LOG_DEBUG ( 2, pattern );
        check_pattern ( it->master_index ( ), correct_patterns, pattern, order );
    }

    LOG_DEBUG ( 1, "=== End Unit square (1-irregular refinement of cell 2) ====" );
}
