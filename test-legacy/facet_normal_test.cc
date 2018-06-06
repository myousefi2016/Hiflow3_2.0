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

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <mpi.h>

#include "test.h"
#include "hiflow.h"

void quad_normal_test ( );
void tria_normal_test ( );
void hexa_normal_test ( );
void tet_normal_test ( );

const int DEBUG_LEVEL = 3;

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

int main ( int argc, char *argv[] )
{
    MPI_Init ( &argc, &argv );
    std::ofstream debug_file ( "facet_normal_test.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );

    quad_normal_test ( );
    tria_normal_test ( );
    hexa_normal_test ( );
    tet_normal_test ( );

    LogKeeper::get_log ( "debug" ).flush ( );
    MPI_Finalize ( );
    return 0;
}

const double quad_vertex_coords[4][2] = {
    { 0.5, 0.5 },
    { 2., 1. },
    { 1.5, 3. },
    { -0.5, 1.5 }
};

MeshPtr build_quad ( )
{
    MeshDbViewBuilder builder ( 2, 2 );

    int vertices[4];

    for ( int k = 0; k < 4; ++k )
    {
        vertices[k] = builder.add_vertex ( std::vector<double>( &quad_vertex_coords[k][0],
                                           &quad_vertex_coords[k][2] ) );
    }

    builder.add_entity ( 2, std::vector<int>( &vertices[0], &vertices[4] ) );
    return builder.build ( );
}

const double tria_vertex_coords[3][2] = {
    { 0.5, 0.5 },
    { 2., 1. },
    { -0.5, 1.5 }
};

MeshPtr build_tria ( )
{
    MeshDbViewBuilder builder ( 2, 2 );

    int vertices[3];

    for ( int k = 0; k < 3; ++k )
    {
        vertices[k] = builder.add_vertex ( std::vector<double>( &tria_vertex_coords[k][0],
                                           &tria_vertex_coords[k][2] ) );
    }

    builder.add_entity ( 2, std::vector<int>( &vertices[0], &vertices[3] ) );
    return builder.build ( );
}

const double hexa_vertex_coords[8][3] = {
    { 0.5, 0.5, 0.0 },
    { 2.0, 1.0, 0.0 },
    { 1.5, 3.0, 0.0 },
    { -0.5, 1.5, 0.0 },
    { 0.5, 0.5, 1.0 },
    { 2.0, 1.0, 1.0 },
    { 1.5, 3.0, 1.0 },
    { -0.5, 1.5, 1.0 }

                                         // // "unit hex"
                                         // {  0.0, 0.0, 0.0 },
                                         // {  1.0, 0.0, 0.0 },
                                         // {  1.0, 1.0, 0.0 },
                                         // {  0.0, 1.0, 0.0 },
                                         // {  0.0, 0.0, 1.0 },
                                         // {  1.0, 0.0, 1.0 },
                                         // {  1.0, 1.0, 1.0 },
                                         // {  0.0, 1.0, 1.0 }

                                         // // "some other hex"
                                         // {  260.871, 189.540, 352.456 },
                                         // {  261.339, 189.836, 352.338 },
                                         // {  261.517, 189.786, 352.000 },
                                         // {  260.904, 189.318, 352.008 },
                                         // {  260.460, 190.311, 352.436 },
                                         // {  261.221, 190.448, 352.285 },
                                         // {  261.427, 190.680, 351.751 },
                                         // {  260.112, 190.637, 351.520 }
};

MeshPtr build_hexa ( )
{
    MeshDbViewBuilder builder ( 3, 3 );

    int vertices[8];

    for ( int k = 0; k < 8; ++k )
    {
        vertices[k] = builder.add_vertex ( std::vector<double>( &hexa_vertex_coords[k][0],
                                           &hexa_vertex_coords[k][3] ) );
    }

    builder.add_entity ( 3, std::vector<int>( &vertices[0], &vertices[8] ) );
    return builder.build ( );
}

const double tet_vertex_coords[4][3] = {
    { 0.1, 0.1, 0.1 },
    { 0.9, 0.2, 0.2 },
    { 0.3, 0.95, 0.3 },
    { 0.5, 0.5, 0.8 }
};

MeshPtr build_tet ( )
{
    MeshDbViewBuilder builder ( 3, 3 );

    int vertices[4];

    for ( int k = 0; k < 4; ++k )
    {
        vertices[k] = builder.add_vertex ( std::vector<double>( &tet_vertex_coords[k][0],
                                           &tet_vertex_coords[k][3] ) );
    }

    builder.add_entity ( 3, std::vector<int>( &vertices[0], &vertices[4] ) );
    return builder.build ( );
}

class QuadAsm : private AssemblyAssistant<2, double>
{
  public:

    QuadAsm ( )
    {
        // will be resetted!
        quadrature.set_quadrature ( "GaussQuadrilateral", 9 );
    }
    Vec<2, double> normal[4];

    void test_normal_on_facet ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussLine", 3 );

        Entity cell = elem.get_cell ( );
        int i = 0;
        for ( IncidentEntityIterator it = cell.begin_incident ( 1 ); it != cell.end_incident ( 1 ); ++it )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), i );
            initialize_for_facet ( elem, quadrature, i );

            // get normal
            for ( int k = 0; k < 2; ++k )
            {
                LOG_DEBUG ( 1, "normal == (" << n ( k )[0] << ", " << n ( k )[1] << ")" );
            }
            normal[i] = n ( 0 );

            std::vector<double> coords;
            it->get_coordinates ( coords );
            LOG_DEBUG ( 1, "normal == (" << normal[i][0] << ", " << normal[i][1] << ")" );

            assert ( coords.size ( ) == 4 );
            double calculated_normal[2];
            // normal == 1/norm * (-(y1-y0), x1-x0)
            calculated_normal[1] = -1. * ( coords[2] - coords[0] );
            calculated_normal[0] = ( coords[3] - coords[1] );
            // normalize
            double res = 0.0;
            for ( int k = 0; k < 2; ++k )
            {
                res += calculated_normal[k] * calculated_normal[k];
            }
            res = std::sqrt ( res );
            calculated_normal[0] /= res;
            calculated_normal[1] /= res;

            LOG_DEBUG ( 1, "calculated normal == (" << calculated_normal[0] << ", " << calculated_normal[1] << ")\n" );

            for ( int j = 0; j < 2; ++j )
            {
                TEST_EQUAL_EPS ( calculated_normal[j], normal[i][j], 10e-4 );
            }
            ++i;
        }
    }

  private:
    Quadrature<double> quadrature;

};

class TriaAsm : private AssemblyAssistant<2, double>
{
  public:

    TriaAsm ( )
    {
        // will be resetted!
        quadrature.set_quadrature ( "GaussTriangle", 4 );
    }
    Vec<2, double> normal[4];

    void test_normal_on_facet ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussLine", 3 );

        Entity cell = elem.get_cell ( );
        int i = 0;
        for ( IncidentEntityIterator it = cell.begin_incident ( 1 ); it != cell.end_incident ( 1 ); ++it )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), i );
            initialize_for_facet ( elem, quadrature, i );

            // get normal
            normal[i] = n ( 0 );
            // get normal
            for ( int k = 0; k < 2; ++k )
            {
                LOG_DEBUG ( 1, "normal == (" << n ( k )[0] << ", " << n ( k )[1] << ")" );
            }
            normal[i] = n ( 0 );

            std::vector<double> coords;
            it->get_coordinates ( coords );
            LOG_DEBUG ( 1, "normal == (" << normal[i][0] << ", " << normal[i][1] << ")" );

            assert ( coords.size ( ) == 4 );
            double calculated_normal[2];
            // normal == 1/norm * (-(y1-y0), x1-x0)
            calculated_normal[1] = -1. * ( coords[2] - coords[0] );
            calculated_normal[0] = ( coords[3] - coords[1] );
            // normalize
            double res = 0.0;
            for ( int k = 0; k < 2; ++k )
            {
                res += calculated_normal[k] * calculated_normal[k];
            }
            res = std::sqrt ( res );
            calculated_normal[0] /= res;
            calculated_normal[1] /= res;

            LOG_DEBUG ( 1, "calculated normal == (" << calculated_normal[0] << ", " << calculated_normal[1] << ")\n" );

            for ( int j = 0; j < 2; ++j )
            {
                TEST_EQUAL_EPS ( calculated_normal[j], normal[i][j], 10e-4 );
            }

            ++i;
        }
    }

  private:
    Quadrature<double> quadrature;

};

class HexaAsm : private AssemblyAssistant<3, double>
{
  public:

    HexaAsm ( )
    {
        // will be resetted!
        quadrature.set_quadrature ( "GaussHexahedron", 27 );
    }
    Vec<3, double> normal[6];

    void test_normal_on_facet ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussQuadrilateral", 9 );

        Entity cell = elem.get_cell ( );
        int i = 0;
        for ( IncidentEntityIterator it = cell.begin_incident ( 2 ); it != cell.end_incident ( 2 ); ++it )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), i );
            initialize_for_facet ( elem, quadrature, i );

            // get normal
            normal[i] = n ( 0 );
            LOG_DEBUG ( 1, "normal == (" << normal[i][0] << ", " << normal[i][1] << ", " << normal[i][2] << ")" );

            std::vector<double> coords;
            it->get_coordinates ( coords );
            assert ( coords.size ( ) == 12 );

            LOG_DEBUG ( 4, "coords == " << string_from_range ( coords.begin ( ), coords.end ( ) ) );

            std::vector<double> coords0 ( 3 );
            coords0[0] = coords[3] - coords[0];
            coords0[1] = coords[4] - coords[1];
            coords0[2] = coords[5] - coords[2];

            std::vector<double> coords1 ( 3 );
            coords1[0] = coords[9] - coords[0];
            coords1[1] = coords[10] - coords[1];
            coords1[2] = coords[11] - coords[2];

            Vec<3, double> coordinates[2];
            coordinates[0] = Vec<3, double>( coords0 );
            coordinates[1] = Vec<3, double>( coords1 );

            Vec<3, double> calculated_normal;
            calculated_normal = cross ( coordinates[0], coordinates[1] );

            // normalize
            double res = 0.0;
            for ( int k = 0; k < 3; ++k )
            {
                res += calculated_normal[k] * calculated_normal[k];
            }
            res = std::sqrt ( res );
            calculated_normal[0] /= res;
            calculated_normal[1] /= res;
            calculated_normal[2] /= res;

            LOG_DEBUG ( 1, "calculated normal == (" << calculated_normal[0] << ", " << calculated_normal[1] << ", " << calculated_normal[2] << ")\n" );

            for ( int j = 0; j < 3; ++j )
            {
                TEST_EQUAL_EPS ( calculated_normal[j], normal[i][j], 10e-4 );
            }

            ++i;
        }
    }

  private:
    Quadrature<double> quadrature;

};

class TetAsm : private AssemblyAssistant<3, double>
{
  public:

    TetAsm ( )
    {
        // will be resetted!
        quadrature.set_quadrature ( "GaussTetrahedron", 24 );
    }
    Vec<3, double> normal[4];

    void test_normal_on_facet ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussTriangle", 4 );

        Entity cell = elem.get_cell ( );
        int i = 0;
        for ( IncidentEntityIterator it = cell.begin_incident ( 2 ); it != cell.end_incident ( 2 ); ++it )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), i );
            initialize_for_facet ( elem, quadrature, i );

            // get normal
            normal[i] = n ( 0 );
            LOG_DEBUG ( 1, "normal == (" << normal[i][0] << ", " << normal[i][1] << ", " << normal[i][2] << ")" );

            std::vector<double> coords;
            it->get_coordinates ( coords );
            assert ( coords.size ( ) == 9 );

            LOG_DEBUG ( 4, "coords == " << string_from_range ( coords.begin ( ), coords.end ( ) ) );

            std::vector<double> coords0 ( 3 );
            coords0[0] = coords[3] - coords[0];
            coords0[1] = coords[4] - coords[1];
            coords0[2] = coords[5] - coords[2];

            std::vector<double> coords1 ( 3 );
            coords1[0] = coords[6] - coords[0];
            coords1[1] = coords[7] - coords[1];
            coords1[2] = coords[8] - coords[2];

            Vec<3, double> coordinates[2];
            coordinates[0] = Vec<3, double>( coords0 );
            coordinates[1] = Vec<3, double>( coords1 );

            Vec<3, double> calculated_normal;
            calculated_normal = cross ( coordinates[0], coordinates[1] );

            // normalize
            double res = 0.0;
            for ( int k = 0; k < 3; ++k )
            {
                res += calculated_normal[k] * calculated_normal[k];
            }
            res = std::sqrt ( res );
            calculated_normal[0] /= res;
            calculated_normal[1] /= res;
            calculated_normal[2] /= res;

            LOG_DEBUG ( 1, "calculated normal == (" << calculated_normal[0] << ", " << calculated_normal[1] << ", " << calculated_normal[2] << ")\n" );

            for ( int j = 0; j < 3; ++j )
            {
                TEST_EQUAL_EPS ( calculated_normal[j], normal[i][j], 10e-4 );
            }

            ++i;
        }
    }

  private:
    Quadrature<double> quadrature;

};

void quad_normal_test ( )
{
    LOG_DEBUG ( 1, "=== Begin Quad Test ====" );

    const int deg = 2;

    MeshPtr quad_mesh = build_quad ( );
    VectorSpace<double> space;
    space.Init ( std::vector<int>( 1, deg ), *quad_mesh );
    Element<double> elem ( space, 0 );

    QuadAsm quad_asm;
    quad_asm.test_normal_on_facet ( elem );

    LOG_DEBUG ( 1, "=== End Quad Test ====" );
}

void tria_normal_test ( )
{
    LOG_DEBUG ( 1, "=== Begin Triangle Test ====" );

    const int deg = 2;

    MeshPtr triangle_mesh = build_tria ( );

    VectorSpace<double> space;
    space.Init ( std::vector<int>( 1, deg ), *triangle_mesh );
    Element<double> elem ( space, 0 );

    TriaAsm tria_asm;
    tria_asm.test_normal_on_facet ( elem );

    LOG_DEBUG ( 1, "=== End Triangle Test ====" );
}

void hexa_normal_test ( )
{
    LOG_DEBUG ( 1, "=== Begin Hexa Test ====" );

    const int deg = 2;

    MeshPtr hexa_mesh = build_hexa ( );
    VectorSpace<double> space;
    space.Init ( std::vector<int>( 1, deg ), *hexa_mesh );
    Element<double> elem ( space, 0 );

    VtkWriter writer;
    writer.write ( "hexa_mesh.vtu", *hexa_mesh );

    HexaAsm hexa_asm;
    hexa_asm.test_normal_on_facet ( elem );

    LOG_DEBUG ( 1, "=== End Hexa Test ====" );
}

void tet_normal_test ( )
{
    LOG_DEBUG ( 1, "=== Begin Tet Test ====" );

    const int deg = 2;

    MeshPtr tet_mesh = build_tet ( );
    VectorSpace<double> space;
    space.Init ( std::vector<int>( 1, deg ), *tet_mesh );
    Element<double> elem ( space, 0 );

    VtkWriter writer;
    writer.write ( "tet_mesh.vtu", *tet_mesh );

    TetAsm tet_asm;
    tet_asm.test_normal_on_facet ( elem );

    LOG_DEBUG ( 1, "=== End Tet Test ====" );
}
