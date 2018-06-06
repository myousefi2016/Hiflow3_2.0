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

/// \author Thomas Gengenbach, Staffan Ronnas

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <mpi.h>

#include "test.h"
#include "hiflow.h"
#include "mesh/cell_type_definitions.cc"

void quad_facet_test ( );
void tria_facet_test ( );
void hexa_facet_test ( );
void tetra_facet_test ( );

double distance ( const int dim, const double* point1, const double* point2 );

const int DEBUG_LEVEL = 3;

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

int main ( int argc, char *argv[] )
{
    MPI_Init ( &argc, &argv );
    std::ofstream debug_file ( "hp_refinement_test.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cerr );

    quad_facet_test ( );
    tria_facet_test ( );
    hexa_facet_test ( );
    tetra_facet_test ( );

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

const double hex_vertex_coords[8][3] = {
                                        // // "unit hex"
                                        // {  0.0, 0.0, 0.0 },
                                        // {  1.0, 0.0, 0.0 },
                                        // {  1.0, 1.0, 0.0 },
                                        // {  0.0, 1.0, 0.0 },
                                        // {  0.0, 0.0, 1.0 },
                                        // {  1.0, 0.0, 1.0 },
                                        // {  1.0, 1.0, 1.0 },
                                        // {  0.0, 1.0, 1.0 }

                                        // // "some hex"
                                        // {  3.2, 0.15, -0.1 },
                                        // {  4.0, 0.15, 0.7 },
                                        // {  3.2, 0.15, 1.5 },
                                        // {  2.4, 0.15, 0.7 },
                                        // {  3.2, 0.85, -0.1 },
                                        // {  4.0, 0.85, 0.7 },
                                        // {  3.2, 0.85, 1.5 },
                                        // {  2.4, 0.85, 0.7 }

                                        // "some other hex"
    { 260.871, 189.540, 352.456 },
    { 261.339, 189.836, 352.338 },
    { 261.517, 189.786, 352.000 },
    { 260.904, 189.318, 352.008 },
    { 260.460, 190.311, 352.436 },
    { 261.221, 190.448, 352.285 },
    { 261.427, 190.680, 351.751 },
    { 260.112, 190.637, 351.520 }

};

MeshPtr build_hexa ( )
{
    MeshDbViewBuilder builder ( 3, 3 );

    int vertices[8];

    for ( int k = 0; k < 8; ++k )
    {
        vertices[k] = builder.add_vertex ( std::vector<double>( &hex_vertex_coords[k][0],
                                           &hex_vertex_coords[k][3] ) );
    }

    builder.add_entity ( 3, std::vector<int>( &vertices[0], &vertices[8] ) );
    return builder.build ( );
}

const double tet_vertex_coords[4][3] = {
                                        // // "unit tet"
                                        // {  0.0, 0.0, 0.0 },
                                        // {  1.0, 0.0, 0.0 },
                                        // {  0.0, 1.0, 0.0 },
                                        // {  0.0, 0.0, 1.0 }

                                        // "some tet"
    { 21.4692, 17.0101, -31.0927 },
    { 21.5796, 17.0332, -31.0916 },
    { 21.5578, 17.0780, -31.0874 },
    { 21.4641, 17.0103, -30.9638 }
};

MeshPtr build_tetra ( )
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
        quadrature.set_quadrature ( "GaussQuadrilateral", 9 );
    }

    void regression_output ( const Element<double>& elem )
    {
        std::ofstream os ( "facet_integral_test_quad_regression.log" );

        initialize_for_element ( elem, quadrature );

        os << "num_dofs_total = " << num_dofs_total ( ) << "\n";

        for ( int s = 0; s < num_dofs ( 0 ); ++s )
        {
            os << "dof(" << s << ") = " << dof_index ( s, 0 ) << "\n";
        }

        for ( int q = 0; q < num_quadrature_points ( ); ++q )
        {
            os << "q = " << q << ": {"
                    << " phi(1) = " << phi ( 1, q )
                    << ", grad_phi(2)[1] = " << grad_phi ( 2, q )[1]
                    << ", detJ = " << detJ ( q )
                    << ", w = " << w ( q ) << "\n";
        }
        os.close ( );
    }

    void test_facet_lengths ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussLine", 3 );

        for ( int f = 0; f < 4; ++f )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), f );

            initialize_for_facet ( elem, quadrature, f );

            const int v0 = f;
            const int v1 = ( f + 1 ) % 4;
            LOG_DEBUG ( 3, "v0 = " << v0 << ", v1 = " << v1 );
            LOG_DEBUG ( 3, "x0 = " << quad_vertex_coords[v0][0] << ", x1 = " << quad_vertex_coords[v1][0] );
            LOG_DEBUG ( 3, "y0 = " << quad_vertex_coords[v0][1] << ", y1 = " << quad_vertex_coords[v1][1] );
            const double correct_length =
                    std::sqrt (
                                std::pow ( ( quad_vertex_coords[v1][0] - quad_vertex_coords[v0][0] ), 2. ) +
                                std::pow ( ( quad_vertex_coords[v1][1] - quad_vertex_coords[v0][1] ), 2. ) );

            LOG_DEBUG ( 2, "Exact length of quad edge " << f << " = " << correct_length );

            double approx_length = 0.;
            for ( int q = 0; q < num_quadrature_points ( ); ++q )
            {
                LOG_DEBUG ( 3, "w(" << q << ") = " << w ( q ) << ", ds(" << q << ") = " << ds ( q ) );
                approx_length += w ( q ) * ds ( q );
            }
            LOG_DEBUG ( 2, "Approx. length of quad edge " << f << " = " << approx_length );

            TEST_EQUAL_EPS ( correct_length, approx_length, 1.e-8 );
        }
    }

    double compute_exact_bilinear_integral ( const double coef[], int f )
    {
        int v0 = f;
        int v1 = ( f + 1 ) % 4;
        double xi = quad_vertex_coords[v0][0];
        double xj = quad_vertex_coords[v1][0];
        double yi = quad_vertex_coords[v0][1];
        double yj = quad_vertex_coords[v1][1];
        double dx = xj - xi;
        double dy = yj - yi;
        double ds = std::sqrt ( std::pow ( dx, 2. ) + std::pow ( dy, 2. ) );
        LOG_DEBUG ( 3, "Exact ds = " << ds );
        double b0 = coef[0] + coef[1] * xi + coef[2] * yi + coef[3] * xi * yi;
        double b1 = coef[1] * dx + coef[2] * dy + coef[3] * ( xi * dy + yi * dx );
        double b2 = coef[3] * dx * dy;
        LOG_DEBUG ( 3, "b0 = " << b0 << " b1 = " << b1 << " b2 = " << b2 );
        return ds * ( b0 + 0.5 * b1 + 1. / 3. * b2 );
    }

    void test_bilinear_polynomial ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussLine", 3 );

        // coefficients a00, a10, a01, a11 for polynomial
        // a00 + a10*x + a01*y + a11*x*y
        const double poly_coeff[4] = { 0.4, -1.2, 4.1, 2.5 };

        LOG_DEBUG ( 2, "Integrating " << poly_coeff[0] << " + " << poly_coeff[1] << "x + "
                    << poly_coeff[2] << "y + " << poly_coeff[3] << "xy\n" )

        for ( int f = 0; f < 4; ++f )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), f );

            initialize_for_facet ( elem, quadrature, f );

            const double exact_value = compute_exact_bilinear_integral ( poly_coeff, f );

            LOG_DEBUG ( 2, "Exact value of integral on edge " << f << " = " << exact_value );

            double approx_value = 0.;
            for ( int q = 0; q < num_quadrature_points ( ); ++q )
            {
                LOG_DEBUG ( 3, "w(" << q << ") = " << w ( q ) << ", ds(" << q << ") = " << ds ( q ) );
                LOG_DEBUG ( 3, "x(" << q << ") = " << x ( q ) );
                approx_value +=
                        w ( q ) * ( poly_coeff[0]
                        + poly_coeff[1] * x ( q )[0]
                        + poly_coeff[2] * x ( q )[1]
                        + poly_coeff[3] * x ( q )[0] * x ( q )[1] )
                        * ds ( q );
            }
            LOG_DEBUG ( 2, "Approx. value of integral on facet " << f << " = " << approx_value );

            TEST_EQUAL_EPS ( exact_value, approx_value, 1.e-8 );
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
        quadrature.set_quadrature ( "GaussTriangle", 4 );
    }

    void test_facet_lengths ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussLine", 3 );

        for ( int f = 0; f < 3; ++f )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), f );

            initialize_for_facet ( elem, quadrature, f );

            const int v0 = f;
            const int v1 = ( f + 1 ) % 3;
            LOG_DEBUG ( 3, "v0 = " << v0 << ", v1 = " << v1 );
            LOG_DEBUG ( 3, "x0 = " << tria_vertex_coords[v0][0] << ", x1 = " << tria_vertex_coords[v1][0] );
            LOG_DEBUG ( 3, "y0 = " << tria_vertex_coords[v0][1] << ", y1 = " << tria_vertex_coords[v1][1] );
            const double correct_length =
                    std::sqrt (
                                std::pow ( ( tria_vertex_coords[v1][0] - tria_vertex_coords[v0][0] ), 2. ) +
                                std::pow ( ( tria_vertex_coords[v1][1] - tria_vertex_coords[v0][1] ), 2. ) );

            LOG_DEBUG ( 2, "Exact length of tria edge " << f << " = " << correct_length );

            double approx_length = 0.;
            for ( int q = 0; q < num_quadrature_points ( ); ++q )
            {
                LOG_DEBUG ( 3, "w(" << q << ") = " << w ( q ) << ", ds(" << q << ") = " << ds ( q ) );
                approx_length += w ( q ) * ds ( q );
            }
            LOG_DEBUG ( 2, "Approx. length of tria edge " << f << " = " << approx_length );

            TEST_EQUAL_EPS ( correct_length, approx_length, 1.e-8 );
        }
    }

  private:
    Quadrature<double> quadrature;

};

class HexAsm : private AssemblyAssistant<3, double>
{
  public:

    HexAsm ( )
    {
        quadrature.set_quadrature ( "GaussHexahedron", 27 );
    }

    void test_facet_surface ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussQuadrilateral", 16 );

        for ( int facet = 0; facet < 6; ++facet )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), facet );

            initialize_for_facet ( elem, quadrature, facet );

            const int v0 = elem.get_cell ( ).cell_type ( ).local_vertices_of_entity ( 2, facet )[0];
            const int v1 = elem.get_cell ( ).cell_type ( ).local_vertices_of_entity ( 2, facet )[1];
            const int v2 = elem.get_cell ( ).cell_type ( ).local_vertices_of_entity ( 2, facet )[2];
            const int v3 = elem.get_cell ( ).cell_type ( ).local_vertices_of_entity ( 2, facet )[3];
            LOG_DEBUG ( 3, "v0 = " << v0 << ", v1 = " << v1 << ", v2 = " << v2 << ", v3 = " << v3 );
            LOG_DEBUG ( 3, "x0 = " << hex_vertex_coords[v0][0] << ", x1 = " << hex_vertex_coords[v1][0] << ", x2 = " << hex_vertex_coords[v2][0] << ", x4 = " << hex_vertex_coords[v3][0] );
            LOG_DEBUG ( 3, "y0 = " << hex_vertex_coords[v0][1] << ", y1 = " << hex_vertex_coords[v1][1] << ", y2 = " << hex_vertex_coords[v2][1] << ", y4 = " << hex_vertex_coords[v3][1] );
            LOG_DEBUG ( 3, "z0 = " << hex_vertex_coords[v0][2] << ", z1 = " << hex_vertex_coords[v1][2] << ", z2 = " << hex_vertex_coords[v2][2] << ", z4 = " << hex_vertex_coords[v3][2] );

            //A=\frac{1}{4}\sqrt{4e^2f^2-\left(b^2+d^2-a^2-c^2\right)^2}
            const double a = distance ( 3, hex_vertex_coords[v0], hex_vertex_coords[v1] );
            LOG_DEBUG ( 4, "a == " << a );
            const double b = distance ( 3, hex_vertex_coords[v1], hex_vertex_coords[v2] );
            LOG_DEBUG ( 4, "b == " << b );
            const double c = distance ( 3, hex_vertex_coords[v2], hex_vertex_coords[v3] );
            LOG_DEBUG ( 4, "c == " << c );
            const double d = distance ( 3, hex_vertex_coords[v3], hex_vertex_coords[v0] );
            LOG_DEBUG ( 4, "d == " << d );
            const double e = distance ( 3, hex_vertex_coords[v0], hex_vertex_coords[v2] );
            LOG_DEBUG ( 4, "e == " << e );
            const double f = distance ( 3, hex_vertex_coords[v1], hex_vertex_coords[v3] );
            LOG_DEBUG ( 4, "f == " << f );
            const double diags = 4. * e * e * f * f;
            const double outsides = b * b + d * d - a * a - c * c;
            const double correct_surface = 0.25 * std::sqrt ( diags - outsides * outsides );

            LOG_DEBUG ( 2, "Exact surface of hexa face " << facet << " = " << correct_surface );

            double approx_surface = 0.;
            for ( int q = 0; q < num_quadrature_points ( ); ++q )
            {
                LOG_DEBUG ( 3, "w(" << q << ") = " << w ( q ) << ", ds(" << q << ") = " << ds ( q ) );
                approx_surface += w ( q ) * ds ( q );
            }
            LOG_DEBUG ( 2, "Approx. surface of hexa face " << facet << " = " << approx_surface );

            TEST_EQUAL_EPS ( correct_surface, approx_surface, 1.e-5 );
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
        quadrature.set_quadrature ( "GaussTetrahedron", 5 );
    }

    void test_facet_surface ( const Element<double>& elem )
    {
        Quadrature<double> base_quadrature;
        base_quadrature.set_quadrature ( "GaussTriangle", 4 );

        for ( int facet = 0; facet < 4; ++facet )
        {
            quadrature.set_facet_quadrature ( base_quadrature, elem.get_cell ( ).cell_type ( ).tag ( ), facet );

            initialize_for_facet ( elem, quadrature, facet );

            const int v0 = elem.get_cell ( ).cell_type ( ).local_vertices_of_entity ( 2, facet )[0];
            const int v1 = elem.get_cell ( ).cell_type ( ).local_vertices_of_entity ( 2, facet )[1];
            const int v2 = elem.get_cell ( ).cell_type ( ).local_vertices_of_entity ( 2, facet )[2];

            LOG_DEBUG ( 3, "v0 = " << v0 << ", v1 = " << v1 << ", v2 = " << v2 );
            LOG_DEBUG ( 3, "x0 = " << tet_vertex_coords[v0][0] << ", x1 = " << tet_vertex_coords[v1][0] << ", x2 = " << tet_vertex_coords[v2][0] );
            LOG_DEBUG ( 3, "y0 = " << tet_vertex_coords[v0][1] << ", y1 = " << tet_vertex_coords[v1][1] << ", y2 = " << tet_vertex_coords[v2][1] );
            LOG_DEBUG ( 3, "z0 = " << tet_vertex_coords[v0][2] << ", z1 = " << tet_vertex_coords[v1][2] << ", z2 = " << tet_vertex_coords[v2][2] );

            // From Wikipedia:
            //
            // The shape of the triangle is determined by the lengths
            // of the sides alone. Therefore the area can also be
            // derived from the lengths of the sides. By Heron's
            // formula:
            //
            // T = \sqrt{s(s-a)(s-b)(s-c)}
            //
            // where s = frac{a+b+c}{2} is the semiperimeter, or half
            // of the triangle's perimeter.
            double a = distance ( 3, tet_vertex_coords[v0], tet_vertex_coords[v1] );
            double b = distance ( 3, tet_vertex_coords[v1], tet_vertex_coords[v2] );
            double c = distance ( 3, tet_vertex_coords[v0], tet_vertex_coords[v2] );
            double s = 0.5 * ( a + b + c );
            const double correct_surface = std::sqrt ( s * ( s - a ) * ( s - b ) * ( s - c ) );

            LOG_DEBUG ( 0, "Exact surface of tetra face " << facet << " = " << correct_surface );

            double approx_surface = 0.;
            for ( int q = 0; q < num_quadrature_points ( ); ++q )
            {
                LOG_DEBUG ( 3, "w(" << q << ") = " << w ( q ) << ", ds(" << q << ") = " << ds ( q ) );
                approx_surface += w ( q ) * ds ( q );
            }
            LOG_DEBUG ( 2, "Approx. surface of tetra face " << facet << " = " << approx_surface );

            TEST_EQUAL_EPS ( correct_surface, approx_surface, 1.e-5 );
        }
    }

  private:
    Quadrature<double> quadrature;

};

// helper to compute distance between two points

double distance ( const int dim, const double* point1, const double* point2 )
{
    double dist = 0.;
    for ( int i = 0; i < dim; ++i )
    {
        dist += ( point1[i] - point2[i] ) * ( point1[i] - point2[i] );
    }
    dist = std::sqrt ( dist );
    LOG_DEBUG ( 5, "dist == " << dist );
    assert ( dist >= 0.0 );
    return dist;
}

void quad_facet_test ( )
{
    LOG_DEBUG ( 1, "\n=== Begin Quad Test ====" );

    const int deg = 2;

    MeshPtr quad_mesh = build_quad ( );
    VectorSpace<double> space;
    space.Init ( std::vector<int>( 1, deg ), *quad_mesh );
    Element<double> elem ( space, 0 );

    QuadAsm quad_asm;
    //quad_asm.regression_output(elem);
    quad_asm.test_facet_lengths ( elem );
    // quad_asm.test_bilinear_polynomial(elem);

    LOG_DEBUG ( 1, "\n=== End Quad Test ====" );
}

void tria_facet_test ( )
{
    LOG_DEBUG ( 1, "\n=== Begin Triangle Test ====" );

    const int deg = 2;

    MeshPtr triangle_mesh = build_tria ( );

    VectorSpace<double> space;
    space.Init ( std::vector<int>( 1, deg ), *triangle_mesh );
    Element<double> elem ( space, 0 );

    TriaAsm tria_asm;
    tria_asm.test_facet_lengths ( elem );

    LOG_DEBUG ( 1, "\n=== End Triangle Test ====" );
}

void hexa_facet_test ( )
{
    LOG_DEBUG ( 1, "\n=== Begin Hexahedron Test ====" );

    const int deg = 3;

    MeshPtr hexahedron_mesh = build_hexa ( );

    mesh::VtkWriter writer;
    std::string filename = "hex_mesh.vtu";
    writer.write ( filename.c_str ( ), *hexahedron_mesh );

    VectorSpace<double> space;
    space.Init ( std::vector<int>( 1, deg ), *hexahedron_mesh );
    Element<double> elem ( space, 0 );

    HexAsm hex_asm;
    hex_asm.test_facet_surface ( elem );

    LOG_DEBUG ( 1, "\n=== End Hexahedron Test ====" );
}

void tetra_facet_test ( )
{
    LOG_DEBUG ( 1, "\n=== Begin Tetrahedron Test ====" );

    const int deg = 3;

    MeshPtr tetrahedron_mesh = build_tetra ( );

    mesh::VtkWriter writer;
    std::string filename = "tet_mesh.vtu";
    writer.write ( filename.c_str ( ), *tetrahedron_mesh );

    VectorSpace<double> space;
    space.Init ( std::vector<int>( 1, deg ), *tetrahedron_mesh );
    Element<double> elem ( space, 0 );

    TetAsm tet_asm;
    tet_asm.test_facet_surface ( elem );

    LOG_DEBUG ( 1, "\n=== End Tetrahedron Test ====" );
}
