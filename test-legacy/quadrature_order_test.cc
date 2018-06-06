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

#ifndef HIFLOW_QUADRATURE_ORDER_TEST
#    define HIFLOW_QUADRATURE_ORDER_TEST

#    include "quadrature/quadrature.h"
#    include <cmath>
#    include <iostream>
#    include <string>
#    include <vector>
#    include "test.h"

using namespace hiflow;

const double INTEGRATION_TOL = 1.e-15;

// Gamma Table: sets gamma_table[i] = Gamma(i), (not zero-based).
const int TABLE_SIZE = 100;
std::vector<double> gamma_table;

void init_gamma_table ( )
{
    gamma_table.resize ( TABLE_SIZE );
    gamma_table.at ( 1 ) = 1;

    for ( int i = 2; i < TABLE_SIZE; ++i )
    {
        gamma_table.at ( i ) = double(i - 1 ) * gamma_table.at ( i - 1 );
    }
}

double beta ( int i, int j )
{
    return double(gamma_table.at ( i ) * gamma_table.at ( j ) ) / double(gamma_table.at ( i + j ) );
}

template<class Integrand>
double integrate ( const std::string& name, int order, Integrand F )
{

    Quadrature<double> quad;
    quad.set_quadrature_by_order ( name, order );

    const int num_q = quad.size ( );

    double I = 0;

    for ( int q = 0; q < num_q; ++q )
    {
        I += quad.w ( q ) * F ( quad, q );
    }
    return I;
}

/// 1d integrand: x^p

struct F1d
{

    F1d ( int p ) : p_ ( p )
    {
    }

    double operator() ( const Quadrature<double>& quad, int q ) const
    {
        return std::pow ( quad.x ( q ), p_ );
    }

    int p_;
};

/// 2d integrand: x^p * y^p

struct F2d
{

    F2d ( int p1, int p2 ) : p1_ ( p1 ), p2_ ( p2 )
    {
    }

    double operator() ( const Quadrature<double>& quad, int q ) const
    {
        return std::pow ( quad.x ( q ), p1_ ) * std::pow ( quad.y ( q ), p2_ );
    }

    int p1_, p2_;
};

/// 3d integrand: x^p * y^p * z^p

struct F3d
{

    F3d ( int p ) : p_ ( p )
    {
    }

    double operator() ( const Quadrature<double>& quad, int q ) const
    {
        return std::pow ( quad.x ( q ), p_ ) * std::pow ( quad.y ( q ), p_ ) * std::pow ( quad.z ( q ), p_ );
    }

    int p_;
};

struct F3dTet
{

    F3dTet ( int p ) : p_ ( p )
    {
    }

    double operator() ( const Quadrature<double>& quad, int q ) const
    {
        return std::pow ( quad.z ( q ), p_ );
    }

    int p_;
};

int main ( int argc, char** argv )
{

    init_gamma_table ( );

    // 1d tests -- \int_0^1{x^p dx} = 1/(p+1)
    std::cout << "Testing gauss on line ...\n";
    for ( int p = 0; p < 30; ++p )
    {
        F1d f ( p );
        const double I = integrate ( "GaussLine", p, f );
        const double I_ref = 1 / double(p + 1 );

        std::cout << "p = " << p << "\n";
        std::cout << "I = " << I << "\n";
        std::cout << "I_ref = " << I_ref << "\n\n";

        TEST_EQUAL_EPS ( I, I_ref, INTEGRATION_TOL );
    }

    // 2d tests -- \int_0^1\int_0^1{x^py^p dx dy} = (1/(p+1))^2
    // Tensor-product Gauss rule
    std::cout << "Testing gauss on quad ...\n";
    for ( int p = 0; p < 80; ++p )
    {
        F2d f ( p, p );
        const double I = integrate ( "GaussQuadrilateral", p, f );
        const double I_ref = std::pow ( 1 / double(p + 1 ), 2 );

        std::cout << "p = " << p << "\n";
        std::cout << "I = " << I << "\n";
        std::cout << "I_ref = " << I_ref << "\n\n";

        TEST_EQUAL_EPS ( I, I_ref, INTEGRATION_TOL );
    }

    // Economical Gauss rule -- DO NOT PASS
    std::cout << "Testing economical gauss on quad ...\n";
    for ( int p = 0; p < 21; ++p )
    {
        F2d f ( p, p );
        const double I = integrate ( "EconomicalGaussQuadrilateral", p, f );
        const double I_ref = std::pow ( 1 / double(p + 1 ), 2 );

        std::cout << "p = " << p << "\n";
        std::cout << "I = " << I << "\n";
        std::cout << "I_ref = " << I_ref << "\n\n";

        //        TEST_EQUAL_EPS(I, I_ref, INTEGRATION_TOL);
    }

    // Triangle -- \int_Tri{x^p1 * y^p2}
    std::cout << "Testing triangle ...\n";
    for ( int p1 = 0; p1 < 20; ++p1 )
    {
        for ( int p2 = 0; p2 < 20 - p1; ++p2 )
        {
            F2d f ( p1, p2 );
            const double I = integrate ( "GaussTriangle", p1 + p2, f );
            const double I_ref = beta ( p2 + 2, p1 + 1 ) / double(p2 + 1 );

            std::cout << "p1 = " << p1 << ", p2 = " << p2 << "\n";
            std::cout << "I = " << I << "\n";
            std::cout << "I_ref = " << I_ref << "\n\n";

            TEST_EQUAL_EPS ( I, I_ref, INTEGRATION_TOL );
        }
    }
    for ( int p2 = 0; p2 < 20; ++p2 )
    {
        for ( int p1 = 0; p1 < 20 - p2; ++p1 )
        {
            F2d f ( p1, p2 );
            const double I = integrate ( "GaussTriangle", p1 + p2, f );
            const double I_ref = beta ( p2 + 2, p1 + 1 ) / double(p2 + 1 );

            std::cout << "p1 = " << p1 << ", p2 = " << p2 << "\n";
            std::cout << "I = " << I << "\n";
            std::cout << "I_ref = " << I_ref << "\n\n";

            TEST_EQUAL_EPS ( I, I_ref, INTEGRATION_TOL );
        }
    }

    // Tetrahedron -- \int_Tet{z^p}
    for ( int p = 0; p < 21; ++p )
    {
        F3dTet f ( p );
        const double I = integrate ( "GaussTetrahedron", p, f );
        const double I_ref = 1. / double(p * p * p + 6. * p * p + 11. * p + 6. );

        std::cout << "p = " << p << "\n";
        std::cout << "I = " << I << "\n";
        std::cout << "I_ref = " << I_ref << "\n\n";

        TEST_EQUAL_EPS ( I, I_ref, INTEGRATION_TOL );
    }

    // Line -- \int_line{x^p}
    std::cout << "Testing line ...\n";
    for ( int p = 0; p < 20; ++p )
    {

        F1d f ( p );
        const double I = integrate ( "GaussLine", p, f );
        const double I_ref = 1.0 / ( p + 1.0 );

        std::cout << "p = " << p << "\n";
        std::cout << "I = " << I << "\n";
        std::cout << "I_ref = " << I_ref << "\n\n";

        TEST_EQUAL_EPS ( I, I_ref, INTEGRATION_TOL );

    }

    return 0;
}

#endif
