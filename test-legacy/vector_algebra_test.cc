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

/// \brief Unit test of classes and functions in vector_algebra.h
/// \TODO Implement automatic testing here.

#include <iostream>

#include "common/vector_algebra.h"

using namespace hiflow;

int main ( int argc, char *argv[] )
{
    double b_val[3] = { 1.5, -2., 4.8 };
    Vec<3, double> b ( &b_val[0] );

    double A_val[9] = { 4.8, -3.2, 1.25, -0.6, 1.5, 4.6, -2.7, -1.7, 6.1 };
    Mat<3, 3, double> A ( &A_val[0] );

    Mat<3, 3, double> Ainv;
    inv ( A, Ainv );

    Vec<3, double> c = A * b;

    for ( int i = 0; i < 3; ++i )
    {
        std::cout << "c(" << i << ") = " << c[i] << std::endl;
    }

    std::cout << "det(A) = " << det ( A ) << std::endl;

    for ( int i = 0; i < 3; ++i )
    {
        for ( int j = 0; j < 3; ++j )
        {
            std::cout << "Ainv(" << i << "," << j << ") = " << Ainv ( i, j ) << "\n";
        }
    }

    double x1_val[3] = { 1., 0., 0. };
    double x2_val[3] = { 0., 1., 0. };
    double x3_val[3] = { 0., 0., 1. };
    Vec<3, double> x1 ( &x1_val[0] );
    Vec<3, double> x2 ( &x2_val[0] );
    Vec<3, double> x3 ( &x3_val[0] );
    std::cout << "orientation = " << dot ( x3, cross ( x1, x2 ) ) << "\n";

    return 0;
}
