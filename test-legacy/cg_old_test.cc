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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-06-30

#include "hiflow.h"
#include "linear_algebra/lmp/solvers/cg.h"
#include "linear_algebra/lmp/lmatrix_dense_cpu.h"
#include <unittestpp.h>
#include <iostream>

using namespace hiflow::la;

class NoPreconditioner
{
};

template <class MatrixType, class VectorType, class PreconType>
struct TestFunctor
{

    void operator() ( ) const
    {
        VectorType b ( 3, "b" );
        b.add_value ( 0, 1.0 );
        b.add_value ( 1, -2.0 );
        b.add_value ( 2, 0.0 );

        VectorType x ( 3, "x" );
        x.Zeros ( );

        //CPUsimple_DENSE_lMatrix<double> A(9, 3, 3, "A"); // If dense matrix is used, the preconditioner build step is crashing.
        MatrixType A ( 9, 3, 3, "A" );
        A.add_value ( 0, 0, 2.0 );
        A.add_value ( 0, 1, -1.0 );
        A.add_value ( 0, 2, 0.0 );
        A.add_value ( 1, 0, -1.0 );
        A.add_value ( 1, 1, 2.0 );
        A.add_value ( 1, 2, -1.0 );
        A.add_value ( 2, 0, 0.0 );
        A.add_value ( 2, 1, -1.0 );
        A.add_value ( 2, 2, 2.0 );

        PreconType precond;
        precond.SetupOperator ( A );
        precond.Init ( x );
        precond.Build ( );

        double rel_eps = 0.0;
        double abs_eps = 0.0;
        int max_iter = 100;
        int print_level = 0;

        int return_value = cg ( &x, &b, &A, rel_eps, abs_eps, max_iter, print_level, &precond );

        CHECK_EQUAL ( 0, return_value );
    }
};

template <class MatrixType, class VectorType>
struct TestFunctor<MatrixType, VectorType, NoPreconditioner>
{

    void operator() ( ) const
    {
        VectorType b ( 3, "b" );
        b.add_value ( 0, 1.0 );
        b.add_value ( 1, -2.0 );
        b.add_value ( 2, 0.0 );

        VectorType x ( 3, "x" );
        x.Zeros ( );

        //CPUsimple_DENSE_lMatrix<double> A(9, 3, 3, "A"); // If dense matrix is used, the preconditioner build step is crashing.
        MatrixType A ( 9, 3, 3, "A" );
        A.add_value ( 0, 0, 2.0 );
        A.add_value ( 0, 1, -1.0 );
        A.add_value ( 0, 2, 0.0 );
        A.add_value ( 1, 0, -1.0 );
        A.add_value ( 1, 1, 2.0 );
        A.add_value ( 1, 2, -1.0 );
        A.add_value ( 2, 0, 0.0 );
        A.add_value ( 2, 1, -1.0 );
        A.add_value ( 2, 2, 2.0 );

        VectorType diagA ( 3, "diagA" );
        diagA.add_value ( 0, 1.0 / 3.0 );
        diagA.add_value ( 1, -1.0 / 2.0 );
        diagA.add_value ( 2, -1.0 );

        double rel_eps = 0.0;
        double abs_eps = 0.0;
        int max_iter = 100;
        int print_level = 1;

        int return_value = cg ( &x, &b, &A, rel_eps, abs_eps, max_iter, print_level );
        CHECK_EQUAL ( 0, return_value );

        // Check result
        int x_size = 3;
        int x_index[3] = { 0, 1, 2 };
        double x_values[3];
        x.GetValues ( &x_index[0], x_size, &x_values[0] );

        //std::cout << x_values[0] << " " << x_values[1] << " " << x_values[2] << std::endl;
        //std::cout << x.buffer[0] << " " << x.buffer[1] << " " << x.buffer[2] << std::endl;

        VectorType res ( 3, "result" );
        res.add_value ( 0, -0.25 );
        res.add_value ( 1, -1.5 );
        res.add_value ( 2, -0.75 );

        //CHECK_EQUAL(res, b);
        CHECK_CLOSE ( res.buffer[0], x.buffer[0], 1e-12 );
        CHECK_CLOSE ( res.buffer[1], x.buffer[1], 1e-12 );
        CHECK_CLOSE ( res.buffer[2], x.buffer[2], 1e-12 );
    }
};

TEST ( cg_simple_dense )
{
    TestFunctor<CPUsimple_DENSE_lMatrix<double>, CPUsimple_lVector<double>, NoPreconditioner >( )( );
}

TEST ( cg_simple_csr )
{
    TestFunctor<CPUsimple_CSR_lMatrix<double>, CPUsimple_lVector<double>, NoPreconditioner >( )( );
}

TEST ( cg_simple_csr_jacobi )
{
    TestFunctor<CPUsimple_CSR_lMatrix<double>, CPUsimple_lVector<double>, lPreconditioner_Jacobi<double> >( )( );
}

#if 0

TEST ( cg_openmp_csr_jacobi )
{
    TestFunctor<CPUopenmp_CSR_lMatrix<double>, CPUopenmp_lVector<double>, lPreconditioner_Jacobi<double> >( )( );
}
#endif

int main ( int argc, char** argv )
{
    return UnitTest::RunAllTests ( );
}
