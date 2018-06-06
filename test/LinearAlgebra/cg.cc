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
#include "linear_algebra/lmp/init_vec_mat.h"
#include "gtest/gtest.h"

using namespace hiflow::la;

struct TestFunctor
{

    void operator() ( const enum PLATFORM &platform,
            const enum IMPLEMENTATION &vector_implementation,
            const enum IMPLEMENTATION &matrix_implementation,
            const enum MATRIX_FORMAT &matrix_format,
            const enum MATRIX_FREE_PRECOND &precond ) const
    {
        lVector<double> *b =
                init_vector<double>( 3, "b", platform, vector_implementation );
        b->add_value ( 0, 1.0 );
        b->add_value ( 1, -2.0 );
        b->add_value ( 2, 0.0 );

        lVector<double> *x =
                init_vector<double>( 3, "x", platform, vector_implementation );
        x->Zeros ( );

        lMatrix<double> *A = init_matrix<double>(
                9, 3, 3, "A", platform, matrix_implementation, matrix_format );

        if ( matrix_format == CSR )
        {
            int num_rows = 3;
            int num_cols = 3;
            int nnz = num_rows * num_cols;

            int *rows = new int[nnz];
            for ( int i = 0; i < num_rows; ++i )
                for ( int j = 0; j < num_rows; ++j ) rows[i * num_rows + j] = i;

            int *cols = new int[nnz];
            for ( int i = 0; i < num_rows; ++i )
                for ( int j = 0; j < num_cols; ++j ) cols[i * num_rows + j] = j;

            A->init_structure ( rows, cols );
        }

        A->add_value ( 0, 0, 2.0 );
        A->add_value ( 0, 1, -1.0 );
        A->add_value ( 0, 2, 0.0 );
        A->add_value ( 1, 0, -1.0 );
        A->add_value ( 1, 1, 2.0 );
        A->add_value ( 1, 2, -1.0 );
        A->add_value ( 2, 0, 0.0 );
        A->add_value ( 2, 1, -1.0 );
        A->add_value ( 2, 2, 2.0 );

        double rel_eps = 0.0;
        double abs_eps = 0.0;
        int max_iter = 100;
        int print_level = -1;

        int return_value = cg ( x, b, A, rel_eps, abs_eps, max_iter, print_level );
        EXPECT_EQ ( 0, return_value );

        // Check result
        int x_size = 3;
        int x_index[3] = { 0, 1, 2 };
        double x_values[3];
        x->GetValues ( &x_index[0], x_size, &x_values[0] );

        EXPECT_NEAR ( -0.25, x_values[0], 1e-12 );
        EXPECT_NEAR ( -1.50, x_values[1], 1e-12 );
        EXPECT_NEAR ( -0.75, x_values[2], 1e-12 );
    }
};

TEST ( CG, simple_dense )
{
    TestFunctor ( )( CPU, NAIVE, NAIVE, DENSE, NOPRECOND );
}

TEST ( CG, simple_csr )
{
    TestFunctor ( )( CPU, NAIVE, NAIVE, CSR, NOPRECOND );
}

#ifdef WITH_OPENMP

TEST ( CG, openmp_csr )
{
    TestFunctor ( )( CPU, OPENMP, OPENMP, CSR, NOPRECOND );
}
#endif

TEST ( CG, simple_csr_jacobi )
{
    TestFunctor ( )( CPU, NAIVE, NAIVE, CSR, JACOBI );
}

#ifdef WITH_OPENMP

TEST ( CG, openmp_csr_jacobi )
{
    TestFunctor ( )( CPU, OPENMP, OPENMP, CSR, JACOBI );
}
#endif

#if 0

TEST ( CG, scalartex_csr )
{
    TestFunctor ( )( GPU, BLAS, SCALAR_TEX, CSR, NOPRECOND );
}
#endif
