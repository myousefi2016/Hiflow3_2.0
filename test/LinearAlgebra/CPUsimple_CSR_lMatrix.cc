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
/// @date 2015-10-15

#include "hiflow.h"
#include "linear_algebra/lmp/solvers/cg.h"
#include "linear_algebra/lmp/lmatrix_dense_cpu.h"
#include "linear_algebra/lmp/init_vec_mat.h"
#include "gtest/gtest.h"

using namespace hiflow::la;

/// Test square matrix

TEST ( CPUsimple_CSR_lMatrix, transpose_2x2 )
{
    CPUsimple_CSR_lMatrix<double> A;
    A.Init ( 4, 2, 2, "A" );
    int row[4] = { 0, 0, 1, 1 };
    int col[4] = { 0, 1, 0, 1 };
    A.init_structure ( row, col );
    A.Zeros ( );
    A.add_value ( 0, 0, 1.0 );
    A.add_value ( 0, 1, 2.0 );
    A.add_value ( 1, 0, 3.0 );
    A.add_value ( 1, 1, 4.0 );

    A.transpose_me ( );

    double value;
    A.get_value ( 0, 0, &value );
    EXPECT_DOUBLE_EQ ( 1.0, value );
    A.get_value ( 0, 1, &value );
    EXPECT_DOUBLE_EQ ( 3.0, value );
    A.get_value ( 1, 0, &value );
    EXPECT_DOUBLE_EQ ( 2.0, value );
    A.get_value ( 1, 1, &value );
    EXPECT_DOUBLE_EQ ( 4.0, value );
}

/// Test rectangular matrix

TEST ( CPUsimple_CSR_lMatrix, transpose_1x3 )
{
    CPUsimple_CSR_lMatrix<double> A;
    A.Init ( 3, 1, 3, "A" );
    int row[3] = { 0, 0, 0 };
    int col[3] = { 0, 1, 2 };
    A.init_structure ( row, col );
    A.Zeros ( );
    A.add_value ( 0, 0, 1.0 );
    A.add_value ( 0, 1, 2.0 );
    A.add_value ( 0, 2, 3.0 );

    A.transpose_me ( );

    double value;
    A.get_value ( 0, 0, &value );
    EXPECT_DOUBLE_EQ ( 1.0, value );
    A.get_value ( 1, 0, &value );
    EXPECT_DOUBLE_EQ ( 2.0, value );
    A.get_value ( 2, 0, &value );
    EXPECT_DOUBLE_EQ ( 3.0, value );
}

TEST ( CPUsimple_CSR_lMatrix, add_values )
{
    int rows[] = { 0, 0, 1, 1, 1, 2, 2 };
    int cols[] = { 0, 1, 0, 1, 2, 1, 2 };
    double values[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

    CPUsimple_CSR_lMatrix<double> A;
    A.Init ( 7, 3, 3, "A" );
    A.init_structure ( rows, cols );
    A.Zeros ( );
    int rows_2[] = { 0, 1, 2 };
    int cols_2[] = { 0, 1, 2 };
    A.add_values ( rows_2, sizeof (rows_2 ) / sizeof (rows_2[0] ), cols_2,
                   sizeof (cols_2 ) / sizeof (cols_2[0] ), values );

    EXPECT_DOUBLE_EQ ( 1.0, A.matrix.val[0] );
    EXPECT_DOUBLE_EQ ( 2.0, A.matrix.val[1] );
    EXPECT_DOUBLE_EQ ( 4.0, A.matrix.val[2] );
    EXPECT_DOUBLE_EQ ( 5.0, A.matrix.val[3] );
    EXPECT_DOUBLE_EQ ( 6.0, A.matrix.val[4] );
    EXPECT_DOUBLE_EQ ( 8.0, A.matrix.val[5] );
    EXPECT_DOUBLE_EQ ( 9.0, A.matrix.val[6] );
}

TEST ( CPUsimple_CSR_lMatrix, MatrixMult )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;

    int A_rows[] = { 0, 0 };
    int A_cols[] = { 0, 1 };
    MatrixType A;
    A.Init ( 2, 1, 2, "A" );
    A.init_structure ( A_rows, A_cols );
    A.matrix.val[0] = 1.0;
    A.matrix.val[1] = 2.0;

    int B_rows[] = { 0, 1 };
    int B_cols[] = { 0, 0 };
    MatrixType B;
    B.Init ( 2, 2, 1, "B" );
    B.init_structure ( B_rows, B_cols );
    B.matrix.val[0] = 3.0;
    B.matrix.val[1] = 4.0;

    MatrixType C = *static_cast < MatrixType* > ( A.MatrixMult ( *( static_cast < const hiflow::la::lMatrix<double>* > ( &B ) ) ) );

    EXPECT_DOUBLE_EQ ( 11.0, C.matrix.val[0] );
}

TEST ( CPUsimple_CSR_lMatrix, MatrixMult2 )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;

    int A_rows[] = { 0, 1 };
    int A_cols[] = { 0, 0 };
    MatrixType A;
    A.Init ( 2, 2, 1, "A" );
    A.init_structure ( A_rows, A_cols );
    A.matrix.val[0] = 1.0;
    A.matrix.val[1] = 2.0;

    int B_rows[] = { 0, 0, 0 };
    int B_cols[] = { 0, 1, 2 };
    MatrixType B;
    B.Init ( 3, 1, 3, "B" );
    B.init_structure ( B_rows, B_cols );
    B.matrix.val[0] = 3.0;
    B.matrix.val[1] = 4.0;
    B.matrix.val[2] = 5.0;

    MatrixType C = *static_cast < MatrixType* > ( A.MatrixMult ( *( static_cast < const hiflow::la::lMatrix<double>* > ( &B ) ) ) );

    EXPECT_DOUBLE_EQ ( 3.0, C.matrix.val[0] );
    EXPECT_DOUBLE_EQ ( 4.0, C.matrix.val[1] );
    EXPECT_DOUBLE_EQ ( 5.0, C.matrix.val[2] );
    EXPECT_DOUBLE_EQ ( 6.0, C.matrix.val[3] );
    EXPECT_DOUBLE_EQ ( 8.0, C.matrix.val[4] );
    EXPECT_DOUBLE_EQ ( 10.0, C.matrix.val[5] );
}

TEST ( CPUsimple_CSR_lMatrix, UnorderedInitStructure_square )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;

    int A_rows[] = { 1, 2, 0, 0, 1 };
    int A_cols[] = { 2, 1, 0, 1, 1 };
    MatrixType A;
    A.Init ( 5, 3, 3, "A" );
    A.init_structure ( A_rows, A_cols );

    EXPECT_EQ ( 0, A.matrix.row[0] );
    EXPECT_EQ ( 2, A.matrix.row[1] );
    EXPECT_EQ ( 4, A.matrix.row[2] );
    EXPECT_EQ ( 5, A.matrix.row[3] );

    EXPECT_EQ ( 0, A.matrix.col[0] );
    EXPECT_EQ ( 1, A.matrix.col[1] );
    EXPECT_EQ ( 1, A.matrix.col[2] );
    EXPECT_EQ ( 2, A.matrix.col[3] );
    EXPECT_EQ ( 1, A.matrix.col[4] );
}

TEST ( CPUsimple_CSR_lMatrix, UnorderedInitStructure_rect )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;

    int A_rows[] = { 2, 0, 0, 1 };
    int A_cols[] = { 1, 0, 1, 1 };
    MatrixType A;
    A.Init ( 4, 3, 2, "A" );
    A.init_structure ( A_rows, A_cols );

    EXPECT_EQ ( 0, A.matrix.row[0] );
    EXPECT_EQ ( 2, A.matrix.row[1] );
    EXPECT_EQ ( 3, A.matrix.row[2] );
    EXPECT_EQ ( 4, A.matrix.row[3] );

    EXPECT_EQ ( 0, A.matrix.col[0] );
    EXPECT_EQ ( 1, A.matrix.col[1] );
    EXPECT_EQ ( 1, A.matrix.col[2] );
    EXPECT_EQ ( 1, A.matrix.col[3] );
}
