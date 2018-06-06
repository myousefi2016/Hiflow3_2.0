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

/// @author Simon Gawlok

#include "hiflow.h"
#include "gtest/gtest.h"
#include "mpi.h"

using namespace hiflow::la;

/// Test adding and getting values to/from HypreMatrix

TEST ( HypreMatrix, add_get_value )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate dimensions
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set up input and output vectors
    std::vector<int> global_offsets ( num_procs + 1 );
    for ( int i = 0; i < num_procs + 1; ++i )
    {
        global_offsets[i] = i * local_dim;
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1, 0 );

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate sparsity structure
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    // first row
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower );
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower + 1 );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        rows_diag.push_back ( i );
        cols_diag.push_back ( i - 1 );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i + 1 );
    }
    // last row
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper - 1 );
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper );

    // Create and initialize HypreMatrix
    HypreMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), rows_diag.size ( ), vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), rows_offdiag.size ( ) );

    // Add values to matrix
    // first row
    matrix.Add ( local_ilower, local_ilower, 2. );
    matrix.Add ( local_ilower, local_ilower + 1, -1. );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        matrix.Add ( i, i - 1, -1. );
        matrix.Add ( i, i, 2. );
        matrix.Add ( i, i + 1, -1. );
    }
    // last row
    matrix.Add ( local_iupper, local_iupper - 1, -1. );
    matrix.Add ( local_iupper, local_iupper, 2. );

    // Get values in matrix
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> val_res ( local_dim*local_dim, 0. );
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        rows.push_back ( i );
        cols.push_back ( i );
    }

    matrix.GetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val_res ) );

    // Check values in val_res
    // first row
    EXPECT_DOUBLE_EQ ( 2., val_res[0] );
    EXPECT_DOUBLE_EQ ( -1., val_res[1] );
    for ( int j = 2; j <= local_iupper - local_ilower; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[j] );
    }

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // first cols
        for ( int j = 0; j <= i - 2; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
        // tridiagonal PARTICULAR
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i - 1] );
        EXPECT_DOUBLE_EQ ( 2., val_res[i * local_dim + i] );
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i + 1] );
        // last cols
        for ( int j = i + 2; j <= local_iupper - local_ilower - 1; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
    }

    // last row
    for ( int j = 0; j <= local_iupper - local_ilower - 2; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + j] );
    }
    EXPECT_DOUBLE_EQ ( -1., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] );
    EXPECT_DOUBLE_EQ ( 2., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] );
}

/// Test setting values in HypreMatrix

TEST ( HypreMatrix, set_value )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate dimensions
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set up input and output vectors
    std::vector<int> global_offsets ( num_procs + 1 );
    for ( int i = 0; i < num_procs + 1; ++i )
    {
        global_offsets[i] = i * local_dim;
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1, 0 );

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate sparsity structure
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    // first row
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower );
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower + 1 );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        rows_diag.push_back ( i );
        cols_diag.push_back ( i - 1 );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i + 1 );
    }
    // last row
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper - 1 );
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper );

    // Create and initialize HypreMatrix
    HypreMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), rows_diag.size ( ), vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), rows_offdiag.size ( ) );

    // Add values to matrix
    // first row
    matrix.SetValue ( local_ilower, local_ilower, 2. );
    matrix.SetValue ( local_ilower, local_ilower + 1, -1. );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        matrix.SetValue ( i, i - 1, -1. );
        matrix.SetValue ( i, i, 2. );
        matrix.SetValue ( i, i + 1, -1. );
    }
    // last row
    matrix.SetValue ( local_iupper, local_iupper - 1, -1. );
    matrix.SetValue ( local_iupper, local_iupper, 2. );

    // Get values in matrix
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> val_res ( local_dim*local_dim, 0. );
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        rows.push_back ( i );
        cols.push_back ( i );
    }

    matrix.GetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val_res ) );

    // Check values in val_res
    // first row
    EXPECT_DOUBLE_EQ ( 2., val_res[0] );
    EXPECT_DOUBLE_EQ ( -1., val_res[1] );
    for ( int j = 2; j <= local_iupper - local_ilower; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[j] );
    }

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // first cols
        for ( int j = 0; j <= i - 2; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
        // tridiagonal PARTICULAR
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i - 1] );
        EXPECT_DOUBLE_EQ ( 2., val_res[i * local_dim + i] );
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i + 1] );
        // last cols
        for ( int j = i + 2; j <= local_iupper - local_ilower - 1; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
    }

    // last row
    for ( int j = 0; j <= local_iupper - local_ilower - 2; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + j] );
    }
    EXPECT_DOUBLE_EQ ( -1., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] );
    EXPECT_DOUBLE_EQ ( 2., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] );
}

/// Test adding values (vectorized) to HypreMatrix

TEST ( HypreMatrix, add_values )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate dimensions
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set up input and output vectors
    std::vector<int> global_offsets ( num_procs + 1 );
    for ( int i = 0; i < num_procs + 1; ++i )
    {
        global_offsets[i] = i * local_dim;
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1, 0 );

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate sparsity structure
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    // first row
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower );
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower + 1 );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        rows_diag.push_back ( i );
        cols_diag.push_back ( i - 1 );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i + 1 );
    }
    // last row
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper - 1 );
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper );

    // Create and initialize HypreMatrix
    HypreMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), rows_diag.size ( ), vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), rows_offdiag.size ( ) );

    // Add values to matrix
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> val ( local_dim*local_dim, 0. );
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        rows.push_back ( i );
        cols.push_back ( i );
    }
    // first row
    val[0] = 2.;
    val[1] = -1.;

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // tridiagonal PARTICULAR
        val[i * local_dim + i - 1] = -1.;
        val[i * local_dim + i] = 2.;
        val[i * local_dim + i + 1] = -1.;
        ;
    }

    // last row
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] = -1.;
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] = 2.;

    matrix.Add ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val ) );

    // Get values in matrix
    std::vector<double> val_res ( local_dim*local_dim, 0. );

    matrix.GetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val_res ) );

    // Check values in val_res
    // first row
    EXPECT_DOUBLE_EQ ( 2., val_res[0] );
    EXPECT_DOUBLE_EQ ( -1., val_res[1] );
    for ( int j = 2; j <= local_iupper - local_ilower; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[j] );
    }

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // first cols
        for ( int j = 0; j <= i - 2; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
        // tridiagonal PARTICULAR
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i - 1] );
        EXPECT_DOUBLE_EQ ( 2., val_res[i * local_dim + i] );
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i + 1] );
        // last cols
        for ( int j = i + 2; j <= local_iupper - local_ilower - 1; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
    }

    // last row
    for ( int j = 0; j <= local_iupper - local_ilower - 2; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + j] );
    }
    EXPECT_DOUBLE_EQ ( -1., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] );
    EXPECT_DOUBLE_EQ ( 2., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] );
}

/// Test setting values (vectorized) to HypreMatrix

TEST ( HypreMatrix, set_values )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate dimensions
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set up input and output vectors
    std::vector<int> global_offsets ( num_procs + 1 );
    for ( int i = 0; i < num_procs + 1; ++i )
    {
        global_offsets[i] = i * local_dim;
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1, 0 );

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate sparsity structure
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    // first row
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower );
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower + 1 );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        rows_diag.push_back ( i );
        cols_diag.push_back ( i - 1 );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i + 1 );
    }
    // last row
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper - 1 );
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper );

    // Create and initialize HypreMatrix
    HypreMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), rows_diag.size ( ), vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), rows_offdiag.size ( ) );

    // Add values to matrix
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> val ( local_dim*local_dim, 0. );
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        rows.push_back ( i );
        cols.push_back ( i );
    }
    // first row
    val[0] = 2.;
    val[1] = -1.;

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // tridiagonal PARTICULAR
        val[i * local_dim + i - 1] = -1.;
        val[i * local_dim + i] = 2.;
        val[i * local_dim + i + 1] = -1.;
        ;
    }

    // last row
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] = -1.;
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] = 2.;

    matrix.SetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val ) );

    // Get values in matrix
    std::vector<double> val_res ( local_dim*local_dim, 0. );

    matrix.GetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val_res ) );

    // Check values in val_res
    // first row
    EXPECT_DOUBLE_EQ ( 2., val_res[0] );
    EXPECT_DOUBLE_EQ ( -1., val_res[1] );
    for ( int j = 2; j <= local_iupper - local_ilower; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[j] );
    }

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // first cols
        for ( int j = 0; j <= i - 2; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
        // tridiagonal PARTICULAR
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i - 1] );
        EXPECT_DOUBLE_EQ ( 2., val_res[i * local_dim + i] );
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i + 1] );
        // last cols
        for ( int j = i + 2; j <= local_iupper - local_ilower - 1; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
    }

    // last row
    for ( int j = 0; j <= local_iupper - local_ilower - 2; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + j] );
    }
    EXPECT_DOUBLE_EQ ( -1., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] );
    EXPECT_DOUBLE_EQ ( 2., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] );
}

/// Test scaling of HypreMatrix

TEST ( HypreMatrix, scale )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate dimensions
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set up input and output vectors
    std::vector<int> global_offsets ( num_procs + 1 );
    for ( int i = 0; i < num_procs + 1; ++i )
    {
        global_offsets[i] = i * local_dim;
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1, 0 );

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate sparsity structure
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    // first row
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower );
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower + 1 );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        rows_diag.push_back ( i );
        cols_diag.push_back ( i - 1 );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i + 1 );
    }
    // last row
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper - 1 );
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper );

    // Create and initialize HypreMatrix
    HypreMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), rows_diag.size ( ), vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), rows_offdiag.size ( ) );

    // Add values to matrix
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> val ( local_dim*local_dim, 0. );
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        rows.push_back ( i );
        cols.push_back ( i );
    }
    // first row
    val[0] = 2.;
    val[1] = -1.;

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // tridiagonal PARTICULAR
        val[i * local_dim + i - 1] = -1.;
        val[i * local_dim + i] = 2.;
        val[i * local_dim + i + 1] = -1.;
        ;
    }

    // last row
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] = -1.;
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] = 2.;

    matrix.SetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val ) );

    // Scale matrix
    const double factor = 0.5;
    matrix.Scale ( factor );

    // Get values in matrix
    std::vector<double> val_res ( local_dim*local_dim, 0. );

    matrix.GetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val_res ) );

    // Check values in val_res
    // first row
    EXPECT_DOUBLE_EQ ( 2. * factor, val_res[0] );
    EXPECT_DOUBLE_EQ ( -1. * factor, val_res[1] );
    for ( int j = 2; j <= local_iupper - local_ilower; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0. * factor, val_res[j] );
    }

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // first cols
        for ( int j = 0; j <= i - 2; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0. * factor, val_res[i * local_dim + j] );
        }
        // tridiagonal PARTICULAR
        EXPECT_DOUBLE_EQ ( -1. * factor, val_res[i * local_dim + i - 1] );
        EXPECT_DOUBLE_EQ ( 2. * factor, val_res[i * local_dim + i] );
        EXPECT_DOUBLE_EQ ( -1. * factor, val_res[i * local_dim + i + 1] );
        // last cols
        for ( int j = i + 2; j <= local_iupper - local_ilower - 1; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0. * factor, val_res[i * local_dim + j] );
        }
    }

    // last row
    for ( int j = 0; j <= local_iupper - local_ilower - 2; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0. * factor, val_res[( local_iupper - local_ilower ) * local_dim + j] );
    }
    EXPECT_DOUBLE_EQ ( -1. * factor, val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] );
    EXPECT_DOUBLE_EQ ( 2. * factor, val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] );
}

/// Test zeroing a HypreMatrix

TEST ( HypreMatrix, zeros )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate dimensions
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set up input and output vectors
    std::vector<int> global_offsets ( num_procs + 1 );
    for ( int i = 0; i < num_procs + 1; ++i )
    {
        global_offsets[i] = i * local_dim;
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1, 0 );

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate sparsity structure
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    // first row
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower );
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower + 1 );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        rows_diag.push_back ( i );
        cols_diag.push_back ( i - 1 );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i + 1 );
    }
    // last row
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper - 1 );
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper );

    // Create and initialize HypreMatrix
    HypreMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), rows_diag.size ( ), vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), rows_offdiag.size ( ) );

    // Add values to matrix
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> val ( local_dim*local_dim, 0. );
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        rows.push_back ( i );
        cols.push_back ( i );
    }
    // first row
    val[0] = 2.;
    val[1] = -1.;

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // tridiagonal PARTICULAR
        val[i * local_dim + i - 1] = -1.;
        val[i * local_dim + i] = 2.;
        val[i * local_dim + i + 1] = -1.;
        ;
    }

    // last row
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] = -1.;
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] = 2.;

    matrix.SetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val ) );

    // Set matrix to zero
    matrix.Zeros ( );

    // Get values in matrix
    std::vector<double> val_res ( local_dim*local_dim, 0. );

    matrix.GetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val_res ) );

    // Check values in val_res
    // first row
    EXPECT_DOUBLE_EQ ( 0., val_res[0] );
    EXPECT_DOUBLE_EQ ( 0., val_res[1] );
    for ( int j = 2; j <= local_iupper - local_ilower; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[j] );
    }

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // first cols
        for ( int j = 0; j <= i - 2; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
        // tridiagonal PARTICULAR
        EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + i - 1] );
        EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + i] );
        EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + i + 1] );
        // last cols
        for ( int j = i + 2; j <= local_iupper - local_ilower - 1; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
    }

    // last row
    for ( int j = 0; j <= local_iupper - local_ilower - 2; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + j] );
    }
    EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] );
    EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] );
}

/// Test Matrix-Vector Multiplication

TEST ( HypreMatrix, vector_mult )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate dimensions
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set up input and output vectors
    std::vector<int> global_offsets ( num_procs + 1 );
    for ( int i = 0; i < num_procs + 1; ++i )
    {
        global_offsets[i] = i * local_dim;
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1, 0 );

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate sparsity structure
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    // first row
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower );
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower + 1 );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        rows_diag.push_back ( i );
        cols_diag.push_back ( i - 1 );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i + 1 );
    }
    // last row
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper - 1 );
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper );

    // Create and initialize HypreMatrix
    HypreMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), rows_diag.size ( ), vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), rows_offdiag.size ( ) );

    // Add values to matrix
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> val ( local_dim*local_dim, 0. );
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        rows.push_back ( i );
        cols.push_back ( i );
    }
    // first row
    val[0] = 2.;
    val[1] = -1.;

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // tridiagonal PARTICULAR
        val[i * local_dim + i - 1] = -1.;
        val[i * local_dim + i] = 2.;
        val[i * local_dim + i + 1] = -1.;
        ;
    }

    // last row
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] = -1.;
    val[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] = 2.;

    matrix.SetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val ) );

    HypreVector<double> in, out;
    in.Init ( comm, laCouplings );
    out.Init ( comm, laCouplings );

    // Set different vector components
    std::vector<int> ind;
    std::vector<double> vals;
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        ind.push_back ( i );
        vals.push_back ( 1. );
    }

    in.Add ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( vals ) );

    matrix.VectorMult ( in, &out );

    // Check for correct result
    EXPECT_DOUBLE_EQ ( 1., out.GetValue ( local_ilower ) );
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        EXPECT_DOUBLE_EQ ( 0., out.GetValue ( i ) );
    }
    EXPECT_DOUBLE_EQ ( 1., out.GetValue ( local_iupper ) );
}

/// Test making some rows of matrix diagonal

TEST ( HypreMatrix, diagonalize_rows )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate dimensions
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set up input and output vectors
    std::vector<int> global_offsets ( num_procs + 1 );
    for ( int i = 0; i < num_procs + 1; ++i )
    {
        global_offsets[i] = i * local_dim;
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1, 0 );

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate sparsity structure
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    // first row
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower );
    rows_diag.push_back ( local_ilower );
    cols_diag.push_back ( local_ilower + 1 );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        rows_diag.push_back ( i );
        cols_diag.push_back ( i - 1 );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i );
        rows_diag.push_back ( i );
        cols_diag.push_back ( i + 1 );
    }
    // last row
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper - 1 );
    rows_diag.push_back ( local_iupper );
    cols_diag.push_back ( local_iupper );

    // Create and initialize HypreMatrix
    HypreMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), rows_diag.size ( ), vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), rows_offdiag.size ( ) );

    // Add values to matrix
    // first row
    matrix.Add ( local_ilower, local_ilower, 2. );
    matrix.Add ( local_ilower, local_ilower + 1, -1. );

    // rows except first and last one
    for ( int i = local_ilower + 1; i <= local_iupper - 1; ++i )
    {
        matrix.Add ( i, i - 1, -1. );
        matrix.Add ( i, i, 2. );
        matrix.Add ( i, i + 1, -1. );
    }
    // last row
    matrix.Add ( local_iupper, local_iupper - 1, -1. );
    matrix.Add ( local_iupper, local_iupper, 2. );

    // Get values in matrix
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> val_res ( local_dim*local_dim, 0. );
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        rows.push_back ( i );
        cols.push_back ( i );
    }

    // diagonalize first and last local row
    std::vector<int> diag_indices;
    diag_indices.push_back ( local_ilower );
    diag_indices.push_back ( local_iupper );

    double diag_val = 10.;
    matrix.diagonalize_rows ( vec2ptr ( diag_indices ), diag_indices.size ( ), diag_val );

    matrix.GetValues ( vec2ptr ( rows ), rows.size ( ), vec2ptr ( cols ), cols.size ( ), vec2ptr ( val_res ) );

    // Check values in val_res
    // first row
    EXPECT_DOUBLE_EQ ( diag_val, val_res[0] );
    EXPECT_DOUBLE_EQ ( 0., val_res[1] );
    for ( int j = 2; j <= local_iupper - local_ilower; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[j] );
    }

    // rows except first and last one
    for ( int i = 1; i <= local_iupper - local_ilower - 1; ++i )
    {
        // first cols
        for ( int j = 0; j <= i - 2; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
        // tridiagonal PARTICULAR
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i - 1] );
        EXPECT_DOUBLE_EQ ( 2., val_res[i * local_dim + i] );
        EXPECT_DOUBLE_EQ ( -1., val_res[i * local_dim + i + 1] );
        // last cols
        for ( int j = i + 2; j <= local_iupper - local_ilower - 1; ++j )
        {
            EXPECT_DOUBLE_EQ ( 0., val_res[i * local_dim + j] );
        }
    }

    // last row
    for ( int j = 0; j <= local_iupper - local_ilower - 2; ++j )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + j] );
    }
    EXPECT_DOUBLE_EQ ( 0., val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower - 1 )] );
    EXPECT_DOUBLE_EQ ( diag_val, val_res[( local_iupper - local_ilower ) * local_dim + ( local_iupper - local_ilower )] );
}
