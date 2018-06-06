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

/// Test setting value in HpreVector

TEST ( HypreVector, set_get_value )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec.SetValue ( i, static_cast < double > ( i ) );
    }

    // Get vector components and check for correct values
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( static_cast < double > ( i ), vec.GetValue ( i ) );
    }
}

/// Test setting values (vectorized) in HpreVector

TEST ( HypreVector, set_get_values )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    std::vector<int> ind;
    std::vector<double> val;
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        ind.push_back ( i );
        val.push_back ( static_cast < double > ( i ) );
    }

    vec.SetValues ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( val ) );

    // Get vector components and check for correct values
    std::vector<double> val_res ( val.size ( ) );
    vec.GetValues ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( val_res ) );

    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( static_cast < double > ( i ), val_res[i - local_ilower] );
    }
}

/// Test adding value in HpreVector

TEST ( HypreVector, add_value )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec.Add ( i, static_cast < double > ( i ) );
    }

    // Get vector components and check for correct values
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( static_cast < double > ( i ), vec.GetValue ( i ) );
    }
}

/// Test adding values (vectorized) in HpreVector

TEST ( HypreVector, add_values )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    std::vector<int> ind;
    std::vector<double> val;
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        ind.push_back ( i );
        val.push_back ( static_cast < double > ( i ) );
    }

    vec.Add ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( val ) );

    // Get vector components and check for correct values
    std::vector<double> val_res ( val.size ( ) );
    vec.GetValues ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( val_res ) );

    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( static_cast < double > ( i ), val_res[i - local_ilower] );
    }
}

/// Test cloning HypreVector without content

TEST ( HypreVector, clone_without_content )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // Clone vector without content
    HypreVector<double> vec2;
    vec2.CloneFromWithoutContent ( vec );

    // Check for correct dimensions
    EXPECT_EQ ( vec.size_local ( ), vec2.size_local ( ) );
    EXPECT_EQ ( vec.size_global ( ), vec2.size_global ( ) );

}

/// Test cloning complete HypreVector

TEST ( HypreVector, clone )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec.Add ( i, static_cast < double > ( i ) );
    }

    // Clone to different HypreVector
    HypreVector<double> vec2;
    vec2.CloneFrom ( vec );

    // Check for correct dimensions
    EXPECT_EQ ( vec.size_local ( ), vec2.size_local ( ) );
    EXPECT_EQ ( vec.size_global ( ), vec2.size_global ( ) );

    // Get vector components and check for correct values
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( static_cast < double > ( i ), vec2.GetValue ( i ) );
    }
}

// Test dot product

TEST ( HypreVector, dot )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec.Add ( i, static_cast < double > ( i ) );
    }

    // Clone structure to different HypreVector
    HypreVector<double> vec2;
    vec2.CloneFromWithoutContent ( vec );

    // Check for correct dimensions
    EXPECT_EQ ( vec.size_local ( ), vec2.size_local ( ) );
    EXPECT_EQ ( vec.size_global ( ), vec2.size_global ( ) );

    // Set components of second vector
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec2.Add ( i, static_cast < double > ( i + 2 ) );
    }

    // Compute dot product
    double dot_res = vec.Dot ( vec2 );

    // Maximum component index
    const double N = static_cast < double > ( global_dim - 1 );
    const double res_expected = ( 1. / 6. ) * N * ( N + 1. ) * ( 2. * N + 1. ) + N * ( N + 1. );

    EXPECT_DOUBLE_EQ ( res_expected, dot_res );
}

// Test axpy

TEST ( HypreVector, axpy )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec.Add ( i, static_cast < double > ( i ) );
    }

    // Clone structure to different HypreVector
    HypreVector<double> vec2;
    vec2.CloneFromWithoutContent ( vec );

    // Check for correct dimensions
    EXPECT_EQ ( vec.size_local ( ), vec2.size_local ( ) );
    EXPECT_EQ ( vec.size_global ( ), vec2.size_global ( ) );

    // Set components of second vector
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec2.Add ( i, static_cast < double > ( i + 2 ) );
    }

    const double factor = 4.0;

    vec.Axpy ( vec2, factor );

    // Get vector components and check for correct values
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( static_cast < double > ( i + factor * ( i + 2 ) ), vec.GetValue ( i ) );
    }

}

/// Test scaling a vector

TEST ( HypreVector, scale )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec.Add ( i, static_cast < double > ( i ) );
    }

    // Scale vector
    const double factor = 0.5;
    vec.Scale ( factor );

    // Get vector components and check for correct values
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( static_cast < double > ( i * factor ), vec.GetValue ( i ) );
    }
}

// Test scaling vector and adding another vector

TEST ( HypreVector, scale_add )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec.Add ( i, static_cast < double > ( i ) );
    }

    // Clone structure to different HypreVector
    HypreVector<double> vec2;
    vec2.CloneFromWithoutContent ( vec );

    // Check for correct dimensions
    EXPECT_EQ ( vec.size_local ( ), vec2.size_local ( ) );
    EXPECT_EQ ( vec.size_global ( ), vec2.size_global ( ) );

    // Set components of second vector
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec2.Add ( i, static_cast < double > ( i + 2 ) );
    }

    const double factor = 4.0;

    vec.ScaleAdd ( vec2, factor );

    // Get vector components and check for correct values
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( static_cast < double > ( i * factor + ( i + 2 ) ), vec.GetValue ( i ) );
    }

}

// Test norm2

TEST ( HypreVector, norm2 )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        vec.Add ( i, static_cast < double > ( i ) );
    }

    // Compute dot product
    double norm_res = vec.Norm2 ( );

    // Maximum component index
    const double N = static_cast < double > ( global_dim - 1 );
    const double res_expected = std::sqrt ( ( 1. / 6. ) * N * ( N + 1. ) * ( 2. * N + 1. ) );

    EXPECT_DOUBLE_EQ ( res_expected, norm_res );

}

/// Test zeroing a HypreVector

TEST ( HypreVector, zeros )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 10;
    int global_dim = num_procs * local_dim;

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

    HypreVector<double> vec;
    vec.Init ( comm, laCouplings );

    // compute first and last component of current process
    int local_ilower = rank * local_dim;
    int local_iupper = ( rank + 1 ) * local_dim - 1;

    // Set different vector components
    std::vector<int> ind;
    std::vector<double> val;
    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        ind.push_back ( i );
        val.push_back ( static_cast < double > ( i ) );
    }

    vec.Add ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( val ) );

    // Set vector to zero
    vec.Zeros ( );

    // Get vector components and check for correct values
    std::vector<double> val_res ( val.size ( ) );
    vec.GetValues ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( val_res ) );

    for ( int i = local_ilower; i <= local_iupper; ++i )
    {
        EXPECT_DOUBLE_EQ ( 0., val_res[i - local_ilower] );
    }
}
