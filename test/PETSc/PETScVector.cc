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
/// @date 2015-11-17

#include "hiflow.h"
#include "gtest/gtest.h"
#include "mpi.h"

using namespace hiflow::la;

/// Test fixture for PETScVectorTest

class PETScVectorTest : public ::testing::Test
{
  protected:

    virtual void SetUp ( )
    {
        // Redirect std::cout into stringstream
        old_buffer = std::cout.rdbuf ( );
        std::cout.rdbuf ( oss.rdbuf ( ) );

        comm = MPI_COMM_WORLD;
        ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &mpi_rank ) ) );
        ASSERT_FALSE ( ( MPI_Comm_size ( comm, &mpi_num_procs ) ) );

        // Setup couplings object
        laCouplings.Init ( comm );

        // Generate global offsets
        local_dim = 2;
        global_dim = mpi_num_procs * local_dim;

        std::vector<int> global_offsets ( mpi_num_procs + 1 );
        global_offsets[0] = 0;
        for ( int i = 0; i < mpi_num_procs; ++i ) global_offsets[i + 1] = global_offsets[i] + local_dim;

        // Generate offdiag offsets
        // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
        std::vector<int> offdiag_offsets ( mpi_num_procs + 1 );
        offdiag_offsets[0] = 0;
        for ( int i = 0; i < mpi_num_procs; ++i ) offdiag_offsets[i + 1] = offdiag_offsets[i] + ( mpi_rank == i ? 0 : local_dim );

        // Generate offdiag cols
        int offdiag_cols_dim = ( mpi_num_procs - 1 ) * local_dim;
        std::vector<int> offdiag_cols ( offdiag_cols_dim );
        for ( int i = 0; i < mpi_rank * local_dim; ++i ) offdiag_cols[i] = i;
        for ( int i = ( mpi_rank + 1 ) * local_dim; i < mpi_num_procs * local_dim; ++i ) offdiag_cols[i - local_dim] = i;

        // Initialize laCouplings
        laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );
    }

    virtual void TearDown ( )
    {
        laCouplings.Clear ( );
        std::cout.rdbuf ( old_buffer );
    }

    void set_values ( PETScVector<double>& v ) const
    {
        for ( int i = 0; i != local_dim; ++i )
        {
            int global_index = 2 * mpi_rank + i;
            v.SetValue ( global_index, global_index + 0.2 );
        }
    }

    MPI_Comm comm;
    int mpi_rank;
    int mpi_num_procs;

    int local_dim;
    int global_dim;

    LaCouplings laCouplings;

    std::ostringstream oss;
    std::streambuf* old_buffer;
};

TEST_F ( PETScVectorTest, default_constructor )
{
    PETScVector<double> v;
}

TEST_F ( PETScVectorTest, init )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
}

TEST_F ( PETScVectorTest, size_local )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    EXPECT_EQ ( local_dim, v.size_local ( ) );
}

TEST_F ( PETScVectorTest, size_global )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    EXPECT_EQ ( global_dim, v.size_global ( ) );
}

TEST_F ( PETScVectorTest, set_value )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    set_values ( v );
}

TEST_F ( PETScVectorTest, zeros )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    for ( int i = 0; i != global_dim; ++i ) v.SetValue ( i, i + 0.2 );
    v.Zeros ( );
    EXPECT_DOUBLE_EQ ( 0.0, v.Norm1 ( ) );
}

TEST_F ( PETScVectorTest, get_value )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    set_values ( v );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        EXPECT_DOUBLE_EQ ( global_index + 0.2, v.GetValue ( global_index ) );
    }
}

TEST_F ( PETScVectorTest, add_value )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    set_values ( v );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        v.Add ( global_index, global_index + 0.2 );
    }
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        EXPECT_DOUBLE_EQ ( 2 * ( global_index + 0.2 ), v.GetValue ( global_index ) );
    }
}

TEST_F ( PETScVectorTest, norm1 )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    double value, norm = 0.0;
    for ( int i = 0; i != global_dim; ++i )
    {
        value = i + 0.2;
        v.SetValue ( i, value );
        norm += std::abs ( value );
    }
    EXPECT_DOUBLE_EQ ( norm, v.Norm1 ( ) );
}

TEST_F ( PETScVectorTest, norm2 )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    double value, norm = 0.0;
    for ( int i = 0; i != global_dim; ++i )
    {
        value = i + 0.2;
        v.SetValue ( i, value );
        norm += value * value;
    }
    norm = std::sqrt ( norm );
    EXPECT_DOUBLE_EQ ( norm, v.Norm2 ( ) );
}

TEST_F ( PETScVectorTest, normMax )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    double value, norm = 0.0;
    for ( int i = 0; i != global_dim; ++i )
    {
        value = i + 0.2;
        v.SetValue ( i, value );
        norm = std::max ( norm, std::abs ( value ) );
    }
    EXPECT_DOUBLE_EQ ( norm, v.NormMax ( ) );
}

TEST_F ( PETScVectorTest, dot )
{
    PETScVector<double> v1;
    v1.Init ( comm, laCouplings );
    v1.SetValue ( 0, 7.2 );
    v1.SetValue ( 1, 3.7 );
    PETScVector<double> v2;
    v2.Init ( comm, laCouplings );
    v2.SetValue ( 0, 2.7 );
    v2.SetValue ( 1, 7.3 );
    double result = 7.2 * 2.7 + 3.7 * 7.3;
    EXPECT_DOUBLE_EQ ( result, v1.Dot ( v2 ) );
    EXPECT_DOUBLE_EQ ( result, v2.Dot ( v1 ) );
}

TEST_F ( PETScVectorTest, scale )
{
    PETScVector<double> v;
    v.Init ( comm, laCouplings );
    for ( int i = 0; i != global_dim; ++i ) v.SetValue ( i, i + 0.2 );
    v.Scale ( 2.0 );
    EXPECT_DOUBLE_EQ ( ( ( mpi_rank * local_dim + 0.2 ) * 2 ), v.GetValue ( mpi_rank * local_dim ) );
    EXPECT_DOUBLE_EQ ( ( ( mpi_rank * local_dim + 1.2 ) * 2 ), v.GetValue ( mpi_rank * local_dim + 1 ) );
}
