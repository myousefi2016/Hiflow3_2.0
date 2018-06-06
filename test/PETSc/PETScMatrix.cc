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
/// @date 2015-11-26

#include "hiflow.h"
#include "gtest/gtest.h"
#include "mpi.h"

using namespace hiflow::la;

/// Test fixture for PETScMatrixTest

class PETScMatrixTest : public ::testing::Test
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

        //std::cout << "laCouplings: \n" << laCouplings << std::endl;

        // Generate COO structure for diagonal matrix
        nnz_diag = local_dim * local_dim;
        rows_diag.resize ( nnz_diag );
        cols_diag.resize ( nnz_diag );

        // Global indices
        for ( int r = 0; r < local_dim; ++r )
        {
            for ( int c = 0; c < local_dim; ++c )
            {
                rows_diag[r * local_dim + c] = r + global_offsets[mpi_rank];
                cols_diag[r * local_dim + c] = c + global_offsets[mpi_rank];
            }
        }

        // Generate COO structure for off-diagonal matrix
        nnz_offdiag = offdiag_cols_dim * local_dim;
        rows_offdiag.resize ( nnz_offdiag );
        cols_offdiag.resize ( nnz_offdiag );

        for ( int r = 0; r < local_dim; ++r )
        {
            for ( int c = 0; c < mpi_rank * local_dim; ++c )
            {
                rows_offdiag[r * offdiag_cols_dim + c] = r + global_offsets[mpi_rank];
                cols_offdiag[r * offdiag_cols_dim + c] = c;
            }
            for ( int c = ( mpi_rank + 1 ) * local_dim; c < mpi_num_procs * local_dim; ++c )
            {
                rows_offdiag[r * offdiag_cols_dim + c - local_dim] = r + global_offsets[mpi_rank];
                cols_offdiag[r * offdiag_cols_dim + c - local_dim] = c;
            }
        }
    }

    virtual void TearDown ( )
    {
        laCouplings.Clear ( );
        std::cout.rdbuf ( old_buffer );
    }

    void set_values ( PETScMatrix<double>& m ) const
    {
        //for (int i = 0; i != global_dim; ++i) m.SetValue(i, i, i + 0.2);
        for ( int i = 0; i != local_dim; ++i )
        {
            int global_index = 2 * mpi_rank + i;
            m.SetValue ( global_index, global_index, global_index + 0.2 );
        }
    }

    MPI_Comm comm;
    int mpi_rank;
    int mpi_num_procs;

    int local_dim;
    int global_dim;

    LaCouplings laCouplings;

    int nnz_diag;
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;

    int nnz_offdiag;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    std::ostringstream oss;
    std::streambuf* old_buffer;
};

TEST_F ( PETScMatrixTest, default_constructor )
{
    PETScMatrix<double> m;
}

TEST_F ( PETScMatrixTest, init )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
}

TEST_F ( PETScMatrixTest, init_structure )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
}

TEST_F ( PETScMatrixTest, set_value )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    set_values ( m );
}

TEST_F ( PETScMatrixTest, get_value )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    set_values ( m );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        EXPECT_DOUBLE_EQ ( global_index + 0.2, m.GetValue ( global_index, global_index ) );
    }
}

TEST_F ( PETScMatrixTest, norm_frobenius )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    set_values ( m );
    double norm = 0.0;
    for ( int i = 0; i != global_dim; ++i ) norm += std::pow ( i + 0.2, 2 );
    norm = std::sqrt ( norm );
    EXPECT_DOUBLE_EQ ( norm, m.NormFrobenius ( ) );
}

TEST_F ( PETScMatrixTest, zeros )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    set_values ( m );
    m.Zeros ( );
    EXPECT_DOUBLE_EQ ( 0.0, m.NormFrobenius ( ) );
}

TEST_F ( PETScMatrixTest, scale )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    set_values ( m );
    m.Scale ( 2.0 );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        EXPECT_DOUBLE_EQ ( 2 * ( global_index + 0.2 ), m.GetValue ( global_index, global_index ) );
    }
}

TEST_F ( PETScMatrixTest, diagonalize_rows )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    set_values ( m );

    double value = 10.0;
    std::vector<int> diag_indices;
    for ( int i = 0; i != local_dim; ++i ) diag_indices.push_back ( 2 * mpi_rank + i );

    m.diagonalize_rows ( vec2ptr ( diag_indices ), diag_indices.size ( ), value );

    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        EXPECT_DOUBLE_EQ ( value, m.GetValue ( global_index, global_index ) );
    }
}

TEST_F ( PETScMatrixTest, mat_vec_mult )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    m.Zeros ( );
    set_values ( m );

    PETScVector<double> v1;
    v1.Init ( comm, laCouplings );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        v1.SetValue ( global_index, global_index + 0.2 );
    }

    PETScVector<double> v2;
    v2.Init ( comm, laCouplings );
    m.VectorMult ( v1, &v2 );

    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        EXPECT_DOUBLE_EQ ( global_index + 0.2, m.GetValue ( global_index, global_index ) );
    }
}
