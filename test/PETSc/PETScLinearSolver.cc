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
/// @date 2015-12-02

#include "hiflow.h"
#include "gtest/gtest.h"
#include "mpi.h"

using namespace hiflow::la;

/// Test fixture for PETScLinearSolverTest

class PETScLinearSolverTest : public ::testing::Test
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

TEST_F ( PETScLinearSolverTest, general_cg )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    m.Zeros ( );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        m.SetValue ( global_index, global_index, global_index + 0.2 );
    }

    PETScVector<double> v1;
    v1.Init ( comm, laCouplings );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        v1.SetValue ( global_index, global_index + 0.2 );
    }

    PETScVector<double> v2;
    v2.Init ( comm, laCouplings );

    PETScGeneralKSP<LADescriptorPETSc> solver ( comm, m, petsc::CG );
    LinearSolverState state = solver.Solve ( v1, &v2 );
    EXPECT_EQ ( kSolverSuccess, state );

    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        EXPECT_DOUBLE_EQ ( 1.0, v2.GetValue ( global_index ) );
    }
}

TEST_F ( PETScLinearSolverTest, general_gmres )
{
    PETScMatrix<double> m;
    m.Init ( comm, laCouplings );
    m.InitStructure ( &rows_diag[0], &cols_diag[0], nnz_diag, &rows_offdiag[0], &cols_offdiag[0], nnz_offdiag );
    m.Zeros ( );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        m.SetValue ( global_index, global_index, global_index + 0.2 );
    }

    PETScVector<double> v1;
    v1.Init ( comm, laCouplings );
    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        v1.SetValue ( global_index, global_index + 0.2 );
    }

    PETScVector<double> v2;
    v2.Init ( comm, laCouplings );

    PETScGeneralKSP<LADescriptorPETSc> solver ( comm, m, petsc::GMRES );
    LinearSolverState state = solver.Solve ( v1, &v2 );
    EXPECT_EQ ( kSolverSuccess, state );

    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = 2 * mpi_rank + i;
        EXPECT_DOUBLE_EQ ( 1.0, v2.GetValue ( global_index ) );
    }
}
