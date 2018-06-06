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
/// @date 2015-12-23

#include "cxx98_prettyprint.h"
#include "hiflow.h"
#include "gtest/gtest.h"
#include "mpi.h"
#include <algorithm>

using namespace hiflow::la;

namespace
{

    struct InputData
    {

        InputData ( std::string const& test_set, petsc::KSPType solver_type )
        : test_set ( test_set ), solver_type ( solver_type )
        {
        }

        std::string test_set;
        petsc::KSPType solver_type;
    };

    void get_sub_coo_structure ( CPUsimple_CSR_lMatrix<double> const& matrix,
                                 std::vector<int> const& rows, std::vector<int> const& cols,
                                 std::vector<int> &rows_indices, std::vector<int> &cols_indices,
                                 std::vector<double> &vals )
    {
        for ( int r = 0; r != matrix.get_num_row ( ); ++r )
        {
            if ( find ( rows.begin ( ), rows.end ( ), r ) == rows.end ( ) ) continue;
            for ( int p = matrix.matrix.row[r]; p != matrix.matrix.row[r + 1]; ++p )
            {
                int c = matrix.matrix.col[p];
                if ( find ( cols.begin ( ), cols.end ( ), c ) == cols.end ( ) ) continue;
                rows_indices.push_back ( r );
                cols_indices.push_back ( c );
                vals.push_back ( matrix.matrix.val[p] );
            }
        }
    }

} // namespace anonymous

/// Test fixture

class PETScLinearSolverTestReadData
: public ::testing::TestWithParam<InputData>
{
};

/// Test body

TEST_P ( PETScLinearSolverTestReadData, PETScGeneralKSP_solver )
{

    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_rank;
    int mpi_num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &mpi_rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &mpi_num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    CPUsimple_CSR_lMatrix<double> simple_matrix;
    simple_matrix.ReadFile ( ( std::string ( DATA_DIR ) + GetParam ( ).test_set + "/matrix.txt" ).c_str ( ) );

    // Generate global offsets
    // Distribute indices equally over all nodes
    int global_dim = simple_matrix.get_num_row ( );
    ASSERT_LE ( mpi_num_procs, global_dim );
    int local_dim = global_dim / mpi_num_procs;
    if ( global_dim % mpi_num_procs > mpi_rank ) ++local_dim;

    std::vector<int> global_offsets ( mpi_num_procs + 1 );
    global_offsets[0] = 0;
    for ( int i = 0; i < mpi_num_procs; ++i )
    {
        global_offsets[i + 1] = global_offsets[i] + global_dim / mpi_num_procs;
        if ( global_dim % mpi_num_procs > i ) ++global_offsets[i + 1];
    }

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are
    // allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( mpi_num_procs + 1 );
    offdiag_offsets[0] = 0;
    for ( int i = 0; i < mpi_num_procs; ++i )
    {
        if ( mpi_rank == i ) offdiag_offsets[i + 1] = offdiag_offsets[i];
        else
        {
            offdiag_offsets[i + 1] = offdiag_offsets[i] + global_dim / mpi_num_procs;
            if ( global_dim % mpi_num_procs > i ) ++offdiag_offsets[i + 1];
        }
    }

    // Generate offdiag cols
    std::vector<int> offdiag_cols ( global_dim - local_dim );
    {
        int j = 0;
        for ( int i = 0; i < global_offsets[mpi_rank]; ++i, ++j )
            offdiag_cols[j] = i;
        for ( int i = global_offsets[mpi_rank + 1]; i < global_dim; ++i, ++j )
            offdiag_cols[j] = i;
    }

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    std::vector<int> rows ( local_dim );
    for ( int i = 0; i < local_dim; ++i ) rows[i] = global_offsets[mpi_rank] + i;

    std::vector<int> diag_cols ( local_dim );
    for ( int i = 0; i < local_dim; ++i ) diag_cols[i] = global_offsets[mpi_rank] + i;

#ifdef PRINT_DEBUG
    std::cout << "global_offsets " << mpi_rank << " = " << global_offsets << std::endl;
    std::cout << "rows " << mpi_rank << " = " << rows << std::endl;
    std::cout << "diag_cols " << mpi_rank << " = " << diag_cols << std::endl;
    std::cout << "offdiag_cols " << mpi_rank << " = " << offdiag_cols << std::endl;
#endif

    // Generate COO structure for diagonal matrix
    std::vector<int> rows_diag, cols_diag;
    std::vector<double> vals_diag;
    get_sub_coo_structure ( simple_matrix, rows, diag_cols, rows_diag, cols_diag, vals_diag );

#ifdef PRINT_DEBUG
    std::cout << "diag rows " << mpi_rank << " = " << rows_diag << std::endl;
    std::cout << "diag cols " << mpi_rank << " = " << cols_diag << std::endl;
    std::cout << "diag vals " << mpi_rank << " = " << vals_diag << std::endl;
#endif

    // Generate COO structure for off-diagonal matrix
    std::vector<int> rows_offdiag, cols_offdiag;
    std::vector<double> vals_offdiag;
    get_sub_coo_structure ( simple_matrix, rows, offdiag_cols, rows_offdiag, cols_offdiag, vals_offdiag );

#ifdef PRINT_DEBUG
    std::cout << "offdiag rows " << mpi_rank << " = " << rows_offdiag << std::endl;
    std::cout << "offdiag cols " << mpi_rank << " = " << cols_offdiag << std::endl;
    std::cout << "offdiag vals " << mpi_rank << " = " << vals_offdiag << std::endl;
#endif

    PETScMatrix<double> matrix;
    matrix.Init ( comm, laCouplings );
    matrix.InitStructure ( &rows_diag[0], &cols_diag[0], rows_diag.size ( ), &rows_offdiag[0],
                           &cols_offdiag[0], rows_offdiag.size ( ) );
    for ( int i = 0; i != vals_diag.size ( ); ++i )
    {
        matrix.SetValue ( rows_diag[i], cols_diag[i], vals_diag[i] );
    }
    for ( int i = 0; i != vals_offdiag.size ( ); ++i )
    {
        matrix.SetValue ( rows_offdiag[i], cols_offdiag[i], vals_offdiag[i] );
    }

    CPUsimple_lVector<double> simple_rhs;
    simple_rhs.ReadFile ( ( std::string ( DATA_DIR ) + GetParam ( ).test_set + "/rhs.txt" ).c_str ( ) );

    PETScVector<double> rhs;
    rhs.Init ( comm, laCouplings );
    rhs.SetValues ( &rows[0], local_dim, simple_rhs.GetBuffer ( ) + global_offsets[mpi_rank] );

    CPUsimple_lVector<double> simple_x0;
    simple_x0.ReadFile ( ( std::string ( DATA_DIR ) + GetParam ( ).test_set + "/x0.txt" ).c_str ( ) );

    PETScVector<double> x;
    x.Init ( comm, laCouplings );
    x.SetValues ( &rows[0], local_dim, simple_x0.GetBuffer ( ) + global_offsets[mpi_rank] );

    EXPECT_EQ ( local_dim, x.size_local ( ) );

    PETScGeneralKSP<LADescriptorPETSc> solver ( comm, matrix, GetParam ( ).solver_type );
    solver.InitControl ( 1000, 1e-12, 1e-12 );
    LinearSolverState state = solver.Solve ( rhs, &x );
    EXPECT_EQ ( kSolverSuccess, state );

    CPUsimple_lVector<double> simple_sol;
    simple_sol.ReadFile ( ( std::string ( DATA_DIR ) + GetParam ( ).test_set + "/sol.txt" ).c_str ( ) );

    for ( int i = 0; i != local_dim; ++i )
    {
        int global_index = global_offsets[mpi_rank] + i;
        EXPECT_NEAR ( simple_sol.GetBuffer ( )[global_index], x.GetValue ( global_index ), 1e-2 );
    }
}

INSTANTIATE_TEST_CASE_P ( AllTests, PETScLinearSolverTestReadData, ::testing::Values (
                                                                                       InputData ( "set_A", petsc::CG ),
                                                                                       InputData ( "set_A", petsc::GMRES ),
                                                                                       InputData ( "set_B", petsc::CG ),
                                                                                       InputData ( "set_B", petsc::GMRES ),
                                                                                       InputData ( "set_C", petsc::CG ),
                                                                                       InputData ( "set_C", petsc::GMRES )
                                                                                       ) );
