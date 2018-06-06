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
/// @date 2015-09-29

#include "hiflow.h"
#include "linear_solver/amg/DirectInterpolation.h"
#include "linear_solver/amg/LevelGenerator.h"
#include "linear_solver/amg/PretendCoarsening.h"
#include "linear_solver/amg/StandardInterpolation.h"
#include "linear_solver/amg/StandardCoarsening.h"
#include "linear_solver/amg/StandardCoarseningInterpolationConnection.h"
#include "linear_solver/amg/TripleProduct.h"
#include "gtest/gtest.h"
#include "mpi.h"

using namespace hiflow::la;
using namespace hiflow::AMG;

TEST ( LevelGeneratorTest, A1_direct_interpolation )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            DirectInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A1.txt" ).c_str ( ) );

    LevelGeneratorType::ResultType level1 = LevelGeneratorType ( )( A );

    // Check
    EXPECT_EQ ( 3, level1.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 3, level1.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 7, level1.ptr_coarse_matrix->get_nnz ( ) );
    EXPECT_NEAR ( 1.0, level1.ptr_coarse_matrix->matrix.val[0], 1e-8 );
    EXPECT_NEAR ( 2.0, level1.ptr_coarse_matrix->matrix.val[1], 1e-8 );
    EXPECT_NEAR ( 3.0, level1.ptr_coarse_matrix->matrix.val[2], 1e-8 );
    EXPECT_NEAR ( 4.0, level1.ptr_coarse_matrix->matrix.val[3], 1e-8 );
    EXPECT_NEAR ( 5.0, level1.ptr_coarse_matrix->matrix.val[4], 1e-8 );
    EXPECT_NEAR ( 6.0, level1.ptr_coarse_matrix->matrix.val[5], 1e-8 );
    EXPECT_NEAR ( 7.0, level1.ptr_coarse_matrix->matrix.val[6], 1e-8 );
}

TEST ( LevelGeneratorTest, A1_std_interpolation )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            StandardInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A1.txt" ).c_str ( ) );

    LevelGeneratorType::ResultType level1 = LevelGeneratorType ( )( A );

    // Check
    EXPECT_EQ ( 3, level1.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 3, level1.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 7, level1.ptr_coarse_matrix->get_nnz ( ) );
    EXPECT_NEAR ( 1.0, level1.ptr_coarse_matrix->matrix.val[0], 1e-8 );
    EXPECT_NEAR ( 2.0, level1.ptr_coarse_matrix->matrix.val[1], 1e-8 );
    EXPECT_NEAR ( 3.0, level1.ptr_coarse_matrix->matrix.val[2], 1e-8 );
    EXPECT_NEAR ( 4.0, level1.ptr_coarse_matrix->matrix.val[3], 1e-8 );
    EXPECT_NEAR ( 5.0, level1.ptr_coarse_matrix->matrix.val[4], 1e-8 );
    EXPECT_NEAR ( 6.0, level1.ptr_coarse_matrix->matrix.val[5], 1e-8 );
    EXPECT_NEAR ( 7.0, level1.ptr_coarse_matrix->matrix.val[6], 1e-8 );
}

TEST ( LevelGeneratorTest, A2_direct_interpolation )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            StandardInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A2.txt" ).c_str ( ) );

    LevelGeneratorType::ResultType level1 = LevelGeneratorType ( )( A );

    // Check dimension
    EXPECT_EQ ( 56, level1.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 56, level1.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 460, level1.ptr_coarse_matrix->get_nnz ( ) );

    LevelGeneratorType::ResultType level2 = LevelGeneratorType ( )( *level1.ptr_coarse_matrix );

    // Check dimension
    EXPECT_EQ ( 23, level2.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 23, level2.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 281, level2.ptr_coarse_matrix->get_nnz ( ) );
}

TEST ( LevelGeneratorTest, A2_standard_interpolation )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            StandardInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A2.txt" ).c_str ( ) );

    LevelGeneratorType::ResultType level1 = LevelGeneratorType ( )( A );

    // Check dimension
    EXPECT_EQ ( 56, level1.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 56, level1.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 460, level1.ptr_coarse_matrix->get_nnz ( ) );

    LevelGeneratorType::ResultType level2 = LevelGeneratorType ( )( *level1.ptr_coarse_matrix );

    // Check dimension
    EXPECT_EQ ( 23, level2.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 23, level2.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 281, level2.ptr_coarse_matrix->get_nnz ( ) );
}

TEST ( LevelGeneratorTest, A3_direct_interpolation )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            DirectInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A3.mtx" ).c_str ( ) );

    LevelGeneratorType::ResultType level1 = LevelGeneratorType ( )( A );

    // Check
    EXPECT_EQ ( 2, level1.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 2, level1.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 4, level1.ptr_coarse_matrix->get_nnz ( ) );
    EXPECT_NEAR ( 1.0, level1.ptr_coarse_matrix->matrix.val[0], 1e-8 );
    EXPECT_NEAR ( -0.5, level1.ptr_coarse_matrix->matrix.val[1], 1e-8 );
    EXPECT_NEAR ( -0.5, level1.ptr_coarse_matrix->matrix.val[2], 1e-8 );
    EXPECT_NEAR ( 1.0, level1.ptr_coarse_matrix->matrix.val[3], 1e-8 );
}

TEST ( LevelGeneratorTest, A3_standard_interpolation )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            StandardInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A3.mtx" ).c_str ( ) );

    LevelGeneratorType::ResultType level1 = LevelGeneratorType ( )( A );

    // Check
    EXPECT_EQ ( 2, level1.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 2, level1.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 4, level1.ptr_coarse_matrix->get_nnz ( ) );
    EXPECT_NEAR ( 7.5, level1.ptr_coarse_matrix->matrix.val[0], 1e-8 );
    EXPECT_NEAR ( 1.5, level1.ptr_coarse_matrix->matrix.val[1], 1e-8 );
    EXPECT_NEAR ( 1.5, level1.ptr_coarse_matrix->matrix.val[2], 1e-8 );
    EXPECT_NEAR ( 7.5, level1.ptr_coarse_matrix->matrix.val[3], 1e-8 );
}

TEST ( LevelGeneratorTest, A4_direct_interpolation )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            DirectInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A4.mtx" ).c_str ( ) );

    LevelGeneratorType::ResultType level1 = LevelGeneratorType ( )( A );

    // Check
    EXPECT_EQ ( 1, level1.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 1, level1.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 1, level1.ptr_coarse_matrix->get_nnz ( ) );
    EXPECT_NEAR ( 0.57380952380952388, level1.ptr_coarse_matrix->matrix.val[0], 1e-8 );
}

TEST ( LevelGeneratorTest, A4_standard_interpolation )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            StandardInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A4.mtx" ).c_str ( ) );

    LevelGeneratorType::ResultType level1 = LevelGeneratorType ( )( A );

#if 0
    std::cout << "P: ";
    std::cout << "row: ";
    for ( size_t i ( 0 ); i != level1.ptr_interpolation_matrix->get_num_row ( ) + 1; ++i ) std::cout << level1.ptr_interpolation_matrix->matrix.row[i] << " ";
    std::cout << std::endl;
    std::cout << "col: ";
    for ( size_t i ( 0 ); i != level1.ptr_interpolation_matrix->get_nnz ( ); ++i ) std::cout << level1.ptr_interpolation_matrix->matrix.col[i] << " ";
    std::cout << std::endl;
    std::cout << "val: ";
    for ( size_t i ( 0 ); i != level1.ptr_interpolation_matrix->get_nnz ( ); ++i ) std::cout << level1.ptr_interpolation_matrix->matrix.val[i] << " ";
    std::cout << std::endl;

    std::cout << "R: ";
    std::cout << "row: ";
    for ( size_t i ( 0 ); i != level1.ptr_restriction_matrix->get_num_row ( ) + 1; ++i ) std::cout << level1.ptr_restriction_matrix->matrix.row[i] << " ";
    std::cout << std::endl;
    std::cout << "col: ";
    for ( size_t i ( 0 ); i != level1.ptr_restriction_matrix->get_nnz ( ); ++i ) std::cout << level1.ptr_restriction_matrix->matrix.col[i] << " ";
    std::cout << std::endl;
    std::cout << "val: ";
    for ( size_t i ( 0 ); i != level1.ptr_restriction_matrix->get_nnz ( ); ++i ) std::cout << level1.ptr_restriction_matrix->matrix.val[i] << " ";
    std::cout << std::endl;
#endif

    // Check
    EXPECT_EQ ( 1, level1.ptr_coarse_matrix->get_num_row ( ) );
    EXPECT_EQ ( 1, level1.ptr_coarse_matrix->get_num_col ( ) );
    EXPECT_EQ ( 1, level1.ptr_coarse_matrix->get_nnz ( ) );
    EXPECT_NEAR ( 5.6, level1.ptr_coarse_matrix->matrix.val[0], 1e-8 );
}

#if 0

TEST ( LevelGeneratorTest, CoupledStandardCoarsening )
{
    typedef CoupledMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            StandardInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int local_dim = 2;
    int global_dim = num_procs * local_dim;

    std::vector<int> global_offsets ( num_procs + 1 );
    global_offsets[0] = 0;
    for ( int i = 0; i < num_procs; ++i ) global_offsets[i + 1] = global_offsets[i] + local_dim;

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1 );
    offdiag_offsets[0] = 0;
    for ( int i = 0; i < num_procs; ++i ) offdiag_offsets[i + 1] = offdiag_offsets[i] + ( rank == i ? 0 : local_dim );

    //for (int i=0; i < num_procs + 1; ++i)
    //  std::cout << "[" << rank << "] offdiag_offsets" << "[" << i << "] = " << offdiag_offsets[i] << std::endl;

    // Generate offdiag cols
    int offdiag_cols_dim = ( num_procs - 1 ) * local_dim;
    std::vector<int> offdiag_cols ( offdiag_cols_dim );
    for ( int i = 0; i < rank * local_dim; ++i ) offdiag_cols[i] = i;
    for ( int i = ( rank + 1 ) * local_dim; i < num_procs * local_dim; ++i ) offdiag_cols[i - local_dim] = i;

    //for (int i=0; i < offdiag_cols_dim; ++i)
    //  std::cout << "[" << rank << "] offdiag_cols" << "[" << i << "] = " << offdiag_cols[i] << std::endl;

    // Initialize laCouplings
    laCouplings.InitializeCouplings ( global_offsets, offdiag_cols, offdiag_offsets );

    // Generate COO structure for diagonal matrix
    int nnz_diag = local_dim * local_dim;
    std::vector<int> rows_diag ( nnz_diag );
    std::vector<int> cols_diag ( nnz_diag );

    // Global indices
    for ( int r = 0; r < local_dim; ++r )
    {
        for ( int c = 0; c < local_dim; ++c )
        {
            rows_diag[r * local_dim + c] = r + global_offsets[rank];
            cols_diag[r * local_dim + c] = c + global_offsets[rank];
        }
    }

    // Generate COO structure for off-diagonal matrix
    int nnz_offdiag = offdiag_cols_dim * local_dim;
    std::vector<int> rows_offdiag ( nnz_offdiag );
    std::vector<int> cols_offdiag ( nnz_offdiag );

    for ( int r = 0; r < local_dim; ++r )
    {
        for ( int c = 0; c < rank * local_dim; ++c )
        {
            rows_offdiag[r * offdiag_cols_dim + c] = r + global_offsets[rank];
            cols_offdiag[r * offdiag_cols_dim + c] = c;
        }
        for ( int c = ( rank + 1 ) * local_dim; c < num_procs * local_dim; ++c )
        {
            rows_offdiag[r * offdiag_cols_dim + c - local_dim] = r + global_offsets[rank];
            cols_offdiag[r * offdiag_cols_dim + c - local_dim] = c;
        }
    }

    MatrixType A;
    A.Init ( comm, laCouplings, CPU, OPENMP, CSR );
    A.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), nnz_diag,
                      vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), nnz_offdiag );

    StandardCoarseningSettings coarsening_settings;
    coarsening_settings.strength_threshold = 0.5;

    LevelGeneratorType::ResultType coarsening_level = LevelGeneratorType ( coarsening_settings )( A );

    //coarsening_level.coarse_matrix->Print();
}
#endif
