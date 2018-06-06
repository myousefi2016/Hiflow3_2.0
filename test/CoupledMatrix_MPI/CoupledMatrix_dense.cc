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
/// @date 2015-09-15

#include "hiflow.h"
#include "gtest/gtest.h"
#include "mpi.h"

using namespace hiflow::la;

/// Test multiplication of a CoupledMatrix with a CoupledVector y = A * x;
/// Structure is build up without FEM mesh.
/// The CoupledMatrix is fully dense occupied.
/// To test an arbitrary number of cores the global dimension is chosen as a multiple of the local dimension.

TEST ( CoupledMatrix_dense, FullGlobalMatrix )
{
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

    CoupledMatrix<double> A;
#ifdef WITH_OPENMP
    A.Init ( comm, laCouplings, CPU, OPENMP, CSR );
#else
    A.Init ( comm, laCouplings, CPU, NAIVE, CSR );
#endif
    A.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), nnz_diag,
                      vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), nnz_offdiag );

    if ( rank == 0 ) A.Add ( 0, 0, 2.0 );

    CoupledVector<double> x;
#ifdef WITH_OPENMP
    x.Init ( comm, laCouplings, CPU, OPENMP );
#else
    x.Init ( comm, laCouplings, CPU, NAIVE );
#endif
    x.InitStructure ( );

    if ( rank == 0 ) x.Add ( 0, 2.0 );

    CoupledVector<double> y;
#ifdef WITH_OPENMP
    y.Init ( comm, laCouplings, CPU, OPENMP );
#else
    y.Init ( comm, laCouplings, CPU, NAIVE );
#endif
    y.InitStructure ( );

    A.VectorMult ( x, &y );

    // Check result
    if ( rank == 0 )
    {
        int indices[global_dim];
        for ( int i = 0; i < global_dim; ++i ) indices[i] = i;
        double values[global_dim];
        y.GetValues ( &indices[0], global_dim, &values[0] );

        EXPECT_NEAR ( 4.0, values[0], 1e-12 );
    }
}
