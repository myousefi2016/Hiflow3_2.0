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
/// @date 2015-09-02

#include "hiflow.h"
#include "gtest/gtest.h"
#include "mpi.h"

using namespace hiflow::la;

/// Test multiplication of a CoupledMatrix with a CoupledVector y = A * x;
/// Structure is build up without FEM mesh.
/// Only working using two MPI workers.
/// The local matrices are dense occupied.

TEST ( CoupledMatrix_mesh, NoOffdiagCoupling )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );
    ASSERT_EQ ( 2, num_procs );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int global_dim = 9;
    int local_dim = rank == 0 ? 6 : 3;

    std::vector<int> global_offsets ( num_procs + 1 );
    global_offsets[0] = 0;
    global_offsets[1] = 6;
    global_offsets[2] = 9;

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1 );
    offdiag_offsets[0] = 0;
    offdiag_offsets[1] = 0;
    offdiag_offsets[2] = 0;

    // Generate offdiag cols
    std::vector<int> offdiag_cols;

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
    int nnz_offdiag = 0;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    CoupledMatrix<double> A;
#ifdef WITH_OPENMP
    A.Init ( comm, laCouplings, CPU, OPENMP, CSR );
#else
    A.Init ( comm, laCouplings, CPU, NAIVE, CSR );
#endif
    A.InitStructure ( vec2ptr ( rows_diag ), vec2ptr ( cols_diag ), nnz_diag,
                      vec2ptr ( rows_offdiag ), vec2ptr ( cols_offdiag ), nnz_offdiag );

    CoupledVector<double> x;
#ifdef WITH_OPENMP
    x.Init ( comm, laCouplings, CPU, OPENMP );
#else
    x.Init ( comm, laCouplings, CPU, NAIVE );
#endif
    x.InitStructure ( );

    CoupledVector<double> y;
#ifdef WITH_OPENMP
    y.Init ( comm, laCouplings, CPU, OPENMP );
#else
    y.Init ( comm, laCouplings, CPU, NAIVE );
#endif
    y.InitStructure ( );

    A.VectorMult ( x, &y );

}

/// Test multiplication of a CoupledMatrix with a CoupledVector y = A * x;
/// Structure is build up without FEM mesh.
/// Only working using two MPI workers.
/// The local matrices are dense occupied.

TEST ( CoupledMatrix_mesh, OffdiagCoupling )
{
    MPI_Comm comm ( MPI_COMM_WORLD );

    int rank, num_procs;
    ASSERT_FALSE ( ( MPI_Comm_rank ( comm, &rank ) ) );
    ASSERT_FALSE ( ( MPI_Comm_size ( comm, &num_procs ) ) );
    ASSERT_EQ ( 2, num_procs );

    // Setup couplings object
    LaCouplings laCouplings;
    laCouplings.Init ( comm );

    // Generate global offsets
    int global_dim = 9;
    int local_dim = rank == 0 ? 6 : 3;

    std::vector<int> global_offsets ( num_procs + 1 );
    global_offsets[0] = 0;
    global_offsets[1] = 6;
    global_offsets[2] = 9;

    // Generate offdiag offsets
    // The first rank don't have off-diagonals since all border vertices are allocated by the lowest rank.
    std::vector<int> offdiag_offsets ( num_procs + 1 );
    offdiag_offsets[0] = 0;
    offdiag_offsets[1] = rank == 0 ? 0 : 3;
    offdiag_offsets[2] = rank == 0 ? 0 : 3;

    // Generate offdiag cols
    std::vector<int> offdiag_cols;
    if ( rank == 1 )
    {
        offdiag_cols.resize ( 3 );
        offdiag_cols[0] = 3;
        offdiag_cols[1] = 4;
        offdiag_cols[2] = 5;
    }

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
    int nnz_offdiag = 0;
    std::vector<int> rows_offdiag;
    std::vector<int> cols_offdiag;

    if ( rank == 1 )
    {
        nnz_offdiag = local_dim * local_dim;
        rows_offdiag.resize ( nnz_offdiag );
        cols_offdiag.resize ( nnz_offdiag );

        for ( int r = 0; r < local_dim; ++r )
        {
            for ( int c = 0; c < local_dim; ++c )
            {
                rows_offdiag[r * local_dim + c] = r + global_offsets[rank];
                cols_offdiag[r * local_dim + c] = c + local_dim;
            }
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

    CoupledVector<double> x;
#ifdef WITH_OPENMP
    x.Init ( comm, laCouplings, CPU, OPENMP );
#else
    x.Init ( comm, laCouplings, CPU, NAIVE );
#endif
    x.InitStructure ( );

    CoupledVector<double> y;
#ifdef WITH_OPENMP
    y.Init ( comm, laCouplings, CPU, OPENMP );
#else
    y.Init ( comm, laCouplings, CPU, NAIVE );
#endif
    y.InitStructure ( );

    A.VectorMult ( x, &y );

}
