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

/// @author Chandramowli Subramanian, Martin Wlotzka

#include "mumps_solver.h"

#include "config.h"
#include "common/log.h"
#include "dof/dof_partition.h"
#include "linear_algebra/la_descriptor.h"
#include "mumps_structure.h"

namespace hiflow
{
    namespace la
    {

        /// standard constructor

        template<class LAD, class MUMPS_STRUCTURE>
        MumpsSolver<LAD, MUMPS_STRUCTURE>::MumpsSolver ( )
        : LinearSolver<LAD>( )
        {
            this->comm_ = MPI_COMM_NULL;

            this->mumps_struc_initialized_ = false;
            this->transposed_ = false;
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Linear solver", "MUMPS" );
            }
        }

        /// destructor

        template<class LAD, class MUMPS_STRUCTURE>
        MumpsSolver<LAD, MUMPS_STRUCTURE>::~MumpsSolver ( )
        {
            this->Clear ( );
            int is_finalized;
            MPI_Finalized ( &is_finalized );
            if ( !is_finalized )
            {
                if ( this->comm_ != MPI_COMM_NULL )
                    MPI_Comm_free ( &this->comm_ );
            }
        }

        /// Inits parameter for Mumps structure.

        template<class LAD, class MUMPS_STRUCTURE>
        void MumpsSolver<LAD, MUMPS_STRUCTURE>::InitParameter ( const MPI_Comm& comm
                                                                /* more params */ )
        {
#ifdef WITH_MUMPS
            this->Clear ( );
            assert ( comm != MPI_COMM_NULL );
            int my_rank = -1;
            int info = MPI_Comm_rank ( comm, &my_rank );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank >= 0 );
            info = MPI_Comm_split ( comm, 0, my_rank, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            this->mumps_struc_.id_.comm_fortran = MPI_Comm_c2f ( comm ); // communicator

            this->mumps_struc_.id_.sym = 0; // unsymmetric problem
            // ( 1 for spd, 2 for general symmetric)
            this->mumps_struc_.id_.par = 1; // host participates in computations

            // initialize MUMPS package
            this->mumps_struc_.id_.job = -1;
            this->mumps_struc_.apply ( );
            this->mumps_struc_initialized_ = true;

            this->mumps_struc_.id_.nrhs = 1; // number of rhs

            // set up default ICNTL parameters (see mumps documentation)
            this->mumps_struc_.id_.ICNTL ( 1 ) = 6; // output stream for error messages (default 6)
            this->mumps_struc_.id_.ICNTL ( 2 ) = 0; // output stream for diagnostic (default 0)
            this->mumps_struc_.id_.ICNTL ( 3 ) = 0; // output stream for global information
            // (default 6)
            this->mumps_struc_.id_.ICNTL ( 4 ) = 2; // level of printing (default 2)
            this->mumps_struc_.id_.ICNTL ( 5 ) = 0; // assembled in coordinate format
            this->mumps_struc_.id_.ICNTL ( 6 ) = 7; // column permutation and/or scaling
            // to get a zero-free diagonal
            // (default 7, i.e. automatic choice)
            this->mumps_struc_.id_.ICNTL ( 7 ) = 7; // matrix ordering
            // (default 7, automatic choice)
            this->mumps_struc_.id_.ICNTL ( 8 ) = 77; // scaling strategy
            // (default 77, automatic choice)
            this->mumps_struc_.id_.ICNTL ( 9 ) = 1; // 1 : Ax=b; other than 1: A^T x=b
            this->mumps_struc_.id_.ICNTL ( 10 ) = 0; // max number of refinement
            // (default 0, no refinement)
            this->mumps_struc_.id_.ICNTL ( 11 ) = 0; // if positive returns statistics (default 0)
            this->mumps_struc_.id_.ICNTL ( 12 ) = 0; // ordering strategy (only for SYM=2)
            this->mumps_struc_.id_.ICNTL ( 13 ) = 0; // if <= 0, use ScaLAPACK for root
            // frontal matrix (default 0)
            //this->mumps_struc_.id_.ICNTL(14) = 20; // percentage of estimated workspace increase (automatic chosen)

            // 15-17 not used

            this->mumps_struc_.id_.ICNTL ( 18 ) = 3; // = 3, distributed assembled matrix input
            this->mumps_struc_.id_.ICNTL ( 19 ) = 0; // return Schur complement (default 0, no return)
            this->mumps_struc_.id_.ICNTL ( 20 ) = 0; // rhs sparse pattern (0 = dense, 1 = sparse)
            this->mumps_struc_.id_.ICNTL ( 21 ) = 0; // solution struc
            // (0 = centralized, 1 = distributed)
            this->mumps_struc_.id_.ICNTL ( 22 ) = 0; // in-core/out-of-core facility
            // (default 0, in-core)
            this->mumps_struc_.id_.ICNTL ( 23 ) = 0; // max size of memory that can be
            // allocated locally
            // (default 0, automatic choice based on
            // estimation)
            this->mumps_struc_.id_.ICNTL ( 24 ) = 0; // detection of null pivot rows
            // (0 nothing done, 1 detection)
            this->mumps_struc_.id_.ICNTL ( 25 ) = 0; // computation of a null space basis
            // (default 0, one possible solution is returned)
            this->mumps_struc_.id_.ICNTL ( 26 ) = 0; // Schur complement options (default 0)

            // 27 experimental

            this->mumps_struc_.id_.ICNTL ( 28 ) = 1; // analyis phase
            // (0 automatic, 1 sequential, 2 parallel)
            this->mumps_struc_.id_.ICNTL ( 29 ) = 0; // ordering tool for parallel analyis

            // now set up customized parameters
            // TODO
#else
            printf ( "MumpsSolver::InitParameter: No MUMPS support.\n" );
            exit ( -1 );
#endif
        }

        /// Inits the operator i.e. sets up the local matrix.

        template<class LAD, class MUMPS_STRUCTURE>
        void MumpsSolver<LAD, MUMPS_STRUCTURE>::SetupOperator ( OperatorType& op )
        {
            this->op_ = &op;
            assert ( this->op_ != NULL );
            assert ( this->op_->is_square_matrix ( ) );

            this->SetModifiedOperator ( true );
        }

        /// Inits the operator i.e. sets up the local matrix.

        template<class LAD, class MUMPS_STRUCTURE>
        void MumpsSolver<LAD, MUMPS_STRUCTURE>::Build ( )
        {
            assert ( this->op_ != NULL );
            assert ( this->op_->is_square_matrix ( ) );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Build Solver", 1 );
            }
            this->CreateLocalMatrix ( );
            this->FactorizeSymbolic ( );
            this->FactorizeNumeric ( );

            this->SetState ( true );
            this->SetModifiedOperator ( false );
        }

        /// Symbolic factorization.

        template<class LAD, class MUMPS_STRUCTURE>
        void MumpsSolver<LAD, MUMPS_STRUCTURE>::FactorizeSymbolic ( )
        {
#ifdef WITH_MUMPS
            // only distributed input matrix
            assert ( this->mumps_struc_.id_.ICNTL ( 18 ) == 3 );
            // host participates in computaions
            assert ( this->mumps_struc_.id_.par == 1 );

            // set up Mumps structure
            this->mumps_struc_.id_.nz_loc = this->op_->nnz_local ( ); // nnz per process
            this->mumps_struc_.id_.n = this->op_->nrows_global ( ); // global size
            // global row indices for local process
            this->mumps_struc_.id_.irn_loc = &( this->row_.front ( ) );
            // global column indices for local process
            this->mumps_struc_.id_.jcn_loc = &( this->col_.front ( ) );
            // values for local process
            this->mumps_struc_.id_.a_loc = &( this->val_.front ( ) );

            // now factorize
            this->mumps_struc_.id_.job = 1;
            this->mumps_struc_.apply ( );

#else
            LOG_ERROR ( "MumpsSolver::FactorizeSymbolic: No MUMPS support." );
            exit ( -1 );
#endif
        }

        /// Numeric factorization.

        template<class LAD, class MUMPS_STRUCTURE>
        void MumpsSolver<LAD, MUMPS_STRUCTURE>::FactorizeNumeric ( )
        {
#ifdef WITH_MUMPS
            // factorize
            this->mumps_struc_.id_.job = 2;
            this->mumps_struc_.apply ( );
#else
            LOG_ERROR ( "MumpsSolver::FactorizeNumeric: No MUMPS support." );
            exit ( -1 );
#endif
        }

        /// Solution step after factorization
        /// @param b right hand side vector
        /// @param x solution vector (needs to be allocated)

        template<class LAD, class MUMPS_STRUCTURE>
        LinearSolverState MumpsSolver<LAD, MUMPS_STRUCTURE>::SolveFactorized ( const VectorType& b,
                                                                               VectorType* x )
        {
            assert ( this->mumps_struc_initialized_ );

#ifdef WITH_MUMPS
            // MUMPS only supports centralized rhs -> gather
            int master_rank = 0;
            DataType* rhs_gathered = new DataType[this->mumps_struc_.id_.n];
            b.Gather ( master_rank, rhs_gathered );

            // size of rhs
            this->mumps_struc_.id_.lrhs = this->mumps_struc_.id_.n;

            this->mumps_struc_.id_.rhs = rhs_gathered; // centralized rhs and sol

            // choose if tranposed system is to be solved
            if ( this->transposed ( ) )
            {
                this->mumps_struc_.id_.ICNTL ( 9 ) = 0;
            }
            else
            {
                this->mumps_struc_.id_.ICNTL ( 9 ) = 1;
            }

            // now solve
            this->mumps_struc_.id_.job = 3;
            this->mumps_struc_.apply ( );

            if ( this->mumps_struc_.id_.INFO ( 1 ) != 0 )
            {
                return kSolverError;
            }

            int my_rank, size_comm;
            int info = MPI_Comm_rank ( this->comm ( ), &my_rank );
            assert ( info == MPI_SUCCESS );
            info = MPI_Comm_size ( this->comm ( ), &size_comm );
            assert ( info == MPI_SUCCESS );
            int tag = 1;
            DataType* sol_recv = new DataType[this->op_->nrows_local ( )];

            // communicate centralized solution
            if ( my_rank == master_rank )
            {

                int send_index_begin;
                int send_index_end = this->op_->row_couplings ( ).nb_dofs ( 0 );

                for ( int id = 1; id < size_comm; ++id )
                {
                    send_index_begin = send_index_end;
                    send_index_end += this->op_->row_couplings ( ).nb_dofs ( id );

                    info = MPI_Send ( &rhs_gathered[send_index_begin],
                                      send_index_end - send_index_begin,
                                      MPI_DOUBLE, id, tag, this->comm ( ) );
                    assert ( info == MPI_SUCCESS );
                }

                // own values
                for ( int i = 0; i<this->op_->nrows_local ( ); ++i )
                {
                    sol_recv[i] = rhs_gathered[i];
                }

            }
            else
            {
                MPI_Status status;
                info = MPI_Recv ( sol_recv, this->op_->nrows_local ( ), MPI_DOUBLE,
                                  master_rank, tag, this->comm ( ), &status );
                assert ( info == MPI_SUCCESS );
            }

            // now set values in solution vector
            x->SetValues ( sol_recv );

            delete[] rhs_gathered;
            delete[] sol_recv;

            // check residuum
#    ifndef NDEBUG
            VectorType res;
            res.CloneFromWithoutContent ( *x );
            res.Zeros ( );
            this->op_->VectorMult ( *x, &res );
            res.ScaleAdd ( b, -1. );
            if ( my_rank == master_rank )
                std::cout << "residuum of linear system: " << res.Norm2 ( ) << std::endl;
#    endif

            return kSolverSuccess;

#else
            LOG_ERROR ( "MumpsSolver::SolveFactorized: No MUMPS support." );
            exit ( -1 );
            return kSolverError;
#endif
        }

        /// Solves the linear system, includes factorization steps
        /// @param b right hand side vector
        /// @param x solution vector (needs to be allocated)

        template<class LAD, class MUMPS_STRUCTURE>
        LinearSolverState MumpsSolver<LAD, MUMPS_STRUCTURE>::Solve ( const VectorType& b,
                                                                     VectorType* x )
        {
#ifdef WITH_MUMPS
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            return this->SolveFactorized ( b, x );
#else
            LOG_ERROR ( "MumpsSolver::Solve: No MUMPS support." );
            exit ( -1 );
            return kSolverError;
#endif
        }

        /// Clears allocated data.

        template<class LAD, class MUMPS_STRUCTURE>
        void MumpsSolver<LAD, MUMPS_STRUCTURE>::Clear ( )
        {
            // deallocate MUMPS package
            if ( this->mumps_struc_initialized_ )
            {
#ifdef WITH_MUMPS
                this->mumps_struc_.id_.job = -2;
                this->mumps_struc_.apply ( );
#else
                LOG_ERROR ( "MumpsSolver::Clear: No MUMPS support." );
                exit ( -1 );
#endif
                this->mumps_struc_initialized_ = false;
            }

            this->row_.clear ( );
            this->col_.clear ( );
            this->val_.clear ( );
            LinearSolver<LAD>::Clear ( );
        }

        /// Creates local matrix in COO format.

        template<class LAD, class MUMPS_STRUCTURE>
        void MumpsSolver<LAD, MUMPS_STRUCTURE>::CreateLocalMatrix ( )
        {
            assert ( this->op_ != NULL );

            this->row_.resize ( this->op_->nnz_local ( ) );
            this->col_.resize ( this->op_->nnz_local ( ) );
            this->val_.resize ( this->op_->nnz_local ( ) );

            // extract CSR matrix
            int* ia = new int[this->op_->nrows_local ( ) + 1];
            this->op_->ExtractCSR ( ia, &( this->col_.front ( ) ), &( this->val_.front ( ) ) );

            // compute row indices
            int row_no = 1; // row_no + i1 treated row (fortran style)
            int ind_row = 0; // index in _row
            int nz_row; // number of nonzeros in row

            for ( int i = 0; i<this->op_->nrows_local ( ); ++i )
            {
                nz_row = ia[row_no] - ia[row_no - 1];

                // insert nz_row times row number
                for ( int k = 0; k < nz_row; ++k )
                {
                    this->row_[ind_row] = row_no + this->op_->row_ownership_begin ( );
                    ++ind_row;
                }
                ++row_no; // next row
            }

            // shift column indices to fortran style
            for ( int i = 0; i<this->op_->nnz_local ( ); ++i )
                this->col_[i]++;

            delete[] ia;
        }

        /// template instantiation
        template class MumpsSolver<LADescriptorCoupledD, MumpsStructureD>;
        template class MumpsSolver<LADescriptorCoupledS, MumpsStructureS>;

    } // namespace la
} // namespace hiflow
