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

/// @author Martin Wlotzka

#ifndef HIFLOW_LINEARSOLVER_MULTIGPUBLOCKASYNCITER_
#    define HIFLOW_LINEARSOLVER_MULTIGPUBLOCKASYNCITER_

#    include "linear_solver/linear_solver.h"
#    include "linear_algebra/lmp/lmatrix_csr_cpu.h"

#    ifdef WITH_CUDA
#        include <cuda_runtime_api.h>
#    endif

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        class AsynchronousIterationGPU : public LinearSolver<LAD>
        {
          public:

            typedef typename LAD::MatrixType MatrixType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            AsynchronousIterationGPU ( const MPI_Comm& comm );

            virtual ~AsynchronousIterationGPU ( );

            void Clear ( );

            void Prepare ( const MatrixType& A,
                           const VectorType& b,
                           VectorType& x,
                           const bool compute_residual_on_dev,
                           const int use_gpus,
                           const int procs_per_node );

            void SetNumIter ( const int n_global, const int n_block = 1, const int n_inner = 5 );
            void SetSolveMode ( const int mode );
            void SetDampingParameter ( const DataType w );
            void SetCudaBlockSize ( const int block_size );

            virtual void SetupOperator ( MatrixType& op )
            {
                A_ = &op;
            }

            /// Sets the relative tolerance.
            /// Needed by Inexact Newton Methods
            /// @param reltol relative tolerance of residual to converge

            void SetRelativeTolerance ( double reltol )
            {
                int maxits = this->control_.maxits ( );
                double atol = this->control_.absolute_tol ( );
                double dtol = this->control_.divergence_tol ( );
                this->control_.Init ( maxits, atol, reltol, dtol );
            }

            virtual LinearSolverState Solve ( const VectorType& b, VectorType* x );

            LinearSolverState SolveNormal ( const VectorType& b, VectorType* x );
            LinearSolverState SolveDamped ( const VectorType& b, VectorType* x );
            LinearSolverState SmoothNormal ( const VectorType& b, VectorType* x );
            LinearSolverState SmoothDamped ( const VectorType& b, VectorType* x );

            void CopyHostToDev_x ( );
            void CopyHostToDev_b_x ( );
            void CopyDevToHost_x ( );

            void Run ( VectorType *x );
            void RunDamped ( VectorType *x );

            void ComputeResidualOnDev ( );

            void BlockUpdate ( );
            void GlobalUpdate ( VectorType *x );

            int n_gpus ( void ) const
            {
                return n_gpus_;
            }

            DataType* dev_b_ptr ( const int dev )
            {
                return b_gpu[dev];
            }

            DataType* host_b_ptr ( const int dev )
            {
                return b_host[dev];
            }

            DataType* dev_x_ptr ( const int dev )
            {
                return x_gpu[dev];
            }

            DataType* host_x_ptr ( const int dev )
            {
                return x_host[dev];
            }

            DataType* dev_xg_ptr ( const int dev )
            {
                return xg_gpu[dev];
            }

            DataType* host_xg_ptr ( const int dev )
            {
                return xg_host[dev];
            }

          protected:

            DataType ComputeBlockSquaredResidualOnDev ( );

            MatrixType *A_;
            VectorType *res_vec_;

            const CPU_CSR_lMatrix<DataType> *lmat_A_diag_, *lmat_A_offdiag_;
            const CPU_lVector<DataType> *lvec_b_;
            CPU_lVector<DataType> *lvec_x_, *lvec_xg_;

            // host pinned memory
            DataType **b_host, **x_host, **xg_host;
            DataType *block_res_host;

            // GPU memory
            DataType **A_diag_gpu, **A_offdiag_gpu;
            DataType **b_gpu, **x_gpu, **xg_gpu;
            DataType **inv_D_gpu, **D_gpu, **block_res_gpu;

#    ifdef WITH_CUDA
            cudaStream_t *stream_;
#    endif

            int n_gpus_;
            int used_gpus_;
            int enable_hyperq_;
            int enable_p2p_;
            int grid_dim_;
            int block_dim_;
            int n_global_iter_;
            int n_block_iter_;
            int n_inner_iter_;
            int solve_mode_;
            int max_threads_per_block_;

            DataType w_;

            MPI_Comm comm_;
            int rank_;
            int nproc_;

            int *row_offs_; // row offsets (diagonal + offdiagonal matrix part + vectors)
            int *nnz_diag_;
            int **row_diag_;
            int **col_diag_;
            int *nnz_offdiag_;
            int **row_offdiag_;
            int **col_offdiag_;

            bool have_offdiag_;
            bool compute_residual_on_dev_;

          private:

            // prohibit copy constructor and assignment operator
            AsynchronousIterationGPU ( const AsynchronousIterationGPU<LAD>& other );
            AsynchronousIterationGPU<LAD>& operator= ( const AsynchronousIterationGPU<LAD>& other );
        };

    } // namespace la
} // namespace hiflow

#endif
