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

/// \author Aksel Alpay, Martin Wlotzka

#include "gmg.h"

#include <signal.h>

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            template<class LAD>
            void SetupJacobiSmoother<LAD>::set_parameters ( const int solve_mode,
                                                            const int global_iter,
                                                            const int inner_iter,
                                                            const DataType w,
                                                            const bool async )
            {
                assert ( solve_mode > 1 );
                assert ( solve_mode < 4 );
                assert ( global_iter > 0 );
                assert ( inner_iter > 0 );
                assert ( w > 0.0 );

                solve_mode_ = solve_mode;
                global_iter_ = global_iter;
                inner_iter_ = inner_iter;
                w_ = w;
                async_ = async;
            }

            template<class LAD>
            void SetupAsyncIterGPUSmoother<LAD>::set_parameters ( const int use_gpus,
                                                                  const int procs_per_node,
                                                                  const int solve_mode,
                                                                  const int global_iter,
                                                                  const int block_iter,
                                                                  const int inner_iter,
                                                                  const DataType w,
                                                                  const int cuda_block_size )
            {
                assert ( use_gpus > 0 );
                assert ( procs_per_node > 0 );
                assert ( solve_mode > 1 );
                assert ( solve_mode < 4 );
                assert ( global_iter > 0 );
                assert ( block_iter > 0 );
                assert ( inner_iter > 0 );
                assert ( w > 0.0 );
                assert ( cuda_block_size > 0 );

                use_gpus_ = use_gpus;
                procs_per_node_ = procs_per_node;
                solve_mode_ = solve_mode;
                global_iter_ = global_iter;
                block_iter_ = block_iter;
                inner_iter_ = inner_iter;
                w_ = w;
                cuda_block_size_ = cuda_block_size;
            }

            template<class LAD>
            void SetupIndividualSmoother<LAD>::SetJacobiParameters ( const int solve_mode,
                                                                     const int global_iter,
                                                                     const int inner_iter,
                                                                     const DataType w,
                                                                     const bool async )
            {
                jacobi_.set_parameters ( solve_mode, global_iter, inner_iter, w, async );
            }

            template<class LAD>
            void SetupIndividualSmoother<LAD>::SetAsyncIterGPUParameters ( const int use_gpus,
                                                                           const int procs_per_node,
                                                                           const int solve_mode,
                                                                           const int global_iter,
                                                                           const int block_iter,
                                                                           const int inner_iter,
                                                                           const DataType w,
                                                                           const int cuda_block_size )
            {
                async_gpu_.set_parameters ( use_gpus, procs_per_node, solve_mode,
                                            global_iter, block_iter, inner_iter, w, cuda_block_size );
            }

            template class SmootherHierarchy<LADescriptorCoupledD>;
            template class SmootherHierarchy<LADescriptorCoupledS>;

            template class SetupJacobiSmoother<LADescriptorCoupledD>;
            template class SetupJacobiSmoother<LADescriptorCoupledS>;

            template class SetupAsyncIterGPUSmoother<LADescriptorCoupledD>;
            template class SetupAsyncIterGPUSmoother<LADescriptorCoupledS>;

            template class SetupIndividualSmoother<LADescriptorCoupledD>;
            template class SetupIndividualSmoother<LADescriptorCoupledS>;

        } // namespace gmg

        template<class LAD>
        LinearSolverState GeometricMultiGrid<LAD>::cycle_process_level ( typename GMGHierarchy::IteratorFromFinest lvl,
                                                                         typename GMGHierarchy::IteratorFromFinest start_lvl,
                                                                         typename GMGHierarchy::IteratorFromFinest stop_lvl,
                                                                         typename gmg::SmootherHierarchy<LAD>::IteratorFromFinest smoother,
                                                                         unsigned num_recursive_calls )
        {
            //   gmg::util::master_ostream mcout(std::cout,0);
            //   unsigned level_idx = hierarchy_->get_levels() - 1 - std::distance(start_lvl, lvl);
            //
            //   std::string position_string = "iter "
            //           + boost::lexical_cast<std::string>(this->iter_)
            //           + " lvl " + boost::lexical_cast<std::string>(level_idx);

            if ( lvl->is_scheduled_to_this_process ( ) )
            {
                if ( lvl != start_lvl )
                    lvl->sol ( )->Zeros ( );

                if ( lvl != stop_lvl )
                {
                    assert ( smoother != smoother_hierarchy_.end_at_coarsest_level ( ) );

                    unsigned num_iter = num_recursive_calls;
                    if ( lvl == start_lvl )
                        num_iter = 1;

                    for ( unsigned i = 0; i != num_iter; ++i )
                    {
                        // Pre-smoothing
                        if ( *smoother )
                        {
                            ( *smoother )->SetupOperator ( *lvl->matrix ( ) );
                            ( *smoother )->Solve ( *lvl->rhs ( ), lvl->sol ( ).get ( ) );
                        }

                        if ( grid_transfer_on_device_ )
                        {
                            // start ghost updates
                            lvl->sol ( )->SendBorder ( );
                            lvl->sol ( )->ReceiveGhost ( );

                            // copy sol and rhs interior parts to device
                            lvl->CopySolToDevice ( );
                            lvl->CopyRhsToDevice ( );

                            // sol ghost update
                            lvl->sol ( )->WaitForRecv ( );
                            lvl->CopySolGhostToDevice ( );

                            // call kernel to compute res = rhs - A * sol
                            lvl->ComputeResidualOnDevice ( );

                            // copy residual to host for ghost update
                            lvl->CopyResToHost ( );
                            lvl->res ( )->SendBorder ( );
                            lvl->res ( )->ReceiveGhost ( );

                            // copy res to tmp on device
                            lvl->CopyResToTmp ( );

                            // res ghost update
                            lvl->res ( )->WaitForRecv ( );
                            lvl->CopyResGhostToDevice ( );

                            // call kernel to compute res = R_diag * tmp + R_offdiag * res_ghost
                            lvl->ComputeRestrictionOnDevice ( );

                            // copy device tmp to host res
                            lvl->CopyResToHost ( );

                            // end sends
                            lvl->sol ( )->WaitForSend ( );
                            lvl->res ( )->WaitForSend ( );
                        }
                        else
                        {
                            // Calculate residual
                            calculate_residual ( lvl );

                            // Restrict residual
                            linear_restriction_ ( *lvl, *lvl->res ( ), lvl->dirichlet_dofs ( ) );
                        }

                        lvl->get_connection_to_next_coarser_grid ( )->transfer_fine_res_to_coarse_rhs ( );

                        cycle_process_level ( lvl + 1, start_lvl, stop_lvl, smoother + 1, num_recursive_calls );

                        if ( use_sleep_ && ( lvl->partial_rank ( ) == 0 ) )
                        {
                            // send signal SIGUSR1 to sleeping procs
                            for ( int p = 1; p != lvl->pids ( ).size ( ); ++p )
                            {
                                kill ( lvl->pids ( )[p], SIGUSR1 );
                            }
                        }

                        // Prolongate the coarse solution
                        lvl->get_connection_to_next_coarser_grid ( )->transfer_coarse_sol_to_fine_res ( );

                        if ( grid_transfer_on_device_ )
                        {
                            lvl->res ( )->SendBorder ( );
                            lvl->res ( )->ReceiveGhost ( );

                            lvl->CopyResToDevice ( );

                            lvl->res ( )->WaitForRecv ( );
                            lvl->CopyResGhostToDevice ( );

                            lvl->ComputeUpdatedSolutionOnDevice ( );

                            lvl->CopySolToHost ( );
                            lvl->res ( )->WaitForSend ( );
                        }
                        else
                        {
                            // prolongate
                            linear_interpolation_ ( *lvl, *lvl->res ( ) );

                            // add correction
                            lvl->sol ( )->Axpy ( *lvl->res ( ), 1.0 );
                        }

                        // Post-smoothing
                        if ( *smoother )
                        {
                            ( *smoother )->SetupOperator ( *lvl->matrix ( ) );
                            ( *smoother )->Solve ( *lvl->rhs ( ), lvl->sol ( ).get ( ) );
                        }
                    }
                }
                else
                {
                    solver_->SetupOperator ( *lvl->matrix ( ) );
                    solver_->Solve ( *lvl->rhs ( ), lvl->sol ( ).get ( ) );
                    total_coarse_grid_iter_ += solver_->iter ( );
                }
            }
            else if ( use_sleep_ )
            {
                pause ( );
            }

            return kSolverSuccess;
        }

        template class GeometricMultiGrid<LADescriptorCoupledD>;
        template class GeometricMultiGrid<LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow
