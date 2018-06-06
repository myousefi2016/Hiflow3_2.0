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

#ifndef GMG_H
#    define GMG_H

#    include "gmg_hierarchy.h"
#    include "gmg_restriction.h"
#    include "linear_solver/jacobi.h"
#    include "linear_solver/gpu-based/async_iter_gpu.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            /// Implements a hierarchy for smoothers. Memory mangement of the stored smoothers
            /// is the responsibility of the user.

            template<class LAD>
            class SmootherHierarchy : public AbstractHierarchy<LinearSolver<LAD>*>
            {
              public:
                typedef LinearSolver<LAD>* SmootherType;

                IMPORT_FROM_BASECLASS ( AbstractHierarchy<SmootherType>, IteratorFromCoarsest );
                IMPORT_FROM_BASECLASS ( AbstractHierarchy<SmootherType>, IteratorFromFinest );
                IMPORT_FROM_BASECLASS ( AbstractHierarchy<SmootherType>, ConstIteratorFromCoarsest );
                IMPORT_FROM_BASECLASS ( AbstractHierarchy<SmootherType>, ConstIteratorFromFinest );

                /// Construct hierarchy
                /// @param n_levels The number of smoothers of the hierarchy

                explicit SmootherHierarchy ( std::size_t n_levels )
                {
                    this->levels_ = std::vector<SmootherType>( n_levels, NULL );
                }

                virtual ~SmootherHierarchy ( )
                {
                }
            };

            /// A functor to be used in conjunction with GeometricMultriGrid::for_each_level()
            /// that sets up Jacobi smoothers on each level. This class also deals with
            /// memory management of the smoothers. It is therefore important that this
            /// object lives at least as long as the GeometricMultiGrid solver.

            template<class LAD>
            class SetupJacobiSmoother
            {
              public:

                TYPE_FROM_CLASS ( LAD, DataType );

                typedef boost::shared_ptr<Jacobi<LAD> > SmootherPtr;

                SetupJacobiSmoother ( )
                : solve_mode_ ( 3 ), global_iter_ ( 3 ), inner_iter_ ( 1 ), w_ ( 0.5 ), async_ ( false )
                {
                }

                SetupJacobiSmoother ( const int solve_mode,
                                      const int global_iter,
                                      const int inner_iter,
                                      const DataType w,
                                      const bool async )
                : solve_mode_ ( solve_mode ), global_iter_ ( global_iter ), inner_iter_ ( inner_iter ), w_ ( w ), async_ ( async )
                {
                }

                void set_parameters ( const int solve_mode,
                                      const int global_iter,
                                      const int inner_iter,
                                      const DataType w,
                                      const bool async );

                /// Sets up a Jacobi Smoother on the supplied level
                /// @param lvl The level on which a Jacobi smoother shall be constructed
                /// @param smoother A reference to the pointer in the smoother hierarchy
                /// where the new smoother will be stored

                template<class LevelType>
                void operator() ( LevelType& lvl, LinearSolver<LAD>*& smoother )
                {
                    std::stringstream s;
                    s << "Creating smoother with parameters:\n"
                            << " | solve_mode  = " << solve_mode_ << "\n"
                            << " | global_iter = " << global_iter_ << "\n"
                            << " | inner_iter  = " << inner_iter_ << "\n"
                            << " | w           = " << w_ << "\n"
                            << " | async       = " << ( async_ ? 1 : 0 );
                    LOG_INFO ( "SetupJacobiSmoother", s.str ( ) );
                    // Setup Jacobi Solver
                    SmootherPtr new_smoother ( new Jacobi<LAD>( ) );

                    new_smoother->Prepare ( *lvl.matrix ( ), *lvl.rhs ( ) );
                    new_smoother->SetSolveMode ( solve_mode_ );
                    new_smoother->SetNumIter ( global_iter_ );
                    new_smoother->SetInnerIter ( inner_iter_ );
                    new_smoother->SetDampingParameter ( w_ );
                    new_smoother->SetAsyncFlag ( async_ );

                    // Store new smoother
                    smoother = new_smoother.get ( );
                    smoothers_.push_back ( new_smoother );
                }

              private:
                std::vector<SmootherPtr> smoothers_;
                int solve_mode_;
                int global_iter_;
                int inner_iter_;
                DataType w_;
                bool async_;
            };

            template<class LAD>
            class SetupAsyncIterGPUSmoother
            {
              public:

                TYPE_FROM_CLASS ( LAD, DataType );

                typedef boost::shared_ptr<AsynchronousIterationGPU<LAD> > SmootherPtr;

                SetupAsyncIterGPUSmoother ( )
                : use_gpus_ ( 1 ), procs_per_node_ ( 1 ),
                solve_mode_ ( 3 ), global_iter_ ( 3 ), block_iter_ ( 1 ), inner_iter_ ( 5 ), w_ ( 0.5 ),
                cuda_block_size_ ( 128 )
                {
                }

                SetupAsyncIterGPUSmoother ( const int use_gpus,
                                            const int procs_per_node,
                                            const int solve_mode,
                                            const int global_iter,
                                            const int block_iter,
                                            const int inner_iter,
                                            const DataType w,
                                            const int cuda_block_size )
                : use_gpus_ ( use_gpus ), procs_per_node_ ( procs_per_node ),
                solve_mode_ ( solve_mode ), global_iter_ ( global_iter ), block_iter_ ( block_iter ),
                inner_iter_ ( inner_iter ), w_ ( w ), cuda_block_size_ ( cuda_block_size )
                {
                }

                void set_parameters ( const int use_gpus,
                                      const int procs_per_node,
                                      const int solve_mode,
                                      const int global_iter,
                                      const int block_iter,
                                      const int inner_iter,
                                      const DataType w,
                                      const int cuda_block_size );

                /// Sets up a Jacobi Smoother on the supplied level
                /// @param lvl The level on which a Jacobi smoother shall be constructed
                /// @param smoother A reference to the pointer in the smoother hierarchy
                /// where the new smoother will be stored

                template<class LevelType>
                void operator() ( LevelType& lvl, LinearSolver<LAD>*& smoother )
                {
                    std::stringstream s;
                    s << "Creating smoother with parameters:\n"
                            << " | solve_mode      = " << solve_mode_ << "\n"
                            << " | global_iter     = " << global_iter_ << "\n"
                            << " | block_iter      = " << block_iter_ << "\n"
                            << " | inner_iter      = " << inner_iter_ << "\n"
                            << " | w               = " << w_ << "\n"
                            << " | CUDA block size = " << cuda_block_size_;
                    LOG_INFO ( "AsynchronousIterationGPU", s.str ( ) );
                    // setup AsynchronousIterationGPU solver
                    SmootherPtr new_smoother ( new AsynchronousIterationGPU<LAD>( lvl.comm ( ) ) );

                    new_smoother->Prepare ( *lvl.matrix ( ),
                                            *lvl.rhs ( ),
                                            *lvl.sol ( ),
                                            false,
                                            use_gpus_,
                                            procs_per_node_ );

                    new_smoother->SetSolveMode ( solve_mode_ );
                    new_smoother->SetNumIter ( global_iter_, block_iter_, inner_iter_ );
                    new_smoother->SetDampingParameter ( w_ );
                    new_smoother->SetCudaBlockSize ( cuda_block_size_ );

                    // Store new smoother
                    smoother = new_smoother.get ( );
                    smoothers_.push_back ( new_smoother );
                }

              private:
                std::vector<SmootherPtr> smoothers_;
                int use_gpus_, procs_per_node_, solve_mode_;
                int global_iter_, block_iter_, inner_iter_;
                int cuda_block_size_;
                DataType w_;
            };

            template<class LAD>
            class SetupIndividualSmoother
            {
              public:

                TYPE_FROM_CLASS ( LAD, DataType );

                SetupIndividualSmoother ( const std::vector<int>& smoother_types )
                : lvl_idx_ ( 0 ), smoother_types_ ( smoother_types ),
                jacobi_ ( ),
                async_gpu_ ( )
                {
                }

                void SetJacobiParameters ( const int solve_mode,
                                           const int global_iter,
                                           const int inner_iter,
                                           const DataType w,
                                           const bool async );

                void SetAsyncIterGPUParameters ( const int use_gpus,
                                                 const int procs_per_node,
                                                 const int solve_mode,
                                                 const int global_iter,
                                                 const int block_iter,
                                                 const int inner_iter,
                                                 const DataType w,
                                                 const int cuda_block_size );

                template<class LevelType>
                void operator() ( LevelType& lvl, LinearSolver<LAD>*& smoother )
                {
                    assert ( lvl_idx_ < smoother_types_.size ( ) );
                    switch ( smoother_types_[lvl_idx_] )
                    {
                        case 0: // Jacobi
                            jacobi_ ( lvl, smoother );
                            break;

                        case 1: // AsynchronousIterationGPU
                            async_gpu_ ( lvl, smoother );
                            break;

                        default:
                            LOG_ERROR ( "unknown smoother type, using Jacobi" );
                            jacobi_ ( lvl, smoother );
                            break;
                    }
                    ++lvl_idx_;
                }

              private:

                int lvl_idx_;
                std::vector<int> smoother_types_;

                SetupJacobiSmoother<LAD> jacobi_;
                SetupAsyncIterGPUSmoother<LAD> async_gpu_;

            };

        } // namespace gmg

        /// Implementation of a geometric multigrid linear solver

        template<class LAD>
        class GeometricMultiGrid : public LinearSolver<LAD>
        {
          public:

            TYPE_FROM_CLASS ( LAD, DataType );
            TYPE_FROM_CLASS ( gmg::HierarchyTypes<LAD>, GMGHierarchy );
            TYPE_FROM_CLASS ( GMGHierarchy, LevelType );
            TYPE_FROM_CLASS ( LevelType, MatrixPtr );
            TYPE_FROM_CLASS ( LevelType, VectorPtr );
            IMPORT_FROM_BASECLASS ( LinearSolver<LAD>, OperatorType );
            IMPORT_FROM_BASECLASS ( LinearSolver<LAD>, VectorType );

            typedef boost::shared_ptr<GMGHierarchy> HierarchyPtr;

            /// Initializes the solver. Collective on the root communicator supplied
            /// to the communicator hierarchy generator.
            /// @param master_rank The rank of the master process relative to the
            /// root communicator.
            /// @param levels How many levels the hierarchy is supposed to have
            /// @param master_mesh The master mesh that is to be refined and distributed.
            /// This parameter only needs to be a valid pointer on the master process.
            /// @param fe_degrees Specifies the fe-degrees. The i-th entry of the
            /// vector is the polynomial degree of the i-th variable of the problem
            /// @param settings The settings that are to be used by the linear algebra
            /// routines
            /// @param comm_generator A pointer to a communicator hierarchy generator.

            GeometricMultiGrid ( unsigned master_rank,
                                 unsigned levels,
                                 const MeshPtr& master_mesh,
                                 const std::vector<int>& fe_degrees,
                                 const gmg::Settings& settings,
                                 gmg::communicator_hierarchy_generators::BasicGenerator* comm_generator )

            : LinearSolver<LAD>( ),
            hierarchy_ ( new GMGHierarchy ( master_rank,
                                            levels,
                                            master_mesh,
                                            fe_degrees,
                                            settings,
                                            comm_generator ) ),
            solver_ ( NULL ), smoother_hierarchy_ ( levels ), use_sleep_ ( false ),
            use_start_cycle_ ( false ), cycle_type_ ( V_CYCLE ), solve_mode_ ( 0 ),
            total_coarse_grid_iter_ ( 0 ), grid_transfer_on_device_ ( false )
            {
            }

            virtual ~GeometricMultiGrid ( )
            {
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

            void SetSolveMode ( const int solve_mode )
            {
                solve_mode_ = solve_mode;

                if ( ( solve_mode_ != 0 ) && ( solve_mode_ != 1 ) )
                {
                    LOG_ERROR ( "GeometricMultiGrid::SetSolveMode called with invalid solve mode "
                                << solve_mode_ << " (must be 0 [solver] or 1 [preconditioner])" );
                    LOG_ERROR ( "GeometricMultiGrid: setting solve mode to 0 [solver]" );
                    solve_mode_ = 0;
                }
            }

            void use_sleep ( const bool use )
            {
                use_sleep_ = use;
            }

            void use_restriction_matrix ( const bool use )
            {
                linear_restriction_.use_restriction_matrix ( use );
            }

            void use_interpolation_matrix ( const bool use )
            {
                linear_interpolation_.use_interpolation_matrix ( use );
            }

            void use_grid_transfer_on_device ( const bool use )
            {
                grid_transfer_on_device_ = use;
            }

            void build_restriction_matrix ( const bool use_transposed_of_interpolation_matrix )
            {
                for ( typename GMGHierarchy::IteratorFromFinest lvl = hierarchy_->begin_from_finest_level ( );
                      lvl != hierarchy_->end_at_coarsest_level ( ) - 1; ++lvl )
                {
                    if ( use_transposed_of_interpolation_matrix )
                    {
                        linear_restriction_.use_transposed_of_interpolation_matrix ( *lvl );
                    }
                    else
                    {
                        linear_restriction_.build_restriction_matrix ( *lvl );
                    }
                }
            }

            void build_interpolation_matrix ( const bool use_transposed_of_restriction_matrix )
            {
                for ( typename GMGHierarchy::IteratorFromFinest lvl = hierarchy_->begin_from_finest_level ( );
                      lvl != hierarchy_->end_at_coarsest_level ( ) - 1; ++lvl )
                {
                    if ( use_transposed_of_restriction_matrix )
                    {
                        linear_interpolation_.use_transposed_of_restriction_matrix ( *lvl );
                    }
                    else
                    {
                        linear_interpolation_.build_interpolation_matrix ( *lvl );
                    }
                }
            }

            void prepare_grid_transfer_on_device ( void )
            {
                typename GMGHierarchy::IteratorFromFinest lvl
                        = hierarchy_->begin_from_finest_level ( );
                typename gmg::SmootherHierarchy<LAD>::IteratorFromFinest smt
                        = smoother_hierarchy_.begin_from_finest_level ( );

                for (; lvl != hierarchy_->end_at_coarsest_level ( ) - 1; ++lvl, ++smt )
                {
                    lvl->prepare_grid_transfer_on_device ( *smt );
                }
                grid_transfer_on_device_ = true;
            }

            /// Assembles all vectors and matrices neccessary to perform the multigrid algorithms
            /// (In more detail: Assembles the right hand side on the finest level and the matrices
            /// of all levels of the hierarchy)
            /// @param vector_assembler The vector assembly function
            /// @param matrix_assembler The matrix assembly funcion.

            void assemble_system ( typename GlobalAssembler<DataType>::VectorAssemblyFunction vector_assembler,
                                   typename GlobalAssembler<DataType>::MatrixAssemblyFunction matrix_assembler )
            {
                hierarchy_->begin_from_finest_level ( )->assemble_system ( vector_assembler, matrix_assembler );

                for ( typename GMGHierarchy::IteratorFromFinest lvl = hierarchy_->begin_from_finest_level ( ) + 1;
                      lvl != hierarchy_->end_at_coarsest_level ( ); ++lvl )
                {
                    lvl->assemble_matrix ( matrix_assembler );

                    if ( lvl->is_scheduled_to_this_process ( ) )
                        lvl->rhs ( )->Zeros ( );
                }
            }

            /// @return A pointer to the solution vector on the finest grid

            VectorPtr sol ( ) const
            {
                return hierarchy_->begin_from_finest_level ( )->sol ( );
            }

            /// @return A pointer to the right hand side vector on the finest grid

            VectorPtr rhs ( ) const
            {
                return hierarchy_->begin_from_finest_level ( )->rhs ( );
            }

            /// @return A pointer to the matrix of the finest grid

            MatrixPtr matrix ( ) const
            {
                return hierarchy_->begin_from_finest_level ( )->matrix ( );
            }

            /// @return Accumulated number of coarse grid solver iterations.

            int get_total_coarse_grid_iter ( ) const
            {
                return total_coarse_grid_iter_;
            }

            enum CycleType
            {
                V_CYCLE,
                W_CYCLE
            };

            /// Set the mode of operation
            /// @param use_full_multigrid_start_cycle whether to start with a full multigrid
            /// iteration (default: no)
            /// @param type The cycle type that shall be used (V- or W-cycles. default: V-cycles)

            void set_mode ( bool use_full_multigrid_start_cycle, CycleType type )
            {
                use_start_cycle_ = use_full_multigrid_start_cycle;
                cycle_type_ = type;
                LOG_INFO ( "GeometricMultiGrid", "set mode: full MG = " << ( use_full_multigrid_start_cycle ? "yes" : "NO" )
                           << ", cycle type = " << ( cycle_type_ == V_CYCLE ? "V" : "W" ) );
            }

            /// Setup Operator.
            /// @param op The matrix on which the solver operates. \c op must be the matrix
            /// return by \c matrix()

            virtual void SetupOperator ( OperatorType& op )
            {
                assert ( matrix ( ).get ( ) == &op );
            }

            /// Solve the system. Requires a previous call to \c setup_solvers(), and, for any
            /// meaningful calculation, a call to \c assemble_system() and \c set_dirichlet_bc()
            /// @param b The right hand side vector. Must be the vector returned by \c rhs()
            /// @param x The solution vector. Must be the vector returned by \c sol()
            /// @return The linear solver state of the multigrid solver

            virtual LinearSolverState Solve ( const VectorType& b, VectorType* x )
            {
                this->iter_ = 0;
                total_coarse_grid_iter_ = 0;
                if ( rhs ( ).get ( ) != &b )
                {
                    rhs ( )->CopyFrom ( b );
                }
                if ( sol ( ).get ( ) != x )
                {
                    sol ( )->CopyFrom ( *x );
                }
                return Solve ( );
            }

            /// Solve the system. Requires a previous call to \c setup_solvers(), and, for any
            /// meaningful calculation, a call to \c assemble_system() and \c set_dirichlet_bc()
            /// @return The linear solver state of the multigrid solver

            virtual LinearSolverState Solve ( )
            {
                switch ( solve_mode_ )
                {
                    case 0:
                        return SolveNormal ( );
                        break;

                    case 1:
                        return Precondition ( );
                        break;

                    default:
                        LOG_ERROR ( "GeometricMultriGrid: unknown solve mode " << solve_mode_
                                    << ", must be 0 [solver] or 1 [preconditioner]" );
                        exit ( -1 );
                        break;
                }
                return kSolverError;
            }

            virtual LinearSolverState SolveNormal ( )
            {
                //     gmg::util::master_ostream mcout(std::cout, 0);

                if ( use_start_cycle_ )
                    full_multigrid ( );

                this->res_ = calculate_finest_level_residual_norm ( );
                LOG_INFO ( "GeometricMultiGrid", "starting with residual norm " << this->res_ );
                for ( this->iter_ = 0;
                      this->control ( ).Check ( this->iter_, this->res_ ) == IterateControl::kIterate;
                      ++( this->iter_ ) )
                {
                    if ( cycle_type_ == V_CYCLE )
                        this->v_cycle ( );
                    else this->w_cycle ( );

                    this->res_ = calculate_finest_level_residual_norm ( );
                    LOG_INFO ( "GeometricMultiGrid", "cycle " << this->iter_ + 1 << " ended with residual norm " << this->res_ );
                }
                IterateControl::State state = this->control ( ).status ( );

                if ( state == IterateControl::kSuccessAbsoluteTol ||
                     state == IterateControl::kSuccessRelativeTol )
                    return kSolverSuccess;
                else if ( state == IterateControl::kFailureMaxitsExceeded )
                    return kSolverExceeded;
                else return kSolverError;
            }

            virtual LinearSolverState Precondition ( )
            {
                if ( cycle_type_ == V_CYCLE )
                    this->v_cycle ( );
                else this->w_cycle ( );

                ++( this->iter_ );

                return kSolverSuccess;
            }

            /// Setup dirichlet boundary conditions
            /// @param eval The dirichlet evaluator that specifies the boundary conditions

            template<class DirichletEvaluator>
            void set_dirichlet_bc ( DirichletEvaluator& eval )
            {
                hierarchy_->begin_from_finest_level ( )->init_boundary_conditions ( eval, false );

                for ( typename GMGHierarchy::IteratorFromFinest lvl = hierarchy_->begin_from_finest_level ( ) + 1;
                      lvl != hierarchy_->end_at_coarsest_level ( ); ++lvl )
                {
                    lvl->init_boundary_conditions ( eval, true );
                }
            }

            /// Sets the coarse solver that shall be used by the multigrid algorithm
            /// @param solver The linear solver that will be used to solve the coarsest level.

            void setup_solver ( LinearSolver<LAD>* solver )
            {
                solver_ = solver;
            }

            /// @return The smoother hierarchy of the multigrid algorithm. It is the users
            /// responsibility to fill the hierarchy with smoothers.

            gmg::SmootherHierarchy<LAD>& get_smoother_hierarchy ( )
            {
                return smoother_hierarchy_;
            }

            /// @return The smoother hierarchy of the multigrid algorithm. It is the users
            /// responsibility to fill the hierarchy with smoothers.

            const gmg::SmootherHierarchy<LAD>& get_smoother_hierarchy ( ) const
            {
                return smoother_hierarchy_;
            }

            /// Perform a V cycle.
            /// @return The linear solver state

            LinearSolverState v_cycle ( )
            {
                return v_cycle ( hierarchy_->end_at_coarsest_level ( ) - 1 );
            }

            /// Perform a V cycle.
            /// @param stop_level An iterator to the level that represents the coarsest level,
            /// i.e. the level at which the V cycle is supposed to stop descending through the
            /// hierarchy.
            /// @return The linear solver state

            LinearSolverState v_cycle ( typename GMGHierarchy::IteratorFromFinest stop_level )
            {
                assert ( solver_ );

                return cycle_process_level ( hierarchy_->begin_from_finest_level ( ),
                                             hierarchy_->begin_from_finest_level ( ),
                                             stop_level,
                                             smoother_hierarchy_.begin_from_finest_level ( ),
                                             1 );
            }

            /// Perform a W cycle.
            /// @return The linear solver state

            LinearSolverState w_cycle ( )
            {
                return w_cycle ( hierarchy_->end_at_coarsest_level ( ) - 1 );
            }

            /// Perform a W cycle
            /// @param stop_level An iterator to the level that represents the coarsest level,
            /// i.e. the level at which the W cycle is supposed to stop descending through the
            /// hierarchy.
            /// @return The linear solver state

            LinearSolverState w_cycle ( typename GMGHierarchy::IteratorFromFinest stop_level )
            {
                assert ( solver_ );

                return cycle_process_level ( hierarchy_->begin_from_finest_level ( ),
                                             hierarchy_->begin_from_finest_level ( ),
                                             stop_level,
                                             smoother_hierarchy_.begin_from_finest_level ( ),
                                             2 );
            }

            /// Perform full multigrid
            /// @return The linear solver state

            LinearSolverState full_multigrid ( )
            {
                typename GMGHierarchy::IteratorFromCoarsest lvl = hierarchy_->begin_from_coarsest_level ( );
                typename GMGHierarchy::IteratorFromFinest lvl_reverse = hierarchy_->end_at_coarsest_level ( ) - 1;

                typename gmg::SmootherHierarchy<LAD>::IteratorFromFinest smoother
                        = smoother_hierarchy_.end_at_coarsest_level ( ) - 1;

                for (; lvl != hierarchy_->end_at_finest_level ( ); ++lvl )
                {
                    // Prolongate solution from the next coarser level
                    if ( lvl != hierarchy_->begin_from_coarsest_level ( ) )
                    {
                        if ( lvl->is_scheduled_to_this_process ( ) )
                        {
                            lvl->get_connection_to_next_coarser_grid ( )->transfer_coarse_sol_to_fine_sol ( );
                            linear_interpolation_ ( *lvl, *lvl->sol ( ) );
                        }
                    }

                    // V cycle starting from each level
                    cycle_process_level ( lvl_reverse, lvl_reverse, hierarchy_->end_at_coarsest_level ( ) - 1, smoother, 1 );

                    --lvl_reverse;
                    --smoother;
                }

                return kSolverSuccess;
            }

            /// Executes a function on each level of the hierarchy. The supplied function
            /// must take a reference to a \c LevelType and a reference
            /// to a pointer to a smoother of type LinearSolver<LAD> as arguments. This provides
            /// for example the means to conveniently setup all smoothers of the smoother hierarchy
            /// if they require access to the matrices of their respective level.
            /// @param f The function that will be executed on each level. It will need to
            /// match the signature \c f(GeometricMultigrid<LAD>::LevelType&, LinearSolver<LAD>*&).
            /// The \c LinearSolver will be the smoother for the \c LevelType.

            template<class Function>
            void for_each_level ( Function& f )
            {
                typename GMGHierarchy::IteratorFromFinest lvl_it = hierarchy_->begin_from_finest_level ( );
                typename gmg::SmootherHierarchy<LAD>::IteratorFromFinest smoother_it
                        = smoother_hierarchy_.begin_from_finest_level ( );

                for (; lvl_it != hierarchy_->end_at_coarsest_level ( ); ++lvl_it, ++smoother_it )
                {
                    assert ( smoother_it != smoother_hierarchy_.end_at_coarsest_level ( ) );
                    if ( lvl_it->is_scheduled_to_this_process ( ) )
                        f ( *lvl_it, *smoother_it );
                }
            }

            template<class Function>
            void for_each_but_coarsest_level ( Function& f )
            {
                typename GMGHierarchy::IteratorFromFinest lvl_it = hierarchy_->begin_from_finest_level ( );
                typename gmg::SmootherHierarchy<LAD>::IteratorFromFinest smoother_it
                        = smoother_hierarchy_.begin_from_finest_level ( );

                for (; lvl_it != hierarchy_->end_at_coarsest_level ( ) - 1; ++lvl_it, ++smoother_it )
                {
                    assert ( smoother_it != smoother_hierarchy_.end_at_coarsest_level ( ) - 1 );
                    if ( lvl_it->is_scheduled_to_this_process ( ) )
                        f ( *lvl_it, *smoother_it );
                }
            }

            /// @return The multilevel hierarchy

            const GMGHierarchy& get_multilevel_hierarchy ( ) const
            {
                assert ( hierarchy_ );
                return *hierarchy_;
            }

            /// @return The multilevel hierarchy

            GMGHierarchy& get_multilevel_hierarchy ( )
            {
                assert ( hierarchy_ );
                return *hierarchy_;
            }
          private:
            /// Calculate the residual vector on a level
            /// @param lvl The level on which the residual shall be calculated

            inline void calculate_residual ( typename GMGHierarchy::IteratorFromFinest lvl )
            {
                if ( lvl->is_scheduled_to_this_process ( ) )
                {
                    lvl->matrix ( )->VectorMult ( *lvl->sol ( ), lvl->res ( ).get ( ) );
                    lvl->res ( )->ScaleAdd ( *lvl->rhs ( ), -1.0 );
                }
            }

            /// Calculate the residual vector on the finest level

            inline void calculate_finest_level_residual ( )
            {
                typename GMGHierarchy::IteratorFromFinest finest_lvl = hierarchy_->begin_from_finest_level ( );

                calculate_residual ( finest_lvl );
            }

            /// Calculate the residual vector on the finest level and return its norm
            /// @return The norm of the residual vector

            inline DataType calculate_finest_level_residual_norm ( )
            {
                calculate_finest_level_residual ( );

                typename GMGHierarchy::IteratorFromFinest finest_lvl = hierarchy_->begin_from_finest_level ( );

                if ( finest_lvl->is_scheduled_to_this_process ( ) )
                    return finest_lvl->res ( )->Norm2 ( );
                else
                    return 0.0;
            }

            inline void visualize_status ( typename GMGHierarchy::IteratorFromFinest lvl,
                                           const std::string& filename )
            {
                lvl->visualize_vector ( lvl->sol ( ), filename + "_solution" );
                lvl->visualize_vector ( lvl->rhs ( ), filename + "_rhs" );
                lvl->visualize_vector ( lvl->res ( ), filename + "_residual" );
            }

            /// Processes one level of the recursive V or W cycle implementation.
            /// @param lvl The current level
            /// @param start_lvl The start level of the cycle
            /// @param stop_lvl The level at which the algorithm shall cease to descend through
            /// the hierarchy.
            /// @param smoother An iterator to the current smoother in the smoother hierarchy
            /// @param num_recursive_calls The number of recursive calls. 1 corresponds to a V cycle,
            /// 2 to a W cycle.
            LinearSolverState cycle_process_level ( typename GMGHierarchy::IteratorFromFinest lvl,
                                                    typename GMGHierarchy::IteratorFromFinest start_lvl,
                                                    typename GMGHierarchy::IteratorFromFinest stop_lvl,
                                                    typename gmg::SmootherHierarchy<LAD>::IteratorFromFinest smoother,
                                                    unsigned num_recursive_calls );

            HierarchyPtr hierarchy_;

            gmg::SmootherHierarchy<LAD> smoother_hierarchy_;
            LinearSolver<LAD>* solver_;

            gmg::LinearInterpolation<LAD> linear_interpolation_;
            gmg::LinearRestriction<LAD> linear_restriction_;

            bool use_start_cycle_;
            bool grid_transfer_on_device_;
            bool use_sleep_;
            CycleType cycle_type_;
            int solve_mode_;
            int total_coarse_grid_iter_;
        };

    } // namespace la
} // namespace hiflow

#endif
