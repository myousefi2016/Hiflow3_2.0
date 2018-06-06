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

#ifndef GMG_LEVEL_H
#    define GMG_LEVEL_H

#    include <time.h>

#    include <string>
#    include <vector>
#    include <mpi.h>
#    include "boost/lexical_cast.hpp"
#    include "boost/shared_ptr.hpp"
#    include "boost/unordered_map.hpp"
#    include <stdexcept>
#    include <cassert>

#    include "mesh/mesh_tools.h"
#    include "assembly/assembly.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_algebra/lmp/lmatrix_csr_cpu.h" // for CPU_CSR_lMatrix
#    include "linear_algebra/lmp/lvector_cpu.h" // for CPU_lVector
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/gpu-based/async_iter_gpu.h"
#    include "space/cell_visualization.h"
#    include "tools/mpi_tools.h" // for mpi_data_type<DataType>::get_type()

#    include "util.h"

#    ifdef WITH_CUDA
#        include "linear_algebra/lmp/cuda/cuda_gmg.h"
#        include "linear_algebra/lmp/cuda/cuda_gmg.h"
#    endif

using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {
            /// Contains basic linear algebra settings

            struct Settings
            {
                PLATFORM platform;
                IMPLEMENTATION implementation;
                MATRIX_FORMAT matrix_format;
            };

            /// This class represents a single level within the MultiLevelHierarchy.
            /// It initializes a mesh, a vectorspace, matrices and vectors
            /// and provides functions to assemble and solve the system on this level.
            /// Note that a BasicSingleLevel object is rather a handle than an actual,
            /// copyable object. Although the object can be copied, it will still
            /// use the same pointers to the underlying, non-copying HiFlow data structures.
            /// A deep copy of a BasicSingleLevel object is not possible.

            template<class LAD>
            class BasicLevel
            {
              public:

                TYPE_FROM_CLASS ( LAD, MatrixType );
                TYPE_FROM_CLASS ( LAD, VectorType );
                TYPE_FROM_CLASS ( LAD, DataType );

                typedef boost::shared_ptr<MatrixType> MatrixPtr;
                typedef boost::shared_ptr<VectorType> VectorPtr;
                typedef boost::shared_ptr<VectorSpace<DataType> > SpacePtr;
                typedef boost::shared_ptr<Couplings<DataType> > CouplingsPtr;

                typedef StandardGlobalAssembler<DataType> GlobalAssemblerType;

                /// Initializes this level. Collective on the supplied communicator.
                /// @param global_comm The communicator that shall be used by this level
                /// @param master_mesh The master mesh that is to be distributed. It
                /// only needs to be a valid mesh pointer on the master rank. The
                /// BasicSingleLevel will then take care of its partitioning
                /// and distribution.
                /// This will be the coarsest mesh of the hierarchy. The MultiLevelHierarchy
                /// will obtain meshes for other levels of the hierarchy by refining this
                /// mesh.
                /// @param fe_degrees Specifies the fe-degrees. The i-th entry of the
                /// vector is the polynomial degree of the i-th variable of the problem
                /// @param global_asm A pointer to a global assembler object
                /// @param Settings An instance of the settings struct that provides
                /// information about the linear algebra routines and data structures
                /// to be used.

                BasicLevel ( MPI_Comm global_comm,
                             MPI_Comm partial_comm,
                             const MeshPtr& local_mesh,
                             int master_rank,
                             const std::vector<int>& fe_degrees,
                             StandardGlobalAssembler<DataType>* global_asm,
                             const Settings& settings );

                virtual ~BasicLevel ( )
                {
                }

                /// @return The process rank of the calling process in the communicator
                /// of this level

                int rank ( void ) const
                {
                    return global_rank_;
                }

                int partial_rank ( void ) const
                {
                    return partial_rank_;
                }

                /// @return The communicator that is used by this level, or MPI_COMM_NULL
                /// if the calling process is contained in the communicator for this level

                MPI_Comm comm ( void ) const
                {
                    return global_comm_;
                }

                MPI_Comm partial_comm ( void ) const
                {
                    return partial_comm_;
                }

                /// @return A pointer to the mesh used by this level, or a NULL pointer
                /// if the calling process is not contained in the communicator for this level.

                MeshPtr mesh ( void ) const
                {
                    return mesh_;
                }

                /// @return A pointer to the vector space used by this level, or a NULL pointer
                /// if the calling process is not contained in the communicator for this level.

                SpacePtr space ( void ) const
                {
                    return space_;
                }

                CouplingsPtr couplings ( void ) const
                {
                    return couplings_;
                }

                /// @return The solution vector of this level, or a NULL pointer
                /// if the calling process is not contained in the communicator for this level.

                VectorPtr sol ( void ) const
                {
                    return sol_;
                }

                /// @return The right hand side vector of this level, or a NULL pointer
                /// if the calling process is not contained in the communicator for this level.

                VectorPtr rhs ( void ) const
                {
                    return rhs_;
                }

                /// @return The matrix of this level, or a NULL pointer if the calling process
                /// is not contained in the communicator for this level.

                MatrixPtr matrix ( void ) const
                {
                    return matrix_;
                }

                MatrixPtr restriction_matrix ( void ) const
                {
                    return restriction_matrix_;
                }

                MatrixPtr interpolation_matrix ( void ) const
                {
                    return interpolation_matrix_;
                }

                VectorPtr tmp_vector ( void ) const
                {
                    return tmp_;
                }

                /// Assembles the Matrix and right hand side vector of this level.
                /// Collective on this level's communicator.
                /// @param v_asm The vector assembling function
                /// @param m_asm The matrix assembling function
                void assemble_system ( typename GlobalAssemblerType::VectorAssemblyFunction v_asm,
                                       typename GlobalAssemblerType::MatrixAssemblyFunction m_asm );

                /// Assembles only the matrix of this level. Collective on this level's
                /// communicator
                /// @param m_asm The matrix assembling function
                void assemble_matrix ( typename GlobalAssemblerType::MatrixAssemblyFunction m_asm );

                /// @return whether the calling process is contained in this level's
                /// communicator, i.e. whether it has been scheduled to do work on this
                /// level.

                inline bool is_scheduled_to_this_process ( ) const
                {
                    return global_comm_ != MPI_COMM_NULL;
                }

                /// Solves the system on this level. Collective on this level's communicator.
                /// @param solver The linear solver that shall be used.

                void solve_system ( LinearSolver<LAD>& solver )
                {
                    if ( is_scheduled_to_this_process ( ) )
                    {
                        solver.SetupOperator ( *matrix_ );
                        solver.Solve ( *rhs_, sol_.get ( ) );
                    }
                }

                /// Writes a visualization of a vector of this level to a file. Collective
                /// on this level's communicator.
                /// The \c num_intervals parameter of the visualization (which determines
                /// how many cells in the visualization correspond to real cells) will
                /// be set to the maximum fe degree of all variables.
                /// @param v A pointer to a vector living on this level.
                /// @param filename The name of the file where the visualization is
                /// to be saved.
                void visualize_vector ( const VectorPtr v, const std::string& filename ) const;

                /// Writes a visualization of the solution vector to a file. Collective
                /// on this level's communicator.
                /// The \c num_intervals parameter of the visualization (which determines
                /// how many cells in the visualization correspond to real cells) will
                /// be set to the maximum fe degree of all variables.
                /// @param filename The name of the file where the visualization is
                /// to be saved.

                void visualize_solution ( const std::string& filename ) const
                {
                    visualize_vector ( sol ( ), filename );
                }

                /// Sets the solution to zero again.
                /// Collective on this level's communicator.

                void reset_solution ( void )
                {
                    if ( is_scheduled_to_this_process ( ) )
                        sol_->Zeros ( );
                }

                /// Sets the matrix, rhs vector and the solution to zero again.
                /// After a call to this function, assemble_system() must be executed
                /// to make the system solvable.
                /// Collective on this level's communicator.

                void reset_system ( void )
                {
                    if ( is_scheduled_to_this_process ( ) )
                    {
                        rhs_->Zeros ( );
                        matrix_->Zeros ( );
                        reset_solution ( );
                    }
                }

                /// Creates a new vector using the vectorspace and settings of this level.

                VectorPtr create_vector ( void ) const;

                MatrixPtr create_matrix ( void ) const;

                MatrixPtr create_matrix ( const LaCouplings& cp, const SparsityStructure& sp ) const;

                void create_restriction_matrix ( void );

                void create_interpolation_matrix ( void );

                /// Initializes Dirichlet boundary conditions
                /// @param eval The Dirichlet Evaluator
                /// @param correct_only_matrix whether dirichlet bc should only be applied to the matrix or
                /// to the right hand side as well.

                template<class DirichletEvaluator>
                void init_boundary_conditions ( DirichletEvaluator& eval, bool correct_only_matrix = false )
                {
                    if ( is_scheduled_to_this_process ( ) )
                    {
                        dirichlet_dofs_.clear ( );
                        dirichlet_values_.clear ( );

                        for ( int var = 0; var < space_->get_nb_var ( ); ++var )
                            compute_dirichlet_dofs_and_values ( eval, *space_, var, dirichlet_dofs_,
                                                                dirichlet_values_ );

                        if ( !dirichlet_dofs_.empty ( ) )
                        {
                            if ( !correct_only_matrix )
                            {
                                rhs_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                                                  vec2ptr ( dirichlet_values_ ) );
                                sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                                                  vec2ptr ( dirichlet_values_ ) );
                            }

                            // Correct Dirichlet dofs.
                            matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
                        }
                    }
                }

                const std::vector<int>& dirichlet_dofs ( void ) const
                {
                    return dirichlet_dofs_;
                }

                const std::vector<int>& pids ( void ) const
                {
                    return pids_;
                }

              protected:

                std::string class_name;

                /// Initializes the linear algebra part of this level
                void init_la ( void );

                /// Initilizes the mesh of this level. Collective on this level's communicator.
                /// The mesh will be partitioned using METIS if WITH_METIS has been defined.
                /// Otherwise, the naive partitioner of HiFlow will be used.
                /// @param master_mesh The master mesh that shall be distributed. Only
                /// needs to be a valid argument on the master process.
                void init_mesh ( const MeshPtr& master_mesh );

                MPI_Comm global_comm_;
                MPI_Comm partial_comm_;
                MeshPtr mesh_;
                SpacePtr space_;
                CouplingsPtr couplings_;
                VectorPtr rhs_, sol_, tmp_;
                MatrixPtr matrix_, restriction_matrix_, interpolation_matrix_;
                StandardGlobalAssembler<DataType>* global_asm_;

                std::vector<int> dirichlet_dofs_;
                std::vector<DataType> dirichlet_values_;

                std::vector<int> fe_degrees_;

                int master_rank_;

                Settings settings_;

                int global_size_;
                int global_rank_;
                int partial_rank_;

                std::vector<int> pids_;

            };

            /// A BasicLevel implementation that contains two pointers
            /// to InterGridConnection objects to the next coarser and next finer grids.

            template<class LAD, class ConnectionType>
            class ConnectedLevel : public BasicLevel<LAD>
            {
              public:

                TYPE_FROM_CLASS ( LAD, DataType );

                typedef ConnectionType Connection;
                /// Initializes this level without setting the pointers to the
                /// InterGridConnection objects.
                /// Collective on the supplied communicator.
                /// @param comm The communicator that shall be used by this level
                /// @param master_mesh The master mesh that is to be distributed. It
                /// only needs to be a valid mesh pointer on the master rank. The
                /// \c BasicSingleLevel will then take care of its partitioning
                /// and distribution.
                /// This will be the coarsest mesh of the hierarchy. The MultiLevelHierarchy
                /// will obtain meshes for other levels of the hierarchy by refining this
                /// mesh.
                /// @param fe_degrees Specifies the fe-degrees. The i-th entry of the
                /// vector is the polynomial degree of the i-th variable of the problem
                /// @param global_asm A pointer to a global assembler object
                /// @param Settings An instance of the settings struct that provides
                /// information about the linear algebra routines and data structures
                /// to be used.

                ConnectedLevel ( MPI_Comm global_comm,
                                 MPI_Comm partial_comm,
                                 const MeshPtr& master_mesh,
                                 unsigned master_rank,
                                 const std::vector<int>& fe_degrees,
                                 StandardGlobalAssembler<DataType>* global_asm,
                                 const Settings& settings )
                : BasicLevel<LAD>( global_comm, partial_comm, master_mesh, master_rank,
                fe_degrees, global_asm, settings )
                {
                }

                /// Initializes this level from a BasicSingleLevel object and pointers
                /// to InterGridConnection objects to the finer and coarser levels.
                /// @param lvl A BasicSingleLevel object from which the ConnectedSingleLevel
                /// will be initialized.
                /// Note that the object created by this constructor will not be
                /// independent of lvl, as the vectors, matrices, vectorspace and mesh
                /// pointers will be shared.
                /// These are non-copyable objects, therefore a deep copy is not possible
                /// (or at least not worth the effort)
                /// @param to_finer A pointer to an InterGridConnection object linking
                /// this level with the next finer level.
                /// @param to_coarser A pointer to an InterGridConnection object linking
                /// this level with the next coarser level.

                ConnectedLevel ( const BasicLevel<LAD>& lvl,
                                 const boost::shared_ptr<Connection>& to_finer,
                                 const boost::shared_ptr<Connection>& to_coarser )
                : BasicLevel<LAD>( lvl ), connection_to_finer_grid_ ( to_finer ),
                connection_to_coarser_grid_ ( to_coarser )
                {
                }

                virtual ~ConnectedLevel ( )
                {
                }

                /// @return A pointer to the InterGridConnection object connecting
                /// this level with the next finer grid.

                boost::shared_ptr<ConnectionType>
                get_connection_to_next_finer_grid ( ) const
                {
                    return connection_to_finer_grid_;
                }

                /// @return A pointer to the InterGridConnection object connecting
                /// this level with the next coarser grid.

                boost::shared_ptr<ConnectionType>
                get_connection_to_next_coarser_grid ( ) const
                {
                    return connection_to_coarser_grid_;
                }

                /// Sets the connections of this level to the next coarser and finer levels.
                /// @param to_coarser A pointer to an InterGridConnection object connecting
                /// this level with the next coarser level.
                /// @param to_finer A pointer to an InterGridConnection object connecting
                /// this level with the next finer level.

                void set_connections ( const boost::shared_ptr<ConnectionType>& to_coarser,
                                       const boost::shared_ptr<ConnectionType>& to_finer )
                {
                    connection_to_coarser_grid_ = to_coarser;
                    connection_to_finer_grid_ = to_finer;
                }

              private:
                boost::shared_ptr<ConnectionType> connection_to_finer_grid_;
                boost::shared_ptr<ConnectionType> connection_to_coarser_grid_;
            };

            /// A connected \c ConnectedLevel that additionally contains a vector to store
            /// the residual, as required by the geometric multigrid algorithm.

            template<class LAD, class ConnectionType>
            class GMGLevel : public ConnectedLevel<LAD, ConnectionType>
            {
              public:

                TYPE_FROM_CLASS ( LAD, DataType );

                typedef ConnectedLevel<LAD, ConnectionType> BaseType;
                IMPORT_FROM_BASECLASS ( BaseType, VectorPtr );

                //   using ConnectedLevel<LAD, ConnectionType>::create_vector;

                /// Initializes this level without setting the pointers to the
                /// InterGridConnection objects.
                /// Collective on the supplied communicator.
                /// @param global_comm The communicator that shall be used by this level
                /// @param master_mesh The master mesh that is to be distributed. It
                /// only needs to be a valid mesh pointer on the master rank. The
                /// \c BasicSingleLevel will then take care of its partitioning
                /// and distribution.
                /// This will be the coarsest mesh of the hierarchy. The MultiLevelHierarchy
                /// will obtain meshes for other levels of the hierarchy by refining this
                /// mesh.
                /// @param fe_degrees Specifies the fe-degrees. The i-th entry of the
                /// vector is the polynomial degree of the i-th variable of the problem
                /// @param global_asm A pointer to a global assembler object
                /// @param Settings An instance of the settings struct that provides
                /// information about the linear algebra routines and data structures
                /// to be used.

                GMGLevel ( MPI_Comm global_comm,
                           MPI_Comm partial_comm,
                           const MeshPtr& master_mesh,
                           unsigned master_rank,
                           const std::vector<int>& fe_degrees,
                           StandardGlobalAssembler<DataType>* global_asm,
                           const Settings& settings )
                : ConnectedLevel<LAD, ConnectionType>( global_comm, partial_comm, master_mesh, master_rank,
                fe_degrees, global_asm, settings )
#    ifdef WITH_CUDA
                ,
                R_diag_ ( 0 ), R_offdiag_ ( 0 ),
                A_diag_ ( 0 ), A_offdiag_ ( 0 ),
                P_diag_ ( 0 ), P_offdiag_ ( 0 ),
                lvec_sol_ ( 0 ), lvec_sol_ghost_ ( 0 ),
                lvec_rhs_ ( 0 ), lvec_rhs_ghost_ ( 0 ),
                lvec_res_ ( 0 ), lvec_res_ghost_ ( 0 ),
                dev_sol_ ( 0 ), dev_sol_ghost_ ( 0 ),
                dev_rhs_ ( 0 ), dev_rhs_ghost_ ( 0 ),
                dev_res_ ( 0 ), dev_res_ghost_ ( 0 ),
                dev_tmp_ ( 0 ),
                host_sol_ ( 0 ), host_sol_ghost_ ( 0 ),
                host_rhs_ ( 0 ), host_rhs_ghost_ ( 0 ),
                host_res_ ( 0 ), host_res_ghost_ ( 0 ),
                row_R_diag_ ( 0 ), row_R_offdiag_ ( 0 ),
                row_A_diag_ ( 0 ), row_A_offdiag_ ( 0 ),
                row_P_diag_ ( 0 ), row_P_offdiag_ ( 0 ),
                col_R_diag_ ( 0 ), col_R_offdiag_ ( 0 ),
                col_A_diag_ ( 0 ), col_A_offdiag_ ( 0 ),
                col_P_diag_ ( 0 ), col_P_offdiag_ ( 0 ),
                nnz_A_diag_ ( 0 ), nnz_A_offdiag_ ( 0 ),
                nnz_R_diag_ ( 0 ), nnz_R_offdiag_ ( 0 ),
                nnz_P_diag_ ( 0 ), nnz_P_offdiag_ ( 0 ),
                nrows_ ( 0 ), block_dim_ ( 128 ),
                grid_transfer_on_device_ ( false ),
                async_iter_gpu_ ( 0 )
#    endif
                {
                    res_ = this->create_vector ( );
                }

                virtual ~GMGLevel ( )
                {
                    clear ( );
                }

                /// @return The residual vector of this level

                VectorPtr res ( void ) const
                {
                    return res_;
                }

                /// Computes the restricted residual from the current solution,
                /// i.e. r = R(b-Ax).
                void prepare_grid_transfer_on_device ( LinearSolver<LAD> * smoother );

                void CopySolToDevice ( )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        if ( async_iter_gpu_ == 0 )
                        {
                            memcpy ( host_sol_, lvec_sol_->buffer, nrows_ * sizeof (DataType ) );
                            cudaSetDevice ( 0 );
                            cudaMemcpyAsync ( dev_sol_, host_sol_,
                                              nrows_ * sizeof (DataType ),
                                              cudaMemcpyHostToDevice, stream_ );
                        }
#    endif
                    }
                }

                void CopySolGhostToDevice ( void )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        if ( nnz_A_offdiag_ > 0 )
                        {
                            memcpy ( host_sol_ghost_, lvec_sol_ghost_->buffer, lvec_sol_ghost_->get_size ( ) * sizeof (DataType ) );
                            cudaSetDevice ( 0 );
                            cudaMemcpyAsync ( dev_sol_ghost_, host_sol_ghost_,
                                              lvec_sol_ghost_->get_size ( ) * sizeof (DataType ),
                                              cudaMemcpyHostToDevice, stream_ );
                        }
#    endif
                    }
                }

                void CopyRhsToDevice ( )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        if ( async_iter_gpu_ == 0 )
                        {
                            memcpy ( host_rhs_, lvec_rhs_->buffer, nrows_ * sizeof (DataType ) );
                            cudaSetDevice ( 0 );
                            cudaMemcpyAsync ( dev_rhs_, host_rhs_,
                                              nrows_ * sizeof (DataType ),
                                              cudaMemcpyHostToDevice, stream_ );
                        }
#    endif
                    }
                }

                void CopyResToDevice ( )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        memcpy ( host_res_, lvec_res_->buffer, nrows_ * sizeof (DataType ) );
                        cudaSetDevice ( 0 );
                        cudaMemcpyAsync ( dev_res_, host_res_,
                                          nrows_ * sizeof (DataType ),
                                          cudaMemcpyHostToDevice, stream_ );
#    endif
                    }
                }

                void CopyResGhostToDevice ( void )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        if ( ( nnz_P_offdiag_ > 0 ) || ( nnz_R_offdiag_ > 0 ) )
                        {
                            memcpy ( host_res_ghost_, lvec_res_ghost_->buffer, lvec_res_ghost_->get_size ( ) * sizeof (DataType ) );
                            cudaSetDevice ( 0 );
                            cudaMemcpyAsync ( dev_res_ghost_, host_res_ghost_,
                                              lvec_res_ghost_->get_size ( ) * sizeof (DataType ),
                                              cudaMemcpyHostToDevice, stream_ );
                        }
#    endif
                    }
                }

                void ComputeResidualOnDevice ( void )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        cudaSetDevice ( 0 );

                        if ( nnz_A_offdiag_ > 0 )
                        {
                            cuda_gmg_compute_residual ( A_diag_, col_A_diag_, row_A_diag_, nrows_,
                                                        A_offdiag_, col_A_offdiag_, row_A_offdiag_,
                                                        dev_sol_, dev_rhs_, dev_res_,
                                                        dev_sol_ghost_,
                                                        grid_dim_, block_dim_, stream_ );
                        }
                        else
                        {
                            cuda_gmg_compute_residual_no_offdiag ( A_diag_, col_A_diag_, row_A_diag_, nrows_,
                                                                   dev_sol_, dev_rhs_, dev_res_,
                                                                   grid_dim_, block_dim_, stream_ );
                        }
#    endif
                    }
                }

                void ComputeRestrictionOnDevice ( void )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        cudaSetDevice ( 0 );

                        if ( nnz_R_offdiag_ > 0 )
                        {
                            cuda_gmg_compute_restriction ( R_diag_, col_R_diag_, row_R_diag_, nrows_,
                                                           R_offdiag_, col_R_offdiag_, row_R_offdiag_,
                                                           dev_res_, dev_res_ghost_, dev_tmp_,
                                                           grid_dim_, block_dim_, stream_ );
                        }
                        else
                        {
                            cuda_gmg_compute_restriction_no_offdiag ( R_diag_, col_R_diag_, row_R_diag_, nrows_,
                                                                      dev_res_, dev_tmp_,
                                                                      grid_dim_, block_dim_, stream_ );
                        }
#    endif
                    }
                }

                void ComputeUpdatedSolutionOnDevice ( void )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        cudaSetDevice ( 0 );

                        if ( nnz_P_offdiag_ > 0 )
                        {
                            cuda_gmg_compute_updated_solution ( P_diag_, col_P_diag_, row_P_diag_, nrows_,
                                                                P_offdiag_, col_P_offdiag_, row_P_offdiag_,
                                                                dev_sol_, dev_res_, dev_res_ghost_,
                                                                grid_dim_, block_dim_, stream_ );
                        }
                        else
                        {
                            cuda_gmg_compute_updated_solution_no_offdiag
                                    ( P_diag_, col_P_diag_, row_P_diag_, nrows_,
                                      dev_sol_, dev_res_,
                                      grid_dim_, block_dim_, stream_ );
                        }
#    endif
                    }
                }

                void CopySolToHost ( void )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        cudaSetDevice ( 0 );
                        cudaMemcpyAsync ( host_sol_, dev_sol_,
                                          nrows_ * sizeof (DataType ),
                                          cudaMemcpyDeviceToHost, stream_ );
                        cudaStreamSynchronize ( stream_ );
                        memcpy ( lvec_sol_->buffer, host_sol_, nrows_ * sizeof (DataType ) );
#    endif
                    }
                }

                void CopyResToHost ( void )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        cudaSetDevice ( 0 );
                        cudaMemcpyAsync ( host_res_, dev_res_,
                                          nrows_ * sizeof (DataType ),
                                          cudaMemcpyDeviceToHost, stream_ );
                        cudaStreamSynchronize ( stream_ );
                        memcpy ( lvec_res_->buffer, host_res_, nrows_ * sizeof (DataType ) );
#    endif
                    }
                }

                void CopyResToTmp ( void )
                {
                    if ( this->is_scheduled_to_this_process ( ) )
                    {
                        assert ( grid_transfer_on_device_ );
#    ifdef WITH_CUDA
                        cudaSetDevice ( 0 );
                        cudaMemcpyAsync ( dev_tmp_, dev_res_,
                                          nrows_ * sizeof (DataType ),
                                          cudaMemcpyDeviceToDevice, stream_ );
#    endif
                    }
                }

                //   void prepare_interpolation_on_device(void);

                void clear ( void )
                {
#    ifdef WITH_CUDA
                    if ( A_diag_ ) cudaFree ( A_diag_ );
                    if ( A_offdiag_ ) cudaFree ( A_offdiag_ );
                    if ( R_diag_ ) cudaFree ( R_diag_ );
                    if ( R_offdiag_ ) cudaFree ( R_offdiag_ );
                    if ( P_diag_ ) cudaFree ( P_diag_ );
                    if ( P_offdiag_ ) cudaFree ( P_offdiag_ );
                    if ( async_iter_gpu_ || dev_sol_ ) cudaFree ( dev_sol_ );
                    if ( async_iter_gpu_ || host_sol_ ) cudaFreeHost ( host_sol_ );
                    if ( async_iter_gpu_ || dev_sol_ghost_ ) cudaFree ( dev_sol_ghost_ );
                    if ( async_iter_gpu_ || host_sol_ghost_ ) cudaFreeHost ( host_sol_ghost_ );
                    if ( async_iter_gpu_ || dev_rhs_ ) cudaFree ( dev_rhs_ );
                    if ( async_iter_gpu_ || host_rhs_ ) cudaFreeHost ( host_rhs_ );
                    if ( dev_rhs_ghost_ ) cudaFree ( dev_rhs_ghost_ );
                    if ( host_rhs_ghost_ ) cudaFreeHost ( host_rhs_ghost_ );
                    if ( dev_res_ ) cudaFree ( dev_res_ );
                    if ( host_res_ ) cudaFreeHost ( host_res_ );
                    if ( dev_res_ghost_ ) cudaFree ( dev_res_ghost_ );
                    if ( host_res_ghost_ ) cudaFreeHost ( host_res_ghost_ );
                    if ( dev_tmp_ ) cudaFree ( dev_tmp_ );
                    if ( row_A_diag_ ) cudaFree ( row_A_diag_ );
                    if ( row_A_offdiag_ ) cudaFree ( row_A_offdiag_ );
                    if ( row_R_diag_ ) cudaFree ( row_R_diag_ );
                    if ( row_R_offdiag_ ) cudaFree ( row_R_offdiag_ );
                    if ( row_P_diag_ ) cudaFree ( row_P_diag_ );
                    if ( row_P_offdiag_ ) cudaFree ( row_P_offdiag_ );
                    if ( col_A_diag_ ) cudaFree ( col_A_diag_ );
                    if ( col_A_offdiag_ ) cudaFree ( col_A_offdiag_ );
                    if ( col_R_diag_ ) cudaFree ( col_R_diag_ );
                    if ( col_R_offdiag_ ) cudaFree ( col_R_offdiag_ );
                    if ( col_P_diag_ ) cudaFree ( col_P_diag_ );
                    if ( col_P_offdiag_ ) cudaFree ( col_P_offdiag_ );
#    endif
                }

              private:
                VectorPtr res_;
                bool grid_transfer_on_device_;

#    ifdef WITH_CUDA
                DataType *A_diag_, *A_offdiag_, *R_diag_, *R_offdiag_, *P_diag_, *P_offdiag_;
                DataType *dev_sol_, *dev_sol_ghost_, *dev_rhs_, *dev_rhs_ghost_, *dev_res_, *dev_res_ghost_, *dev_tmp_;
                DataType *host_sol_, *host_sol_ghost_, *host_rhs_, *host_rhs_ghost_, *host_res_, *host_res_ghost_;
                int *row_A_diag_, *row_A_offdiag_, *col_A_diag_, *col_A_offdiag_;
                int *row_R_diag_, *row_R_offdiag_, *col_R_diag_, *col_R_offdiag_;
                int *row_P_diag_, *row_P_offdiag_, *col_P_diag_, *col_P_offdiag_;
                int nnz_A_diag_, nnz_A_offdiag_, nnz_R_diag_, nnz_R_offdiag_, nnz_P_diag_, nnz_P_offdiag_;
                int nrows_;
                int grid_dim_, block_dim_;

                CPU_lVector<DataType> *lvec_sol_, *lvec_sol_ghost_, *lvec_rhs_, *lvec_rhs_ghost_, *lvec_res_, *lvec_res_ghost_;
                AsynchronousIterationGPU<LAD> *async_iter_gpu_;
                cudaStream_t stream_;
#    endif
            };

        } // namespace gmg
    } // namespace la
} // namespace hiflow

#endif
