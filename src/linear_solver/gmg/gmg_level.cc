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

#include <signal.h>

#ifndef NDEBUG
#    include <fstream>
#endif

#include "gmg_level.h"
#include "mesh/communication.h" // for SharedVertexTable

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            // function to be registered as signal handler for SIGUSR1, does nothing;
            // default signal handler might terminate the process, but it shall just wake up

            void do_nothing ( int unused )
            {
            }

            template<class LAD>
            BasicLevel<LAD>::BasicLevel ( MPI_Comm global_comm,
                                          MPI_Comm partial_comm,
                                          const MeshPtr& local_mesh,
                                          int master_rank,
                                          const std::vector<int>& fe_degrees,
                                          StandardGlobalAssembler<DataType>* global_asm,
                                          const Settings& settings )
            : global_comm_ ( global_comm ), partial_comm_ ( partial_comm ),
            global_size_ ( 0 ), global_rank_ ( -1 ), partial_rank_ ( -1 ), global_asm_ ( global_asm ),
            master_rank_ ( master_rank ), fe_degrees_ ( fe_degrees ), settings_ ( settings )
            {
                if ( global_comm_ != MPI_COMM_NULL )
                {
                    assert ( partial_comm_ != MPI_COMM_NULL );

                    MPI_Comm_rank ( global_comm_, &global_rank_ );
                    MPI_Comm_rank ( partial_comm_, &partial_rank_ );
                    MPI_Comm_size ( global_comm_, &global_size_ );

                    LOG_INFO ( "BasicLevel", "Creating level with " << global_size_ << " processes." );

                    SharedVertexTable shared_verts;
                    mesh_ = compute_ghost_cells ( *local_mesh, global_comm_, shared_verts );

                    init_la ( );

                    // share pid of all procs of the partial comm
                    int partial_size = 0;
                    MPI_Comm_size ( partial_comm_, &partial_size );
                    assert ( partial_size > 0 );
                    int my_pid = static_cast < int > ( getpid ( ) );
                    pids_ = std::vector<int>( partial_size, 0 );
                    MPI_Allgather ( &my_pid, 1, MPI_INT,
                                    util::raw_array ( pids_ ), 1, MPI_INT,
                                    partial_comm_ );
                }

                // register do_nothing function for signal SIGUSR1
                signal ( SIGUSR1, do_nothing );
            }

            template<class LAD>
            void BasicLevel<LAD>::assemble_system ( typename GlobalAssemblerType::VectorAssemblyFunction v_asm,
                                                    typename GlobalAssemblerType::MatrixAssemblyFunction m_asm )
            {
                if ( is_scheduled_to_this_process ( ) )
                {
                    assemble_matrix ( m_asm );

                    global_asm_->assemble_vector ( *space_, v_asm, *rhs_ );
                }
            }

            template<class LAD>
            void BasicLevel<LAD>::assemble_matrix ( typename GlobalAssemblerType::MatrixAssemblyFunction m_asm )
            {
                if ( is_scheduled_to_this_process ( ) )
                {
                    global_asm_->assemble_matrix ( *space_, m_asm, *matrix_ );
                    matrix_->Compress ( );
                }
            }

            template<class LAD>
            void BasicLevel<LAD>::visualize_vector ( const VectorPtr v, const std::string& filename ) const
            {
                if ( is_scheduled_to_this_process ( ) )
                {
                    v->Update ( );

                    int max_fe_degree = *std::max_element ( fe_degrees_.begin ( ),
                                                            fe_degrees_.end ( ) );

                    int num_intervals = max_fe_degree;
                    ParallelCellVisualization<DataType> visu ( *space_, num_intervals,
                                                               global_comm_, master_rank_ );

                    std::vector<DataType> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
                    std::vector<DataType> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
                    std::vector<DataType> material_number ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

                    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
                          it != mesh_->end ( mesh_->tdim ( ) );
                          ++it )
                    {
                        int temp1, temp2;
                        //       mesh_->get_attribute_value("_remote_index_", mesh_->tdim(),
                        //                                   it->index(),
                        //                                   &temp1);
                        mesh_->get_attribute_value ( "_sub_domain_", mesh_->tdim ( ),
                                                     it->index ( ),
                                                     &temp2 );
                        remote_index.at ( it->index ( ) ) = temp1;
                        sub_domain.at ( it->index ( ) ) = temp2;
                        material_number.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
                    }

                    for ( int var = 0; var < space_->get_nb_var ( ); ++var )
                    {
                        std::string identifier = "var";
                        identifier += boost::lexical_cast<std::string>( var );

                        visu.visualize ( EvalFeFunction<LAD>( *space_, *v, var ), identifier );
                    }

                    //     visu.visualize_cell_data(remote_index, "_remote_index_");
                    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
                    visu.visualize_cell_data ( material_number, "Material Id" );
                    visu.write ( filename );
                }
            }

            template<class LAD>
            typename BasicLevel<LAD>::VectorPtr BasicLevel<LAD>::create_vector ( void ) const
            {
                if ( is_scheduled_to_this_process ( ) )
                {
                    VectorPtr vec = VectorPtr ( new VectorType ( ) );
                    vec->Init ( global_comm_, *couplings_, settings_.platform, settings_.implementation );
                    vec->InitStructure ( );
                    vec->Zeros ( );
                    return vec;
                }
                else return VectorPtr ( );
            }

            template<class LAD>
            typename BasicLevel<LAD>::MatrixPtr BasicLevel<LAD>::create_matrix ( void ) const
            {
                if ( is_scheduled_to_this_process ( ) )
                {
                    MatrixPtr mat = MatrixPtr ( new MatrixType ( ) );
                    mat->Init ( global_comm_, *couplings_,
                                settings_.platform, settings_.implementation, settings_.matrix_format );

                    SparsityStructure sparsity;
                    global_asm_->compute_sparsity_structure ( *space_, sparsity );

                    mat->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                                         vec2ptr ( sparsity.diagonal_cols ),
                                         sparsity.diagonal_rows.size ( ),
                                         vec2ptr ( sparsity.off_diagonal_rows ),
                                         vec2ptr ( sparsity.off_diagonal_cols ),
                                         sparsity.off_diagonal_rows.size ( ) );
                    mat->Zeros ( );
                    return mat;
                }
                else return MatrixPtr ( );
            }

            template<class LAD>
            typename BasicLevel<LAD>::MatrixPtr BasicLevel<LAD>::create_matrix ( const LaCouplings& cp, const SparsityStructure& sp ) const
            {
                if ( is_scheduled_to_this_process ( ) )
                {
                    MatrixPtr mat = MatrixPtr ( new MatrixType ( ) );
                    mat->Init ( global_comm_, cp,
                                settings_.platform, settings_.implementation, settings_.matrix_format );
                    mat->InitStructure ( vec2ptr ( sp.diagonal_rows ),
                                         vec2ptr ( sp.diagonal_cols ),
                                         sp.diagonal_rows.size ( ),
                                         vec2ptr ( sp.off_diagonal_rows ),
                                         vec2ptr ( sp.off_diagonal_cols ),
                                         sp.off_diagonal_rows.size ( ) );
                    mat->Zeros ( );
                    return mat;
                }
                else return MatrixPtr ( );
            }

            template<class LAD>
            void BasicLevel<LAD>::create_restriction_matrix ( void )
            {
                restriction_matrix_ = create_matrix ( );
                if ( !tmp_ )
                    tmp_ = create_vector ( );
            }

            template<class LAD>
            void BasicLevel<LAD>::create_interpolation_matrix ( void )
            {
                interpolation_matrix_ = create_matrix ( );
                if ( !tmp_ )
                    tmp_ = create_vector ( );
            }

            template<class LAD>
            void BasicLevel<LAD>::init_la ( )
            {
                space_ = SpacePtr ( new VectorSpace<DataType>( global_comm_ ) );
                couplings_ = CouplingsPtr ( new Couplings<DataType>( ) );

                space_->Init ( fe_degrees_, *mesh_ );

                couplings_->Init ( global_comm_, space_->dof ( ) );

                SparsityStructure sparsity;
                global_asm_->compute_sparsity_structure ( *space_, sparsity );

                couplings_->InitializeCouplings ( sparsity.off_diagonal_rows,
                                                  sparsity.off_diagonal_cols );

                matrix_ = MatrixPtr ( new MatrixType ( ) );

                matrix_->Init ( global_comm_, *couplings_, settings_.platform, settings_.implementation,
                                settings_.matrix_format );

                matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                                         vec2ptr ( sparsity.diagonal_cols ),
                                         sparsity.diagonal_rows.size ( ),
                                         vec2ptr ( sparsity.off_diagonal_rows ),
                                         vec2ptr ( sparsity.off_diagonal_cols ),
                                         sparsity.off_diagonal_rows.size ( ) );

                matrix_->Zeros ( );

                rhs_ = create_vector ( );
                sol_ = create_vector ( );
            }

            // template<class LAD, class ConnectionType>
            // void GMGLevel<LAD, ConnectionType>::prepare_interpolation_on_device(void)
            // {
            //   LOG_ERROR("GMGLevel::prepare_interpolation_on_device() not yet implemented!");
            //   exit(-1);
            // }

            template class BasicLevel<LADescriptorCoupledD>;
            template class BasicLevel<LADescriptorCoupledS>;

            template<class ConnectionType> class ConnectedLevel<LADescriptorCoupledD, ConnectionType>;
            template<class ConnectionType> class ConnectedLevel<LADescriptorCoupledS, ConnectionType>;

            template<class ConnectionType> class GMGLevel<LADescriptorCoupledD, ConnectionType>;
            template<class ConnectionType> class GMGLevel<LADescriptorCoupledS, ConnectionType>;

        } // namespace gmg
    } // namespace la
} // namespace hiflow
