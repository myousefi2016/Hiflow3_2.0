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

/// \author Philipp Gerstner

#include "parallel_adaptive_jump_term_test.h"

const int PRINT_PROC = 1;
//#ifndef MESHES_DATADIR
#define MESHES_DATADIR "../examples/data/"
//#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class JumpTest
{
  public:

    JumpTest ( const std::string& param_filename,
               const std::string& path_mesh )
    : path_mesh ( path_mesh ),
    comm_ ( MPI_COMM_WORLD ),
    rank_ ( -1 ),
    num_partitions_ ( -1 ),
    is_done_ ( false ),
    refinement_level_ ( 0 ),
    sol_ ( new CoupledVector<Scalar> )
    {
        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );

#ifndef USE_MESH_P4EST
        // Application only works sequentially, not parallel!
        if ( num_partitions_ > 1 )
        {
            std::cerr
                    << "Test can only be run sequentially!"
                    << std::endl;
            // exit ( -1 );
        }
#endif
    }

    // Main algorithm

    void run ( )
    {
        // Construct / read in the initial mesh.
        build_initial_mesh ( );

        // local refinement
        adapt ( );

        // Initialize space and linear algebra.
        prepare_system ( );

        // assemble jump terms
        assemble_jump_terms ( );

        visualize ( );
    }

    ~JumpTest ( )
    {
        delete sol_;
    }

  private:
    // Member functions

    // Read and distribute mesh.
    std::string path_mesh;
    void build_initial_mesh ( );

    // Setup space, linear algebra, and compute Dirichlet values.
    void prepare_system ( );

    // Compute the matrix and rhs.
    void assemble_jump_terms ( );

    // Adapt the space (mesh and/or degree).
    void adapt ( );

    void visualize ( );

    // member variables
    // MPI communicator.
    MPI_Comm comm_;
    // Local process rank and number of processes.
    int rank_, num_partitions_;

    // Local mesh and mesh on master process.
    MeshPtr mesh_, master_mesh_;
    // Solution space.
    VectorSpace<double> space_;

    // Linear algebra couplings helper object.
    Couplings<double> couplings_;
    // Vectors for solution and load vector.
    CoupledVector<Scalar>* sol_;

    // Global assembler.
    // HpFem instead of StandardAssembler to allow different polynomial degrees.
    HpFemAssembler<double> global_asm_;
    DGGlobalAssembler<double> jump_term_asm_;

    // Flag for stopping adaptive loop.
    bool is_done_;
    // Current refinement level.
    int refinement_level_;

    // Vectors for error norms and error estimators
    std::vector<Scalar> rho_jump_;

}; // end class PoissonAdaptive

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    // Set default parameter file
    std::string param_filename ( "" );
    std::string path_mesh;
    // If set take parameter file specified on console
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }
    // If set take mesh following path specified on console
    if ( argc > 2 )
    {
        path_mesh = std::string ( argv[2] );
    }
    try
    {
        // Create log files for INFO and DEBUG output
        //    std::ofstream info_log("poisson_tutorial_info_log");
        LogKeeper::get_log ( "info" ).set_target ( &( std::cout ) );
        //std::ofstream debug_log("poisson_tutorial_debug_log");
        LogKeeper::get_log ( "debug" ).set_target ( &( std::cout ) );
        //     std::ofstream error_log ( "poisson_tutorial_error_log" );
        LogKeeper::get_log ( "error" ).set_target ( &( std::cout ) );

        // Create application object and run it
        JumpTest app ( param_filename, path_mesh );
#ifdef USE_MESH_P4EST
        if ( DIMENSION == 1 )
        {
            std::cout << "Combination of MeshP4est and DIMENSION==1 does not work" << std::endl;
            return 0;
        }
#endif
        app.run ( );

    }
    catch ( std::exception& e )
    {
        std::cerr << "\nProgram ended with uncaught exception.\n";
        std::cerr << e.what ( ) << "\n";
        return -1;
    }
    MPI_Finalize ( );
    return 0;
}

void JumpTest::build_initial_mesh ( )
{
    // Read in the mesh.
    //The mesh is chosen according to the dimension of the problem.
    std::string mesh_name;

    mesh_name = "parallel_adaptive_jump_term_test_mesh.inp";

    std::string mesh_filename;
    if ( path_mesh.empty ( ) )
    {
        mesh_filename = std::string ( DATADIR ) + mesh_name;
    }
    else
    {
        mesh_filename = path_mesh + mesh_name;
    }
    std::vector<MasterSlave> period ( 0, MasterSlave ( 0., 0., 0., 0 ) );

    //#ifdef USE_MESH_P4EST
    SharedVertexTable shared_verts;
    if ( rank_ == MASTER_RANK )
    {
#ifdef USE_MESH_P4EST  
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period, mesh::IMPL_P4EST );
#else
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period, mesh::IMPL_DBVIEW );
#endif
        // Refine the mesh until the initial refinement level is reached.
        const int initial_ref_lvl = 0;
        if ( initial_ref_lvl > 0 )
        {
            master_mesh_ = master_mesh_->refine_uniform_seq ( initial_ref_lvl );
        }
        refinement_level_ = initial_ref_lvl;
    }
    // Broadcast information from master to slaves.
    MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
    int uniform_ref_steps = 0;
#ifdef USE_MESH_P4EST
    mesh_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &uniform_ref_steps, mesh::IMPL_P4EST );
    mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_P4EST );
#else
    mesh_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &uniform_ref_steps, mesh::IMPL_DBVIEW );
    mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_DBVIEW );
#endif
    std::vector<int> local_index;
    for ( int l = 0; l < mesh_->num_entities ( DIMENSION ); ++l )
    {
        local_index.push_back ( l );
    }
    AttributePtr local_index_attr ( new IntAttribute ( local_index ) );
    mesh_->add_attribute ( "_local_index_", DIMENSION, local_index_attr );

    PVtkWriter writer ( comm_ );
    std::ostringstream name;
    name << "jump_term_mesh_" << refinement_level_ << ".pvtu";
    std::string output_file = name.str ( );

    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

void JumpTest::prepare_system ( )
{
    MeshPtr mesh_ptr = mesh_;

    if ( sol_ != NULL )
    {
        delete sol_;
    }

    sol_ = new CoupledVector<Scalar>;

    QuadratureSelection q_sel ( 3 );
    global_asm_.set_quadrature_selection_function ( q_sel );

    // Assign degrees to each element.
    const int fe_degree = 1;
    std::vector< int > degrees ( 1, fe_degree );

    // Initialize the VectorSpace object.
    space_.Clear ( );
    space_.Init ( degrees, *mesh_ptr );

    // Setup couplings object.
    couplings_.Clear ( );
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    sol_->Init ( comm_, couplings_, CPU, NAIVE );
    sol_->InitStructure ( );
    sol_->Zeros ( );

    LOG_DEBUG ( 2, "[" << rank_ << "] sparsity.diagonal_cols = " << string_from_range ( sparsity.diagonal_cols.begin ( ), sparsity.diagonal_cols.end ( ) ) );
    LOG_DEBUG ( 2, "[" << rank_ << "] sparsity.diagonal_rows = " << string_from_range ( sparsity.diagonal_rows.begin ( ), sparsity.diagonal_rows.end ( ) ) );
    LOG_DEBUG ( 2, "[" << rank_ << "] sparsity.off_diagonal_cols = " << string_from_range ( sparsity.off_diagonal_cols.begin ( ), sparsity.off_diagonal_cols.end ( ) ) );
    LOG_DEBUG ( 2, "[" << rank_ << "] sparsity.off_diagonal_rows = " << string_from_range ( sparsity.off_diagonal_rows.begin ( ), sparsity.off_diagonal_rows.end ( ) ) );
    LOG_DEBUG ( 2, "[" << rank_ << "] local dofs = " << sol_->size_local ( ) << ", ghost size = " << sol_->size_local_ghost ( ) );
}

void JumpTest::assemble_jump_terms ( )
{
    MeshPtr mesh_ptr = mesh_;

    // ***************************************************************
    // set dofs
    for ( mesh::EntityIterator it = mesh_ptr->begin ( DIMENSION ), end_it = mesh_ptr->end ( DIMENSION ); it != end_it; ++it )
    {
        std::vector<int> global_dof_ids;
        space_.GetDofIndices ( 0, *it, &global_dof_ids );
        int num_dofs = global_dof_ids.size ( );
        std::vector<double> values;
        values.resize ( num_dofs, 0. );

        std::vector< Coord > coords;
        space_.dof ( ).get_coord_on_cell ( 0, it->index ( ), coords );
        for ( int i = 0; i < num_dofs; i++ )
        {
            if ( rank_ == space_.dof ( ).owner_of_dof ( global_dof_ids.at ( i ) ) )
            {
                double x = coords[i][0];
                double y = coords[i][1];
                double val = 0.;
                if ( x <= 1. )
                {
                    val = x;
                }
                else
                {
                    val = 2. - x;
                }
                sol_->SetValues ( &global_dof_ids.at ( i ), 1, &val );

                LOG_DEBUG ( 2, "dof " << global_dof_ids.at ( i ) << " has coord " << x << ", " << y );
            }
        }
    }

    LOG_DEBUG ( 2, "[" << rank_ << "] Vector contains pp data before updatecouplings()? " << sol_->HasPpData ( ) );
    sol_->UpdateCouplings ( );
    LOG_DEBUG ( 2, "[" << rank_ << "] Vector contains pp data after updatecouplings()? " << sol_->HasPpData ( ) );

    if ( sol_->pp_data ( ).dof_ids.size ( ) > 0 )
    {
        LOG_DEBUG ( 2, "[" << rank_ << "] pp ids " << string_from_range ( sol_->pp_data ( ).dof_ids.begin ( ), sol_->pp_data ( ).dof_ids.end ( ) ) );
        LOG_DEBUG ( 2, "[" << rank_ << "] pp val " << string_from_range ( sol_->pp_data ( ).values.begin ( ), sol_->pp_data ( ).values.end ( ) ) );
    }
    else
    {
        LOG_DEBUG ( 2, "[" << rank_ << "] empty pp data !" );
    }
    for ( mesh::EntityIterator it = mesh_ptr->begin ( DIMENSION ), end_it = mesh_ptr->end ( DIMENSION ); it != end_it; ++it )
    {
        std::vector<int> global_dof_ids;
        space_.GetDofIndices ( 0, *it, &global_dof_ids );
        int num_dofs = global_dof_ids.size ( );
        std::vector<double> values;
        values.resize ( num_dofs, 0. );

        std::vector< Coord > coords;
        space_.dof ( ).get_coord_on_cell ( 0, it->index ( ), coords );
        for ( int i = 0; i < num_dofs; i++ )
        {
            double x = coords[i][0];
            double y = coords[i][1];
            double val = 0.;

            sol_->GetValues ( &global_dof_ids.at ( i ), 1, &val );

            double val_ref = 0.;
            if ( x <= 1. )
            {
                val_ref = x;
            }
            else
            {
                val_ref = 2. - x;
            }
            if ( val != val_ref )
            {
                LOG_DEBUG ( 2, "[" << rank_ << "]: coord = " << x << ", " << y << ",  u = " << val << " BUG !!!" );
            }
        }
    }

    // ***************************************************************
    // get number of interior facets and boundary facets
    InterfaceList if_list = InterfaceList::create ( mesh_ptr );

    int Number_Boundary = 0;
    for ( InterfaceList::const_iterator it = if_list.begin ( ),
          end_it = if_list.end ( );
          it != end_it;
          ++it )
    {

        int remote_index_master = -10;
        mesh_ptr->get_attribute_value
                ( "_remote_index_",
                  mesh_ptr->tdim ( ),
                  it->master_index ( ),
                  &remote_index_master
                  );

        const int num_slaves = it->num_slaves ( );
        if ( remote_index_master == -1 )
        {
            if ( num_slaves == 0 )
            {
                Number_Boundary += 1;
            }
        }
    }

    // ***************************************************************
    // rho_E : jump terms
    rho_jump_.clear ( );
    JumpTermAssembler rho_jump_int ( *sol_ );
    jump_term_asm_.assemble_interface_scalar_cells
            (
              space_,
              rho_jump_int,
              rho_jump_
              );

    // Create attribute with inner edge jump term for output.
    AttributePtr rho_jump_attr ( new DoubleAttribute ( rho_jump_ ) );
    mesh_ptr->add_attribute ( "rho_E", DIMENSION, rho_jump_attr );
    double total_rho_jump = std::accumulate
            (
              rho_jump_.begin ( ),
              rho_jump_.end ( ),
              0.
              );
    double global_rho_jump = 0.;
    MPI_Allreduce ( &total_rho_jump,
                    &global_rho_jump,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_
                    );
    int local_size_E = if_list.size ( ) - Number_Boundary;
    int global_size_E = 0;
    MPI_Allreduce ( &local_size_E,
                    &global_size_E,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_
                    );

    LOG_INFO ( "rho_E", "partition " << rank_ << ": size = " << local_size_E << ", local value = " << std::sqrt ( total_rho_jump ) );
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "rho_E", "global size = " << global_size_E << ", global value = " << std::sqrt ( global_rho_jump ) );
    }
}

void JumpTest::visualize ( )
{
    MeshPtr mesh_ptr = mesh_;

    sol_->UpdateCouplings ( );

    // Setup visualization object.
    const int num_intervals = 1;
    ParallelCellVisualization<double> visu
            ( space_,
              num_intervals,
              comm_,
              MASTER_RANK );

    // Generate filename.
    std::stringstream name;
    name << "solution" << refinement_level_;

    std::vector<double> remote_index
            ( mesh_ptr->num_entities ( mesh_ptr->tdim ( ) ),
              0. );
    std::vector<double> sub_domain
            ( mesh_ptr->num_entities ( mesh_ptr->tdim ( ) ),
              0. );
    std::vector<double> material_number
            ( mesh_ptr->num_entities ( mesh_ptr->tdim ( ) ),
              0. );

    std::vector<double> local_index
            ( mesh_ptr->num_entities ( mesh_ptr->tdim ( ) ),
              -1. );

    for ( mesh::EntityIterator it = mesh_ptr->begin
          ( mesh_ptr->tdim ( ) );
          it != mesh_ptr->end ( mesh_ptr->tdim ( ) );
          ++it )
    {
        int temp1, temp2, temp3;
        mesh_ptr->get_attribute_value ( "_remote_index_",
                                        mesh_ptr->tdim ( ),
                                        it->index ( ),
                                        &temp1 );
        mesh_ptr->get_attribute_value ( "_sub_domain_",
                                        mesh_ptr->tdim ( ),
                                        it->index ( ),
                                        &temp2 );

        mesh_ptr->get_attribute_value ( "_local_index_",
                                        mesh_ptr->tdim ( ),
                                        it->index ( ),
                                        &temp3 );

        remote_index.at ( it->index ( ) ) = temp1;
        sub_domain.at ( it->index ( ) ) = temp2;
        local_index.at ( it->index ( ) ) = temp3;
        material_number.at ( it->index ( ) ) = mesh_ptr->get_material_number
                ( mesh_ptr->tdim ( ),
                  it->index ( ) );
    }

    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ) ), "u" );

    // visualize error measures
    visu.visualize_cell_data ( rho_jump_, "rho_E" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.visualize_cell_data ( local_index, "_local_index_" );
    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.write ( name.str ( ) );
}

void JumpTest::adapt ( )
{

    std::vector<int> refinemnts ( mesh_->num_entities ( DIMENSION ), 0 );

    if ( rank_ == 0 )
    {
        refinemnts[0] = 1;
        //		refinemnts[1] = 1;
    }
    if ( rank_ == 1 )
    {
        refinemnts[0] = 1;
    }

    SharedVertexTable shared_verts;
    mesh_ = mesh_->refine ( refinemnts );
#ifdef USE_MESH_P4EST    
    mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_P4EST );
#else
    mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_DBVIEW );
#endif
    refinement_level_++;
    std::vector<int> local_index;
    for ( int l = 0; l < mesh_->num_entities ( DIMENSION ); ++l )
    {
        local_index.push_back ( l );
    }
    AttributePtr local_index_attr ( new IntAttribute ( local_index ) );
    mesh_->add_attribute ( "_local_index_", DIMENSION, local_index_attr );

    PVtkWriter writer ( comm_ );
    std::ostringstream name;
    name << "jump_term_mesh_" << refinement_level_ << ".pvtu";
    std::string output_file = name.str ( );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

