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

/// \author Staffan Ronnas<br>Julian Kraemer

#include "poisson_tutorial.h"

static const char* PARAM_FILENAME = "poisson_tutorial.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class PoissonTutorial
{
  public:

    PoissonTutorial ( const std::string& param_filename, const std::string& path_mesh )
    : path_mesh ( path_mesh ),
    comm_ ( MPI_COMM_WORLD ),
    rank_ ( -1 ),
    num_partitions_ ( -1 ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    rhs_ ( 0 ), sol_ ( 0 ), matrix_ ( 0 ),
    is_done_ ( false ),
    refinement_level_ ( 0 )
    {
        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );
    }

    // Main algorithm

    void run ( )
    {
        // Construct / read in the initial mesh.
        build_initial_mesh ( );
        // Main adaptation loop.
        while ( !is_done_ )
        {
            // Initialize space and linear algebra.
            prepare_system ( );
            // Compute the stiffness matrix and right-hand side.
            assemble_system ( );
            // Solve the linear system.
            solve_system ( );
            // Compute the error to the exact solution.
            compute_error ( );
            // Visualize the solution and the errors.
            visualize ( );
            adapt ( );
        }
    }

    ~PoissonTutorial ( )
    {
        delete matrix_;
        delete sol_;
        delete rhs_;
    }

  private:
    // Member functions

    // Read and distribute mesh.
    std::string path_mesh;
    void build_initial_mesh ( );
    // Setup space, linear algebra, and compute Dirichlet values.
    void prepare_system ( );

    // Compute the matrix and rhs.
    void assemble_system ( );
    // Compute solution x.
    void solve_system ( );

    // Compute errors compared to exact solution.
    void compute_error ( );
    // Visualize the results.
    void visualize ( );
    // Adapt the space (mesh and/or degree).
    void adapt ( );

    // member variables

    // MPI communicator.
    MPI_Comm comm_;
    // Local process rank and number of processes.
    int rank_, num_partitions_;

    // Parameter data read in from file.
    PropertyTree params_;

    // Local mesh and mesh on master process.
    MeshPtr mesh_, mesh_without_ghost_, master_mesh_;
    // Solution space.
    VectorSpace<double> space_;

    // Linear algebra couplings helper object.
    Couplings<double> couplings_;
    // Vectors for solution and load vector.
    CoupledVector<Scalar>* rhs_, *sol_;
    // System matrix.
    CoupledMatrix<Scalar>* matrix_;

    // Global assembler.
    StandardGlobalAssembler<double> global_asm_;

    // Flag for stopping adaptive loop.
    bool is_done_;
    // Current refinement level.
    int refinement_level_;

    // Vectors for error norms
    std::vector<Scalar> L2_err_, H1_err_;

    // Dof id:s for Dirichlet boundary conditions.
    std::vector<int> dirichlet_dofs_;
    // Dof values for Dirichlet boundary conditions.
    std::vector<Scalar> dirichlet_values_;
}; // end class PoissonTutorial

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    // set default parameter file
    std::string param_filename ( PARAM_FILENAME );
    std::string path_mesh;
    // if set take parameter file specified on console
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }
    // if set take mesh following path specified on console
    if ( argc > 2 )
    {
        path_mesh = std::string ( argv[2] );
    }
    try
    {
        // Create log files for INFO and DEBUG output
        //std::ofstream info_log("poisson_tutorial_info_log");
        LogKeeper::get_log ( "info" ).set_target ( &( std::cout ) );
        //std::ofstream debug_log("poisson_tutorial_debug_log");
        LogKeeper::get_log ( "debug" ).set_target ( &( std::cout ) );
        std::ofstream error_log ( "poisson_tutorial_error_log" );
        LogKeeper::get_log ( "error" ).set_target ( &( std::cout ) );

        // Create application object and run it
        PoissonTutorial app ( param_filename, path_mesh );
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

//////////////// PoissonTutorial implementation //////////////

void PoissonTutorial::build_initial_mesh ( )
{
#ifdef USE_MESH_P4EST
    mesh::IMPL impl = mesh::IMPL_P4EST;
#else
    mesh::IMPL impl = mesh::IMPL_DBVIEW;
#endif

    // Read in the mesh on the master process. The mesh is chosen according to the dimension of the problem.
    if ( rank_ == MASTER_RANK )
    {
        std::string mesh_name;

        switch ( DIMENSION )
        {
            case 1:
            {
                mesh_name =
                        params_["Mesh"]["Filename1"].get<std::string>( "unit_line.inp" );
                break;
            }
            case 2:
            {
                mesh_name =
                        params_["Mesh"]["Filename2"].get<std::string>( "unit_square.inp" );
                break;
            }
            case 3:
            {
                mesh_name =
                        params_["Mesh"]["Filename3"].get<std::string>( "unit_cube.inp" );
                break;
            }

            default: assert ( 0 );
        }
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

        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period, impl );

        // Refine the mesh until the initial refinement level is reached.
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( 3 );

        if ( initial_ref_lvl > 0 )
        {
            master_mesh_ = master_mesh_->refine_uniform_seq ( initial_ref_lvl );
        }
        refinement_level_ = initial_ref_lvl;
    }

    // 1D parallel execution is not yet implemented.
    if ( DIMENSION == 1 )
    {
        mesh_ = master_mesh_;
    }
    else
    {
        MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
        SharedVertexTable shared_verts;

        int uniform_ref_steps;
        mesh_without_ghost_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &uniform_ref_steps, impl );

        assert ( mesh_without_ghost_ != 0 );
        mesh_ = compute_ghost_cells ( *mesh_without_ghost_, comm_, shared_verts, impl );

        // Write out mesh of initial refinement level      
        PVtkWriter writer ( comm_ );
        std::ostringstream name;
        name << "poisson_tutorial_mesh_" << refinement_level_ << ".pvtu";
        std::string output_file = name.str ( );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );
    }
}

void PoissonTutorial::prepare_system ( )
{
    // Assign degrees to each element.

    const int fe_degree = params_["Mesh"]["FeDegree"].get<int>( 1 );
    std::vector< int > degrees ( 1, fe_degree );

    // Initialize the VectorSpace object.
    space_.Init ( degrees, *mesh_ );

    // Setup couplings object.
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    // Setup linear algebra objects.
    delete matrix_;
    delete rhs_;
    delete sol_;

    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( "CoupledMatrix" ) )->
            params ( params_["LinearAlgebra"] );
    matrix_->Init ( comm_, couplings_ );
    CoupledVectorFactory<Scalar> CoupVecFact;
    rhs_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( "CoupledVector" ) )->
            params ( params_["LinearAlgebra"] );
    sol_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( "CoupledVector" ) )->
            params ( params_["LinearAlgebra"] );
    rhs_->Init ( comm_, couplings_ );
    sol_->Init ( comm_, couplings_ );

    // Initialize structure of LA objects.
    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );

    rhs_->InitStructure ( );
    sol_->InitStructure ( );

    // Zero all linear algebra objects.
    matrix_->Zeros ( );
    rhs_->Zeros ( );
    sol_->Zeros ( );

    // Compute Dirichlet BC dofs and values using known exact solution.
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    DirichletZero zero;

    // The function compute_dirichlet_dofs_and_values des not yet work for 1D.
    if ( DIMENSION == 1 )
    {
        dirichlet_values_.resize ( 2, 0.0 );

        // Loop over all cells.
        for ( EntityIterator facet_it = mesh_->begin ( DIMENSION - 1 ), facet_end = mesh_->end ( DIMENSION - 1 );
              facet_it != facet_end; ++facet_it )
        {

            // Returns the number of neighbors for each cell, to check if it is on the facet.
            const EntityCount num_cell_neighbors = facet_it->num_incident_entities ( DIMENSION );

            // If it lies on the facet, the corresponding DOF is a Dirichlet DOF and is added to dirichlet_dofs_.
            if ( num_cell_neighbors == 1 )
            {
                std::vector<int> dof_number_;
                space_.dof ( ).get_dofs_on_subentity ( 0, facet_it->begin_incident ( DIMENSION )->index ( ), 0, facet_it->index ( ), dof_number_ );
                dirichlet_dofs_.push_back ( dof_number_[0] );
            }
        }
    }
    else
    {
        compute_dirichlet_dofs_and_values ( zero, space_, 0, dirichlet_dofs_,
                                            dirichlet_values_ );
    }
}

void PoissonTutorial::assemble_system ( )
{
    // Assemble matrix and right-hand-side vector.
    LocalPoissonAssembler<ExactSol> local_asm;
    global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, *rhs_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
}

void PoissonTutorial::solve_system ( )
{
    LinearSolver<LAD>* solver_;
    LinearSolverFactory<LAD> SolFact;
    solver_ = SolFact.Get (
                            params_["LinearSolver"]["Name"].get<std::string>( "CG" ) )->
            params ( params_["LinearSolver"] );
    solver_->SetupOperator ( *matrix_ );
    solver_->Solve ( *rhs_, sol_ );
    delete solver_;
}

void PoissonTutorial::compute_error ( )
{
    // prepare sol_ for post processing
    sol_->UpdateCouplings ( );

    L2_err_.clear ( );
    H1_err_.clear ( );

    // Compute square of the L2 error on each element, putting the
    // values into L2_err_.
    L2ErrorIntegrator<ExactSol> L2_int ( *( sol_ ) );
    global_asm_.assemble_scalar ( space_, L2_int, L2_err_ );

    // Create attribute with L2 error for output.
    AttributePtr L2_err_attr ( new DoubleAttribute ( L2_err_ ) );
    mesh_->add_attribute ( "L2 error", DIMENSION, L2_err_attr );
    double total_L2_err = std::accumulate ( L2_err_.begin ( ), L2_err_.end ( ), 0. );
    double global_L2_err = 0.;
    MPI_Reduce ( &total_L2_err, &global_L2_err, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK,
                 comm_ );
    LOG_INFO ( "error", "Local L2 error on partition " << rank_ << " = "
               << std::sqrt ( total_L2_err ) );

    // Compute square of the H1 error on each element, putting the
    // values into H1_err_.

    H1ErrorIntegrator<ExactSol> H1_int ( *( sol_ ) );
    global_asm_.assemble_scalar ( space_, H1_int, H1_err_ );

    // Create attribute with H1 error for output.
    AttributePtr H1_err_attr ( new DoubleAttribute ( H1_err_ ) );
    mesh_->add_attribute ( "H1 error", DIMENSION, H1_err_attr );
    double total_H1_err = std::accumulate ( H1_err_.begin ( ), H1_err_.end ( ), 0. );
    double global_H1_err = 0.;
    MPI_Reduce ( &total_H1_err, &global_H1_err, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK,
                 comm_ );
    LOG_INFO ( "error", "Local H1 error on partition " << rank_ << " = "
               << std::sqrt ( total_H1_err ) );

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "Global L2 error = " << std::sqrt ( global_L2_err ) << "\n"
                << "Global H1 error = " << std::sqrt ( global_H1_err ) << "\n";
    }
}

void PoissonTutorial::visualize ( )
{
    // Setup visualization object.
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    // Generate filename.
    std::stringstream name;
    name << "solution" << refinement_level_;

    std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> material_number ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
          it != mesh_->end ( mesh_->tdim ( ) );
          ++it )
    {
        int temp1, temp2;
        if ( DIMENSION > 1 )
        {
            mesh_->get_attribute_value ( "_remote_index_", mesh_->tdim ( ),
                                         it->index ( ),
                                         &temp1 );
            mesh_->get_attribute_value ( "_sub_domain_", mesh_->tdim ( ),
                                         it->index ( ),
                                         &temp2 );
            remote_index.at ( it->index ( ) ) = temp1;
            sub_domain.at ( it->index ( ) ) = temp2;
        }
        material_number.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
    }

    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ) ), "u" );

    // visualize error measures
    visu.visualize_cell_data ( L2_err_, "L2" );
    visu.visualize_cell_data ( H1_err_, "H1" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.write ( name.str ( ) );
}

void PoissonTutorial::adapt ( )
{
    // Refine mesh on master process. 1D parallel execution is not yet implemented.
    if ( DIMENSION == 1 )
    {
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
            if ( refinement_level_ >= final_ref_level )
            {
                is_done_ = true;
            }
            else
            {
                mesh_ = mesh_->refine ( );
                ++refinement_level_;
            }
        }
    }
    else
    {
#ifdef USE_MESH_P4EST
        const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
        if ( refinement_level_ >= final_ref_level )
        {
            is_done_ = true;
        }
        else
        {
            mesh_without_ghost_ = mesh_without_ghost_->refine ( );
            ++refinement_level_;
        }
        if ( !is_done_ )
        {
            SharedVertexTable shared_verts;
            std::vector<MasterSlave> period ( 0, MasterSlave ( 0., 0., 0., 0 ) );
            mesh_ = compute_ghost_cells ( *mesh_without_ghost_, comm_, shared_verts, mesh::IMPL_P4EST );
        }
#else
#    ifdef WITH_PARMETIS
        const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
        if ( refinement_level_ >= final_ref_level )
        {
            is_done_ = true;
        }
        else
        {
            mesh_without_ghost_ = mesh_without_ghost_->refine ( );
            ++refinement_level_;
        }

        if ( !is_done_ )
        {
            // Repartition the new mesh.
            ParMetisGraphPartitioner parmetis_partitioner;
            MeshPtr local_mesh = repartition_mesh ( mesh_without_ghost_, comm_, &parmetis_partitioner );
            assert ( local_mesh != 0 );
            mesh_.reset ( );
            SharedVertexTable shared_verts;
            mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
            local_mesh.reset ( );
        }
#    else
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
            if ( refinement_level_ >= final_ref_level )
            {
                is_done_ = true;
            }
            else
            {
                master_mesh_ = master_mesh_->refine ( );
                ++refinement_level_;
            }
        }
        // Broadcast information from master to slaves.
        MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
        MPI_Bcast ( &is_done_, 1, MPI_CHAR, MASTER_RANK, comm_ );

        if ( !is_done_ )
        {
            // Distribute the new mesh.
            MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_ );
            assert ( local_mesh != 0 );
            SharedVertexTable shared_verts;
            mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
        }
#    endif
#endif
        PVtkWriter writer ( comm_ );
        std::ostringstream name;
        name << "poisson_tutorial_mesh_" << refinement_level_ << ".pvtu";
        std::string output_file = name.str ( );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );
    }
}
