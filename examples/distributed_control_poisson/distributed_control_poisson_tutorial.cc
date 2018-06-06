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

/// \author Julian Kraemer

#include "distributed_control_poisson_tutorial.h"

static const char* PARAM_FILENAME = "distributed_control_poisson_tutorial.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class DistributedControlPoisson
{
  public:

    DistributedControlPoisson ( const std::string& param_filename,
                                const std::string& path_mesh )
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
            // Modify the space through refinement. Set is_done_ = true when finished.
            adapt ( );
        }
    }

    ~DistributedControlPoisson ( )
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
    const MPI_Comm comm_;
    // Local process rank and number of processes.
    int rank_, num_partitions_;

    // Parameter data read in from file.
    PropertyTree params_;

    // Local mesh and mesh on master process.
    MeshPtr mesh_, master_mesh_;
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

    // Vectors for error norm
    std::vector<Scalar> L2_err_y_, L2_err_p_, L2_err_u_;

    // Dof id:s for Dirichlet boundary conditions.
    std::vector<int> dirichlet_dofs_y_;
    std::vector<int> dirichlet_dofs_p_;

    // Dof values for Dirichlet boundary conditions.
    std::vector<Scalar> dirichlet_values_;
#ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> ilupp_;
    bool use_ilupp_;
#endif
}; // end class distributedControlPoisson

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
        std::ofstream info_log ( "distributed_control_poisson_tutorial_info_log" );
        LogKeeper::get_log ( "info" ).set_target ( &info_log );
        std::ofstream debug_log ( "distributed_control_poisson_tutorial_debug_log" );
        LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

        // Create application object and run it
        DistributedControlPoisson app ( param_filename, path_mesh );
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

//////////////// distributedControlPoisson implementation //////////////

void DistributedControlPoisson::build_initial_mesh ( )
{
    // Read in the mesh on the master process. The mesh is chosen according to
    // the dimension of the problem.
    if ( rank_ == MASTER_RANK )
    {
        std::string mesh_name;
        switch ( DIMENSION )
        {
            case 1:
            {
                mesh_name =
                        params_["Mesh"]["Filename1"].get<std::string>( );
                break;
            }
            case 2:
            {
                mesh_name =
                        params_["Mesh"]["Filename2"].get<std::string>( );
                break;
            }
            case 3:
            {
                mesh_name =
                        params_["Mesh"]["Filename3"].get<std::string>( );
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
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION,
                                             DIMENSION, 0 );

        // Refine the mesh until the initial refinement level is reached.
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );
        for ( int r = 0; r < initial_ref_lvl; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
            ++refinement_level_;
        }
    }

    // 1D parallel execution is not yet implemented.
    if ( DIMENSION == 1 )
    {
        mesh_ = master_mesh_;
    }
    else
    {
        MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );

        // Distribute mesh over all processes, and compute ghost cells
        MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK,
                                                        comm_ );
        assert ( local_mesh != 0 );
        SharedVertexTable shared_verts;
        mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );

        // Write out mesh of initial refinement level
        PVtkWriter writer ( comm_ );
        std::string output_file =
                std::string ( "distributed_control_poisson_tutorial_mesh.pvtu" );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );
    }
}

void DistributedControlPoisson::prepare_system ( )
{
    // Assign degrees to each element for each variable.
    const int fe_degree_p = params_["Mesh"]["FeDegreep"].get<int>( );
    const int fe_degree_y = params_["Mesh"]["FeDegreey"].get<int>( );
    const int fe_degree_u = params_["Mesh"]["FeDegreeu"].get<int>( );
    std::vector< int > degrees ( 3, 0 );
    degrees[0] = fe_degree_p;
    degrees[1] = fe_degree_y;
    degrees[2] = fe_degree_u;

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
    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    matrix_->Init ( comm_, couplings_ );
    CoupledVectorFactory<Scalar> CoupVecFact;
    rhs_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    sol_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
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
    dirichlet_dofs_y_.clear ( );
    dirichlet_dofs_p_.clear ( );
    dirichlet_values_.clear ( );

    DirichletZero zero;
    // The function compute_dirichlet_dofs_and_values des not yet work for 1D.
    if ( DIMENSION == 1 )
    {
        dirichlet_values_.resize ( 2, 0.0 );
        dirichlet_dofs_y_.resize ( 0 );
        dirichlet_dofs_p_.resize ( 0 );

        // Loop over all cells.
        for ( EntityIterator facet_it = mesh_->begin ( DIMENSION - 1 ), facet_end =
              mesh_->end ( DIMENSION - 1 ); facet_it != facet_end; ++facet_it )
        {

            // Returns the number of neighbors for each cell, to check if it is on
            // the facet.
            const EntityCount num_cell_neighbors =
                    facet_it->num_incident_entities ( DIMENSION );

            if ( num_cell_neighbors == 1 )
            {
                // If it lies on the facet, the corresponding DOF is a Dirichlet DOF
                // and is added to dirichlet_dofs_y_ or dirichlet_dofs_p, respectively.
                std::vector<int> dof_number_p;
                std::vector<int> dof_number_y;
                space_.dof ( ).get_dofs_on_subentity ( 0, facet_it->begin_incident
                                                       ( DIMENSION )->index ( ), 0, facet_it->index ( ), dof_number_p );
                space_.dof ( ).get_dofs_on_subentity ( 1, facet_it->begin_incident
                                                       ( DIMENSION )->index ( ), 0, facet_it->index ( ), dof_number_y );
                dirichlet_dofs_y_.push_back ( dof_number_p[0] );
                dirichlet_dofs_p_.push_back ( dof_number_y[0] );
            }
        }
    }
    else
    {
        compute_dirichlet_dofs_and_values ( zero, space_, 0, dirichlet_dofs_p_,
                                            dirichlet_values_ );
        compute_dirichlet_dofs_and_values ( zero, space_, 1, dirichlet_dofs_y_,
                                            dirichlet_values_ );
    }
}

void DistributedControlPoisson::assemble_system ( )
{
    // Assemble matrix and right-hand-side vector.
    const double lambda = params_["General"]["Lambda"].get<double>( );
    LocalDistributedControlPoissonAssembler local_asm ( lambda );
    global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, *rhs_ );

    if ( !dirichlet_dofs_p_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_p_ ), dirichlet_dofs_p_.size ( ), 1.0 );
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_p_ ), dirichlet_dofs_p_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_p_ ), dirichlet_dofs_p_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }

    if ( !dirichlet_dofs_y_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_y_ ), dirichlet_dofs_y_.size ( ), 1.0 );
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_y_ ), dirichlet_dofs_y_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_y_ ), dirichlet_dofs_y_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }

    // Communicate modified values outside if-block
    rhs_->UpdateCouplings ( );
    sol_->UpdateCouplings ( );
}

void DistributedControlPoisson::solve_system ( )
{
    LinearSolver<LAD>* solver_;
    LinearSolverFactory<LAD> SolFact;
    solver_ = SolFact.Get (
                            params_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( params_["LinearSolver"] );

#ifdef WITH_ILUPP
    // prepare preconditioner
    ilupp_.InitParameter ( params_["ILUPP"]["PreprocessingType"].get<int>( ),
                           params_["ILUPP"]["PreconditionerNumber"].get<int>( ),
                           params_["ILUPP"]["MaxMultilevels"].get<int>( ),
                           params_["ILUPP"]["MemFactor"].get<double>( ),
                           params_["ILUPP"]["PivotThreshold"].get<double>( ),
                           params_["ILUPP"]["MinPivot"].get<double>( ) );
    ilupp_.SetupOperator ( *matrix_ );
    solver_->SetupPreconditioner ( ilupp_ );

#endif

    solver_->SetupOperator ( *matrix_ );
    solver_->Solve ( *rhs_, sol_ );
    delete solver_;
}

void DistributedControlPoisson::compute_error ( )
{
    L2_err_p_.clear ( );
    L2_err_y_.clear ( );

    sol_->UpdateCouplings ( );
    // Compute square of the L2 error on each element, putting the
    // values into L2_err_.
    const double lambda = params_["General"]["Lambda"].get<double>( );

    ExactSol_p sol_p ( lambda );
    L2ErrorIntegrator_p<ExactSol_p> L2_int_p ( *( sol_ ), sol_p );
    global_asm_.assemble_scalar ( space_, L2_int_p, L2_err_p_ );

    ExactSol_y sol_y;
    L2ErrorIntegrator_y<ExactSol_y> L2_int_y ( *( sol_ ), sol_y );
    global_asm_.assemble_scalar ( space_, L2_int_y, L2_err_y_ );

    ExactSol_u sol_u;
    L2ErrorIntegrator_u<ExactSol_u> L2_int_u ( *( sol_ ), sol_u );
    global_asm_.assemble_scalar ( space_, L2_int_u, L2_err_u_ );

    // Create attribute with L2 error for output.
    AttributePtr L2_err_attr_p ( new DoubleAttribute ( L2_err_p_ ) );
    mesh_->add_attribute ( "L2 error p", DIMENSION, L2_err_attr_p );
    double total_L2_err_p = std::accumulate ( L2_err_p_.begin ( ), L2_err_p_.end ( ), 0. );
    double global_L2_err_p = 0.;
    MPI_Reduce ( &total_L2_err_p, &global_L2_err_p, 1, MPI_DOUBLE, MPI_SUM, 0,
                 comm_ );

    AttributePtr L2_err_attr_y ( new DoubleAttribute ( L2_err_y_ ) );
    mesh_->add_attribute ( "L2 error y", DIMENSION, L2_err_attr_y );
    double total_L2_err_y = std::accumulate ( L2_err_y_.begin ( ), L2_err_y_.end ( ), 0. );
    double global_L2_err_y = 0.;
    MPI_Reduce ( &total_L2_err_y, &global_L2_err_y, 1, MPI_DOUBLE, MPI_SUM, 0,
                 comm_ );

    AttributePtr L2_err_attr_u ( new DoubleAttribute ( L2_err_u_ ) );
    mesh_->add_attribute ( "L2 error u", DIMENSION, L2_err_attr_u );
    double total_L2_err_u = std::accumulate ( L2_err_u_.begin ( ), L2_err_u_.end ( ), 0. );
    double global_L2_err_u = 0.;
    MPI_Reduce ( &total_L2_err_u, &global_L2_err_u, 1, MPI_DOUBLE, MPI_SUM, 0,
                 comm_ );

    // Sum of the 3 Errors
    double error_sum = 0.;
    error_sum = global_L2_err_y + global_L2_err_p + global_L2_err_u;
    LOG_INFO ( "error", "Local L2 error on partition " << rank_ << " = "
               << std::sqrt ( error_sum ) );

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "Global L2 error = " << std::sqrt ( error_sum ) << "\n";
        std::cout << "u error = " << std::sqrt ( global_L2_err_u ) << "\n";
        std::cout << "p error = " << std::sqrt ( global_L2_err_p ) << "\n";
        std::cout << "y error = " << std::sqrt ( global_L2_err_y ) << "\n";
    }
}

void DistributedControlPoisson::visualize ( )
{
    // Setup visualization object.
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );
    std::vector<std::string> names ( 3 );
    names[0] = "p";
    names[1] = "y";
    names[2] = "u";

    // Generate filename.
    std::stringstream input;
    input << "solution" << refinement_level_;

    std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> material_number ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
          it != mesh_->end ( mesh_->tdim ( ) );
          ++it )
    {
        int temp1, temp2;
        mesh_->get_attribute_value ( "_remote_index_", mesh_->tdim ( ),
                                     it->index ( ),
                                     &temp1 );
        mesh_->get_attribute_value ( "_sub_domain_", mesh_->tdim ( ),
                                     it->index ( ),
                                     &temp2 );
        material_number.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
        remote_index.at ( it->index ( ) ) = temp1;
        sub_domain.at ( it->index ( ) ) = temp2;
    }

    // Visualize all variables.

    //UpdateCouplings not nescessary here because of update in compute error
    //sol_->UpdateCouplings();
    for ( int i = 0; i < 3; i++ )
    {
        visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), i ), names[i] );
    }

    //visualize attributes and error data
    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );

    visu.visualize_cell_data ( L2_err_u_, "L2 error u" );
    visu.visualize_cell_data ( L2_err_p_, "L2 error p" );
    visu.visualize_cell_data ( L2_err_y_, "L2 error y" );
    visu.write ( input.str ( ) );
}

void DistributedControlPoisson::adapt ( )
{
    // Refine mesh on master process. 1D parallel execution is not yet implemented.
    if ( DIMENSION == 1 )
    {
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( );
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
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( );
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
            MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK,
                                                            comm_ );
            assert ( local_mesh != 0 );
            SharedVertexTable shared_verts;
            mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
        }
    }
}
