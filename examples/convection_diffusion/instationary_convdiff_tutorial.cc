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

#include "instationary_convdiff_tutorial.h"

#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif

// if you work with extra makefile BASEDIR is defined
// static const char* DATADIR = BASEDIR "/meshes/";
static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;
static const char* PARAM_FILENAME = "instationary_convdiff_tutorial.xml";

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Main program.

int main ( int argc, char** argv )
{
    // Initialize MPI.
    MPI_Init ( &argc, &argv );

    // Set default parameter file
    std::string param_filename ( PARAM_FILENAME );
    std::string path_mesh;
    // Specify parameter file on console
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }
    // Specify path to geometry mesh on console
    if ( argc > 2 )
    {
        path_mesh = std::string ( argv[2] );
    }
    try
    {
        // Run ConvDiff application.
        int rank = -1;
        MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
        if ( rank == MASTER_RANK )
        {
            std::cout << "==============================================="
                    << "=================================\n"
                    << "Running Convection-Diffusion tutorial...\n";
        }
        ConvectionDiffusionApp app ( param_filename, path_mesh );
        app.run ( );
    }
    catch ( std::exception& e )
    {
        std::cerr << "\nProgram ended with uncaught exception.\n";
        std::cerr << e.what ( ) << "\n";
        return -1;
    }

    // Finalize MPI.
    MPI_Finalize ( );
    return 0;
}

// -------------------------------------------------------------------

// Conctructor of ConvectionDiffusionApplication

ConvectionDiffusionApp::ConvectionDiffusionApp ( const std::string&
                                                 param_filename, const std::string& path_mesh )
: params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
comm_ ( MPI_COMM_WORLD ),
rank_ ( -1 ),
num_partitions_ ( -1 ),
master_rank_ ( MASTER_RANK ),
path_mesh ( path_mesh )
{

    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

}

// -------------------------------------------------------------------

// Destructor of ConvectionDiffusionApplication

ConvectionDiffusionApp::~ConvectionDiffusionApp ( )
{
#ifndef RUNGE_KUTTA
    delete sol_;
    delete sol_prev_;
    delete rhs_;
    delete matrix_;
    delete solver_;
#endif
}

// -------------------------------------------------------------------

// Run the Convection-Diffusion Application

void ConvectionDiffusionApp::run ( )
{
#ifndef RUNGE_KUTTA
    // Create log files for INFO and DEBUG output
    std::ofstream info_log ( params_["Output"]["LogFilename"].
                             get<std::string>( ).c_str ( ) );
    LogKeeper::get_log ( "info" ).set_target ( &info_log );
    std::ofstream debug_log ( params_["Output"]["DebugFilename"].
                              get<std::string>( ).c_str ( ) );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

    LOG_INFO ( "MPI Processes", num_partitions_ );

    // Setup mesh and distribute it if run in parallel
    read_and_distribute_mesh ( );

    // Initialize space and linear algebra
    prepare_system ( );

    // Time discretization
    int i = 0;
    time_ = 0.0;
    method = params_["Instationary"]["Method"].get<std::string>( );
    if ( method == "CrankNicolson" )
    {
        theta_ = 0.5;
    }
    else if ( method == "ImplicitEuler" )
    {
        theta_ = 1.;
    }
    else if ( method == "ExplicitEuler" )
    {
        theta_ = 0.;
    }
    else
    {
        theta_ = 0.5;
        std::cout << "=================================\n"
                << "Unknown method for solving instationary problem\n"
                << "=================================\n"
                << "Default method: Crank Nicolson method\n";
    }

    // Number of time steps to be calculated
    Tmax_ = params_["Instationary"]["MaxIts"].get<int>( );
    // Time stepp size
    delta_t_ = params_["Instationary"]["DeltaT"].get<double>( );
    // Coefficient of diffusion
    nu_ = params_["Application"]["Diffusion"].get<double>( );
    for ( i = 0; i < Tmax_; ++i )
    {
        if ( i > 0 )
        {
            // Compute the stiffness matrix and right-hand side
            assemble_system ( );

            // Solve the linear system
            solve_system ( );
        }
        // Visualize the solution and the errors.
        visualize ( i );
        time_ += delta_t_;
    }
#else
    // Read parameters.
    read_parameters ( );

    // Setup mesh and distribute it if run in parallel.
    read_and_distribute_mesh ( );

    // Initialize space and linear algebra.
    prepare_system ( );

    // Compute the stiffness matrix and right-hand side.
    assemble_system ( );

    // Implement time discretization.
    int i = 0;
    time_ = 0.0;

    for ( i = 0; i < Tmax_; ++i )
    {
        if ( i > 0 )
        {
            // Solve the linear system.
            solve_system ( );
        }

        // Visualize the solution and the errors.
        visualize ( i );

        time_ += delta_t_;
    }
#endif
}

// -------------------------------------------------------------------

void ConvectionDiffusionApp::read_parameters ( )
{
    // Parameters for application.
    refinement_level_ = params_["Mesh"]["RefinementLevel"].get<int>( );
    fe_degree_ = params_["FiniteElements"]["Degree"].get<int>( );
    nu_ = params_["Application"]["Diffusion"].get<double>( );

    // Parameters for linear solver.
    basis_size_ = params_["LinearSolver"]["SizeBasis"].get<int>( );
    maxits_ = params_["LinearSolver"]["MaxIterations"].get<int>( );
    abstol_ = params_["LinearSolver"]["AbsTolerance"].get<double>( );
    reltol_ = params_["LinearSolver"]["RelTolerance"].get<double>( );
    divtol_ = params_["LinearSolver"]["DivTolerance"].get<double>( );

    // Parameters for time-stepping
    Tmax_ = params_["Instationary"]["MaxIts"].get<int>( ); // Number of time steps to be calculated
    delta_t_ = params_["Instationary"]["DeltaT"].get<double>( ); // Time stepping size
}
// -------------------------------------------------------------------

// Homogenous Dirichlet boundary conditions

struct DirichletBC
{
    // Callback function to evaluate boundary values

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        return std::vector<double>( coords_on_face.size ( ), 0.0 );
    }
};

// -------------------------------------------------------------------

// Read in mesh and distribute it, if run in parallel

void ConvectionDiffusionApp::read_and_distribute_mesh ( )
{
    MeshPtr master_mesh ( 0 );

    refinement_level_ = params_["Mesh"]["RefinementLevel"].get<int>( );
    // Only master rank reads in mesh
    if ( rank_ == master_rank_ )
    {
        const std::string mesh_name =
                params_["Mesh"]["Filename"].get<std::string>( );
        std::string mesh_filename;
        if ( path_mesh.empty ( ) )
        {
            mesh_filename = std::string ( DATADIR ) + mesh_name;
        }
        else
        {
            mesh_filename = path_mesh + mesh_name;
        }
        master_mesh = read_mesh_from_file ( mesh_filename, DIM, DIM, 0 );
        assert ( master_mesh != 0 );

        // Global refinement of mesh
        for ( int r = 0; r < refinement_level_; ++r )
        {
            master_mesh = master_mesh->refine ( );
        }
        std::cout << "========================================\n"
                << "Mesh refinement level: "
                << refinement_level_ << "\n";
    }

    // Distribute mesh
    MeshPtr local_mesh = partition_and_distribute ( master_mesh, master_rank_, comm_ );
    assert ( local_mesh != 0 );

    // Compute ghost cells
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
}

// -------------------------------------------------------------------

void ConvectionDiffusionApp::prepare_system ( )
{
#ifndef RUNGE_KUTTA
    // Assign degrees to each element
    fe_degree_ = params_["FiniteElements"]["Degree"].get<int>( );
    std::vector<int> degrees ( 1, fe_degree_ );

    // Initialize the VectorSpace object
    space_.Init ( degrees, *mesh_ );

    // Setup couplings object
    couplings_.Clear ( );
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    // Initialize system matrix, solution vector and vector of right-hand side
    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    matrix_->Init ( comm_, couplings_ );
    CoupledVectorFactory<Scalar> CoupVecFact;
    sol_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    sol_->Init ( comm_, couplings_ );
    sol_prev_ = CoupVecFact.Get (
                                  params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    sol_prev_->Init ( comm_, couplings_ );

    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );
    matrix_->Zeros ( );

    sol_->InitStructure ( );
    sol_->Zeros ( );

    sol_prev_->InitStructure ( );
    sol_prev_->Zeros ( );

    rhs_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    rhs_->Init ( comm_, couplings_ );
    rhs_->InitStructure ( );
    rhs_->Zeros ( );

    // Prepare boundary conditions
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    DirichletBC bc;
    compute_dirichlet_dofs_and_values ( bc, space_, 0, dirichlet_dofs_,
                                        dirichlet_values_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct solution with dirichlet BC
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
    // Prepare linear system
    solver_ = SolFact.Get (
                            params_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( params_["LinearSolver"] );
#else
    // Assign degrees to each element.
    std::vector<int> degrees ( 1, fe_degree_ );

    // Initialize the VectorSpace object.
    space_.Init ( degrees, *mesh_ );

    // Setup couplings object.
    couplings_.Clear ( );
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    // Initialize matrices and vectors
    mass_.Init ( comm_, couplings_, la_platform, la_impl, la_matrix_format );
    mass_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ), vec2ptr ( sparsity.diagonal_cols ), sparsity.diagonal_rows.size ( ), vec2ptr ( sparsity.off_diagonal_rows ), vec2ptr ( sparsity.off_diagonal_cols ), sparsity.off_diagonal_rows.size ( ) );

    stiff_.Init ( comm_, couplings_, la_platform, la_impl, la_matrix_format );
    stiff_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ), vec2ptr ( sparsity.diagonal_cols ), sparsity.diagonal_rows.size ( ), vec2ptr ( sparsity.off_diagonal_rows ), vec2ptr ( sparsity.off_diagonal_cols ), sparsity.off_diagonal_rows.size ( ) );

    conv_.Init ( comm_, couplings_, la_platform, la_impl, la_matrix_format );
    conv_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ), vec2ptr ( sparsity.diagonal_cols ), sparsity.diagonal_rows.size ( ), vec2ptr ( sparsity.off_diagonal_rows ), vec2ptr ( sparsity.off_diagonal_cols ), sparsity.off_diagonal_rows.size ( ) );

    source_.Init ( comm_, couplings_, la_platform, la_impl );
    source_.InitStructure ( );
    source_.Zeros ( );

    k1.Init ( comm_, couplings_, la_platform, la_impl );
    k1.InitStructure ( );
    k1.Zeros ( );

    k2.Init ( comm_, couplings_, la_platform, la_impl );
    k2.InitStructure ( );
    k2.Zeros ( );

    k3.Init ( comm_, couplings_, la_platform, la_impl );
    k3.InitStructure ( );
    k3.Zeros ( );

    k4.Init ( comm_, couplings_, la_platform, la_impl );
    k4.InitStructure ( );
    k4.Zeros ( );

    help1.Init ( comm_, couplings_, la_platform, la_impl );
    help1.InitStructure ( );
    help1.Zeros ( );

    help2.Init ( comm_, couplings_, la_platform, la_impl );
    help2.InitStructure ( );
    help2.Zeros ( );

    help3.Init ( comm_, couplings_, la_platform, la_impl );
    help3.InitStructure ( );
    help3.Zeros ( );

    help4.Init ( comm_, couplings_, la_platform, la_impl );
    help4.InitStructure ( );
    help4.Zeros ( );

    s.Init ( comm_, couplings_, la_platform, la_impl );
    s.InitStructure ( );
    s.Zeros ( );

    sol_.Init ( comm_, couplings_, la_platform, la_impl );
    sol_.InitStructure ( );
    sol_.Zeros ( );

    rhs_.Init ( comm_, couplings_, la_platform, la_impl );
    rhs_.InitStructure ( );
    rhs_.Zeros ( );

    // Prepare boundary conditions
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    DirichletBC bc;
    compute_dirichlet_dofs_and_values ( bc, space_, 0, dirichlet_dofs_, dirichlet_values_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }

    // Setup linear solver
    solver_.InitControl ( maxits_, abstol_, reltol_, divtol_ );
    solver_.InitParameter ( "NoPreconditioning" );
    // set the matrix to be used as the operator
    solver_.SetupOperator ( mass_ );
#endif
}

// -------------------------------------------------------------------

// Function to assemble system matrix and right-hand side

void ConvectionDiffusionApp::assemble_system ( )
{
#ifndef RUNGE_KUTTA
    // Assemble matrix and right-hand side vector
    ConvectionDiffusionAssembler local_asm ( nu_ );

    // Initialize timestepping
    local_asm.set_timestep_parameters ( theta_, delta_t_ );
    local_asm.set_time ( time_ );
    local_asm.set_prev_solution ( sol_prev_ );

    if ( time_ == delta_t_ )
    {
        global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );
    }
    global_asm_.assemble_vector ( space_, local_asm, *rhs_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ),
                                    dirichlet_dofs_.size ( ), 1.0 );
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
#else
    //assemble mass matrix
    MassAssembler mass_asm;
    global_asm_.assemble_matrix ( space_, mass_asm, mass_ );

    //stiffness mass matrix
    StiffnessAssembler stiff_asm;
    global_asm_.assemble_matrix ( space_, stiff_asm, stiff_ );

    //assemble convection matrix
    ConvectionAssembler conv_asm;
    global_asm_.assemble_matrix ( space_, conv_asm, conv_ );

    //assemble source vector
    SourceAssembler source_asm;
    global_asm_.assemble_vector ( space_, source_asm, source_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        mass_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
        stiff_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
        conv_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );

        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }
#endif
}

// -------------------------------------------------------------------

// Function to solve system

void ConvectionDiffusionApp::solve_system ( )
{
#ifndef RUNGE_KUTTA
    // Solve linear system
    solver_->SetupOperator ( *matrix_ );
    solver_->Solve ( *rhs_, sol_ );

    // update the ghost DoFs, since in the next time step
    // they are needed in the assembly
    sol_->UpdateCouplings ( );
    // store old time solution (this also clones the ghost DoFs)
    sol_prev_->CloneFrom ( *sol_ );

#else

    // ------------------------------------------
    // setup rhs for k1
    u_.CloneFrom ( sol_ );
    double t = time_ - delta_t_;
    help1.Zeros ( );
    stiff_.VectorMult ( u_, &help1 );

    help2.Zeros ( );
    conv_.VectorMult ( u_, &help2 );

    help3.Zeros ( );
    help3.CloneFrom ( help2 );
    help3.Axpy ( help1, nu_ );

    s.Zeros ( );
    s.Axpy ( source_, std::exp ( -std::pow ( t, 10.0 ) ) );

    rhs_.Zeros ( );
    rhs_.CloneFrom ( s );
    rhs_.Axpy ( help3, -1.0 );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
        // correct solution with dirichlet BC
        k1.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }

    // Solve linear system for k1.
    solver_.Solve ( rhs_, &k1 );

    // ------------------------------------------
    // setup rhs for k2
    u_.CloneFrom ( sol_ );
    u_.Axpy ( k1, 0.5 * delta_t_ );
    t = time_ - 0.5 * delta_t_;

    help1.Zeros ( );
    stiff_.VectorMult ( u_, &help1 );

    help2.Zeros ( );
    conv_.VectorMult ( u_, &help2 );

    help3.Zeros ( );
    help3.CloneFrom ( help2 );
    help3.Axpy ( help1, nu_ );

    s.Zeros ( );
    s.Axpy ( source_, std::exp ( -std::pow ( t, 10.0 ) ) );

    rhs_.Zeros ( );
    rhs_.CloneFrom ( s );
    rhs_.Axpy ( help3, -1.0 );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
        // correct solution with dirichlet BC
        k2.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }

    // Solve linear system for k2.
    solver_.Solve ( rhs_, &k2 );

    // ------------------------------------------
    // setup rhs for k3
    u_.CloneFrom ( sol_ );
    u_.Axpy ( k2, 0.5 * delta_t_ );
    t = time_ - 0.5 * delta_t_;

    help1.Zeros ( );
    stiff_.VectorMult ( u_, &help1 );

    help2.Zeros ( );
    conv_.VectorMult ( u_, &help2 );

    help3.Zeros ( );
    help3.CloneFrom ( help2 );
    help3.Axpy ( help1, nu_ );

    s.Zeros ( );
    s.Axpy ( source_, std::exp ( -std::pow ( t, 10.0 ) ) );

    rhs_.Zeros ( );
    rhs_.CloneFrom ( s );
    rhs_.Axpy ( help3, -1.0 );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
        // correct solution with dirichlet BC
        k3.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }

    // Solve linear system for k3.
    solver_.Solve ( rhs_, &k3 );

    // ------------------------------------------
    // setup rhs for k4
    u_.CloneFrom ( sol_ );
    u_.Axpy ( k3, delta_t_ );
    t = time_;

    help1.Zeros ( );
    stiff_.VectorMult ( u_, &help1 );

    help2.Zeros ( );
    conv_.VectorMult ( u_, &help2 );

    help3.Zeros ( );
    help3.CloneFrom ( help2 );
    help3.Axpy ( help1, nu_ );

    s.Zeros ( );
    s.Axpy ( source_, std::exp ( -std::pow ( t, 10.0 ) ) );

    rhs_.Zeros ( );
    rhs_.CloneFrom ( s );
    rhs_.Axpy ( help3, -1.0 );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
        // correct solution with dirichlet BC
        k4.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }

    // Solve linear system for k4.
    solver_.Solve ( rhs_, &k4 );

    // ------------------------------------------
    // setup rhs for sol
    help4.Zeros ( );
    help4.Axpy ( k1, 1.0 );
    help4.Axpy ( k2, 2.0 );
    help4.Axpy ( k3, 2.0 );
    help4.Axpy ( k4, 1.0 );

    rhs_.Zeros ( );
    rhs_.CloneFrom ( sol_ );

    rhs_.Axpy ( help4, delta_t_ / 6.0 );

    sol_.CloneFrom ( rhs_ );
#endif
}

// -------------------------------------------------------------------

void ConvectionDiffusionApp::visualize ( int time )
{

    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    // Generate filename
    std::stringstream input;
    input << num_partitions_ << "_solution_" << refinement_level_ << "_" << time;

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

#ifndef RUNGE_KUTTA
    sol_->UpdateCouplings ( );
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ) ), "u" );
#else
    sol_.UpdateCouplings ( );
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_ ), "u" );
#endif
    // Visualize the source and convection terms and all attributes
    visu.visualize ( EvalSourceFunction ( space_, time_ ), "s" );
    visu.visualize ( EvalConvectionFunction ( space_, 0 ), "a_x" );
    visu.visualize ( EvalConvectionFunction ( space_, 1 ), "a_y" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );

#ifdef WITH_HDF5
    //XDMF Visualization
    static bool first_time = true;
    XdmfVisualization<double> xdmf_visu ( comm_, "xdmf_visualization", !first_time );
    //write grid
    if ( first_time )
    {
        int xdmf_num_intervals = 1;
        std::string grid_name = "mesh";
        xdmf_visu.add_to_view ( &space_, 0, "_sol_", "u" );
        xdmf_visu.write_view ( xdmf_num_intervals, grid_name );
        first_time = false;
    }

    xdmf_visu.add_timestep ( time_ );
#    ifndef RUNGE_KUTTA
    xdmf_visu.write_solution_vector ( sol_, "_sol_" );
#    else
    xdmf_visu.write_solution_vector ( &sol_, "_sol_" );
#    endif

    xdmf_visu.print_xdmf ( );
#endif
}

// -------------------------------------------------------------------
