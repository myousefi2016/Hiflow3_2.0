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

#include "stabilized_convdiff_tutorial.h"

#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif

// if you work with extra makefile BASEDIR is defined
// static const char* DATADIR = BASEDIR "/meshes/";
static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;
static const char* PARAM_FILENAME = "stabilized_convdiff_tutorial.xml";

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
#ifdef SUPG
        std::cout << "...with SUPG.\n";
#endif
#ifdef GLS
        std::cout << "...with GLS.\n";
#endif
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
    delete sol_;
    delete rhs_;
    delete matrix_;
}

// -------------------------------------------------------------------

// Run the Convection-Diffusion Application

void ConvectionDiffusionApp::run ( )
{

    // Create log files for INFO and DEBUG output
    std::ofstream info_log ( params_["Output"]["LogFilename"].get<std::string>( ).c_str ( ) );
    LogKeeper::get_log ( "info" ).set_target ( &info_log );
    std::ofstream debug_log ( params_["Output"]["DebugFilename"].get<std::string>( ).c_str ( ) );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

    LOG_INFO ( "MPI Processes", num_partitions_ );

    // Setup mesh and distribute it if run in parallel
    read_and_distribute_mesh ( );

    // Initialize space and linear algebra
    prepare_system ( );

    // Compute the stiffness matrix and right-hand side
    assemble_system ( );

    // Solve the linear system
    solve_system ( );

    // Visualize the solution and the errors.
    visualize ( );
}

// -------------------------------------------------------------------

// Setup mesh and distribute it, if run in parallel

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

// Dirichlet boundary conditions

struct DirichletBC
{
    // Callback function to evaluate boundary values

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        // Array with Dirichlet values for dof:s on boundary face
        std::vector<Coordinate> coordinates;
        face.get_coordinates ( coordinates );

        // Dirichlet BC on (x,y=0)
        if ( ( coordinates[1] == 0 ) && ( coordinates[1 + DIM] == 0 ) )
        {
            return std::vector<double>( coords_on_face.size ( ), 0.0 );
        }

        // Dirichlet BC on (x = 0, y)
        if ( ( coordinates[0] == 0 ) && ( coordinates[DIM] == 0 ) )
        {
            std::vector<double> val ( coords_on_face.size ( ), 0.0 );

            // y >= 0.2
            if ( ( coordinates[1] >= 0.2 ) && ( coordinates[1 + DIM] >= 0.2 ) )
            {
                return std::vector<double>( coords_on_face.size ( ), 1.0 );
            }

            if ( coordinates[1] >= 0.2 )
            {
                val[0] = 1.0;
            }
            if ( coordinates[1 + DIM] >= 0.2 )
            {
                val[0] = 1.0;
            }
            return val;
        }
        return std::vector<double>( 0, 0.0 );
    }
};

// -------------------------------------------------------------------

void ConvectionDiffusionApp::prepare_system ( )
{
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

    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );
    matrix_->Zeros ( );

    sol_->InitStructure ( );
    sol_->Zeros ( );

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
}

// -------------------------------------------------------------------

// Function to assemble system matrix and right-hand side

void ConvectionDiffusionApp::assemble_system ( )
{
    // Assemble matrix and right-hand side vector

    // Coefficient of diffusion
    nu_ = params_["Application"]["Diffusion"].get<double>( );

    a_.resize ( DIM );
    a_[0] = cos ( M_PI / 6 );
    a_[1] = sin ( M_PI / 6 );

    double norm_a = 0.0;
    for ( int var = 0; var < DIM; ++var )
    {
        norm_a += a_[var] * a_[var];
    }
    norm_a = std::sqrt ( norm_a );

    // Compute Peclet-number
    Scalar pec = ( std::pow ( 0.5, refinement_level_ ) ) * norm_a / ( 2 * nu_ );
    std::cout << "       Peclet number     " << pec << ".\n";
    // Assemble matrix and rhs
    ConvectionDiffusionAssembler local_asm ( a_, nu_, pec );

    global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, *rhs_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ),
                                    dirichlet_dofs_.size ( ), 1.0 );
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
}

// -------------------------------------------------------------------

// Function to prepare linear solver and solve system

void ConvectionDiffusionApp::solve_system ( )
{
    // Prepare linear system
    solver_ = SolFact.Get (
                            params_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( params_["LinearSolver"] );
    // Solve linear system
    solver_->SetupOperator ( *matrix_ );
    solver_->Solve ( *rhs_, sol_ );

    delete solver_;
}

// -------------------------------------------------------------------

void ConvectionDiffusionApp::visualize ( )
{
    // Setup visualization object
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    // Generate filename
    std::stringstream input;
    input << num_partitions_ << "_solution_" << refinement_level_;

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

    sol_->UpdateCouplings ( );
    // Visualize the solution and all attributes
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ) ), "u" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );

    visu.write ( input.str ( ) );
}

// -------------------------------------------------------------------
