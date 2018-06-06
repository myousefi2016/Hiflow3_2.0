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

/// \author Thomas Gengenbach

#include "convdiff_dg.h"

#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "hiflow.h"

static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;
static const int DEBUG_LEVEL = 0;
const TDim tdim = DIM;
const GDim gdim = DIM;

/// Program main loop: setup MPI, read parameters and start the application

int main ( int argc, char** argv )
{

    MPI_Init ( &argc, &argv );

    Timer main_timer;

    // Read parameters
    std::string fName;
    fName = std::string ( DATADIR ) + "../../examples/convection_diffusion/convdiff_dg.xml";
    PropertyTree param ( fName, MASTER_RANK, MPI_COMM_WORLD );

    // Start Logging
    std::ofstream info_file ( param["Output"]["LogFilename"].get<std::string>( ).c_str ( ) );
    LogKeeper::get_log ( "info" ).set_target ( &info_file );

    std::ofstream debug_file ( param["Output"]["DebugFilename"].get<std::string>( ).c_str ( ) );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_file );

    //    LOG_INFO("Running on System", std::string(popen("uname -a", "r")));
    int num_procs = -1;
    MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );
    LOG_INFO ( "MPI Processes", num_procs );

    // Run ConvDiff apllication
    std::vector<double> beta ( tdim, 1.0 );
    double nu = 1. / 121.;
    ConvectionDiffusionApp app ( param, beta, nu );
    app.run ( );

    // flush log here to avoid problems
    LogKeeper::get_log ( "info" ).flush ( );

    main_timer.stop ( );
    LOG_INFO ( "Total program run time [without MPI Init and Finalize]", main_timer );

    MPI_Finalize ( );
    return 0;
}

///////////// ConvectionDiffusionAssembler /////////////

void ConvectionDiffusionAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
        LocalMatrix& lm )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );

    Vec<DIM, double> beta;
    for ( int var = 0; var < DIM; ++var )
    {
        beta[var] = beta_[var];
    }

    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    // loop q
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );

        // assemble a1(u,v) = \int \nu * {\grad(u) : \grad(v)}
        for ( int i = 0; i < num_dofs ( 0 ); ++i )
        {
            for ( int j = 0; j < num_dofs ( 0 ); ++j )
            {
                lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                        wq * ( nu_ * dot ( grad_phi ( j, q, 0 ), grad_phi ( i, q, 0 ) ) ) * dJ;
            }
        }

        // assemble a2(u,v) = \int { (beta*\grad{u})*v }
        for ( int i = 0; i < num_dofs ( 0 ); ++i )
        {
            for ( int j = 0; j < num_dofs ( 0 ); ++j )
            {
                lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                        wq * ( dot ( beta, grad_phi ( j, q, 0 ) ) * phi ( i, q, 0 ) ) * dJ;
            }
        }
    }
}

void ConvectionDiffusionAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );
        const double gamma = 1.0 / nu_;

        const double x1 = x ( q )[0];
        const double x2 = x ( q )[1];
        // l0(v) = \int( dot(u_n - u_k, v))
        for ( int i = 0; i < num_dofs ( 0 ); ++i )
        {
            lv[dof_index ( i, 0 )] +=
                    wq * ( gamma * ( -nu_ * ( gamma * exp ( ( x1 - 1 ) * gamma )
                    * ( x1 * gamma + 2 )
                    * sin ( M_PI * x2 ) + M_PI * M_PI * x1
                    * ( -1.0 * exp ( -1.0 * gamma ) )
                    * ( exp ( x1 * gamma ) - exp ( gamma ) )
                    * sin ( M_PI * x2 ) )
                    + ( beta_[0] ) * ( M_PI * x1 * exp ( -gamma )
                    * ( exp ( x1 * gamma ) - exp ( gamma ) )
                    * cos ( M_PI * x2 ) )
                    + ( beta_[1] ) * ( sin ( M_PI * x2 )
                    * ( exp ( ( x1 - 1 ) * gamma )
                    * ( x1 * gamma + 1 ) - 1 ) ) )
                    * phi ( i, q, 0 ) ) * dJ;
        }
    }
}

///////////// ConvectionDiffusionApplication /////////////

ConvectionDiffusionApp::ConvectionDiffusionApp ( PropertyTree &param, const std::vector<double>& beta, double nu )
: param_ ( &param ),
comm_ ( MPI_COMM_WORLD ),
rank_ ( -1 ),
num_partitions_ ( -1 ),
master_rank_ ( 0 ),
nu_ ( nu ),
beta_ ( beta )

{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );
    refinement_level_ = ( *param_ )["Mesh"]["RefinementLevel"].get<int>( );
    mesh_filename_ = std::string ( DATADIR ) + ( *param_ )["Mesh"]["Filename"].get<std::string>( );
    std::cout << "Filename: " << mesh_filename_ << "\n";
    read_and_distribute_mesh ( mesh_filename_ );

    // setup linear algebra platform
    la_sys_.Platform = APP_PLATFORM;
    init_platform ( la_sys_ );

    local_asm_ = new ConvectionDiffusionAssembler ( beta_, nu_ );
}

ConvectionDiffusionApp::~ConvectionDiffusionApp ( )
{
    delete local_asm_;
}

void ConvectionDiffusionApp::run ( )
{
    prepare_space ( );
    prepare_bc ( );
    prepare_lin_alg_structures ( );
    prepare_linear_solver ( );
    assemble_system ( );
    solve_system ( );
    visualize_solution ( );
}

void ConvectionDiffusionApp::prepare_space ( )
{
    // setup space Q2 elements
    std::vector<int> degrees ( 1, 2 );
    std::vector<bool> is_cg ( 1, false );
    space_.Init ( degrees, *mesh_, is_cg );
}

void ConvectionDiffusionApp::prepare_lin_alg_structures ( )
{
    // Initialize linear algebra structures
    couplings_.Init ( comm_, space_.dof ( ) );

    // compute matrix graph to build matrix structure

    std::vector < std::vector<bool> > coupling_vars ( 1, std::vector<bool>( 1, true ) );

    SparsityStructure sparsity;

    DGGlobalAssembler<double> global_asm;
    global_asm.compute_sparsity_structure ( space_, sparsity, &coupling_vars );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows, sparsity.off_diagonal_cols );

    // Initialize matrices and vectors
    matrix_.Init ( comm_, couplings_, la_sys_.Platform,
                   ConvectionDiffusionApp::APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    matrix_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                            vec2ptr ( sparsity.diagonal_cols ),
                            sparsity.diagonal_rows.size ( ),
                            vec2ptr ( sparsity.off_diagonal_rows ),
                            vec2ptr ( sparsity.off_diagonal_cols ),
                            sparsity.off_diagonal_rows.size ( ) );

    sol_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    sol_.InitStructure ( );
    sol_.Zeros ( );

    rhs_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    rhs_.InitStructure ( );
    rhs_.Zeros ( );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }
}

void ConvectionDiffusionApp::prepare_linear_solver ( )
{
    int maxits = ( *param_ )["LinearSolver"]["GMRES"]["MaximumIterations"].get<int>( );
    double abstol = ( *param_ )["LinearSolver"]["GMRES"]["AbsoluteTolerance"].get<double>( );
    double reltol = ( *param_ )["LinearSolver"]["GMRES"]["RelativeTolerance"].get<double>( );
    double divtol = ( *param_ )["LinearSolver"]["GMRES"]["DivergenceLimit"].get<double>( );
    int basis_size = ( *param_ )["LinearSolver"]["GMRES"]["BasisSize"].get<int>( );
    solver_.InitControl ( maxits, abstol, reltol, divtol );
    solver_.InitParameter ( basis_size, "NoPreconditioning" );
    // set the matrix to be used as the operator
    solver_.SetupOperator ( matrix_ );
}

void ConvectionDiffusionApp::assemble_system ( )
{

    // Get DG Parameters
    dg_theta_ = ( *param_ )["DGMethod"]["Theta"].get<double>( );
    dg_gamma_ = ( *param_ )["DGMethod"]["Gamma"].get<double>( );

    DGGlobalAssembler<double> global_asm_dg;
    global_asm_dg.assemble_matrix ( space_, *local_asm_, matrix_ );
    global_asm_dg.assemble_vector ( space_, *local_asm_, rhs_ );
    DGJumpMatrixAssembler jump_asm ( beta_, nu_, dg_theta_, dg_gamma_ );
    global_asm_dg.assemble_interface_matrix ( space_, jump_asm, matrix_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        matrix_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                                   static_cast < LAD::DataType > ( 1. ) );
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }

    rhs_.UpdateCouplings ( );
    sol_.UpdateCouplings ( );
}

void ConvectionDiffusionApp::solve_system ( )
{
    // solve linear system
    solver_.SetupOperator ( matrix_ );
    solver_.Solve ( rhs_, &sol_ );
}

// Homogenous Dirichlet boundary conditions

struct Dirichlet_0_BC
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        //homogenous Dirichlet boundary
        return std::vector<double>( coords_on_face.size ( ), 0. );
    }
};

void ConvectionDiffusionApp::prepare_bc ( )
{
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    Dirichlet_0_BC bc;
    compute_dirichlet_dofs_and_values ( bc, space_, 0,
                                        dirichlet_dofs_, dirichlet_values_ );
}

void ConvectionDiffusionApp::read_and_distribute_mesh ( const std::string& filename )
{
    MeshPtr master_mesh ( 0 );

    if ( rank_ == master_rank_ )
    {
        master_mesh = read_mesh_from_file ( filename, DIM, DIM, 0 );
        assert ( master_mesh != 0 );

        for ( int r = 0; r < refinement_level_; ++r )
        {
            master_mesh = master_mesh->refine ( );
        }
    }

    // Distribute mesh
    MeshPtr local_mesh = partition_and_distribute ( master_mesh, master_rank_, MPI_COMM_WORLD );
    assert ( local_mesh != 0 );

    // Compute ghost cells
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, MPI_COMM_WORLD, shared_verts );
}

void ConvectionDiffusionApp::visualize_solution ( )
{
    sol_.UpdateCouplings ( );
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

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

    visu.visualize ( EvalFeFunction<LAD>( space_, sol_ ), "u" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );
}
