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

#include "chorin.h"

static const char* DATADIR = MESHES_DATADIR;

/// Program main loop: setup MPI, read parameters and start the application

int main ( int argc, char** argv )
{

    MPI_Init ( &argc, &argv );

    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    // Set info log file
    std::ofstream info_log ( "chorin_instationary_navier_stokes_info_log" );
    LogKeeper::get_log ( "info" ).set_target ( &info_log );

    // Run application
    if ( rank == 0 )
        std::cout << "================================================================================\n"
            << "Solve the instationary incompressible Navier Stokes equations\n"
            << "with Chorin:s method.\n"
            << "================================================================================\n";

    InstationaryFlowApp application;
    application.run ( );

    MPI_Finalize ( );
    return 0;
}

void InstationaryFlowApp::EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F )
{
    // compute the residual vector
    compute_residual ( u, F );
    F->Scale ( -1. );

    // correct BC -- set Dirichlet dofs to desired value
    if ( !intermediate_dirichlet_dofs_.empty ( ) )
    {
        // no slip at upper and lower boundary
        std::vector<LAD::DataType> zeros ( intermediate_dirichlet_dofs_.size ( ), 0. );
        F->SetValues ( vec2ptr ( intermediate_dirichlet_dofs_ ), intermediate_dirichlet_dofs_.size ( ), vec2ptr ( zeros ) );
    }
    F->UpdateCouplings ( );
}

void InstationaryFlowApp::EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
    // assemble the matrix for the linearized system
    compute_matrix ( u, DF );
    // correct BC -- set Dirichlet rows to identity
    if ( !intermediate_dirichlet_dofs_.empty ( ) )
    {
        DF->diagonalize_rows ( vec2ptr ( intermediate_dirichlet_dofs_ ), intermediate_dirichlet_dofs_.size ( ), 1. );
    }

#ifdef WITH_ILUPP
    ilupp_->SetupOperator ( intermediate_matrix_ );
#endif
}

/// Constructor: setup application

InstationaryFlowApp::InstationaryFlowApp ( )
: comm_ ( MPI_COMM_WORLD ),
rank_ ( -1 ),
num_partitions_ ( -1 ),
master_rank_ ( 0 ),
refinement_level_ ( 2 ),
nls_ ( new Newton<LAD>( &res_, &intermediate_matrix_ ) ),
#ifdef WITH_ILUPP
ilupp_ ( new PreconditionerIlupp<LAD>( ) ),
#endif
intermediate_linear_solver_ ( new GMRES<LAD>( ) ),
pressure_linear_solver_ ( new CG<LAD>( ) ),
velocity_linear_solver_ ( new GMRES<LAD>( ) )

{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    // Create directories for outputs of visualization
    bool ret = system ( "mkdir -p out; mkdir -p out_sd" );
    if ( ret != 0 ) LOG_ERROR ( "ERROR in line " << __LINE__ << ", error code = " << ret << std::endl; );

    std::string mesh_filename;
    if ( DIM == 2 )
        mesh_filename = std::string ( DATADIR ) + std::string ( "dfg_bench2d.inp" );
        //	mesh_filename = std::string(DATADIR) + std::string("channel_2d.inp");
        //	mesh_filename = std::string(DATADIR) + std::string("2d_lung_4g.inp");
        //	mesh_filename = std::string(DATADIR) + std::string("dfg_bench_2d.vtu");

    else
        mesh_filename = std::string ( DATADIR ) + std::string ( "dfg_bench3d_rect.inp" );

    read_and_distribute_mesh ( mesh_filename );

    VtkWriter writer;
    std::stringstream sstr;
    sstr << "out_sd/instationary_flow_mesh_sd_" << rank_ << ".vtu";
    writer.write ( sstr.str ( ).c_str ( ), *mesh_ );

    MeshPtr bdy_mesh = MeshPtr ( mesh_->extract_boundary_mesh ( ) );
    VtkWriter bdy_writer;
    bdy_writer.add_attribute ( "_mesh_facet_index_", DIM - 1 );
    std::stringstream bdy_sstr;
    bdy_sstr << "out_sd/instationary_flow_bdy_sd_" << rank_ << ".vtu";
    bdy_writer.write ( bdy_sstr.str ( ).c_str ( ), *bdy_mesh );

    // setup linear algebra platform
    la_sys_.Platform = APP_PLATFORM;
    init_platform ( la_sys_ );

    // setup names for visualization
    visualization_names_.push_back ( "u" );
    visualization_names_.push_back ( "v" );
    if ( DIM == 3 ) visualization_names_.push_back ( "w" );
    visualization_names_.push_back ( "p" );

    // setup instationary Navier Stokes weak form
    /// Density rho
    rho_ = 1.;

    /// Kinematic viscosity nu
    nu_ = 1.e-3;

    intermediate_asm_ = new IntermediateVelocityAssembler ( nu_, rho_ );
    pressure_asm_ = new PressureCorrectionAssembler ( rho_ );
    velocity_asm_ = new VelocityCorrectionAssembler ( rho_ );

    /// channel height [m]
    H_ = 0.41;
}

/// Destructor: clean up after application

InstationaryFlowApp::~InstationaryFlowApp ( )
{
    delete intermediate_asm_;
    delete pressure_asm_;
    delete velocity_asm_;
    stop_platform ( la_sys_ );
}

/// Mesh handling

void InstationaryFlowApp::read_and_distribute_mesh ( const std::string& filename )
{
    MeshPtr master_mesh ( 0 );

    if ( rank_ == master_rank_ )
    {
        master_mesh = read_mesh_from_file ( filename, DIM, DIM, 0 );
        assert ( master_mesh != 0 );

        // Refine mesh
        for ( int r = 0; r < refinement_level_; ++r )
        {
            master_mesh = master_mesh->refine ( );
        }
    }

    // Distribute mesh
    MeshPtr local_mesh = partition_and_distribute ( master_mesh, master_rank_, comm_ );
    assert ( local_mesh != 0 );

    // Compute ghost cells
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
}

/// Main loop of application

void InstationaryFlowApp::run ( )
{
    prepare_space ( );
    prepare_bc ( );
    prepare_lin_alg_structures ( ); // needs only be done once, since space does not change
    prepare_linear_solver ( );
    prepare_nls ( );

    // set time-stepping method
    //method_ = "ImplicitEuler";
    method_ = "CrankNicolson";
    //method_ = "ExplicitEuler";
    if ( rank_ == master_rank_ )
        std::cout << "Time-stepping method: " << method_ << "\n";
    LOG_INFO ( "Time-stepping method", method_ );

    int Tmax = 9999; // Number of time steps to be calculated
    LOG_INFO ( "Maximum number of time steps", Tmax );
    double delta_t = 0.01; // Time stepping size
    LOG_INFO ( "delta_t", delta_t );
    LOG_INFO ( "Simulated time", delta_t * Tmax );

    std::vector<double> alphas ( 5, 0.0 );
    int time_step = 0;

    if ( rank_ == master_rank_ )
    {
        std::cout << "--->  Visualizing solution at time index " << time_step << " with current time: " << delta_t * time_step << std::endl;
        std::cout << "================================================================================" << std::endl << std::endl;
    }
    // vars to find the right dofs in visualization.
    std::vector<int> var_p ( 1, DIM );
    std::vector<int> vars ( DIM );
    for ( int i = 0; i < DIM; i++ )
        vars[i] = i;
    visualize ( intermediate_sol_, intermediate_space_, time_step, "intermediate", vars );
    visualize ( pressure_sol_, pressure_space_, time_step, "pressure", var_p );
    visualize ( velocity_sol_, velocity_space_, time_step, "velocity", vars );

    // Begin time-stepping loop
    while ( time_step < Tmax )
    {
        ++time_step;
        compute_alphas ( 0, delta_t, &alphas );
        intermediate_asm_->set_timestep_parameters ( alphas );
        intermediate_asm_->set_prev_solution ( sol_prev_ );

        // intermediate velocity
        if ( rank_ == master_rank_ )
            std::cout << "\n----------------------------------------\n"
                << "Solving Viscous Burgers equation...\n"
                << "----------------------------------------\n";

        solve_nlp ( );
        visualize ( intermediate_sol_, intermediate_space_, time_step, "intermediate", vars );

        // solve pressure correction
        pressure_asm_->set_timestep ( delta_t );
        pressure_asm_->set_intermediate_velocity_solution ( intermediate_sol_ );

        if ( rank_ == master_rank_ )
            std::cout << "\n----------------------------------------\n"
                << "Solving pressure correction...\n"
                << "----------------------------------------\n";

        solve_pressure_correction ( );
        visualize ( pressure_sol_, pressure_space_, time_step, "pressure", var_p );

        // solve velocity correction
        velocity_asm_->set_timestep ( delta_t );
        velocity_asm_->set_pressure_solution ( pressure_sol_ );
        velocity_asm_->set_intermediate_velocity_solution ( intermediate_sol_ );

        if ( rank_ == master_rank_ )
            std::cout << "\n----------------------------------------\n"
                << "Solving velocity correction...\n"
                << "----------------------------------------\n";

        solve_velocity_correction ( );
        visualize ( velocity_sol_, velocity_space_, time_step, "velocity", vars );
        sol_prev_.CloneFrom ( velocity_sol_ );

        if ( rank_ == master_rank_ )
        {
            std::cout << "--->  Visualizing solution at time index " << time_step << " with current time: " << delta_t * time_step << std::endl;
            std::cout << "================================================================================" << std::endl << std::endl;
        }
        prepare_bc ( time_step * delta_t );
    }
}

/// Prepare the vector space

void InstationaryFlowApp::prepare_space ( )
{
    // setup space for intermediate velocity
    std::vector<int> intermediate_degrees ( DIM, 2 );
    intermediate_degrees.push_back ( 1 );
    intermediate_space_.Init ( intermediate_degrees, *mesh_ );

    std::vector<int> degrees ( DIM + 1, 2 );
    degrees[DIM] = 1;

    // setup space for pressure correction
    pressure_space_.Init ( degrees, *mesh_ );

    // setup space for velocity correction
    velocity_space_.Init ( degrees, *mesh_ );
}

/// Prepare the linear solvers

void InstationaryFlowApp::prepare_linear_solver ( )
{
    // default parameters: maxits, abstol, reltol, divtol
    intermediate_linear_solver_->InitControl ( 5000, 1.e-14, 1.e-12, 1.e6 );
    pressure_linear_solver_->InitControl ( 5000, 1.e-12, 1.e-12, 1.e6 );
    velocity_linear_solver_->InitControl ( 5000, 1.e-12, 1.e-12, 1.e6 );

    // default GMRES parameters: size of basis, preconditioning method
#ifdef WITH_ILUPP
    intermediate_linear_solver_->InitParameter ( 300, "RightPreconditioning" );
    //    ilupp_->InitParameter(0, 1000, 10, 0.8, 2.75, 0.01);
    ilupp_->InitParameter ( 0, 1010, 50, 0.8, 5., 0.01 );
    intermediate_linear_solver_->SetupPreconditioner ( *ilupp_ );
#else
    intermediate_linear_solver_->InitParameter ( 300, "NoPreconditioning" );
#endif
    pressure_linear_solver_->InitParameter ( "NoPreconditioning" );
    velocity_linear_solver_->InitParameter ( 300, "NoPreconditioning" );

    // set the matrix to be used as the operator
    intermediate_linear_solver_->SetupOperator ( intermediate_matrix_ );
    pressure_linear_solver_->SetupOperator ( pressure_matrix_ );
    velocity_linear_solver_->SetupOperator ( velocity_matrix_ );
}

/// Prepare the non-zero structures for the linear algebra module

void InstationaryFlowApp::prepare_lin_alg_structures ( )
{
    std::vector<int> diagonal_rows, diagonal_cols, off_diagonal_rows, off_diagonal_cols;
    // Initialize linear algebra structures

    std::vector < std::vector<bool> > coupling_vars;
    coupling_vars.resize ( DIM + 1 );
    for ( int i = 0; i < DIM; ++i )
    {
        for ( int j = 0; j < DIM; ++j )
        {
            coupling_vars[i].push_back ( true );
        }
        coupling_vars[i].push_back ( false );
    }
    for ( int i = 0; i < DIM; ++i )
    {
        coupling_vars[DIM].push_back ( false );
    }
    coupling_vars[DIM].push_back ( true );

    // Intermediate solution
    SparsityStructure intermediate_sparsity;
    global_asm_.compute_sparsity_structure ( intermediate_space_, intermediate_sparsity, &coupling_vars );

    intermediate_couplings_.Clear ( );
    intermediate_couplings_.Init ( comm_, intermediate_space_.dof ( ) );

    intermediate_couplings_.InitializeCouplings ( intermediate_sparsity.off_diagonal_rows,
                                                  intermediate_sparsity.off_diagonal_cols );

    intermediate_matrix_.Init ( comm_, intermediate_couplings_, la_sys_.Platform,
                                APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    intermediate_matrix_.InitStructure ( vec2ptr ( intermediate_sparsity.diagonal_rows ), vec2ptr ( intermediate_sparsity.diagonal_cols ),
                                         intermediate_sparsity.diagonal_rows.size ( ), vec2ptr ( intermediate_sparsity.off_diagonal_rows ),
                                         vec2ptr ( intermediate_sparsity.off_diagonal_cols ), intermediate_sparsity.off_diagonal_rows.size ( ) );

    intermediate_sol_.Init ( comm_, intermediate_couplings_, la_sys_.Platform,
                             APP_LINALG_IMPLEMENTATION );
    intermediate_sol_.InitStructure ( );
    intermediate_sol_.Zeros ( );

    // correct solution with dirichlet BC for visualization
    if ( !intermediate_dirichlet_dofs_.empty ( ) )
    {
        intermediate_sol_.SetValues ( vec2ptr ( intermediate_dirichlet_dofs_ ),
                                      intermediate_dirichlet_dofs_.size ( ),
                                      vec2ptr ( intermediate_dirichlet_values_ ) );
    }
    intermediate_sol_.UpdateCouplings ( );

    sol_prev_.Init ( comm_, intermediate_couplings_, la_sys_.Platform,
                     APP_LINALG_IMPLEMENTATION );
    sol_prev_.InitStructure ( );
    sol_prev_.Zeros ( );

    intermediate_rhs_.Init ( comm_, intermediate_couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    intermediate_rhs_.InitStructure ( );
    intermediate_rhs_.Zeros ( );

    // Pressure correction
    SparsityStructure pressure_sparsity;
    global_asm_.compute_sparsity_structure ( pressure_space_, pressure_sparsity, &coupling_vars );

    pressure_couplings_.Clear ( );
    pressure_couplings_.Init ( comm_, pressure_space_.dof ( ) );

    pressure_couplings_.InitializeCouplings ( pressure_sparsity.off_diagonal_rows,
                                              pressure_sparsity.off_diagonal_cols );

    pressure_matrix_.Init ( comm_, pressure_couplings_, la_sys_.Platform,
                            APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    pressure_matrix_.InitStructure ( vec2ptr ( pressure_sparsity.diagonal_rows ), vec2ptr ( pressure_sparsity.diagonal_cols ),
                                     pressure_sparsity.diagonal_rows.size ( ), vec2ptr ( pressure_sparsity.off_diagonal_rows ),
                                     vec2ptr ( pressure_sparsity.off_diagonal_cols ), pressure_sparsity.off_diagonal_rows.size ( ) );

    pressure_sol_.Init ( comm_, pressure_couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    pressure_sol_.InitStructure ( );
    pressure_sol_.Zeros ( );

    // correct solution with dirichlet BC for visualization
    if ( !pressure_dirichlet_dofs_.empty ( ) )
    {
        pressure_sol_.SetValues ( vec2ptr ( pressure_dirichlet_dofs_ ),
                                  pressure_dirichlet_dofs_.size ( ),
                                  vec2ptr ( pressure_dirichlet_values_ ) );
    }
    pressure_sol_.UpdateCouplings ( );

    pressure_rhs_.Init ( comm_, pressure_couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    pressure_rhs_.InitStructure ( );
    pressure_rhs_.Zeros ( );

    // Velocity correction
    SparsityStructure velocity_sparsity;
    global_asm_.compute_sparsity_structure ( velocity_space_, velocity_sparsity, &coupling_vars );

    velocity_couplings_.Clear ( );
    velocity_couplings_.Init ( comm_, velocity_space_.dof ( ) );

    velocity_couplings_.InitializeCouplings ( velocity_sparsity.off_diagonal_rows,
                                              velocity_sparsity.off_diagonal_cols );

    velocity_matrix_.Init ( comm_, velocity_couplings_, la_sys_.Platform,
                            APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    velocity_matrix_.InitStructure ( vec2ptr ( velocity_sparsity.diagonal_rows ), vec2ptr ( velocity_sparsity.diagonal_cols ),
                                     velocity_sparsity.diagonal_rows.size ( ), vec2ptr ( velocity_sparsity.off_diagonal_rows ),
                                     vec2ptr ( velocity_sparsity.off_diagonal_cols ), velocity_sparsity.off_diagonal_rows.size ( ) );

    velocity_sol_.Init ( comm_, velocity_couplings_, la_sys_.Platform,
                         APP_LINALG_IMPLEMENTATION );
    velocity_sol_.InitStructure ( );
    velocity_sol_.Zeros ( );

    velocity_rhs_.Init ( comm_, velocity_couplings_, la_sys_.Platform,
                         APP_LINALG_IMPLEMENTATION );
    velocity_rhs_.InitStructure ( );
    velocity_rhs_.Zeros ( );
}

/// Prepare the non-linear solver

void InstationaryFlowApp::prepare_nls ( )
{
    // default parameters: maxits, abstol, reltol, divtol
    nls_->InitControl ( 1000, 1.e-9, 1.e-9, 1.e6 );

    nls_->SetLinearSolver ( *intermediate_linear_solver_ );

    nls_->SetOperator ( *this );
    // we use our own initial solution -- this needs to be indicated
    // to the Newton-solver
    nls_->InitParameter ( Newton<LAD>::NewtonInitialSolutionOwn );
}

/// Solve non-linear problem

void InstationaryFlowApp::solve_nlp ( )
{
    // intermediate velocity
    intermediate_rhs_.UpdateCouplings ( );
    intermediate_sol_.UpdateCouplings ( );

    ArmijoDamping<LAD> armijo ( 1., 1.e-6, 0.5, 1.e-6, 15 );
    nls_->SetDampingStrategy ( armijo );

    NonlinearSolverState state = nls_->Solve ( intermediate_rhs_, &intermediate_sol_ );
    intermediate_rhs_.UpdateCouplings ( );
    intermediate_sol_.UpdateCouplings ( );
    if ( rank_ == master_rank_ )
        std::cout << "Nonlinear solver ended with state " << state
            << " and residual norm " << nls_->GetResidual ( )
        << " after " << nls_->iter ( ) << " iterations\n";
    if ( state != 0 )
    {
        std::cerr << "Program terminated because nonlinear solver did not converge!\n";
        exit ( -1 );
    }
}

/// Solve poisson problem

void InstationaryFlowApp::solve_pressure_correction ( )
{
    // pressure correction

    global_asm_.assemble_vector ( pressure_space_, *pressure_asm_, pressure_rhs_ );
    global_asm_.assemble_matrix ( pressure_space_, *pressure_asm_, pressure_matrix_ );

    // correct BC -- set Dirichlet rows to identity
    if ( !pressure_dirichlet_dofs_.empty ( ) )
    {
        pressure_matrix_.diagonalize_rows ( vec2ptr ( pressure_dirichlet_dofs_ ), pressure_dirichlet_dofs_.size ( ), 1. );
        pressure_rhs_.SetValues ( vec2ptr ( pressure_dirichlet_dofs_ ), pressure_dirichlet_dofs_.size ( ), vec2ptr ( pressure_dirichlet_values_ ) );
        pressure_sol_.SetValues ( vec2ptr ( pressure_dirichlet_dofs_ ), pressure_dirichlet_dofs_.size ( ), vec2ptr ( pressure_dirichlet_values_ ) );
    }
    pressure_rhs_.UpdateCouplings ( );
    pressure_sol_.UpdateCouplings ( );

    //    pressure_linear_solver_->SetupOperator(pressure_matrix_);
    pressure_linear_solver_->Solve ( pressure_rhs_, &pressure_sol_ );

    pressure_sol_.UpdateCouplings ( );
}

/// Project velocity on subspace of solenoidal functions

void InstationaryFlowApp::solve_velocity_correction ( )
{
    // velocity correction
    StandardGlobalAssembler<double> global_asm;
    global_asm.assemble_vector ( velocity_space_, *velocity_asm_, velocity_rhs_ );
    global_asm.assemble_matrix ( velocity_space_, *velocity_asm_, velocity_matrix_ );

    // correct BC -- set Dirichlet rows to identity
    if ( !velocity_dirichlet_dofs_.empty ( ) )
    {
        velocity_matrix_.diagonalize_rows ( vec2ptr ( velocity_dirichlet_dofs_ ), velocity_dirichlet_dofs_.size ( ), 1. );
        // no slip at upper and lower boundary
        velocity_rhs_.SetValues ( vec2ptr ( velocity_dirichlet_dofs_ ), velocity_dirichlet_dofs_.size ( ), vec2ptr ( velocity_dirichlet_values_ ) );
        velocity_sol_.SetValues ( vec2ptr ( velocity_dirichlet_dofs_ ), velocity_dirichlet_dofs_.size ( ), vec2ptr ( velocity_dirichlet_values_ ) );
    }
    velocity_rhs_.UpdateCouplings ( );
    velocity_sol_.UpdateCouplings ( );

    //    velocity_linear_solver_->SetupOperator(velocity_matrix_);
    velocity_linear_solver_->Solve ( velocity_rhs_, &velocity_sol_ );

    velocity_sol_.UpdateCouplings ( );
}

struct PressureDirichletBC2d
{
    // Material numbers:
    // 10 -> left edge (inflow)
    // 11 -> bottom edge (u = 0)
    // 12 -> right edge (outflow)
    // 13 -> top edge (u = 0)
    // 14 -> cylinder (u = 0)

    PressureDirichletBC2d ( int var, double H, double time )
    : var_ ( var ), H_ ( H ), time_ ( time )
    {
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );

        const bool inflow = ( material_num == 10 );
        const bool outflow = ( material_num == 12 );
        //	const bool cylinder = (material_num == 14);

        if ( outflow )
        {
            values.resize ( coords_on_face.size ( ) );
            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                values[i] = 0.;
            }
        }
        /*        if (cylinder) {
                  values.resize(coords_on_face.size());
                  // loop over points on the face
                  for (int i = 0; i < coords_on_face.size(); ++i) {
                  values[i] = 0.;
                  }
                  }
         */
        if ( inflow )
        {
            values.resize ( coords_on_face.size ( ) );
            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                values[i] = 1.; //sin(3.0 * time_);
                //const Coord& pt = coords_on_face[i];
                //values[i] = 4. * 0.64 * pt[1] * (H_ - pt[1]) / (H_ * H_);

            }
        }
        return values;
    }

    const int var_;
    const double H_;
    const double time_;
};

struct PressureDirichletBC3d
{
    // Material numbers:
    // 10 -> left edge (inflow)
    // 11 -> bottom edge (u = 0)
    // 12 -> right edge (outflow)
    // 13 -> top edge (u = 0)
    // 14 -> cylinder (u = 0)

    PressureDirichletBC3d ( int var, double H, double time )
    : var_ ( var ), H_ ( H ), time_ ( time )
    {
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );

        const bool inflow = ( material_num == 10 );
        const bool outflow = ( material_num == 12 );
        //	const bool cylinder = (material_num == 14);

        if ( outflow )
        {
            values.resize ( coords_on_face.size ( ) );
            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                values[i] = 0.;
            }
        }
        /*        if (cylinder) {
                  values.resize(coords_on_face.size());
                  // loop over points on the face
                  for (int i = 0; i < coords_on_face.size(); ++i) {
                  values[i] = 0.;
                  }
                  }
         */
        if ( inflow )
        {
            values.resize ( coords_on_face.size ( ) );
            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                values[i] = 1.; //sin(3.0 * time_);

                // 		const Coord& pt = coords_on_face[i];
                // 		const double y = pt[1];
                // 		const double z = pt[2];
                // 		values[i] = 16. * Um_ * y * z * (H_ - y) * (H_ - z) / (std::pow(H_, 4.));
            }
        }
        return values;
    }

    const int var_;
    const double H_;
    const double time_;
};

struct NoSlipBoundaryCondition
{
    // Material numbers:
    // 10 -> left edge (inflow)
    // 11 -> bottom edge (u = 0)
    // 12 -> right edge (outflow)
    // 13 -> top edge (u = 0)
    // 14 -> cylinder (u = 0)

    NoSlipBoundaryCondition ( )
    {
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );

        const bool inflow = ( material_num == 10 );
        const bool outflow = ( material_num == 12 );
        //	const bool cylinder = (material_num == 14);

        if ( !outflow && !inflow )
        {
            values.resize ( coords_on_face.size ( ) );
            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                values[i] = 0.;
            }
        }
        return values;
    }
};

void InstationaryFlowApp::prepare_bc ( double time )
{
    intermediate_dirichlet_dofs_.clear ( );
    intermediate_dirichlet_values_.clear ( );

    pressure_dirichlet_dofs_.clear ( );
    pressure_dirichlet_values_.clear ( );

    velocity_dirichlet_dofs_.clear ( );
    velocity_dirichlet_values_.clear ( );

    if ( DIM == 2 )
    {
        // pressure boundary conditions
        PressureDirichletBC2d pressure_bc = PressureDirichletBC2d ( DIM, H_, time );

        compute_dirichlet_dofs_and_values ( pressure_bc, pressure_space_, DIM,
                                            pressure_dirichlet_dofs_,
                                            pressure_dirichlet_values_ );

        // no slip boundary condition for channel walls
        NoSlipBoundaryCondition no_slip_bc[2] = { NoSlipBoundaryCondition ( ),
                                                 NoSlipBoundaryCondition ( ) };
        for ( int var = 0; var < DIM; ++var )
        {
            compute_dirichlet_dofs_and_values ( no_slip_bc[var], velocity_space_, var,
                                                velocity_dirichlet_dofs_,
                                                velocity_dirichlet_values_ );
            compute_dirichlet_dofs_and_values ( no_slip_bc[var], intermediate_space_, var,
                                                intermediate_dirichlet_dofs_,
                                                intermediate_dirichlet_values_ );
        }

    }
    else
    {
        assert ( DIM == 3 );
        PressureDirichletBC3d pressure_bc3d = PressureDirichletBC3d ( DIM, H_, time );

        compute_dirichlet_dofs_and_values ( pressure_bc3d, pressure_space_, DIM,
                                            pressure_dirichlet_dofs_,
                                            pressure_dirichlet_values_ );

        // no slip boundary condition for channel walls
        NoSlipBoundaryCondition no_slip_bc3d[3] = { NoSlipBoundaryCondition ( ),
                                                   NoSlipBoundaryCondition ( ),
                                                   NoSlipBoundaryCondition ( ) };
        for ( int var = 0; var < DIM; ++var )
        {
            compute_dirichlet_dofs_and_values ( no_slip_bc3d[var], velocity_space_, var,
                                                velocity_dirichlet_dofs_,
                                                velocity_dirichlet_values_ );
            compute_dirichlet_dofs_and_values ( no_slip_bc3d[var], intermediate_space_, var,
                                                intermediate_dirichlet_dofs_,
                                                intermediate_dirichlet_values_ );
        }
    }
}

/// Visualization of the solution

void InstationaryFlowApp::visualize ( LAD::VectorType& sol, VectorSpace<double>& space, int step, std::string name, std::vector<int> vars )
{

    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space, num_intervals, comm_, master_rank_ );

    std::stringstream input;

    if ( step < 10 )
        input << "instationary_flow_values_" << name << "_000" << step;
    else if ( step < 100 )
        input << "instationary_flow_values_" << name << "_00" << step;
    else if ( step < 1000 )
        input << "instationary_flow_values_" << name << "_0" << step;
    else
        input << "instationary_flow_values_" << name << "_" << step;

    std::vector<double> material_number ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
          it != mesh_->end ( mesh_->tdim ( ) );
          ++it )
    {
        material_number.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
    }

    sol.UpdateCouplings ( );

    for ( int i = 0; i < static_cast < int > ( vars.size ( ) ); ++i )
    {
        visu.visualize ( EvalFeFunction<LAD>( space, sol, vars[i] ), visualization_names_[vars[i]] );
    }

    visu.visualize_cell_data ( material_number, "Material Id" );

    visu.write ( "out/" + input.str ( ) );
}

/// Set parameters for the different time-stepping methods

void InstationaryFlowApp::compute_alphas ( int sub_step, double delta_t, std::vector<double>* alphas )
{
    assert ( alphas->size ( ) == 5 );
    if ( method_ == "ImplicitEuler" )
    {
        alphas->at ( 0 ) = delta_t;
        alphas->at ( 1 ) = delta_t;
        alphas->at ( 2 ) = 0.;
        alphas->at ( 3 ) = 0.;
        alphas->at ( 4 ) = delta_t;
    }
    else if ( method_ == "CrankNicolson" )
    {
        alphas->at ( 0 ) = 0.5 * delta_t;
        alphas->at ( 1 ) = delta_t;
        alphas->at ( 2 ) = 0.5 * delta_t;
        alphas->at ( 3 ) = 0.5 * delta_t;
        alphas->at ( 4 ) = 0.5 * delta_t;
    }
    else if ( method_ == "ExplicitEuler" )
    {
        alphas->at ( 0 ) = 0.;
        alphas->at ( 1 ) = delta_t;
        alphas->at ( 2 ) = delta_t;
        alphas->at ( 3 ) = delta_t;
        alphas->at ( 4 ) = 0.;
    }
    else
        assert ( 0 );
}

/// Compute the residual for non-linear solver

void InstationaryFlowApp::compute_residual ( const LAD::VectorType& u, LAD::VectorType* F )
{
    intermediate_asm_->set_newton_solution ( u );
    global_asm_.assemble_vector ( intermediate_space_, *intermediate_asm_, *F );
    F->UpdateCouplings ( );
}

/// Compute the matrix for the non-linear solver

void InstationaryFlowApp::compute_matrix ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
    intermediate_asm_->set_newton_solution ( u );
    global_asm_.assemble_matrix ( intermediate_space_, *intermediate_asm_, *DF );
}

//////////////////////////////////////////////////////////////////////
// IntermediateVelocityAssembler to solve Viscous Burgers equation
//////////////////////////////////////////////////////////////////////

void IntermediateVelocityAssembler::set_newton_solution ( const LAD::VectorType& newton_sol )
{
    newton_sol_ = &newton_sol;
}

void IntermediateVelocityAssembler::set_prev_solution ( const LAD::VectorType& prev_sol )
{
    prev_sol_ = &prev_sol;
}

void IntermediateVelocityAssembler::set_timestep_parameters ( const std::vector<double>& alphas )
{
    assert ( alphas.size ( ) == 5 );
    alpha1_ = alphas[0];
    alpha2_ = alphas[1];
    alpha3_ = alphas[2];
    alpha4_ = alphas[3];
    alpha5_ = alphas[4];
}

void IntermediateVelocityAssembler::initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    // recompute previous solution values
    for ( int v = 0; v < DIM; ++v )
    {
        prev_ns_vel_[v].clear ( );
        prev_ts_vel_[v].clear ( );
        grad_prev_ns_vel_[v].clear ( );
        grad_prev_ts_vel_[v].clear ( );
        evaluate_fe_function ( *newton_sol_, v, prev_ns_vel_[v] );
        evaluate_fe_function ( *prev_sol_, v, prev_ts_vel_[v] );
        evaluate_fe_function_gradients ( *newton_sol_, v, grad_prev_ns_vel_[v] );
        evaluate_fe_function_gradients ( *prev_sol_, v, grad_prev_ts_vel_[v] );
    }
}

void IntermediateVelocityAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
        LocalMatrix& lm )
{
    initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    // loop q
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );

        // get previous solution in vector form
        Vec<DIM, double> vel_k;
        for ( int var = 0; var < DIM; ++var )
        {
            vel_k[var] = prev_ns_vel_[var][q];
        }

        // assemble a0(u,v) = \int {dot(u,v)}
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < num_dofs ( u_var ); ++j )
                {
                    lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                            wq * ( phi ( j, q, u_var ) * phi ( i, q, u_var ) ) * dJ;
                }
            }
        }

        // assemble a1(u,v) = \int \alpha_1 * {\grad(u) : \grad(v)}
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < num_dofs ( u_var ); ++j )
                {
                    lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                            wq * alpha1_ * ( nu_ * dot ( grad_phi ( j, q, u_var ), grad_phi ( i, q, u_var ) ) ) * dJ;
                }
            }
        }

        // assemble a2(u,v) = \int \alpha_1 * { (vel_k*\grad{u})*v }
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < num_dofs ( u_var ); ++j )
                {
                    lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                            wq * alpha1_ * ( dot ( vel_k, grad_phi ( j, q, u_var ) ) * phi ( i, q, u_var ) ) * dJ;
                }
            }
        }

        // assemble a3(u,v) = \int \alpha_1 * { (u\grad{u_k}*v }
        for ( int test_var = 0; test_var < DIM; ++test_var )
        {
            for ( int trial_var = 0; trial_var < DIM; ++trial_var )
            {
                for ( int i = 0; i < num_dofs ( test_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( trial_var ); ++j )
                    {
                        lm ( dof_index ( i, test_var ), dof_index ( j, trial_var ) ) +=
                                wq * alpha1_ * ( grad_prev_ns_vel_[test_var][q][trial_var] *
                                phi ( j, q, trial_var ) *
                                phi ( i, q, test_var ) ) * dJ;
                    }
                }
            }
        }

        // assemble mass matrix for pressure component
        for ( int i = 0; i < num_dofs ( DIM ); ++i )
        {
            for ( int j = 0; j < num_dofs ( DIM ); ++j )
            {
                lm ( dof_index ( i, DIM ), dof_index ( j, DIM ) ) +=
                        wq * phi ( j, q, DIM ) * phi ( i, q, DIM ) * dJ;
            }
        }
    }
}

void IntermediateVelocityAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
{
    initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );

        // get previous newton step solution in vector form
        Vec<DIM, double> vel_k;
        for ( int var = 0; var < DIM; ++var )
        {
            vel_k[var] = prev_ns_vel_[var][q];
        }

        // get previous time step solution in vector form
        Vec<DIM, double> vel_n;
        for ( int var = 0; var < DIM; ++var )
        {
            vel_n[var] = prev_ts_vel_[var][q];
        }

        // l0(v) = \int( dot(u_n - u_k, v))
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] +=
                        wq * ( ( vel_n[v_var] - vel_k[v_var] ) * phi ( i, q, v_var ) ) * dJ;
            }
        }

        // l1n(v) = -\alpha_3 * \nu * \int( \grad{u_n} : \grad{v} )
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] +=
                        -wq * alpha3_ * ( nu_ * dot ( grad_phi ( i, q, v_var ), grad_prev_ts_vel_[v_var][q] ) ) * dJ;
            }
        }

        // l1k(v) = -\alpha_1 * \nu * \int( \grad{u_k} : \grad{v} )
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] +=
                        -wq * alpha1_ * ( nu_ * dot ( grad_phi ( i, q, v_var ), grad_prev_ns_vel_[v_var][q] ) ) * dJ;
            }
        }

        // l2n(v) = -\alpha_3 * \int(u_n*\grad{u_n}*v)
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] +=
                        -wq * alpha3_ * ( dot ( grad_prev_ts_vel_[v_var][q], vel_n ) * phi ( i, q, v_var ) ) * dJ;
            }
        }

        // l2k(v) = -\alpha_1 * \int(u_k*\grad{u_k}*v)
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] +=
                        -wq * alpha1_ * ( dot ( grad_prev_ns_vel_[v_var][q], vel_k ) * phi ( i, q, v_var ) ) * dJ;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// PressureCorrectionAssembler to solve Poisson equation with div(u)
// on the right hand side
//////////////////////////////////////////////////////////////////////

void PressureCorrectionAssembler::set_intermediate_velocity_solution ( const LAD::VectorType& intermediate_sol )
{
    intermediate_sol_ = &intermediate_sol;
}

void PressureCorrectionAssembler::set_timestep ( double delta_t )
{
    delta_t_ = delta_t;
}

void PressureCorrectionAssembler::initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    // recompute intermediate velocity solution values
    for ( int v = 0; v < DIM; ++v )
    {
        intermediate_velocity_gradient_[v].clear ( );
        evaluate_fe_function_gradients ( *intermediate_sol_, v, intermediate_velocity_gradient_[v] );
    }
}

void PressureCorrectionAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
{
    initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    // matrix set up with only one variable
    const int pressure_dofs = num_dofs ( DIM );
    const int total_dofs = num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    // loop q
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );

        // solve laplace equation for pressure with known intermediate
        // velocity
        //
        // assemble a(p,q) = \int {dot(p, q)}
        for ( int i = 0; i < pressure_dofs; ++i )
        {
            for ( int j = 0; j < pressure_dofs; ++j )
            {
                lm ( dof_index ( i, DIM ), dof_index ( j, DIM ) ) +=
                        wq * dot ( grad_phi ( j, q, DIM ), grad_phi ( i, q, DIM ) ) * dJ;
            }
        }
    }
}

void PressureCorrectionAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
{
    initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    // matrix set up with only one variable
    const int pressure_dofs = num_dofs ( DIM );
    const int total_dofs = num_dofs_total ( );

    lv.resize ( total_dofs, 0. );

    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = detJ ( q );
        double div_u0 = 0.0;
        // calculate div(u_intermediate)
        for ( int i = 0; i < DIM; ++i )
        {
            div_u0 += intermediate_velocity_gradient_[i][q][i];
        }

        // l(q) = rho/dt * \int{div(u_inter) * qx}
        for ( int i = 0; i < pressure_dofs; ++i )
        {
            lv[dof_index ( i, DIM )] += -wq * rho_ / delta_t_ * div_u0 * phi ( i, q, DIM ) * dJ;
        }
    }
}

//////////////////////////////////////////////////////////////////////
// VelocityCorrectionAssembler to project intermediate velocity on the
// subspace of solenoidal functions
//////////////////////////////////////////////////////////////////////

void VelocityCorrectionAssembler::set_intermediate_velocity_solution ( const LAD::VectorType& intermediate_sol )
{
    intermediate_sol_ = &intermediate_sol;
}

void VelocityCorrectionAssembler::set_pressure_solution ( const LAD::VectorType& pressure_sol )
{
    pressure_sol_ = &pressure_sol;
}

void VelocityCorrectionAssembler::set_timestep ( double delta_t )
{
    delta_t_ = delta_t;
}

void VelocityCorrectionAssembler::initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    // recompute intermediate velocity solution values
    for ( int v = 0; v < DIM; ++v )
    {
        intermediate_velocity_[v].clear ( );
        evaluate_fe_function ( *intermediate_sol_, v, intermediate_velocity_[v] );
    }

    // recompute pressure solution values
    pressure_gradient_.clear ( );
    evaluate_fe_function_gradients ( *pressure_sol_, DIM, pressure_gradient_ );
}

void VelocityCorrectionAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
{
    initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    int velocity_dofs = 0;
    for ( int i = 0; i < DIM; ++i )
        velocity_dofs += num_dofs ( i );

    const int total_dofs = num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    // loop q
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );

        // assemble a(u,v) = \int {dot(u, v)}
        for ( int var = 0; var < DIM; ++var )
        {
            for ( int i = 0; i < num_dofs ( var ); ++i )
            {
                for ( int j = 0; j < num_dofs ( var ); ++j )
                {
                    lm ( dof_index ( i, var ), dof_index ( j, var ) ) +=
                            wq * phi ( j, q, var ) * phi ( i, q, var ) * dJ;
                }
            }
        }
    }
}

void VelocityCorrectionAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
{
    initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    int velocity_dofs = 0;
    for ( int i = 0; i < DIM; ++i )
        velocity_dofs += num_dofs ( i );

    const int total_dofs = num_dofs_total ( );
    lv.resize ( total_dofs, 0. );

    // loop q
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );

        // assemble l(v) = \int {dot(u_inter, v) - dt/rho * dot(grad(p), v)}
        for ( int var = 0; var < DIM; ++var )
        {
            for ( int i = 0; i < num_dofs ( var ); ++i )
            {
                lv[dof_index ( i, var )] +=
                        wq * ( intermediate_velocity_[var][q] * phi ( i, q, var )
                        - delta_t_ / rho_ * pressure_gradient_[q][var] * phi ( i, q, var ) ) * dJ;
            }
        }
    }
}
