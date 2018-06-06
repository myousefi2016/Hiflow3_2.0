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

#include "instationary_flow_benchmark.h"

static const char* DATADIR = MESHES_DATADIR;

/// Program main loop: setup MPI, read parameters and start the application

int main ( int argc, char** argv )
{

    MPI_Init ( &argc, &argv );

    // Read parameters

    // Run application
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

    // correct BC -- set Dirichlet dofs to 0
    if ( !dirichlet_dofs_.empty ( ) )
    {
        std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
        F->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( zeros ) );
    }

}

void InstationaryFlowApp::EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
    // assemble the matrix for the linearized system
    compute_matrix ( u, DF );

    // correct BC -- set Dirichlet rows to identity
    if ( !dirichlet_dofs_.empty ( ) )
    {
        DF->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
    }

#ifdef WITH_MUMPS
    linear_solver_->SetupOperator ( matrix_ );
    //    utils::Timer  timer;
    //    timer.start();
    linear_solver_->Build ( );
    //    timer.stop();
    //    std::cout << "Time for factorizing matrix: " << timer.get_duration() << std::endl;
#endif

}

/// Constructor: setup application

InstationaryFlowApp::InstationaryFlowApp ( )
: comm_ ( MPI_COMM_WORLD ),
rank_ ( -1 ),
num_partitions_ ( -1 ),
master_rank_ ( 0 ),
visu_p2s_ ( false ),
nls_ ( new Newton<LAD>( &res_, &matrix_ ) ),
#ifdef WITH_MUMPS
linear_solver_ ( new MumpsSolver<LAD, MumpsStructureD>( ) ),
#else
linear_solver_ ( new GMRES<LAD>( ) ),
#endif
refinement_level_ ( 3 )
{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    // Create directories for outputs of visualisation
    bool ret = system ( "mkdir -p out; mkdir -p out_seq; mkdir -p out_sd;" );
    if ( ret != 0 ) LOG_ERROR ( "ERROR in line " << __LINE__ << ", error code = " << ret << std::endl; );

    std::string mesh_filename;
    if ( DIM == 2 )
        mesh_filename = std::string ( DATADIR ) + std::string ( "dfg_bench2d.inp" );
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

    // setup Instationary NavierStokes weak form
    rho_ = 1.; // density
    nu_ = 1.e-3; // kinematic viscosity

    local_asm_ = new InstationaryFlowAssembler ( nu_, rho_ );

    H_ = 0.41; // channel height [m]
    if ( DIM == 2 )
        Um_ = 1.5; // mean inflow speed [m/s]
    else
        Um_ = 2.25; // mean inflow speed [m/s]

}

/// Destructor: clean up after application

InstationaryFlowApp::~InstationaryFlowApp ( )
{
    delete local_asm_;
    stop_platform ( la_sys_ );
}

void InstationaryFlowApp::read_and_distribute_mesh ( const std::string& filename )
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

    // Distribute mesh using Naive partitioner
    MeshPtr local_mesh = partition_and_distribute ( master_mesh, master_rank_,
                                                    comm_ );
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

    // setup time stepping method
    //method_ = "ImplicitEuler";
    //method_ = "FractionalStepThetaScheme";
    method_ = "CrankNicolson";

    int Tmax = 50; // Number of time steps to be calculated
    double delta_t = 0.1; // Time stepping size

    std::vector<double> alphas ( 5, 0.0 );
    int time_step = 0;

    if ( rank_ == master_rank_ )
        std::cout << "Solving instationary Navier Stokes with " << method_ << std::endl;

    if ( rank_ == master_rank_ )
    {
        std::cout << "--->  Visualizing solution at time index " << time_step << " with current time: " << delta_t * time_step << std::endl;
        std::cout << "==============================" << std::endl << std::endl;
    }

    visualize_solution ( 0 );

    // Beginn Time Stepping Loop

    while ( time_step < Tmax )
    {
        if ( method_ == "FractionalStepThetaScheme" )
        {

            int sub_step = 0;

            while ( sub_step < 3 )
            {
                if ( rank_ == master_rank_ )
                {
                    std::cout << "Calculating sub step #" << sub_step << std::endl;
                    std::cout << "----------------------" << std::endl;
                }

                compute_alphas ( sub_step, delta_t, &alphas );

                local_asm_->set_timestep_parameters ( alphas );
                local_asm_->set_prev_solution ( sol_prev_ );

                solve_nlp ( );

                sol_prev_.CloneFrom ( sol_ );

                ++sub_step;

                if ( rank_ == master_rank_ )
                    std::cout << "----------------------" << std::endl;
            }
        }
        else
        {
            compute_alphas ( 0, delta_t, &alphas );

            local_asm_->set_timestep_parameters ( alphas );
            local_asm_->set_prev_solution ( sol_prev_ );

            solve_nlp ( );

            sol_prev_.CloneFrom ( sol_ );
        }

        ++time_step;

        if ( rank_ == master_rank_ )
        {
            std::cout << "--->  Visualizing solution at time index " << time_step << " with current time: " << delta_t * time_step << std::endl;
            std::cout << "==============================" << std::endl << std::endl;
        }

        visualize_solution ( time_step );
    }

#ifdef WITH_MUMPS
    linear_solver_->Clear ( );
#endif
}

void InstationaryFlowApp::prepare_space ( )
{
    // setup space Q2-Q1 elements
    std::vector<int> degrees ( DIM + 1, 2 );
    degrees[DIM] = 1;
    space_.Init ( degrees, *mesh_ );
}

void InstationaryFlowApp::prepare_linear_solver ( )
{
#ifdef WITH_MUMPS
    // linear_solver_->InitParameter(MPI_COMM_WORLD);
    linear_solver_->InitParameter ( comm_ );
#else
    // default parameters: maxits, abstol, reltol, divtol
    linear_solver_->InitControl ( 5000, 1.e-12, 1.e-12, 1.e6 );

    // default GMRES parameter: size of basis, preconditioning method
    linear_solver_->InitParameter ( 500, "NoPreconditioning" );

    // set the matrix to be used as the operator
    linear_solver_->SetupOperator ( matrix_ );
#endif
}

void InstationaryFlowApp::prepare_lin_alg_structures ( )
{
    // Initialize linear algebra structures
    couplings_.Init ( comm_, space_.dof ( ) );

    // compute matrix graph to build matrix structure
    std::vector<int> diagonal_rows, diagonal_cols, off_diagonal_rows, off_diagonal_cols;
    std::vector < std::vector<bool> > coupling_vars;
    coupling_vars.resize ( DIM + 1 );
    for ( int i = 0; i < DIM; ++i )
    {
        for ( int j = 0; j < DIM + 1; ++j )
        {
            coupling_vars[i].push_back ( true );
        }
    }
    for ( int i = 0; i < DIM; ++i )
    {
        coupling_vars[DIM].push_back ( true );
    }
    coupling_vars[DIM].push_back ( false );

    InitStructure ( space_, &diagonal_rows, &diagonal_cols,
                    &off_diagonal_rows, &off_diagonal_cols, &coupling_vars );
    couplings_.InitializeCouplings ( off_diagonal_rows, off_diagonal_cols );

    // Initialize matrices and vectors
    matrix_.Init ( comm_, couplings_, la_sys_.Platform,
                   InstationaryFlowApp::APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    matrix_.InitStructure ( vec2ptr ( diagonal_rows ),
                            vec2ptr ( diagonal_cols ),
                            diagonal_rows.size ( ),
                            vec2ptr ( off_diagonal_rows ),
                            vec2ptr ( off_diagonal_cols ),
                            off_diagonal_rows.size ( ) );

    sol_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    sol_.InitStructure ( );
    sol_.Zeros ( );

    sol_prev_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    sol_prev_.InitStructure ( );
    sol_prev_.Zeros ( );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
        //sol_prev_.SetValues(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(), vec2ptr(dirichlet_values_));
    }

    rhs_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    rhs_.InitStructure ( );
    rhs_.Zeros ( );
}

void InstationaryFlowApp::prepare_nls ( )
{
    // default parameters: maxits, abstol, reltol, divtol
    nls_->InitControl ( 100, 1.e-9, 1.e-9, 1.e6 );

    nls_->SetLinearSolver ( *linear_solver_ );

    nls_->SetOperator ( *this );
    // we use our own initial solution -- this needs to be indicated
    // to the Newton-solver
    nls_->InitParameter ( Newton<LAD>::NewtonInitialSolutionOwn );
}

void InstationaryFlowApp::solve_nlp ( )
{
    NonlinearSolverState state = nls_->Solve ( rhs_, &sol_ );
    if ( rank_ == master_rank_ )
        std::cout << "Nonlinear solver ended with state " << state
            << " and residual norm " << nls_->GetResidual ( )
        << " after " << nls_->iter ( ) << " iterations\n";
}

struct DirichletBC2d
{
    // Material numbers:
    // 10 -> left edge (inflow)
    // 11 -> bottom edge (u = 0)
    // 12 -> right edge (outflow)
    // 13 -> top edge (u = 0)
    // 14 -> cylinder (u = 0)

    DirichletBC2d ( int var, double H, double Um )
    : var_ ( var ), H_ ( H ), Um_ ( Um )
    {
        assert ( var_ == 0 || var_ == 1 );

    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );

        const bool inflow = ( material_num == 10 );
        const bool outflow = ( material_num == 12 );

        if ( !outflow )
        {
            values.resize ( coords_on_face.size ( ) );

            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {

                // evaluate dirichlet function at each point
                const Coord& pt = coords_on_face[i];
                if ( inflow )
                {
                    if ( var_ == 0 )
                    { // x-component
                        values[i] = 4. * Um_ * pt[1] * ( H_ - pt[1] ) / ( H_ * H_ );
                    }
                    else if ( var_ == 1 )
                    { // y-component
                        values[i] = 0.;
                    }
                    else
                    {
                        assert ( false );
                    }
                }
                else
                {
                    // not inflow: u = 0
                    values[i] = 0.;
                }
            }

            // outflow: do-nothing
        }
        return values;
    }

    const int var_;
    const double H_;
    const double Um_;
};

struct DirichletBC3d
{
    // Material numbers:
    // 10 -> left edge (inflow)
    // 11 -> bottom edge (u = 0)
    // 12 -> right edge (outflow)
    // 13 -> top edge (u = 0)
    // 14 -> cylinder (u = 0)

    DirichletBC3d ( int var, double H, double Um )
    : var_ ( var ), H_ ( H ), Um_ ( Um )
    {
        assert ( var_ == 0 || var_ == 1 || var_ == 2 );

    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );
        assert ( material_num >= 10 && material_num <= 14 );

        const bool inflow = ( material_num == 10 );
        const bool outflow = ( material_num == 12 );

        if ( !outflow )
        {
            values.resize ( coords_on_face.size ( ) );

            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {

                // evaluate dirichlet function at each point
                const Coord& pt = coords_on_face[i];
                if ( inflow )
                {
                    if ( var_ == 0 )
                    { // x-component
                        const double y = pt[1];
                        const double z = pt[2];
                        values[i] = 16. * Um_ * y * z * ( H_ - y ) * ( H_ - z ) / ( std::pow ( H_, 4. ) );
                    }
                    else if ( var_ == 1 || var_ == 2 )
                    { // y- and z-components
                        values[i] = 0.;
                    }
                    else
                    {
                        assert ( false );
                    }
                }
                else
                {
                    // not inflow: u = 0
                    values[i] = 0.;
                }
            }

            // outflow: do-nothing
        }
        return values;
    }

    const int var_;
    const double H_;
    const double Um_;
};

void InstationaryFlowApp::prepare_bc ( )
{
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    if ( DIM == 2 )
    {
        DirichletBC2d bc[2] = { DirichletBC2d ( 0, H_, Um_ ),
                               DirichletBC2d ( 1, H_, Um_ ) };

        for ( int var = 0; var < DIM; ++var )
        {
            compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                                dirichlet_dofs_, dirichlet_values_ );
        }
    }
    else
    {
        assert ( DIM == 3 );
        DirichletBC3d bc[3] = { DirichletBC3d ( 0, H_, Um_ ),
                               DirichletBC3d ( 1, H_, Um_ ),
                               DirichletBC3d ( 2, H_, Um_ ) };
        for ( int var = 0; var < DIM; ++var )
        {
            compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                                dirichlet_dofs_, dirichlet_values_ );
        }
    }
    //std::cerr << "Num dirichlet dofs = " << dirichlet_dofs_.size() << "\n";
}

void InstationaryFlowApp::visualize_solution ( int step )
{
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, master_rank_ );

    sol_.UpdateCouplings ( );

    std::stringstream input;

    if ( step < 10 )
        input << "instationary_flow_values000" << step;
    else if ( step < 100 )
        input << "instationary_flow_values00" << step;
    else if ( step < 1000 )
        input << "instationary_flow_values0" << step;
    else
        input << "instationary_flow_values" << step;
    std::vector<double> material_number ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
          it != mesh_->end ( mesh_->tdim ( ) );
          ++it )
    {
        material_number.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
    }

    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 0 ), "u" );
    if ( DIM >= 2 )
    {
        visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 1 ), "v" );
    }
    if ( DIM == 3 )
    {
        visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 2 ), "w" );
    }
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, DIM ), "p" );

    visu.visualize_cell_data ( material_number, "Material Id" );

    visu.write ( "out/" + input.str ( ) );
}

void InstationaryFlowApp::compute_alphas ( int sub_step, double delta_t, std::vector<double>* alphas )
{
    assert ( alphas->size ( ) == 5 );

    if ( method_ == "FractionalStepThetaScheme" )
    {
        double theta = 1. - ( sqrt ( 2. ) / 2. );
        double theta_p = 1. - ( 2. * theta );
        double beta = theta_p / ( 1. - theta );
        double gamma = 1. - beta;

        if ( sub_step == 0 )
        {
            alphas->at ( 0 ) = beta * theta*delta_t;
            alphas->at ( 1 ) = theta*delta_t;
            alphas->at ( 2 ) = gamma * theta*delta_t;
            alphas->at ( 3 ) = theta*delta_t;
            alphas->at ( 4 ) = 0.;
        }
        else if ( sub_step == 1 )
        {
            alphas->at ( 0 ) = gamma * theta_p*delta_t;
            alphas->at ( 1 ) = theta_p*delta_t;
            alphas->at ( 2 ) = beta * theta_p*delta_t;
            alphas->at ( 3 ) = 0.;
            alphas->at ( 4 ) = theta_p*delta_t;
        }
        else if ( sub_step == 2 )
        {
            alphas->at ( 0 ) = beta * theta*delta_t;
            alphas->at ( 1 ) = theta*delta_t;
            alphas->at ( 2 ) = gamma * theta*delta_t;
            alphas->at ( 3 ) = theta*delta_t;
            alphas->at ( 4 ) = 0.;
        }
        else
            assert ( 0 );
    }
    else if ( method_ == "ImplicitEuler" )
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
    else
        assert ( 0 );
}

void convert_to_std_vector ( const LAD::VectorType& in, std::vector<double>& out )
{
    const int size = in.size_global ( );
    std::vector<int> ind ( size );
    for ( int i = 0; i < size; ++i )
    {
        ind[i] = i;
    }
    out.resize ( size );
    in.GetValues ( vec2ptr ( ind ), size, vec2ptr ( out ) );
}

void InstationaryFlowAssembler::set_newton_solution ( const LAD::VectorType& newton_sol )
{
    newton_sol_ = &newton_sol;
}

void InstationaryFlowAssembler::set_prev_solution ( const LAD::VectorType& prev_sol )
{
    prev_sol_ = &prev_sol;
}

void InstationaryFlowAssembler::set_timestep_parameters ( const std::vector<double>& alphas )
{
    assert ( alphas.size ( ) == 5 );
    alpha1_ = alphas[0];
    alpha2_ = alphas[1];
    alpha3_ = alphas[2];
    alpha4_ = alphas[3];
    alpha5_ = alphas[4];
}

void InstationaryFlowAssembler::initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
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
    pressure_k_.clear ( );
    evaluate_fe_function ( *newton_sol_, DIM, pressure_k_ );
}

void InstationaryFlowAssembler::initialize_for_facet ( const Element<double>& element, const Quadrature<double>& quadrature, int facet_number )
{
    AssemblyAssistant<DIM, double>::initialize_for_facet ( element, quadrature, facet_number );
}

void InstationaryFlowAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
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

        // assemble b(p, v) = - \alpha_2 * \int{p div{v}}
        const int p_var = DIM;
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                for ( int j = 0; j < num_dofs ( p_var ); ++j )
                {
                    lm ( dof_index ( i, v_var ), dof_index ( j, p_var ) ) +=
                            -wq * alpha2_ * ( inv_rho_ * phi ( j, q, p_var ) *
                            grad_phi ( i, q, v_var )[v_var] ) * dJ;
                }
            }
        }

        // assemble bT(u, q) = - \int{q div(u)}
        const int q_var = DIM;
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < num_dofs ( q_var ); ++i )
            {
                for ( int j = 0; j < num_dofs ( u_var ); ++j )
                {
                    lm ( dof_index ( i, q_var ), dof_index ( j, u_var ) ) +=
                            wq * ( inv_rho_ * phi ( i, q, q_var ) *
                            grad_phi ( j, q, u_var )[u_var] ) * dJ;
                }
            }
        }
    }
}

void InstationaryFlowAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
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

        // l3(v) = 1/rho*\alpha_2*\int(p_k*div(v))
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] +=
                        wq * alpha2_ * ( inv_rho_ * pressure_k_[q] * grad_phi ( i, q, v_var )[v_var] ) * dJ;
            }
        }

        // l4(q) = -\int(q * div(u_k))
        const int q_var = DIM;
        double div_u_k = 0.;
        for ( int d = 0; d < DIM; ++d )
        {
            div_u_k += grad_prev_ns_vel_[d][q][d];
        }

        for ( int i = 0; i < num_dofs ( q_var ); ++i )
        {
            lv[dof_index ( i, q_var )] +=
                    -wq * ( inv_rho_ * div_u_k * phi ( i, q, q_var ) ) * dJ;
        }
    }
}

void InstationaryFlowApp::compute_residual ( const LAD::VectorType& u, LAD::VectorType* F )
{
    local_asm_->set_newton_solution ( u );
    StandardGlobalAssembler<double> global_asm;
    global_asm.assemble_vector ( space_, *local_asm_, *F );
}

void InstationaryFlowApp::compute_matrix ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
    local_asm_->set_newton_solution ( u );
    StandardGlobalAssembler<double> global_asm;
    global_asm.assemble_matrix ( space_, *local_asm_, *DF );

}
