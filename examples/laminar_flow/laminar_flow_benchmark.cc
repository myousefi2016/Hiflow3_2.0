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

#include "laminar_flow_benchmark.h"
#include <iostream>

static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;

/// Program main loop: setup MPI, read parameters and start the application

int main ( int argc, char** argv )
{

    MPI_Init ( &argc, &argv );

    // Read parameters
    std::string fName ( "laminar_flow_benchmark.xml" );
    PropertyTree config ( fName, MASTER_RANK, MPI_COMM_WORLD );

    // Start Logging
    std::ofstream info_file ( config["Output"]["InfoFilename"].get<std::string>( ).c_str ( ) );
    //LogKeeper::get_log("info").set_target(&info_file);
    LogKeeper::get_log ( "info" ).set_target ( &std::cout );
    LogKeeper::get_log ( "error" ).set_target ( &std::cout );
    int size = -1;
    MPI_Comm_size ( MPI_COMM_WORLD, &size );
    LOG_INFO ( "MPI Processes", size );

    // Run application
    LaminarFlowApp application ( config );
    application.run ( );

    // flush log here to avoid problems
    LogKeeper::get_log ( "info" ).flush ( );

    MPI_Finalize ( );
    return 0;
}

LaminarFlowApp::LaminarFlowApp ( PropertyTree &config )
: comm_ ( MPI_COMM_WORLD ),
num_partitions_ ( -1 ),
rank_ ( -1 ),
master_mesh_ ( 0 ),
damping_ ( false ),
forcing_ ( false ),
#ifdef WITH_MUMPS
linear_solver_ ( new MumpsSolver<LAD, MumpsStructureD>( ) ),
#endif
refinement_level_ ( 1 )
{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    assert ( DIM == 2 || DIM == 3 );
    LOG_INFO ( "Problem dimension", DIM );
    refinement_level_ = config["Mesh"]["RefinementLevel"].get<int>( );
    LOG_INFO ( "Refinement level", refinement_level_ );

    if ( rank_ == MASTER_RANK )
    {
        // read mesh (sequential, dimension DIM
        std::string mesh_filename = std::string ( DATADIR ) + config["Mesh"]["Filename"].get<std::string>( );
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIM, DIM, 0 );

        for ( int r = 0; r < refinement_level_; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
        }
    }

    MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, MPI_COMM_WORLD );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, MPI_COMM_WORLD, shared_verts );

    std::ostringstream rank_str;
    rank_str << rank_;
    PVtkWriter writer ( MPI_COMM_WORLD );
    std::string output_file = config["Output"]["MeshFilename"].get<std::string>( );
    writer.write ( output_file.c_str ( ), *mesh_ );

    // setup linear algebra platform
    la_sys_.Platform = APP_PLATFORM;
    init_platform ( la_sys_ );

    // setup NavierStokes weak form
    config["Application"]["Density"].read<double>( rho_ ); // density
    config["Application"]["Viscosity"].read<double>( nu_ ); // kinematic viscosity
    LOG_INFO ( "Density", "rho = " << rho_ );
    LOG_INFO ( "Kinematic viscosity", "nu = " << nu_ );

    config["Application"]["ChannelHeight"].read<double>( H_ );
    config["Application"]["MeanInflowSpeed"].read<double>( Um_ );
    LOG_INFO ( "Mean inflow speed [m/s]", Um_ );

    time_measurement_res_ = 0.0;
    res_count_ = 0;
    time_measurement_mat_ = 0.0;
    mat_count_ = 0;

    // setup quadrature rule to use
    QuadratureFactory<double> QuadFact;
    quadratures_.resize ( 1 );
    quadratures_[0] = QuadFact.Get ( "Default" )->params ( config["Quadrature"] );

    // setup solvers
#ifndef WITH_MUMPS
    LinearSolverFactory<LAD> LinSolFact;
    linear_solver_ = LinSolFact.Get ( config["LinearSolver"]["Name"].get<std::string>( ) )->params ( config["LinearSolver"] );
#endif

    NonlinearSolverFactory<LAD> NLSFact;
    nls_ = NLSFact.Get ( config["NonlinearSolver"]["Name"].get<std::string>( ) )->params ( &res_, &matrix_, config["NonlinearSolver"] );
    // we use our own initial solution -- this needs to be indicated
    // to the Newton-solver
    if ( config["NonlinearSolver"]["Name"].get<std::string>( ) == "Newton" )
    {
        ( ( Newton<LAD>* )nls_ )->InitParameter ( Newton<LAD>::NewtonInitialSolutionOwn );
        // Armijo Damping
        if ( config["NonlinearSolver"]["Damping"].get<std::string>( ) == "Armijo" )
        {
            damping_strategy_ = new ArmijoDamping<LAD>;
            damping_ = true;
            ( ( Newton<LAD>* )nls_ )->SetDampingStrategy ( *damping_strategy_ );
        }
        if ( config["NonlinearSolver"]["Forcing"].get<std::string>( ) == "EisenstatWalker1" )
        {
            forcing_strategy_ = new EWForcing<LAD>( 1 );
            forcing_ = true;
            ( ( Newton<LAD>* )nls_ )->SetForcingStrategy ( *forcing_strategy_ );
        }
        if ( config["NonlinearSolver"]["Forcing"].get<std::string>( ) == "EisenstatWalker2" )
        {
            forcing_strategy_ = new EWForcing<LAD>( 2 );
            forcing_ = true;
            ( ( Newton<LAD>* )nls_ )->SetForcingStrategy ( *forcing_strategy_ );
        }
    }
}

/// Destructor: clean up after application

LaminarFlowApp::~LaminarFlowApp ( )
{
#ifndef WITH_MUMPS
    if ( linear_solver_ )
    {
        delete linear_solver_;
    }
    if ( damping_ )
    {
        delete damping_strategy_;
    }
    if ( forcing_ )
    {
        delete forcing_strategy_;
    }
#endif
    stop_platform ( la_sys_ );
}

/// Main loop of application

void LaminarFlowApp::run ( )
{
    prepare_space ( );
    prepare_bc ( );
    prepare_lin_alg_structures ( ); // needs only be done once, since space does not change
    prepare_linear_solver ( );
    prepare_nls ( );

    solve_nls ( );

    visualize_solution ( );

    std::cout << "Mean matrix assemble time: " << time_measurement_mat_ / mat_count_ << std::endl;
    std::cout << "Mean residual assemble time: " << time_measurement_res_ / res_count_ << std::endl;

#ifdef WITH_MUMPS
    linear_solver_->Clear ( );
#endif
}

void LaminarFlowApp::EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F )
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

void LaminarFlowApp::EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
    // assemble the matrix for the linearized system
    compute_matrix ( u, DF );

    // correct BC -- set Dirichlet rows to identity
    if ( !dirichlet_dofs_.empty ( ) )
    {
        DF->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
    }
#ifdef WITH_MUMPS
    linear_solver_->Clear ( );
    linear_solver_->InitParameter ( comm_ );
    linear_solver_->SetupOperator ( matrix_ );
#endif

    // DEBUG
    //    DF->diagonal().WriteFile("yalaplace_matrix.mat");
}

void LaminarFlowApp::prepare_space ( )
{
    // setup space Q2-Q1 elements
    std::vector<int> degrees ( DIM + 1, 2 );
    degrees[DIM] = 1;
    space_.Init ( degrees, *mesh_ );
}

void LaminarFlowApp::prepare_linear_solver ( )
{
#ifdef WITH_MUMPS
    // linear_solver_->InitParameter(MPI_COMM_WORLD);
    linear_solver_->InitParameter ( comm_ );
#else
    // set the matrix to be used as the operator
    linear_solver_->SetupOperator ( matrix_ );
#endif
}

void LaminarFlowApp::prepare_lin_alg_structures ( )
{
    // Initialize linear algebra structures
    couplings_.Clear ( );
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
                   LaminarFlowApp::APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    matrix_.InitStructure ( vec2ptr ( diagonal_rows ),
                            vec2ptr ( diagonal_cols ),
                            diagonal_rows.size ( ),
                            vec2ptr ( off_diagonal_rows ),
                            vec2ptr ( off_diagonal_cols ),
                            off_diagonal_rows.size ( ) );

    sol_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    sol_.InitStructure ( );
    sol_.Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }

    rhs_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    rhs_.InitStructure ( );
    rhs_.Zeros ( );
}

void LaminarFlowApp::prepare_nls ( )
{
    nls_->SetOperator ( *this );
    nls_->SetLinearSolver ( *linear_solver_ );
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

        const Coord& pt0 = coords_on_face[0];
        const bool outflow = std::abs ( pt0[0] - 2.2 ) < 1.e-9;
        const bool inflow = std::abs ( pt0[0] - 0.0 ) < 1.e-9;

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
        const double COORD_TOL = 1.e-9;

        const Coord& pt0 = coords_on_face[0];
        const bool outflow = std::abs ( pt0[0] - 2.5 ) < COORD_TOL;
        const bool inflow = std::abs ( pt0[0] - 0.0 ) < COORD_TOL;

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

        }

        return values;
    }

    const int var_;
    const double H_;
    const double Um_;
};

void LaminarFlowApp::prepare_bc ( )
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

void LaminarFlowApp::solve_nls ( )
{
    //    std::vector<double> out;
    //    convert_to_std_vector(sol_, out);
    //    std::cout << "Initial solution = " << string_from_range(out.begin(), out.end()) << "\n";

    NonlinearSolverState state = nls_->Solve ( rhs_, &sol_ );
    std::cout << "Nonlinear solver ended with state " << state
            << " and residual norm " << nls_->GetResidual ( )
            << " after " << nls_->iter ( ) << " iterations\n";

    LOG_INFO ( "Nonlinear solver residual", nls_->GetResidual ( ) );
    LOG_INFO ( "Nonlinear solver steps", nls_->iter ( ) );

    //    out.clear();
    //    convert_to_std_vector(sol_, out);
    //    std::cout << "Final solution = " << string_from_range(out.begin(), out.end()) << "\n";
}

void LaminarFlowApp::visualize_solution ( )
{
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    sol_.UpdateCouplings ( );

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

    visu.write ( "laminar_flow_values" );
}

// LaminarFlowAssembler ////////////////

class LaminarFlowAssembler : private AssemblyAssistant<LaminarFlowApp::DIM, double>
{
  public:
    static const int DIM = LaminarFlowApp::DIM;

    LaminarFlowAssembler ( const LAD::VectorType& solution, double nu, double rho )
    : solution_ ( solution ), nu_ ( nu ), inv_rho_ ( 1. / rho )
    {
    }

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
    {
        AssemblyAssistant<LaminarFlowApp::DIM, double>::initialize_for_element ( element, quadrature );

        // recompute previous solution values
        for ( int v = 0; v < DIM; ++v )
        {
            prev_vel_[v].clear ( );
            grad_prev_vel_[v].clear ( );
            evaluate_fe_function ( solution_, v, prev_vel_[v] );
            evaluate_fe_function_gradients ( solution_, v, grad_prev_vel_[v] );
        }
        pressure_k_.clear ( );
        evaluate_fe_function ( solution_, DIM, pressure_k_ );
    }

    void initialize_for_facet ( const Element<double>& element, const Quadrature<double>& quadrature, int facet_number )
    {
        AssemblyAssistant<DIM, double>::initialize_for_facet ( element, quadrature, facet_number );
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
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
                vel_k[var] = prev_vel_[var][q];
            }

            // assemble a1(u,v) = \int {\grad(u) : \grad(v)}
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * ( nu_ * dot ( grad_phi ( j, q, u_var ), grad_phi ( i, q, u_var ) ) ) * dJ;
                    }
                }
            }

            // assemble a2(u,v) = \int { (vel_k*\grad{u})*v }
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * ( dot ( vel_k, grad_phi ( j, q, u_var ) ) * phi ( i, q, u_var ) ) * dJ;
                    }
                }
            }

            // assemble a3(u,v) = \int { (u\grad{u_k}*v }
            for ( int test_var = 0; test_var < DIM; ++test_var )
            {
                for ( int trial_var = 0; trial_var < DIM; ++trial_var )
                {
                    for ( int i = 0; i < num_dofs ( test_var ); ++i )
                    {
                        for ( int j = 0; j < num_dofs ( trial_var ); ++j )
                        {
                            lm ( dof_index ( i, test_var ), dof_index ( j, trial_var ) ) +=
                                    wq * ( grad_prev_vel_[test_var][q][trial_var] *
                                    phi ( j, q, trial_var ) *
                                    phi ( i, q, test_var ) ) * dJ;
                        }
                    }
                }
            }

            // assemble b(p, v) = - \int{p div{v}}
            const int p_var = DIM;
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( p_var ); ++j )
                    {
                        lm ( dof_index ( i, v_var ), dof_index ( j, p_var ) ) +=
                                -wq * ( inv_rho_ * phi ( j, q, p_var ) *
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
#if 0
        if ( element.get_cell_index ( ) == 0 )
        {
            std::cerr << "local matrix =\n";
            for ( int i = 0; i < total_dofs; ++i )
            {
                for ( int j = 0; j < total_dofs; ++j )
                {
                    std::cerr << lm ( i, j ) << " ";
                }
                std::cerr << "\n";
            }
        }
#endif
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
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

            // get previous solution in vector form
            Vec<DIM, double> vel_k;
            for ( int var = 0; var < DIM; ++var )
            {
                vel_k[var] = prev_vel_[var][q];
            }

            // l1(v) = -\nu * \int( \grad{u_k} : \grad{v} )
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            -wq * ( nu_ * dot ( grad_phi ( i, q, v_var ), grad_prev_vel_[v_var][q] ) ) * dJ;
                }
            }

            // l2(v) = -\int(u_k*\grad{u_k}*v)
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            -wq * ( dot ( grad_prev_vel_[v_var][q], vel_k ) * phi ( i, q, v_var ) ) * dJ;
                }
            }

            // l3(v) = 1/rho*\int(p_k*div(v))
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            wq * ( inv_rho_ * pressure_k_[q] * grad_phi ( i, q, v_var )[v_var] ) * dJ;
                }
            }

            // l4(q) = -\int(q * div(u_k))
            const int q_var = DIM;
            double div_u_k = 0.;
            for ( int d = 0; d < DIM; ++d )
            {
                div_u_k += grad_prev_vel_[d][q][d];
            }

            for ( int i = 0; i < num_dofs ( q_var ); ++i )
            {
                lv[dof_index ( i, q_var )] +=
                        -wq * ( inv_rho_ * div_u_k * phi ( i, q, q_var ) ) * dJ;
            }
        }
#if 0
        if ( element.get_cell_index ( ) == 1 )
        {
            std::cerr << "local vector(1) = " << string_from_range ( lv.begin ( ), lv.end ( ) ) << "\n";
        }
#endif
    }

    void assemble_local_bc_vector ( const Element<double>& element, LocalVector& lv, int facet_number ) const
    {
    }

  private:
    const LAD::VectorType& solution_;
    double nu_, inv_rho_;
    FunctionValues<double> prev_vel_[DIM];
    FunctionValues<double> pressure_k_;
    FunctionValues< Vec<DIM, double> > grad_prev_vel_[DIM];
};

void LaminarFlowApp::compute_residual ( const LAD::VectorType& u, LAD::VectorType* F )
{
    Timer res_timer;

    LaminarFlowAssembler local_asm ( u, nu_, rho_ );
    StandardGlobalAssembler<double> global_asm;
    global_asm.assemble_vector ( space_, local_asm, *F );

    res_timer.stop ( );
    // std::cerr << "Residual computation took " << res_timer << "\n";
    time_measurement_res_ += res_timer.get_duration ( );
    ++res_count_;

#if 0
    std::vector<double> dbg;
    convert_to_std_vector ( *F, dbg );
    std::cerr << "res = " << string_from_range ( dbg.begin ( ), dbg.end ( ) ) << "\n";
#endif
}

void LaminarFlowApp::compute_matrix ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
    Timer asm_timer;
    LaminarFlowAssembler local_asm ( u, nu_, rho_ );
    StandardGlobalAssembler<double> global_asm;
    global_asm.assemble_matrix ( space_, local_asm, *DF );

    asm_timer.stop ( );
    // std::cerr << "Matrix assembly took " << asm_timer << "\n";
    time_measurement_mat_ += asm_timer.get_duration ( );
    ++mat_count_;
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

//// DIRICHLET BC based on material id -- to be used when material id works correctly ////

#if 0

struct DirichletBC
{
    // Material numbers:
    // 10 -> left edge (inflow)
    // 11 -> bottom edge (u = 0)
    // 12 -> right edge (outflow)
    // 13 -> top edge (u = 0)
    // 14 -> cylinder (u = 0)

    DirichletBC ( int var, double H, double Um )
    : var_ ( var ), H_ ( H ), Um_ ( Um )
    {
        assert ( var_ == 0 || var_ == 1 );

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
            for ( int i = 0; i < coords_on_face.size ( ); ++i )
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
#endif
