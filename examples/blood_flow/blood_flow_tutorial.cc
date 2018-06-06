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

/// \author Jonas Kratzke, Katrin Mang

#include "blood_flow_tutorial.h"

#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif

// Mesh directory, main MPI process and parameter filename
static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;
static const char* PARAM_FILENAME = "blood_flow_tutorial.xml";

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    // Set default parameter file
    std::string param_filename ( PARAM_FILENAME );
    if ( argc > 1 )
        param_filename = argv[1];
    try
    {
        BloodFlowTutorial app ( param_filename );
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

BloodFlowTutorial::BloodFlowTutorial ( const std::string& param_filename )
: comm_ ( MPI_COMM_WORLD ),
params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
matrix_ ( 0 ), sol_ ( 0 ), prev_sol_ ( 0 ), res_ ( 0 ), rhs_ ( 0 ), results_writer_ ( "results.csv" )
{
}

// Destructor

BloodFlowTutorial::~BloodFlowTutorial ( )
{
    delete sol_;
    delete prev_sol_;
    delete res_;
    delete rhs_;
    delete matrix_;
    delete nls_;
    delete linsolver_;
    if ( nls_forcing_ ) delete nls_forcing_;
    if ( nls_damping_ ) delete nls_damping_;
}

void BloodFlowTutorial::run ( )
{
    simul_name_ = params_["OutputPrefix"].get<std::string>( );

    // Log Output on console
    LogKeeper::get_log ( "info" ).set_target ( &std::cout );
    // std::ofstream info_log((simul_name_ + "_info_log").c_str());
    // LogKeeper::get_log("info").set_target(&info_log);
    LogKeeper::get_log ( "debug" ).set_target ( &std::cout );
    // std::ofstream debug_log((simul_name_ + "_debug_log").c_str());
    // LogKeeper::get_log("debug").set_target(&debug_log);
    LogKeeper::get_log ( "error" ).set_target ( &std::cout );
    // std::ofstream error_log((simul_name_ + "_error_log").c_str());
    // LogKeeper::get_log("error").set_target(&error_log);

    // Output parameters for debugging
    LOG_INFO ( "parameters", params_ );

    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    read_mesh ( );

    prepare_parameters ( );

    prepare_space ( );

    prepare_solver ( );

    run_time_loop ( );
}

void BloodFlowTutorial::EvalFunc ( const CVector& newton_sol, CVector* res )
{

    // Compute the residual vector
    InstationaryFlowAssembler local_asm ( nu_, rho_, newton_sol, *prev_sol_, alpha1_, dt_, alpha2_, stab_fac_ );
    global_asm_.assemble_vector ( space_, local_asm, *res );
    WindkesselBoundaryAssembler windkessel_local_asm ( windkessel_materials_, rho_, dt_, eval_windkessel_flowrates_, eval_windkessel_pressures_, decay_, resistance_, compliance_ );
    global_asm_.should_reset_assembly_target ( false );
    global_asm_.assemble_vector_boundary ( space_, windkessel_local_asm, *res );
    global_asm_.should_reset_assembly_target ( true );

    // Correct BC -- set Dirichlet dofs to 0
    if ( !dirichlet_dofs_.empty ( ) )
    {
        std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
        res->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                         vec2ptr ( zeros ) );
    }
}

void BloodFlowTutorial::EvalGrad ( const CVector& newton_sol, CMatrix* jac )
{

    // Assemble the matrix for the linearized system
    InstationaryFlowAssembler local_asm ( nu_, rho_, newton_sol, *prev_sol_, alpha1_, dt_, alpha2_, stab_fac_ );
    global_asm_.assemble_matrix ( space_, local_asm, *jac );

    // Correct BC -- set Dirichlet rows to identity
    if ( !dirichlet_dofs_.empty ( ) )
    {
        jac->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
    }

    // Set preconditioning operator if enabled
#ifdef WITH_ILUPP
    {
        if ( use_ilupp_ )
        {
            ilupp_.SetupOperator ( *jac );
        }
    }
#endif
}

void BloodFlowTutorial::read_mesh ( )
{
    LOG_INFO ( "mesh", "reading in mesh" );

    MeshPtr master_mesh;

    if ( rank_ == MASTER_RANK )
    {
        const std::string mesh_name =
                params_["Mesh"]["Filename"].get<std::string>( );
        std::string mesh_filename = std::string ( DATADIR ) + mesh_name;

        master_mesh = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0 );

        const int ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );

        for ( int r = 0; r < ref_lvl; ++r )
        {
            master_mesh = master_mesh->refine ( );
        }
        LOG_INFO ( "mesh", "Refinement level = " << ref_lvl );
    }

    MeshPtr local_mesh = partition_and_distribute ( master_mesh, MASTER_RANK, comm_ );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );

    std::ostringstream rank_str;
    rank_str << rank_;

    PVtkWriter writer ( comm_ );
    std::string output_file = std::string ( "mesh_local.pvtu" );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

void BloodFlowTutorial::prepare_parameters ( )
{
    LOG_INFO ( "Prepare", "Parameters" );

    // Prepare timestep
    ts_ = 0;
    dt_ = params_["Instationary"]["Timestep"].get<Scalar>( );
    smooth_start_up_time_ = params_["Instationary"]["SmoothStartTime"].get<Scalar>( );
    period_ = params_["Instationary"]["Period"].get<Scalar>( );
    visu_interval_ = params_["Instationary"]["VisualizationInterval"].get<int>( );
    LOG_INFO ( "Prepare", "dt: " << dt_ << ", period: " << period_ << ", smooth start: " << smooth_start_up_time_ );

    // Set the alpha coefficients for the Crank-Nicolson method
    Scalar theta = params_["Instationary"]["Theta"].get<Scalar>( );
    alpha1_ = theta * dt_;
    alpha2_ = ( 1. - theta ) * dt_;

    // Prepare the problem parameters
    rho_ = params_["FlowModel"]["Density"].get<Scalar>( );
    nu_ = params_["FlowModel"]["Viscosity"].get<Scalar>( );
    stab_fac_ = params_["FlowModel"]["StabilizationFactor"].get<Scalar>( );
    LOG_INFO ( "Prepare", "Rho: " << rho_ << ", viscosity: " << nu_ << ", stabilization: " << stab_fac_ );

    // Prepare BC parameters for Dirichlet boundary
    std::stringstream ss;
    wall_bdy_ = params_["Boundary"]["WallMaterial"].get<int>( );

    // Resize all parameters for Dirichlet/Poiseuille boundaries
    const int num_poiseuille_bcs = params_["Boundary"]["NumberOfPoiseuilleBoundaries"].get<int>( );
    LOG_INFO ( "Prepare", "Num poiseuille bcs: " << num_poiseuille_bcs );
    poiseuille_materials_.resize ( num_poiseuille_bcs );
    eval_poiseuille_flowrates_.resize ( poiseuille_materials_.size ( ) );
    eval_poiseuille_pressures_.resize ( poiseuille_materials_.size ( ) );
    radius_.resize ( num_poiseuille_bcs );
    exponent_.resize ( num_poiseuille_bcs );
    normal_.resize ( num_poiseuille_bcs );
    timestamps_.resize ( num_poiseuille_bcs );
    flow_rate_.resize ( num_poiseuille_bcs );
    center_.resize ( num_poiseuille_bcs );

    std::vector<std::string> results_names;
    for ( int i = 0; i < num_poiseuille_bcs; ++i )
    {
        std::stringstream ss;
        ss << "Poiseuille_" << i;
        // Read in the Poiseuille boundary parameters
        poiseuille_materials_[i] = params_["Boundary"][ss.str ( )]["Material"].get<int>( );
        radius_[i] = params_["Boundary"][ss.str ( )]["Radius"].get<Scalar>( );
        exponent_[i] = params_["Boundary"][ss.str ( )]["Exponent"].get<int>( );
        params_["Boundary"][ss.str ( )]["Normal"].read<Scalar>( normal_[i] );
        params_["Boundary"][ss.str ( )]["TimeStamps"].read<Scalar>( timestamps_[i] );
        params_["Boundary"][ss.str ( )]["FlowRates"].read<Scalar>( flow_rate_[i] );
        params_["Boundary"][ss.str ( )]["Center"].read<Scalar>( center_[i] );

        // Evaluate spline for each Poiseuille Boundary -> vector with three interpolation values
        calculate_spline ( i );
        if ( rank_ == MASTER_RANK )
        {
            std::stringstream name_f;
            name_f << "flow_bdy_" << poiseuille_materials_[i];
            results_names.push_back ( name_f.str ( ) );
            std::stringstream name_p;
            name_p << "press_bdy_" << poiseuille_materials_[i];
            results_names.push_back ( name_p.str ( ) );
        }
    }

    // Prepare BC parameters for Neumann/Windkessel boundary
    const int num_windkessel_bcs = params_["Boundary"]["NumberOfWindkesselBoundaries"].get<int>( );
    windkessel_materials_.resize ( num_windkessel_bcs );
    eval_windkessel_flowrates_.resize ( windkessel_materials_.size ( ), 0 );
    eval_windkessel_pressures_.resize ( windkessel_materials_.size ( ), 0 );
    resistance_.resize ( num_windkessel_bcs );
    compliance_.resize ( num_windkessel_bcs );
    decay_.resize ( num_windkessel_bcs );

    for ( int i = 0; i < num_windkessel_bcs; ++i )
    {
        std::stringstream ss;
        ss << "Windkessel_" << i;
        windkessel_materials_[i] = params_["Boundary"][ss.str ( )]["Material"].get<int>( );
        resistance_[i] = params_["Boundary"][ss.str ( )]["Resistance"].get<Scalar>( );
        compliance_[i] = params_["Boundary"][ss.str ( )]["Compliance"].get<Scalar>( );
        decay_[i] = params_["Boundary"][ss.str ( )]["Decay"].get<Scalar>( );

        if ( rank_ == MASTER_RANK )
        {
            std::stringstream name_f;
            name_f << "flow_bdy_" << windkessel_materials_[i];
            results_names.push_back ( name_f.str ( ) );
            std::stringstream name_p;
            name_p << "press_bdy_" << windkessel_materials_[i];
            results_names.push_back ( name_p.str ( ) );
        }
    }
    results_writer_.Init ( results_names );
}

void BloodFlowTutorial::prepare_space ( )
{
    LOG_INFO ( "Prepare", "FEM space" );

    // Prepare space
    std::vector< int > degrees ( DIMENSION + 1 );
    const int u_deg = params_["FiniteElements"]["VelocityDegree"].get<int>( );
    const int p_deg = params_["FiniteElements"]["PressureDegree"].get<int>( );
    for ( int c = 0; c < DIMENSION; ++c )
    {
        degrees.at ( c ) = u_deg;
    }
    degrees.at ( DIMENSION ) = p_deg;

    space_.Init ( degrees, *mesh_ );

    // Prepare linear algebra structures
    couplings_.Clear ( );
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute matrix graph

    std::vector < std::vector<bool> > coupling_vars;
    coupling_vars.resize ( DIMENSION + 1 );
    for ( int i = 0; i < DIMENSION; ++i )
    {
        for ( int j = 0; j < DIMENSION + 1; ++j )
        {
            coupling_vars[i].push_back ( true );
        }
    }
    for ( int i = 0; i < DIMENSION; ++i )
    {
        coupling_vars[DIMENSION].push_back ( true );
    }
    coupling_vars[DIMENSION].push_back ( false );

    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity, &coupling_vars );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    CoupledMatrixFactory<Scalar> CoupMaFact;
    CoupledVectorFactory<Scalar> CoupVecFact;

    // Prepare linear equation system: matrix, right hand side, solution, previous solution
    matrix_ = CoupMaFact.Get ( params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->params ( params_["LinearAlgebra"] );
    matrix_->Init ( comm_, couplings_ );
    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );
    matrix_->Zeros ( );

    // Prepare solution vector
    sol_ = CoupVecFact.Get ( params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->params ( params_["LinearAlgebra"] );
    sol_->Init ( comm_, couplings_ );
    sol_->InitStructure ( );
    sol_->Zeros ( );

    // Prepare solution vector of last step
    prev_sol_ = CoupVecFact.Get ( params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->params ( params_["LinearAlgebra"] );
    prev_sol_->Init ( comm_, couplings_ );
    prev_sol_->InitStructure ( );
    prev_sol_->Zeros ( );

    // Prepare result vector
    res_ = CoupVecFact.Get ( params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->params ( params_["LinearAlgebra"] );
    res_->Init ( comm_, couplings_ );
    res_->InitStructure ( );
    res_->Zeros ( );

    // Prepare Right hand side vector
    rhs_ = CoupVecFact.Get ( params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->params ( params_["LinearAlgebra"] );
    rhs_->Init ( comm_, couplings_ );
    rhs_->InitStructure ( );
    rhs_->Zeros ( );
}

void BloodFlowTutorial::prepare_solver ( )
{
    LOG_INFO ( "Prepare", "Solver" );

    // Setup linear solver
    LinearSolverFactory<LAD> LinSolFact;
    linsolver_ = LinSolFact.Get (
                                  params_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( params_["LinearSolver"] );
    linsolver_->SetupOperator ( *matrix_ );
    int basis_size = params_["LinearSolver"]["SizeBasis"].get<int>( );
    ( ( GMRES<LAD>* )linsolver_ )->InitParameter ( basis_size, "NoPreconditioning" );

#ifdef WITH_ILUPP
    // Prepare ILU preconditioner
    std::string precon_type;
    params_["LinearSolver"]["Preconditioning"].read<std::string>( precon_type );
    use_ilupp_ = ( precon_type == "ILUPP" );
    if ( use_ilupp_ )
    {
        ilupp_.InitParameter ( params_["LinearSolver"]["ILUPP"]["PreprocessingType"].get<int>( ),
                               params_["LinearSolver"]["ILUPP"]["PreconditionerNumber"].get<int>( ),
                               params_["LinearSolver"]["ILUPP"]["MaxMultilevels"].get<int>( ),
                               params_["LinearSolver"]["ILUPP"]["MemFactor"].get<double>( ),
                               params_["LinearSolver"]["ILUPP"]["PivotThreshold"].get<double>( ),
                               params_["LinearSolver"]["ILUPP"]["MinPivot"].get<double>( ) );
        linsolver_->SetupPreconditioner ( ilupp_ );
        ( ( GMRES<LAD>* )linsolver_ )->InitParameter ( basis_size, "RightPreconditioning" );
    }
#endif

    // Setup nonlinear solver parameters
    NonlinearSolverFactory<LAD> NLSFact;
    nls_ = NLSFact.Get ( params_["NonlinearSolver"]["Name"].get<std::string>( ) )->
            params ( res_, matrix_, params_["NonlinearSolver"] );

    // We use our own initial solution - this needs to be indicated to the Newton-solver
    if ( params_["NonlinearSolver"]["Name"].get<std::string>( ) == "Newton" )
        ( ( Newton<LAD>* )nls_ )->InitParameter ( Newton<LAD>::NewtonInitialSolutionOwn );

    nls_->SetOperator ( *this );
    nls_->SetLinearSolver ( *linsolver_ );

    if ( params_["NonlinearSolver"]["ForcingStrategy"].get<std::string>( "" ) == "EisenstatWalker1" )
    {
        nls_forcing_ = new EWForcing<LAD> ( params_["NonlinearSolver"]["InitialForcingValue"].get<Scalar>( ),
                params_["NonlinearSolver"]["MaxForcingValue"].get<Scalar>( ),
                1 );
        ( ( Newton<LAD>* )nls_ )->SetForcingStrategy ( *nls_forcing_ );
    }
    else if ( params_["NonlinearSolver"]["ForcingStrategy"].get<std::string>( "" ) == "EisenstatWalker2" )
    {
        nls_forcing_ = new EWForcing<LAD> ( params_["NonlinearSolver"]["InitialForcingValue"].get<Scalar>( ),
                params_["NonlinearSolver"]["MaxForcingValue"].get<Scalar>( ),
                2,
                params_["NonlinearSolver"]["EW2Gamma"].get<Scalar>( ),
                params_["NonlinearSolver"]["EW2Alpha"].get<Scalar>( ) );
        ( ( Newton<LAD>* )nls_ )->SetForcingStrategy ( *nls_forcing_ );
    }

    if ( params_["NonlinearSolver"]["DampingStrategy"].get<std::string>( "" ) == "Armijo" )
    {
        nls_damping_ = new ArmijoDamping<LAD> ( params_["NonlinearSolver"]["InitialDampingValue"].get<Scalar>( ),
                params_["NonlinearSolver"]["MinimalDampingValue"].get<Scalar>( ),
                params_["NonlinearSolver"]["ArmijoDecrease"].get<Scalar>( ),
                params_["NonlinearSolver"]["SufficientDecrease"].get<Scalar>( ),
                params_["NonlinearSolver"]["MaxDampingIterations"].get<int>( ) );
        ( ( Newton<LAD>* )nls_ )->SetDampingStrategy ( *nls_damping_ );
    }
}

void BloodFlowTutorial::prepare_bc ( )
{

    // Empty the dirichlet vectors
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    // Prepare Dirichlet/Poiseuille boundaries
    std::vector<Scalar> temporal_flow_rate ( poiseuille_materials_.size ( ) );
    for ( int j = 0; j < poiseuille_materials_.size ( ); j++ )
    {
        temporal_flow_rate[j] = evaluate_spline ( j, dt_ * ts_ );
        LOG_INFO ( "Prepare_bc", "Flowrate at bdy " << poiseuille_materials_[j] << ": " << temporal_flow_rate[j] );
    }
    PoiseuilleBC bc[3] = { PoiseuilleBC ( 0, poiseuille_materials_, exponent_, radius_, center_, normal_, temporal_flow_rate ),
                          PoiseuilleBC ( 1, poiseuille_materials_, exponent_, radius_, center_, normal_, temporal_flow_rate ),
                          PoiseuilleBC ( 2, poiseuille_materials_, exponent_, radius_, center_, normal_, temporal_flow_rate ) };
    for ( int i = 0; i < 3; i++ )
    {
        compute_dirichlet_dofs_and_values ( bc[i], space_, i, dirichlet_dofs_, dirichlet_values_ );
    }

    WallBC wall_bc ( wall_bdy_ );

    for ( int i = 0; i < 3; i++ )
    {
        compute_dirichlet_dofs_and_values ( wall_bc, space_, i, dirichlet_dofs_, dirichlet_values_ );
    }

    // Apply BC to initial solution
    LOG_INFO ( "Prepare_bc", "Num dir dofs = " << dirichlet_dofs_.size ( ) );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct solution with dirichlet BC
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
}

void BloodFlowTutorial::calculate_spline ( int bdy_id )
{

    // Prepare flowrates in scaling them with respect to the peak velocity for each Poiseuille profile
    // Scaling the flowrate from volume flow rate to maximal velocity of the Poiseuille profile
    const int N = timestamps_[bdy_id].size ( ) - 1;
    std::vector<Scalar> h ( N, 0. );
    Scalar velo_scaling = ( exponent_[bdy_id] + 2. ) / ( exponent_[bdy_id] * M_PI * radius_[bdy_id] * radius_[bdy_id] );
    for ( int i = 0; i < N; ++i )
    {
        h[i] = timestamps_[bdy_id][i + 1] - timestamps_[bdy_id][i];
        flow_rate_[bdy_id][i] *= velo_scaling;
    }
    flow_rate_[bdy_id][N] *= velo_scaling;

    std::vector<Scalar> w ( N, 0. );
    std::vector<Scalar> q ( N, 0. );
    w[0] = 2. * ( h[0] + h[N - 1] );
    q[0] = 3. * ( flow_rate_[bdy_id][1] / h[0] - flow_rate_[bdy_id][0]*( 1. / h[0] + 1. / h[N - 1] ) + flow_rate_[bdy_id][N - 1] / h[N - 1] );

    for ( int i = 1; i < N; ++i )
    {
        w[i] = 2. * ( h[i] + h[i - 1] );
        q[i] = 3. * ( flow_rate_[bdy_id][i + 1] / h[i] - flow_rate_[bdy_id][i]*( 1. / h[i] + 1. / h[i - 1] ) + flow_rate_[bdy_id][i - 1] / h[i - 1] );
    }

    la::SeqDenseMatrix<Scalar> matrix;
    matrix.Resize ( N, N );
    for ( int i = 0; i < N; ++i )
    {
        matrix ( i, i ) = w[i];
        matrix ( ( i + 1 ) % N, i ) = h[i];
        matrix ( i, ( i + 1 ) % N ) = h[i];
    }

    // Resize the vectors to number of Poiseuille materials
    c_.resize ( poiseuille_materials_.size ( ) );
    b_.resize ( poiseuille_materials_.size ( ) );
    d_.resize ( poiseuille_materials_.size ( ) );

    c_[bdy_id].resize ( N, 0. );
    matrix.Solve ( q, c_[bdy_id] );

    b_[bdy_id].resize ( N, 0. );
    d_[bdy_id].resize ( N, 0. );
    for ( int i = 0; i < N; ++i )
    {
        d_[bdy_id][i] = ( c_[bdy_id][( i + 1 ) % N] - c_[bdy_id][i] ) / ( 3. * h[i] );
        b_[bdy_id][i] = ( flow_rate_[bdy_id][i + 1] - flow_rate_[bdy_id][i] ) / h[i] - h[i]*( c_[bdy_id][( i + 1 ) % N] + 2. * c_[bdy_id][i] ) / 3.;
    }
    if ( smooth_start_up_time_ > 0 )
    {
        smooth_initial_c_.resize ( poiseuille_materials_.size ( ) );
        smooth_initial_d_.resize ( poiseuille_materials_.size ( ) );
        const Scalar t = smooth_start_up_time_;
        const Scalar val = flow_rate_[bdy_id][0];
        const Scalar der = b_[bdy_id][0];
        smooth_initial_c_[bdy_id] = ( 3. * val - t * der ) / ( t * t );
        smooth_initial_d_[bdy_id] = ( t * der - 2. * val ) / ( t * t * t );
    }
}

Scalar BloodFlowTutorial::evaluate_spline ( int bdy_id, Scalar current_time )
{
    Scalar value = 0;
    if ( current_time < smooth_start_up_time_ )
    {
        value = smooth_initial_c_[bdy_id] * current_time * current_time + smooth_initial_d_[bdy_id] * current_time * current_time*current_time;
    }
    else
    {
        const int N = timestamps_[bdy_id].size ( ) - 1;
        assert ( N == flow_rate_[bdy_id].size ( ) - 1 && N == b_[bdy_id].size ( ) && N == c_[bdy_id].size ( ) && N == d_[bdy_id].size ( ) );
        current_time -= smooth_start_up_time_;
        current_time -= std::floor ( current_time / period_ ) * period_;
        assert ( current_time >= 0. );
        assert ( current_time <= period_ );
        int index = 0;
        while ( current_time >= timestamps_[bdy_id][index + 1] && current_time < period_ && index < N )
        {
            ++index;
        }
        current_time -= timestamps_[bdy_id][index];
        value = flow_rate_[bdy_id][index] + ( b_[bdy_id][index] + ( c_[bdy_id][index] + d_[bdy_id][index] * current_time ) * current_time ) * current_time;
    }
    return value;
}

void BloodFlowTutorial::run_time_loop ( )
{
    LOG_INFO ( "timestep", "starting runtime loop" );

    // Visualize initial solution
    visualize ( );
    const Scalar end_time = params_["Instationary"]["Endtime"].get<Scalar>( );

    LOG_INFO ( "timestep", "End time = " << end_time );
    LOG_INFO ( "timestep", "Step length = " << dt_ );

    // Crank-Nicolson time-stepping method. At the beginning of each
    // time-step, the solution from the previous time-step is stored
    // in prev_sol_, which is used in InstationaryFlowAssembler. The
    // variable ts_ is used to keep track of the current
    // time-step. The solution is visualized at the end of each
    // time-step, in order to be able to animate it in Paraview.

    ++ts_;
    while ( ts_ * dt_ <= end_time )
    {
        LOG_INFO ( "timestep", "Solving time step " << ts_ << " (time " << ts_ * dt_ << ")" );

        // Prepare BCs: Poiseuille and Windkessel BC
        prepare_bc ( );

        // Solve nonlinear problem
        solve ( );

        // Get the previous solution
        prev_sol_->CloneFrom ( *sol_ );

        // Visualize results
        if ( ts_ % visu_interval_ == 0 )
        {
            LOG_INFO ( "timestep", "Visualizing solution at time " << ts_ * dt_
                       << " (time-step " << ts_ << ")" );
            visualize ( );
        }

        // Evaluate results
        evaluate ( );
        ++ts_;
    }
}

void BloodFlowTutorial::solve ( )
{

    NonlinearSolverState state = nls_->Solve ( *rhs_, sol_ );

    // Information state, residual norm and number of iterations
    std::cout << "Nonlinear solver ended with state " << state
            << " and residual norm " << nls_->GetResidual ( )
            << " after " << nls_->iter ( ) << " iterations\n";

    LOG_INFO ( "Nonlinear solver residual", nls_->GetResidual ( ) );
    LOG_INFO ( "Nonlinear solver steps", nls_->iter ( ) );
    sol_->UpdateCouplings ( );

}

void BloodFlowTutorial::evaluate ( )
{

    std::vector<Scalar> values_csv;

    // Evaluate pressure and flowrate at all material numbers.
    for ( int m = 0; m < poiseuille_materials_.size ( ); ++m )
    {
        eval_poiseuille_flowrates_[m] = evaluate_bdy_flowrate ( poiseuille_materials_[m] );
        LOG_INFO ( "results poiseuille flowrates", eval_poiseuille_flowrates_[m] );

        eval_poiseuille_pressures_[m] = evaluate_bdy_pressure ( poiseuille_materials_[m] );
        LOG_INFO ( "results poiseuille pressures", eval_poiseuille_pressures_[m] );

        if ( rank_ == MASTER_RANK )
        {
            values_csv.push_back ( eval_poiseuille_flowrates_[m] );
            values_csv.push_back ( eval_poiseuille_pressures_[m] );
        }
    }
    // Iterate over all Windkessel boundaries
    for ( int m = 0; m < windkessel_materials_.size ( ); ++m )
    {
        eval_windkessel_flowrates_[m] = evaluate_bdy_flowrate ( windkessel_materials_[m] );
        LOG_INFO ( "results windkessel flowrates", eval_windkessel_flowrates_[m] );

        eval_windkessel_pressures_[m] = evaluate_bdy_pressure ( windkessel_materials_[m] );
        LOG_INFO ( "results windkessel pressures", eval_windkessel_pressures_[m] );

        if ( rank_ == MASTER_RANK )
        {
            values_csv.push_back ( eval_windkessel_flowrates_[m] );
            values_csv.push_back ( eval_windkessel_pressures_[m] );
        }
    }
    // Write evaluated results in a csv-file.
    if ( rank_ == MASTER_RANK )
    {
        results_writer_.write ( values_csv );
    }
}

Scalar BloodFlowTutorial::evaluate_bdy_flowrate ( const int bdy_material )
{

    double result_local;
    FlowrateIntegrator integrator ( *sol_, bdy_material );
    global_asm_.integrate_scalar_boundary ( space_, integrator, result_local );
    double result;
    MPI_Allreduce ( &result_local, &result, 1, MPI_DOUBLE, MPI_SUM, space_.get_mpi_comm ( ) );

    return result;
}

Scalar BloodFlowTutorial::evaluate_bdy_pressure ( const int bdy_material )
{

    double results_local[2];
    BoundaryAreaIntegrator area_integrator ( bdy_material );
    global_asm_.integrate_scalar_boundary ( space_, area_integrator, results_local[0] );

    BoundaryValueIntegrator integrator ( DIMENSION, *sol_, bdy_material );
    global_asm_.integrate_scalar_boundary ( space_, integrator, results_local[1] );
    double results[2];

    MPI_Allreduce ( &results_local, &results, 2, MPI_DOUBLE, MPI_SUM, space_.get_mpi_comm ( ) );

    assert ( results[0] != 0 );
    return results[1] / results[0];
}

void BloodFlowTutorial::visualize ( )
{

    int num_intervals = 1;
    ParallelCellVisualization<Scalar> visu ( space_, num_intervals, comm_, MASTER_RANK );

    std::stringstream input;

    input << simul_name_ << "_solution" << std::setw ( 5 ) << std::setfill ( '0' ) << ts_;
    std::vector<Scalar> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
          it != mesh_->end ( mesh_->tdim ( ) );
          ++it )
    {
        int temp;
        mesh_->get_attribute_value ( "_sub_domain_", mesh_->tdim ( ),
                                     it->index ( ),
                                     &temp );
        sub_domain.at ( it->index ( ) ) = temp;
    }

    // Create visualization from post-processing vector.
    // Three velocity values, pressure p and Wall Shear Stress WSS
    // Sub domains for parallel computing
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 0 ), "u" );
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 1 ), "v" );
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 2 ), "w" );
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), DIMENSION ), "p" );
    visu.visualize ( EvalWSSmagn<LAD>( space_, *( sol_ ), wall_bdy_, rho_ * nu_ ), "WSS" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );
}
