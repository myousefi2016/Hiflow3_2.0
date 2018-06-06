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

/// \author Jonas Kratzke

#include "fsi_benchmark.h"

static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );
    std::string param_filename;
    if ( argc > 1 )
    {
        param_filename = argv[1];
    }
    else
    {
        std::cout << "Pass XMl parameter file as argument!" << std::endl;
        std::cout << "Usage: ./fsi_benchmark <xml-filename>" << std::endl;
        exit ( -1 );
    }

    // Run application
    try
    {
        PropertyTree config ( param_filename, MASTER_RANK, MPI_COMM_WORLD );
        // Start Logging
        std::ofstream info_file ( ( config["OutputPrefix"].get<std::string>( ) + "_info_log" ).c_str ( ) );
        if ( config["OutputType"].get<std::string>( ) == "infofile" )
        {
            LogKeeper::get_log ( "info" ).set_target ( &info_file );
        }
        else
        {
            LogKeeper::get_log ( "info" ).set_target ( &std::cout );
        }
        std::ofstream debug_log ( ( config["OutputPrefix"].get<std::string>( ) + "_debug_log" ).c_str ( ) );
        LogKeeper::get_log ( "debug" ).set_target ( &debug_log );
        FSIBenchmark application ( config );
        application.run ( );
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

FSIBenchmark::FSIBenchmark ( PropertyTree& config )
: comm_ ( MPI_COMM_WORLD ),
type_ ( NOT_SET ),
config_ ( config ),
simul_name_ ( config_["OutputPrefix"].get<std::string>( ) ),
mesh_ ( 0 ),
refinement_level_ ( 0 ),
newton_forcing_ ( false ),
newton_damping_ ( false ),
matrix_ ( 0 ), sol_ ( 0 ), prev_sol_ ( 0 ), res_ ( 0 ), rhs_ ( 0 ),
results_writer_ ( simul_name_ + "_results.csv" )
{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    initialize_parameters ( );
    initialize_mesh ( );
    initialize_space ( );
    initialize_la_structures ( );
    initialize_postprocessing ( );
}

FSIBenchmark::~FSIBenchmark ( )
{
    if ( newton_forcing_ ) delete ew_forcing_;
    if ( newton_damping_ ) delete armijo_damping_;

    delete sol_;
    delete prev_sol_;
    delete res_;
    delete rhs_;
    delete matrix_;

#ifdef WITH_MUMPS
    linear_solver_->Clear ( );
#endif
    delete linear_solver_;
    delete nls_;
}

void FSIBenchmark::initialize_mesh ( )
{
    assert ( DIM == 2 || DIM == 3 );
    LOG_INFO ( "Problem dimension", DIM );
    refinement_level_ = config_["Mesh"]["RefinementLevel"].get<int>( );
    LOG_INFO ( "Refinement level", refinement_level_ );

    MeshPtr master_mesh ( 0 );
    if ( rank_ == MASTER_RANK )
    {
        // read mesh (sequential, dimension DIM)
        std::string mesh_filename = std::string ( DATADIR ) + config_["Mesh"]["Filename"].get<std::string>( );
        master_mesh = read_mesh_from_file ( mesh_filename, DIM, DIM, 0 );

        for ( int r = 0; r < refinement_level_; ++r )
        {
            master_mesh = master_mesh->refine ( );
        }
        LOG_INFO ( "mesh", "Initial refinement level = " << refinement_level_ );
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

void FSIBenchmark::initialize_parameters ( )
{
    Scalar young, poisson;
    if ( config_["Application"]["Type"].get<std::string>( ) == "Channel" )
    {
        type_ = Channel;
    }
    else if ( config_["Application"]["Type"].get<std::string>( ) == "Couette" )
    {
        type_ = Couette;
    }
    else
    {
        LOG_INFO ( "Parameters", "Benchmark type not defined: Exiting!" );
        exit ( -1 );
    }
    rho_ = config_["Application"]["Density"].get<double>( );
    nu_ = config_["Application"]["Viscosity"].get<double>( );
    Vm_ = config_["Application"]["Speed"].get<double>( );
    mesh_diff_ = config_["Application"]["MeshDiffFac"].get<double>( );
    linear_ = ( config_["Application"]["LinearElast"].get<int>( ) == 1 );
    LOG_INFO ( "Parameters", "Linear: " << linear_ );

    poisson = config_["Application"]["Poisson"].get<double>( );
    young = config_["Application"]["Young"].get<double>( );
    mu_ = 0.5 * young / ( 1. + poisson );
    lambda_ = 2. * mu_ * poisson / ( 1. - 2. * poisson );

    if ( type_ == Channel )
    {
        H_ = config_["Application"]["InflowHeight"].get<double>( );
    }
    else if ( type_ == Couette )
    {
        config_["Application"]["Radii"].read<double>( radii_ );
        assert ( radii_.size ( ) == 3 );
        // Computation of the analytical solution of the displacement at the interface
        LOG_INFO ( "Parameters", "Vm_: " << Vm_ );
        LOG_INFO ( "Parameters", "rho_: " << rho_ );
        LOG_INFO ( "Parameters", "nu_: " << nu_ );
        LOG_INFO ( "Parameters", "radii_: " << radii_[0] << ", " << radii_[1] << ", " << radii_[2] );
        LOG_INFO ( "Parameters", "mu_: " << mu_ );

        Um_ = Vm_ * rho_ * nu_ * radii_[0] / ( radii_[1] * radii_[1] - radii_[0] * radii_[0] );
        Um_ *= radii_[1] * ( radii_[1] * radii_[1] - radii_[2] * radii_[2] ) / ( mu_ * radii_[2] * radii_[2] );
        LOG_INFO ( "Parameters", "Analytical interface displacement: " << Um_ );
    }

    fluid_mat_ = config_["Domain"]["FluidMaterial"].get<int>( );
    solid_mat_ = config_["Domain"]["SolidMaterial"].get<int>( );
    inflow_mat_ = config_["Boundary"]["InflowMaterial"].get<int>( );
    outflow_mat_ = config_["Boundary"]["OutflowMaterial"].get<int>( );
    obstacle_mat_ = config_["Boundary"]["ObstacleMaterial"].get<int>( );
    interface_mat_ = config_["Boundary"]["InterfaceMaterial"].get<int>( );
    side_mat_ = config_["Boundary"]["SideMaterial"].get<int>( );

    ts_ = 1;
    Dt_ = config_["Instationary"]["Timestep"].get<double>( );
    dt_ = Dt_ / 10.;
    t_ = dt_;
    start_up_time_ = config_["Instationary"]["SmoothStartup"].get<double>( );
    theta_ = config_["Instationary"]["Theta"].get<double>( );
    visu_intervall_ = config_["Instationary"]["VisuIntervall"].get<int>( );
}

void FSIBenchmark::initialize_space ( )
{
    // setup space P2-P1-P2 elements
    std::vector<int> degrees ( 2 * DIM + 1, 2 );
    degrees[DIM] = 1;
    // The handling of pressure dofs on the structure domain requires
    // linear ansatz functions for the pressure!
    assert ( degrees[DIM] == 1 );
    space_.Init ( degrees, *mesh_ );

    LOG_INFO ( "dofs", "Total number: " << space_.dof ( ).ndofs_global ( ) << ", local: " << space_.dof ( ).ndofs_on_sd ( rank_ ) );
}

void FSIBenchmark::initialize_la_structures ( )
{
    // Initialize linear algebra structures
    couplings_.Clear ( );
    couplings_.Init ( comm_, space_.dof ( ) );

    // compute matrix graph to build matrix structure
    std::vector<int> diagonal_rows, diagonal_cols,
            off_diagonal_rows, off_diagonal_cols;

    std::vector < std::vector<bool> > coupling_vars;
    // First entry (row): Test variable
    // Second entry (column): Trial variable
    coupling_vars.resize ( 2 * DIM + 1 );
    // NSE couplings:
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
    coupling_vars[DIM].push_back ( true );
    // Displacement couplings:
    for ( int i = 0; i < DIM; ++i )
    {
        for ( int j = 0; j < DIM; ++j )
        {
            coupling_vars[i].push_back ( true );
        }
    }
    for ( int i = 0; i < DIM; ++i )
    {
        coupling_vars[DIM].push_back ( true );
    }
    for ( int i = DIM + 1; i < 2 * DIM + 1; ++i )
    {
        for ( int j = 0; j < DIM; ++j )
        {
            coupling_vars[i].push_back ( true );
        }
        coupling_vars[i].push_back ( false );
        for ( int j = 0; j < DIM; ++j )
        {
            coupling_vars[i].push_back ( true );
        }
    }

    InitStructure ( space_, &diagonal_rows, &diagonal_cols,
                    &off_diagonal_rows, &off_diagonal_cols, &coupling_vars );
    couplings_.InitializeCouplings ( off_diagonal_rows, off_diagonal_cols );

    // Initialize matrices and vectors
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity, &coupling_vars );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               config_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( config_["LinearAlgebra"] );
    matrix_->Init ( comm_, couplings_ );
    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );
    matrix_->Zeros ( );

    CoupledVectorFactory<Scalar> CoupVecFact;
    sol_ = CoupVecFact.Get (
                             config_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( config_["LinearAlgebra"] );
    sol_->Init ( comm_, couplings_ );
    sol_->InitStructure ( );
    sol_->Zeros ( );
    prev_sol_ = CoupVecFact.Get (
                                  config_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( config_["LinearAlgebra"] );
    prev_sol_->Init ( comm_, couplings_ );
    prev_sol_->InitStructure ( );
    prev_sol_->Zeros ( );
    res_ = CoupVecFact.Get (
                             config_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( config_["LinearAlgebra"] );
    res_->Init ( comm_, couplings_ );
    res_->InitStructure ( );
    res_->Zeros ( );
    rhs_ = CoupVecFact.Get (
                             config_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( config_["LinearAlgebra"] );
    rhs_->Init ( comm_, couplings_ );
    rhs_->InitStructure ( );
    rhs_->Zeros ( );
}

void FSIBenchmark::initialize_linear_solver ( )
{
    // Direct solver if compiled
#ifdef WITH_MUMPS
    linear_solver_ = new MumpsSolver<LAD, MumpsStructureD>( );
    linear_solver_->InitParameter ( comm_ );
#else

    LinearSolverFactory<LAD> LinSolFact;
    linear_solver_ = LinSolFact.Get (
                                      config_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( config_["LinearSolver"] );
    linear_solver_->SetPrintLevel ( 1 );
#    ifdef WITH_ILUPP
    use_ilupp_ = ( config_["LinearSolver"]["Preconditioner"].get<std::string>( ) == "ilupp" );
    if ( use_ilupp_ )
    {
        LOG_INFO ( "solver", "using ILUPP" );
        ilupp_.InitParameter ( config_["LinearSolver"]["ILUPP"]["PreprocessingType"].get<int>( ),
                               config_["LinearSolver"]["ILUPP"]["PreconditionerNumber"].get<int>( ),
                               config_["LinearSolver"]["ILUPP"]["MaxMultilevels"].get<int>( ),
                               config_["LinearSolver"]["ILUPP"]["MemFactor"].get<double>( ),
                               config_["LinearSolver"]["ILUPP"]["PivotThreshold"].get<double>( ),
                               config_["LinearSolver"]["ILUPP"]["MinPivot"].get<double>( ) );
        linear_solver_->SetupPreconditioner ( ilupp_ );
        ( ( GMRES<LAD>* )linear_solver_ )->InitParameter ( ( ( GMRES<LAD>* )linear_solver_ )->size_basis ( ), "RightPreconditioning" );
    }
#    else
    ( ( GMRES<LAD>* )linear_solver_ )->InitParameter ( ( ( GMRES<LAD>* )linear_solver_ )->size_basis ( ), "NoPreconditioning" );
#    endif
    linear_solver_->SetupOperator ( *matrix_ );
#endif
}

void FSIBenchmark::initialize_nonlinear_solver ( )
{
    NonlinearSolverFactory<LAD> NLSFact;
    nls_ = NLSFact.Get ( config_["NonlinearSolver"]["Name"].get<std::string>( ) )->
            params ( res_, matrix_, config_["NonlinearSolver"] );

    // we use our own initial solution -- this needs to be indicated
    // to the Newton-solver
    if ( config_["NonlinearSolver"]["Name"].get<std::string>( ) == "Newton" )
        ( ( Newton<LAD>* )nls_ )->InitParameter ( Newton<LAD>::NewtonInitialSolutionOwn );

    nls_->SetOperator ( *this );
    nls_->SetLinearSolver ( *linear_solver_ );

    // get forcing strategy parameters from param file
    forcing_strategy_ = config_["NonlinearSolver"]["ForcingStrategy"].get<std::string>( );
    eta_initial_ = config_["NonlinearSolver"]["InitialValueForcingTerm"].get<double>( );
    eta_max_ = config_["NonlinearSolver"]["MaxValueForcingTerm"].get<double>( );
    gamma_EW2_ = config_["NonlinearSolver"]["GammaParameterEW2"].get<double>( );
    alpha_EW2_ = config_["NonlinearSolver"]["AlphaParameterEW2"].get<double>( );

    // setup forcing strategy object within the nonlinear solver
    if ( forcing_strategy_ == "EisenstatWalker1" )
    {
        ew_forcing_ = new EWForcing<LAD> ( eta_initial_, eta_max_, 1 );
        ( ( Newton<LAD>* )nls_ )->SetForcingStrategy ( *ew_forcing_ );
        newton_forcing_ = true;
    }
    else if ( forcing_strategy_ == "EisenstatWalker2" )
    {
        ew_forcing_ = new EWForcing<LAD> ( eta_initial_, eta_max_, 2, gamma_EW2_, alpha_EW2_ );
        ( ( Newton<LAD>* )nls_ )->SetForcingStrategy ( *ew_forcing_ );
        newton_forcing_ = true;
    }

    //get damping strategy parameters from param file
    damping_strategy_ = config_["NonlinearSolver"]["DampingStrategy"].get<std::string>( );
    theta_initial_ = config_["NonlinearSolver"]["ThetaInitial"].get<double>( );
    theta_min_ = config_["NonlinearSolver"]["ThetaMinimal"].get<double>( );
    armijo_dec_ = config_["NonlinearSolver"]["ArmijoDecrease"].get<double>( );
    suff_dec_ = config_["NonlinearSolver"]["SufficientDecrease"].get<double>( );
    max_armijo_ite_ = config_["NonlinearSolver"]["MaxArmijoIteration"].get<int>( );
    if ( damping_strategy_ == "Armijo" )
    {
        armijo_damping_ = new ArmijoDamping<LAD>( theta_initial_, theta_min_, armijo_dec_, suff_dec_, max_armijo_ite_ );
        ( ( Newton<LAD>* )nls_ )->SetDampingStrategy ( *armijo_damping_ );
    }
    else
    {
        ( ( Newton<LAD>* )nls_ )->InitParameter ( Newton<LAD>::NewtonDampingStrategyNone );
    }
}

void FSIBenchmark::initialize_postprocessing ( )
{
    config_["Application"]["EvaluationPoint"].read<double>( eval_point_ );
    std::vector<std::string> results_names;
    results_names.push_back ( "time" );
    results_names.push_back ( "displacement_X" );
    results_names.push_back ( "displacement_Y" );
    results_names.push_back ( "obstacle_force_X" );
    results_names.push_back ( "obstacle_force_Y" );
    if ( type_ == Couette )
    {
        results_names.push_back ( "force_magn_err" );
        results_names.push_back ( "vel_L2_err" );
        results_names.push_back ( "press_L2_err" );
        results_names.push_back ( "disp_L2_err" );
        results_names.push_back ( "diffusion_L2_err" );
    }
    results_writer_.Init ( results_names );
}

/// \brief Definition of Dirichlet velocity boundary condition
///        for the FSI benchmark scenario
///
/// \param Vm       peak inflow velocity
/// \param H        height of the channel
/// \param factor   velocity scaling factor, e.g. for the smooth startup
/// \param var      FEM variable
/// \param inflow   inflow material number
/// \param outflow  outflow material number

struct FSIBenchmarkVeloDirichletBC2d
{

    FSIBenchmarkVeloDirichletBC2d (
                                    const Scalar Vm,
                                    const Scalar H,
                                    const Scalar factor,
                                    const int var,
                                    const int inflow,
                                    const int outflow )
    : Vm_ ( Vm ), H_ ( H ), factor_ ( factor ), var_ ( var ), inflow_ ( inflow ), outflow_ ( outflow )
    {
        assert ( var_ == 0 || var_ == 1 );
    }

    std::vector<Scalar> evaluate ( const Entity& face, const std::vector<Coord >& coords_on_face ) const
    {
        std::vector<Scalar> values;

        const int material_num = face.get_material_number ( );

        Scalar pi = acos ( 0.0 )*2.0;

        if ( material_num != outflow_ )
        {
            values.resize ( coords_on_face.size ( ) );

            // loop over points on the face
            for ( int i = 0; i < coords_on_face.size ( ); ++i )
            {

                // evaluate dirichlet function at each point
                const Coord pt = coords_on_face[i];
                if ( material_num == inflow_ )
                {
                    if ( var_ == 0 )
                    { // x-component
                        values[i] = factor_ * 4. * Vm_ * pt[1] * ( H_ - pt[1] ) / ( H_ * H_ );
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
    const Scalar Vm_, H_, factor_;
    const int var_, inflow_, outflow_;
};

/// \brief Definition of Dirichlet velocity boundary condition
///        for the FSI Couette scenario
///
/// \param Vm       peak inflow velocity
/// \param factor   velocity scaling factor, e.g. for the smooth startup
/// \param radii    inner, interface and outer radius of the circular geometry
/// \param var      FEM variable
/// \param inflow   inflow material number

struct FSICouetteVeloDirichletBC2d
{

    FSICouetteVeloDirichletBC2d (
                                  const Scalar Vm,
                                  const Scalar factor,
                                  const std::vector<Scalar>& radii,
                                  const int var,
                                  const int inflow )
    : Vm_ ( Vm ), factor_ ( factor ), radii_ ( radii ), var_ ( var ), inflow_ ( inflow )
    {
        assert ( var_ == 0 || var_ == 1 );
    }

    std::vector<Scalar> evaluate ( const Entity& face, const std::vector<Coord >& coords_on_face ) const
    {
        std::vector<Scalar> values ( coords_on_face.size ( ), 0. );

        const int material_num = face.get_material_number ( );

        if ( material_num == inflow_ )
        {
            // loop over points on the face
            for ( int i = 0; i < coords_on_face.size ( ); ++i )
            {
                // evaluate dirichlet function at each point
                const Coord pt = coords_on_face[i];
                Scalar radius = pt[0] * pt[0] + pt[1] * pt[1];
                radius = sqrt ( radius );
                if ( var_ == 0 )
                { // x-component
                    values[i] = -factor_ * Vm_ * pt[1] / radius;
                    values[i] *= radii_[0] / ( radii_[1] * radii_[1] - radii_[0] * radii_[0] )*( radii_[1] * radii_[1] / radius - radius );
                }
                else if ( var_ == 1 )
                { // y-component
                    values[i] = factor_ * Vm_ * pt[0] / radius;
                    values[i] *= radii_[0] / ( radii_[1] * radii_[1] - radii_[0] * radii_[0] )*( radii_[1] * radii_[1] / radius - radius );
                }
                else
                {
                    assert ( false );
                }
            }
        }
        return values;
    }
    const Scalar Vm_, factor_;
    const std::vector<Scalar>& radii_;
    const int var_, inflow_;
};

/// \brief Definition of Dirichlet displacement boundary condition
///        for the FSI Couette scenario
///
/// \param Um       exact interface displacement
/// \param radii    inner, interface and outer radius of the circular geometry
/// \param var      FEM variable
/// \param inflow   material number of outer boundary

struct FSIBenchmarkExactDispDirichletBC2d
{

    FSIBenchmarkExactDispDirichletBC2d (
                                         const Scalar Um,
                                         const std::vector<double>& radii,
                                         const int var,
                                         const int outer_mat )
    : Um_ ( Um ), radii_ ( radii ), var_ ( var ), outer_mat_ ( outer_mat )
    {
        assert ( var_ == 3 || var_ == 4 );
    }

    std::vector<Scalar> evaluate ( const Entity& face, const std::vector<Coord >& coords_on_face ) const
    {
        std::vector<Scalar> values ( coords_on_face.size ( ), 0. );

        const int material_num = face.get_material_number ( );

        if ( material_num == outer_mat_ )
        {
            // loop over points on the face
            for ( int i = 0; i < coords_on_face.size ( ); ++i )
            {
                // evaluate dirichlet function at each point
                const Coord pt = coords_on_face[i];
                Scalar radius = pt[0] * pt[0] + pt[1] * pt[1];
                radius = std::sqrt ( radius );
                double exact_disp = -Um_ * radii_[1] * ( radius - radii_[2] * radii_[2] / radius ) / ( radii_[1] * radii_[1] - radii_[2] * radii_[2] );

                // compute angle coordinate
                Scalar phi = std::atan2 ( pt[1], pt[0] );
                // add the radial displacement
                phi += exact_disp / radius;

                if ( var_ == 3 )
                { // x-component
                    values[i] = radius * std::cos ( phi ) - pt[0];
                }
                else if ( var_ == 4 )
                { // y-component
                    values[i] = radius * std::sin ( phi ) - pt[1];
                }
                else
                {
                    assert ( false );
                }
            }
        }
        return values;
    }
    const Scalar Um_;
    const std::vector<double>& radii_;
    const int var_, outer_mat_;
};

/// \brief Definition of Dirichlet displacement boundary condition

struct FSIBenchmarkDispDirichletBC2d
{

    FSIBenchmarkDispDirichletBC2d ( )
    {
    }

    std::vector<Scalar> evaluate ( const Entity& face, const std::vector<Coord >& coords_on_face ) const
    {
        std::vector<Scalar> values ( coords_on_face.size ( ), 0. );
        return values;
    }
};

void FSIBenchmark::compute_bc ( )
{
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    // Smooth start up factor
    Scalar factor = 1.;
    if ( t_ < start_up_time_ )
    {
        factor = ( 1. - cos ( t_ * acos ( 0.0 )*2.0 / start_up_time_ ) ) * 0.5;
    }
    LOG_INFO ( "timestep", "Smooth start factor = " << factor );

    if ( type_ == Channel )
    {
        FSIBenchmarkVeloDirichletBC2d vbc[2] = { FSIBenchmarkVeloDirichletBC2d ( Vm_, H_, factor, 0, inflow_mat_, outflow_mat_ ),
                                                FSIBenchmarkVeloDirichletBC2d ( Vm_, H_, factor, 1, inflow_mat_, outflow_mat_ ) };

        for ( int v = 0; v < DIM; ++v )
        {
            compute_dirichlet_dofs_and_values ( vbc[v], space_, v,
                                                dirichlet_dofs_, dirichlet_values_ );
        }

        FSIBenchmarkDispDirichletBC2d dbc;

        for ( int v = DIM + 1; v < 2 * DIM + 1; ++v )
        {
            compute_dirichlet_dofs_and_values ( dbc, space_, v,
                                                dirichlet_dofs_, dirichlet_values_ );
        }
    }
    else if ( type_ == Couette )
    {
        FSICouetteVeloDirichletBC2d vbc[2] = { FSICouetteVeloDirichletBC2d ( Vm_, factor, radii_, 0, inflow_mat_ ),
                                              FSICouetteVeloDirichletBC2d ( Vm_, factor, radii_, 1, inflow_mat_ ) };

        for ( int v = 0; v < DIM; ++v )
        {
            compute_dirichlet_dofs_and_values ( vbc[v], space_, v,
                                                dirichlet_dofs_, dirichlet_values_ );
        }

        FSIBenchmarkExactDispDirichletBC2d dbc[2] = { FSIBenchmarkExactDispDirichletBC2d ( Um_, radii_, 3, side_mat_ ),
                                                     FSIBenchmarkExactDispDirichletBC2d ( Um_, radii_, 4, side_mat_ ) };

        for ( int v = DIM + 1; v < 2 * DIM + 1; ++v )
        {
            compute_dirichlet_dofs_and_values ( dbc[v - DIM - 1], space_, v,
                                                dirichlet_dofs_, dirichlet_values_ );
        }
    }

    fix_solid_pressure_dofs ( dirichlet_dofs_, dirichlet_values_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ),
                          dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
}

void FSIBenchmark::fix_solid_pressure_dofs ( std::vector<int>& dirichlet_dofs,
                                             std::vector<LAD::DataType>& dirichlet_values )
{

    // This procedure requires linear ansatz functions for the pressure!!!

    // Get dofs and coords
    const DofPartition<Scalar>& dof = space_.dof ( );

    int fixed_pressure_dofs = 0;
    // Iterate vertices
    for ( EntityIterator it_v = mesh_->begin ( 0 );
          it_v != mesh_->end ( 0 );
          ++it_v )
    {

        // iterate incident cells
        bool all_inc_cells_solid = true;
        for ( IncidentEntityIterator it_inc_cell = it_v->begin_incident ( DIM );
              it_inc_cell != it_v->end_incident ( DIM );
              ++it_inc_cell )
        {
            if ( fluid_mat_ == mesh_->get_material_number ( DIM, it_inc_cell->index ( ) ) )
            {
                all_inc_cells_solid = false;
            }
        }

        // if all cells are solid cells
        if ( all_inc_cells_solid )
        {
            Entity inc_cell = *( it_v->begin_incident ( DIM ) );
            // find vertex number w.r.t. cell
            int vertex_number = 0;
            for ( IncidentEntityIterator it_vertex = inc_cell.begin_incident ( 0 );
                  it_vertex->index ( ) != it_v->index ( );
                  ++it_vertex, ++vertex_number )
            {
                assert ( vertex_number < inc_cell.num_vertices ( ) );
            }
            // get dof on subentity
            std::vector<int> dofs_on_current_vertex;
            dof.get_dofs_on_subentity ( DIM, inc_cell.index ( ), 0, vertex_number, dofs_on_current_vertex );
            // add dof to diri dofs and values
            if ( dof.owner_of_dof ( dofs_on_current_vertex[0] ) == rank_ )
            {
                dirichlet_dofs.push_back ( dofs_on_current_vertex[0] );
                Scalar value = 0;
                dirichlet_values.push_back ( value );
                ++fixed_pressure_dofs;
            }
        }
    }
    LOG_INFO ( "dof_fix", "Fixed pressure dofs in structure domain: " << fixed_pressure_dofs );
}

void FSIBenchmark::run ( )
{
    compute_bc ( );
    initialize_linear_solver ( );
    initialize_nonlinear_solver ( );

    Scalar end_time = config_["Instationary"]["Endtime"].get<double>( );

    LOG_INFO ( "timestep", "End time = " << end_time );
    LOG_INFO ( "timestep", "Step length, initial: " << dt_ << ", final: " << Dt_ );

    Timer timerT;
    timerT.start ( );

    while ( t_ <= end_time )
    {
        LOG_INFO ( "timestep", "Solving time step " << ts_ << " at time " << t_ );
        compute_bc ( );

        Timer timerNL;
        timerNL.start ( );
        solve_nonlinear_system ( );
        timerNL.stop ( );
        LOG_INFO ( "nonlinear", "Solving nonlinear system time: " << timerNL.get_duration ( ) );

        sol_->UpdateCouplings ( );
        prev_sol_->CloneFrom ( *sol_ );

        LOG_INFO ( "timestep", "Postprocessing solution at time " << t_
                   << " (time-step " << ts_ << ")" );

        postprocessing ( );

        ++ts_;
        if ( dt_ < Dt_ )
        {
            dt_ += dt_;
        }
        else
        {
            dt_ = Dt_;
        }
        t_ += dt_;

        LogKeeper::get_log ( "info" ).flush ( );
    }
    timerT.stop ( );
    LOG_INFO ( "timestep", "Instationary solver time: " << timerT.get_duration ( ) );
}

void FSIBenchmark::EvalFunc ( const CVector& u, CVector* F )
{
    assemble_residual ( u, F );
}

void FSIBenchmark::EvalGrad ( const CVector& u, CMatrix* DF )
{
    assemble_matrix ( u, DF );
}

void FSIBenchmark::assemble_matrix ( const CVector& u, CMatrix* DF )
{

    FSIBenchmarkAssembler local_asm_ ( rho_, nu_, lambda_, mu_, mesh_diff_, linear_ );
    local_asm_.SetNewtonSolution ( &u );
    local_asm_.SetTimeSolution ( prev_sol_ );
    local_asm_.SetMaterials ( fluid_mat_, solid_mat_, outflow_mat_ );
    local_asm_.SetTimeStepping ( theta_, dt_ );

    global_asm_.assemble_matrix ( space_, local_asm_, *DF );

    if ( type_ == Channel )
    {
        global_asm_.should_reset_assembly_target ( false );
        global_asm_.assemble_matrix_boundary ( space_, local_asm_, *DF );
        global_asm_.should_reset_assembly_target ( true );
    }

    // correct BC -- set Dirichlet rows
    if ( !dirichlet_dofs_.empty ( ) )
    {
        DF->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
    }

#ifdef WITH_MUMPS
    {
        linear_solver_->Clear ( );
        linear_solver_->InitParameter ( comm_ );
        linear_solver_->SetupOperator ( *DF );
    }
#else
#    ifdef WITH_ILUPP
    {
        if ( use_ilupp_ )
        {
            ilupp_.SetupOperator ( *DF );
        }
    }
#    endif
#endif
}

void FSIBenchmark::assemble_residual ( const CVector& u, CVector* F )
{

    FSIBenchmarkAssembler local_asm_ ( rho_, nu_, lambda_, mu_, mesh_diff_, linear_ );
    local_asm_.SetNewtonSolution ( &u );
    local_asm_.SetTimeSolution ( prev_sol_ );
    local_asm_.SetMaterials ( fluid_mat_, solid_mat_, outflow_mat_ );
    local_asm_.SetTimeStepping ( theta_, dt_ );

    global_asm_.assemble_vector ( space_, local_asm_, *F );

    if ( type_ == Channel )
    {
        global_asm_.should_reset_assembly_target ( false );
        global_asm_.assemble_vector_boundary ( space_, local_asm_, *F );
        global_asm_.should_reset_assembly_target ( true );
    }

    // correct BC -- set Dirichlet rows
    if ( !dirichlet_dofs_.empty ( ) )
    {
        std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
        F->SetValues ( &dirichlet_dofs_.front ( ), dirichlet_dofs_.size ( ), &zeros.front ( ) );
    }
}

void FSIBenchmark::solve_nonlinear_system ( )
{
    NonlinearSolverState state = nls_->Solve ( *rhs_, sol_ );
    LOG_INFO ( "nonlinear", "Nonlinear solver ended with state " << state
               << " and residual norm " << nls_->GetResidual ( )
               << " after " << nls_->iter ( ) << " iterations" );
}

void FSIBenchmark::postprocessing ( )
{
    // Visualization
    if ( ts_ % visu_intervall_ == 0 )
    {
        std::stringstream input;
        input << simul_name_ << std::setw ( 5 ) << std::setfill ( '0' ) << ts_;
        visualize_solution ( *sol_, input.str ( ) );
    }

    // Calculation of tip position
    std::vector<Scalar> values ( 1 + 2 * DIM, 0 );
    values[0] = t_;
    PointEvaluator<Scalar> evaluator ( space_ );
    for ( int v = 0; v < DIM; ++v )
    {
        EvalFeFunction<LAD> fe_eval ( space_, *sol_, v + DIM + 1 );
        evaluator.evaluate_fun_global ( fe_eval, eval_point_, values[v + 1], comm_ );
        LOG_INFO ( "results", "Tip displacement, dir" << v << ": " << values[v + 1] );
    }

    // Calculation of interface forces
    DGGlobalAssembler<Scalar> dg_global_asm;
    dg_global_asm.should_reset_assembly_target ( true );

    Scalar local_force[2] = { 0, 0 };
    for ( int v = 0; v < DIM; ++v )
    {
        std::vector<Scalar> force_values;

        FSIForceIntegrator force_integrator_bar ( *sol_, interface_mat_, fluid_mat_, rho_*nu_, v );
        dg_global_asm.assemble_interface_scalar ( space_, force_integrator_bar, force_values );

        local_force[v] = std::accumulate ( force_values.begin ( ), force_values.end ( ), 0. );

        if ( type_ == Channel )
        {
            force_values.clear ( );
            FSIForceIntegrator force_integrator_obstacle ( *sol_, obstacle_mat_, fluid_mat_, rho_*nu_, v );
            dg_global_asm.assemble_interface_scalar ( space_, force_integrator_obstacle, force_values );

            local_force[v] += std::accumulate ( force_values.begin ( ), force_values.end ( ), 0. );
        }
        LOG_INFO ( "results", "Local force, dir" << v << ": " << local_force[v] );
    }

    // Sum up local values from each process
    MPI_Reduce ( &local_force[0], &values[DIM + 1], 2, MPI_DOUBLE, MPI_SUM, 0, space_.get_mpi_comm ( ) );

    if ( type_ == Couette )
    {
        std::vector<Scalar> force_values;

        FSIForceMagnIntegrator force_integrator ( *sol_, interface_mat_, fluid_mat_, rho_ * nu_ );
        dg_global_asm.assemble_interface_scalar ( space_, force_integrator, force_values );

        Scalar local_f = std::accumulate ( force_values.begin ( ), force_values.end ( ), 0. );
        Scalar global_f = 0;
        // Sum up local values from each process
        MPI_Reduce ( &local_f, &global_f, 1, MPI_DOUBLE, MPI_SUM, 0, space_.get_mpi_comm ( ) );
        // analytical calculated shear stress
        Scalar exact_force = 4. * acos ( 0.0 )*2.0 * rho_ * nu_ * Vm_ * radii_[1] * radii_[0] / ( radii_[1] * radii_[1] - radii_[0] * radii_[0] );
        // write out the difference
        // Since the DG assembler integrates every facet twice, we have to divide by 2
        values.push_back ( std::abs ( exact_force - global_f / 2. ) );

        L2_err_vel_.clear ( );
        Scalar local_flow_error, global_flow_error;
        flowL2error flow_error_integrator ( *sol_, fluid_mat_, Vm_, radii_ );
        global_asm_.assemble_scalar ( space_, flow_error_integrator, L2_err_vel_ );
        local_flow_error = std::accumulate ( L2_err_vel_.begin ( ), L2_err_vel_.end ( ), 0. );
        MPI_Reduce ( &local_flow_error, &global_flow_error, 1, MPI_DOUBLE, MPI_SUM, 0, space_.get_mpi_comm ( ) );
        values.push_back ( global_flow_error );

        L2_err_press_.clear ( );
        Scalar local_press_error, global_press_error;
        pressureL2error press_error_integrator ( *sol_, fluid_mat_, Vm_, rho_, radii_ );
        global_asm_.assemble_scalar ( space_, press_error_integrator, L2_err_press_ );
        local_press_error = std::accumulate ( L2_err_press_.begin ( ), L2_err_press_.end ( ), 0. );
        MPI_Reduce ( &local_press_error, &global_press_error, 1, MPI_DOUBLE, MPI_SUM, 0, space_.get_mpi_comm ( ) );
        values.push_back ( global_press_error );

        L2_err_disp_.clear ( );
        Scalar local_disp_error, global_disp_error;
        dispL2error disp_error_integrator ( *sol_, solid_mat_, Um_, radii_ );
        global_asm_.assemble_scalar ( space_, disp_error_integrator, L2_err_disp_ );
        local_disp_error = std::accumulate ( L2_err_disp_.begin ( ), L2_err_disp_.end ( ), 0. );
        MPI_Reduce ( &local_disp_error, &global_disp_error, 1, MPI_DOUBLE, MPI_SUM, 0, space_.get_mpi_comm ( ) );
        values.push_back ( global_disp_error );

        L2_err_diff_.clear ( );
        Scalar local_diffusion_error, global_diffusion_error;
        diffusionL2error diffusion_error_integrator ( *sol_, fluid_mat_, Um_, radii_ );
        global_asm_.assemble_scalar ( space_, diffusion_error_integrator, L2_err_diff_ );
        local_diffusion_error = std::accumulate ( L2_err_diff_.begin ( ), L2_err_diff_.end ( ), 0. );
        MPI_Reduce ( &local_diffusion_error, &global_diffusion_error, 1, MPI_DOUBLE, MPI_SUM, 0, space_.get_mpi_comm ( ) );
        values.push_back ( global_diffusion_error );
    }

    if ( rank_ == MASTER_RANK )
    {
        results_writer_.write ( values );
    }
}

void FSIBenchmark::visualize_solution ( CVector& u, std::string const& filename ) const
{
    int num_intervals = 2;
    ParallelCellVisualization<Scalar> visu ( space_, num_intervals, comm_, 0 );

    std::vector<Scalar> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<Scalar> mat_num ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
          it != mesh_->end ( mesh_->tdim ( ) );
          ++it )
    {
        int temp;
        mesh_->get_attribute_value ( "_sub_domain_", mesh_->tdim ( ),
                                     it->index ( ),
                                     &temp );
        sub_domain.at ( it->index ( ) ) = temp;
        mat_num.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
    }

    std::stringstream input;
    input << filename;

    if ( num_partitions_ > 1 )
        input << ".pvtu";
    else
        input << ".vtu";

    visu.visualize ( EvalFeFunction<LAD>( space_, u, 0 ), "Vx" );
    visu.visualize ( EvalFeFunction<LAD>( space_, u, 1 ), "Vy" );
    visu.visualize ( EvalFeFunction<LAD>( space_, u, DIM ), "p" );
    visu.visualize ( EvalFeFunction<LAD>( space_, u, DIM + 1 ), "Ux" );
    visu.visualize ( EvalFeFunction<LAD>( space_, u, DIM + 2 ), "Uy" );

    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.visualize_cell_data ( mat_num, "material" );

    if ( type_ == Couette )
    {
        visu.visualize_cell_data ( L2_err_vel_, "L2-vel" );
        visu.visualize_cell_data ( L2_err_disp_, "L2-disp" );
        visu.visualize_cell_data ( L2_err_press_, "L2-press" );
        visu.visualize_cell_data ( L2_err_diff_, "L2-diff" );
    }

    visu.write ( input.str ( ) );
}

void FSIBenchmarkAssembler::operator() ( const Element<Scalar>& element, const Quadrature<Scalar>& quadrature, LocalMatrix& lm )
{
    initialize_for_element ( element, quadrature );

    const int mat = element.get_cell ( ).get_material_number ( );
    assert ( ( fluid_mat_ == mat ) || ( solid_mat_ == mat ) );

    const int num_q = num_quadrature_points ( );

    if ( fluid_mat_ == mat )
    {

        // recompute previous solution values
        for ( int d = 0; d < DIM; ++d )
        {
            vel_ns_[d].clear ( );
            vel_ts_[d].clear ( );
            grad_vel_ns_[d].clear ( );

            evaluate_fe_function ( *u_, d, vel_ns_[d] );
            evaluate_fe_function ( *prev_sol_, d, vel_ts_[d] );
            evaluate_fe_function_gradients ( *u_, d, grad_vel_ns_[d] );

            disp_ns_[d].clear ( );
            disp_ts_[d].clear ( );
            grad_disp_ns_[d].clear ( );
            grad_disp_ts_[d].clear ( );

            evaluate_fe_function ( *u_, d + 1 + DIM, disp_ns_[d] );
            evaluate_fe_function ( *prev_sol_, d + 1 + DIM, disp_ts_[d] );
            evaluate_fe_function_gradients ( *u_, d + 1 + DIM, grad_disp_ns_[d] );
            evaluate_fe_function_gradients ( *prev_sol_, d + 1 + DIM, grad_disp_ts_[d] );
        }
        // recompute pressure
        p_ns_.clear ( );
        evaluate_fe_function ( *u_, DIM, p_ns_ );

        // Allocation of helper matrices and vectors
        Vec<DIM, Scalar> vel_ts, disp_ts, vel_ns, disp_ns, dummy_vec;
        Mat<DIM, DIM, Scalar> grad_vel_ns, grad_disp_ns, grad_disp_ts, dummy_mat, dummy_mat_T;
        Mat<DIM, DIM, Scalar> F, F_ts, F_inv, F_inv_T;
        Mat<DIM, DIM, Scalar> I;
        for ( int d = 0; d < DIM; ++d )
        {
            I ( d, d ) = 1.;
        }
        assert ( det ( I ) == 1. );
        Mat<DIM, DIM, Scalar> grad_trial_phi, grad_test_phi;

        // loop quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const Scalar wq_dJ = w ( q ) * std::abs ( detJ ( q ) );

            // previous pressure
            const Scalar p_ns = p_ns_[q];

            // Vectors with values of previous newton and time step
            for ( int var = 0; var < DIM; ++var )
            {
                vel_ts[var] = vel_ts_[var][q];
                vel_ns[var] = vel_ns_[var][q];
                disp_ts[var] = disp_ts_[var][q];
                disp_ns[var] = disp_ns_[var][q];
            }

            // Matrices with values of previous newton and time step
            // The derivatives have the format
            // [\nabla F]_ij = dF_i/dx_j
            // Thus, we actually deal with the Jacobian / transposed vector gradient
            for ( int var_c = 0; var_c < DIM; ++var_c )
            { // variable for the component of the vector field
                for ( int var_d = 0; var_d < DIM; ++var_d )
                { // variable for the derivative index
                    grad_vel_ns ( var_c, var_d ) = grad_vel_ns_[var_c][q][var_d];
                    grad_disp_ns ( var_c, var_d ) = grad_disp_ns_[var_c][q][var_d];
                    grad_disp_ts ( var_c, var_d ) = grad_disp_ts_[var_c][q][var_d];
                }
            }

            F = grad_disp_ns + I;
            F_ts = grad_disp_ts + I;

            inv ( F, F_inv );
            trans ( F_inv, F_inv_T );

            const Scalar J = det ( F );
            const Scalar J_ts = det ( F_ts );

            const Scalar alpha = mesh_diff_ / J;
            const Scalar J_theta = theta_ * J + ( 1 - theta_ ) * J_ts;

            // Time derivative terms: Velocity trial variables
            for ( int var_vel = 0; var_vel < DIM; ++var_vel )
            {
                const int n_dofs = num_dofs ( var_vel );
                for ( int i = 0; i < n_dofs; ++i )
                {
                    for ( int j = 0; j < n_dofs; ++j )
                    {
                        lm ( dof_index ( i, var_vel ), dof_index ( j, var_vel ) ) +=
                                wq_dJ * J_theta / dt_ * phi ( j, q, var_vel ) * phi ( i, q, var_vel );
                    }
                }
            }

            // Time derivative terms: Velocity trial variables
            for ( int trial_var_vel = 0; trial_var_vel < DIM; ++trial_var_vel )
            {
                grad_trial_phi *= 0.;
                for ( int j = 0; j < num_dofs ( trial_var_vel ); ++j )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_vel, d ) = grad_phi ( j, q, trial_var_vel )[d];
                    }
                    dummy_vec = grad_trial_phi * ( F_inv * ( disp_ns - disp_ts ) );
                    dummy_vec *= -J / dt_;
                    for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                    {
                        for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                        {
                            lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_vel ) ) +=
                                    wq_dJ * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                        }
                    }
                }
            }

            // Time derivative terms: Displacement trial variables
            for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
            {
                grad_trial_phi *= 0.;
                for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
                {
                    dummy_vec *= 0.;
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                    }
                    dummy_vec[trial_var_disp - 1 - DIM] = phi ( j, q, trial_var_disp );
                    dummy_vec = -J / dt_ * grad_vel_ns * ( F_inv * dummy_vec );
                    dummy_vec += theta_ * J / dt_ * trace ( F_inv * grad_trial_phi ) * ( vel_ns - vel_ts );
                    dummy_vec += -J / dt_ * trace ( F_inv * grad_trial_phi ) * ( grad_vel_ns * ( F_inv * ( disp_ns - disp_ts ) ) );
                    dummy_vec += J / dt_ * ( grad_vel_ns * ( F_inv * grad_trial_phi * F_inv * ( disp_ns - disp_ts ) ) );
                    for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                    {
                        for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                        {
                            lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_disp ) ) +=
                                    wq_dJ * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                        }
                    }
                }
            }

            // Convection terms: Velocity trial variables
            for ( int trial_var_vel = 0; trial_var_vel < DIM; ++trial_var_vel )
            {
                grad_trial_phi *= 0.;
                for ( int j = 0; j < num_dofs ( trial_var_vel ); ++j )
                {
                    dummy_vec *= 0.;
                    dummy_vec[trial_var_vel] = phi ( j, q, trial_var_vel );
                    dummy_vec = theta_ * J * grad_vel_ns * F_inv*dummy_vec;
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_vel, d ) = grad_phi ( j, q, trial_var_vel )[d];
                    }
                    dummy_vec += theta_ * J * ( grad_trial_phi * ( F_inv * vel_ns ) );
                    for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                    {
                        for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                        {
                            lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_vel ) ) +=
                                    wq_dJ * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                        }
                    }
                }
            }

            // Convection terms: Displacement trial variables
            for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
            {
                grad_trial_phi *= 0.;
                for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                    }
                    dummy_vec = theta_ * J * trace ( F_inv * grad_trial_phi ) * grad_vel_ns * ( F_inv * vel_ns );
                    dummy_vec += -theta_ * J * grad_vel_ns * ( F_inv * grad_trial_phi * F_inv * vel_ns );
                    for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                    {
                        for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                        {
                            lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_disp ) ) +=
                                    wq_dJ * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                        }
                    }
                }
            }

            // Fluid stress terms: Velocity trial variables
            for ( int trial_var_vel = 0; trial_var_vel < DIM; ++trial_var_vel )
            {
                for ( int j = 0; j < num_dofs ( trial_var_vel ); ++j )
                {
                    dummy_mat *= 0.;
                    for ( int d = 0; d < DIM; ++d )
                    {
                        dummy_mat ( trial_var_vel, d ) = grad_phi ( j, q, trial_var_vel )[d];
                    }
                    dummy_mat *= F_inv;
                    trans ( dummy_mat, dummy_mat_T );
                    dummy_mat += dummy_mat_T;
                    dummy_mat *= F_inv_T;
                    dummy_mat *= theta_ * J * nu_;
                    for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                    {
                        grad_test_phi *= 0.;
                        for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                        {
                            for ( int d = 0; d < DIM; ++d )
                            {
                                grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                            }
                            lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_vel ) ) +=
                                    wq_dJ * frob ( dummy_mat, grad_test_phi );
                        }
                    }
                }
            }

            // Fluid stress terms: Displacement trial variables
            trans ( grad_vel_ns*F_inv, dummy_mat_T );
            const Mat<DIM, DIM, Scalar> stress ( nu_ * ( grad_vel_ns * F_inv + dummy_mat_T ) * F_inv_T );
            for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
            {
                grad_trial_phi *= 0.;
                for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                    }
                    const Scalar fac_stress = theta_ * J * trace ( F_inv * grad_trial_phi );
                    dummy_mat = grad_vel_ns * F_inv * grad_trial_phi*F_inv;
                    trans ( dummy_mat, dummy_mat_T );
                    dummy_mat += dummy_mat_T;
                    dummy_mat *= -theta_ * J * nu_ * F_inv_T;
                    trans ( grad_trial_phi, dummy_mat_T );
                    dummy_mat += -theta_ * J * stress * dummy_mat_T*F_inv_T;
                    for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                    {
                        grad_test_phi *= 0.;
                        for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                        {
                            for ( int d = 0; d < DIM; ++d )
                            {
                                grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                            }
                            lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_disp ) ) +=
                                    wq_dJ * (
                                    fac_stress * frob ( stress, grad_test_phi )
                                    + frob ( dummy_mat, grad_test_phi )
                                    );
                        }
                    }
                }
            }

            // Mesh regularity
            for ( int var_disp = DIM + 1; var_disp < 2 * DIM + 1; ++var_disp )
            {
                const int n_dofs = num_dofs ( var_disp );
                for ( int i = 0; i < n_dofs; ++i )
                {
                    for ( int j = 0; j < n_dofs; ++j )
                    {
                        lm ( dof_index ( i, var_disp ), dof_index ( j, var_disp ) ) +=
                                wq_dJ *
                                alpha * dot ( grad_phi ( j, q, var_disp ), grad_phi ( i, q, var_disp ) );
                    }
                }
            }

            // Mesh regularity: Linearization of the factor alpha
            for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
            {
                grad_trial_phi *= 0.;
                for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                    }
                    const Scalar fac_reg = -alpha * trace ( F_inv * grad_trial_phi );
                    for ( int test_var_disp = DIM + 1; test_var_disp < 2 * DIM + 1; ++test_var_disp )
                    {
                        grad_test_phi *= 0.;
                        for ( int i = 0; i < num_dofs ( test_var_disp ); ++i )
                        {
                            for ( int d = 0; d < DIM; ++d )
                            {
                                grad_test_phi ( test_var_disp - 1 - DIM, d ) = grad_phi ( i, q, test_var_disp )[d];
                            }
                            lm ( dof_index ( i, test_var_disp ), dof_index ( j, trial_var_disp ) ) +=
                                    wq_dJ *
                                    fac_reg * frob ( grad_disp_ns, grad_test_phi );
                        }
                    }
                }
            }

            // Divergence term: Velocity trial variables
            for ( int trial_var_vel = 0; trial_var_vel < DIM; ++trial_var_vel )
            {
                grad_trial_phi *= 0.;
                for ( int j = 0; j < num_dofs ( trial_var_vel ); ++j )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_vel, d ) = grad_phi ( j, q, trial_var_vel )[d];
                    }
                    const Scalar fac_div_1 = J * trace ( grad_trial_phi * F_inv );
                    const int test_var_press = DIM;
                    for ( int i = 0; i < num_dofs ( test_var_press ); ++i )
                    {
                        lm ( dof_index ( i, test_var_press ), dof_index ( j, trial_var_vel ) ) +=
                                wq_dJ * fac_div_1 * phi ( i, q, test_var_press );
                    }
                }
            }

            // Divergence term: Displacement trial variables
            for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
            {
                grad_trial_phi *= 0.;
                for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                    }
                    Scalar fac_div_2 = J * trace ( grad_vel_ns * F_inv ) * trace ( F_inv * grad_trial_phi );
                    fac_div_2 += -J * trace ( grad_vel_ns * F_inv * grad_trial_phi * F_inv );
                    const int test_var_press = DIM;
                    for ( int i = 0; i < num_dofs ( test_var_press ); ++i )
                    {
                        lm ( dof_index ( i, test_var_press ), dof_index ( j, trial_var_disp ) ) +=
                                wq_dJ * fac_div_2 * phi ( i, q, test_var_press );
                    }
                }
            }

            // Pressure terms: Pressure trial variable
            const int trial_var_press = DIM;
            for ( int j = 0; j < num_dofs ( trial_var_press ); ++j )
            {
                const Scalar fac_press = -J / rho_ * phi ( j, q, trial_var_press );
                for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                {
                    grad_test_phi *= 0.;
                    for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                    {
                        for ( int d = 0; d < DIM; ++d )
                        {
                            grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                        }
                        lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_press ) ) +=
                                wq_dJ * fac_press * frob ( F_inv_T, grad_test_phi );
                    }
                }
            }

            // Pressure terms: Displacement trial variables
            for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
            {
                for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
                {
                    grad_trial_phi *= 0.;
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                    }
                    const Scalar fac_press_2 = -J / rho_ * p_ns * trace ( F_inv * grad_trial_phi );
                    trans ( grad_trial_phi, dummy_mat_T );
                    grad_trial_phi = F_inv_T * dummy_mat_T*F_inv_T;
                    grad_trial_phi *= p_ns * J / rho_;
                    for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                    {
                        grad_test_phi *= 0.;
                        for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                        {
                            for ( int d = 0; d < DIM; ++d )
                            {
                                grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                            }
                            lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_disp ) ) +=
                                    wq_dJ *
                                    ( fac_press_2 * frob ( F_inv_T, grad_test_phi )
                                    + frob ( grad_trial_phi, grad_test_phi )
                                    );
                        }
                    }
                }
            }
        } // end loop q
    }
    else
    { // Elasticity assembly
        // recompute previous solution values
        for ( int d = 0; d < DIM; ++d )
        {
            grad_disp_ns_[d].clear ( );
            evaluate_fe_function_gradients ( *u_, d + 1 + DIM, grad_disp_ns_[d] );
        }

        // Allocation of helper matrices and vectors
        Mat<DIM, DIM, Scalar> I;
        for ( int d = 0; d < DIM; ++d )
        {
            I ( d, d ) = 1.;
        }
        assert ( det ( I ) == 1. );
        Mat<DIM, DIM, Scalar> grad_disp_ns, F;
        Mat<DIM, DIM, Scalar> grad_trial_phi, grad_test_phi, dummy_mat, dummy_mat_T;

        // loop quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const Scalar wq_dJ = w ( q ) * std::abs ( detJ ( q ) );

            // Matrices with values of previous newton and time step
            // The derivatives have the format
            // [\nabla F]_ij = dF_i/dx_j
            // Thus, we actually deal with the Jacobian / transposed vector gradient
            for ( int var_c = 0; var_c < DIM; ++var_c )
            { // variable for the component of the vector field
                for ( int var_d = 0; var_d < DIM; ++var_d )
                { // variable for the derivative index
                    grad_disp_ns ( var_c, var_d ) = grad_disp_ns_[var_c][q][var_d];
                }
            }
            F = grad_disp_ns + I;

            // Velocity time derivative
            const Scalar inv_dt = 1. / dt_;
            for ( int var_vel = 0; var_vel < DIM; ++var_vel )
            {
                const int n_dofs = num_dofs ( var_vel );
                for ( int i = 0; i < n_dofs; ++i )
                {
                    for ( int j = 0; j < n_dofs; ++j )
                    {
                        lm ( dof_index ( i, var_vel ), dof_index ( j, var_vel ) ) +=
                                wq_dJ * inv_dt * phi ( j, q, var_vel ) * phi ( i, q, var_vel );
                    }
                }
            }

            // Displacement time derivative
            for ( int var_disp = DIM + 1; var_disp < 2 * DIM + 1; ++var_disp )
            {
                const int n_dofs = num_dofs ( var_disp );
                for ( int i = 0; i < n_dofs; ++i )
                {
                    for ( int j = 0; j < n_dofs; ++j )
                    {
                        lm ( dof_index ( i, var_disp ), dof_index ( j, var_disp ) ) +=
                                wq_dJ * inv_dt * phi ( j, q, var_disp ) * phi ( i, q, var_disp );
                    }
                }
            }

            // Velocity tested with displacement functions
            for ( int var_vel = 0; var_vel < DIM; ++var_vel )
            {
                int var_disp = var_vel + 1 + DIM;
                for ( int j = 0; j < num_dofs ( var_vel ); ++j )
                {
                    for ( int i = 0; i < num_dofs ( var_disp ); ++i )
                    {
                        lm ( dof_index ( i, var_disp ), dof_index ( j, var_vel ) ) +=
                                -wq_dJ * theta_ * phi ( j, q, var_vel ) * phi ( i, q, var_disp );
                    }
                }
            }

            // Structure stress terms
            if ( linear_ )
            {
                for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
                {
                    grad_trial_phi *= 0.;
                    for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
                    {
                        for ( int d = 0; d < DIM; ++d )
                        {
                            grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                        }
                        trans ( grad_trial_phi, dummy_mat );
                        dummy_mat += grad_trial_phi;
                        dummy_mat = lambda_ * 0.5 * trace ( dummy_mat ) * I + mu_*dummy_mat;
                        dummy_mat *= theta_ / rho_;
                        for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                        {
                            grad_test_phi *= 0.;
                            for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                            {
                                for ( int d = 0; d < DIM; ++d )
                                {
                                    grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                                }
                                lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_disp ) ) +=
                                        wq_dJ * frob ( dummy_mat, grad_test_phi );
                            }
                        }
                    }
                }
            }
            else
            {
                Mat<DIM, DIM, Scalar> stress;
                trans ( F, stress );
                stress *= F;
                stress -= I;
                stress *= 0.5;
                stress = ( lambda_ * trace ( stress ) * I ) + ( 2. * mu_ * stress );
                for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
                {
                    grad_trial_phi *= 0.;
                    for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
                    {
                        for ( int d = 0; d < DIM; ++d )
                        {
                            grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                        }
                        trans ( grad_trial_phi, dummy_mat );
                        dummy_mat *= F;
                        trans ( dummy_mat, dummy_mat_T );
                        dummy_mat += dummy_mat_T;
                        const Scalar tr = 0.5 * lambda_ * trace ( dummy_mat );
                        dummy_mat *= mu_;
                        for ( int d = 0; d < DIM; ++d )
                        {
                            dummy_mat ( d, d ) = dummy_mat ( d, d ) + tr;
                        }
                        dummy_mat = theta_ / rho_ * F*dummy_mat;
                        dummy_mat += theta_ / rho_ * grad_trial_phi * stress;
                        for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                        {
                            grad_test_phi *= 0.;
                            for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                            {
                                for ( int d = 0; d < DIM; ++d )
                                {
                                    grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                                }
                                lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_disp ) ) +=
                                        wq_dJ * frob ( dummy_mat, grad_test_phi );
                            }
                        }
                    }
                }
            }
        } // end loop q
    } // end if material
}

void FSIBenchmarkAssembler::operator() ( const Element<Scalar>& element, const Quadrature<Scalar>& quadrature, LocalVector& lv )
{
    initialize_for_element ( element, quadrature );

    const int mat = element.get_cell ( ).get_material_number ( );
    assert ( ( fluid_mat_ == mat ) || ( solid_mat_ == mat ) );

    // indices j -> trial variable, i -> test variable
    // basis functions \phi -> velocity components, \eta -> pressure, \psi -> displacement

    const int num_q = num_quadrature_points ( );

    if ( fluid_mat_ == mat )
    {
        // recompute previous solution values
        for ( int d = 0; d < DIM; ++d )
        {
            vel_ns_[d].clear ( );
            vel_ts_[d].clear ( );
            grad_vel_ns_[d].clear ( );
            grad_vel_ts_[d].clear ( );

            evaluate_fe_function ( *u_, d, vel_ns_[d] );
            evaluate_fe_function ( *prev_sol_, d, vel_ts_[d] );
            evaluate_fe_function_gradients ( *u_, d, grad_vel_ns_[d] );
            evaluate_fe_function_gradients ( *prev_sol_, d, grad_vel_ts_[d] );

            disp_ns_[d].clear ( );
            disp_ts_[d].clear ( );
            grad_disp_ns_[d].clear ( );
            grad_disp_ts_[d].clear ( );

            evaluate_fe_function ( *u_, d + 1 + DIM, disp_ns_[d] );
            evaluate_fe_function ( *prev_sol_, d + 1 + DIM, disp_ts_[d] );
            evaluate_fe_function_gradients ( *u_, d + 1 + DIM, grad_disp_ns_[d] );
            evaluate_fe_function_gradients ( *prev_sol_, d + 1 + DIM, grad_disp_ts_[d] );
        }

        // recompute pressure
        p_ns_.clear ( );
        evaluate_fe_function ( *u_, DIM, p_ns_ );

        // Allocation of helper matrices and vectors
        Vec<DIM, Scalar> vel_ts, disp_ts, vel_ns, disp_ns;
        Mat<DIM, DIM, Scalar> grad_vel_ns, grad_vel_ts, grad_disp_ns, grad_disp_ts;
        Mat<DIM, DIM, Scalar> I;
        for ( int d = 0; d < DIM; ++d )
        {
            I ( d, d ) = 1.;
        }
        assert ( det ( I ) == 1. );
        Mat<DIM, DIM, Scalar> F, F_ts, F_inv, F_inv_T, F_inv_ts, F_inv_T_ts;

        Vec<DIM, Scalar> dummy_vec;
        Mat<DIM, DIM, Scalar> grad_test_phi, dummy_mat, dummy_mat_T;

        // loop quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const Scalar wq_dJ = w ( q ) * std::abs ( detJ ( q ) );

            // previous pressure
            const Scalar p_ns = p_ns_[q];

            // Vectors with values of previous newton and time step
            for ( int var = 0; var < DIM; ++var )
            {
                vel_ts[var] = vel_ts_[var][q];
                vel_ns[var] = vel_ns_[var][q];
                disp_ts[var] = disp_ts_[var][q];
                disp_ns[var] = disp_ns_[var][q];
            }

            // Matrices with values of previous newton and time step
            // The derivatives have the format
            // [\nabla F]_ij = dF_i/dx_j
            // Thus, we actually deal with the Jacobian / transposed vector gradient
            for ( int var_c = 0; var_c < DIM; ++var_c )
            { // variable for the component of the vector field
                for ( int var_d = 0; var_d < DIM; ++var_d )
                { // variable for the derivative index
                    grad_vel_ns ( var_c, var_d ) = grad_vel_ns_[var_c][q][var_d];
                    grad_vel_ts ( var_c, var_d ) = grad_vel_ts_[var_c][q][var_d];
                    grad_disp_ns ( var_c, var_d ) = grad_disp_ns_[var_c][q][var_d];
                    grad_disp_ts ( var_c, var_d ) = grad_disp_ts_[var_c][q][var_d];
                }
            }

            F = grad_disp_ns + I;
            F_ts = grad_disp_ts + I;

            inv ( F, F_inv );
            trans ( F_inv, F_inv_T );
            inv ( F_ts, F_inv_ts );
            trans ( F_inv_ts, F_inv_T_ts );

            const Scalar J = det ( F );
            const Scalar J_ts = det ( F_ts );

            const Scalar alpha = mesh_diff_ / J;
            const Scalar J_theta = theta_ * J + ( 1 - theta_ ) * J_ts;

            // Time derivative terms
            dummy_vec = J_theta / dt_ * ( vel_ns - vel_ts );
            dummy_vec += -J / dt_ * ( grad_vel_ns * ( F_inv * ( disp_ns - disp_ts ) ) );
            for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
            {
                for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                {
                    lv[dof_index ( i, test_var_vel )] +=
                            wq_dJ * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                }
            }

            // Convection terms
            dummy_vec = theta_ * J * ( grad_vel_ns * ( F_inv * vel_ns ) );
            dummy_vec += ( 1 - theta_ ) * J_ts * ( grad_vel_ts * ( F_inv_ts * vel_ts ) );
            for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
            {
                for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                {
                    lv[dof_index ( i, test_var_vel )] +=
                            wq_dJ * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                }
            }

            // Fluid stress terms
            trans ( grad_vel_ns*F_inv, dummy_mat_T );
            dummy_mat = theta_ * nu_ * J * ( ( grad_vel_ns * F_inv + dummy_mat_T ) * F_inv_T );
            trans ( grad_vel_ts*F_inv_ts, dummy_mat_T );
            dummy_mat += ( 1 - theta_ ) * nu_ * J_ts * ( ( grad_vel_ts * F_inv_ts + dummy_mat_T ) * F_inv_T_ts );
            for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
            {
                grad_test_phi *= 0.;
                for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                    }
                    lv[dof_index ( i, test_var_vel )] +=
                            wq_dJ * frob ( dummy_mat, grad_test_phi );
                }
            }

            // Mesh regularity
            for ( int test_var_disp = DIM + 1; test_var_disp < 2 * DIM + 1; ++test_var_disp )
            {
                grad_test_phi *= 0.;
                for ( int i = 0; i < num_dofs ( test_var_disp ); ++i )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_test_phi ( test_var_disp - 1 - DIM, d ) = grad_phi ( i, q, test_var_disp )[d];
                    }
                    lv[dof_index ( i, test_var_disp )] +=
                            wq_dJ *
                            alpha * frob ( grad_disp_ns, grad_test_phi );
                }
            }

            // Divergence term
            const Scalar div = J * trace ( grad_vel_ns * F_inv );
            const int test_var_press = DIM;
            for ( int i = 0; i < num_dofs ( test_var_press ); ++i )
            {
                lv[dof_index ( i, test_var_press )] +=
                        wq_dJ *
                        div * phi ( i, q, test_var_press );
            }

            // Pressure term
            const Scalar press_ns = -p_ns * J / rho_;
            for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
            {
                grad_test_phi *= 0.;
                for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                    }
                    lv[dof_index ( i, test_var_vel )] +=
                            wq_dJ *
                            press_ns * frob ( F_inv_T, grad_test_phi );
                }
            }
        } // end quadrature loop

    }
    else
    { // Elasticity assembly

        // recompute previous solution values
        for ( int d = 0; d < DIM; ++d )
        {
            vel_ns_[d].clear ( );
            vel_ts_[d].clear ( );

            evaluate_fe_function ( *u_, d, vel_ns_[d] );
            evaluate_fe_function ( *prev_sol_, d, vel_ts_[d] );

            disp_ns_[d].clear ( );
            disp_ts_[d].clear ( );
            grad_disp_ns_[d].clear ( );
            grad_disp_ts_[d].clear ( );

            evaluate_fe_function ( *u_, d + 1 + DIM, disp_ns_[d] );
            evaluate_fe_function ( *prev_sol_, d + 1 + DIM, disp_ts_[d] );
            evaluate_fe_function_gradients ( *u_, d + 1 + DIM, grad_disp_ns_[d] );
            evaluate_fe_function_gradients ( *prev_sol_, d + 1 + DIM, grad_disp_ts_[d] );
        }

        // Allocation of helper matrices and vectors
        Vec<DIM, Scalar> vel_ts, disp_ts, vel_ns, disp_ns;
        Mat<DIM, DIM, Scalar> grad_disp_ns, grad_disp_ts;
        Mat<DIM, DIM, Scalar> I;
        for ( int d = 0; d < DIM; ++d )
        {
            I ( d, d ) = 1.;
        }
        assert ( det ( I ) == 1. );

        Vec<DIM, Scalar> dummy_vec;
        Mat<DIM, DIM, Scalar> grad_test_phi, dummy_mat, dummy_mat_2;

        // loop quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const Scalar wq_dJ = w ( q ) * std::abs ( detJ ( q ) );

            // Vectors with values of previous newton and time step
            for ( int var = 0; var < DIM; ++var )
            {
                vel_ts[var] = vel_ts_[var][q];
                vel_ns[var] = vel_ns_[var][q];
                disp_ts[var] = disp_ts_[var][q];
                disp_ns[var] = disp_ns_[var][q];
            }

            // Matrices with values of previous newton and time step
            // The derivatives have the format
            // [\nabla F]_ij = dF_i/dx_j
            // Thus, we actually deal with the Jacobian / transposed vector gradient
            for ( int var_c = 0; var_c < DIM; ++var_c )
            { // variable for the component of the vector field
                for ( int var_d = 0; var_d < DIM; ++var_d )
                { // variable for the derivative index
                    grad_disp_ns ( var_c, var_d ) = grad_disp_ns_[var_c][q][var_d];
                    grad_disp_ts ( var_c, var_d ) = grad_disp_ts_[var_c][q][var_d];
                }
            }

            // Velocity time derivative
            dummy_vec = 1. / dt_ * ( vel_ns - vel_ts );
            for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
            {
                for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                {
                    lv[dof_index ( i, test_var_vel )] +=
                            wq_dJ * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                }
            }

            // Displacement time derivative
            dummy_vec = 1. / dt_ * ( disp_ns - disp_ts );
            dummy_vec += ( -theta_ * vel_ns + ( theta_ - 1 ) * vel_ts );
            for ( int test_var_disp = DIM + 1; test_var_disp < 2 * DIM + 1; ++test_var_disp )
            {
                for ( int i = 0; i < num_dofs ( test_var_disp ); ++i )
                {
                    lv[dof_index ( i, test_var_disp )] +=
                            wq_dJ * dummy_vec[test_var_disp - 1 - DIM] * phi ( i, q, test_var_disp );
                }
            }

            // Structure stress terms
            if ( linear_ )
            {
                trans ( grad_disp_ns, dummy_mat );
                dummy_mat = 0.5 * ( dummy_mat + grad_disp_ns );
                dummy_mat = theta_ / rho_ * ( lambda_ * trace ( dummy_mat ) * I + 2. * mu_ * dummy_mat );
                trans ( grad_disp_ts, dummy_mat_2 );
                dummy_mat_2 = 0.5 * ( dummy_mat_2 + grad_disp_ts );
                dummy_mat += ( ( 1 - theta_ ) / rho_ * ( lambda_ * trace ( dummy_mat_2 ) * I + 2. * mu_ * dummy_mat_2 ) );
            }
            else
            {
                trans ( grad_disp_ns + I, dummy_mat );
                trans ( grad_disp_ts + I, dummy_mat_2 );
                dummy_mat = 0.5 * ( dummy_mat * ( grad_disp_ns + I ) - I );
                dummy_mat = ( theta_ / rho_ * ( grad_disp_ns + I )*( lambda_ * trace ( dummy_mat ) * I + 2. * mu_ * dummy_mat ) );
                dummy_mat_2 = 0.5 * ( dummy_mat_2 * ( grad_disp_ts + I ) - I );
                dummy_mat += ( ( 1 - theta_ ) / rho_ * ( grad_disp_ts + I )*( lambda_ * trace ( dummy_mat_2 ) * I + 2. * mu_ * dummy_mat_2 ) );
            }
            for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
            {
                grad_test_phi *= 0.;
                for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_test_phi ( test_var_vel, d ) = grad_phi ( i, q, test_var_vel )[d];
                    }
                    lv[dof_index ( i, test_var_vel )] +=
                            wq_dJ * frob ( dummy_mat, grad_test_phi );
                }
            }
        } // end quadrature loop
    } // end material
}

void FSIBenchmarkAssembler::operator() ( const Element<Scalar>& element,
        const int facet_number,
        const Quadrature<Scalar>& quadrature,
        LocalMatrix& lm )
{
    // Procedure to get the facet entity
    mesh::IncidentEntityIterator facet = element.get_cell ( ).begin_incident ( DIM - 1 );
    for ( int i = 0; i < facet_number; ++i, ++facet )
    {
    }

    // Check if it belongs to the boundary of the fictitious domain
    if ( facet->get_material_number ( ) != outflow_mat_ ) return;

    // Initialize the quadrature for integration over the facet
    initialize_for_facet ( element, quadrature, facet_number );

    // recompute previous solution values
    for ( int d = 0; d < DIM; ++d )
    {
        grad_vel_ns_[d].clear ( );
        evaluate_fe_function_gradients ( *u_, d, grad_vel_ns_[d] );

        grad_disp_ns_[d].clear ( );
        evaluate_fe_function_gradients ( *u_, d + 1 + DIM, grad_disp_ns_[d] );
    }

    // indices j -> trial variable, i -> test variable
    // basis functions \phi -> velocity components, \eta -> pressure, \psi -> displacement

    const int num_q = num_quadrature_points ( );

    // Allocation of helper matrices and vectors
    Mat<DIM, DIM, Scalar> grad_vel_ns, grad_disp_ns;
    Mat<DIM, DIM, Scalar> I;
    for ( int d = 0; d < DIM; ++d )
    {
        I ( d, d ) = 1.;
    }
    assert ( det ( I ) == 1. );
    Mat<DIM, DIM, Scalar> F, F_inv, F_inv_T;

    Vec<DIM, Scalar> dummy_vec;
    Mat<DIM, DIM, Scalar> grad_vel_ns_T, grad_trial_phi, grad_trial_phi_T;

    // loop quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const Scalar wq_dsurf = w ( q ) * ds ( q );
        const Vec<DIM, Scalar> normal = n ( q );

        // Matrices with values of previous newton and time step
        // The derivatives have the format
        // [\nabla F]_ij = dF_i/dx_j
        // Thus, we actually deal with the Jacobian / transposed vector gradient
        for ( int var_c = 0; var_c < DIM; ++var_c )
        { // variable for the component of the vector field
            for ( int var_d = 0; var_d < DIM; ++var_d )
            { // variable for the derivative index
                grad_vel_ns ( var_c, var_d ) = grad_vel_ns_[var_c][q][var_d];
                grad_disp_ns ( var_c, var_d ) = grad_disp_ns_[var_c][q][var_d];
            }
        }

        F = grad_disp_ns + I;
        inv ( F, F_inv );
        trans ( F_inv, F_inv_T );

        const Scalar J = det ( F );

        // Velocity correction terms on outflow boundary
        Mat<DIM, DIM, Scalar> g;
        trans ( grad_vel_ns, grad_vel_ns_T );
        g = nu_ * ( F_inv_T * grad_vel_ns_T );

        // Velocity correction: Displacement trial variables
        for ( int trial_var_disp = DIM + 1; trial_var_disp < 2 * DIM + 1; ++trial_var_disp )
        {
            grad_trial_phi *= 0.;
            for ( int j = 0; j < num_dofs ( trial_var_disp ); ++j )
            {
                for ( int d = 0; d < DIM; ++d )
                {
                    grad_trial_phi ( trial_var_disp - 1 - DIM, d ) = grad_phi ( j, q, trial_var_disp )[d];
                }
                trans ( grad_trial_phi, grad_trial_phi_T );
                dummy_vec = -theta_ * J * trace ( F_inv * grad_trial_phi ) * ( g * F_inv_T * normal );
                dummy_vec += theta_ * J * nu_ * ( F_inv_T * grad_trial_phi_T * F_inv_T * grad_vel_ns_T * F_inv_T * normal );
                dummy_vec += theta_ * J * ( g * F_inv_T * grad_trial_phi_T * F_inv_T * normal );
                for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                {
                    for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                    {
                        lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_disp ) ) +=
                                wq_dsurf * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                    }
                }
            }
        }

        // Velocity correction terms on outflow boundary: Velocity trial variables
        for ( int trial_var_vel = 0; trial_var_vel < DIM; ++trial_var_vel )
        {
            grad_trial_phi *= 0.;
            for ( int j = 0; j < num_dofs ( trial_var_vel ); ++j )
            {
                for ( int d = 0; d < DIM; ++d )
                {
                    grad_trial_phi ( trial_var_vel, d ) = grad_phi ( j, q, trial_var_vel )[d];
                }
                trans ( grad_trial_phi, grad_trial_phi_T );
                dummy_vec = -theta_ * J * nu_ * ( ( F_inv_T * grad_trial_phi_T * F_inv_T ) * normal );
                for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
                {
                    for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
                    {
                        lm ( dof_index ( i, test_var_vel ), dof_index ( j, trial_var_vel ) ) +=
                                wq_dsurf * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
                    }
                }
            }
        }
    }
}

void FSIBenchmarkAssembler::operator() ( const Element<Scalar>& element,
        const int facet_number,
        const Quadrature<Scalar>& quadrature,
        LocalVector& lv )
{
    // Procedure to get the facet entity
    mesh::IncidentEntityIterator facet = element.get_cell ( ).begin_incident ( DIM - 1 );
    for ( int i = 0; i < facet_number; ++i, ++facet )
    {
    }

    // Check if it belongs to the boundary of the fictitious domain
    if ( facet->get_material_number ( ) != outflow_mat_ ) return;

    // Initialize the quadrature for integration over the facet
    initialize_for_facet ( element, quadrature, facet_number );

    // recompute previous solution values
    for ( int d = 0; d < DIM; ++d )
    {
        grad_vel_ns_[d].clear ( );
        grad_vel_ts_[d].clear ( );
        evaluate_fe_function_gradients ( *u_, d, grad_vel_ns_[d] );
        evaluate_fe_function_gradients ( *prev_sol_, d, grad_vel_ts_[d] );

        grad_disp_ns_[d].clear ( );
        grad_disp_ts_[d].clear ( );
        evaluate_fe_function_gradients ( *u_, d + 1 + DIM, grad_disp_ns_[d] );
        evaluate_fe_function_gradients ( *prev_sol_, d + 1 + DIM, grad_disp_ts_[d] );
    }

    // indices j -> trial variable, i -> test variable
    // basis functions \phi -> velocity components, \eta -> pressure, \psi -> displacement

    // Allocation of helper matrices and vectors
    Mat<DIM, DIM, Scalar> grad_vel_ns, grad_vel_ts, grad_disp_ns, grad_disp_ts;
    Mat<DIM, DIM, Scalar> F, F_inv, F_inv_T, F_ts, F_inv_ts, F_inv_T_ts;

    Mat<DIM, DIM, Scalar> I;
    for ( int d = 0; d < DIM; ++d )
    {
        I ( d, d ) = 1.;
    }
    assert ( det ( I ) == 1. );

    Vec<DIM, Scalar> dummy_vec;
    Mat<DIM, DIM, Scalar> dummy_mat;

    const int num_q = num_quadrature_points ( );

    // loop quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const Scalar wq_dsurf = w ( q ) * ds ( q );
        const Vec<DIM, Scalar> normal = n ( q );

        // Matrices with values of previous newton and time step
        // The derivatives have the format
        // [\nabla F]_ij = dF_i/dx_j
        // Thus, we actually deal with the Jacobian / transposed vector gradient
        for ( int var_c = 0; var_c < DIM; ++var_c )
        { // variable for the component of the vector field
            for ( int var_d = 0; var_d < DIM; ++var_d )
            { // variable for the derivative index
                grad_vel_ns ( var_c, var_d ) = grad_vel_ns_[var_c][q][var_d];
                grad_vel_ts ( var_c, var_d ) = grad_vel_ts_[var_c][q][var_d];
                grad_disp_ns ( var_c, var_d ) = grad_disp_ns_[var_c][q][var_d];
                grad_disp_ts ( var_c, var_d ) = grad_disp_ts_[var_c][q][var_d];
            }
        }

        F = grad_disp_ns + I;
        F_ts = grad_disp_ts + I;

        inv ( F, F_inv );
        trans ( F_inv, F_inv_T );
        inv ( F_ts, F_inv_ts );
        trans ( F_inv_ts, F_inv_T_ts );

        const Scalar J = det ( F );
        const Scalar J_ts = det ( F_ts );

        // Velocity correction terms on outflow boundary
        trans ( grad_vel_ns, dummy_mat );
        dummy_vec = -theta_ * J * nu_ * ( ( F_inv_T * dummy_mat * F_inv_T ) * normal );
        trans ( grad_vel_ts, dummy_mat );
        dummy_vec += ( 1 - theta_ ) * J_ts * nu_ * ( ( F_inv_T_ts * dummy_mat * F_inv_T_ts ) * normal );
        for ( int test_var_vel = 0; test_var_vel < DIM; ++test_var_vel )
        {
            for ( int i = 0; i < num_dofs ( test_var_vel ); ++i )
            {
                lv[dof_index ( i, test_var_vel )] +=
                        wq_dsurf * dummy_vec[test_var_vel] * phi ( i, q, test_var_vel );
            }
        }
    }
}

void FSIForceIntegrator::operator() ( const Element<Scalar>& left_elem, const Element<Scalar>& right_elem,
        const Quadrature<Scalar>& left_quad, const Quadrature<Scalar>& right_quad,
        int left_facet_number, int right_facet_number,
        InterfaceSide left_if_side, InterfaceSide right_if_side,
        Scalar& force )
{

    if ( left_elem.get_cell ( ).get_material_number ( ) != fluid_material_ ) return;
    // check the material number of the facet
    mesh::IncidentEntityIterator fac_it = left_elem.get_cell ( ).begin_incident ( DIM - 1 );
    for ( int i = 0; i < left_facet_number; ++i, ++fac_it )
    {
    }
    if ( fac_it->get_material_number ( ) != mat_num_ ) return;

    initialize_for_interface ( left_elem, right_elem,
                               left_quad, right_quad,
                               left_facet_number, right_facet_number,
                               left_if_side, right_if_side );

    pressure_.clear ( );
    trial ( ).evaluate_fe_function ( sol_, DIM, pressure_ );

    for ( int d = 0; d < DIM; ++d )
    {
        grad_vel_[d].clear ( );
        grad_disp_[d].clear ( );
        trial ( ).evaluate_fe_function_gradients ( sol_, d, grad_vel_[d] );
        trial ( ).evaluate_fe_function_gradients ( sol_, d + DIM + 1, grad_disp_[d] );
    }

    // Allocation of helper matrices and vectors
    Mat<DIM, DIM, Scalar> grad_vel, grad_disp;

    Mat<DIM, DIM, Scalar> I;
    for ( int d = 0; d < DIM; ++d )
    {
        I ( d, d ) = 1.;
    }
    assert ( det ( I ) == 1. );

    Mat<DIM, DIM, Scalar> F, F_inv, F_inv_T, stress;

    for ( int q = 0; q < num_quadrature_points ( ); ++q )
    {
        const Vec<DIM, Scalar> normal = trial ( ).n ( q );
        const Scalar w_ds_q = w ( q ) * ds ( q );

        for ( int var_c = 0; var_c < DIM; ++var_c )
        { // variable for the component of the vector field
            for ( int var_d = 0; var_d < DIM; ++var_d )
            { // variable for the derivative index
                grad_vel ( var_c, var_d ) = grad_vel_[var_c][q][var_d];
                grad_disp ( var_c, var_d ) = grad_disp_[var_c][q][var_d];
            }
        }

        F = grad_disp + I;
        inv ( F, F_inv );
        trans ( F_inv, F_inv_T );

        const Scalar J = det ( F );

        trans ( grad_vel*F_inv, stress );
        stress = -dyn_visc_ * ( grad_vel * F_inv + stress );

        stress += pressure_[q] * I;
        stress = J * stress*F_inv_T;

        // multiply with normal vector
        for ( int d = 0; d < DIM; ++d )
        {
            force += w_ds_q * stress ( var_, d ) * normal[d];
        }
    }
}

void FSIForceMagnIntegrator::operator() ( const Element<Scalar>& left_elem, const Element<Scalar>& right_elem,
        const Quadrature<Scalar>& left_quad, const Quadrature<Scalar>& right_quad,
        int left_facet_number, int right_facet_number,
        InterfaceSide left_if_side, InterfaceSide right_if_side,
        Scalar& force )
{

    if ( left_elem.get_cell ( ).get_material_number ( ) != fluid_material_ ) return;
    // check the material number of the facet
    mesh::IncidentEntityIterator fac_it = left_elem.get_cell ( ).begin_incident ( DIM - 1 );
    for ( int i = 0; i < left_facet_number; ++i, ++fac_it )
    {
    }
    if ( fac_it->get_material_number ( ) != mat_num_ ) return;

    initialize_for_interface ( left_elem, right_elem,
                               left_quad, right_quad,
                               left_facet_number, right_facet_number,
                               left_if_side, right_if_side );

    for ( int d = 0; d < DIM; ++d )
    {
        grad_vel_[d].clear ( );
        grad_disp_[d].clear ( );
        trial ( ).evaluate_fe_function_gradients ( sol_, d, grad_vel_[d] );
        trial ( ).evaluate_fe_function_gradients ( sol_, d + DIM + 1, grad_disp_[d] );
    }

    // Allocation of helper matrices and vectors
    Mat<DIM, DIM, Scalar> grad_vel, grad_disp;

    Mat<DIM, DIM, Scalar> I;
    for ( int d = 0; d < DIM; ++d )
    {
        I ( d, d ) = 1.;
    }
    assert ( det ( I ) == 1. );

    Mat<DIM, DIM, Scalar> F, F_inv, F_inv_T, stress;

    for ( int q = 0; q < num_quadrature_points ( ); ++q )
    {
        const Vec<DIM, Scalar> normal = trial ( ).n ( q );
        const Scalar w_ds_q = w ( q ) * ds ( q );

        for ( int var_c = 0; var_c < DIM; ++var_c )
        { // variable for the component of the vector field
            for ( int var_d = 0; var_d < DIM; ++var_d )
            { // variable for the derivative index
                grad_vel ( var_c, var_d ) = grad_vel_[var_c][q][var_d];
                grad_disp ( var_c, var_d ) = grad_disp_[var_c][q][var_d];
            }
        }

        F = grad_disp + I;
        inv ( F, F_inv );
        trans ( F_inv, F_inv_T );

        const Scalar J = det ( F );

        trans ( grad_vel*F_inv, stress );
        stress = -dyn_visc_ * ( grad_vel * F_inv + stress );
        stress = J * stress*F_inv_T;

        // multiply with normal vector
        force += w_ds_q * norm ( stress * normal );
    }
}

Scalar flowL2error::exact_vel ( int var, const Vec<DIM, Scalar>& x )
{
    const Scalar radius = norm ( x );

    Scalar vel = ( ( var == 1 ) ? -1. : 1. ) * x[( var + 1 ) % 2] / radius;
    vel *= vel0_ * radii_[0] / ( radii_[0] * radii_[0] - radii_[1] * radii_[1] )*( radii_[1] * radii_[1] / radius - radius );

    return vel;
}

void flowL2error::operator() ( const Element<Scalar>& element,
        const Quadrature<Scalar>& quadrature, Scalar& error )
{

    if ( element.get_cell ( ).get_material_number ( ) != fluid_material_ ) return;

    initialize_for_element ( element, quadrature );

    // recompute solution values
    for ( int d = 0; d < DIM; ++d )
    {
        vel_[d].clear ( );
        disp_[d].clear ( );
        evaluate_fe_function ( sol_, d, vel_[d] );
        evaluate_fe_function ( sol_, DIM + 1 + d, disp_[d] );
    }

    const int num_q = num_quadrature_points ( );

    // loop quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const Scalar wq_dJ = w ( q ) * std::abs ( detJ ( q ) );
        Vec<DIM, Scalar> xq = x ( q );
        for ( int d = 0; d < DIM; ++d )
        {
            xq[d] += disp_[d][q];
        }

        for ( int var = 0; var < DIM; ++var )
        {
            Scalar diff = vel_[var][q] - exact_vel ( var, xq );
            error += wq_dJ * diff * diff;
        }
    }
}

Scalar pressureL2error::exact_press ( const Vec<DIM, Scalar>& x )
{
    Scalar press = 0;
    const Scalar radius = norm ( x );

    press = vel0_ * radii_[0] / ( radii_[0] * radii_[0] - radii_[1] * radii_[1] );
    press *= rho_ * press;
    press *= 0.5 * radius * radius
            - 2. * radii_[1] * radii_[1] * log ( radius )
            - 0.5 * radii_[1] * radii_[1] * radii_[1] * radii_[1] / ( radius * radius )
            + 2. * radii_[1] * radii_[1] * log ( radii_[1] );

    return press;
}

void pressureL2error::operator() ( const Element<Scalar>& element,
        const Quadrature<Scalar>& quadrature, Scalar& error )
{

    if ( element.get_cell ( ).get_material_number ( ) != fluid_material_ ) return;

    initialize_for_element ( element, quadrature );

    // recompute solution values
    press_.clear ( );
    evaluate_fe_function ( sol_, DIM, press_ );
    for ( int d = 0; d < DIM; ++d )
    {
        disp_[d].clear ( );
        evaluate_fe_function ( sol_, DIM + 1 + d, disp_[d] );
    }

    const int num_q = num_quadrature_points ( );

    // loop quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const Scalar wq_dJ = w ( q ) * std::abs ( detJ ( q ) );
        Vec<DIM, Scalar> xq = x ( q );
        for ( int d = 0; d < DIM; ++d )
        {
            xq[d] += disp_[d][q];
        }

        Scalar diff = press_[q] - exact_press ( xq );
        error += wq_dJ * diff * diff;
    }
}

Scalar dispL2error::exact_disp ( const int var, const Vec<DIM, Scalar>& x )
{
    const Scalar radius = norm ( x );

    Scalar disp = -disp1_ * radii_[1] * ( radius - radii_[2] * radii_[2] / radius ) / ( radii_[1] * radii_[1] - radii_[2] * radii_[2] );

    // compute angle coordinate
    Scalar phi = std::atan2 ( x[1], x[0] );
    // add the radial displacement
    phi += disp / radius;
    if ( var == 0 )
    { // x-component
        disp = radius * std::cos ( phi ) - x[0];
    }
    else if ( var == 1 )
    { // y-component
        disp = radius * std::sin ( phi ) - x[1];
    }
    else
    {
        assert ( false );
    }

    return disp;
}

void dispL2error::operator() ( const Element<Scalar>& element,
        const Quadrature<Scalar>& quadrature, Scalar& error )
{

    if ( element.get_cell ( ).get_material_number ( ) != solid_material_ ) return;

    initialize_for_element ( element, quadrature );

    // recompute solution values
    for ( int d = 0; d < DIM; ++d )
    {
        disp_[d].clear ( );
        evaluate_fe_function ( sol_, DIM + 1 + d, disp_[d] );
    }

    const int num_q = num_quadrature_points ( );

    // loop quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const Scalar wq_dJ = w ( q ) * std::abs ( detJ ( q ) );
        Vec<DIM, Scalar> xq = x ( q );

        for ( int var = 0; var < DIM; ++var )
        {
            Scalar diff = disp_[var][q] - exact_disp ( var, xq );
            error += wq_dJ * diff * diff;
        }
    }
}

Scalar diffusionL2error::exact_diff ( const int var, const Vec<DIM, Scalar>& x )
{
    const Scalar radius = norm ( x );

    Scalar disp = ( ( var == 1 ) ? -1. : 1. ) * disp1_ * x[( var + 1 ) % 2] / radius;
    disp /= log ( radii_[1] ) - log ( radii_[0] );
    disp *= log ( radius ) - log ( radii_[0] );

    return disp;
}

void diffusionL2error::operator() ( const Element<Scalar>& element,
        const Quadrature<Scalar>& quadrature, Scalar& error )
{

    if ( element.get_cell ( ).get_material_number ( ) != fluid_material_ ) return;

    initialize_for_element ( element, quadrature );

    // recompute solution values
    for ( int d = 0; d < DIM; ++d )
    {
        disp_[d].clear ( );
        evaluate_fe_function ( sol_, DIM + 1 + d, disp_[d] );
    }

    const int num_q = num_quadrature_points ( );

    // loop quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const Scalar wq_dJ = w ( q ) * std::abs ( detJ ( q ) );
        Vec<DIM, Scalar> xq = x ( q );

        for ( int var = 0; var < DIM; ++var )
        {
            Scalar diff = disp_[var][q] - exact_diff ( var, xq );
            error += wq_dJ * diff * diff;
        }
    }
}
