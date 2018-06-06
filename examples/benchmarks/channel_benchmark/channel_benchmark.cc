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

#include "channel_benchmark.h"

#include <iomanip>
#include <iterator>
#include <algorithm>
#include <vector>

namespace
{
    static const char* DATADIR = MESHES_DATADIR;
    static const int MASTER_RANK = 0;
    static const char* PARAM_FILENAME = "dfg_bench_param.xml";

    static bool CONSOLE_OUTPUT_ACTIVE = true;
    static const int CONSOLE_THRESHOLD_LEVEL = 3;
}
#define CONSOLE_OUTPUT(lvl, x) { if (CONSOLE_OUTPUT_ACTIVE &&            \
                                    lvl <= CONSOLE_THRESHOLD_LEVEL) {   \
            for (int i = 0; i < lvl; ++i) { std::cout << "  "; }        \
            std::cout << x << "\n"; }}

//#if !(DIMENSION == 3)
//#error "The channel benchmark only works in 3d!"
//#endif

struct TimingData
{
    double time_elapsed;
};

class TimingScope
{
  public:

    TimingScope ( const std::string& name )
    {
        if ( report_ )
        {
            report_->begin_section ( name );
        }
    }

    TimingScope ( int iteration )
    {
        if ( report_ )
        {
            std::stringstream sstr;
            sstr << "Iteration " << iteration;
            report_->begin_section ( sstr.str ( ) );
            timer_.reset ( );
            timer_.start ( );
        }
    }

    ~TimingScope ( )
    {
        timer_.stop ( );
        if ( report_ )
        {
            TimingData* data = report_->end_section ( );
            data->time_elapsed = timer_.get_duration ( );
        }
    }

    static void set_report ( HierarchicalReport<TimingData>* report )
    {
        report_ = report;
    }

  private:
    static HierarchicalReport<TimingData>* report_;
    Timer timer_;
};

HierarchicalReport<TimingData>* TimingScope::report_ = 0;

class TimingReportOutputVisitor
{
  public:

    TimingReportOutputVisitor ( std::ostream& os )
    : os_ ( os ), level_ ( 0 )
    {
    }

    void enter ( const std::string& name, TimingData* data )
    {
        if ( name == "root" )
        {
            os_ << "+++ Timing Report +++\n\n";
        }
        else
        {
            for ( int l = 0; l < level_; ++l )
            {
                os_ << "  ";
            }
            os_ << name << " took " << data->time_elapsed << " s.\n";
            ++level_;
        }
    }

    void exit ( const std::string& name, TimingData* data )
    {
        if ( name == "root" )
        {
            os_ << "\n+++ End Timing Report +++\n\n";
        }
        else
        {
            --level_;
        }
    }

  private:
    std::ostream& os_;
    int level_;
};

class ChannelBenchmark : public NonlinearProblem<LAD>
{
  public:

    ChannelBenchmark ( const std::string& param_filename )
    : comm_ ( MPI_COMM_WORLD ),
    params_ ( param_filename.c_str ( ), MASTER_RANK, MPI_COMM_WORLD ),
    use_pressure_filter_ ( false ),
    refinement_level_ ( 0 ),
    is_done_ ( false )

    {
    }

    virtual ~ChannelBenchmark ( )
    {
    }

    virtual void run ( )
    {
        simul_name_ = params_["OutputPrefix"].get<std::string>( );

        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );

        // Turn off INFO log except on master proc.
        if ( rank_ != MASTER_RANK )
        {
            INFO = false;
            CONSOLE_OUTPUT_ACTIVE = false;
        }

        std::ofstream info_log ( ( simul_name_ + "_info_log" ).c_str ( ) );
        LogKeeper::get_log ( "info" ).set_target ( &( std::cout ) );
        std::ofstream debug_log ( ( simul_name_ + "_debug_log" ).c_str ( ) );
        LogKeeper::get_log ( "debug" ).set_target ( &( std::cout ) );

        CONSOLE_OUTPUT ( 0, "============================================================" );
        CONSOLE_OUTPUT ( 0, "==== ChannelBenchmark                                    ===" );
        CONSOLE_OUTPUT ( 0, "====    built using HiFlow3.                             ===" );
        CONSOLE_OUTPUT ( 0, "====                                                     ===" );
        CONSOLE_OUTPUT ( 0, "==== Engineering Mathematics and Computing Lab (EMCL)    ===" );
        CONSOLE_OUTPUT ( 0, "============================================================" );
        CONSOLE_OUTPUT ( 0, "" );

        // output parameters for debugging
        LOG_INFO ( "parameters", params_ );

        // setup timing report
        TimingScope::set_report ( &time_report_ );

        {
            TimingScope tscope ( "Setup" );

            setup_linear_algebra ( );

            read_mesh ( );

            // The simulation has two modes: stationary and
            // instationary. Which one is used, depends on the parameter
            // Instationary.SolveInstationary. In stationary mode, solve()
            // is only called once, whereas in instationary mode, it is
            // called several times via the time-stepping method implemented in run_time_loop() .
            solve_instationary_ = params_["Instationary"]["SolveInstationary"].get<bool>( );

            prepare ( );
        }

        if ( solve_instationary_ )
        {
            LOG_INFO ( "simulation", "Solving instationary problem" );
            run_time_loop ( );
        }
        else
        {
            LOG_INFO ( "simulation", "Solving stationary problem" );
            if ( use_hiflow_newton_ )
            {

                Timer timer;

                newton_.Solve ( &sol_ );
                CONSOLE_OUTPUT ( 1, "Newton ended with residual norm " << newton_.GetResidual ( )
                                 << " after " << newton_.iter ( ) << " iterations." );

                timer.stop ( );

                CONSOLE_OUTPUT ( 0, "" );
                CONSOLE_OUTPUT ( 1, "Measured time of interest " << ( int ) timer.get_duration ( ) / 60 << "m"
                                 << ( timer.get_duration ( ) - ( ( int ) timer.get_duration ( ) / 60 )*60 ) << "s" );
                CONSOLE_OUTPUT ( 0, "" );
                CONSOLE_OUTPUT ( 1, "Measured time in seconds " << timer.get_duration ( ) );

            }
            else
            {
                solve ( );
            }

            visualize ( );
            if ( compute_bench_quantities_ )
            {
                compute_dfg_benchmark_quantities ( );
                output_norms ( );
            }
        }

        CONSOLE_OUTPUT ( 0, "" );

        if ( rank_ == MASTER_RANK )
        {
            // Output time report
            TimingReportOutputVisitor visitor ( std::cout );
            time_report_.traverse_depth_first ( visitor );

            // Output results table
            std::vector< std::string > column_names;
            column_names.push_back ( "Time" );
            column_names.push_back ( "|u|_L2" );
            column_names.push_back ( "|p|_L2" );
            column_names.push_back ( "|u|_H1" );
            column_names.push_back ( "|p|_H1" );
            column_names.push_back ( "Fd" );
            column_names.push_back ( "Cd" );
            column_names.push_back ( "Fl" );
            column_names.push_back ( "Cl" );
            column_names.push_back ( "delta-P" );
            //  results_table_.print(std::cout, column_names);

            //  std::ofstream file((simul_name_ + "_results.csv").c_str());
            //  results_table_.print_csv(file, column_names);
        }

        LogKeeper::get_log ( "info" ).flush ( );
        LogKeeper::get_log ( "debug" ).flush ( );
        LogKeeper::get_log ( "info" ).set_target ( 0 );
        LogKeeper::get_log ( "debug" ).set_target ( 0 );

        CONSOLE_OUTPUT ( 0, "============================================================" );
    }

  private:

    const MPI_Comm& communicator ( )
    {
        return comm_;
    }

    int rank ( )
    {
        return rank_;
    }

    int num_partitions ( )
    {
        return num_partitions_;
    }

    PLATFORM la_platform ( ) const
    {
        return la_sys_.Platform;
    }

    IMPLEMENTATION la_implementation ( ) const
    {
        return la_impl_;
    }

    MATRIX_FORMAT la_matrix_format ( ) const
    {
        return la_matrix_format_;
    }

    // Read, refine and partition mesh.
    void read_mesh ( );

    // Set up datastructures and read in some parameters.
    void prepare ( );

    // Set up boundary conditions
    void prepare_bc ( );

    // Solve nonlinear problem
    void solve ( );

    // Update step in nonlinear solver
    void update_solution ( );
    void update_solution_naive ( );
    void update_solution_armijo ( );

    // Compute forcing term for inexact Newton method.
    void choose_forcing_term ( int iter );

    // Computate instationary solution by time-stepping method.
    void run_time_loop ( );

    // Compute exact solution for one constructed example, see docu of function
    // at end of file (where implemented)
    void eval_exact_sol ( );

    // Visualize the solution in a file. In stationary mode, the
    // filename contains 'stationary', in instationary mode, it contains the current time-step ts_.
    void visualize ( );

    // Compute L2 or H1 norm of variable defined through vars on the
    // master process.
    double compute_norm ( int norm_type, const std::vector<int>& vars );

    // Output various norms of the solution.
    void output_norms ( );

    virtual void Reinit ( );

    // Helper functions for nonlinear solver
    void ApplyFilter ( LAD::VectorType& u );
    virtual void EvalFunc ( const LAD::VectorType& in, LAD::VectorType* out );
    void compute_residual ( const LAD::VectorType& in, LAD::VectorType* out ); // updates res_ with the residual
    void compute_stationary_residual ( const LAD::VectorType& in, LAD::VectorType* out ); // residual computation in stationary mode
    void compute_instationary_residual ( const LAD::VectorType& in, LAD::VectorType* out ); // residual computation in instationary mode

    virtual void EvalGrad ( const LAD::VectorType& in, LAD::MatrixType* out );
    void compute_jacobian ( const LAD::VectorType& in, LAD::MatrixType* out ); // updates matrix_ with the jacobian matrix
    void compute_stationary_matrix ( const LAD::VectorType& in, LAD::MatrixType* out ); // jacobi matrix computation in stationary mode
    void compute_instationary_matrix ( const LAD::VectorType& in, LAD::MatrixType* out ); // jacobi matrix computation in instationary mode

    // Pressure filter: substracts the mean of the pressure from each
    // pressure dof in sol_ .
    void filter_pressure ( );

    // Linear algebra set up
    void setup_linear_algebra ( );

    // compute L2-Error and H1semi-Error
    void compute_errors ( );

    // compute difference between solution last and penultmate timestep
    void compute_difference ( );

    void compute_dfg_benchmark_quantities ( );
    double compute_drag_force ( );
    void compute_forces ( double& drag_force, double& lift_force );
    void compute_force_coefficients ( double drag_force, double lift_force, double& drag_coef, double& lift_coef ) const;

    void find_cylinder_boundary_dofs ( int cylinder_mat_num,
                                       int var,
                                       std::vector<int>& bdy_dofs );

    // MPI stuff
    MPI_Comm comm_;
    int rank_, num_partitions_;

    // Linear algebra stuff
    SYSTEM la_sys_;
    IMPLEMENTATION la_impl_;
    MATRIX_FORMAT la_matrix_format_;

    // Parameter data read in from file.
    PropertyTree params_;

    std::string simul_name_; // parameter 'OutputPrefix': prefix for output files

    // Time-stepping variables
    int ts_;
    double dt_;
    double alpha1_, alpha2_, alpha3_;

    // Flow model variables
    double Um_, H_, W_, rho_, nu_;
    double U_mean_, diam_; // (DFG-specific) mean velocity and diameter of cylinder.

    // Flag for pressure filter -- parameter 'UsePressureFilter'
    bool use_pressure_filter_;

    // Meshes
    MeshPtr mesh_;
    int refinement_level_;

    VectorSpace<double> space_;

    // linear algebra objects
    Couplings<double> couplings_;
    CMatrix matrix_;
    CVector sol_, prev_sol_, cor_, res_, pressure_correction_, exact_sol_, error_;

    // linear solver parameters
    int lin_max_iter;
    double lin_abs_tol;
    double lin_rel_tol;
    double lin_div_tol;
    int basis_size;

    // linear solver
    GMRES<LAD> gmres_;

    // nonlinear solver parameters
    int nls_max_iter;
    double nls_abs_tol;
    double nls_rel_tol;
    double nls_div_tol;
    double eta_; // forcing term
    std::vector<double> residual_history_, forcing_term_history_;
    bool do_armijo_update_;
    std::string forcing_strategy_;
    bool use_forcing_strategy_;

    //damping strategy paramters
    double theta_initial;
    double theta_min;
    double armijo_dec;
    double suff_dec;
    int max_armijo_ite;

    //forcing strategy parameters
    double eta_initial;
    double eta_max;
    double gamma_EW2;
    double alpha_EW2;

    // nonlinear solver
    Newton<LAD> newton_;

    StandardGlobalAssembler<double> global_asm_;

    bool is_done_, solve_instationary_, convergence_test_;

    std::vector<int> dirichlet_dofs_;
    std::vector<double> dirichlet_values_;

#ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> ilupp_;
    bool use_ilupp_;
#endif

    bool is_dfg_benchmark_;
    bool compute_bench_quantities_;
    bool use_hiflow_newton_;

    HierarchicalReport<TimingData> time_report_;
    Table results_table_;

    // CSV output
    CSVWriter<double> bench_quantity_writer_;
    std::vector<std::string> bench_names_;
    std::vector<double> bench_quantities_;
};

// program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    std::string param_filename ( PARAM_FILENAME );
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }

    try
    {
        ChannelBenchmark app ( param_filename );
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

void ChannelBenchmark::read_mesh ( )
{
    TimingScope tscope ( "read_mesh" );

    MeshPtr master_mesh;

    refinement_level_ = 0;
    const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );

    const bool use_bdd = params_["UseBoundaryDomainDescriptor"].get<bool>( false );

    Coordinate radius = 0.05;
    std::vector<Coordinate> center ( 2 );
#if DIMENSION == 2
    center[0] = 0.2;
#else
    center[0] = 0.5;
#endif
    center[1] = 0.2;
    CylinderDescriptor cyl ( radius, center );
#ifdef WITH_PARMETIS

    if ( rank ( ) == MASTER_RANK )
    {
#    if DIMENSION == 2
        const std::string mesh_name = params_["Mesh"]["Filename1"].get<std::string>( );
#    else
        const std::string mesh_name = params_["Mesh"]["Filename2"].get<std::string>( );
#    endif
        std::string mesh_filename = std::string ( DATADIR ) + mesh_name;

        master_mesh = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0 );

        CONSOLE_OUTPUT ( 1, "Read mesh with " << master_mesh->num_entities ( DIMENSION ) << " cells." );

        while ( refinement_level_ < initial_ref_lvl && master_mesh->num_entities ( DIMENSION ) < 8 * this->num_partitions_ )
        {
            master_mesh = master_mesh->refine ( );
#    if DIMENSION ==3
            if ( use_bdd )
                adapt_boundary_to_function ( master_mesh, cyl );
#    endif
            ++refinement_level_;
        }
        LOG_INFO ( "mesh", "Initial refinement level = " << refinement_level_ );
        CONSOLE_OUTPUT ( 1, "Refined mesh (level " << refinement_level_ << ") has "
                         << master_mesh->num_entities ( DIMENSION ) << " cells." );
    }

    MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
    // Distribute mesh using METIS
    ParMetisGraphPartitioner metis;

    MeshPtr local_mesh;
    if ( num_partitions_ <= 1 )
    {
        NaiveGraphPartitioner partitioner;
        std::cerr << "METIS not used, because number of partitions is "
                << num_partitions_ << ".\n";
        local_mesh = partition_and_distribute ( master_mesh,
                                                MASTER_RANK,
                                                comm_,
                                                &partitioner );
        if ( rank_ == MASTER_RANK )
        {
            LOG_INFO ( "Partitioner", "NAIVE" );
        }
    }
    else
    {
        local_mesh = partition_and_distribute ( master_mesh,
                                                MASTER_RANK,
                                                comm_,
                                                &metis );
        if ( rank_ == MASTER_RANK )
        {
            LOG_INFO ( "Partitioner", "METIS" );
        }
    }
    assert ( local_mesh != 0 );

    master_mesh.reset ( );

    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "Mesh", "Partitioned and distributed" );
    }

    while ( refinement_level_ < initial_ref_lvl )
    {
        local_mesh = local_mesh->refine ( );
        ++refinement_level_;
    }

    assert ( local_mesh != 0 );

    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "Mesh refinement level", refinement_level_ );
    }

    // Compute ghost cells.
    SharedVertexTable shared_verts;

    MeshPtr repartitioned_mesh;
    ParMetisGraphPartitioner parmetis_partitioner;
    repartitioned_mesh = repartition_mesh ( local_mesh, comm_, &parmetis_partitioner );
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "Repartitioning of mesh", "done" );
    }
    mesh_ = compute_ghost_cells ( *repartitioned_mesh, comm_, shared_verts );

    local_mesh.reset ( );

    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "Ghost cell computation", "done" );
    }

#else
    if ( rank ( ) == MASTER_RANK )
    {
#    if DIMENSION == 2
        const std::string mesh_name = params_["Mesh"]["Filename1"].get<std::string>( );
#    else
        const std::string mesh_name = params_["Mesh"]["Filename2"].get<std::string>( );
#    endif
        std::string mesh_filename = std::string ( DATADIR ) + mesh_name;

        master_mesh = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0 );

        CONSOLE_OUTPUT ( 1, "Read mesh with " << master_mesh->num_entities ( DIMENSION ) << " cells." );

        for ( int r = 0; r < initial_ref_lvl; ++r )
        {
            master_mesh = master_mesh->refine ( );
#    if DIMENSION ==3
            if ( use_bdd )
                adapt_boundary_to_function ( master_mesh, cyl );
#    endif
            ++refinement_level_;
        }
        LOG_INFO ( "mesh", "Initial refinement level = " << refinement_level_ );
        CONSOLE_OUTPUT ( 1, "Refined mesh (level " << refinement_level_ << ") has "
                         << master_mesh->num_entities ( DIMENSION ) << " cells." );
    }

    MeshPtr local_mesh = partition_and_distribute ( master_mesh, MASTER_RANK, comm_ );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
#endif
}

void ChannelBenchmark::prepare ( )
{
    TimingScope tscope ( "prepare" );

    // prepare timestep
    ts_ = 0;
    dt_ = params_["Instationary"]["Timestep"].get<double>( );

    // set the alpha coefficients correctly for the
    // Crank-Nicolson method.
    alpha1_ = 0.5 * dt_;
    alpha2_ = dt_;
    alpha3_ = 0.5 * dt_;

#ifdef WITH_ILUPP
    // prepare preconditioner
    use_ilupp_ = params_["LinearSolver"]["Preconditioning"].get<bool>( );
    if ( use_ilupp_ )
    {
        ilupp_.InitParameter ( params_["ILUPP"]["PreprocessingType"].get<int>( ),
                               params_["ILUPP"]["PreconditionerNumber"].get<int>( ),
                               params_["ILUPP"]["MaxMultilevels"].get<int>( ),
                               params_["ILUPP"]["MemFactor"].get<double>( ),
                               params_["ILUPP"]["PivotThreshold"].get<double>( ),
                               params_["ILUPP"]["MinPivot"].get<double>( ) );
    }
#endif

    // prepare problem parameters
    rho_ = params_["FlowModel"]["Density"].get<double>( );
    nu_ = params_["FlowModel"]["Viscosity"].get<double>( );

    Um_ = params_["FlowModel"]["InflowSpeed"].get<double>( );
    H_ = params_["FlowModel"]["InflowHeight"].get<double>( );
    W_ = params_["FlowModel"]["InflowWidth"].get<double>( );

    // prepare space
    std::vector< int > degrees ( DIMENSION + 1 );
    const int u_deg = params_["FiniteElements"]["VelocityDegree"].get<int>( );
    const int p_deg = params_["FiniteElements"]["PressureDegree"].get<int>( );
    for ( int c = 0; c < DIMENSION; ++c )
    {
        degrees.at ( c ) = u_deg;
    }
    degrees.at ( DIMENSION ) = p_deg;

    std::vector<bool> is_cg ( DIMENSION + 1, true );

    space_.Init ( degrees, *mesh_, is_cg, HIFLOW_CLASSIC );

    CONSOLE_OUTPUT ( 1, "Total number of dofs = " << space_.dof ( ).ndofs_global ( ) );

    for ( int p = 0; p < num_partitions ( ); ++p )
    {
        CONSOLE_OUTPUT ( 2, "Num dofs on process " << p << " = " << space_.dof ( ).ndofs_on_sd ( p ) );
    }

    // pressure filter
    use_pressure_filter_ = params_["UsePressureFilter"].get<bool>( );

    // prepare linear algebra structures
    couplings_.Clear ( );
    couplings_.Init ( communicator ( ), space_.dof ( ) );

    // prepare global assembler
    QuadratureSelection q_sel ( params_["QuadratureOrder"].get<int>( ) );
    global_asm_.set_quadrature_selection_function ( q_sel );

    // set DFG benchmark flag
    is_dfg_benchmark_ = params_["DFGbenchmark"].get<bool>( );
    compute_bench_quantities_ = params_["BenchQuantities"].get<bool>( );
    use_hiflow_newton_ = params_["NonlinearSolver"]["UseHiFlowNewton"].get<bool>( );

    if ( is_dfg_benchmark_ )
    {
#if DIMENSION == 2
        U_mean_ = 2. / 3. * Um_;
#elif DIMENSION == 3
        U_mean_ = 4. / 9. * Um_;
#endif
        diam_ = 0.1;
        CONSOLE_OUTPUT ( 1, "Reynolds number = " << U_mean_ * diam_ / nu_ );
    }

    // compute matrix graph

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

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows, sparsity.off_diagonal_cols );

    matrix_.Init ( communicator ( ), couplings_, la_platform ( ), la_implementation ( ), la_matrix_format ( ) );
    sol_.Init ( communicator ( ), couplings_, la_platform ( ), la_implementation ( ) );
    prev_sol_.Init ( communicator ( ), couplings_, la_platform ( ), la_implementation ( ) );
    cor_.Init ( communicator ( ), couplings_, la_platform ( ), la_implementation ( ) );
    res_.Init ( communicator ( ), couplings_, la_platform ( ), la_implementation ( ) );

    matrix_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                            vec2ptr ( sparsity.diagonal_cols ),
                            sparsity.diagonal_rows.size ( ),
                            vec2ptr ( sparsity.off_diagonal_rows ),
                            vec2ptr ( sparsity.off_diagonal_cols ),
                            sparsity.off_diagonal_rows.size ( ) );
    matrix_.Zeros ( );

    sol_.InitStructure ( );
    sol_.Zeros ( );
#if EXACTSOL == 1
    // debugging example with constructed right hand side
    eval_exact_sol ( );
#endif

    prev_sol_.InitStructure ( );
    prev_sol_.Zeros ( );

    cor_.InitStructure ( );
    cor_.Zeros ( );

    res_.InitStructure ( );
    res_.Zeros ( );

    // setup linear solver
    lin_max_iter = params_["LinearSolver"]["MaximumIterations"].get<int>( );
    lin_abs_tol = params_["LinearSolver"]["AbsoluteTolerance"].get<double>( );
    lin_rel_tol = params_["LinearSolver"]["RelativeTolerance"].get<double>( );
    lin_div_tol = params_["LinearSolver"]["DivergenceLimit"].get<double>( );
    basis_size = params_["LinearSolver"]["BasisSize"].get<int>( );

#ifdef WITH_ILUPP
    if ( use_ilupp_ )
    {
        gmres_.SetupPreconditioner ( ilupp_ );
        gmres_.InitParameter ( basis_size, "RightPreconditioning" );
    }
    else
    {
        gmres_.InitParameter ( basis_size, "NoPreconditioning" );
    }
#else
    gmres_.InitParameter ( basis_size, "NoPreconditioning" );
#endif

    gmres_.InitControl ( lin_max_iter, lin_abs_tol, lin_rel_tol, lin_div_tol );
    gmres_.SetupOperator ( matrix_ );

    // get nonlinear solver parameters from param file
    nls_max_iter = params_["NonlinearSolver"]["MaximumIterations"].get<int>( );
    nls_abs_tol = params_["NonlinearSolver"]["AbsoluteTolerance"].get<double>( );
    nls_rel_tol = params_["NonlinearSolver"]["RelativeTolerance"].get<double>( );
    nls_div_tol = params_["NonlinearSolver"]["DivergenceLimit"].get<double>( );
    do_armijo_update_ = params_["NonlinearSolver"]["ArmijoUpdate"].get<bool>( );
    forcing_strategy_ = params_["NonlinearSolver"]["ForcingStrategy"].get<std::string>( );
    use_forcing_strategy_ = ( forcing_strategy_ != "None" );
    eta_ = 1.e-4; // initial value of forcing term

    //get damping strategy parameters from param file
    theta_initial = params_["NonlinearSolver"]["ThetaInitial"].get<double>( );
    theta_min = params_["NonlinearSolver"]["ThetaMinimal"].get<double>( );
    armijo_dec = params_["NonlinearSolver"]["ArmijoDecrease"].get<double>( );
    suff_dec = params_["NonlinearSolver"]["SufficientDecrease"].get<double>( );
    max_armijo_ite = params_["NonlinearSolver"]["MaxArmijoIteration"].get<int>( );

    //get forcing strategy parameters from param file
    eta_initial = params_["NonlinearSolver"]["InitialValueForcingTerm"].get<double>( );
    eta_max = params_["NonlinearSolver"]["MaxValueForcingTerm"].get<double>( );
    gamma_EW2 = params_["NonlinearSolver"]["GammaParameterEW2"].get<double>( );
    alpha_EW2 = params_["NonlinearSolver"]["AlphaParameterEW2"].get<double>( );

    // setup nonlinear solver
    newton_.InitParameter ( &res_, &matrix_ );
    newton_.InitParameter ( Newton<LAD>::NewtonInitialSolutionOwn );
    newton_.InitControl ( nls_max_iter, nls_abs_tol, nls_rel_tol, nls_div_tol );
    newton_.SetOperator ( *this );
    newton_.SetLinearSolver ( gmres_ );

    //Damping strategy object
    if ( do_armijo_update_ )
    {
        ArmijoDamping<LAD>* Armijo_Damping = new ArmijoDamping<LAD>( theta_initial, theta_min, armijo_dec, suff_dec, max_armijo_ite );
        newton_.SetDampingStrategy ( *Armijo_Damping );
    }

    //Forcing strategy object
    if ( forcing_strategy_ == "EisenstatWalker1" )
    {
        EWForcing<LAD>* EW_Forcing = new EWForcing<LAD>( eta_initial, eta_max, 1 );
        newton_.SetForcingStrategy ( *EW_Forcing );
    }
    else if ( forcing_strategy_ == "EisenstatWalker2" )
    {
        EWForcing<LAD>* EW_Forcing = new EWForcing<LAD>( eta_initial, eta_max, 2, gamma_EW2, alpha_EW2 );
        newton_.SetForcingStrategy ( *EW_Forcing );
    }

    // prepare dirichlet BC
    prepare_bc ( );
}

void ChannelBenchmark::prepare_bc ( )
{
    TimingScope tscope ( "prepare_bc" );

    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    const int inflow_bdy = params_["Boundary"]["InflowMaterial"].get<int>( );
    const int outflow_bdy = params_["Boundary"]["OutflowMaterial"].get<int>( );
#if DIMENSION == 3
#    if EXACTSOL == 1
    // debugging example with constructed right hand side
    ExactSolChannelFlowBC3d bc[3] = { ExactSolChannelFlowBC3d ( 0 ),
                                     ExactSolChannelFlowBC3d ( 1 ),
                                     ExactSolChannelFlowBC3d ( 2 ) };
#    else
    ChannelFlowBC3d bc[3] = { ChannelFlowBC3d ( 0, W_, H_, Um_, inflow_bdy, outflow_bdy ),
                             ChannelFlowBC3d ( 1, W_, H_, Um_, inflow_bdy, outflow_bdy ),
                             ChannelFlowBC3d ( 2, W_, H_, Um_, inflow_bdy, outflow_bdy ) };
#    endif
    for ( int var = 0; var < DIMENSION; ++var )
    {
        compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                            dirichlet_dofs_, dirichlet_values_ );
    }
#else
    ChannelFlowBC2d bc[2] = { ChannelFlowBC2d ( 0, H_, Um_, inflow_bdy, outflow_bdy ),
                             ChannelFlowBC2d ( 1, H_, Um_, inflow_bdy, outflow_bdy ) };

    for ( int var = 0; var < DIMENSION; ++var )
    {
        compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                            dirichlet_dofs_, dirichlet_values_ );
    }
#endif

    // apply BC to initial solution
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }
}

void ChannelBenchmark::solve ( )
{
    TimingScope tscope ( "Nonlinear solver" );

    // Newton's method is used to solve the nonlinear problem.  The
    // functions compute_residual(...) and compute_jacobian(...)
    // update the variables matrix_ and res_, respectively.
    // The vector cor_ is set up to be used for the correction, and
    // the solution state is stored in sol_ .
    //

    residual_history_.clear ( );
    forcing_term_history_.clear ( );

    // the computation of the residual needs updated ghost DoFs
    sol_.UpdateCouplings ( );
    compute_residual ( sol_, &res_ );

    int iter = 0;
    const double initial_res_norm = res_.Norm2 ( );
    residual_history_.push_back ( initial_res_norm );

    LOG_INFO ( "nonlinear", "Nonlinear solver starts with residual norm " << initial_res_norm );
    double res_norm = initial_res_norm;

    CONSOLE_OUTPUT ( 2, "Newton solver starts with residual norm = " << res_norm );

    // initial value of forcing term

    forcing_term_history_.push_back ( eta_ );

    while ( iter < nls_max_iter
            && res_norm > nls_abs_tol
            && res_norm > nls_rel_tol * initial_res_norm )
    {
        TimingScope tscope_iter ( iter );

        // Solve DF * cor = res

        // the computation of the Jacobian needs updated ghost DoFs,
        // but at this point of the code, the input vector is always
        // up-to-date
        compute_jacobian ( sol_, &matrix_ );

        {
            TimingScope tscope_solve ( "Linear solve" );

            if ( use_forcing_strategy_ )
            {
                // Set linear solver parameters based on forcing term.
                choose_forcing_term ( iter );
                CONSOLE_OUTPUT ( 4, "Linear tolerance = " << eta_ * res_norm );
                gmres_.InitControl ( lin_max_iter, eta_ * res_norm, 1.e-20, lin_div_tol );
            }

            // Compute correction cor
            cor_.Zeros ( );
            gmres_.Solve ( res_, &cor_ );
            CONSOLE_OUTPUT ( 3, "Linear solver computed correction in " << gmres_.iter ( ) << " iterations." );
            CONSOLE_OUTPUT ( 3, "Residual norm for correction = " << gmres_.res ( ) );
        }

        // Update solution
        update_solution ( );

        if ( use_pressure_filter_ )
        {
            ApplyFilter ( sol_ );
        }

        // the computation of the residual needs updated ghost DoFs
        sol_.UpdateCouplings ( );
        compute_residual ( sol_, &res_ );
        res_norm = res_.Norm2 ( );
        residual_history_.push_back ( res_norm );

        ++iter;
        CONSOLE_OUTPUT ( 2, "Newton step " << iter << " ends with residual norm = " << res_norm );
    }
    CONSOLE_OUTPUT ( 3, "Forcing term history = "
                     << string_from_range ( forcing_term_history_.begin ( ), forcing_term_history_.end ( ) ) );
    CONSOLE_OUTPUT ( 3, "Residual history = "
                     << string_from_range ( residual_history_.begin ( ), residual_history_.end ( ) ) );

    LOG_INFO ( "Relative residual stop criterium", nls_rel_tol * initial_res_norm );
    LOG_INFO ( "Nonlinear solver residual", res_norm );
    LOG_INFO ( "Nonlinear solver steps", iter );

    CONSOLE_OUTPUT ( 2, "Newton solver finished after " << iter << " iterations." );
}

void ChannelBenchmark::update_solution ( )
{
    TimingScope tscope_update ( "Update" );

    if ( do_armijo_update_ )
    {
        update_solution_armijo ( );
    }
    else
    {
        update_solution_naive ( );
    }
}

void ChannelBenchmark::update_solution_naive ( )
{
    sol_.Axpy ( cor_, -1. );
}

void ChannelBenchmark::update_solution_armijo ( )
{
    // update solution using simple armijo update
    const double initial_res = res_.Norm2 ( );
    double theta = 1.; // initial theta
    const double suff_dec = 1.e-4; // sufficient decrease
    const int max_armijo_iter = 10;
    const double armijo_dec = 0.5; // armijo decrease factor
    const double theta_min = 0.1;

    double res_norm = initial_res;

    for ( int i = 0; i < max_armijo_iter; ++i )
    {
        // update solution
        sol_.Axpy ( cor_, -theta );

        // the computation of the residual needs updated ghost DoFs
        sol_.UpdateCouplings ( );
        compute_residual ( sol_, &res_ );
        res_norm = res_.Norm2 ( );

        CONSOLE_OUTPUT ( 4, "\t\t\tRes norm = " << res_norm
                         << ", criterion = " << ( ( 1 - suff_dec * ( 1 - eta_ ) ) * initial_res ) );

        if ( res_norm <= ( ( 1 - suff_dec * ( 1 - eta_ ) ) * initial_res ) )
        {
            // condition satisfied -- we are done
            break;
        }

        // condition not satisfied

        // restore solution
        sol_.Axpy ( cor_, theta );

        theta *= armijo_dec;
        eta_ = 1 - theta * ( 1 - eta_ );

        // break if theta becomes too small
        if ( theta < theta_min )
        {
            break;
        }
    }
    CONSOLE_OUTPUT ( 3, "Armijo update last theta = " << theta );
}

void ChannelBenchmark::choose_forcing_term ( int iter )
{
    if ( forcing_strategy_ == "EisenstatWalker2" )
    {
        // Eisenstat-Walker Choice 2

        if ( iter == 0 )
        {
            // Simply set initial forcing term on first iteration.
            eta_ = 0.9;
        }
        else
        {
            const double golden_ratio = 0.5 * ( 1 + std::sqrt ( 5.0 ) );
            const double alpha = 2.;
            const double gamma = 0.9;
            const double eta_max = 0.9;

            const double last_res = residual_history_.at ( iter - 1 );
            const double res_norm = residual_history_.at ( iter );

            const double last_eta = forcing_term_history_.at ( iter - 1 );

            CONSOLE_OUTPUT ( 4, "last_res = " << last_res << ", curr_res = " << res_norm );

            eta_ = gamma * std::pow ( res_norm / last_res, alpha );

            const double threshold = std::pow ( last_eta, golden_ratio );
            if ( threshold > 0.1 )
            {
                eta_ = std::max ( eta_, threshold );
            }

            // safeguard
            eta_ = std::min ( eta_, eta_max );
        }
    }
    else if ( forcing_strategy_ == "Constant" )
    {
        eta_ = params_["NonlinearSolver"]["ConstantForcingTerm"].get<double>( );
    }
    else
    {
        std::cerr << "Unknown forcing strategy " << forcing_strategy_ << "\n";
        exit ( 1 );
    }

    forcing_term_history_.push_back ( eta_ );
    CONSOLE_OUTPUT ( 3, "Forcing term eta = " << eta_ );
}

void ChannelBenchmark::run_time_loop ( )
{
    TimingScope tscope ( "Timestepping loop" );

#ifdef WITH_HDF5
    if ( params_["Backup"]["Restore"].get<bool>( ) )
    {
        ts_ = params_["Backup"]["LastTimeStep"].get<int>( );
        std::stringstream filename;
        filename << this->num_partitions_ << "_" << params_["Backup"]["Filename"].get<std::string>( );
        const std::string backup_name = filename.str ( );
        std::ostringstream vec_name;
        vec_name << "sol_" << ts_;
        prev_sol_.ReadHDF5 ( backup_name, "backup", vec_name.str ( ) );
        prev_sol_.UpdateCouplings ( );
        CONSOLE_OUTPUT ( 1, "Restarting from backup in file " << backup_name << " after timestep " << ts_ );
    }
#endif

    // Visualize initial solution.
    if ( ts_ == 0 )
    {
        if ( rank ( ) == MASTER_RANK )
        {
            results_table_.insert ( "Time", 0. );

            results_table_.insert ( "Fd", 0. );
            results_table_.insert ( "Cd", 0. );
            results_table_.insert ( "Fl", 0. );
            results_table_.insert ( "Cl", 0. );
            results_table_.insert ( "delta-P", 0. );
        }

        visualize ( );

        if ( compute_bench_quantities_ )
        {
            output_norms ( );
        }

#ifdef WITH_HDF5
        std::stringstream filename;
        filename << this->num_partitions_ << "_" << params_["Backup"]["Filename"].get<std::string>( );
        const std::string backup_name = filename.str ( );
        std::ostringstream vec_name;
        vec_name << "sol_" << ts_;
        prev_sol_.WriteHDF5 ( backup_name, "backup", vec_name.str ( ) );
#endif
    }

    const double end_time = params_["Instationary"]["Endtime"].get<double>( );
    LOG_INFO ( "timestep", "End time = " << end_time );
    LOG_INFO ( "timestep", "Step length = " << dt_ );

    // Set up CSV output of benchmarking quantities
    bench_names_.push_back ( "Timestep" );
    bench_names_.push_back ( "Iterations (Newton)" );
    bench_names_.push_back ( "Time to compute time-step [s]" );
    bench_names_.push_back ( "Time to compute residuals [s]" );
    bench_names_.push_back ( "Time to compute Jacobians and setup preconditioner [s]" );
    bench_names_.push_back ( "(Total) Number of GMRES iterations per time-step" );

    bench_quantities_.resize ( bench_names_.size ( ), 0. );

    // Initialize CSV Writer
    std::stringstream bench_file;
    bench_file
            << num_partitions_
            << "_"
            << dt_
            << "_"
            << "benchmarking_quantities.csv";

    if ( this->rank_ == MASTER_RANK )
    {
        LOG_INFO ( "Benchmarking quantities file", bench_file.str ( ) );
    }
    bench_quantity_writer_.InitFilename ( bench_file.str ( ) );

    if ( params_["Backup"]["Restore"].get<bool>( ) && ts_ != 0 )
    {
        std::vector< std::vector < double > > stored_bench;
        std::vector< std::vector < double > > preserve_bench;
        bench_quantity_writer_.read ( stored_bench );

        for ( int i = 0; i < static_cast < int > ( stored_bench.size ( ) ); ++i )
        {
            if ( stored_bench[i][0] <= ts_ )
            {
                preserve_bench.push_back ( stored_bench[i] );
            }
        }

        // keep only bench quantities of former timesteps
        if ( rank_ == MASTER_RANK )
        {
            bench_quantity_writer_.Init ( bench_names_ );

            for ( int i = 0; i < static_cast < int > ( preserve_bench.size ( ) ); ++i )
            {
                bench_quantity_writer_.write ( preserve_bench[i] );
            }
        }
    }
    else
    {
        if ( this->rank_ == MASTER_RANK )
        {
            bench_quantity_writer_.Init ( bench_names_ );
        }
    }

    // Set first timestep to compute.

    CONSOLE_OUTPUT ( 1, "Starting time loop from t = " << ts_ * dt_ << " to " << end_time << " with timestep " << dt_ );
    ++ts_;

    // Crank-Nicolson time-stepping method. At the beginning of each
    // time-step, the solution from the previous time-step is stored
    // in prev_sol_, which is used in InstationaryFlowAssembler. The
    // variable ts_ is used to keep track of the current
    // time-step. The solution is visualized at the end of each
    // time-step, in order to be able to animate it in Paraview.
    while ( ts_ * dt_ <= end_time )
    {
        TimingScope tscope ( ts_ );
        CONSOLE_OUTPUT ( 1, "Solving timestep " << ts_ << " (for t = " << ts_ * dt_ << ")" );
        LOG_INFO ( "timestep", "Solving time step " << ts_ );

        bench_quantities_[1] = bench_quantities_[2] = bench_quantities_[3] = bench_quantities_[4] = bench_quantities_[5] = 0.;

        // check benchmarking quantities
        bench_quantities_[0] = ts_;

        Timer time_step_timer;
        time_step_timer.start ( );

        if ( use_hiflow_newton_ )
        {
            newton_.Solve ( &sol_ );
            CONSOLE_OUTPUT ( 1, "Newton ended with residual norm " << newton_.GetResidual ( )
                             << " after " << newton_.iter ( ) << " iterations." );
            bench_quantities_[1] = newton_.iter ( );
            bench_quantities_[5] += gmres_.iter ( );
        }
        else
        {
            solve ( );
        }

        time_step_timer.stop ( );
        bench_quantities_[2] = time_step_timer.get_duration ( );

        if ( this->rank_ == MASTER_RANK )
        {
            bench_quantity_writer_.write ( bench_quantities_ );
        }

        prev_sol_.CloneFrom ( sol_ );
        // this also clones the ghost DoFs, so we don't need to call
        // UpdateCouplings()

#ifdef WITH_HDF5
        std::stringstream filename;
        filename << this->num_partitions_ << "_" << params_["Backup"]["Filename"].get<std::string>( );
        const std::string backup_name = filename.str ( );
        std::ostringstream vec_name;
        vec_name << "sol_" << ts_;
        prev_sol_.WriteHDF5 ( backup_name, "backup", vec_name.str ( ) );
#endif

        results_table_.insert ( "Time", ts_ * dt_ );

        LOG_INFO ( "timestep", "Visualizing solution at time " << ts_ * dt_ << " (time-step " << ts_ << ")" );
        visualize ( );

        if ( compute_bench_quantities_ )
        {
            compute_dfg_benchmark_quantities ( );
            output_norms ( );
        }

        ++ts_;
    }
}

void ChannelBenchmark::visualize ( )
{
    TimingScope tscope ( "Visualization" );

    // Setup visualization object.
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    std::stringstream input;

    input << simul_name_ << "_solution";

    if ( solve_instationary_ )
    {
        if ( ts_ < 10 )
            input << "000" << ts_;
        else if ( ts_ < 100 )
            input << "00" << ts_;
        else if ( ts_ < 1000 )
            input << "0" << ts_;
        else
            input << "" << ts_;
    }
    else
    {
        input << "_stationary";
    }

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

    sol_.UpdateCouplings ( );

    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 0 ), "u" );
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 1 ), "v" );
#if DIMENSION == 2
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 2 ), "p" );
#elif DIMENSION == 3
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 2 ), "w" );
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 3 ), "p" );
#endif

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );
}

double ChannelBenchmark::compute_norm ( int norm_type, const std::vector<int>& vars )
{
    double local_norm = -1.e30, global_norm = 0.;
    sol_.UpdateCouplings ( );
    switch ( norm_type )
    {
        case 0: // L2-norm
        {
            L2NormIntegratorPp L2_int ( sol_, vars );
            global_asm_.integrate_scalar ( space_, L2_int, local_norm );
            break;
        }
        case 1: // H1-seminorm
        {
            H1semiNormIntegratorPp H1_int ( sol_, vars );
            global_asm_.integrate_scalar ( space_, H1_int, local_norm );
            break;
        }
        default:
            std::cerr << "unknown type of norm!\n";
            assert ( false );
    };

    // NB: global value will only be returned on master proc -- others will return 0.

    MPI_Reduce ( &local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK, comm_ );

    return std::sqrt ( global_norm );
}

void ChannelBenchmark::output_norms ( )
{
    TimingScope tscope ( "Norm computation" );
    std::vector<int> vel_vars, p_var;
    vel_vars.push_back ( 0 );
    vel_vars.push_back ( 1 );
    vel_vars.push_back ( 2 );
    p_var.push_back ( 3 );

    const double L2_vel_norm = compute_norm ( 0, vel_vars );
    const double L2_p_norm = compute_norm ( 0, p_var );
    const double H1_vel_norm = compute_norm ( 1, vel_vars );
    const double H1_p_norm = compute_norm ( 1, p_var );

    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "L2-norm of velocity", L2_vel_norm );
        LOG_INFO ( "L2-norm of pressure", L2_p_norm );
        LOG_INFO ( "H1-seminorm of velocity", H1_vel_norm );
        LOG_INFO ( "H1-seminorm of pressure", H1_p_norm );

        results_table_.insert ( "|u|_L2", L2_vel_norm );
        results_table_.insert ( "|p|_L2", L2_p_norm );
        results_table_.insert ( "|u|_H1", H1_vel_norm );
        results_table_.insert ( "|p|_H1", H1_p_norm );
    }
}

void ChannelBenchmark::Reinit ( )
{
    //  prev_sol_.CloneFrom(sol_);
}

void ChannelBenchmark::EvalFunc ( const LAD::VectorType& in,
                                  LAD::VectorType* out )
{
    Timer assembly_timer;
    assembly_timer.start ( );
    compute_residual ( in, out );
    //    out->Scale(-1.0);
    assembly_timer.stop ( );

    bench_quantities_[3] += assembly_timer.get_duration ( );
}

void ChannelBenchmark::compute_residual ( const LAD::VectorType& in,
                                          LAD::VectorType* out )
{
    TimingScope tscope ( "Compute Residual" );

    // the evaluation of the residual needs updated ghost DoFs,
    // so make sure you call UpdateCouplings() on the input vector
    // before you call this function!

    if ( solve_instationary_ )
    {
        compute_instationary_residual ( in, out );
    }
    else
    {
        compute_stationary_residual ( in, out );
    }

    // correct BC -- set Dirichlet dofs to 0
    if ( !dirichlet_dofs_.empty ( ) )
    {
        std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
        out->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                         vec2ptr ( zeros ) );
    }
}

void ChannelBenchmark::compute_stationary_residual ( const LAD::VectorType& in,
                                                     LAD::VectorType* out )
{
    StationaryFlowAssembler local_asm ( in, nu_, rho_ );
    global_asm_.assemble_vector ( space_, local_asm, *out );
}

void ChannelBenchmark::compute_instationary_residual ( const LAD::VectorType& in,
                                                       LAD::VectorType* out )
{
    InstationaryFlowAssembler local_asm ( nu_, rho_ );
    local_asm.set_newton_solution ( &in );
    // the computation of the instationary residual also needs updated
    // ghost DoFs of the previous solution, so make sure you call
    // UpdateCouplings() on prev_sol_ before calling this function
    local_asm.set_time_solution ( &prev_sol_ );
    local_asm.set_time_stepping_weights ( alpha1_, alpha2_, alpha3_ );

    global_asm_.assemble_vector ( space_, local_asm, *out );
}

void ChannelBenchmark::EvalGrad ( const LAD::VectorType& in,
                                  LAD::MatrixType* out )
{
    bench_quantities_[5] += gmres_.iter ( );

    Timer assembly_timer;
    assembly_timer.start ( );

    compute_jacobian ( in, out );

    assembly_timer.stop ( );
    bench_quantities_[4] += assembly_timer.get_duration ( );
}

void ChannelBenchmark::compute_jacobian ( const LAD::VectorType& in,
                                          LAD::MatrixType* out )
{
    {
        TimingScope tscope ( "Compute Jacobian" );

        // the computation of the Jacobian needs updated ghost DoFs,
        // so make sure you call UpdateCouplings() on the input vector
        // before calling this function!

        if ( solve_instationary_ )
        {
            compute_instationary_matrix ( in, out );
        }
        else
        {
            compute_stationary_matrix ( in, out );
        }

        // correct BC -- set Dirichlet rows to identity
        if ( !dirichlet_dofs_.empty ( ) )
        {
            out->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
        }
    }

#ifdef WITH_ILUPP
    {
        TimingScope tscope ( "ILU++ Factorization" );
        if ( use_ilupp_ )
        {
            ilupp_.SetupOperator ( *out );
        }
    }
#endif
}

void ChannelBenchmark::compute_stationary_matrix ( const LAD::VectorType& in,
                                                   LAD::MatrixType* out )
{
    StationaryFlowAssembler local_asm ( in, nu_, rho_ );
    global_asm_.assemble_matrix ( space_, local_asm, *out );
}

void ChannelBenchmark::compute_instationary_matrix ( const LAD::VectorType& in,
                                                     LAD::MatrixType* out )
{
    InstationaryFlowAssembler local_asm ( nu_, rho_ );

    local_asm.set_newton_solution ( &in );
    // the computation of the Jacobian in the instationary case also
    // also needs updated ghost DoFs of the previous solution, so
    // make sure you call UpdateCouplings() on prev_sol_ before
    // calling this function
    local_asm.set_time_solution ( &prev_sol_ );
    local_asm.set_time_stepping_weights ( alpha1_, alpha2_, alpha3_ );

    global_asm_.assemble_matrix ( space_, local_asm, *out );
}

//////////////// Pressure Filtering ////////////////

struct PressureIntegral : private AssemblyAssistant<DIMENSION, Scalar>
{

    PressureIntegral ( const CoupledVector<Scalar>& sol ) : sol_ ( sol )
    {
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            double& pressure )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        evaluate_fe_function ( sol_, DIMENSION, p_ );

        const int num_q = num_quadrature_points ( );

        // loop over quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            pressure += wq * p_[q] * dJ;
        }
    }

    const CoupledVector<Scalar>& sol_;
    FunctionValues<double> p_;
};

struct VolumeIntegral : private AssemblyAssistant<DIMENSION, double>
{

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, double& vol )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        const int num_q = num_quadrature_points ( );

        // loop over quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );
            vol += wq * dJ;
        }
    }
};

void ChannelBenchmark::ApplyFilter ( LAD::VectorType& u )
{
    if ( !use_pressure_filter_ ) return;
    double recv;
    u.UpdateCouplings ( );
    PressureIntegral int_p ( u );
    double total_pressure;
    global_asm_.integrate_scalar ( space_, int_p, total_pressure );

    MPI_Allreduce ( &total_pressure, &recv, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    total_pressure = recv;

    double integrated_vol;
    VolumeIntegral vol_int;
    global_asm_.integrate_scalar ( space_, vol_int, integrated_vol );

    MPI_Allreduce ( &integrated_vol, &recv, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    integrated_vol = recv;

    const double average_pressure = total_pressure / integrated_vol;

    LOG_INFO ( "pressure_filter", "Average pressure before filter = " << average_pressure );

    pressure_correction_.CloneFromWithoutContent ( u );
    pressure_correction_.Zeros ( );

    // set value for pressure dofs to average pressure
    std::vector<int> cell_p_dofs;
    std::vector<int> local_p_dofs;
    for ( EntityIterator it = mesh_->begin ( DIMENSION ), end = mesh_->end ( DIMENSION ); it != end; ++it )
    {
        cell_p_dofs.clear ( );
        space_.GetDofIndices ( DIMENSION, *it, &cell_p_dofs );
        for ( int i = 0, sz = cell_p_dofs.size ( ); i < sz; ++i )
        {
            if ( space_.dof ( ).is_dof_on_sd ( cell_p_dofs[i] ) )
            {
                local_p_dofs.push_back ( cell_p_dofs[i] );
            }
        }
    }

    std::sort ( local_p_dofs.begin ( ), local_p_dofs.end ( ) );
    std::unique ( local_p_dofs.begin ( ), local_p_dofs.end ( ) );

    // remove average pressure from solution
    std::vector<double> p_correction_values ( local_p_dofs.size ( ) );
    std::fill ( p_correction_values.begin ( ), p_correction_values.end ( ), average_pressure );

    pressure_correction_.SetValues ( vec2ptr ( local_p_dofs ), local_p_dofs.size ( ), vec2ptr ( p_correction_values ) );

    u.Axpy ( pressure_correction_, -1. );

    u.UpdateCouplings ( );
    PressureIntegral int_p_check ( u );
    global_asm_.integrate_scalar ( space_, int_p_check, total_pressure );
    MPI_Allreduce ( &total_pressure, &recv, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    total_pressure = recv;
    LOG_INFO ( "pressure_filter", "Average pressure after filter = " << total_pressure / integrated_vol );
}

class ForceIntegral : private AssemblyAssistant<DIMENSION, double>
{
  public:

    enum FORCE_TYPE
    {
        DRAG = 0, LIFT = 1
    };

    ForceIntegral ( double nu, const CoupledVector<Scalar>* sol, const std::vector<int>& bdy_dofs, FORCE_TYPE type )
    : nu_ ( nu ), x_var_ ( type == DRAG ? 0 : 1 ), bdy_dofs_ ( bdy_dofs ), sol_ ( sol )
    {
    } // rho assumed to be one

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, double& val )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        std::vector<int> dofs;
        element.get_dof_indices ( dofs );

        std::sort ( dofs.begin ( ), dofs.end ( ) );

        std::vector<int> bdofs; // cylinder boundary dofs on the cell

        // find dofs on cell that also lie on cylinder boundary
        std::set_intersection ( dofs.begin ( ), dofs.end ( ), bdy_dofs_.begin ( ), bdy_dofs_.end ( ),
                                std::back_inserter ( bdofs ) );

        if ( !bdofs.empty ( ) )
        {
            dofs.clear ( );
            element.get_dof_indices ( dofs );

            // We compute function values here only, since otherwise they will not be needed.
            recompute_function_values ( );

            for ( std::vector<int>::const_iterator d = bdofs.begin ( ), d_end = bdofs.end ( );
                  d != d_end; ++d )
            {
                // Find local dof number for *d
                std::vector<int>::iterator i_it =
                        std::find ( dofs.begin ( ) + dof_index ( 0, x_var_ ),
                                    dofs.begin ( ) + dof_index ( 0, x_var_ ) + num_dofs ( x_var_ ),
                                    *d );
                const int i = std::distance ( dofs.begin ( ) + dof_index ( 0, x_var_ ), i_it );
                assert ( i >= 0 );
                assert ( i < num_dofs ( x_var_ ) );

                const int num_q = num_quadrature_points ( );

                for ( int q = 0; q < num_q; ++q )
                {
                    const double wq = w ( q );
                    const double dJ = std::abs ( detJ ( q ) );

                    val -= wq *
                            ( nu_ * dot ( grad_phi ( i, q, x_var_ ), grad_u_[x_var_][q] )
                            - p_[q] * grad_phi ( i, q, x_var_ )[x_var_] )
                            * dJ;

                    // non-linear term involves all components of u_ and grad_u_.
                    for ( int v = 0; v < DIMENSION; ++v )
                    {
                        val -= wq * ( u_[v][q] * grad_u_[v][q][x_var_] * phi ( i, q, x_var_ ) ) * dJ;
                    }
                }
            }
        }
    }

  private:

    void recompute_function_values ( )
    {
        for ( int d = 0; d < DIMENSION; ++d )
        {
            u_[d].clear ( );
            grad_u_[d].clear ( );

            evaluate_fe_function ( *sol_, d, u_[d] );
            evaluate_fe_function_gradients ( *sol_, d, grad_u_[d] );
        }
        p_.clear ( );
        evaluate_fe_function ( *sol_, DIMENSION, p_ );
    }

    double nu_;
    int surface_material_num_, x_var_;
    FunctionValues<double> u_[DIMENSION], p_;
    FunctionValues< Vec<DIMENSION, double> > grad_u_[DIMENSION];
    const std::vector<int>& bdy_dofs_;
    const CoupledVector<Scalar>* sol_;

};

void ChannelBenchmark::compute_forces ( double& drag_force, double& lift_force )
{
    std::vector<int> dof_ids;
    std::vector<double> dof_values; // dummy argument -- not used.

    const int cylinder_material = params_["Boundary"]["CylinderMaterial"].get<int>( );
    double local_force[2] = { 0., 0. };

    find_cylinder_boundary_dofs ( cylinder_material, 0, dof_ids );
    LOG_INFO ( "dfg_benchmark", "Total number of dof indices on cylinder bdy 0 = " << dof_ids.size ( ) );

    ForceIntegral int_F0 ( nu_, &sol_, dof_ids, ForceIntegral::DRAG );
    global_asm_.integrate_scalar ( space_, int_F0, local_force[0] );

    find_cylinder_boundary_dofs ( cylinder_material, 1, dof_ids );
    LOG_INFO ( "dfg_benchmark", "Total number of dof indices on cylinder bdy 1 = " << dof_ids.size ( ) );

    ForceIntegral int_F1 ( nu_, &sol_, dof_ids, ForceIntegral::LIFT );
    global_asm_.integrate_scalar ( space_, int_F1, local_force[1] );

    LOG_INFO ( "dfg_benchmark", "On process " << rank ( ) << ", Fd = " << local_force[0]
               << " and " << " Fl = " << local_force[1] );

    double recv[2] = { 0., 0. };
    MPI_Allreduce ( &local_force[0], &recv[0], 2, MPI_DOUBLE, MPI_SUM, comm_ );

    drag_force = recv[0];
    lift_force = recv[1];
}

void ChannelBenchmark::compute_force_coefficients ( double drag_force, double lift_force, double& drag_coef, double& lift_coef ) const
{
#if DIMENSION == 2
    drag_coef = 2. * drag_force / ( U_mean_ * U_mean_ * diam_ );
    lift_coef = 2. * lift_force / ( U_mean_ * U_mean_ * diam_ );
#elif DIMENSION == 3
    drag_coef = 2. * drag_force / ( U_mean_ * U_mean_ * diam_ * H_ );
    lift_coef = 2. * lift_force / ( U_mean_ * U_mean_ * diam_ * H_ );
#endif
}

void ChannelBenchmark::compute_dfg_benchmark_quantities ( )
{
    TimingScope tscope ( "Compute DFG Benchmark Quantities" );

    // compute force coefficients
    // NB: rho is assumed to be one
    double Fd, Fl, Cd, Cl;
    sol_.UpdateCouplings ( );
    compute_forces ( Fd, Fl );
    compute_force_coefficients ( Fd, Fl, Cd, Cl );

    if ( rank ( ) == MASTER_RANK )
    {
        LOG_INFO ( "dfg_benchmark", "Drag force = " << Fd );
        LOG_INFO ( "dfg_benchmark", "Drag coefficient = " << Cd );

        LOG_INFO ( "dfg_benchmark", "Lift force = " << Fl );
        LOG_INFO ( "dfg_benchmark", "Lift coefficient = " << Cl );

        results_table_.insert ( "Fd", Fd );
        results_table_.insert ( "Cd", Cd );
        results_table_.insert ( "Fl", Fl );
        results_table_.insert ( "Cl", Cl );
    }

    // compute pressure drop
    std::vector< std::vector<double> > p ( 2, std::vector<double>( DIMENSION ) );

#if DIMENSION == 2
    p[0][0] = 0.15;
    p[0][1] = 0.2;
    p[1][0] = 0.25;
    p[1][1] = 0.2;
#elif DIMENSION == 3
    p[0][0] = 0.448;
    p[0][1] = 0.2;
    p[0][2] = 0.205;
    p[1][0] = 0.5502;
    p[1][1] = 0.2;
    p[1][2] = 0.205;
#endif
    PointEvaluator<double> p_eval ( space_ );

    double recv_p[2];

    for ( int i = 0; i < 2; ++i )
    {
        p_eval.evaluate_fun_global ( EvalFeFunction<LAD>( space_, sol_, DIMENSION ), p[i], recv_p[i], communicator ( ) );
        LOG_INFO ( "dfg_benchmark", "Computed pressure " << i << " on process "
                   << rank_ << ": " << recv_p[i] );
    }

    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "dfg_benchmark", "Pressure difference delta-p = " << recv_p[0] - recv_p[1] );
        results_table_.insert ( "delta-P", recv_p[0] - recv_p[1] );
    }
}

void ChannelBenchmark::find_cylinder_boundary_dofs ( int cylinder_mat_num,
                                                     int var,
                                                     std::vector<int>& bdy_dofs )
{
    bdy_dofs.clear ( );

    const Mesh& mesh = space_.mesh ( );
    const TDim tdim = mesh.tdim ( );

    MeshPtr boundary_mesh = mesh.extract_boundary_mesh ( );
    const bool is_sequential = ( num_partitions ( ) == 1 );
    if ( !is_sequential )
    {
        assert ( mesh.has_attribute ( "_sub_domain_", tdim ) );
    }

    std::vector<doffem::DofID> dofs_on_face;
    for ( EntityIterator it_boundary = boundary_mesh->begin ( tdim - 1 );
          it_boundary != boundary_mesh->end ( tdim - 1 );
          ++it_boundary )
    {
        // get id of boundary face
        const Id boundary_id = it_boundary->id ( );

        // check if the boundary face exists and get the location
        // where the entity number should be stored
        int face_number;
        const bool check = mesh.find_entity ( tdim - 1, boundary_id, &face_number );
        assert ( check );

        // Get the face to be able to access to the data associated with the face
        Entity face = mesh.get_entity ( tdim - 1, face_number );
        if ( face.get_material_number ( ) != cylinder_mat_num )
        {
            continue;
        }

        IncidentEntityIterator cell = face.begin_incident ( tdim );

        // loop over all faces of the cell to get the local face index for identifying the dofs
        int local_face_number = 0;
        for ( IncidentEntityIterator global_face = cell->begin_incident ( tdim - 1 );
              global_face != cell->end_incident ( tdim - 1 );
              ++global_face )
        {
            // if the global face id equals the boundary id the local face index is found
            if ( global_face->id ( ) == boundary_id )
            {
                break;
            }
            else
            {
                local_face_number++;
            }
        }

        dofs_on_face.clear ( );
        space_.dof ( ).get_dofs_on_subentity ( var,
                                               cell->index ( ),
                                               tdim - 1,
                                               local_face_number,
                                               dofs_on_face );
        bdy_dofs.insert ( bdy_dofs.end ( ), dofs_on_face.begin ( ), dofs_on_face.end ( ) );
    }

    std::sort ( bdy_dofs.begin ( ), bdy_dofs.end ( ) );
    std::vector<int>::iterator new_end = std::unique ( bdy_dofs.begin ( ), bdy_dofs.end ( ) );
    bdy_dofs.resize ( std::distance ( bdy_dofs.begin ( ), new_end ) );
}

void ChannelBenchmark::setup_linear_algebra ( )
{
    TimingScope tscope ( "setup_linear_algebra" );
    const std::string platform_str = params_["LinearAlgebra"]["Platform"].get<std::string>( );
    if ( platform_str == "CPU" )
    {
        la_sys_.Platform = CPU;
    }
    else if ( platform_str == "GPU" )
    {
        la_sys_.Platform = GPU;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.Platform", platform_str );
    }
    init_platform ( la_sys_ );

    const std::string impl_str = params_["LinearAlgebra"]["Implementation"].get<std::string>( );
    if ( impl_str == "Naive" )
    {
        la_impl_ = NAIVE;
    }
    else if ( impl_str == "BLAS" )
    {
        la_impl_ = BLAS;
    }
    else if ( impl_str == "MKL" )
    {
        la_impl_ = MKL;
    }
    else if ( impl_str == "OPENMP" )
    {
        la_impl_ = OPENMP;
    }
    else if ( impl_str == "SCALAR" )
    {
        la_impl_ = SCALAR;
    }
    else if ( impl_str == "SCALAR_TEX" )
    {
        la_impl_ = SCALAR_TEX;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.Implementation", impl_str );
    }

    const std::string matrix_str = params_["LinearAlgebra"]["MatrixFormat"].get<std::string>( );
    if ( matrix_str == "CSR" )
    {
        la_matrix_format_ = CSR;
    }
    else if ( matrix_str == "COO" )
    {
        la_matrix_format_ = COO;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.MatrixFormat", impl_str );
    }
}

// function to evaluate given exact solution u = (z-y, x-z, y-x) and
// p = x + 2y + 3z (divergence free, do nothing boundary condition does not hold
// hence complete Dirichlet boundary must be used.
// Used for debugging purpose to exclude errors in the mesh and dof numbering
// Use this as start solution for stationary case, if no bug, residual of newton
// method should fulfill stop criteria from the beginning.

void ChannelBenchmark::eval_exact_sol ( )
{
    ExactSol3D exact_sol;
    std::vector<double> sol_tmp;
    std::vector<int> sol_dofs;

    for ( mesh::EntityIterator cell = mesh_->begin ( DIMENSION ); cell != mesh_->end ( DIMENSION ); ++cell )
    {
        std::vector<doffem::DofID> dofs_on_cell;
        std::vector<Coord> coords_on_cell;

        for ( int var = 0; var < DIMENSION + 1; ++var )
        {
            space_.GetDofIndices ( var, mesh_->get_entity ( DIMENSION, cell->index ( ) ), &dofs_on_cell );
            space_.dof ( ).get_coord_on_cell ( var, cell->index ( ), coords_on_cell );
            for ( int i = 0; i < static_cast < int > ( dofs_on_cell.size ( ) ); ++i )
            {
                // check if dof is owned by subdomain
                if ( space_.dof ( ).is_dof_on_sd ( dofs_on_cell[i] ) )
                {
                    // save value of dof and dof number
                    sol_tmp.push_back ( exact_sol ( static_cast < Vec<DIMENSION, double> > ( coords_on_cell[i] ), var ) );
                    sol_dofs.push_back ( dofs_on_cell[i] );
                }
            }
        }
    }
    // set values to sol_ as startvalues
    sol_.SetValues ( vec2ptr ( sol_dofs ), sol_dofs.size ( ), vec2ptr ( sol_tmp ) );
}
