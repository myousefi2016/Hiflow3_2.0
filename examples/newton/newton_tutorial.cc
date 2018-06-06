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

#include "newton_tutorial.h"

#include <iomanip>
#include <iterator>
#include <algorithm>
#include <vector>

namespace
{
    static const char* DATADIR = MESHES_DATADIR;
    static const int MASTER_RANK = 0;
    static const char* PARAM_FILENAME = "newton_tutorial.xml";

    static bool CONSOLE_OUTPUT_ACTIVE = true;
    static const int CONSOLE_THRESHOLD_LEVEL = 3;
}
#define CONSOLE_OUTPUT(lvl, x) { if (CONSOLE_OUTPUT_ACTIVE &&            \
                                    lvl <= CONSOLE_THRESHOLD_LEVEL) {   \
            for (int i = 0; i < lvl; ++i) { std::cout << "  "; }        \
            std::cout << x << "\n"; }}

// time measurement of the different calculation steps

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

class NewtonTutorial : public NonlinearProblem<LAD>
{
  public:

    NewtonTutorial ( const std::string& param_filename )
    : comm_ ( MPI_COMM_WORLD ),
    params_ ( param_filename.c_str ( ), MASTER_RANK, MPI_COMM_WORLD ),
    refinement_level_ ( 0 ),
    is_done_ ( false )

    {
    }

    virtual ~NewtonTutorial ( )
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
        LogKeeper::get_log ( "info" ).set_target ( &info_log );
        std::ofstream debug_log ( ( simul_name_ + "_debug_log" ).c_str ( ) );
        LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

        CONSOLE_OUTPUT ( 0, "============================================================" );
        CONSOLE_OUTPUT ( 0, "==== NewtonTutorial                                    ===" );
        CONSOLE_OUTPUT ( 0, "====    built using HiFlow3.                             ===" );
        CONSOLE_OUTPUT ( 0, "====                                                     ===" );
        CONSOLE_OUTPUT ( 0, "==== Engineering Mathematics and Computing Lab (EMCL)    ===" );
        CONSOLE_OUTPUT ( 0, "============================================================" );
        CONSOLE_OUTPUT ( 0, "" );

        // Output parameters for debugging
        LOG_INFO ( "parameters", params_ );

        // Setup timing report
        TimingScope::set_report ( &time_report_ );

        {
            TimingScope tscope ( "Setup" );
            setup_linear_algebra ( );
            read_mesh ( );
            prepare ( );
        }

        LOG_INFO ( "simulation", "Solving stationary problem" );

        // Measurement of the global calculation time
        Timer timer;

        newton_.Solve ( &sol_ );
        CONSOLE_OUTPUT ( 1, "Newton ended with residual norm " << newton_.GetResidual ( )
                         << " after " << newton_.iter ( ) << " iterations." );

        timer.stop ( );
        CONSOLE_OUTPUT ( 0, "" );
        CONSOLE_OUTPUT ( 1, "Measured time of interest " << static_cast < int > ( timer.get_duration ( ) / 60 ) << "m"
                         << ( timer.get_duration ( ) - ( static_cast < int > ( timer.get_duration ( ) / 60 ) )*60 ) << "s" );

        visualize ( );

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

    // Visualize the solution in a file.
    void visualize ( );

    // Helper functions for nonlinear solver
    virtual void EvalFunc ( const LAD::VectorType& in, LAD::VectorType* out );
    void compute_residual ( const LAD::VectorType& in, LAD::VectorType* out ); // updates res_ with the residual
    void compute_stationary_residual ( const LAD::VectorType& in, LAD::VectorType* out ); // residual computation in stationary mode

    virtual void EvalGrad ( const LAD::VectorType& in, LAD::MatrixType* out );
    void compute_jacobian ( const LAD::VectorType& in, LAD::MatrixType* out ); // updates matrix_ with the jacobian matrix
    void compute_stationary_matrix ( const LAD::VectorType& in, LAD::MatrixType* out ); // jacobi matrix computation in stationary mode

    // Linear algebra set up
    void setup_linear_algebra ( );

    // MPI stuff
    const MPI_Comm comm_;
    int rank_, num_partitions_;

    // Linear algebra stuff
    SYSTEM la_sys_;
    IMPLEMENTATION la_impl_;
    MATRIX_FORMAT la_matrix_format_;

    // Parameter data read in from file.
    PropertyTree params_;
    std::string simul_name_; // parameter 'OutputPrefix': prefix for output files

    // Flow model variables
    double Um_, H_, W_, rho_, nu_;
    double U_mean_, diam_; // (DFG-specific) mean velocity and diameter of cylinder.

    // Meshes
    MeshPtr master_mesh_, mesh_;
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
    std::vector<double> residual_history_, forcing_term_history_;
    std::string forcing_strategy_;

    // forcing strategy parameters
    double eta_initial;
    double eta_max;
    double gamma_EW2;
    double alpha_EW2;

    // nonlinear solver
    Newton<LAD> newton_;

    StandardGlobalAssembler<double> global_asm_;

    bool is_done_, convergence_test_;

    std::vector<int> dirichlet_dofs_;
    std::vector<Scalar> dirichlet_values_;

#ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> ilupp_;
    bool use_ilupp_;
#endif

    HierarchicalReport<TimingData> time_report_;
    Table results_table_;
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
        NewtonTutorial app ( param_filename );
        app.run ( );
    }
    catch ( const std::exception& e )
    {
        std::cerr << "\nProgram ended with uncaught exception.\n";
        std::cerr << e.what ( ) << "\n";
        return -1;
    }
    MPI_Finalize ( );

    return 0;
}

void NewtonTutorial::read_mesh ( )
{
    TimingScope tscope ( "read_mesh" );
    if ( rank ( ) == MASTER_RANK )
    {
        const std::string mesh_name = params_["Mesh"]["Filename"].get<std::string>( );
        std::string mesh_filename = std::string ( DATADIR ) + mesh_name;

        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0 );

        CONSOLE_OUTPUT ( 1, "Read mesh with " << master_mesh_->num_entities ( DIMENSION ) << " cells." );

        refinement_level_ = 0;
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );

        for ( int r = 0; r < initial_ref_lvl; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
            ++refinement_level_;
        }
        LOG_INFO ( "mesh", "Initial refinement level = " << refinement_level_ );
        CONSOLE_OUTPUT ( 1, "Refined mesh (level " << refinement_level_ << ") has "
                         << master_mesh_->num_entities ( DIMENSION ) << " cells." );
    }

    MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_ );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );

    std::ostringstream rank_str;
    rank_str << rank ( );

    PVtkWriter writer ( comm_ );
    std::string output_file = std::string ( "mesh_local.pvtu" );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

void NewtonTutorial::prepare ( )
{
    TimingScope tscope ( "prepare" );

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

    space_.Init ( degrees, *mesh_ );

    CONSOLE_OUTPUT ( 1, "Total number of dofs = " << space_.dof ( ).ndofs_global ( ) );

    for ( int p = 0; p < num_partitions ( ); ++p )
    {
        CONSOLE_OUTPUT ( 2, "Num dofs on process " << p << " = " << space_.dof ( ).ndofs_on_sd ( p ) );
    }

    // prepare linear algebra structures
    couplings_.Clear ( );
    couplings_.Init ( communicator ( ), space_.dof ( ) );

    // prepare global assembler
    QuadratureSelection q_sel ( params_["QuadratureOrder"].get<int>( ) );
    global_asm_.set_quadrature_selection_function ( q_sel );

    // compute matrix graph
    SparsityStructure sparsity;

    global_asm_.compute_sparsity_structure ( space_, sparsity );

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

    // get forcing strategy parameters from param file
    forcing_strategy_ = params_["NonlinearSolver"]["ForcingStrategy"].get<std::string>( );
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

    // setup forcing strategy object within the nonlinear solver
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

void NewtonTutorial::prepare_bc ( )
{
    TimingScope tscope ( "prepare_bc" );

    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    const int inflow_bdy = params_["Boundary"]["InflowMaterial"].get<int>( );
    const int outflow_bdy = params_["Boundary"]["OutflowMaterial"].get<int>( );

    ChannelFlowBC2d bc[2] = { ChannelFlowBC2d ( 0, H_, Um_, inflow_bdy, outflow_bdy ),
                             ChannelFlowBC2d ( 1, H_, Um_, inflow_bdy, outflow_bdy ) };

    for ( int var = 0; var < DIMENSION; ++var )
    {
        compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                            dirichlet_dofs_, dirichlet_values_ );
    }

    // apply BC to initial solution
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }
}

void NewtonTutorial::visualize ( )
{
    TimingScope tscope ( "Visualization" );

    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    std::stringstream input;

    input << simul_name_ << "_solution";

    input << "_stationary";

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

    // Create visualization vector from post-processing vector.
    sol_.UpdateCouplings ( );

    // Visualize the solution vector and all attributes
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 0 ), "u" );
#if DIMENSION >= 2
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 1 ), "v" );
#endif
#if DIMENION == 3
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 2 ), "w" );
#endif
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, DIMENSION ), "p" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );
}

// Member function of the nonlinear solver newton_
// Used to launch the residual vector computation at each newton iteration

void NewtonTutorial::EvalFunc ( const LAD::VectorType& in,
                                LAD::VectorType* out )
{
    compute_residual ( in, out );
}

// Computes residual vector F(sol)

void NewtonTutorial::compute_residual ( const LAD::VectorType& in,
                                        LAD::VectorType* out )
{
    TimingScope tscope ( "Compute Residual" );

    compute_stationary_residual ( in, out );
    // correct BC -- set Dirichlet dofs to 0
    if ( !dirichlet_dofs_.empty ( ) )
    {
        std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
        out->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                         vec2ptr ( zeros ) );
    }
}

// Computes the residual vector on the stationary mode

void NewtonTutorial::compute_stationary_residual ( const LAD::VectorType& in,
                                                   LAD::VectorType* out )
{
    StationaryFlowAssembler local_asm ( in, nu_, rho_ );
    global_asm_.assemble_vector ( space_, local_asm, *out );
}

// Member function of the nonlinear solver newton_
// Used to launch the jacobian matrix computation at each newton iteration

void NewtonTutorial::EvalGrad ( const LAD::VectorType& in,
                                LAD::MatrixType* out )
{
    compute_jacobian ( in, out );
}

// Returns the jacobian matrix J of nonlinear problem F at point x

void NewtonTutorial::compute_jacobian ( const LAD::VectorType& in,
                                        LAD::MatrixType* out )
{
    {
        TimingScope tscope ( "Compute Jacobian" );

        compute_stationary_matrix ( in, out );
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

// Computes the jacobian matrix J on the stationary mode

void NewtonTutorial::compute_stationary_matrix ( const LAD::VectorType& in,
                                                 LAD::MatrixType* out )
{
    StationaryFlowAssembler local_asm ( in, nu_, rho_ );
    global_asm_.assemble_matrix ( space_, local_asm, *out );
}

void NewtonTutorial::setup_linear_algebra ( )
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
