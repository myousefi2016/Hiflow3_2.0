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

#include "conv_diff.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>

/// Structure for the construction of dirichlet boundaries

struct ConvDiffDirichletBC
{

    ConvDiffDirichletBC ( int var, std::vector<int> dirichlet_bdy, std::vector<double> dirichlet_val, int T_var )
    : var_ ( var ),
    t_var_ ( T_var ), // quadrant of the temperature
    d_bdy_ ( dirichlet_bdy ),
    d_val_ ( dirichlet_val )
    {
        assert ( ( var_ == 0 ) );
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;
        if ( var_ != t_var_ )
        {
            return values;
        }

        const int material_num = face.get_material_number ( );

        for ( int l = 0; l < d_bdy_.size ( ); ++l )
        {
            if ( material_num == d_bdy_[l] )
            {
                std::cout << " MATERIAL " << material_num << std::endl;
                values.resize ( coords_on_face.size ( ), d_val_[l] );
                return values;
            }
        }

        return values;
    }

    std::vector<int> d_bdy_;
    std::vector<double> d_val_;
    const int var_;
    int t_var_;
};

/// ***************************************************************************
/// INITIALIZATION
/// ***************************************************************************

MetFlowConvDiffApp::MetFlowConvDiffApp ( const std::string& root_path, const std::string& param_filename, const std::string& base_param_filename, bool resumed_run )
: MetFlowApp ( root_path, param_filename, base_param_filename, resumed_run ),
computed_T_mass_ ( false ),
computed_T_mass_dual_ ( false ),

#ifdef RUN_MODE2
T_linear_solver_dual_ ( new FGMRES<LAD>( ) ),
#endif
T_linear_solver_ ( new FGMRES<LAD>( ) )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "Create MetFlowConvDiffApp" << std::endl;

    MetFlowApp::local_asm_primal_ = &this->T_local_asm_;
    MetFlowApp::local_asm_dual_ = &this->T_local_asm_dual_;
    MetFlowApp::local_asm_est_ = &this->T_local_asm_est_;
    MetFlowApp::local_asm_prod_ = &this->T_local_asm_prod_;
    MetFlowApp::linear_solver_ = this->T_linear_solver_;
    MetFlowApp::linear_solver_dual_ = this->T_linear_solver_dual_;

    if ( COSYSTEM == 1 )
        MetFlowApp::cyl_coord_ = true;
    else
        MetFlowApp::cyl_coord_ = false;

    MetFlowApp::filename_nonlinear_solver_ = "/log/NonlinearSolver";
    MetFlowApp::filename_linear_solver_ = "/log/LinearSolver";

    // let the compiler select wheter HDF5 should be used
#ifdef WITH_HDF5
    use_hdf5_ = true;
    if ( rank ( ) == master_rank ( ) ) std::cout << "HDF5 I/O activated within MetFlowBoussinesqTEHD" << std::endl;
#endif

    this->vel_var_.resize ( DIM, 0 );
    this->vel_var_[1] = 1;
    this->vel_var_[DIM - 1] = DIM - 1;

    this->t_var_ = 0;
    this->num_vars_ = 1;
    this->num_eq_ = 1;
    this->is_linear_problem_ = true;
}

MetFlowConvDiffApp::~MetFlowConvDiffApp ( )
{
}

// prepare application

void MetFlowConvDiffApp::initial_prepare ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare ConvDiff problem " << std::endl;
    Timer timer;
    timer.start ( );

    MetFlowApp::prepare_periodicity ( );

    // Build initial mesh
#ifndef PARALLEL_PARTITIONING
    MetFlowApp::build_initial_mesh ( this->adapt_counter_ );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );
#else
    MetFlowApp::build_initial_mesh_parallel ( this->adapt_counter_ );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );
#endif

    MetFlowApp::compute_comm_pattern ( "/log/CommPattern.txt" );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // read in parameters
    MetFlowConvDiffApp::prepare_parameters ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // **********************************************
    // Primal problem
    // **********************************************

    // setup FE space
    MetFlowConvDiffApp::prepare_space ( *this->space_, this->coupling_vars_, this->mesh_, 1 );
    MetFlowConvDiffApp::prepare_space_convection ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // setup peridoidc BC
    MetFlowApp::periodify_space ( *this->space_, this->mesh_ );
    MetFlowApp::periodify_space ( this->V_space_, this->mesh_ );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // Setup LA
    MetFlowConvDiffApp::prepare_lin_alg_structures_convection ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // Setup Time discretization
    MetFlowApp::build_initial_time_mesh ( this->adapt_counter_ );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    MetFlowApp::prepare_time_method ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // Setup Assmebler
    MetFlowConvDiffApp::prepare_assembler ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // Compute characteristic quantities
    MetFlowConvDiffApp::compute_char_quant ( );

    // **********************************************
    // Dual problem
    // **********************************************

#ifdef RUN_MODE2
    // read in parameters
    MetFlowConvDiffApp::prepare_goal_functional ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // setup FE space
    MetFlowConvDiffApp::prepare_space ( *this->space_dual_, this->coupling_vars_dual_, this->mesh_, -1 );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // setup peridoidc BC
    MetFlowApp::periodify_space ( *this->space_dual_, this->mesh_ );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // Setup Time discretization
    MetFlowApp::prepare_time_method_dual ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // Setup Assmebler
    MetFlowConvDiffApp::prepare_assembler_dual ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );
#endif
}

void MetFlowConvDiffApp::dwr_loop_prepare ( )
{
    MetFlowApp::dwr_loop_prepare ( );
    T_local_asm_.set_goal_functional ( this->j_ );
    T_local_asm_dual_.set_goal_functional ( this->j_ );
    T_local_asm_est_.set_goal_functional ( this->j_ );
}

// Prepare application-specific parameters

void MetFlowConvDiffApp::prepare_parameters ( )
{
    MetFlowApp::prepare_parameters ( );

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare ConvDiff parameters\n";

    kappa_ = params_["Fluid"]["ThermalDiffusivity"].get<double>( );
    gamma_ = params_["Fluid"]["ThermalConvection"].get<double>( );
    lambda_ = params_["Fluid"]["ThermalReaction"].get<double>( );
    start_T_ = params_["ExperimentSetup"]["StartTemperature"].get<double>( );

    std::string base_bdy = "DirichletBdy";
    std::string base_val = "DirichletVal";
    dirichlet_val_.clear ( );
    dirichlet_bdy_.clear ( );
    for ( int l = 1; l <= 4; ++l )
    {
        int bdy = params_["ExperimentSetup"]["BC"][base_bdy + static_cast < ostringstream* > ( &( ostringstream ( ) << l ) )->str ( )].get<int>( );
        double val = params_["ExperimentSetup"]["BC"][base_val + static_cast < ostringstream* > ( &( ostringstream ( ) << l ) )->str ( )].get<double>( );

        if ( bdy >= 0 )
        {
            dirichlet_bdy_.push_back ( bdy );
            dirichlet_val_.push_back ( val );
        }
    }

    // output parameters for debugging
    LOG_INFO ( "parameters", params_ );
}

/// Prepare FE space

void MetFlowConvDiffApp::prepare_space ( VectorSpace<double>& space, std::vector< std::vector< bool> >& coupling_vars, MeshPtr mesh, int mode )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare FE space" << std::endl;

    // setup finite element ansatz
    const int velocity_deg = params_["FiniteElements"]["VelocityDegree"] .get<int>( );
    const int temperature_deg = params_["FiniteElements"]["TemperatureDegree"].get<int>( );

    std::vector<int> degrees ( this->num_vars_, velocity_deg );
    std::vector<bool> is_cg ( this->num_vars_, true );

    degrees[t_var_] = temperature_deg;

    int high_deg = velocity_deg + 2 * temperature_deg;
    int quad_order = params_["FiniteElements"]["QuadratureOrder"].get<int>( );

    // set up finite element space

    std::string ordering = this->base_params_["LinearAlgebra"]["DOFOrdering"].get<std::string>( );
    if ( ordering == "Classic" )
        space.Init ( degrees, *mesh, is_cg, hiflow::doffem::HIFLOW_CLASSIC );
    if ( ordering == "Cuthill" )
        space.Init ( degrees, *mesh, is_cg, hiflow::doffem::CUTHILL_MCKEE );
    if ( ordering == "King" )
        space.Init ( degrees, *mesh, is_cg, hiflow::doffem::KING );

    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "  Total number of dofs:     " << space.dof ( ).ndofs_global ( ) << std::endl;
        std::cout << "  Number of variables:      " << space.get_nb_var ( ) << std::endl;
        std::cout << "  FE degree of velocity:    " << velocity_deg << std::endl;
        std::cout << "  FE degree of temperature: " << temperature_deg << std::endl;
        std::cout << "  Highest polnomial degree: " << high_deg << std::endl;
        std::cout << "  Quadrature order:         " << quad_order << std::endl;
    }

    // quadrature
    QuadratureSelection qsel ( quad_order );
    this->global_asm_.set_quadrature_selection_function ( qsel );
    if ( params_["FiniteElements"]["QuadratureOrder"].get<int>( ) < high_deg )
    {
        if ( rank ( ) == master_rank ( ) ) std::cout << "Quadrature Order too low? " << std::endl;
    }

    // Variable couplings vel , press, temp, pot,(aug_press)
    coupling_vars.clear ( );
    coupling_vars.resize ( this->num_vars_ );
    coupling_vars[0].push_back ( true );
}

void MetFlowConvDiffApp::prepare_space_convection ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare convection space" << std::endl;

    // setup finite element ansatz
    const int velocity_deg = params_["FiniteElements"]["VelocityDegree"] .get<int>( );
    const int temperature_deg = params_["FiniteElements"]["TemperatureDegree"].get<int>( );

    std::vector<int> degrees ( this->num_vars_, velocity_deg );
    std::vector<bool> is_cg ( this->num_vars_, true );

    degrees[t_var_] = temperature_deg;

    int high_deg = velocity_deg + 2 * temperature_deg;
    int quad_order = params_["FiniteElements"]["QuadratureOrder"].get<int>( );

    // setup space for convection field
    std::vector<int> conv_degrees ( DIM, velocity_deg );
    std::vector<bool> conv_is_cg ( DIM, true );

    std::string ordering = this->base_params_["LinearAlgebra"]["DOFOrdering"].get<std::string>( );
    if ( ordering == "Classic" )
        V_space_.Init ( conv_degrees, *mesh_, conv_is_cg, hiflow::doffem::HIFLOW_CLASSIC );
    if ( ordering == "Cuthill" )
        V_space_.Init ( conv_degrees, *mesh_, conv_is_cg, hiflow::doffem::CUTHILL_MCKEE );
    if ( ordering == "King" )
        V_space_.Init ( conv_degrees, *mesh_, conv_is_cg, hiflow::doffem::KING );

    this->V_coupling_vars_.resize ( DIM );
    for ( int l = 0; l < DIM; ++l )
    {
        this->V_coupling_vars_[l].resize ( DIM, true );
    }
}

void MetFlowConvDiffApp::clear_problem ( )
{
    // clear general stuff
    MetFlowApp::clear_problem ( );

    if ( rank ( ) == master_rank ( ) )
        std::cout << "> Clear ConvDiff Problem" << std::endl;

    // clear problem specific stuff
    this->V_coupling_vars_.clear ( );
    this->V_space_.Clear ( );
    this->V_couplings_.Clear ( );
    this->V_.Clear ( );
    this->V_prev_.Clear ( );
    this->V_next_.Clear ( );

    this->computed_T_mass_ = false;
    this->computed_T_mass_dual_ = false;
    T_mass_matrix_.Clear ( );
    T_mass_matrix_dual_.Clear ( );

    if ( COSYSTEM == 1 )
        MetFlowApp::cyl_coord_ = true;
    else
        MetFlowApp::cyl_coord_ = false;

    MetFlowApp::filename_nonlinear_solver_ = "/log/NonlinearSolver";
    MetFlowApp::filename_linear_solver_ = "/log/LinearSolver";

#ifdef WITH_HDF5
    use_hdf5_ = true;
#endif

    this->vel_var_.resize ( DIM, 0 );
    this->vel_var_[1] = 1;
    this->vel_var_[DIM - 1] = DIM - 1;

    this->t_var_ = 0;
    this->num_vars_ = 1;
    this->num_eq_ = 1;
    this->is_linear_problem_ = true;
}

void MetFlowConvDiffApp::prepare_lin_alg_structures_convection ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare LA structures for primal problem " << std::endl;

    SparsityStructure sparsity;
    this->global_asm_.compute_sparsity_structure ( V_space_, sparsity, &this->V_coupling_vars_ );

    // Initialize linear algebra structures
    V_couplings_.Init ( comm_, V_space_.dof ( ) );
    V_couplings_.InitializeCouplings ( sparsity.off_diagonal_rows, sparsity.off_diagonal_cols );

    // Initialize vectors
#ifdef USE_HYPRE
    V_.Init ( comm_, V_couplings_ );
    V_prev_.Init ( comm_, V_couplings_ );
    V_next_.Init ( comm_, V_couplings_ );
#else
    V_.Init ( comm_, V_couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    V_.InitStructure ( );
    V_prev_.Init ( comm_, V_couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    V_prev_.InitStructure ( );
    V_next_.Init ( comm_, V_couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    V_next_.InitStructure ( );
#endif

    V_.Zeros ( );
    V_prev_.Zeros ( );
    V_next_.Zeros ( );
}

// prepare local assembler

void MetFlowConvDiffApp::prepare_assembler ( )
{
    MetFlowApp::prepare_assembler ( );

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare local primal TEHD assembler " << std::endl;

    T_local_asm_.set_kappa ( kappa_ );
    T_local_asm_.set_lambda ( lambda_ );
    T_local_asm_.set_gamma ( gamma_ );

    // set perturbation
    int temp_type = params_["Perturbation"]["TempType"].get<int>( );
    int perturb_primal = this->base_params_["Perturbation"]["PerturbEF"].get<int>( );

    if ( perturb_primal > 0 )
    {
        if ( temp_type == 0 ) T_local_asm_.set_perturb_type ( this->t_var_, true, false );
        if ( temp_type == 1 ) T_local_asm_.set_perturb_type ( this->t_var_, false, true );
        if ( temp_type == 2 ) T_local_asm_.set_perturb_type ( this->t_var_, true, true );
    }

    /*
    T_local_asm_.set_conv        (V_);
    T_local_asm_.set_conv_prev   (V_prev_);
    T_local_asm_.set_conv_next   (V_next_);
    T_local_asm_.set_conv_space  (V_space_);
     */
    T_local_asm_.setup_fe_evaluators ( );

    double width = params_["ExperimentSetup"]["LocalSource"]["Width"].get<double>( );
    double ampl = params_["ExperimentSetup"]["LocalSource"]["Amplitude"].get<double>( );
    double x = params_["ExperimentSetup"]["LocalSource"]["XPosition"].get<double>( );
    double y = params_["ExperimentSetup"]["LocalSource"]["YPosition"].get<double>( );
    double z = params_["ExperimentSetup"]["LocalSource"]["ZPosition"].get<double>( );
    double t0 = params_["ExperimentSetup"]["LocalSource"]["StartAt"].get<double>( );
    double t1 = params_["ExperimentSetup"]["LocalSource"]["StopAt"].get<double>( );

    Vec<DIM, double> pos;
    pos[0] = x;
    if ( DIM > 1 )
        pos[1] = y;
    if ( DIM > 2 )
        pos[2] = z;

    LocalSourceTerm<DIM, double>* source_struct = new LocalSourceTerm<DIM, double> ( t_var_, pos, width, ampl, t0, t1 );
    T_local_asm_.set_source_term ( source_struct );
}

// prepare local assembler
// TODO weitere parameter notwendig? z.B. stabilisierung?

void MetFlowConvDiffApp::prepare_assembler_dual ( )
{
    MetFlowApp::prepare_assembler_dual ( );

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare local dual TEHD assembler " << std::endl;

    T_local_asm_dual_.set_kappa ( kappa_ );
    T_local_asm_dual_.set_lambda ( lambda_ );
    T_local_asm_dual_.set_gamma ( gamma_ );

    /*
    T_local_asm_dual_.set_conv        (V_);
    T_local_asm_dual_.set_conv_prev   (V_prev_);
    T_local_asm_dual_.set_conv_next   (V_next_);
    T_local_asm_dual_.set_conv_space  (V_space_);
     */
    T_local_asm_dual_.setup_fe_evaluators ( );

}

// prepare local assembler

void MetFlowConvDiffApp::prepare_assembler_est ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare local estimator assembler " << std::endl;

    // call general setup routines
    MetFlowApp::prepare_assembler_est ( );

    std::vector<int> order_p ( this->num_vars_, 1 );
    std::vector<int> order_d ( this->num_vars_, 1 );

    if ( this->method_dual_ == "GalerkinDC" )
    {
        order_d[0] = 0;
    }
    if ( this->is_stationary_ )
    {
        order_p[0] = 0;
        order_d[0] = 0;
    }
    T_local_asm_est_.set_time_fem_order ( order_p, order_d );

    T_local_asm_est_.set_kappa ( kappa_ );
    T_local_asm_est_.set_lambda ( lambda_ );
    T_local_asm_est_.set_gamma ( gamma_ );

    /*
    T_local_asm_est_.set_conv        (V_);
    T_local_asm_est_.set_conv_prev   (V_prev_);
    T_local_asm_est_.set_conv_next   (V_next_);
    T_local_asm_est_.set_conv_space  (V_space_);
     */
    //T_local_asm_est_.setup_fe_evaluators();

    double width = params_["ExperimentSetup"]["LocalSource"]["Width"].get<double>( );
    double ampl = params_["ExperimentSetup"]["LocalSource"]["Amplitude"].get<double>( );
    double x = params_["ExperimentSetup"]["LocalSource"]["XPosition"].get<double>( );
    double y = params_["ExperimentSetup"]["LocalSource"]["YPosition"].get<double>( );
    double z = params_["ExperimentSetup"]["LocalSource"]["ZPosition"].get<double>( );
    double t0 = params_["ExperimentSetup"]["LocalSource"]["StartAt"].get<double>( );
    double t1 = params_["ExperimentSetup"]["LocalSource"]["StopAt"].get<double>( );
    Vec<DIM, double> pos;
    pos[0] = x;
    if ( DIM > 1 )
        pos[1] = y;
    if ( DIM > 2 )
        pos[2] = z;

    LocalSourceTerm<DIM, double>* source_struct = new LocalSourceTerm<DIM, double> ( t_var_, pos, width, ampl, t0, t1 );
    T_local_asm_est_.set_source_term ( source_struct );
}

/// prepare dirichlet boundary conditions

void MetFlowConvDiffApp::prepare_bc ( )
{
    this->prepare_bc ( 0., this->cur_time_ );
}

/// prepare dirichlet boundary conditions

void MetFlowConvDiffApp::prepare_bc ( double rotation, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare ConvDiff boundary conditions at time " << time << std::endl;
    this->dirichlet_dofs_.clear ( );
    this->dirichlet_values_.clear ( );

    ConvDiffDirichletBC bc[1] = { ConvDiffDirichletBC ( 0, this->dirichlet_bdy_, this->dirichlet_val_, this->t_var_ ) };

    for ( int var = 0; var < this->num_vars_; ++var )
    {
        compute_dirichlet_dofs_and_values ( bc[var],
                                            *this->space_,
                                            var,
                                            this->dirichlet_dofs_,
                                            this->dirichlet_values_ );
    }
}

/// prepare dirichlet boundary conditions for dual problem -> same structure as primal problem, but homogenous

void MetFlowConvDiffApp::prepare_bc_dual ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare ConvDiff dual boundary conditions" << std::endl;
    this->dirichlet_dofs_dual_.clear ( );
    this->dirichlet_values_dual_.clear ( );
    std::vector<double> zeros ( dirichlet_val_.size ( ), 0. );

    ConvDiffDirichletBC bc[1] = { ConvDiffDirichletBC ( 0, this->dirichlet_bdy_, zeros, this->t_var_ ) };

    for ( int var = 0; var < this->num_vars_; ++var )
    {
        compute_dirichlet_dofs_and_values ( bc[var],
                                            *this->space_dual_,
                                            var,
                                            this->dirichlet_dofs_dual_,
                                            this->dirichlet_values_dual_ );
    }
}

/// prepare goal functional

void MetFlowConvDiffApp::prepare_goal_functional ( )
{
    std::string type = params_["Adaptivity"]["GoalFunctional"]["Type"].get<std::string>( );
    int tmp_ic = this->base_params_["Adaptivity"]["GoalFunctional"]["ActiveFinalTime"].get<int>( );
    int tmp_force = this->base_params_["Adaptivity"]["GoalFunctional"]["ActiveInterval"].get<int>( );

    bool active_force = false;
    bool active_ic = false;

    if ( tmp_ic != 0 )
        active_ic = true;

    if ( tmp_force != 0 )
        active_force = true;

    double y_min = this->base_params_["Adaptivity"]["GoalFunctional"]["yMin"].get<double>( );
    double y_max = this->base_params_["Adaptivity"]["GoalFunctional"]["yMax"].get<double>( );
    double x_min = this->base_params_["Adaptivity"]["GoalFunctional"]["xMin"].get<double>( );
    double x_max = this->base_params_["Adaptivity"]["GoalFunctional"]["xMax"].get<double>( );
    double z_min = this->base_params_["Adaptivity"]["GoalFunctional"]["zMin"].get<double>( );
    double z_max = this->base_params_["Adaptivity"]["GoalFunctional"]["zMax"].get<double>( );
    double t_min = this->base_params_["Adaptivity"]["GoalFunctional"]["tMin"].get<double>( );
    double t_max = this->base_params_["Adaptivity"]["GoalFunctional"]["tMax"].get<double>( );
    double scale = this->base_params_["Adaptivity"]["GoalFunctional"]["Scale"].get<double>( );

    if ( type == "VarOnSub" )
    {
        j_var_on_sub_.set_variable ( this->t_var_ );
        j_ = &j_var_on_sub_;
    }
    else
    {
        if ( rank ( ) == master_rank ( ) ) std::cout << "> Type of goal functional not defined!" << std::endl;
        quit_program ( );
    }
    j_->set_active_parts ( active_force, active_ic );
    j_->set_cart_box ( t_min, t_max, x_min, x_max, y_min, y_max, z_min, z_max );
    j_->set_scale_factor ( scale );

    T_local_asm_.set_goal_functional ( this->j_ );
    T_local_asm_dual_.set_goal_functional ( this->j_ );
    T_local_asm_est_.set_goal_functional ( this->j_ );
}
