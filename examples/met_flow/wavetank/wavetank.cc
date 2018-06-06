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

#include "wavetank.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>

/// Structure for the construction of dirichlet boundaries

struct DirichletBC3D
{

    DirichletBC3D ( int var, double Omega, double warm_T, double cold_T, int warm_bdy, int cold_bdy, int T_var )
    : var_ ( var ),
    T_var_ ( T_var ), // quadrant of the temperature
    omega_ ( Omega ), // in rad/s
    w_T_ ( warm_T ),
    w_bdy_ ( warm_bdy ),
    c_bdy_ ( cold_bdy ),
    c_T_ ( cold_T )
    {
        assert ( ( var_ >= 0 ) && ( var_ <= 4 ) );
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );

        const bool cold = ( material_num == c_bdy_ );
        const bool warm = ( material_num == w_bdy_ );
        const bool bottom = ( material_num == 13 );
        const bool top = ( material_num == 23 );
        bool frame = ( warm || cold || bottom );

        // velocity
        if ( frame && ( var_ < DIM ) )
        {
            values.resize ( coords_on_face.size ( ), 0. );
#ifndef ROTATING_FOR
            if ( var_ == 0 )
            {
                // rotation in phi direction
                for ( int i = 0; i < coords_on_face.size ( ); ++i )
                {
                    const Coord& pt = coords_on_face[i];
                    values[i] = omega_ * pt[1];
                }
            }
#endif
            return values;
        }

        // no outflow in z-direction
        if ( top && ( var_ == 2 ) )
        {
            values.resize ( coords_on_face.size ( ), 0. );
            return values;
        }

        // temperature
        if ( warm && ( var_ == T_var_ ) )
        {
            values.resize ( coords_on_face.size ( ), w_T_ );
            return values;
        }
        else if ( cold && ( var_ == T_var_ ) )
        {
            values.resize ( coords_on_face.size ( ), c_T_ );
            return values;
        }
    }

    const double w_T_;
    const double c_T_;
    const int var_;
    const int c_bdy_;
    const int w_bdy_;
    int T_var_;
    double omega_;
};

MetFlowBousCyl3dApp::MetFlowBousCyl3dApp ( const std::string& root_path, const std::string& param_filename, const std::string base_param_filename, bool resume_run )
: MetFlowApp ( root_path, param_filename, base_param_filename, resume_run ),
T_var_ ( DIM + 1 ),
temp_supg_mode_ ( 0 ),
#ifdef SCHUR_SOLVER
PVT_linear_solver_ ( new FGMRES<LAD>( ) )
#else
PVT_linear_solver_ ( new GMRES<LAD>( ) )
#endif
{
    std::cout << "Create MetFlowBous3DApp" << std::endl;
    MetFlowApp::local_asm_primal_ = &this->PVT_local_asm_;
    MetFlowApp::linear_solver_ = this->PVT_linear_solver_;

    MetFlowApp::cyl_coord_ = true;

    gravity_.resize ( DIM, 0.0 );

    this->num_vars_ = DIM + 2;
    // let the compiler select wheter HDF5 should be used
#ifdef WITH_HDF5
    use_hdf5_ = true;
    if ( rank ( ) == master_rank ( ) ) std::cout << "HDF5 I/O activated within MetFlowBoussinesq" << std::endl;
#endif
}

MetFlowBousCyl3dApp::~MetFlowBousCyl3dApp ( )
{

}

// prepare application
// prepare application

void MetFlowBousCyl3dApp::initial_prepare ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare Wavetank problem " << std::endl;
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
    MetFlowBousCyl3dApp::prepare_parameters ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // **********************************************
    // Primal problem
    // **********************************************

    // setup FE space
    MetFlowBousCyl3dApp::prepare_space ( *this->space_, this->coupling_vars_, this->mesh_, 1 );

    // setup peridoidc BC
    MetFlowBousCyl3dApp::periodify_space ( *this->space_, this->mesh_ );

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
    MetFlowBousCyl3dApp::prepare_assembler ( );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  " << timer.get_duration ( ) << " sec" << std::endl;
    timer.reset ( );
    timer.start ( );

    // Compute characteristic quantities
    MetFlowBousCyl3dApp::compute_char_quant ( );
}

void MetFlowBousCyl3dApp::dwr_loop_prepare ( )
{
    MetFlowApp::dwr_loop_prepare ( );
}

void MetFlowBousCyl3dApp::prepare_assembler ( )
{
    MetFlowApp::prepare_assembler ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare local assembler " << std::endl;

    PVT_local_asm_.set_rho ( this->rho_ );
    PVT_local_asm_.set_nu ( this->nu_ );
    PVT_local_asm_.set_kappa ( this->kappa_ );
    PVT_local_asm_.set_gravity ( gravity_ );
    PVT_local_asm_.set_alpha_g ( alpha_g_ );
    PVT_local_asm_.set_ref_T ( ref_T_ );
    PVT_local_asm_.set_start_T ( start_T_ );
    PVT_local_asm_.set_omega ( omega_ );

    if ( temp_supg_mode_ == 0 ) PVT_local_asm_.set_temp_supg_mode_to_OFF ( );
    if ( temp_supg_mode_ == 1 ) PVT_local_asm_.set_temp_supg_mode_to_STD ( );
    PVT_local_asm_.set_temp_supg ( temp_supg_gamma_, 0., 0., 0. );

    if ( graddiv_mode_ == 0 ) PVT_local_asm_.set_graddiv_mode_to_OFF ( );
    if ( graddiv_mode_ == 1 ) PVT_local_asm_.set_graddiv_mode_to_CONST ( );
    if ( graddiv_mode_ == 2 ) PVT_local_asm_.set_graddiv_mode_to_VMS ( );
    if ( graddiv_mode_ == 3 ) PVT_local_asm_.set_graddiv_mode_to_SUPG ( );
    PVT_local_asm_.set_graddiv ( graddiv_gamma_, graddiv_C_ );

    if ( skew_mode_ == 0 ) PVT_local_asm_.set_skew_mode_to_OFF ( );
    if ( skew_mode_ == 1 ) PVT_local_asm_.set_skew_mode_to_ON ( );

    PVT_local_asm_.set_nusselt_area ( nusselt_r_min_, nusselt_r_max_, nusselt_surface_id_ );

    if ( rank ( ) == master_rank ( ) ) this->PVT_local_asm_.print_parameters ( );
}

/// Prepare application-specific parameters

void MetFlowBousCyl3dApp::prepare_parameters ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "\n> Prepare parameters\n";

    rad_i_ = params_["Mesh"]["InnerRadius"].get<double>( );
    rad_o_ = params_["Mesh"]["OuterRadius"].get<double>( );
    height_ = params_["Mesh"]["Height"].get<double>( );

    rho_ = params_["Fluid"]["Density"].get<double>( );
    nu_ = params_["Fluid"]["Viscosity"].get<double>( );
    kappa_ = params_["Fluid"]["ThermalDiffusivity"].get<double>( );
    alpha_g_ = params_["Fluid"]["ThermalExpansionCoefficient"].get<double>( );

    gravity_[0] = params_["ExperimentSetup"]["GravityX"].get<double>( );
    gravity_[1] = params_["ExperimentSetup"]["GravityY"].get<double>( );
    gravity_[2] = params_["ExperimentSetup"]["GravityZ"].get<double>( );

    ref_T_ = params_["ExperimentSetup"]["ReferenceTemperature"].get<double>( );
    start_T_ = params_["ExperimentSetup"]["StartTemperature"].get<double>( );
    warm_T_ = params_["ExperimentSetup"]["WarmTemperature"].get<double>( );
    cold_T_ = params_["ExperimentSetup"]["ColdTemperature"].get<double>( );
    omega_ = params_["ExperimentSetup"]["RotationSpeed"].get<double>( ); // in rpm
    omega_ = ( omega_ / 60. )*2. * pi; // in rad/s = rps*2pi -> later multiplication with radius

    warm_material_ = params_["ExperimentSetup"]["WarmBoundary"].get<int>( );
    cold_material_ = params_["ExperimentSetup"]["ColdBoundary"].get<int>( );

    temp_supg_mode_ = params_["Stabilization"]["TempSUPG"]["Mode"].get<int>( );
    temp_supg_gamma_ = params_["Stabilization"]["TempSUPG"]["Gamma"].get<double>( );

    nusselt_r_min_ = rad_i_ + params_["PostProcessing"]["NusseltNumber"]["InnerRadius"].get<double>( ) * ( rad_o_ - rad_i_ );
    nusselt_r_max_ = rad_i_ + params_["PostProcessing"]["NusseltNumber"]["OuterRadius"].get<double>( ) * ( rad_o_ - rad_i_ );
    nusselt_surface_id_ = params_["PostProcessing"]["NusseltNumber"]["Surface"].get<int>( );
}

void MetFlowBousCyl3dApp::prepare_ic ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare initial conditions " << std::endl;

    int perturb_ic = this->base_params_["Perturbation"]["PerturbIC"].get<int>( 0 );

    std::string starting_type = this->base_params_["InitialCondition"]["Type"].get<std::string>( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  StartingType: " << starting_type.c_str ( ) << std::endl;

    // enable the counting of jacobian assemblies for efficient preconditioning
    this->ilupp_time_step_ = 0;

    // **********************************************
    // Load IC
    // **********************************************
    if ( starting_type == "Load" )
    {
        // Load start solution and return
        const int velocity_deg = params_["FiniteElements"]["VelocityDegree"] .get<int>( );
        const int pressure_deg = params_["FiniteElements"]["PressureDegree"] .get<int>( );
        const int temperature_deg = params_["FiniteElements"]["TemperatureDegree"].get<int>( );

        std::string prefix = this->base_params_["InitialCondition"]["FilePrefix"].get<std::string>( );
        std::stringstream ss;
        std::string mesh_name = this->base_params_["Mesh"]["Filename"].get<std::string>( );
        std::string fluid_name = params_["Fluid"]["Name"].get<std::string>( );
        mesh_name.erase ( 0, 3 );
        int ref_level = this->base_params_["Mesh"]["InitialRefLevel"].get<int>( );
        std::string filename_start = this->root_ + "/" + prefix + "-" + mesh_name + "-" + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << ref_level ) )->str ( ) + "-"
                + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->num_partitions_ ) )->str ( ) + "-"
                + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << velocity_deg ) )->str ( )
                + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << pressure_deg ) )->str ( )
                + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << temperature_deg ) )->str ( ) + "-"
#ifdef AUGMENT_PRESS
                + "AUG-"
#endif
                + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->ref_T_ ) )->str ( ) + "-"
                + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->cold_T_ ) )->str ( ) + "-"
                + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->warm_T_ ) )->str ( ) + "-"
                + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->omega_ ) )->str ( ) + "-"
                + fluid_name
                + ".h5";
        std::string groupname_start = "start";
        std::string prefix_start = "start";
        ss << prefix_start;

        double cur_time;
        this->read_file ( *this->solP_, filename_start, groupname_start, ss.str ( ) );

        this->base_->CloneFrom ( *this->solP_ );

        // perturb initial solution
        if ( perturb_ic > 0 )
        {
            this->perturb_ic ( );
        }
        return;
    }

    // Compute start solution
    double initial_rotation_speed = omega_;
    if ( params_["InitialCondition"].contains ( "InitialRotationSpeed" ) )
    {
        initial_rotation_speed = params_["InitialCondition"]["InitialRotationSpeed"].get<double>( );
        initial_rotation_speed = ( initial_rotation_speed / 60. )*2. * pi;
    }

    // ************************************************************************
    // set initial velocity field

    if ( rank ( ) == master_rank ( ) ) std::cout << " > Set up velocity and temperature field" << std::endl;

    for ( int var = 0; var < DIM + 2; ++var )
    {
        const int tdim = this->space_->mesh ( ).tdim ( );
        for ( mesh::EntityIterator it = this->space_->mesh ( ).begin ( tdim ), end_it = this->space_->mesh ( ).end ( tdim ); it != end_it; ++it )
        {
            std::vector<int> global_dof_ids;
            this->space_->GetDofIndices ( var, *it, &global_dof_ids );
            int num_dofs = global_dof_ids.size ( );
            std::vector<double> values;
            values.resize ( num_dofs, 0. );

#ifndef ROTATING_FOR
            if ( var == 0 )
            {
                std::vector< Coord > coords;
                this->space_->dof ( ).get_coord_on_cell ( var, it->index ( ), coords );
                for ( int i = 0; i < num_dofs; i++ )
                {
                    if ( rank_ == this->space_->dof ( ).owner_of_dof ( global_dof_ids.at ( i ) ) )
                    {
                        values[i] = initial_rotation_speed * coords[i][1];
                        this->solP_->SetValues ( &global_dof_ids.at ( i ), 1, &values.at ( i ) );
                        this->solP_prev_->SetValues ( &global_dof_ids.at ( i ), 1, &values.at ( i ) );
                    }
                }
            }
#endif
            if ( var == T_var_ )
            {
                std::vector< Coord > coords;
                this->space_->dof ( ).get_coord_on_cell ( var, it->index ( ), coords );
                for ( int i = 0; i < num_dofs; i++ )
                {
                    if ( rank_ == this->space_->dof ( ).owner_of_dof ( global_dof_ids.at ( i ) ) )
                    {
                        this->solP_->SetValues ( &global_dof_ids.at ( i ), 1, &start_T_ );
                        this->solP_prev_->SetValues ( &global_dof_ids.at ( i ), 1, &start_T_ );
                    }
                }
            }
        }
    }

    // reset BC for stationary run
    if ( rank ( ) == master_rank ( ) ) std::cout << " > Preparing Dirichlet dofs for stationary problem (initial_rotation_speed=" << initial_rotation_speed << ")" << std::endl;

#ifdef ROTATING_FOR
    PVT_local_asm_.set_omega ( initial_rotation_speed );
    MetFlowBousCyl3dApp::prepare_bc ( 0.0, 0. );
#else
    MetFlowBousCyl3dApp::prepare_bc ( initial_rotation_speed, 0. );
#endif

    // apply BC to initial solution
    if ( !this->dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        this->solP_->SetValues ( vec2ptr ( this->dirichlet_dofs_ ), this->dirichlet_dofs_.size ( ), vec2ptr ( this->dirichlet_values_ ) );
        this->solP_prev_->SetValues ( vec2ptr ( this->dirichlet_dofs_ ), this->dirichlet_dofs_.size ( ), vec2ptr ( this->dirichlet_values_ ) );
    }
    this->solP_->Update ( );
    this->solP_prev_->Update ( );

    if ( starting_type == "Zero" )
    {
        this->base_->CloneFrom ( *this->solP_ );
        // perturb initial solution
        if ( perturb_ic > 0 )
        {
            this->perturb_ic ( );
        }
        return;
    }
    // ************************************************************************
    // Solve stationary problem

    if ( rank ( ) == master_rank ( ) ) std::cout << " > Solve stationary problem" << std::endl;

    // prepare local assembler
    PVT_local_asm_.set_time_discretization_off ( );

    // number of nudge steps, 0 -> no nudging
    int num_nudge_steps = this->base_params_["PrimalNonlinearSolver"]["Stationary"]["NumberNudgeSteps"].get<int>( );

    // viscosities are lowered by this 1/factor in each nudge step
    double factor = this->base_params_["PrimalNonlinearSolver"]["Stationary"]["NudgeFactor"].get<double>( 10. );

    // viscosities of stationary problem are larger by this factor compared to the instationary parameters (given in the param)
    double stationary_instationary_factor = this->base_params_["PrimalNonlinearSolver"]["Stationary"]["StationaryInstationaryNudgeFactor"].get<double>( 1. );

    MetFlowApp::prepare_nls ( true );
    for ( int i = 0; i <= num_nudge_steps; ++i )
    {
        double effective_thermal_diffusion = this->kappa_ * pow ( factor, ( num_nudge_steps - i ) ) * stationary_instationary_factor;
        double effective_nu = this->nu_ * pow ( factor, ( num_nudge_steps - i ) ) * stationary_instationary_factor;
        PVT_local_asm_.set_kappa ( effective_thermal_diffusion );
        PVT_local_asm_.set_nu ( effective_nu );

        if ( rank ( ) == master_rank ( ) )
        {
            if ( num_nudge_steps > 0 )
            {
                std::cout << "Nudging:" << std::endl;
                std::cout << "-> effective thermal diffusion: " << effective_thermal_diffusion << std::endl;
                std::cout << "-> effective nu:                " << effective_nu << std::endl;
            }
            else
            {
                std::cout << " > Nudging deactivated." << std::endl;
            }
        }

        PVT_local_asm_.set_solP_prev ( *this->solP_prev_ );

        // solve
        bool success;
        success = solve_nlp ( 1 );
        this->solP_prev_->CloneFrom ( *this->solP_ );

        this->visualize_solution ( *this->solP_, this->space_, "out/nudged_start_solutions_step_", i );
    }

    // set correct parameters
    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << " > Preparing the instationary problem:" << std::endl;
        std::cout << "   thermal diffusion: " << this->kappa_ << std::endl;
        std::cout << "   nu:                " << this->nu_ << std::endl;
    }

    PVT_local_asm_.set_kappa ( this->kappa_ );
    PVT_local_asm_.set_nu ( this->nu_ );

    // reset local assembler
    if ( method_ == "GalerkinCD" ) PVT_local_asm_.set_time_discretization_to_galerkin_cd ( mod_galerkin_ );
    if ( method_ == "Theta" ) PVT_local_asm_.set_time_discretization_to_theta ( this->theta_ );

    // reset BC for instationary run
    if ( rank ( ) == master_rank ( ) ) std::cout << " > Preparing Dirichlet dofs for instationary problem (rotation_speed=" << omega_ << ")" << std::endl;
    this->prepare_bc ( );
#ifdef ROTATING_FOR
    PVT_local_asm_.set_omega ( omega_ );
#endif

    this->base_->CloneFrom ( *this->solP_ );

    // add perturbation
    if ( perturb_ic > 0 )
    {
        this->perturb_ic ( );
    }

    // enable the counting of jacobian assemblies for efficient preconditioning
    this->ilupp_time_step_ = 0;
}

void MetFlowBousCyl3dApp::prepare_linear_solver ( bool ic )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare primal linear solver" << std::endl;
    const int max_iter = this->base_params_["PrimalLinearSolver"]["MaximumIterations"].get<int>( );
    const double abs_tol = this->base_params_["PrimalLinearSolver"]["AbsoluteTolerance"].get<double>( );
    const double rel_tol = this->base_params_["PrimalLinearSolver"]["RelativeTolerance"].get<double>( );
    const double div_tol = this->base_params_["PrimalLinearSolver"]["DivergenceLimit"].get<double>( );
    const int basis_size = this->base_params_["PrimalLinearSolver"]["BasisSize"].get<int>( );

    PVT_linear_solver_->InitControl ( max_iter, abs_tol, rel_tol, div_tol );

#ifdef WITH_ILUPP

    PVT_linear_solver_->InitParameter ( basis_size, "RightPreconditioning" );

    // prepare ILUPP
    const int prepro_type = this->base_params_["ILUPP"]["PreprocessingType"].get<int>( );
    const int precond_no = this->base_params_["ILUPP"]["PreconditionerNumber"].get<int>( );
    const int max_levels = this->base_params_["ILUPP"]["MaxMultilevels"].get<int>( );
    const double mem_factor = this->base_params_["ILUPP"]["MemFactor"].get<double>( );
    const double threshold = this->base_params_["ILUPP"]["PivotThreshold"].get<double>( );
    const double min_pivot = this->base_params_["ILUPP"]["MinPivot"].get<double>( );

    PVT_precond_.InitParameter ( prepro_type, precond_no, max_levels, mem_factor, threshold, min_pivot );
    PVT_linear_solver_->SetupPreconditioner ( this->PVT_precond_ );
#else
    PVT_linear_solver_->InitParameter ( basis_size, "NoPreconditioning" );
#endif
    // set the matrix to be used as the operator
    PVT_linear_solver_->SetupOperator ( *this->matrix_ );

    PVT_linear_solver_->SetReuseBasis ( false );
}

void MetFlowBousCyl3dApp::perturb_ic ( )
{

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Perturb initial condition " << std::endl;

    int vel_var = params_["Perturbation"]["InitialCondition"]["Velocity"]["Component"].get<int>( );
    double vel_ampl = params_["Perturbation"]["InitialCondition"]["Velocity"]["Amplitude"].get<double>( );
    double vel_offset = params_["Perturbation"]["InitialCondition"]["Velocity"]["Offset"].get<double>( );
    double vel_wn_phi = params_["Perturbation"]["InitialCondition"]["Velocity"]["WaveNumberPhi"].get<double>( );
    double vel_wn_rad = params_["Perturbation"]["InitialCondition"]["Velocity"]["WaveNumberRad"].get<double>( );
    double vel_wn_z = params_["Perturbation"]["InitialCondition"]["Velocity"]["WaveNumberZ"].get<double>( );

    double temp_ampl = params_["Perturbation"]["InitialCondition"]["Temperature"]["Amplitude"].get<double>( );
    double temp_offset = params_["Perturbation"]["InitialCondition"]["Temperature"]["Offset"].get<double>( );
    double temp_wn_phi = params_["Perturbation"]["InitialCondition"]["Temperature"]["WaveNumberPhi"].get<double>( );
    double temp_wn_rad = params_["Perturbation"]["InitialCondition"]["Temperature"]["WaveNumberRad"].get<double>( );
    double temp_wn_z = params_["Perturbation"]["InitialCondition"]["Temperature"]["WaveNumberZ"].get<double>( );

    double r_min = rad_i_ + params_["Perturbation"]["InitialCondition"]["BoundingBox"]["rMin"].get<double>( ) * ( rad_o_ - rad_i_ );
    double r_max = rad_i_ + params_["Perturbation"]["InitialCondition"]["BoundingBox"]["rMax"].get<double>( ) * ( rad_o_ - rad_i_ );
    double z_min = params_["Perturbation"]["InitialCondition"]["BoundingBox"]["zMin"].get<double>( ) * height_;
    double z_max = params_["Perturbation"]["InitialCondition"]["BoundingBox"]["zMax"].get<double>( ) * height_;
    double phi_min = params_["Perturbation"]["InitialCondition"]["BoundingBox"]["phiMin"].get<double>( ) * pi;
    double phi_max = params_["Perturbation"]["InitialCondition"]["BoundingBox"]["phiMax"].get<double>( ) * pi;

    double d_r = r_max - r_min;
    double d_z = z_max - z_min;
    double d_phi = phi_max - phi_min;

    double perturb_scale = this->base_params_["Perturbation"]["InitialCondition"]["Scale"].get<double>( 1. );

    std::string type = this->base_params_["Perturbation"]["InitialCondition"]["Type"].get<std::string>( );

    if ( type == "Create" )
    {
        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "  Velocity ------------ " << std::endl;
            std::cout << "  Component:            " << vel_var << std::endl;
            std::cout << "  Amplitude:            " << vel_ampl << std::endl;
            std::cout << "  Azimuthal wavenumber: " << vel_wn_phi << std::endl;
            std::cout << "  Radial wavenumber:    " << vel_wn_rad << std::endl;
            std::cout << "  Axial wavenumber:     " << vel_wn_z << std::endl;
            std::cout << "  Temperature---------- " << std::endl;
            std::cout << "  Amplitude:            " << temp_ampl << std::endl;
            std::cout << "  Azimuthal wavenumber: " << temp_wn_phi << std::endl;
            std::cout << "  Radial wavenumber:    " << temp_wn_rad << std::endl;
            std::cout << "  Axial wavenumber:     " << temp_wn_z << std::endl;
        }

        this->perturb_->Zeros ( );

        std::vector<int> vars;
        std::vector<double> ampl;
        std::vector<double> wn_phi;
        std::vector<double> wn_rad;
        std::vector<double> wn_z;
        std::vector<double> offset;
        vars.resize ( 2, 0 );
        ampl.resize ( 2, 0. );
        wn_phi.resize ( 2, 0. );
        wn_rad.resize ( 2, 0. );
        wn_z.resize ( 2, 0. );
        offset.resize ( 2, 0. );

        vars[0] = vel_var;
        vars[1] = DIM + 1;
        ampl[0] = vel_ampl;
        ampl[1] = temp_ampl;
        wn_phi[0] = vel_wn_phi;
        wn_phi[1] = temp_wn_phi;
        wn_rad[0] = vel_wn_rad;
        wn_rad[1] = temp_wn_rad;
        wn_z[0] = vel_wn_z;
        wn_z[1] = temp_wn_z;
        offset[0] = vel_offset;
        offset[1] = temp_offset;
        int tdim = DIM;

        for ( mesh::EntityIterator it = this->space_->mesh ( ).begin ( tdim ), end_it = this->space_->mesh ( ).end ( tdim ); it != end_it; ++it )
        {
            for ( int v = 0; v < vars.size ( ); ++v )
            {
                int var = vars[v];
                std::vector<int> global_dof_ids;
                this->space_->GetDofIndices ( var, *it, &global_dof_ids );

                int num_dofs = global_dof_ids.size ( );
                std::vector<double> values;
                values.resize ( num_dofs, 0. );

                std::vector<Coord> coords;
                this->space_->dof ( ).get_coord_on_cell ( var, it->index ( ), coords );
                for ( int i = 0; i < num_dofs; i++ )
                {
                    if ( rank_ == this->space_->dof ( ).owner_of_dof ( global_dof_ids.at ( i ) ) )
                    {
                        double phi = coords.at ( i )[0];
                        double r = coords.at ( i )[1];
                        double z = coords.at ( i )[2];

                        if ( r < r_min || r > r_max )
                            continue;
                        if ( z < z_min || z > z_max )
                            continue;
                        if ( phi < phi_min || phi > phi_max )
                            continue;

                        phi = ( phi - phi_min ) / d_phi;
                        r = ( r - r_min ) / d_r;
                        z = ( z - z_min ) / d_z;

                        double perturbation = ampl[v]
                                * ( offset[v]
                                + cos ( wn_phi[v] * 2.0 * pi * phi )
                                * cos ( wn_z[v] * 2. * pi * z )
                                * cos ( wn_rad[v] * 2. * pi * r ) );

                        // ensure continuity
                        perturbation *= sin ( pi * z ) * sin ( pi * r );

                        this->perturb_->Add ( global_dof_ids.at ( i ), perturbation );
                    }
                }
            }
        }

        this->perturb_->Scale ( perturb_scale );
        this->perturb_->Update ( );
        // store
        std::string filename_perturb = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsIn"].get<std::string>( );
        std::string groupname_perturb = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsGroup"].get<std::string>( );
        std::string prefix_perturb = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsPrefix"].get<std::string>( );
        std::stringstream ss_perturb;
        ss_perturb << prefix_perturb << 0;
        this->write_file ( *this->perturb_, filename_perturb, groupname_perturb, ss_perturb.str ( ) );
    }
    if ( type == "Load" )
    {
        // load
        std::string filename_perturb = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsIn"].get<std::string>( );
        std::string groupname_perturb = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsGroup"].get<std::string>( );
        std::string prefix_perturb = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsPrefix"].get<std::string>( );
        std::stringstream ss_perturb;
        ss_perturb << prefix_perturb << 0;
        MetFlowApp::read_file ( *this->perturb_, filename_perturb, groupname_perturb, ss_perturb.str ( ) );
    }

    if ( type == "Dual" )
    {
        // load dual
        std::string filename_perturb = this->base_params_["Perturbation"]["Test"]["DualSolution"]["SnapshotsIn"].get<std::string>( );
        std::string groupname_perturb = this->base_params_["Perturbation"]["Test"]["DualSolution"]["SnapshotsGroup"].get<std::string>( );
        std::string prefix_perturb = this->base_params_["Perturbation"]["Test"]["DualSolution"]["SnapshotsPrefix"].get<std::string>( );
        std::stringstream ss_perturb;
        ss_perturb << prefix_perturb << 0;
        MetFlowApp::read_file ( *this->perturb_, filename_perturb, groupname_perturb, ss_perturb.str ( ) );
        this->perturb_->Scale ( perturb_scale );
        this->perturb_->Update ( );
    }

    this->solP_->Axpy ( *this->perturb_, -1. );

    this->solP_->Update ( );
}

void MetFlowBousCyl3dApp::prepare_space ( VectorSpace<double>& space, std::vector< std::vector< bool> >& coupling_vars, MeshPtr mesh, int mode )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare coupled space" << std::endl;

    // setup finite element ansatz
    const int velocity_deg = params_["FiniteElements"]["VelocityDegree"] .get<int>( );
    const int pressure_deg = params_["FiniteElements"]["PressureDegree"] .get<int>( );
    const int temperature_deg = params_["FiniteElements"]["TemperatureDegree"].get<int>( );
    int add_press_var = 0;
    int quad_order = params_["FiniteElements"]["QuadratureOrder"].get<int>( );

    std::vector<int> degrees ( DIM + 2, velocity_deg );
    std::vector<bool> is_cg ( DIM + 2, true );
    assert ( DIM + 1 == T_var_ );
    degrees[DIM] = pressure_deg;
    degrees[DIM + 1] = temperature_deg;

#ifdef AUGMENT_PRESS
    degrees.push_back ( 0 );
    is_cg.push_back ( false );
    add_press_var = 1;
#endif

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
#ifdef AUGMENT_PRESS
        std::cout << "  FE degree of pressure:    " << pressure_deg << " + " << degrees[DIM + 2] << std::endl;
#else
        std::cout << "  FE degree of pressure:    " << pressure_deg << std::endl;
#endif
        std::cout << "  FE degree of temperature: " << temperature_deg << std::endl;
        std::cout << "  Quadrature order:         " << quad_order << std::endl;
    }
    // quadrature
    QuadratureSelection qsel ( quad_order );
    this->global_asm_.set_quadrature_selection_function ( qsel );

    // Variable couplings vel , press, temp,(aug_press)
    coupling_vars.resize ( DIM + 2 + add_press_var );

    //    x x x x x (x)
    //  x x x x x (x)
    //  x x x x x (x)
    //  x x x 0 0 (0)
    //  x x x 0 x (0)
    // (x x x 0 0  0)

    // momentum equation
    for ( int i = 0; i < DIM; ++i )
    {
        for ( int j = 0; j < DIM + 2 + add_press_var; ++j )
        {
            coupling_vars[i].push_back ( true );
        }
    }

    // incompressibility constraint
    for ( int i = 0; i < DIM; ++i )
        coupling_vars[DIM].push_back ( true );
    for ( int i = DIM; i < DIM + 2 + add_press_var; ++i )
        coupling_vars[DIM].push_back ( false );
#ifdef AUGMENT_PRESS
    for ( int i = 0; i < DIM; ++i )
        coupling_vars[DIM + 2].push_back ( true );
    for ( int i = DIM; i < DIM + 2 + add_press_var; ++i )
        coupling_vars[DIM + 2].push_back ( false );
#endif

    // heat equation
    for ( int i = 0; i < DIM; ++i )
        coupling_vars[DIM + 1].push_back ( true );

    coupling_vars[DIM + 1].push_back ( false );
    coupling_vars[DIM + 1].push_back ( true );
#ifdef AUGMENT_PRESS
    coupling_vars[DIM + 1].push_back ( false );
#endif
}

void MetFlowBousCyl3dApp::update_assembler ( )
{
}

void MetFlowBousCyl3dApp::filter_solution ( )
{
}

void MetFlowBousCyl3dApp::create_log_files ( int offset )
{
    this->create_log_file ( offset, "/log/FlowCharacteristics", "txt" );
    this->create_log_file ( offset, "/log/NonlinearSolver", "txt" );
    this->create_log_file ( offset, "/log/LinearSolver", "txt" );
    this->create_log_file ( offset, "/log/Norm2", "txt" );
}

/// prepare dirichlet boundary conditions

void MetFlowBousCyl3dApp::prepare_bc ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare coupled boundary conditions" << std::endl;
    // use omega_ as rotation speed in boundary conditions
    this->prepare_bc ( omega_, this->cur_time_ );
}

/// prepare dirichlet boundary conditions

void MetFlowBousCyl3dApp::prepare_bc ( double rotation, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare boundary conditions in cylindrical coordinates" << std::endl;
    this->dirichlet_dofs_.clear ( );
    this->dirichlet_values_.clear ( );

    DirichletBC3D bc[5] = { DirichletBC3D ( 0, rotation, warm_T_, cold_T_, warm_material_, cold_material_, T_var_ ),
                           DirichletBC3D ( 1, rotation, warm_T_, cold_T_, warm_material_, cold_material_, T_var_ ),
                           DirichletBC3D ( 2, rotation, warm_T_, cold_T_, warm_material_, cold_material_, T_var_ ),
                           DirichletBC3D ( 3, rotation, warm_T_, cold_T_, warm_material_, cold_material_, T_var_ ),
                           DirichletBC3D ( 4, rotation, warm_T_, cold_T_, warm_material_, cold_material_, T_var_ ) };

    for ( int var = 0; var < DIM + 2; ++var )
    {
        if ( var < DIM || var == T_var_ )
            compute_dirichlet_dofs_and_values ( bc[var],
                                                *this->space_,
                                                var,
                                                this->dirichlet_dofs_,
                                                this->dirichlet_values_ );
    }

    if ( rank ( ) == master_rank ( ) )
        std::cout << "   Number of Dirichlet dofs total: " << this->dirichlet_dofs_.size ( ) << std::endl;
}

void MetFlowBousCyl3dApp::update_preconditioner ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
#ifdef SCHUR_SOLVER
    this->update_schur_preconditioner ( u, *this->solP_prev_ );
#else
#    ifdef WITH_ILUPP
    this->PVT_precond_.SetupOperator ( *DF );
#    endif
#endif
}

void MetFlowBousCyl3dApp::write_ic ( )
{
    const int velocity_deg = params_["FiniteElements"]["VelocityDegree"] .get<int>( );
    const int pressure_deg = params_["FiniteElements"]["PressureDegree"] .get<int>( );
    const int temperature_deg = params_["FiniteElements"]["TemperatureDegree"].get<int>( );

    std::stringstream ss;
    std::string prefix = this->base_params_["InitialCondition"]["FilePrefix"].get<std::string>( );
    std::string mesh_name = this->base_params_["Mesh"]["Filename"].get<std::string>( );
    std::string fluid_name = params_["Fluid"]["Name"].get<std::string>( );
    mesh_name.erase ( 0, 3 );
    int ref_level = this->base_params_["Mesh"]["InitialRefLevel"].get<int>( );
    std::string filename_start = this->root_ + "/" + prefix + "-" + mesh_name + "-" + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << ref_level ) )->str ( ) + "-"
            + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->num_partitions_ ) )->str ( ) + "-"
            + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << velocity_deg ) )->str ( )
            + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << pressure_deg ) )->str ( )
            + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << temperature_deg ) )->str ( ) + "-"
#ifdef AUGMENT_PRESS
            + "AUG-"
#endif
            + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->ref_T_ ) )->str ( ) + "-"
            + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->cold_T_ ) )->str ( ) + "-"
            + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->warm_T_ ) )->str ( ) + "-"
            + static_cast < std::ostringstream* > ( &( std::ostringstream ( ) << this->omega_ ) )->str ( ) + "-"
            + fluid_name
            + ".h5";
    std::string groupname_start = "start";
    std::string prefix_start = "start";
    ss << prefix_start;

    double time = 0.;
    MetFlowApp::write_file ( *this->solP_prev_, filename_start, groupname_start, ss.str ( ) );
}

void MetFlowBousCyl3dApp::visualize_solution ( int step )
{
    std::string prefix = "out/boussinesq";
    this->visualize_solution ( *this->solP_, this->space_, prefix, step );
}

void MetFlowBousCyl3dApp::visualize_solution ( LAD::VectorType& sol, const VectorSpace<double>* space, std::string const& prefix, int time_step )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Visualize solution in cylindrical coordinates at time step " << time_step << std::endl;
    std::stringstream input;
    if ( time_step < 10 )
        input << prefix << ".000" << time_step;
    else if ( time_step < 100 )
        input << prefix << ".00" << time_step;
    else if ( time_step < 1000 )
        input << prefix << ".0" << time_step;
    else
        input << prefix << "." << time_step;

    if ( num_partitions_ > 1 )
        input << ".pvtu";
    else
        input << ".vtu";

    std::string visu_filename = input.str ( );

    CylCellVisualization visu ( *this->space_, 1, comm_, master_rank_ );

    this->solP_->Update ( );
    visu.visualize ( EvalCylFeFunction ( *this->space_, sol, 0, rank_ ), "vel", DIM );
    visu.visualize ( EvalFeFunction<LAD>( *this->space_, sol, DIM ), "P", 1 );
    visu.visualize ( EvalFeFunction<LAD>( *this->space_, sol, DIM + 1 ), "T", 1 );

    // parallel statistics
    if ( false )
    {
        std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
        std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

        AttributePtr sub = mesh_->get_attribute ( "_sub_domain_", mesh_->tdim ( ) );
        AttributePtr remote = mesh_->get_attribute ( "_remote_index_", mesh_->tdim ( ) );
        for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
              it != mesh_->end ( mesh_->tdim ( ) );
              ++it )
        {
            remote_index.at ( it->index ( ) ) = remote->get_int_value ( it->index ( ) );
            sub_domain.at ( it->index ( ) ) = sub->get_int_value ( it->index ( ) );
        }
        visu.visualize_cell_data ( remote_index, "_remote_index_" );
        visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    }

    // write file
    visu.write ( visu_filename );
}

void MetFlowBousCyl3dApp::post_processing ( )
{
    int compute_quant = params_["PostProcessing"]["FlowQuantities"].get<int>( );
    int compute_norm = params_["PostProcessing"]["Norm"].get<int>( );

    if ( compute_quant == 1 )
        MetFlowBousCyl3dApp::compute_flow_quantities ( "/log/FlowCharacteristics.txt", this->solP_, this->cur_time_ );
    if ( compute_norm == 1 )
        MetFlowBousCyl3dApp::compute_solution_norm ( "/log/Norm2.txt", this->solP_, this->solP_prev_, this->cur_time_ );
}

void MetFlowBousCyl3dApp::compute_quant_vector ( std::string filename, double time )
{
    // compute integrated quantities
    int incomp_scalars = 24;
    int bous_scalars = 5;
    int num_scalars = incomp_scalars + bous_scalars;
    std::vector<double> local_quant;
    std::vector<double> global_quant;
    local_quant.resize ( num_scalars, 0. );
    global_quant.resize ( num_scalars, 0. );

    PVT_local_asm_.set_solP ( *this->solP_ );
    PVT_local_asm_.set_solP_prev ( *this->solP_prev_ );
    PVT_local_asm_.set_vector_mode_to_quant ( incomp_scalars, bous_scalars );
    global_asm_.integrate_multiple_scalar ( *this->space_, boost::ref ( PVT_local_asm_ ), num_scalars, local_quant );
    PVT_local_asm_.set_vector_mode_to_std ( );

    for ( int l = 0; l < num_scalars; ++l )
        MPI_Reduce ( &local_quant[l], &global_quant[l], 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );

    // take square roots if necessary
    global_quant[16] = std::sqrt ( global_quant[16] );
    global_quant[17] = std::sqrt ( global_quant[17] );
    global_quant[18] = std::sqrt ( global_quant[18] );
    global_quant[19] = std::sqrt ( global_quant[19] );
    global_quant[20] = std::sqrt ( global_quant[20] );
    global_quant[21] = std::sqrt ( global_quant[21] );
    global_quant[22] = std::sqrt ( global_quant[22] );
    global_quant[23] = std::sqrt ( global_quant[23] );
    global_quant[25] = std::sqrt ( global_quant[25] );
    global_quant[26] = std::sqrt ( global_quant[26] );
    global_quant[28] = std::sqrt ( global_quant[28] );

    // assign corresponding labels
    std::vector< std::string > labels;
    labels.resize ( num_scalars );
    std::vector< int > types;
    types.resize ( num_scalars, 0 );

    labels[0] = "KinEnergy       ";
    types[0] = 0;
    labels[1] = "azimKinEnergy   ";
    types[1] = 1;
    labels[2] = "radialKinEnergy ";
    types[2] = 2;
    labels[3] = "axialKinEnergy  ";
    types[3] = 3;
    labels[4] = "DissEnergy      ";
    types[4] = 0;
    labels[5] = "azimDissEnergy  ";
    types[5] = 1;
    labels[6] = "radialDissEnergy";
    types[6] = 2;
    labels[7] = "axialDissEnergy ";
    types[7] = 3;
    labels[8] = "azimPressure    ";
    types[8] = 1;
    labels[9] = "radialPressure  ";
    types[9] = 2;
    labels[10] = "axialPressure   ";
    types[10] = 3;
    labels[11] = "azimGradDiv     ";
    types[11] = 1;
    labels[12] = "radialGradDiv   ";
    types[12] = 2;
    labels[13] = "axialGradDiv    ";
    types[13] = 3;
    labels[14] = "azi2radConv     ";
    types[14] = 1;
    labels[15] = "rad2aziConv     ";
    types[15] = 2;
    labels[16] = "L2_Divergence   ";
    types[16] = 0;
    labels[17] = "L2_azimVort     ";
    types[17] = 0;
    labels[18] = "L2_radialVort   ";
    types[18] = 0;
    labels[19] = "L2_axialVort    ";
    types[19] = 0;
    labels[20] = "L2_azimBaseDiff ";
    types[20] = 4;
    labels[21] = "L2_radBaseDiff  ";
    types[21] = 4;
    labels[22] = "L2_axiBaseDiff  ";
    types[22] = 4;
    labels[23] = "L2_pressBaseDiff";
    types[23] = 4;

    labels[24] = "axialBuoyancy   ";
    types[24] = 3;
    labels[25] = "L2_Temperature  ";
    types[25] = 0;
    labels[26] = "L2_azimTempGrad ";
    types[26] = 0;
    labels[27] = "Nusselt         ";
    types[27] = -1;
    labels[28] = "L2_tempBaseDiff ";
    types[28] = 4;

    int nusselt_id = 27;

    // compute cellwise quantities
    // Max / min of mean div(vel) in each cell
    std::vector<double> div_mean;
    div_mean.clear ( );
    PVT_local_asm_.set_scalar_mode_to_div_mean ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_ ), div_mean );

    double lmax_div_mean = -1e6;
    double lmin_div_mean = 1e6;

    for ( int i = 0; i < div_mean.size ( ); i++ )
    {
        if ( div_mean[i] > lmax_div_mean ) lmax_div_mean = div_mean[i];
        if ( div_mean[i] < lmin_div_mean ) lmin_div_mean = div_mean[i];
    }

    double max_div_mean = -1e6;
    double min_div_mean = 1e6;

    MPI_Reduce ( &lmax_div_mean, &max_div_mean, 1, MPI_DOUBLE, MPI_MAX, master_rank ( ), comm_ );
    MPI_Reduce ( &lmin_div_mean, &min_div_mean, 1, MPI_DOUBLE, MPI_MIN, master_rank ( ), comm_ );
    min_div_mean = std::abs ( min_div_mean );
    if ( min_div_mean > max_div_mean )
        max_div_mean = min_div_mean;

    global_quant.push_back ( max_div_mean );
    labels.push_back ( "MaxMeanDiv      " );
    types.push_back ( 0 );

    // Max of div(vel) in each cell
    std::vector<double> div_max;
    div_max.clear ( );
    PVT_local_asm_.set_scalar_mode_to_div_max ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_ ), div_max );

    double lmax_div = -1e6;

    for ( int i = 0; i < div_max.size ( ); i++ )
    {
        if ( div_max[i] > lmax_div ) lmax_div = div_max[i];
    }

    double max_div = -1e6;

    MPI_Reduce ( &lmax_div, &max_div, 1, MPI_DOUBLE, MPI_MAX, master_rank ( ), comm_ );

    global_quant.push_back ( max_div );
    labels.push_back ( "MaxDiv          " );
    types.push_back ( 0 );

    // Max / min of mean temp in each cell
    std::vector<double> temp_mean;
    temp_mean.clear ( );
    PVT_local_asm_.set_scalar_mode_to_temp_mean ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_ ), temp_mean );

    double lmax_temp_mean = -1e6;
    double lmin_temp_mean = 1e6;

    for ( int i = 0; i < temp_mean.size ( ); i++ )
    {
        if ( temp_mean[i] > lmax_temp_mean ) lmax_temp_mean = temp_mean[i];
        if ( temp_mean[i] < lmin_temp_mean ) lmin_temp_mean = temp_mean[i];
    }

    double max_temp_mean = -1e6;
    double min_temp_mean = 1e6;

    MPI_Reduce ( &lmax_temp_mean, &max_temp_mean, 1, MPI_DOUBLE, MPI_MAX, master_rank ( ), comm_ );
    MPI_Reduce ( &lmin_temp_mean, &min_temp_mean, 1, MPI_DOUBLE, MPI_MIN, master_rank ( ), comm_ );
    global_quant.push_back ( max_temp_mean );
    global_quant.push_back ( min_temp_mean );
    labels.push_back ( "MaxTemp         " );
    labels.push_back ( "MinTemp         " );
    types.push_back ( 0 );
    types.push_back ( 0 );

    // correct Nusselt number
    double fac = 1.;

    if ( params_["PostProcessing"]["NusseltNumber"]["Type"].get<int>( ) == 1 )
    {
        fac = 1. / ( nusselt_r_max_ - nusselt_r_min_ );
        global_quant[nusselt_id] *= fac;
    }
    else
    {
        std::vector<double> rad_heat;
        rad_heat.clear ( );
        PVT_local_asm_.set_scalar_mode_to_rad_heat_surface ( );
        pp_asm_.assemble_scalar_boundary ( *this->space_, boost::ref ( PVT_local_asm_ ), rad_heat );
        double global_rad_heat = 0.;
        double total_rad_heat = std::accumulate ( rad_heat.begin ( ), rad_heat.end ( ), 0. );
        MPI_Reduce ( &total_rad_heat, &global_rad_heat, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
        global_quant[nusselt_id] = global_rad_heat;
    }

    if ( time == 0 )
    {
        Nu0_ = global_quant[nusselt_id];
        Nu_ = 1.;
        global_quant[nusselt_id] = Nu_;
        if ( rank ( ) == master_rank ( ) )
        {
            std::ofstream out;
            std::string path = this->root_ + "/log/Nusselt0.txt";
            out.open ( path.c_str ( ), std::ios::out );

            out.precision ( 16 );
            out << std::scientific;

            out << Nu0_ << "\n ";
            out.close ( );
        }
    }
    else
    {
        std::ifstream in;
        std::string path = this->root_ + "/log/Nusselt0.txt";
        in.open ( path.c_str ( ), std::ios::in );

        in >> Nu0_;
        in.close ( );

        Nu_ = global_quant[nusselt_id] / Nu0_;
        global_quant[nusselt_id] = Nu_;
    }

    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "  " << labels[0] << ": " << global_quant[0] << std::endl
                << "  " << labels[4] << ": " << global_quant[4] << std::endl
                << "  " << labels[17] << ": " << global_quant[17] << std::endl
                << "  " << labels[24] << ": " << global_quant[24] << std::endl
                << "  " << labels[25] << ": " << global_quant[25] << std::endl
                << "  " << labels[20] << ": " << global_quant[20] << std::endl
                << "  ---------------------" << std::endl
                << "  " << labels[1] << ": " << global_quant[1] << std::endl
                << "  " << labels[5] << ": " << global_quant[5] << std::endl
                << "  " << labels[8] << ": " << global_quant[8] << std::endl
                << "  " << labels[12] << ": " << global_quant[12] << std::endl
                << "  " << labels[15] << ": " << global_quant[15] << std::endl
                //              << "  " << labels[24] << ": " << global_quant[24] << std::endl
                << "  " << labels[21] << ": " << global_quant[21] << std::endl
                << "  ---------------------" << std::endl
                << "  " << labels[2] << ": " << global_quant[2] << std::endl
                << "  " << labels[6] << ": " << global_quant[6] << std::endl
                << "  " << labels[9] << ": " << global_quant[9] << std::endl
                << "  " << labels[13] << ": " << global_quant[13] << std::endl
                << "  " << labels[16] << ": " << global_quant[16] << std::endl
                //              << "  " << labels[25] << ": " << global_quant[25] << std::endl
                << "  " << labels[22] << ": " << global_quant[22] << std::endl
                << "  ---------------------" << std::endl
                << "  " << labels[3] << ": " << global_quant[3] << std::endl
                << "  " << labels[7] << ": " << global_quant[7] << std::endl
                << "  " << labels[10] << ": " << global_quant[10] << std::endl
                << "  " << labels[14] << ": " << global_quant[14] << std::endl
                << "  " << labels[11] << ": " << global_quant[11] << std::endl
                //              << "  " << labels[26] << ": " << global_quant[26] << std::endl
                << "  " << labels[23] << ": " << global_quant[23] << std::endl
                << "  ---------------------" << std::endl
                << "  " << labels[18] << ": " << global_quant[18] << std::endl
                << "  " << labels[19] << ": " << global_quant[19] << std::endl
                << "  " << labels[26] << ": " << global_quant[26] << std::endl
                << "  " << labels[27] << ": " << global_quant[27] << std::endl;

        std::string path = this->root_ + filename;
        std::ofstream out;
        out.open ( path.c_str ( ), std::ios::out | std::ios::app );

        out.precision ( 10 );
        out << std::scientific;

        out << time << " ";
        for ( int l = 0; l < num_scalars + 4; ++l )
            out << global_quant[l] << " ";

        out << "\n";
        out.close ( );

        if ( time == 0. )
        {
            path = this->root_ + "/log/Labels.txt";
            out.open ( path.c_str ( ), std::ios::out );
            out << "Time,";
            for ( int l = 0; l < num_scalars + 3; ++l )
                out << labels[l] << ",";

            out << "\n";
            out.close ( );

            path = this->root_ + "/log/Types.txt";
            out.open ( path.c_str ( ), std::ios::out );
            out << "-1 ";
            for ( int l = 0; l < num_scalars + 3; ++l )
                out << types[l] << " ";

            out << "\n";
            out.close ( );
        }
    }
}

/// Compute some flow characteristics

void MetFlowBousCyl3dApp::compute_flow_quantities ( std::string filename, LAD::VectorType* sol, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Compute characteristics of flow " << std::endl;
    Timer timer;
    timer.start ( );

    // Set solution vector w.r.t. which quantities should be computed
    PVT_local_asm_.set_solP ( *sol );

    //  compute_quant_scalar();
    compute_quant_vector ( filename, time );

    double vol_int = MetFlowApp::compute_volume_int ( );
    double press_int = MetFlowApp::compute_pressure_int ( *sol, DIM );
#ifdef AUGMENT_PRESS
    press_int += MetFlowApp::compute_pressure_int ( *sol, this->aug_p_var_ );
#endif
    double av_press = press_int / vol_int;

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) )
        std::cout << "  took " << timer.get_duration ( ) << " sec" << std::endl;

    // Reset solution vector
    PVT_local_asm_.set_solP ( *this->solP_ );
}

/// Compute some flow characteristics

void MetFlowBousCyl3dApp::compute_solution_norm ( std::string filename, LAD::VectorType* sol, LAD::VectorType* sol_prev, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Compute norm of solution vector " << std::endl;

    Timer timer;
    timer.start ( );

    std::vector<int> vel_var ( DIM, 0 );
    for ( int i = 0; i < DIM; ++i )
        vel_var[i] = i;
    int p_var = DIM;

    // Set soultion vector w.r.t. which quantities should be computed
    PVT_local_asm_prod_.set_left_vector ( *sol );
    PVT_local_asm_prod_.set_right_vector ( *sol );
    PVT_local_asm_prod_.mark_vector_field ( vel_var );

    std::vector<bool> no_L2, L2_vel, L2_press, L2_temp;
    std::vector<bool> no_H1, H1_vel, H1_temp;
    std::vector<bool> no_H2;

    no_L2.resize ( this->num_vars_, false );
    no_H1.resize ( this->num_vars_, false );
    no_H2.resize ( this->num_vars_, false );

    // vel
    L2_vel.resize ( this->num_vars_, false );
    H1_vel.resize ( this->num_vars_, false );
    for ( int i = 0; i < vel_var.size ( ); ++i )
    {
        L2_vel[i] = true;
        H1_vel[i] = true;
    }

    // press
    L2_press.resize ( this->num_vars_, false );
    L2_press[p_var] = true;

    // temp
    L2_temp.resize ( this->num_vars_, false );
    H1_temp.resize ( this->num_vars_, false );
    L2_temp[T_var_] = true;
    H1_temp[T_var_] = true;

    // ********************************************************************************
    // norms
    // velocity
    // L2
    PVT_local_asm_prod_.set_mode ( L2_vel, no_H1, no_H2 );
    std::vector<double> cell_vals;
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    double local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    double global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double vel_l2 = std::sqrt ( global_val );

    // H1
    PVT_local_asm_prod_.set_mode ( no_L2, H1_vel, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double vel_h1 = std::sqrt ( global_val );
    double vel_w12 = std::sqrt ( vel_l2 * vel_l2 + vel_h1 * vel_h1 );

    // pressure
    // L2
    PVT_local_asm_prod_.set_mode ( L2_press, no_H1, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double press_l2 = std::sqrt ( global_val );

    // temperature
    // L2
    PVT_local_asm_prod_.set_mode ( L2_temp, no_H1, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double temp_l2 = std::sqrt ( global_val );

    // H1
    PVT_local_asm_prod_.set_mode ( no_L2, H1_temp, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double temp_h1 = std::sqrt ( global_val );
    double temp_w12 = std::sqrt ( temp_l2 * temp_l2 + temp_h1 * temp_h1 );

    // **************************************************************
    // difference
    LAD::VectorType diff;
    diff.CloneFrom ( *sol );
    diff.Axpy ( *sol_prev, -1. );

    diff.Update ( );

    PVT_local_asm_prod_.set_left_vector ( diff );
    PVT_local_asm_prod_.set_right_vector ( diff );

    // velocity
    // L2
    PVT_local_asm_prod_.set_mode ( L2_vel, no_H1, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double diff_vel_l2 = std::sqrt ( global_val );

    // H1
    PVT_local_asm_prod_.set_mode ( no_L2, H1_vel, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double diff_vel_h1 = std::sqrt ( global_val );
    double diff_vel_w12 = std::sqrt ( diff_vel_l2 * diff_vel_l2 + diff_vel_h1 * diff_vel_h1 );

    // pressure
    // L2
    PVT_local_asm_prod_.set_mode ( L2_press, no_H1, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double diff_press_l2 = std::sqrt ( global_val );

    // temperature
    // L2
    PVT_local_asm_prod_.set_mode ( L2_temp, no_H1, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double diff_temp_l2 = std::sqrt ( global_val );

    // H1
    PVT_local_asm_prod_.set_mode ( no_L2, H1_temp, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( PVT_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double diff_temp_h1 = std::sqrt ( global_val );
    double diff_temp_w12 = std::sqrt ( diff_temp_l2 * diff_temp_l2 + diff_temp_h1 * diff_temp_h1 );

    // file output
    if ( rank ( ) == master_rank ( ) )
    {
        string path = this->root_ + filename;
        ofstream out;
        out.open ( path.c_str ( ), ios::out | ios::app );

        out.precision ( 10 );
        out << std::scientific;
        out << time << " " << vel_l2 << " " << vel_h1 << " " << press_l2 << " " << temp_l2 << " " << temp_h1
                << " " << diff_vel_l2 << " " << diff_vel_h1 << " " << diff_press_l2 << " " << diff_temp_l2 << " " << diff_temp_h1 << "\n";
        out.close ( );

        std::cout << "  L2 (vel_c):             " << vel_l2 << std::endl
                << "  H1 (vel_c):             " << vel_h1 << std::endl
                << "  L2 (vel_c - vel_p):     " << diff_vel_l2 << std::endl
                << "  H1 (vel_c - vel_p):     " << diff_vel_h1 << std::endl
                << "  rel W_12 diff:          " << diff_vel_w12 / vel_w12 << std::endl
                << " ---------------------    " << std::endl
                << "  L2 (press_c):           " << press_l2 << std::endl
                << "  L2 (press_c - press_p): " << diff_press_l2 << std::endl
                << "  rel L_2 diff:           " << diff_press_l2 / press_l2 << std::endl
                << " ---------------------    " << std::endl
                << "  L2 (temp_c):            " << temp_l2 << std::endl
                << "  H1 (temp_c):            " << temp_h1 << std::endl
                << "  L2 (temp_c - temp_p):   " << diff_temp_l2 << std::endl
                << "  H1 (temp_c - temp_p):   " << diff_temp_h1 << std::endl
                << "  rel W_12 diff:          " << diff_temp_w12 / temp_w12 << std::endl
                << " ---------------------    " << std::endl;
    }

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) )
        std::cout << "  took " << timer.get_duration ( ) << " sec" << std::endl;
}
