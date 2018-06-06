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

#include "met_flow_main.h"
#include "common/timer.h"
#include "adaptivity/refinement_strategies.h"

/// Console output macro

static bool CONSOLE_OUTPUT_ACTIVE = true;
static const int CONSOLE_THRESHOLD_LEVEL = 3;

#define CONSOLE_OUTPUT(lvl, x)                                  \
{                                                               \
  if (CONSOLE_OUTPUT_ACTIVE && lvl <= CONSOLE_THRESHOLD_LEVEL)  \
  {                                                             \
    for (int i = 0; i < lvl; ++i)                               \
      std::cout << "  ";                                        \
    std::cout << x << "\n";                                     \
  }                                                             \
}

// For pressure filtering

struct PressureIntegral : private AssemblyAssistant<DIM, double>
{

    PressureIntegral ( const LAD::VectorType& sol, bool cyl_coord, int var ) : sol_ ( sol ), cyl_coord_ ( cyl_coord ), var_ ( var )
    {
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            double& pressure )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
        evaluate_fe_function ( sol_, var_, p_ );

        const int num_q = num_quadrature_points ( );

        // loop over quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            double dJ = 0.0;
            if ( cyl_coord_ )
            {
                const double r = this->x ( q )[1];
                dJ = r * std::fabs ( this->detJ ( q ) );
            }
            else
            {
                dJ = std::abs ( detJ ( q ) );
            }
            pressure += wq * p_[q] * dJ;
        }
    }
    bool cyl_coord_;
    const LAD::VectorType& sol_;
    FunctionValues<double> p_;
    int var_;
};

struct VolumeIntegral : private AssemblyAssistant<DIM, double>
{

    VolumeIntegral ( bool cyl_coord ) : cyl_coord_ ( cyl_coord )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, double& vol )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
        const int num_q = num_quadrature_points ( );

        // loop over quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = this->w ( q );
            double dJ = 0.0;
            if ( cyl_coord_ )
            {
                const double r = this->x ( q )[1];
                dJ = r * std::fabs ( this->detJ ( q ) );
            }
            else
            {
                dJ = std::abs ( detJ ( q ) );
            }
            vol += wq * dJ;
        }
    }
    int cyl_coord_;
};

// *****************************************************************************
// MetFlowApp member function definitions
// *****************************************************************************

/// Constructor: setup application

MetFlowApp::MetFlowApp ( const std::string& root_path, const std::string& param_filename, const std::string base_param_filename, bool resumed_run )
: DynamicMeshProblem<LAD, MESHIMPL, DIM>( MPI_COMM_WORLD ),
comm_ ( MPI_COMM_WORLD ),
rank_ ( -1 ),
num_partitions_ ( -1 ),
master_rank_ ( 0 ),
refinement_level_ ( 0 ),
root_ ( root_path ),
adapt_counter_ ( -1 ),
cur_time_ ( 0. ),
time_step_ ( 0 ),
aug_p_var_ ( DIM + 1 ),
resumed_run_ ( resumed_run ),
num_vars_ ( 0 ),
num_eq_ ( 0 ),
update_every_newton_step_ ( 1 ),
update_time_step_ ( 1 ),
is_linear_problem_ ( false ),
matrix_assembled_ ( false ),
dual_matrix_assembled_ ( false ),
is_stationary_ ( false ),
base_params_ ( base_param_filename.c_str ( ), master_rank ( ), MPI_COMM_WORLD ),
params_ ( param_filename.c_str ( ), master_rank ( ), MPI_COMM_WORLD ),
use_pressure_filter_ ( false ),
use_pressure_filter_dual_ ( false ),
use_hdf5_ ( false ),
ilupp_time_step_ ( 0 ),
ilupp_time_step_dual_ ( 0 ),
problem_mode_ ( 1 ),
pp_step_ ( 1 ),
test_step_ ( 1 ),
lp_div_tol_ ( 1.e6 ),
min_step_primal_ ( 0 ),
max_step_primal_ ( 0 ),
min_step_dual_ ( 0 ),
max_step_dual_ ( 0 )
{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    // Setup Parallel Output / Logging
    if ( rank_ == master_rank_ )
        INFO = true;
    else
        INFO = false;

    // setup linear algebra platform
    la_sys_.Platform = APP_PLATFORM;
    init_platform ( la_sys_ );

    //is_pod_active_ = false;

    if ( !this->resumed_run_ )
    {
        this->create_log_file ( 0, "/log/LastBackupStep", "txt", " ", false );
    }

    MetFlowApp::space_ = new VectorSpace<double> ( );
    MetFlowApp::matrix_ = new LAD::MatrixType;
    MetFlowApp::perturb_ = new LAD::VectorType;
    MetFlowApp::perturb_prev_ = new LAD::VectorType;
    MetFlowApp::solP_ = new LAD::VectorType;
    MetFlowApp::solP_prev_ = new LAD::VectorType;
    MetFlowApp::solP_next_ = new LAD::VectorType;
    MetFlowApp::rhs_ = new LAD::VectorType;
    MetFlowApp::res_ = new LAD::VectorType;
    MetFlowApp::base_ = new LAD::VectorType;
    MetFlowApp::nls_ = new MyNewton<LAD>( this->res_, this->matrix_ );

    MetFlowApp::space_dual_ = new VectorSpace<double> ( );
    MetFlowApp::matrix_dual_ = new LAD::MatrixType;
    MetFlowApp::solD_prev_ = new LAD::VectorType;
    MetFlowApp::solD_next_ = new LAD::VectorType;
    MetFlowApp::solD_ = new LAD::VectorType;
    MetFlowApp::rhs_dual_ = new LAD::VectorType;
    MetFlowApp::res_dual_ = new LAD::VectorType;
}

/// Destructor: clean up after application
/// TODO destroy properly

MetFlowApp::~MetFlowApp ( )
{
}

/// ***************************************************************************
/// MAIN LOOPS
/// ***************************************************************************

void MetFlowApp::run ( )
{
    this->adapt_counter_ = base_params_["Mesh"]["RestartAt"].get<int>( -1 );
    std::string adapt_type = base_params_["Adaptivity"]["DynamicMesh"]["Type"].get<std::string>( );

    // initialize dynamic mesh handler
    this->init_dmh ( adapt_type );
    int num_mesh = this->dmh_->num_mesh_change ( ) + 1;

    // load last spatial and temporal meshes to dmh
    if ( this->adapt_counter_ >= 0 )
    {
        std::string name_space = this->root_ + "/mesh/space_mesh";
        std::string name_time = this->root_ + "/mesh/time_mesh";

        this->t_mesh_.load_all ( name_time, this->adapt_counter_ + 1 );
        this->dmh_->load_mesh_from_file ( this->adapt_counter_, num_mesh, name_space );
    }

    // Prepare application
    // if adapt_counter = -1: build initial spatial mesh
    this->initial_prepare ( );

    // copy initial meshes
    // in case of no restart, new mesh is created -> make copies and store in mesh list
    if ( this->adapt_counter_ < 0 )
    {
        std::string name_space = this->root_ + "/mesh/space_mesh";
        std::string name_time = this->root_ + "/mesh/time_mesh";

        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "> Set initial space mesh " << std::endl;
        }
        this->dmh_->set_initial_mesh ( this->mesh_ );

#ifdef RUN_MODE2
        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "> Save initial space mesh " << std::endl;
        }

        this->dmh_->save_mesh_to_file ( 0, name_space );

        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "> Save initial time mesh " << std::endl;
            this->t_mesh_.save_all ( name_time );
        }
#endif
    }

    // setup fe spaces corresponding to initial meshes
    if ( rank ( ) == master_rank ( ) )
        std::cout << "> Init FE spaces in dmh object " << std::endl;

    this->dmh_->init_fe_spaces ( );

    // setup initial linear algebra objects
    this->prepare_lin_alg_structures ( *this->dmh_->get_couplings_by_index ( 0, 1 ) );

#ifdef RUN_MODE2
    this->prepare_lin_alg_structures_dual ( *this->dmh_->get_couplings_by_index ( 0, -1 ) );
#endif

    if ( this->adapt_counter_ < 0 )
    {
        this->adapt_counter_ = 0;
    }

    this->create_log_file ( this->adapt_counter_, "/log/estimator", "txt", " ", false );
    this->create_log_file ( this->adapt_counter_ * 2, "/log/time_steps", "csv", ", ", false );
    this->create_log_file ( this->adapt_counter_ * 2, "/log/time_estimator", "csv", ", ", false );
    this->create_log_file ( this->adapt_counter_, "/log/problem_size", "csv", ", ", false );

#ifdef RUN_MODE0
    this->dwr_run ( );
#endif
#ifdef RUN_MODE1
    this->dwr_run ( );
#endif
#ifdef RUN_MODE2
    this->dwr_run ( );
#endif
#ifdef RUN_MODE3
    this->perturb_test_run ( );
#endif
#ifdef RUN_MODE4
    this->qoi_run ( );
#endif
#ifdef RUN_MODE5
    this->pp_run ( );
#endif
}

/// Adaptive main loop -> RUN_MODE2

void MetFlowApp::dwr_run ( )
{
    int num_adaptations = base_params_["Adaptivity"]["NumAdaptionCycles"].get<int>( );
    int basis_length_dual = 0;
    int backup_step_primal = 1;

    // ************************************************************
    // ************************************************************
    // ************************************************************
    // Adaption loop
    while ( this->adapt_counter_ < num_adaptations )
    {
        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "Number of adaptions " << num_adaptations << "\n";
            std::cout << "ADAPTION STEP " << this->adapt_counter_ << "\n";
        }

        int time_step = 0;
        this->time_step_ = time_step;

        this->dwr_loop_prepare ( );

        // Log data
        std::map<std::string, double> log_data;

        // ************************************************************
        // PRIMAL PROBLEM
        // ************************************************************
        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "NOW THE PRIMAL PROBLEM (AC=" << this->adapt_counter_ << ")" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }

        Timer timer;
        timer.start ( );

        std::string primal_mode = base_params_["BackUp"]["PrimalMode"].get<std::string>( );

        if ( primal_mode == "Calculate" )
        {
            // ************************************************************
            // ************************************************************
            // Calculate primal solution

            if ( rank ( ) == master_rank ( ) )
                std::cout << "> PrimalMode = " << primal_mode << std::endl;
            time_step = 0;

            // initialized or resume simulation
            if ( rank ( ) == master_rank ( ) )
                std::cout << "> Initializing Simulation" << std::endl;

            // ************************************************************
            // init primal initial conditions
#ifdef RUN_MODE0
            time_step = this->init_run ( true );
            return;
#else
            time_step = this->init_run ( false );
#endif
            this->time_step_ = time_step;

            // ************************************************************
            // Solve primal problem

            time_step = this->solve_primal ( this->adapt_counter_, time_step );

            //this->num_time_steps_ = time_step;
            this->min_step_primal_ = 0;
            if ( !this->is_stationary_ )
            {
                this->max_step_primal_ = this->get_num_intervals ( );
                this->max_step_dual_ = this->get_num_intervals ( );
            }
            else
            {
                this->max_step_primal_ = 0;
                this->max_step_dual_ = 0;
            }
            this->cur_time_ = this->duration_;
            this->time_step_ = this->get_num_intervals ( );

        }
        else if ( primal_mode == "Load" )
        {
            // ************************************************************
            // ************************************************************
            // Do nothing with primal solution

            if ( rank ( ) == master_rank ( ) )
                std::cout << "-> PrimalMode = " << primal_mode << std::endl;

            time_step = this->get_num_intervals ( );
            this->time_step_ = time_step;
        }

        timer.stop ( );

        if ( rank ( ) == master_rank ( ) )
            std::cout << "Primal problem solved." << std::endl;

#ifdef RUN_MODE1
        return;
#endif
        // ************************************************************
        // DUAL PROBLEM
        // ************************************************************
        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "NOW THE DUAL PROBLEM (AC=" << this->adapt_counter_ << ")" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }

        timer.reset ( );
        timer.start ( );
        std::string dual_mode = base_params_["BackUp"]["DualMode"].get<std::string>( );

        if ( dual_mode == "Calculate" )
        {
            // ************************************************************
            // ************************************************************
            // Calculate dual solution

            if ( rank ( ) == master_rank ( ) )
                std::cout << "-> DualMode = " << dual_mode << std::endl;

            // ************************************************************
            // Init initial condition for dual solution
            time_step = this->init_run_dual ( );
            this->time_step_ = time_step;

            // ************************************************************
            // Solve dual problem
            time_step = this->solve_dual ( this->adapt_counter_, time_step );
            this->time_step_ = time_step;
        }
        else if ( dual_mode == "Load" )
        {
            // ************************************************************
            // ************************************************************
            // Do nothing with dual solution

            if ( rank ( ) == master_rank ( ) )
                std::cout << "-> DualMode = " << dual_mode << std::endl;
        }

        timer.stop ( );
        //log_data["Dual Problem, Time"] = timer.get_duration();

        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "DUAL PROBLEM SOLVED ... (AC=" << this->adapt_counter_ << ")" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }

        // ************************************************************
        // ERROR ESTIMATION & MESH ADAPTION
        // ************************************************************
        // compute error indicators
        this->estimate_error ( this->adapt_counter_, this->get_num_intervals ( ) );

        // determine points in time when to change the mesh
        this->adapt_mesh_change_list ( );

        // compute reduced estimators -> combine estimators for several time steps to one estimation
        this->compute_reduced_estimators ( );

        // adapt space and time mesh
        this->adapt_temporal_mesh ( );
        this->adapt_spatial_mesh ( );
        this->dmh_->init_fe_spaces ( );

        // write out some statistics
        int num_intervals = this->get_num_intervals ( );
        std::vector<int> num_cells ( this->dmh_->num_mesh ( ) );
        for ( int m = 0; m<this->dmh_->num_mesh ( ); ++m )
        {
            num_cells[m] = this->dmh_->get_mesh_by_index ( m )->num_global_cells ( this->comm_ );
        }
        if ( rank_ == master_rank_ )
        {
            std::stringstream log_pre;
            log_pre << this->root_ << "/log/problem_size.csv";
            std::string log_filename = log_pre.str ( );

            ofstream logfile;
            logfile.open ( log_filename.c_str ( ), std::ios_base::app );
            logfile << this->adapt_counter_ + 1 << ", ";

            for ( int m = 0; m<this->dmh_->num_mesh ( ); ++m )
            {
                logfile << num_cells[m] << ", ";
            }
            logfile << num_intervals << " \n";
            logfile.close ( );
        }
        // reset problem
        this->clear_problem ( );

        this->adapt_counter_++;
    }
}

/// Post Process HDF5 data -> RUN_MODE5

void MetFlowApp::perturb_test_run ( )
{
    // Initialize problem
    this->initial_prepare ( );
    this->create_log_file ( 0, "/log/dJ_dp_ef", "txt" );
    this->create_log_file ( 0, "/log/dJ_dp_ic", "txt" );

    this->local_asm_primal_->set_goal_functional ( this->j_ );

    int frequence = this->base_params_["Perturbation"]["Test"]["TimeStep"].get<int>( 10 );
    int final_step = this->base_params_["Perturbation"]["Test"]["FinalStep"].get<int>( 1000 );

    // For HDF5
    std::string filename_dual = this->base_params_["Perturbation"]["Test"]["DualSolution"]["SnapshotsIn"].get<std::string>( );
    std::string groupname_dual = this->base_params_["Perturbation"]["Test"]["DualSolution"]["SnapshotsGroup"].get<std::string>( );
    std::string prefix_dual = this->base_params_["Perturbation"]["Test"]["DualSolution"]["SnapshotsPrefix"].get<std::string>( );
    filename_dual = this->root_ + "/" + filename_dual;

    std::string filename_ef = this->base_params_["Perturbation"]["ExternalForce"]["SnapshotsIn"].get<std::string>( );
    std::string groupname_ef = this->base_params_["Perturbation"]["ExternalForce"]["SnapshotsGroup"].get<std::string>( );
    std::string prefix_ef = this->base_params_["Perturbation"]["ExternalForce"]["SnapshotsPrefix"].get<std::string>( );
    filename_ef = this->root_ + "/" + filename_ef;

    std::string filename_ic = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsIn"].get<std::string>( );
    std::string groupname_ic = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsGroup"].get<std::string>( );
    std::string prefix_ic = this->base_params_["Perturbation"]["InitialCondition"]["SnapshotsPrefix"].get<std::string>( );
    filename_ic = this->root_ + "/" + filename_ic;

    this->cur_time_ = 0.;

    int perturb_ic = this->base_params_["Perturbation"]["PerturbIC"].get<int>( );
    int perturb_ef = this->base_params_["Perturbation"]["PerturbEF"].get<int>( );

    if ( perturb_ic != 0 )
    {
        std::stringstream ss_ic;
        ss_ic << prefix_ic << 0;
        MetFlowApp::read_file ( *this->solP_, filename_ic, groupname_ic, ss_ic.str ( ) );

        std::stringstream ss_dual;
        ss_dual << prefix_dual << 0;
        MetFlowApp::read_file ( *this->solP_prev_, filename_dual, groupname_dual, ss_dual.str ( ) );

        this->solP_->Update ( );
        this->solP_prev_->Update ( );

        this->local_asm_primal_->set_solP ( *this->solP_ );
        this->local_asm_primal_->set_solP_prev ( *this->solP_prev_ );
        this->local_asm_primal_->set_solP_next ( *this->solP_next_ );

        MetFlowApp::evaluate_inner_prod ( "/log/dJ_dp_ic", this->solP_, this->solP_prev_, this->cur_time_ );
    }
    if ( perturb_ef == 0 )
        final_step = 0;

    for ( int l = 1; l <= final_step; ++l )
    {
        this->delta_t_ = this->get_delta_t ( l );
        this->cur_time_ += this->delta_t_;

        if ( l % frequence == 0 )
        {
            std::stringstream ss_ef;
            ss_ef << prefix_ef << l;
            MetFlowApp::read_file ( *this->solP_, filename_ef, groupname_ef, ss_ef.str ( ) );

            std::stringstream ss_dual;
            ss_dual << prefix_dual << l;
            MetFlowApp::read_file ( *this->solP_prev_, filename_dual, groupname_dual, ss_dual.str ( ) );

            this->solP_->Update ( );
            this->solP_prev_->Update ( );

            this->local_asm_primal_->set_solP ( *this->solP_ );
            this->local_asm_primal_->set_solP_prev ( *this->solP_prev_ );
            this->local_asm_primal_->set_solP_next ( *this->solP_next_ );

            MetFlowApp::evaluate_inner_prod ( "/log/dJ_dp_ef", this->solP_, this->solP_prev_, this->cur_time_ );

            if ( rank ( ) == master_rank ( ) )
                std::cout << "  " << l << " / " << final_step << " done " << std::endl;
        }
    }
}

/// Post Process HDF5 data -> RUN_MODE4

void MetFlowApp::qoi_run ( )
{
    this->prepare_goal_functional ( );
    this->local_asm_primal_->set_goal_functional ( this->j_ );

    int frequence = this->base_params_["Perturbation"]["Test"]["TimeStep"].get<int>( 10 );
    int final_step = this->base_params_["Perturbation"]["Test"]["FinalStep"].get<int>( this->get_num_intervals ( ) );

    // For HDF5
    std::string filename = this->base_params_["Perturbation"]["Test"]["PerturbedSolution"]["SnapshotsIn"].get<std::string>( );
    std::string groupname = this->base_params_["Perturbation"]["Test"]["PerturbedSolution"]["SnapshotsGroup"].get<std::string>( );
    std::string prefix = this->base_params_["Perturbation"]["Test"]["PerturbedSolution"]["SnapshotsPrefix"].get<std::string>( );
    filename = this->root_ + "/" + filename;

    this->cur_time_ = 0.;

    for ( int l = 0; l <= final_step; ++l )
    {
        this->delta_t_ = this->get_delta_t ( l );
        this->local_asm_primal_->set_dT_pc ( this->delta_t_ );
        this->local_asm_primal_->set_time ( this->cur_time_ );

        if ( l % frequence == 0 )
        {
            std::stringstream ss;
            ss << prefix << l;
            MetFlowApp::read_file ( *this->solP_, filename, groupname, ss.str ( ) );

            this->solP_->Update ( );

            local_asm_primal_->set_solP ( *this->solP_ );

            MetFlowApp::evaluate_qoi_int ( "/log/QoI_int", this->solP_, this->cur_time_ );

            if ( l == final_step )
                MetFlowApp::evaluate_qoi_fin ( "/log/QoI_fin", this->solP_, this->cur_time_ );

            if ( rank ( ) == master_rank ( ) )
                std::cout << "  " << l << " / " << final_step << " done " << std::endl;
        }
        this->cur_time_ += this->delta_t_;
    }
}

/// ***************************************************************************
/// NONLINEAR PROBLEM
/// ***************************************************************************

void MetFlowApp::EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F )
{
    // compute the residual vector
    if ( this->problem_mode_ == 1 )
    {
        compute_residual ( &u, F, 1 );
    }
    else if ( this->problem_mode_ == -1 )
    {
        compute_residual ( NULL, F, -1 );
    }
}

void MetFlowApp::EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
    // assemble the matrix for the linearized system
    if ( this->problem_mode_ == 1 )
    {
        this->compute_matrix ( &u, DF, 1 );
    }
    else if ( this->problem_mode_ == -1 )
    {
        this->compute_matrix ( NULL, DF, -1 );
    }
}

void MetFlowApp::compute_residual ( const LAD::VectorType* u, LAD::VectorType* F, int mode )
{
    Timer timer;
    timer.reset ( );
    timer.start ( );

    if ( u != NULL && mode == 1 )
    {
        local_asm_primal_->set_newton_solution ( *u );
    }
    if ( mode == 1 )
    {
        global_asm_.assemble_vector ( *space_, boost::ref ( *local_asm_primal_ ), *F );
        if ( !dirichlet_dofs_.empty ( ) )
        {
            std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
            F->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( zeros ) );
        }
    }
    else if ( mode == -1 )
    {
        global_asm_.assemble_vector ( *space_dual_, boost::ref ( *local_asm_dual_ ), *F );
        if ( !dirichlet_dofs_dual_.empty ( ) )
        {
            std::vector<LAD::DataType> zeros ( dirichlet_dofs_dual_.size ( ), 0. );
            F->SetValues ( vec2ptr ( dirichlet_dofs_dual_ ), dirichlet_dofs_dual_.size ( ), vec2ptr ( zeros ) );
        }
    }

    timer.stop ( );

    if ( rank ( ) == master_rank ( ) ) std::cout << "  Vector assembly time: " << timer.get_duration ( ) << "s, " << std::endl;

    //  Note: Hiflow convention: J dx = res, x<-x-dx
    F->Update ( );
}

void MetFlowApp::compute_matrix ( const LAD::VectorType* u, LAD::MatrixType* DF, int mode )
{
    assert ( mode == 1 || mode == -1 );

    bool write_time_statistics = false;
    Timer timer;
    timer.reset ( );

    // ASSEMBLY
    bool constant_matrix = false;
    bool matrix_assembled = false;
    if ( mode == 1 )
    {
        constant_matrix = this->local_asm_primal_->is_matrix_constant ( );
        matrix_assembled = this->matrix_assembled_;
    }
    else if ( mode == -1 )
    {
        constant_matrix = this->local_asm_dual_->is_matrix_constant ( );
        matrix_assembled = this->dual_matrix_assembled_;
    }

    if ( !constant_matrix || !matrix_assembled )
    {
        DF->Zeros ( );
        if ( u != NULL && mode == 1 )
        {
            local_asm_primal_->set_newton_solution ( *u );
        }
        timer.start ( );

        if ( mode == 1 )
        {
            global_asm_.assemble_matrix ( *space_, boost::ref ( *local_asm_primal_ ), *DF );
            if ( !dirichlet_dofs_.empty ( ) )
            {
                DF->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
            }
            this->matrix_assembled_ = true;
        }
        else if ( mode == -1 )
        {
            global_asm_.assemble_matrix ( *space_, boost::ref ( *local_asm_dual_ ), *DF );
            if ( !dirichlet_dofs_dual_.empty ( ) )
            {
                DF->diagonalize_rows ( vec2ptr ( dirichlet_dofs_dual_ ), dirichlet_dofs_dual_.size ( ), 1. );
            }
            this->dual_matrix_assembled_ = true;
        }

        timer.stop ( );
        if ( rank ( ) == master_rank ( ) ) std::cout << "  Matrix assembly time: " << timer.get_duration ( ) << "s, " << std::endl;

        // PRECONDITIONING
        int iter = 0;
        if ( mode == 1 )
        {
            iter = nls_->iter ( );
        }
        else if ( mode == -1 )
        {
            //iter = nls_dual_->iter();
            iter = 0;
        }

        timer.reset ( );

        if ( update_every_newton_step_ != 0 && update_every_newton_step_ != 1 )
        {
            // not setting the parameter results in "0"
            std::cout << "Preconditioner|UpdateEveryNewtonStep must be either 0 or 1, but is: " << update_every_newton_step_ << std::endl;
            quit_program ( );
        }

        // Update of Preconditioner
        if ( ( ( update_every_newton_step_ == 1 ) || ( iter == 1 ) ) && ( ilupp_time_step_ % update_time_step_ == 0 ) )
        {
            timer.start ( );
            if ( mode == 1 )
            {
                this->update_preconditioner ( *u, DF );
            }
            else if ( mode == -1 )
            {
                this->update_preconditioner_dual ( DF );
            }
            timer.stop ( );
        }
        ++ilupp_time_step_;

        if ( rank ( ) == master_rank ( ) )
            std::cout << "  Precond setup time: " << timer.get_duration ( ) << " s, " << std::endl;
        if ( write_time_statistics )
            std::cout << std::endl;
    }
}

/// ***************************************************************************
/// INITIALIZATION
/// ***************************************************************************

void MetFlowApp::clear_problem ( )
{
    if ( rank ( ) == master_rank ( ) )
        std::cout << "> Clear MetFlow Problem" << std::endl;

    // clear FE stuff
    this->coupling_vars_.clear ( );
    this->coupling_vars_dual_.clear ( );

    this->matrix_assembled_ = false;
    this->dual_matrix_assembled_ = false;

    this->fineP_.Clear ( );
    this->fineP_prev_.Clear ( );
    this->fineP_next_.Clear ( );
    this->fineD_.Clear ( );
    this->fineD_prev_.Clear ( );
    this->fineD_next_.Clear ( );

    // clear assemblers
    this->local_asm_primal_->clear ( );
    this->local_asm_dual_->clear ( );
    this->local_asm_est_->clear ( );

    // Clear interpolation
    this->patch_interpolation_.clear ( );
    this->patch_interpolation_dual_.clear ( );

    // Clear dirichlet boundary conditions
    this->dirichlet_dofs_.clear ( );
    this->dirichlet_values_.clear ( );
    this->dirichlet_dofs_dual_.clear ( );
    this->dirichlet_values_dual_.clear ( );

    // Clear estimators
    this->time_indicator_.clear ( );
    this->space_indicator_.clear ( );
    this->reduced_space_indicator_.clear ( );
    this->reduced_time_indicator_.clear ( );

    this->local_equation_estimator_.clear ( );
    this->local_time_estimator_ = 0.;
    this->local_space_estimator_ = 0.;
    this->global_equation_estimator_.clear ( );
    this->global_time_estimator_ = 0.;
    this->global_space_estimator_ = 0.;
    this->global_estimator_ = 0.;

    // other stuff
    this->period_.clear ( );
    this->cur_time_ = 0.;
    this->time_step_ = 0;
    this->min_step_primal_ = 0;
    this->max_step_primal_ = 0;
    this->min_step_dual_ = 0;
    this->max_step_dual_ = 0;

    this->cur_time_ = 0.;
    this->time_step_ = 0;
    this->problem_mode_ = 1;
    this->num_vars_ = 0;
    this->num_eq_ = 0;
    this->update_every_newton_step_ = 1;
    this->update_time_step_ = 1;
    this->is_linear_problem_ = false;
    this->is_stationary_ = false;
    this->use_pressure_filter_ = false;
    this->use_pressure_filter_dual_ = false;
    this->ilupp_time_step_ = 0;
    this->ilupp_time_step_dual_ = 0;
    this->pp_step_ = 1;
    this->test_step_ = 1;
    this->lp_div_tol_ = 1.e6;

    this->active_mesh_index_ = -1;
}

void MetFlowApp::dwr_loop_prepare ( )
{
    this->prepare_assembler ( );
#ifdef RUN_MODE2
    this->prepare_assembler_dual ( );
#endif
}

void MetFlowApp::prepare_periodicity ( )
{
    period_.clear ( );
    if ( base_params_["Mesh"].contains ( "xPeriodicity" ) )
    {
        double h, master, slave;
        h = base_params_["Mesh"]["xPeriodicity"]["h"] .get<double>( );
        master = base_params_["Mesh"]["xPeriodicity"]["Master"] .get<double>( );
        slave = base_params_["Mesh"]["xPeriodicity"]["Slave"] .get<double>( );
        period_.push_back ( MasterSlave ( master, slave, h, 0 ) );
    }
    if ( base_params_["Mesh"].contains ( "yPeriodicity" ) )
    {
        double h, master, slave;
        h = base_params_["Mesh"]["yPeriodicity"]["h"] .get<double>( );
        master = base_params_["Mesh"]["yPeriodicity"]["Master"] .get<double>( );
        slave = base_params_["Mesh"]["yPeriodicity"]["Slave"] .get<double>( );
        period_.push_back ( MasterSlave ( master, slave, h, 1 ) );
    }
    if ( base_params_["Mesh"].contains ( "zPeriodicity" ) )
    {
        double h, master, slave;
        h = base_params_["Mesh"]["zPeriodicity"]["h"] .get<double>( );
        master = base_params_["Mesh"]["zPeriodicity"]["Master"] .get<double>( );
        slave = base_params_["Mesh"]["zPeriodicity"]["Slave"] .get<double>( );
        period_.push_back ( MasterSlave ( master, slave, h, 2 ) );
    }
#ifdef INFINITE_CYLINDER
    double h, master, slave;
    h = 0.5 * base_params_["Mesh"]["Height"] .get<double>( );
    master = 0.;
    slave = base_params_["Mesh"]["Height"] .get<double>( );
    period_.push_back ( MasterSlave ( master, slave, h, 2 ) );
#endif
}

void MetFlowApp::build_initial_mesh ( int adapt_counter )
{
    if ( adapt_counter >= 0 )
    {
        return;
    }

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Build sequential initial mesh " << std::endl;

    // read in mesh
    const std::string mesh_filename = this->root_ + "/" + base_params_["Mesh"]["Filename"].get<std::string>( );

#ifdef USE_PXEST
    this->build_initial_mesh_p4est ( adapt_counter );
#else
    // read in and refine
    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << ">  read mesh from file " << mesh_filename << std::endl;
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIM, DIM, 0, period_, mesh::IMPL_DBVIEW );
        int init_ref_lvl = base_params_["Mesh"]["InitialRefLevel"].get<int>( 0 );

        std::cout << ">  refine sequential mesh " << init_ref_lvl << " times " << std::endl;
        if ( init_ref_lvl > 0 )
        {
            master_mesh_ = master_mesh_->refine_uniform_seq ( init_ref_lvl );
        }
        refinement_level_ = init_ref_lvl;
    }
    MPI_Bcast ( &refinement_level_, 1, MPI_INT, master_rank ( ), comm_ );

    // partition
#    ifdef WITH_METIS
    MetisGraphPartitioner partitioner;
    if ( rank ( ) == master_rank ( ) ) std::cout << ">  Partition initial mesh with metis " << std::endl;
#    else
    NaiveGraphPartitioner partitioner;
    if ( rank ( ) == master_rank ( ) ) std::cout << ">  Partition initial mesh with naive partitioner " << std::endl;
#    endif
    //const GraphPartitioner* p = (num_partitions() == 1 ? 0 : &partitioner);
    const GraphPartitioner* p = &partitioner;

    MeshPtr local_mesh = partition_and_distribute ( master_mesh_, master_rank ( ), comm_, p );
    int num_local_cells = local_mesh->num_entities ( DIM );
    assert ( local_mesh != 0 );

    // compute ghost cells
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts, mesh::IMPL_DBVIEW, 1 );

    int num_ghost_cells = mesh_->num_entities ( DIM ) - num_local_cells;
    std::cout << "   #local cells: " << num_local_cells
            << ", #ghost cells: " << num_ghost_cells
            << ", #sum: " << mesh_->num_entities ( DIM )
            << std::endl;
#endif

    PVtkWriter writer ( comm_ );
    std::string output_file = this->root_ + "/mesh/space_mesh.0.pvtu";
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

void MetFlowApp::build_initial_mesh_parallel ( int adapt_counter )
{
    if ( adapt_counter >= 0 )
    {
        return;
    }

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Build parallel initial mesh " << std::endl;

    int init_ref_lvl = base_params_["Mesh"]["InitialRefLevel"].get<int>( 0 );

    // Read in mesh on master
    if ( rank ( ) == master_rank ( ) ) std::cout << ">  Read in mesh on master " << std::endl;
    const std::string mesh_filename = this->root_ + "/" + base_params_["Mesh"]["Filename"].get<std::string>( );

#ifdef USE_PXEST
    this->build_initial_mesh_p4est ( adapt_counter );
#else
#    ifdef WITH_PARMETIS
    if ( rank ( ) == master_rank ( ) )
    {
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIM, DIM, 0, period_, mesh::IMPL_DBVIEW );
        // uniform sequential refinement
        int seq_ref_lvl = base_params_["Mesh"]["SequentialRefLevel"].get<int>( 0 );
        if ( seq_ref_lvl > 0 )
        {
            master_mesh_ = master_mesh_->refine_uniform_seq ( seq_ref_lvl );
        }
        refinement_level_ = seq_ref_lvl;
    }
    MPI_Bcast ( &refinement_level_, 1, MPI_INT, master_rank ( ), this->comm_ );

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Partition initial mesh with parmetis " << std::endl;

    ParMetisGraphPartitioner par1;
    const GraphPartitioner* p1 = &par1;
    int uniform_ref_steps = 0;

    MeshPtr local_mesh_tmp = partition_and_distribute ( this->master_mesh_, master_rank ( ), this->comm_, &uniform_ref_steps, mesh::IMPL_DBVIEW );
    assert ( local_mesh_tmp != 0 );

    // Refine mesh locally
    if ( rank ( ) == master_rank ( ) ) std::cout << ">  Initial refinement of mesh " << std::endl;

    for ( int r = refinement_level_; r < init_ref_lvl; ++r )
    {
        local_mesh_tmp = local_mesh_tmp->refine ( );
    }
    refinement_level_ = init_ref_lvl;

    // Repartition Mesh
    if ( rank ( ) == master_rank ( ) ) std::cout << ">  Repartition mesh " << std::endl;
    ParMetisGraphPartitioner par2;

    MeshPtr local_mesh = repartition_mesh ( local_mesh_tmp, this->comm_, &par2 );
    assert ( local_mesh != 0 );

    int num_local_cells = local_mesh->num_entities ( DIM );

    // Compute ghosts
    MPI_Barrier ( this->comm_ );
    if ( rank ( ) == master_rank ( ) ) std::cout << ">  Compute ghosts " << std::endl;

    SharedVertexTable shared_verts;
    this->mesh_ = compute_ghost_cells ( *local_mesh, this->comm_, shared_verts, mesh::IMPL_DBVIEW, 1 );

    int num_ghost_cells = mesh_->num_entities ( DIM ) - num_local_cells;
    std::cout << "   #local cells: " << num_local_cells
            << ", #ghost cells: " << num_ghost_cells
            << ", #sum: " << mesh_->num_entities ( DIM )
            << std::endl;

    MPI_Barrier ( this->comm_ );
    if ( rank ( ) == master_rank ( ) ) std::cout << ">  Write out mesh " << std::endl;
    PVtkWriter writer ( this->comm_ );
    std::string output_file = this->root_ + "/mesh/space_mesh.0.pvtu";
    writer.add_all_attributes ( *this->mesh_, true );
    writer.write ( output_file.c_str ( ), *this->mesh_ );

#    else
    std::cout << "Parallel refinenement needs parmetis or p4est !! ";
#    endif
#endif
}

void MetFlowApp::build_initial_mesh_p4est ( int adapt_counter )
{
    if ( adapt_counter >= 0 )
    {
        return;
    }

    int init_ref_lvl = base_params_["Mesh"]["InitialRefLevel"].get<int>( 0 );
    const std::string mesh_filename = this->root_ + "/" + base_params_["Mesh"]["Filename"].get<std::string>( );

#ifdef USE_PXEST
    if ( rank ( ) == master_rank ( ) )
    {
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIM, DIM, 0, period_, mesh::IMPL_P4EST );

        // uniform sequential refinement
        int seq_ref_lvl = base_params_["Mesh"]["SequentialRefLevel"].get<int>( 0 );
        std::cout << ">  refine sequential mesh " << seq_ref_lvl << " times " << std::endl;
        if ( seq_ref_lvl > 0 )
        {
            master_mesh_ = master_mesh_->refine_uniform_seq ( seq_ref_lvl );
        }
        refinement_level_ = seq_ref_lvl;
    }
    MPI_Bcast ( &refinement_level_, 1, MPI_INT, master_rank ( ), this->comm_ );

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Partition initial mesh with p4est " << std::endl;

    int uniform_ref_steps;
    SharedVertexTable shared_verts;

    mesh_ = partition_and_distribute ( master_mesh_, master_rank ( ), comm_, &uniform_ref_steps, mesh::IMPL_P4EST );
    boost::intrusive_ptr<MeshPXest> mesh_pXest = boost::static_pointer_cast<MeshPXest> ( mesh_ );

    //  mesh_ = compute_ghost_cells      ( *mesh_, comm_, shared_verts, mesh::IMPL_P4EST, 2 );

    refinement_level_ += uniform_ref_steps;
    for ( int r = refinement_level_; r < init_ref_lvl; ++r )
    {
        mesh_ = mesh_->refine ( );
    }

    /*
    std::vector<int> ref(mesh_->num_entities(DIM), 0);
    for (int l=0; l<mesh_->num_entities(DIM); ++l)
    {
        if (l%3==0)
        {
            ref[l] = 1;
        }
    }
    mesh_pXest->set_patch_mode   (true);
    mesh_pXest->set_connection_mode (2);
    mesh_ = mesh_->refine ( ref );
     */
    mesh_ = compute_ghost_cells ( *mesh_, this->comm_, shared_verts, mesh::IMPL_P4EST, GHOST_LAYER_WIDTH );

    refinement_level_ = init_ref_lvl;

    std::cout << "  #local cells: " << mesh_->num_local_cells ( )
            << ", #ghost cells: " << mesh_->num_ghost_cells ( )
            << ", #sum: " << mesh_->num_entities ( DIM )
            << std::endl;

#endif
}

void MetFlowApp::build_initial_time_mesh ( int adapt_counter )
{
    const double eps = 1e-8;
    if ( adapt_counter >= 0 )
    {
        return;
    }

    this->t_mesh_.clear ( );

    this->delta_t_ = base_params_["TimeDiscretization"]["DeltaT"].get<double>( );
    double time = 0.;

    while ( time < base_params_["TimeDiscretization"]["EndTime"] .get<double>( ) - eps )
    {
        time += this->delta_t_;
        this->t_mesh_.add_time ( time, 0 );

    }
    this->num_time_steps_ = this->t_mesh_.num_intervals ( ) + 1;
    this->t_mesh_.make_regular ( 0 );

    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "> Prepare time mesh with " << this->num_time_steps_ - 1 << " intervals " << std::endl;

        std::string prefix = "time_mesh";

        std::stringstream pre;
        pre << this->root_ << "/mesh/" << prefix << "." << 0 << ".csv";
        std::string visu_filename = pre.str ( );
        std::cout << "> Save time mesh as " << visu_filename << std::endl;

        std::stringstream log_pre;
        log_pre << this->root_ << "/log/time_steps.csv";
        std::string log_filename = log_pre.str ( );

        ofstream myfile;
        myfile.open ( visu_filename.c_str ( ) );
        ofstream logfile;
        logfile.open ( log_filename.c_str ( ), std::ios_base::app );

        for ( int t = 0; t<this->t_mesh_.num_intervals ( ) + 1; ++t )
        {
            myfile << this->t_mesh_.time ( t, 0 ) << ", ";
            logfile << this->t_mesh_.time ( t, 0 ) << ", ";
        }
        myfile << "\n";
        myfile.close ( );
        logfile << "\n";
        for ( int t = 0; t<this->t_mesh_.num_intervals ( 0 ); ++t )
        {
            logfile << this->t_mesh_.time ( t + 1, 0 ) - this->t_mesh_.time ( t, 0 ) << ", ";
        }
        logfile << "\n";
        logfile.close ( );
    }
}

void MetFlowApp::prepare_parameters ( )
{
    duration_ = this->base_params_["TimeDiscretization"]["EndTime"].get<double>( );
    visual_step_ = this->base_params_["Visualization"]["VisTimestep"].get<int>( );
    if ( visual_step_ <= 0 ) visual_step_ = 1; // default value

    backup_step_ = this->base_params_["BackUp"]["BackupTimestep"].get<int>( );
    if ( backup_step_ <= 0 ) backup_step_ = 1; // default value

    this->pp_step_ = this->base_params_["PostProcessing"]["PPTimeStep"].get<int>( );
    this->print_level_ = this->base_params_["PrimalLinearSolver"]["PrintLevel"].get<int>( );
    this->test_step_ = this->base_params_["DualLinearSolver"]["TestStep"].get<int>( );
    this->lp_div_tol_ = this->base_params_["DualLinearSolver"]["DivergenceLimit"].get<int>( );

    // output parameters for debugging
    LOG_INFO ( "base_parameters", base_params_ );
}

void MetFlowApp::prepare_time_method ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare time discretization" << std::endl;

    method_ = base_params_["TimeDiscretization"]["Method"].get<std::string>( );
    theta_ = base_params_["TimeDiscretization"]["Theta"] .get<double>( );
    int mod_gal = base_params_["TimeDiscretization"]["ModifiedGalerkin"] .get<int>( );
    if ( mod_gal == 0 )
        mod_galerkin_ = false;
    else
        mod_galerkin_ = true;

    interminable_assert ( method_ == "GalerkinCD" || method_ == "Theta" || method_ == "Stationary" );

    if ( method_ == "GalerkinCD" )
    {
        this->local_asm_primal_->set_time_discretization_to_galerkin_cd ( mod_galerkin_ );
    }
    if ( method_ == "Theta" )
    {
        this->local_asm_primal_->set_time_discretization_to_theta ( this->theta_ );
    }
    if ( method_ == "Stationary" )
    {
        this->is_stationary_ = true;
        this->local_asm_primal_->set_time_discretization_off ( );
    }
}

void MetFlowApp::prepare_time_method_dual ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare dual time discretization" << std::endl;
    method_dual_ = base_params_["TimeDiscretization"]["DualMethod"].get<std::string>( );
    theta_dual_ = base_params_["TimeDiscretization"]["DualTheta"] .get<double>( );

    interminable_assert ( method_dual_ == "GalerkinDC" || method_dual_ == "GalerkinCD" || method_dual_ == "Theta" || method_dual_ == "Simple" || method_dual_ == "Stationary" );

    if ( method_dual_ == "GalerkinCD" )
    {
        this->local_asm_dual_->set_time_discretization_to_galerkin_cd ( false );
    }
    if ( method_dual_ == "GalerkinDC" )
    {
        this->local_asm_dual_->set_time_discretization_to_galerkin_dc ( false );
    }
    if ( method_dual_ == "Theta" )
    {
        this->local_asm_dual_->set_time_discretization_to_theta ( this->theta_dual_ );
    }
    if ( method_dual_ == "Simple" )
    {
        this->local_asm_dual_->set_time_discretization_to_simple ( );
    }
    if ( method_dual_ == "Stationary" )
    {
        this->is_stationary_ = true;
        this->local_asm_dual_->set_time_discretization_off ( );
    }
    /*
    if (method_dual_ != "GalerkinDC")
    {
        this->solP_prev_->Zeros();
        this->solP_prev_->Update();
    }
     * */
}

void MetFlowApp::periodify_space ( VectorSpace<DATATYPE>& space, MeshPtr mesh )
{
    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "> Periodify space" << std::endl;
    }

    assert ( mesh != 0 );
    const int tdim = mesh->tdim ( );

    // loop over all cells and init the cell transformations
    for ( mesh::EntityIterator it = mesh->begin ( tdim ); it != mesh->end ( tdim ); ++it )
    {
        Coord coord_vtx;
        it->get_coordinates ( coord_vtx );

        std::vector<MasterSlave> period = mesh->get_period ( );
        coord_vtx = unperiodify ( coord_vtx, DIM, period );
        space.fe_manager ( ).get_cell_transformation ( it->index ( ) )->reinit ( coord_vtx );
    }
}

void MetFlowApp::prepare_lin_alg_structures ( Couplings<double>& couplings )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare LA structures for primal problem " << std::endl;
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( *this->space_, sparsity, &this->coupling_vars_ );

    // Initialize linear algebra structures
    couplings.Init ( comm_, this->space_->dof ( ) );
    couplings.InitializeCouplings ( sparsity.off_diagonal_rows, sparsity.off_diagonal_cols );

    // Initialize matrices and vectors
#ifdef USE_HYPRE
    matrix_->Init ( comm_, couplings );
#else
    matrix_->Init ( comm_, couplings, la_sys_.Platform, MetFlowApp::APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
#endif

    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );

#ifdef USE_HYPRE
    solP_ ->Init ( comm_, couplings );
    solP_prev_ ->Init ( comm_, couplings );
    solP_next_ ->Init ( comm_, couplings );
    rhs_ ->Init ( comm_, couplings );
    perturb_ ->Init ( comm_, couplings );
    perturb_prev_->Init ( comm_, couplings );
    base_ ->Init ( comm_, couplings );
#else
    solP_ ->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    solP_ ->InitStructure ( );
    solP_prev_ ->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    solP_prev_ ->InitStructure ( );
    solP_next_ ->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    solP_next_ ->InitStructure ( );
    rhs_ ->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    rhs_ ->InitStructure ( );
    base_ ->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    base_ ->InitStructure ( );
    perturb_ ->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    perturb_ ->InitStructure ( );
    perturb_prev_->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    perturb_prev_->InitStructure ( );
#endif

    solP_ ->Zeros ( );
    solP_prev_ ->Zeros ( );
    solP_next_ ->Zeros ( );
    rhs_ ->Zeros ( );
    perturb_ ->Zeros ( );
    perturb_prev_->Zeros ( );
    base_ ->Zeros ( );

    // prepare dirichlet BC
    this->prepare_bc ( );
    this->set_bc ( 1 );
}

void MetFlowApp::prepare_lin_alg_structures_dual ( Couplings<double>& couplings )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare LA structures for dual problem " << std::endl;
    SparsityStructure sparsity_dual;
    global_asm_.compute_sparsity_structure ( *this->space_dual_, sparsity_dual, &this->coupling_vars_dual_ );

    // Initialize linear algebra structures
    couplings.Init ( comm_, this->space_dual_->dof ( ) );
    couplings.InitializeCouplings ( sparsity_dual.off_diagonal_rows, sparsity_dual.off_diagonal_cols );

    // Initialize matrices and vectors
#ifdef USE_HYPRE
    matrix_dual_->Init ( comm_, couplings );
#else
    matrix_dual_->Init ( comm_, couplings, la_sys_.Platform, MetFlowApp::APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
#endif
    matrix_dual_->InitStructure ( vec2ptr ( sparsity_dual.diagonal_rows ),
                                  vec2ptr ( sparsity_dual.diagonal_cols ),
                                  sparsity_dual.diagonal_rows.size ( ),
                                  vec2ptr ( sparsity_dual.off_diagonal_rows ),
                                  vec2ptr ( sparsity_dual.off_diagonal_cols ),
                                  sparsity_dual.off_diagonal_rows.size ( ) );

#ifdef USE_HYPRE
    solD_ ->Init ( comm_, couplings );
    solD_next_->Init ( comm_, couplings );
    solD_prev_->Init ( comm_, couplings );
    rhs_dual_ ->Init ( comm_, couplings );
#else
    solD_ ->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    solD_ ->InitStructure ( );
    solD_next_->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    solD_next_->InitStructure ( );
    solD_prev_->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    solD_prev_->InitStructure ( );
    rhs_dual_ ->Init ( comm_, couplings, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    rhs_dual_ ->InitStructure ( );
#endif

    solD_ ->Zeros ( );
    solD_next_->Zeros ( );
    solD_prev_->Zeros ( );
    rhs_dual_ ->Zeros ( );

    // prepare dirichlet BC
    this->prepare_bc_dual ( );
    this->set_bc ( -1 );
}

void MetFlowApp::set_bc ( int mode )
{
    if ( mode == 1 )
    {
        if ( !dirichlet_dofs_.empty ( ) )
        {
            // correct solution with dirichlet BC
            solP_ ->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
            //solP_prev_->SetValues(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(), vec2ptr(dirichlet_values_));
            //solP_next_->SetValues(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(), vec2ptr(dirichlet_values_));
        }
        solP_ ->Update ( );
        //solP_prev_->Update();
        //solP_next_->Update();
    }
    else if ( mode == -1 )
    {
        if ( !dirichlet_dofs_dual_.empty ( ) )
        {
            // correct solution with dirichlet BC
            solD_ ->SetValues ( vec2ptr ( dirichlet_dofs_dual_ ), dirichlet_dofs_dual_.size ( ), vec2ptr ( dirichlet_values_dual_ ) );
            //solD_next_->SetValues(vec2ptr(dirichlet_dofs_dual_), dirichlet_dofs_dual_.size(), vec2ptr(dirichlet_values_dual_));
            //solD_prev_->SetValues(vec2ptr(dirichlet_dofs_dual_), dirichlet_dofs_dual_.size(), vec2ptr(dirichlet_values_dual_));
        }

        solD_ ->Update ( );
        //solD_next_->Update();
        //solD_prev_->Update();
    }
}

void MetFlowApp::prepare_nls ( bool stationary )
{
    if ( rank ( ) == master_rank ( ) )
    {
        if ( stationary )
            std::cout << "> Prepare nonlinear solver for stationary problem" << std::endl;
        else
            std::cout << "> Prepare nonlinear solver for instationary problem" << std::endl;
    }
    int maxits;
    double abstol;
    double reltol;
    double divtol = base_params_["PrimalNonlinearSolver"]["DivergenceLimit"].get<double>( );
    std::string damping_strat;

    if ( stationary )
    {
        maxits = base_params_["PrimalNonlinearSolver"]["Stationary"]["MaximumIterations"].get<int>( );
        abstol = base_params_["PrimalNonlinearSolver"]["Stationary"]["AbsoluteTolerance"].get<double>( );
        reltol = base_params_["PrimalNonlinearSolver"]["Stationary"]["RelativeTolerance"].get<double>( );
        use_pressure_filter_ = base_params_["PrimalNonlinearSolver"]["Stationary"]["PressureFilter"].get<bool>( );
        damping_strat = base_params_["PrimalNonlinearSolver"]["Stationary"]["DampingStrategy"].get<std::string>( );
    }
    else
    {
        maxits = base_params_["PrimalNonlinearSolver"]["TimeLoop"]["MaximumIterations"].get<int>( );
        abstol = base_params_["PrimalNonlinearSolver"]["TimeLoop"]["AbsoluteTolerance"].get<double>( );
        reltol = base_params_["PrimalNonlinearSolver"]["TimeLoop"]["RelativeTolerance"].get<double>( );
        use_pressure_filter_ = base_params_["PrimalNonlinearSolver"]["TimeLoop"]["PressureFilter"].get<bool>( );
        damping_strat = base_params_["PrimalNonlinearSolver"]["TimeLoop"]["DampingStrategy"].get<std::string>( );
    }

    std::string forcing_strat = base_params_["PrimalNonlinearSolver"]["ForcingStrategy"].get<std::string>( );
    int ew_type = base_params_["PrimalNonlinearSolver"]["EisenstatWalker"]["Type"].get<int>( );
    double max_forcing_term = base_params_["PrimalNonlinearSolver"]["EisenstatWalker"]["MaxForcingTerm"].get<double>( );
    double forcing_term = base_params_["PrimalNonlinearSolver"]["EisenstatWalker"]["InitialForcingTerm"].get<double>( );
    double alpha = base_params_["PrimalNonlinearSolver"]["EisenstatWalker"]["Alpha"].get<double>( );
    double gamma = base_params_["PrimalNonlinearSolver"]["EisenstatWalker"]["Gamma"].get<double>( );

    double initial = base_params_["PrimalNonlinearSolver"]["ArmijoUpdate"]["Initial"].get<double>( );
    double minimal = base_params_["PrimalNonlinearSolver"]["ArmijoUpdate"]["Minimal"].get<double>( );
    double decrease = base_params_["PrimalNonlinearSolver"]["ArmijoUpdate"]["Decrease"].get<double>( );
    double suff_dec = base_params_["PrimalNonlinearSolver"]["ArmijoUpdate"]["SuffDec"].get<double>( );
    int max_loop = base_params_["PrimalNonlinearSolver"]["ArmijoUpdate"]["MaxLoop"].get<int>( );

    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "> ForcingStrategy:      " << forcing_strat.c_str ( ) << ew_type << std::endl;
        std::cout << "  max forcing term:     " << max_forcing_term << std::endl;
        std::cout << "  initial forcing term: " << forcing_term << std::endl;
        std::cout << "  gamma:                " << gamma << std::endl;
        std::cout << "  alpha:                " << alpha << std::endl;

        std::cout << "> DampingStrategy: " << damping_strat.c_str ( ) << std::endl;
        std::cout << "  initial step:    " << initial << std::endl;
        std::cout << "  minimal step:    " << minimal << std::endl;
        std::cout << "  sufficient dec:  " << suff_dec << std::endl;
        std::cout << "  max loop:        " << max_loop << std::endl;
    }

    if ( forcing_strat == "EisenstatWalker" )
    {
        if ( ew_type == 1 )
        {
            this->forcing_strategy_ = new EWForcing<LAD>( forcing_term, max_forcing_term, 1 );
            this->nls_->SetForcingStrategy ( *this->forcing_strategy_ );
        }
        else if ( ew_type == 2 )
        {
            this->forcing_strategy_ = new EWForcing<LAD>( forcing_term, max_forcing_term, 2, gamma, alpha );
            this->nls_->SetForcingStrategy ( *this->forcing_strategy_ );
        }
    }

    if ( damping_strat == "Armijo" )
    {
        this->damping_strategy_ = new ArmijoDamping<LAD>( initial, minimal, decrease, suff_dec, max_loop );
        this->nls_->SetDampingStrategy ( *this->damping_strategy_ );
    }

    nls_->InitControl ( maxits, abstol, reltol, divtol );
    nls_->SetLinearSolver ( *linear_solver_ );
    nls_->SetOperator ( *this );

    std::string path = this->root_ + this->filename_linear_solver_;
    nls_->SetOutputFilename ( path );

    // we use our own initial solution -- this needs to be indicated to the Newton-solver
    nls_->MyNewton<LAD>::InitParameter ( MyNewton<LAD>::NewtonInitialSolutionOwn );
    //  nls_->MyNewton<LAD>::InitParameter(MyNewton<LAD>::NewtonInitialSolution0);
}

int MetFlowApp::init_run ( bool run_mode0 )
{
    int time_step = 0;

    // For HDF5
    std::string filename_primal = base_params_["BackUp"]["PrimalSnapshotsIn"].get<std::string>( );
    std::string groupname_primal = base_params_["BackUp"]["PrimalSnapshotsGroupIn"].get<std::string>( );
    std::string prefix_primal = base_params_["BackUp"]["PrimalSnapshotsPrefixIn"].get<std::string>( );
    int prefix_primal_with_timestep = base_params_["BackUp"]["PrimalSnapshotsTimeStepIn"].get<int>( 1 );
    filename_primal = this->root_ + "/" + filename_primal;
    this->problem_mode_ = 1;

    if ( ( base_params_["BackUp"].contains ( "PrimalOffset" ) || this->resumed_run_ ) && ( run_mode0 == false ) )
    {
        // resume time loop
        int offset = 0;

        // get last backup step
        if ( this->resumed_run_ )
        {
            ifstream in;
            string path = this->root_ + "/log/LastBackupStep.txt";
            in.open ( path.c_str ( ), ios::in );
            in >> offset;
            in.close ( );
        }
        else
        {
            offset = base_params_["BackUp"]["PrimalOffset"].get<int>( );
        }

        interminable_assert ( offset >= 0 );

        // PrimalOffset is the index of the last existing file
        time_step = 1 + offset;

        // setup correct fe space and mesh
        this->dmh_->update ( time_step, 1, true );
        this->prepare_linear_solver ( true );

        // read last primal solution from file
        std::stringstream ss;
        ss << prefix_primal;
        if ( prefix_primal_with_timestep == 1 )
        {
            ss << offset;
        }

        // last solution belongs to same fe space as current one
        this->read_file ( *this->solP_prev_, offset, filename_primal, groupname_primal, ss.str ( ) );

        // read base state from file
        std::stringstream ss_base;
        std::string groupname_start = "start";
        ss_base << "start";
        this->create_ic_hdf5_filename ( );

        //this->read_file(*this->base_, this->filename_base_, groupname_start, ss_base.str());
        this->base_->Zeros ( );

        this->local_asm_primal_->set_base ( *this->base_ );

        this->cur_time_ = 0.;
        for ( int l = 0; l <= offset; ++l )
        {
            this->cur_time_ += this->get_delta_t ( l );
        }

        this->solP_->CloneFrom ( *this->solP_prev_ );
        this->local_asm_primal_->set_time ( this->cur_time_ );

        if ( rank ( ) == master_rank ( ) ) std::cout << "> Resuming Simulation at step / time " << time_step << " / " << this->cur_time_ << std::endl;

        this->delta_t_ = this->get_delta_t ( time_step );
        this->cur_time_ += this->delta_t_;
        this->time_step_ = time_step;
        if ( rank ( ) == master_rank ( ) )
        {
            this->create_log_files ( offset + 1 );
        }

        //this->visualize_solution(9999);
        //exit(-1);
    }
    else
    {
        // first time step
        time_step = 0;
        this->time_step_ = time_step;
        this->cur_time_ = 0.0;
        this->local_asm_primal_->set_time ( this->cur_time_ );

        // create empty files for flow characteristics
        if ( rank ( ) == master_rank ( ) )
        {
            this->create_log_files ( 0 );
        }

        std::string starting_type = base_params_["InitialCondition"]["Type"].get<std::string>( );

        if ( rank ( ) == master_rank ( ) ) std::cout << "> Starting Simulation at step " << time_step << std::endl;

        // setup correct fe space
        this->dmh_->update ( time_step, 1, true );
        this->prepare_linear_solver ( true );

        // Compute initial condition
        if ( !this->is_stationary_ )
        {
            this->prepare_ic ( );
        }
        else
        {
            this->solP_->Zeros ( );
        }

        this->solP_prev_->CloneFrom ( *this->solP_ );

        this->solP_prev_->Update ( );
        this->solP_->Update ( );

        this->local_asm_primal_->set_base ( *this->base_ );
        this->local_asm_primal_->set_solP ( *this->solP_ );
        this->local_asm_primal_->set_solP_prev ( *this->solP_prev_ );

        if ( !this->is_stationary_ )
        {
            this->visualize_solution ( time_step );
        }

        // Write solution
        if ( starting_type != "Load" )
        {
            this->write_ic ( );
        }

        // Post Processing
        this->post_processing ( );

        time_step = 1;
        this->time_step_ = time_step;
        this->delta_t_ = this->get_delta_t ( time_step );
        this->cur_time_ = this->delta_t_;
    }
    return time_step;
}

int MetFlowApp::init_run_dual ( )
{
    int time_step = 0;
    this->problem_mode_ = -1;

    // For HDF5
    std::string filename_primal = base_params_["BackUp"]["PrimalSnapshotsIn"].get<std::string>( );
    std::string groupname_primal = base_params_["BackUp"]["PrimalSnapshotsGroupIn"].get<std::string>( );
    std::string prefix_primal = base_params_["BackUp"]["PrimalSnapshotsPrefixIn"].get<std::string>( );
    int prefix_primal_with_timestep = base_params_["BackUp"]["PrimalSnapshotsTimeStepIn"].get<int>( 1 );

    filename_primal = this->root_ + "/" + filename_primal;

    std::string filename_dual = base_params_["BackUp"]["DualSnapshotsIn"].get<std::string>( );
    std::string groupname_dual = base_params_["BackUp"]["DualSnapshotsGroupIn"].get<std::string>( );
    std::string prefix_dual = base_params_["BackUp"]["DualSnapshotsPrefixIn"].get<std::string>( );
    int prefix_dual_with_timestep = base_params_["BackUp"]["DualSnapshotsTimeStepIn"].get<int>( 1 );

    std::vector<std::string> datasetname_dual;
    filename_dual = this->root_ + "/" + filename_dual;

    if ( base_params_["BackUp"].contains ( "DualOffset" ) || this->resumed_run_ )
    {
        // ************************************************************
        // Read in most recent dual solution
        int offset = 0;
        if ( this->resumed_run_ )
        {
            ifstream in;
            string path = this->root_ + "/log/LastBackupStep.txt";
            in.open ( path.c_str ( ), ios::in );
            in >> offset;
            in.close ( );
        }
        else
        {
            offset = base_params_["BackUp"]["DualOffset"].get<int>( );
        }
        interminable_assert ( offset >= 0 );

        // PrimalOffset is the index of the last existing file
        time_step = offset - 1;

        // setup correct mesh and fe space
        this->dmh_->update ( time_step, -1, true );
        this->prepare_linear_solver_dual ( );

        // read last primal and dual solution from file
        std::stringstream ss;
        ss << prefix_dual;
        if ( prefix_dual_with_timestep == 1 )
        {
            ss << offset;
        }

        this->read_file ( *this->solD_, offset, filename_dual, groupname_dual, ss.str ( ) );

        this->create_log_files_dual ( offset - 1 );

        // Read last primal solution
        if ( method_dual_ == "GalerkinDC" )
        {
            double time;
            std::stringstream ss_primal;
            ss_primal << prefix_primal;
            if ( prefix_primal_with_timestep == 1 )
            {
                ss_primal << offset;
            }
            if ( rank ( ) == master_rank ( ) )
                std::cout << "  Reading primal solution at time step " << offset << " to initialize dual solution.\n";

            // last solution belongs to same fe space as current one
            this->read_file ( *this->solP_, offset, filename_primal, groupname_primal, ss_primal.str ( ) );

            if ( offset >= 1 )
            {
                std::stringstream ss_primal_prev;
                ss_primal_prev << prefix_primal;
                if ( prefix_primal_with_timestep == 1 )
                {
                    ss_primal_prev << offset - 1;
                }

                if ( rank ( ) == master_rank ( ) )
                    std::cout << "  Reading primal solution at time step " << offset - 1 << " to initialize dual solution.\n";

                this->read_file ( *this->solP_prev_, filename_primal, groupname_primal, ss_primal_prev.str ( ) );
            }
        }
        if ( method_dual_ == "GalerkinCD" || method_dual_ == "Theta" )
        {
            double time;
            std::stringstream ss_primal;
            ss_primal << prefix_primal;
            if ( prefix_primal_with_timestep == 1 )
            {
                ss_primal << offset;
            }

            if ( rank ( ) == master_rank ( ) )
                std::cout << "  Reading primal solution at time step " << offset << " to initialize dual solution.\n";

            // last solution belongs to same fe space as current one
            this->read_file ( *this->solP_, offset, filename_primal, groupname_primal, ss_primal.str ( ) );
            this->solP_prev_->Zeros ( );
        }

        this->cur_time_ = 0.;
        for ( int l = 0; l <= offset; ++l )
        {
            this->cur_time_ += this->get_delta_t ( l );
        }

        this->local_asm_dual_->set_time ( this->cur_time_ );

        if ( rank ( ) == master_rank ( ) )
            std::cout << "> Resuming Dual Simulation after step / time " << time_step + 1 << " / " << this->cur_time_ << std::endl;

        // possibly change delta_t_
        this->delta_t_ = this->get_delta_t ( time_step );
        this->delta_t_next_ = this->get_delta_t ( time_step + 1 );
        this->cur_time_ -= this->delta_t_;
        this->time_step_ = time_step;
    }
    else
    {
        // ************************************************************
        // Compute initial dual solution

        time_step = this->get_num_intervals ( );
        this->cur_time_ = this->duration_;
        this->local_asm_dual_->set_time ( this->cur_time_ );

        std::string starting_type = base_params_["InitialCondition"]["Type"].get<std::string>( );
        this->create_log_files_dual ( 0 );

        if ( rank ( ) == master_rank ( ) ) std::cout << "> Starting Dual Simulation at step " << time_step << std::endl;

        // setup correct mesh and fe space
        this->dmh_->update ( time_step, -1, true );
        this->prepare_linear_solver_dual ( );

        // Read last primal solution
        if ( method_dual_ == "GalerkinDC" )
        {
            double time;
            std::stringstream ss_primal_prev;
            ss_primal_prev << prefix_primal;
            if ( prefix_primal_with_timestep == 1 )
            {
                ss_primal_prev << this->get_num_intervals ( );
            }

            if ( rank ( ) == master_rank ( ) )
                std::cout << "  Reading primal solution at time step " << this->get_num_intervals ( ) << " to initialize dual solution.\n";

            this->read_file ( *this->solP_prev_, filename_primal, groupname_primal, ss_primal_prev.str ( ) );
            this->solP_->Zeros ( );
        }
        if ( method_dual_ == "GalerkinCD" || method_dual_ == "Theta" || method_dual_ == "Stationary" )
        {
            double time;
            int offset = 0;
            if ( this->is_stationary_ )
            {
                offset = -1;
            }
            std::stringstream ss_primal;
            ss_primal << prefix_primal;
            if ( prefix_primal_with_timestep == 1 )
            {
                ss_primal << this->get_num_intervals ( ) + offset;
            }

            if ( rank ( ) == master_rank ( ) )
                std::cout << "  Reading primal solution at time step " << this->get_num_intervals ( ) << " to initialize dual solution.\n";

            this->read_file ( *this->solP_, this->get_num_intervals ( ) + offset, filename_primal, groupname_primal, ss_primal.str ( ), 1 );
        }

        // Compute IC
        int off = 0;
        if ( method_dual_ == "GalerkinCD" || method_dual_ == "Theta" )
        {
            this->prepare_ic_dual ( );
            this->solD_->CloneFrom ( *this->solD_next_ );
        }
        if ( method_dual_ == "GalerkinDC" )
        {
            off = 1;
            this->solD_->CloneFrom ( *this->solP_prev_ );
            this->solD_->Scale ( 0. );

            this->solD_->Update ( );
            this->solP_->Update ( );
        }
        if ( method_dual_ == "Stationary" )
        {
            this->solD_next_->Zeros ( );
            this->solD_->Zeros ( );
        }

        // Post Processing
        this->post_processing_dual ( );

        // Visualization
        if ( !this->is_stationary_ )
        {
            this->visualize_solution_dual ( time_step + off );
        }

        std::stringstream ss;
        ss << prefix_dual;
        if ( prefix_dual_with_timestep == 1 )
        {
            ss << time_step + off;
        }

        if ( !this->is_stationary_ )
        {
            this->write_file ( *this->solD_, filename_dual, groupname_dual, ss.str ( ) );
        }

        time_step = time_step - 1 + off;
        this->delta_t_ = this->get_delta_t ( time_step );

        if ( method_dual_ == "GalerkinDC" )
            this->delta_t_next_ = 0.;
        else
            this->delta_t_next_ = this->get_delta_t ( time_step + 1 );

        this->local_asm_dual_->set_time ( this->cur_time_ );
        this->local_asm_dual_->set_dT_pc ( this->delta_t_ );
        this->local_asm_dual_->set_dT_cn ( this->delta_t_next_ );
    }

    if ( method_dual_ == "GalerkinDC" )
    {
        this->max_step_dual_ += 1;
    }
    return time_step;
}

void MetFlowApp::prepare_assembler ( )
{
    this->local_asm_primal_->set_mode_to_primal ( );
    this->local_asm_primal_->set_solP ( *this->solP_ );
    this->local_asm_primal_->set_solP_prev ( *this->solP_prev_ );
    this->local_asm_primal_->set_perturb ( *this->perturb_ );
    this->local_asm_primal_->set_perturb_prev ( *this->perturb_prev_ );
    this->local_asm_primal_->set_base ( *this->base_ );

    this->local_asm_primal_->set_dT_pc ( delta_t_ );

    if ( this->is_stationary_ )
    {
        this->local_asm_primal_->set_time_discretization_off ( );
    }
    else
    {
        if ( method_ == "GalerkinCD" )this->local_asm_primal_->set_time_discretization_to_galerkin_cd ( mod_galerkin_ );
        if ( method_ == "Theta" ) this->local_asm_primal_->set_time_discretization_to_theta ( this->theta_ );
    }
}

void MetFlowApp::prepare_assembler_dual ( )
{
    this->local_asm_dual_->set_solD ( *this->solD_ );
    this->local_asm_dual_->set_solD_prev ( *this->solD_prev_ );
    this->local_asm_dual_->set_solD_next ( *this->solD_next_ );
    this->local_asm_dual_->set_solP ( *this->solP_ );
    this->local_asm_dual_->set_solP_prev ( *this->solP_prev_ );
    this->local_asm_dual_->set_solP_next ( *this->solP_next_ );

    this->local_asm_dual_->set_mode_to_dual ( );

    this->local_asm_dual_->set_dT_pc ( delta_t_ );
    this->local_asm_dual_->set_dT_cn ( delta_t_ );

    if ( this->is_stationary_ == 1 )
    {
        this->local_asm_dual_->set_time_discretization_off ( );
    }
    else
    {
        if ( method_dual_ == "GalerkinCD" ) this->local_asm_dual_->set_time_discretization_to_galerkin_cd ( false );
        if ( method_dual_ == "GalerkinDC" ) this->local_asm_dual_->set_time_discretization_to_galerkin_dc ( false );
        if ( method_dual_ == "Theta" ) this->local_asm_dual_->set_time_discretization_to_theta ( this->theta_dual_ );
        if ( method_dual_ == "Simple" ) this->local_asm_dual_->set_time_discretization_to_simple ( );
    }

}

void MetFlowApp::prepare_assembler_est ( )
{
    this->local_asm_est_->set_solP ( *this->solP_ );
    this->local_asm_est_->set_solP_prev ( *this->solP_prev_ );
    this->local_asm_est_->set_solP_next ( *this->solP_next_ );

    this->local_asm_est_->set_solD ( *this->solD_ );
    this->local_asm_est_->set_solD_next ( *this->solD_next_ );
    this->local_asm_est_->set_solD_prev ( *this->solD_prev_ );

    this->local_asm_est_->set_fineP ( this->fineP_ );
    this->local_asm_est_->set_fineP_prev ( this->fineP_prev_ );
    this->local_asm_est_->set_fineP_next ( this->fineP_next_ );

    this->local_asm_est_->set_fineD ( this->fineD_ );
    this->local_asm_est_->set_fineD_prev ( this->fineD_prev_ );
    this->local_asm_est_->set_fineD_next ( this->fineD_next_ );

    this->local_asm_est_->set_fineP_space ( *this->fine_space_ );
    this->local_asm_est_->set_fineP_space_next ( *this->fine_space_ );
    this->local_asm_est_->set_fineP_space_prev ( *this->fine_space_ );

    this->local_asm_est_->set_fineD_space ( *this->fine_space_dual_ );
    this->local_asm_est_->set_fineD_space_next ( *this->fine_space_dual_ );
    this->local_asm_est_->set_fineD_space_prev ( *this->fine_space_dual_ );

    int quad_order;
    if ( this->is_stationary_ )
    {
        quad_order = 0;
    }
    else
    {
        quad_order = base_params_["Adaptivity"]["Estimator"]["TimeQuadOrder"].get<int>( 1 );
    }
    this->local_asm_est_->setup_time_quadrature ( quad_order );

    bool primal = base_params_["Adaptivity"]["Estimator"]["PrimalIndicator"].get<int>( 1 );
    bool dual = base_params_["Adaptivity"]["Estimator"]["DualIndicator"].get<int>( 1 );
    bool temporal = base_params_["Adaptivity"]["Estimator"]["TemporalIndicator"].get<int>( 1 );
    bool spatial = base_params_["Adaptivity"]["Estimator"]["SpatialIndicator"].get<int>( 1 );
    bool csi = base_params_["Adaptivity"]["Estimator"]["UseCauchySchwarz"].get<int>( 0 );
    bool no_dwr = base_params_["Adaptivity"]["Estimator"]["NoDWR"].get<int>( 0 );
    this->local_asm_est_->set_indicators ( primal, dual, temporal, spatial, csi );

    if ( no_dwr )
    {
        this->local_asm_est_->use_dwr ( false );
    }
}

/// ***************************************************************************
/// DynamicMeshProblem
/// ***************************************************************************

void MetFlowApp::set_active_space ( VectorSpace<DATATYPE>* space, int mode )
{
    if ( mode == 1 )
    {
        this->space_ = space;
    }
    else if ( mode == -1 )
    {
        this->space_dual_ = space;
    }
}

void MetFlowApp::set_active_mesh ( MeshPtr mesh )
{
    this->mesh_ = mesh;
}

void MetFlowApp::init_mesh_change_list ( )
{
    std::vector<DATATYPE> mesh_change_times;
    if ( this->adapt_counter_ < 0 )
    {
        std::string start_type = base_params_["Adaptivity"]["DynamicMesh"]["StartChanges"].get<std::string>( );
        if ( start_type == "Load" )
        {
            this->read_mesh_change_list ( -1, mesh_change_times );
        }
        this->write_mesh_change_list ( 0, mesh_change_times );
    }
    else
    {
        this->read_mesh_change_list ( this->adapt_counter_, mesh_change_times );
    }

    this->dmh_->set_initial_mesh_change_times ( mesh_change_times );
}

void MetFlowApp::set_update_vectors ( )
{
    std::vector<VectorType*> p_vectors;
    std::vector<VectorType*> d_p_vectors;
    std::vector<VectorType*> d_vectors;

    p_vectors.push_back ( this->solP_prev_ );
    d_p_vectors.push_back ( this->solP_prev_ );
    d_p_vectors.push_back ( this->solP_next_ );
    d_vectors.push_back ( this->solD_prev_ );
    d_vectors.push_back ( this->solD_next_ );

    this->dmh_->set_update_vectors ( 1, 1, p_vectors );
    this->dmh_->set_update_vectors ( -1, 1, d_p_vectors );
    this->dmh_->set_update_vectors ( -1, -1, d_vectors );
    this->dmh_->set_update_vectors ( 0, 1, d_p_vectors );
    this->dmh_->set_update_vectors ( 0, -1, d_vectors );
}

void MetFlowApp::setup_LA_primal ( VectorSpace<DATATYPE>& space, VectorSpace<DATATYPE>& space_dual,
                                   std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                   Couplings<DATATYPE>& couplings, Couplings<DATATYPE>& couplings_dual )
{
    if ( rank ( ) == master_rank ( ) )
        std::cout << "> Setup primal LA objects in dmh->update() " << std::endl;

    this->space_ = &space;
    this->coupling_vars_ = coupling_vars_;
    this->prepare_lin_alg_structures ( couplings );

}

void MetFlowApp::setup_LA_dual ( VectorSpace<DATATYPE>& space, VectorSpace<DATATYPE>& space_dual,
                                 std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                 Couplings<DATATYPE>& couplings, Couplings<DATATYPE>& couplings_dual )
{
    if ( rank ( ) == master_rank ( ) )
        std::cout << "> Setup dual LA objects in dmh->update() " << std::endl;

    this->space_ = &space;
    this->space_dual_ = &space_dual;
    this->coupling_vars_ = coupling_vars_;
    this->coupling_vars_dual_ = coupling_vars_dual_;
    this->prepare_lin_alg_structures ( couplings );
    this->prepare_lin_alg_structures_dual ( couplings_dual );
}

void MetFlowApp::setup_LA_est ( VectorSpace<DATATYPE>& space, VectorSpace<DATATYPE>& space_dual,
                                std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                Couplings<DATATYPE>& couplings, Couplings<DATATYPE>& couplings_dual )
{
    if ( rank ( ) == master_rank ( ) )
        std::cout << "> Setup estimator LA objects in dmh->update() " << std::endl;

    this->space_ = &space;
    this->space_dual_ = &space_dual;
    this->coupling_vars_ = coupling_vars_;
    this->coupling_vars_dual_ = coupling_vars_dual_;
    this->prepare_lin_alg_structures ( couplings );
    this->prepare_lin_alg_structures_dual ( couplings_dual );

}

void MetFlowApp::setup_space ( VectorSpace<DATATYPE>& space, std::vector< std::vector< bool> >& coupling_vars, MeshPtr mesh, int mode )
{
    if ( rank ( ) == master_rank ( ) )
        std::cout << "> Setup space in dmh->init_fe_spaces() for mode " << mode << std::endl;

    this->prepare_space ( space, coupling_vars, mesh, mode );
    this->periodify_space ( space, mesh );
}

/// ***************************************************************************
/// SOLVER
/// ***************************************************************************

int MetFlowApp::solve_primal ( int adapt, int time_step )
{
    std::string filename_primal = base_params_["BackUp"]["PrimalSnapshotsOut"].get<std::string>( );
    std::string groupname_primal = base_params_["BackUp"]["PrimalSnapshotsGroupOut"].get<std::string>( );
    std::string prefix_primal = base_params_["BackUp"]["PrimalSnapshotsPrefixOut"].get<std::string>( );
    int prefix_primal_with_timestep = base_params_["BackUp"]["PrimalSnapshotsTimeStepOut"].get<int>( 1 );

    filename_primal = this->root_ + "/" + filename_primal;

    // setup correct fe space (including linear sovler)
    bool updated = this->dmh_->update ( time_step, 1, false );
    this->prepare_linear_solver ( false );

    this->matrix_assembled_ = false;

    int off = 0;
    if ( this->is_stationary_ )
    {
        off = -1;
    }

    if ( time_step == 1 && !this->is_stationary_ )
    {
        std::stringstream ss;
        ss << prefix_primal;
        if ( prefix_primal_with_timestep == 1 )
        {
            ss << time_step - 1;
        }
        this->write_file ( *this->solP_prev_, filename_primal, groupname_primal, ss.str ( ) );
    }

    // setup nonlinear solver
    if ( !this->is_linear_problem_ )
    {
        this->prepare_nls ( false );
    }

    // ********************************************************
    // Time loop
    // ********************************************************
    if ( rank ( ) == master_rank ( ) )
    {
        this->local_asm_primal_->print_parameters ( );
    }
    while ( time_step <= this->get_num_intervals ( ) )
    {
        this->cur_time_ = this->get_time ( time_step );
        this->delta_t_ = this->get_delta_t ( time_step );
        int current_mesh_index = this->dmh_->mesh_index ( time_step );

        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << std::endl << "======================================================" << std::endl;
            std::cout << "> Adaption Cycle " << this->adapt_counter_ << std::endl;
            std::cout << "> Time step      " << time_step << " / " << this->get_num_intervals ( ) << std::endl;
            std::cout << "> Time           " << this->get_time ( time_step ) << " / " << duration_ << std::endl;
            std::cout << "> delta_T        " << this->delta_t_ << std::endl;
            std::cout << "> Current Mesh index " << current_mesh_index << std::endl;
            std::cout << "> Active mesh index  " << this->active_mesh_index_ << std::endl;
        }

        // update mesh and fe space
        updated = this->dmh_->update ( time_step, 1, false );
        if ( updated )
        {
            this->prepare_linear_solver ( false );
        }

        // General update of assembler
        this->local_asm_primal_->set_dT_pc ( this->delta_t_ );
        this->local_asm_primal_->set_time ( this->cur_time_ );

        // Problem specific update of assembler
        this->update_assembler ( );

        // Update BOundary conditions in case of varying voltage
        this->prepare_bc ( );
        this->set_bc ( 1 );

        // perturb solution
        this->perturb_solution ( time_step );

        // solve stationary problem
        bool success;
        if ( rank ( ) == master_rank ( ) ) std::cout << "> Solve stationary problem " << std::endl;

        if ( !this->is_linear_problem_ )
        {
            success = solve_nlp ( 1 );
        }
        else
        {
            success = solve_lp ( 1 );
        }

        // apply filter to solution
        this->filter_solution ( );

        if ( time_step % pp_step_ == 0 || this->is_stationary_ )
        {
            this->post_processing ( );
        }
        if ( time_step % visual_step_ == 0 || this->is_stationary_ )
        {
            this->visualize_solution ( time_step + off );
        }

        // Write solution every n-th time step
        if ( ( time_step % backup_step_ == 0 ) || ( time_step == this->get_num_intervals ( ) ) || this->is_stationary_ )
        {
            std::stringstream ss;
            ss << prefix_primal;
            if ( prefix_primal_with_timestep == 1 )
            {
                ss << time_step + off;
            }
            this->write_file ( *this->solP_, filename_primal, groupname_primal, ss.str ( ) );

            ofstream out;
            string path = this->root_ + "/log/LastBackupStep.txt";
            out.open ( path.c_str ( ), ios::out );
            out << time_step << "\n ";
            out.close ( );
        }

        if ( success == false )
        {
            std::cerr << "No solution found for nonlinear problem. \nProgram terminates now.\n";
            interminable_assert ( 0 );
        }

        time_step++;
        this->time_step_ = time_step;
        this->solP_prev_ ->CloneFrom ( *this->solP_ );
    }
    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "END OF PRIMAL TIME LOOP" << std::endl;
    }
    return time_step - 1;
}

int MetFlowApp::solve_dual ( int adapt, int time_step )
{
    // For HDF5
    std::string filename_primal = base_params_["BackUp"]["PrimalSnapshotsIn"].get<std::string>( );
    std::string groupname_primal = base_params_["BackUp"]["PrimalSnapshotsGroupIn"].get<std::string>( );
    std::string prefix_primal = base_params_["BackUp"]["PrimalSnapshotsPrefixIn"].get<std::string>( );
    int prefix_primal_with_timestep = base_params_["BackUp"]["PrimalSnapshotsTimeStepIn"].get<int>( 1 );
    filename_primal = this->root_ + "/" + filename_primal;

    std::string filename_dual = base_params_["BackUp"]["DualSnapshotsOut"].get<std::string>( );
    std::string groupname_dual = base_params_["BackUp"]["DualSnapshotsGroupOut"].get<std::string>( );
    std::string prefix_dual = base_params_["BackUp"]["DualSnapshotsPrefixOut"].get<std::string>( );
    int prefix_dual_with_timestep = base_params_["BackUp"]["DualSnapshotsTimeStepOut"].get<int>( 1 );

    std::vector<std::string> datasetname_dual;
    filename_dual = this->root_ + "/" + filename_dual;

    int step_offset = 0;
    this->min_step_dual_ = 0;
    if ( method_dual_ == "GalerkinDC" )
    {
        this->min_step_dual_ = 1;
        step_offset = 0;
    }

    // ************************************************************
    // setup correct fe space (including linear solver)
    bool updated = this->dmh_->update ( time_step - step_offset, -1, false );
    this->prepare_linear_solver_dual ( );

    // ************************************************************
    // Dual Time Stepping Loop
    // NOTE: Theta and cGdG-> time_steps refers to specific point in time, dGcG-> time_step refers to interval
    // NOTE: delta_t_         = t_n - t_{n-1}
    // NOTE: delta_t_next     = t_{n+1} -t_n
    // NOTE: cGdG: t_n = t_{time_step}, use dt = delta_t_next only
    // NOTE: dGcG: t_n = t_{time_step}, use dt = delta_t and dt = delta_t_next
    // NOTE: cur_time = t_n

    this->dual_matrix_assembled_ = false;

    if ( rank ( ) == master_rank ( ) )
    {
        this->local_asm_dual_->print_parameters ( );
    }
    while ( time_step >= this->min_step_dual_ )
    {
        this->delta_t_next_ = this->get_delta_t ( time_step + 1 );
        this->delta_t_ = this->get_delta_t ( time_step );
        this->cur_time_ = this->get_time ( time_step );
        int current_mesh_index = this->dmh_->mesh_index ( time_step - step_offset );

        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << std::endl << "======================================================" << std::endl;
            std::cout << "> Adaption Cycle " << this->adapt_counter_ << std::endl;
            if ( method_dual_ == "GalerkinCD" || method_dual_ == "Theta" )
            {
                std::cout << "> Time step       " << time_step << " / " << this->get_num_intervals ( ) << std::endl;
                std::cout << "> Time            " << this->cur_time_ << " / " << duration_ << std::endl;
            }
            if ( method_dual_ == "GalerkinDC" )
            {
                std::cout << "> Time step         " << time_step << " / " << this->get_num_intervals ( ) << std::endl;
                std::cout << "> Time interval     " << time_step << ":  (" << this->cur_time_ - this->delta_t_ << ", " << this->cur_time_ << ")" << std::endl;
                std::cout << "> dT_pc             " << this->delta_t_ << ", dT_cn " << this->delta_t_next_ << std::endl;
            }
            std::cout << "> Current mesh index " << current_mesh_index << std::endl;
            std::cout << "> Active mesh index  " << this->active_mesh_index_ << std::endl;
        }

        // pass already computed time steps to assembler
        this->solD_next_->CloneFrom ( *this->solD_ );
        this->solP_next_->CloneFrom ( *this->solP_ );

        // update mesh and fe space
        updated = this->dmh_->update ( time_step - step_offset, -1, false );
        if ( updated )
        {
            this->prepare_linear_solver_dual ( );
        }

        if ( method_dual_ == "GalerkinDC" )
        {
            this->solP_->CloneFrom ( *this->solP_prev_ );

            if ( time_step >= 1 )
            {
                std::stringstream ss_primal;
                ss_primal << prefix_primal;
                if ( prefix_primal_with_timestep == 1 )
                {
                    ss_primal << time_step - 1;
                }
                double time;
                this->read_file ( *this->solP_prev_, time_step - 1, filename_primal, groupname_primal, ss_primal.str ( ), 1 );
            }

            if ( time_step == this->get_num_intervals ( ) )
            {
                // In this case, the final condition for the dual problem is part of the weak formulation
                this->local_asm_dual_->set_final_time_goal_contrib ( true );
            }
            else
            {
                this->local_asm_dual_->set_final_time_goal_contrib ( false );
            }
        }
        if ( method_dual_ == "GalerkinCD" || method_dual_ == "Theta" || this->is_stationary_ )
        {
            std::stringstream ss_primal;
            ss_primal << prefix_primal;
            if ( prefix_primal_with_timestep == 1 )
            {
                ss_primal << time_step;
            }

            this->read_file ( *this->solP_, time_step, filename_primal, groupname_primal, ss_primal.str ( ), 1 );
            this->solP_prev_->Zeros ( );
        }

        // general update of assembler
        this->local_asm_dual_->set_dT_pc ( this->delta_t_ );
        this->local_asm_dual_->set_dT_cn ( this->delta_t_next_ );
        this->local_asm_dual_->set_time ( this->cur_time_ );

        // problem specific update of assembler
        this->update_assembler_dual ( );

        bool success;
        if ( rank ( ) == master_rank ( ) ) std::cout << "> Solve dual linear problem " << std::endl;
        success = MetFlowApp::solve_lp ( -1 );

        // Post processing
        if ( time_step % pp_step_ == 0 || this->is_stationary_ )
        {
            this->post_processing_dual ( );
        }

        // Visualize solution
        if ( ( time_step ) % visual_step_ == 0 || this->is_stationary_ )
            visualize_solution_dual ( time_step );

        // Write solution
        if ( ( time_step % backup_step_ == 0 ) || ( time_step == 0 ) || this->is_stationary_ )
        {
            std::stringstream ss;
            ss << prefix_dual;
            if ( prefix_dual_with_timestep == 1 )
            {
                ss << time_step;
            }
            this->write_file ( *this->solD_, filename_dual, groupname_dual, ss.str ( ) );

            ofstream out;
            string path = this->root_ + "/log/LastBackupStep.txt";
            out.open ( path.c_str ( ), ios::out );
            out << time_step << "\n ";
            out.close ( );
        }

        if ( success == false )
        {
            std::cerr << "No solution found for linear problem. \nProgram terminates now.\n";
            interminable_assert ( 0 );
        }

        // possibly change delta_t_
        --time_step;
        this->time_step_ = time_step;
    }
    return 0;
}

bool MetFlowApp::solve_nlp ( int mode )
{
    assert ( mode == 1 );

    Timer timer;
    timer.start ( );

    this->problem_mode_ = mode;
    NonlinearSolverState state;
    int iter = 0;
    double res = 0.;
    double duration = 0.;

    if ( mode == 1 )
    {
        state = nls_->Solve ( *rhs_, solP_, *space_ );
        iter = nls_->iter ( );
        res = nls_->GetResidual ( );
    }
    /*
    else if (mode == -1)
    {
        state = nls_dual_->Solve(*rhs_dual_, solD_, *space_dual_);
        iter = nls_dual_->iter();
        res = nls_dual_->GetResidual();
    }
     * */
    timer.stop ( );
    duration = timer.get_duration ( );

    if ( rank_ == master_rank_ )
    {
        std::cout << "  Nonlinear solver ended with state " << state
                << " and residual norm " << res
                << " after " << iter << " iterations\n";

        string path = this->root_ + this->filename_nonlinear_solver_;
        ofstream out;
        out.open ( path.c_str ( ), ios::out | ios::app );
        out.precision ( 6 );
        out << std::scientific;
        out << this->cur_time_ << " " << iter << " " << duration << " " << res << "\n";
        out.close ( );
    }

    // interminable_assert(state == 0);
    if ( rank_ == master_rank_ )
    {
        std::cout << "  Time for solution of nonlinear problem: " << duration << std::endl << std::endl;
    }

    // Resetting forcing term and solution vector
    if ( mode == 1 )
    {
        this->local_asm_primal_->set_solP ( *this->solP_ );

        if ( base_params_["PrimalNonlinearSolver"]["ForcingStrategy"].get<std::string>( ) == "EisenstatWalker" )
        {
            double forcing_term = ( base_params_["PrimalNonlinearSolver"]["EisenstatWalker"]["InitialForcingTerm"].get<double>( ) );
            nls_->SetForcingTerm ( forcing_term );
            // Note: Changed from "nls_->ResetForcingTerm(forcing_term)"
        }
    }
    /*
    else if (mode == -1)
    {
        if (base_params_["DualNonlinearSolver"]["ForcingStrategy"].get<std::string>("")=="EisenstatWalker1")
        {
            double forcing_term = (base_params_["DualNonlinearSolver"]["ConstantForcingTerm"].get<double>(1e-3));
            nls_dual_->SetForcingTerm(forcing_term);
            // Note: Changed from "nls_->ResetForcingTerm(forcing_term)"
        }
    }
     */

    if ( state != 0 )
        return false;
    else
        return true;
}

bool MetFlowApp::solve_lp ( int mode )
{

    bool state = true;
    Timer timer;
    timer.start ( );

    LAD::MatrixType* matrix;
    LAD::VectorType* rhs;
    LAD::VectorType* sol;
    VectorSpace<double>* space;
    LinearSolver<LAD>* linear_solver;

    bool use_pressure_filter;

    this->solP_->Update ( );
    this->solP_prev_->Update ( );

    if ( mode == 1 )
    {
        matrix = this->matrix_;
        rhs = this->rhs_;
        sol = this->solP_;
        linear_solver = this->linear_solver_;
        space = this->space_;
        use_pressure_filter_ = base_params_["PrimalLinearSolver"]["OuterPressureFilter"].get<int>( 0 );
        use_pressure_filter = use_pressure_filter_;
    }
    else if ( mode == -1 )
    {
        matrix = this->matrix_dual_;
        rhs = this->rhs_dual_;
        sol = this->solD_;
        linear_solver = this->linear_solver_dual_;
        space = this->space_dual_;
        use_pressure_filter_dual_ = base_params_["DualLinearSolver"]["OuterPressureFilter"].get<int>( 0 );
        use_pressure_filter = use_pressure_filter_dual_;
        this->solP_next_->Update ( );
        this->solD_->Update ( );
        this->solD_prev_->Update ( );
        this->solD_next_->Update ( );
    }

    // assemble A
    this->compute_matrix ( NULL, matrix, mode );

    // assemble rhs =  A*sol_old + b(sol_prev)
    //              =: A*sol_old + b
    rhs->Zeros ( );
    this->compute_residual ( NULL, rhs, mode );
    rhs->Update ( );

    double norm_old_rhs = rhs->Norm2 ( );

    // prepare solution process
    LAD::VectorType cor;
    cor.CloneFromWithoutContent ( *rhs );
    cor.Zeros ( );

    // solve  A * cor = rhs + res
    linear_solver->SetupOperator ( *matrix );
    linear_solver->Solve ( *rhs, &cor );
    cor.Update ( );
    interpolate_constrained_vector ( *space, cor );
    cor.Update ( );

    double norm_cor = cor.Norm2 ( );

    // update solution sol = sol_old - cor
    // -> rhs_new     = A * sol + b
    //                = A * (sol_old - cor) + b
    //                = A * (sol_old - A^{-1} (A*sol_old + b) - A^{-1} * res) + b
    //                = A * (-A^{-1} * b - A^{-1} * res) + b
    //                = - res
    sol->Axpy ( cor, -1. );
    sol->Update ( );

    if ( use_pressure_filter )
    {
        this->ApplyFilter ( *sol );
        sol->Update ( );
    }

    // testing each nth calculation
    double norm_new_rhs = 0.;

    if ( this->time_step_ % this->test_step_ == 0 )
    {
        // evaluate residual
        rhs->Zeros ( );

        // rhs_new = A * sol + b
        compute_residual ( NULL, rhs, mode );
        rhs->Update ( );
        norm_new_rhs = rhs->Norm2 ( );

        // handle diverging situation
        if ( norm_new_rhs >= this->lp_div_tol_ )
        {
            //    visualize_solution_dual(rhs_dual_, "out/residual_dual", 0);
            state = false; // diverging solution
        }
        if ( rank_ == master_rank_ )
        {
            std::cout << "  COR_NORM  " << norm_cor << std::endl;
            std::cout << "  RHS_NORM  " << norm_new_rhs << std::endl;
            std::cout << "  INIT_RHS_NORM  " << norm_old_rhs << std::endl;
            std::cout << "  These should be equal:  Final res comp: " << norm_new_rhs << ", solver abs res: " << linear_solver->res ( ) << std::endl;
        }
    }

    timer.stop ( );

    // interminable_assert(state == 0);
    if ( rank_ == master_rank_ )
    {
        std::cout << "  Time for solution of linear problem: " << timer.get_duration ( ) << " s with " << linear_solver->iter ( )
                << "  iterations " << std::endl
                << "  final absolute residual " << linear_solver->res ( ) << std::endl
                << "  final relative residual " << linear_solver->res ( ) / norm_old_rhs << std::endl << std::endl;
    }
    if ( state == false )
    {
        MPI_Barrier ( MPI_COMM_WORLD );
        if ( rank_ == master_rank_ )
            std::cout << "state = FALSE -> I better should quit the simulation" << std::endl;
    }
    return state;
}

/// ***************************************************************************
/// IO
/// ***************************************************************************

void MetFlowApp::compute_comm_pattern ( std::string const& filename )
{
    std::vector<int> connected_to_proc ( num_partitions_, 0 );
    for ( int l = 0; l < mesh_->num_entities ( DIM ); ++l )
    {
        int subdomain;
        mesh_->get_attribute_value ( "_sub_domain_", DIM, l, &subdomain );
        connected_to_proc[subdomain] = 1;
    }

    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "> Compute communication pattern " << std::endl;
        std::vector< std::vector<int> > comm_pattern ( num_partitions_ );
        comm_pattern[0] = connected_to_proc;

        for ( int p = 1; p < num_partitions_; ++p )
        {
            comm_pattern[p].resize ( num_partitions_, 0 );
            MPI_Status state;
            MPI_Recv ( &comm_pattern[p][0], num_partitions_, MPI_INT, p, p, this->comm_, &state );
        }

        ofstream out;
        std::string path = this->root_ + filename;
        out.open ( path.c_str ( ), ios::out );
        for ( int p = 0; p < num_partitions_; ++p )
        {
            for ( int q = 0; q < num_partitions_; ++q )
            {
                out << comm_pattern[p][q] << " ";
            }
            out << "\n";
        }
        out.close ( );
    }
    else
    {
        MPI_Send ( &connected_to_proc[0], num_partitions_, MPI_INT, 0, rank ( ), this->comm_ );
    }
}

void MetFlowApp::visualize_function ( const LAD::VectorType& sol, const VectorSpace<double>* space, std::string const& prefix, int time_step, std::vector<std::string>& var_names )
{
    // Filename
    std::stringstream input;
    if ( time_step < 10 )
        input << prefix << ".000" << time_step;
    else if ( time_step < 100 )
        input << prefix << ".00" << time_step;
    else if ( time_step < 1000 )
        input << prefix << ".0" << time_step;
    else
        input << prefix << "." << time_step;

    std::string tmp_filename = input.str ( );

    if ( this->num_partitions_ > 1 )
        input << ".pvtu";
    else
        input << ".vtu";

    std::string visu_filename = input.str ( );

    // CellVisualization
    ParallelCellVisualization<double> visu ( *space, 1, comm_, master_rank ( ) );

    for ( int i = 0; i < space_->get_nb_var ( ); i++ )
    {
        visu.visualize ( EvalFeFunction<LAD>( *space, sol, i ), var_names[i] );
    }

    visu.write ( visu_filename );
}

void MetFlowApp::visualize_error_indicators ( int time_step )
{
    // Filename
    std::stringstream input;
    std::string prefix = "ind";

    std::stringstream pre;
    pre << this->root_ << "/indicators/" << prefix << "." << this->adapt_counter_;

    if ( time_step < 10 )
        input << ".000" << time_step;
    else if ( time_step < 100 )
        input << ".00" << time_step;
    else if ( time_step < 1000 )
        input << ".0" << time_step;
    else
        input << "." << time_step;

    if ( this->num_partitions_ > 1 )
        input << ".pvtu";
    else
        input << ".vtu";

    std::string visu_filename = pre.str ( ) + input.str ( );

    // CellVisualization
    ParallelCellVisualization<double> visu ( *space_, 1, comm_, master_rank ( ) );

    std::vector< std::string> est_type ( 4 );
    est_type[0] = "res_h";
    est_type[1] = "weight_h";
    est_type[2] = "res_tau";
    est_type[3] = "weight_tau";

    for ( int q = 0; q<this->num_eq_; ++q )
    {
        // cell indicators
        for ( int e = 0; e < 4; ++e )
        {
            std::vector<double> attr_val ( mesh_->num_entities ( DIM ), 0. );
            std::vector<double> attr_val_dual ( mesh_->num_entities ( DIM ), 0. );
            std::string attr_name = "cell_primal_" + static_cast < ostringstream* > ( &( ostringstream ( ) << q ) )->str ( ) + "_" + est_type[e];
            std::string attr_name_dual = "cell_dual_" + static_cast < ostringstream* > ( &( ostringstream ( ) << q ) )->str ( ) + "_" + est_type[e];
            AttributePtr attr_ptr = mesh_->get_attribute ( attr_name, DIM );
            AttributePtr attr_ptr_dual = mesh_->get_attribute ( attr_name_dual, DIM );

            for ( mesh::EntityIterator it = mesh_->begin ( DIM ); it != mesh_->end ( DIM ); ++it )
            {
                attr_val .at ( it->index ( ) ) = attr_ptr ->get_double_value ( it->index ( ) );
                attr_val_dual.at ( it->index ( ) ) = attr_ptr_dual->get_double_value ( it->index ( ) );
            }
            visu.visualize_cell_data ( attr_val, attr_name );
            visu.visualize_cell_data ( attr_val_dual, attr_name_dual );
        }

        // facet indicators
        for ( int e = 0; e < 4; ++e )
        {
            std::vector<double> attr_val ( mesh_->num_entities ( DIM ), 0. );
            std::vector<double> attr_val_dual ( mesh_->num_entities ( DIM ), 0. );
            std::string attr_name = "jump_primal_" + static_cast < ostringstream* > ( &( ostringstream ( ) << q ) )->str ( ) + "_" + est_type[e];
            std::string attr_name_dual = "jump_dual_" + static_cast < ostringstream* > ( &( ostringstream ( ) << q ) )->str ( ) + "_" + est_type[e];
            AttributePtr attr_ptr = mesh_->get_attribute ( attr_name, DIM );
            AttributePtr attr_ptr_dual = mesh_->get_attribute ( attr_name_dual, DIM );

            for ( mesh::EntityIterator it = mesh_->begin ( DIM ); it != mesh_->end ( DIM ); ++it )
            {
                attr_val .at ( it->index ( ) ) = attr_ptr ->get_double_value ( it->index ( ) );
                attr_val_dual.at ( it->index ( ) ) = attr_ptr_dual->get_double_value ( it->index ( ) );
            }
            visu.visualize_cell_data ( attr_val, attr_name );
            visu.visualize_cell_data ( attr_val_dual, attr_name_dual );
        }
    }

    visu.write ( visu_filename );
}

void MetFlowApp::visualize_reduced_error_indicators ( int time_step, int mesh_index )
{
    VectorSpace<DATATYPE>* active_space = this->dmh_->get_space_by_index ( mesh_index, 1 );

    // Filename
    std::stringstream input;
    std::string prefix = "red_ind";

    std::stringstream pre;
    pre << this->root_ << "/indicators/" << prefix << "." << this->adapt_counter_;

    if ( time_step < 10 )
        input << ".000" << time_step;
    else if ( time_step < 100 )
        input << ".00" << time_step;
    else if ( time_step < 1000 )
        input << ".0" << time_step;
    else
        input << "." << time_step;

    if ( this->num_partitions_ > 1 )
        input << ".pvtu";
    else
        input << ".vtu";

    std::string visu_filename = pre.str ( ) + input.str ( );

    // CellVisualization
    ParallelCellVisualization<double> visu ( *active_space, 1, comm_, master_rank ( ) );

    visu.visualize_cell_data ( this->space_indicator_[time_step], "space_ind" );
    visu.visualize_cell_data ( this->time_indicator_[time_step], "time_ind" );

    visu.write ( visu_filename );
}

void MetFlowApp::visualize_reduced_space_indicators ( )
{
    int num_mesh = this->dmh_->num_mesh ( );

    // Filename
    std::string prefix = "red_space_ind";

    for ( int m = 0; m < num_mesh; ++m )
    {
        std::stringstream input;
        std::stringstream pre;
        pre << this->root_ << "/indicators/" << prefix << "." << this->adapt_counter_ << "." << m;

        if ( this->num_partitions_ > 1 )
            input << ".pvtu";
        else
            input << ".vtu";

        std::string visu_filename = pre.str ( ) + input.str ( );

        // CellVisualization
        VectorSpace<DATATYPE>* active_space = this->dmh_->get_space_by_index ( m, 2 );
        ParallelCellVisualization<double> visu ( *active_space, 1, comm_, master_rank ( ) );

        visu.visualize_cell_data ( this->reduced_space_indicator_[m], "red_space_ind" );
        visu.write ( visu_filename );
    }
}

void MetFlowApp::visualize_reduced_time_indicators ( )
{
    // Filename
    std::stringstream input;
    std::string prefix = "red_time_ind";

    std::stringstream pre;
    pre << this->root_ << "/indicators/" << prefix << "." << this->adapt_counter_ << ".csv";
    std::string visu_filename = pre.str ( );

    std::stringstream log_pre;
    log_pre << this->root_ << "/log/time_estimator.csv";
    std::string log_filename = log_pre.str ( );

    ofstream myfile;
    myfile.open ( visu_filename.c_str ( ) );

    ofstream logfile;
    logfile.open ( log_filename.c_str ( ), std::ios_base::app );

    for ( int t = 0; t<this->reduced_time_indicator_.size ( ); ++t )
    {
        myfile << this->get_time ( t + 1 ) << ", ";
        logfile << this->get_time ( t + 1 ) << ", ";
    }
    myfile << "\n";
    logfile << "\n";
    for ( int t = 0; t<this->reduced_time_indicator_.size ( ); ++t )
    {
        myfile << this->reduced_time_indicator_[t] << ", ";
        logfile << this->reduced_time_indicator_[t] << ", ";
    }
    myfile << "\n";
    myfile.close ( );
    logfile << "\n";
    logfile.close ( );
}

void MetFlowApp::visualize_patch_interpolation ( int time_step )
{
    Timer timer;
    timer.start ( );

    std::stringstream pre;
    pre << this->root_ << "/out_patch/primal_patch." << this->adapt_counter_;

    std::vector<std::string> var_names;
    for ( int i = 0; i<this->num_vars_; ++i )
    {
        var_names.push_back ( "patch_" + static_cast < ostringstream* > ( &( ostringstream ( ) << i ) )->str ( ) );
    }
    this->visualize_function ( fineP_prev_, fine_space_, pre.str ( ), time_step, var_names );
    this->visualize_function ( fineP_, fine_space_, pre.str ( ), time_step + 1, var_names );

    std::stringstream pre_dual;
    pre_dual << this->root_ << "/out_patch/dual_patch." << this->adapt_counter_;

    std::vector<std::string> var_names_dual;
    for ( int i = 0; i<this->num_vars_; ++i )
    {
        var_names_dual.push_back ( "dual_patch_" + static_cast < ostringstream* > ( &( ostringstream ( ) << i ) )->str ( ) );
    }
    this->visualize_function ( fineD_prev_, fine_space_dual_, pre_dual.str ( ), time_step, var_names_dual );
    this->visualize_function ( fineD_, fine_space_dual_, pre_dual.str ( ), time_step + 1, var_names_dual );

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  took " << timer.get_duration ( ) << "sec" << std::endl;

}

void MetFlowApp::read_file ( LAD::VectorType& sol,
                             std::string const& filename,
                             std::string const& groupname,
                             std::string const& datasetname )
{
    Timer timer;
    timer.start ( );
    if ( use_hdf5_ == true )
    {
        std::string new_filename = filename + "." + static_cast < ostringstream* > ( &( ostringstream ( ) << this->adapt_counter_ ) )->str ( ) + ".h5";

        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "> Read HDF5 file " << new_filename << ", group " << groupname << ", data set " << datasetname << std::endl;
        }
        MPI_Barrier ( this->comm_ );

#ifdef WITH_HDF5
        sol.ReadHDF5 ( new_filename, groupname, datasetname );
        /*
                double* buffer;

                H5FilePtr      file_ptr    ( new H5File   ( filename, "r", this->comm_ ) );
                H5GroupPtr      group_ptr    ( new H5Group  ( file_ptr, "time", "r" ) );
                H5DatasetPtr dataset_ptr( new H5Dataset( group_ptr, num_partitions_, datasetname, "r", buffer ) );
                dataset_ptr->read( 1, 0, buffer );
                time = *buffer;
         */

#else
        if ( rank ( ) == master_rank ( ) )
            std::cout << "HDF5 not set up!" << std::endl;
        quit_program ( );
#endif
    }
}

void MetFlowApp::write_file ( LAD::VectorType & sol,
                              std::string const& filename,
                              std::string const& groupname,
                              std::string const& datasetname )
{
    Timer timer;
    timer.start ( );
    if ( use_hdf5_ == true )
    {
        std::string new_filename = filename + "." + static_cast < ostringstream* > ( &( ostringstream ( ) << this->adapt_counter_ ) )->str ( ) + ".h5";
        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "> Write HDF5 file " << new_filename << ", group " << groupname << ", data set " << datasetname << std::endl;
        }
        MPI_Barrier ( this->comm_ );
#ifdef WITH_HDF5
        sol.WriteHDF5 ( new_filename, groupname, datasetname );
        /*
                double* data;
                data = &time;

                H5FilePtr      file_ptr    ( new H5File   ( filename, "w", this->comm_ ) );
                H5GroupPtr      group_ptr    ( new H5Group  ( file_ptr, "time", "w" ) );
                H5DatasetPtr dataset_ptr( new H5Dataset( group_ptr, num_partitions_, datasetname, "w", data ) );
                dataset_ptr->write( 1, rank_, data );
         */
#else
        std::cout << "HDF5 not set up!" << std::endl;
        quit_program ( );
#endif
    }
    else
    {
        if ( rank ( ) > 1 )
        {
            std::cout << "Implement me: This function is not parallel-ready." << std::endl;
            quit_program ( );
        }

        std::cout << "> Write ASCII file " << filename << std::endl;

        // filename    -> used as filename
        // groupname   -> ignored
        // datasetname -> used as file-extension
        // => ASCII filename: [filename].[datasetname].dat.ascii

        //assert(groupname=="NoHDF5");
        //assert(datasetname=="NoHDF5");

        // convert CoupledVector to std::vector

        std::vector<LAD::DataType> backup_vec ( sol.size_global ( ), 1.e20 );

        std::vector<int> dof_ids;
        std::vector<double> values;

        //PpVector<LAD> pp_sol;
        //pp_sol.Init(MPI_COMM_WORLD, space_.dof(), sol);
        //pp_sol.InitStructure();
        //pp_sol.UpdateValues();
        //pp_sol.GetDofsAndValues(dof_ids, values);

        for ( int i = 0; i < values.size ( ); ++i )
            backup_vec.at ( dof_ids[i] ) = values.at ( i );

        // store values to file

        std::stringstream ss;
        ss << filename << "." << datasetname << ".dat.ascii";
        std::ofstream out_file;
        out_file.open ( ss.str ( ).c_str ( ) );

        // write number of entries in the very first line
        out_file << values.size ( ) << std::endl;

        // write the data, one value per line
        for ( int i = 0; i < backup_vec.size ( ); ++i )
            out_file << std::setfill ( ' ' ) << std::setprecision ( 20 ) << std::setw ( 25 ) << backup_vec.at ( i ) << "\n";

        out_file.close ( );
    }
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) )
        std::cout << "  took " << timer.get_duration ( ) << " sec " << std::endl;

}

void MetFlowApp::read_file ( LAD::VectorType& sol,
                             int time_step,
                             std::string const& filename,
                             std::string const& groupname,
                             std::string const& datasetname,
                             int mode )
{
    // last solution belongs to a different fe space as current one
    int read_mesh_index = this->dmh_->mesh_index ( time_step );
    if ( read_mesh_index != this->active_mesh_index_ )
    {
        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "  Read solution with mesh index " << read_mesh_index
                    << " and interpolate to mesh index " << this->active_mesh_index_ << std::endl;
        }

        // create temporary vector which fits to read_space
LAD:
        VectorType tmp;
#ifdef USE_HYPRE
        tmp.Init ( comm_, *this->dmh_->get_couplings_by_index ( read_mesh_index, mode ) );
#else
        tmp.Init ( comm_, *this->dmh_->get_couplings_by_index ( read_mesh_index, mode ), la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
        tmp.InitStructure ( );
#endif
        // read in solution
        this->read_file ( tmp, filename, groupname, datasetname );

        // interpolate solution to new space
        this->interpolate_vector ( tmp, read_mesh_index, sol, this->active_mesh_index_, mode );
    }
    else
    {
        this->read_file ( sol, filename, groupname, datasetname );
    }
}

/// Create empty log files for flow characteristics

void MetFlowApp::create_log_file ( int offset, std::string filename, std::string suffix, std::string seperator, bool print_counter )
{
    if ( offset == 0 )
    {
        ofstream out;
        if ( print_counter )
        {
            std::string path = this->root_ + filename + "." + static_cast < ostringstream* > ( &( ostringstream ( ) << this->adapt_counter_ ) )->str ( ) + "." + suffix;
            out.open ( path.c_str ( ), ios::out );
        }
        else
        {
            std::string path = this->root_ + filename + "." + suffix;
            out.open ( path.c_str ( ), ios::out );
        }
        out.close ( );
    }
    else
    {
        ifstream in;
        if ( print_counter )
        {
            std::string path = this->root_ + filename + "." + static_cast < ostringstream* > ( &( ostringstream ( ) << this->adapt_counter_ ) )->str ( ) + "." + suffix;
            in.open ( path.c_str ( ) );
        }
        else
        {
            std::string path = this->root_ + filename + "." + suffix;
            in.open ( path.c_str ( ) );
        }
        std::vector < std::vector<double> > data;

        const std::string& delim = seperator;

        string line;
        string strnum;

        while ( getline ( in, line ) )
        {
            data.push_back ( vector<double>( ) );

            for ( string::const_iterator i = line.begin ( ); i != line.end ( ); ++i )
            {
                // If i is not a delim, then append it to strnum
                if ( delim.find ( *i ) == string::npos )
                {
                    strnum += *i;
                    if ( i + 1 != line.end ( ) ) // If it's the last char, do not continue
                        continue;
                }

                // if strnum is still empty, it means the previous char is also a
                // delim (several delims appear together). Ignore this char.
                if ( strnum.empty ( ) )
                    continue;

                // If we reach here, we got a number. Convert it to double.
                double number;

                istringstream ( strnum ) >> number;
                data.back ( ).push_back ( number );

                strnum.clear ( );
            }
        }
        in.close ( );

        ofstream out;
        std::string path = this->root_ + filename;
        out.open ( path.c_str ( ) );
        out.precision ( 10 );
        out << std::scientific;

        for ( int i = 0; i < data.size ( ); ++i )
        {
            for ( int j = 0; j < data[i].size ( ); ++j )
            {
                out << data[i][j] << " ";

            }
            out << "\n";
        }
        out.close ( );
        data.clear ( );
    }
}

/// ***************************************************************************
/// POSTPROCESSING
/// ***************************************************************************
/// ***************************************************************************
/// TOOLS
/// ***************************************************************************

double MetFlowApp::compute_pressure_int ( LAD::VectorType& u, int var )
{
    double recv;
    double total_pressure;
    u.Update ( );

    PressureIntegral int_p ( u, cyl_coord_, var );
    global_asm_.integrate_scalar ( *space_, int_p, total_pressure );

    MPI_Allreduce ( &total_pressure, &recv, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    return recv;
}

double MetFlowApp::compute_volume_int ( )
{
    double integrated_vol;
    double recv;
    VolumeIntegral vol_int ( cyl_coord_ );
    global_asm_.integrate_scalar ( *space_, vol_int, integrated_vol );

    MPI_Allreduce ( &integrated_vol, &recv, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    return recv;
}

void MetFlowApp::ApplyFilter ( LAD::VectorType& u )
{
    if ( !this->use_pressure_filter ( ) )
    {
        return;
    }

    this->ApplyFilter ( u, DIM );
#ifdef AUGMENT_PRESS
    this->ApplyFilter ( u, this->aug_p_var_ );
#endif
}

void MetFlowApp::ApplyFilter ( LAD::VectorType& u, int var )
{
    double total_pressure = compute_pressure_int ( u, var );
    double integrated_vol = compute_volume_int ( );
    const double average_pressure = total_pressure / integrated_vol;

    LOG_INFO ( "pressure_filter", "Average pressure before filter = " << average_pressure );

    LAD::VectorType pressure_correction;
    pressure_correction.CloneFromWithoutContent ( u );
    pressure_correction.Zeros ( );

    // set value for pressure dofs to average pressure
    std::vector<int> cell_p_dofs;
    std::vector<int> local_p_dofs;
    for ( EntityIterator it = mesh_->begin ( DIM ), end = mesh_->end ( DIM ); it != end; ++it )
    {
        cell_p_dofs.clear ( );
        space_->GetDofIndices ( var, *it, &cell_p_dofs );
        for ( int i = 0, sz = cell_p_dofs.size ( ); i < sz; ++i )
        {
            if ( space_->dof ( ).is_dof_on_sd ( cell_p_dofs[i] ) )
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

    pressure_correction.SetValues ( vec2ptr ( local_p_dofs ), local_p_dofs.size ( ), vec2ptr ( p_correction_values ) );

    u.Axpy ( pressure_correction, -1. );

    total_pressure = compute_pressure_int ( u, var );
    LOG_INFO ( "pressure_filter", "Average pressure after filter = " << total_pressure / integrated_vol );
}

int MetFlowApp::rank ( ) const
{
    return rank_;
}

int MetFlowApp::master_rank ( ) const
{
    return master_rank_;
}

int MetFlowApp::num_partitions ( ) const
{
    return num_partitions_;
}

bool MetFlowApp::use_pressure_filter ( )
{
    if ( this->problem_mode_ == 1 )
    {
        return use_pressure_filter_;
    }
    if ( this->problem_mode_ == -1 )
    {
        return use_pressure_filter_dual_;
    }
}

void MetFlowApp::evaluate_qoi_int ( std::string filename, LAD::VectorType* sol, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Evaluate integral part of goal functional " << std::endl;
    Timer timer;
    timer.start ( );

    // Set soultion vector w.r.t. which quantities should be compute
    local_asm_primal_->set_goal_functional ( this->j_ );
    local_asm_primal_->set_solP ( *sol );

    // compute norm
    std::vector<double> qoi;
    qoi.clear ( );
    local_asm_primal_->set_scalar_mode_to_goal_int ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( *local_asm_primal_ ), qoi );

    double total_qoi = std::accumulate ( qoi.begin ( ), qoi.end ( ), 0. );
    double global_qoi = 0.;
    MPI_Reduce ( &total_qoi, &global_qoi, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );

    // file output
    if ( rank ( ) == master_rank ( ) )
    {
        std::string path = this->root_ + filename + "." + static_cast < ostringstream* > ( &( ostringstream ( ) << this->adapt_counter_ ) )->str ( ) + ".txt";
        ofstream out;
        out.open ( path.c_str ( ), ios::out | ios::app );

        out.precision ( 10 );
        out << std::scientific;
        out << time << " " << global_qoi << "\n";
        out.close ( );
    }

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) )
        std::cout << "  took " << timer.get_duration ( ) << " sec" << std::endl;

    // Reset solution vector
    local_asm_primal_->set_solP ( *this->solP_ );
}

void MetFlowApp::evaluate_qoi_fin ( std::string filename, LAD::VectorType* sol, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Evaluate final step part of goal functional " << std::endl;
    Timer timer;
    timer.start ( );

    // Set soultion vector w.r.t. which quantities should be computed
    this->local_asm_primal_->set_goal_functional ( this->j_ );
    this->local_asm_primal_->set_solP ( *sol );

    // compute norm
    std::vector<double> qoi;
    qoi.clear ( );
    this->local_asm_primal_->set_scalar_mode_to_goal_fin ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( *this->local_asm_primal_ ), qoi );

    double total_qoi = std::accumulate ( qoi.begin ( ), qoi.end ( ), 0. );
    double global_qoi = 0.;
    MPI_Reduce ( &total_qoi, &global_qoi, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );

    // file output
    if ( rank ( ) == master_rank ( ) )
    {
        std::string path = this->root_ + filename + "." + static_cast < ostringstream* > ( &( ostringstream ( ) << this->adapt_counter_ ) )->str ( ) + ".txt";
        ofstream out;
        out.open ( path.c_str ( ), ios::out | ios::app );

        out.precision ( 10 );
        out << std::scientific;
        out << time << " " << global_qoi << "\n";
        out.close ( );
    }

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) )
        std::cout << "  took " << timer.get_duration ( ) << " sec" << std::endl;

    // Reset solution vector
    this->local_asm_primal_->set_solP ( *this->solP_ );
}

void MetFlowApp::evaluate_inner_prod ( std::string filename, LAD::VectorType* solA, LAD::VectorType* solB, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Evaluate inner product in W_{1,2}(Omega) " << std::endl;
    Timer timer;
    timer.start ( );

    // Set soultion vector w.r.t. which quantities should be computed
    this->local_asm_prod_->set_left_vector ( *solA );
    this->local_asm_prod_->set_right_vector ( *solB );
    std::vector<bool> L2, H1, H2;

    L2 = this->local_asm_primal_->get_L2_perturb_flags ( );
    H1 = this->local_asm_primal_->get_H1_perturb_flags ( );
    H2.resize ( this->num_vars_, false );

    this->local_asm_prod_->set_mode ( L2, H1, H2 );
    this->local_asm_prod_->mark_vector_field ( this->vel_var_ );

    // compute norm
    std::vector<double> ip;
    ip.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( *this->local_asm_prod_ ), ip );

    double total_ip = std::accumulate ( ip.begin ( ), ip.end ( ), 0. );
    double global_ip = 0.;
    MPI_Reduce ( &total_ip, &global_ip, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );

    // file output
    if ( rank ( ) == master_rank ( ) )
    {
        std::string path = this->root_ + filename + "." + static_cast < ostringstream* > ( &( ostringstream ( ) << this->adapt_counter_ ) )->str ( ) + ".txt";
        ofstream out;
        out.open ( path.c_str ( ), ios::out | ios::app );

        out.precision ( 10 );
        out << std::scientific;
        out << time << " " << global_ip << "\n";
        out.close ( );
    }

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) )
        std::cout << "  took " << timer.get_duration ( ) << " sec" << std::endl;
}

void MetFlowApp::perturb_solution ( int time_step )
{
    int perturb_primal = this->base_params_["Perturbation"]["PerturbEF"].get<int>( 0 );
    double perturb_scale = this->base_params_["Perturbation"]["ExternalForce"]["Scale"].get<double>( 0. );
    std::string filename_perturb = this->base_params_["Perturbation"]["ExternalForce"]["SnapshotsIn"].get<std::string>( );
    std::string groupname_perturb = this->base_params_["Perturbation"]["ExternalForce"]["SnapshotsGroup"].get<std::string>( );
    std::string prefix_perturb = this->base_params_["Perturbation"]["ExternalForce"]["SnapshotsPrefix"].get<std::string>( );

    // Set perturbation
    if ( perturb_primal > 0 )
    {
        std::stringstream ss_perturb;
        ss_perturb << prefix_perturb << time_step;
        MetFlowApp::read_file ( *this->perturb_, filename_perturb, groupname_perturb, ss_perturb.str ( ) );

        if ( time_step >= 2 )
        {
            std::stringstream ss_perturb_prev;
            ss_perturb_prev << prefix_perturb << time_step - 1;
            MetFlowApp::read_file ( *this->perturb_prev_, filename_perturb, groupname_perturb, ss_perturb_prev.str ( ) );
        }
        else
        {
            this->perturb_prev_->Zeros ( );
        }

        this->perturb_->Scale ( perturb_scale );
        this->perturb_prev_->Scale ( perturb_scale );

        this->perturb_->Update ( );
        this->perturb_prev_->Update ( );

        this->local_asm_primal_->set_perturb ( *this->perturb_ );
        this->local_asm_primal_->set_perturb_prev ( *this->perturb_prev_ );
        this->local_asm_primal_->set_perturb_scale ( 1. );
    }
}

void MetFlowApp::compute_cell_diameters ( std::vector<double>& hK )
{
    hK.clear ( );
    global_asm_.assemble_scalar ( *space_, boost::ref ( *local_asm_est_ ), hK );
}

/// ***************************************************************************
/// ADAPTIVE REFINEMENT
/// ***************************************************************************

double MetFlowApp::get_delta_t ( int time_step )
{
    if ( this->is_stationary_ )
    {
        return this->duration_;
    }

    return this->t_mesh_.delta_t ( time_step, this->adapt_counter_ );
}

int MetFlowApp::get_num_intervals ( )
{
    return this->t_mesh_.num_intervals ( this->adapt_counter_ );
}

double MetFlowApp::get_time ( int time_step )
{
    return this->t_mesh_.time ( time_step, this->adapt_counter_ );
}

void MetFlowApp::write_mesh_change_list ( int adapt_counter, std::vector<DATATYPE>& mesh_change_times )
{
    ofstream out;
    std::string path = this->root_ + "/mesh/mesh_change." + static_cast < ostringstream* > ( &( ostringstream ( ) << adapt_counter ) )->str ( ) + ".txt";
    out.open ( path.c_str ( ) );

    for ( int l = 0; l < mesh_change_times.size ( ); ++l )
    {
        out << mesh_change_times[l] << "\n";
    }
    out.close ( );
}

void MetFlowApp::read_mesh_change_list ( int adapt_counter, std::vector<DATATYPE>& mesh_change_times )
{
    // read in points in time when to change the mesh
    ifstream in;
    string path;
    if ( adapt_counter < 0 )
    {
        path = this->root_ + "/in/mesh_change.txt";
    }
    else
    {
        path = this->root_ + "/mesh/mesh_change." + static_cast < ostringstream* > ( &( ostringstream ( ) << adapt_counter ) )->str ( ) + ".txt";
    }
    in.open ( path.c_str ( ) );
    std::vector < std::vector<double> > data;

    const std::string& delim = " ";

    string line;
    string strnum;

    while ( getline ( in, line ) )
    {
        data.push_back ( vector<double>( ) );

        for ( string::const_iterator i = line.begin ( ); i != line.end ( ); ++i )
        {
            // If i is not a delim, then append it to strnum
            if ( delim.find ( *i ) == string::npos )
            {
                strnum += *i;
                if ( i + 1 != line.end ( ) ) // If it's the last char, do not continue
                    continue;
            }

            // if strnum is still empty, it means the previous char is also a
            // delim (several delims appear together). Ignore this char.
            if ( strnum.empty ( ) )
                continue;

            // If we reach here, we got a number. Convert it to double.
            double number;

            istringstream ( strnum ) >> number;
            data.back ( ).push_back ( number );

            strnum.clear ( );
        }
    }
    in.close ( );

    int num_changes = data.size ( );
    mesh_change_times.clear ( );
    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "> Mesh change times: " << std::endl;
    }
    for ( int t = 0; t < num_changes; ++t )
    {
        if ( data[t][0] >= 0. )
        {
            mesh_change_times.push_back ( data[t][0] );
            if ( rank ( ) == master_rank ( ) )
            {
                std::cout << mesh_change_times[t] << " ";
            }
        }
    }
    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << std::endl;
    }
}

void MetFlowApp::estimate_error ( int adaption_counter, int num_time_steps )
{
    Timer timer;

    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "***************************************************" << std::endl;
        std::cout << "NOW THE ERROR ESTIMATION (AC=" << adaption_counter << ", TS=" << num_time_steps << ")" << std::endl;
        std::cout << "***************************************************" << std::endl << std::endl;
        std::cout << "> Min primal step: " << min_step_primal_ << std::endl;
        std::cout << "> Max primal step: " << max_step_primal_ << std::endl;
        std::cout << "> Min dual step  : " << min_step_dual_ << std::endl;
        std::cout << "> Max dual step  : " << max_step_dual_ << std::endl;
    }

    timer.reset ( );
    timer.start ( );

    // For HDF5
    std::string filename_primal = base_params_["BackUp"]["PrimalSnapshotsIn"].get<std::string>( );
    std::string groupname_primal = base_params_["BackUp"]["PrimalSnapshotsGroupIn"].get<std::string>( );
    std::string prefix_primal = base_params_["BackUp"]["PrimalSnapshotsPrefixIn"].get<std::string>( );
    int prefix_primal_with_timestep = base_params_["BackUp"]["PrimalSnapshotsTimeStepIn"].get<int>( 1 );
    filename_primal = this->root_ + "/" + filename_primal;

    std::string filename_dual = base_params_["BackUp"]["DualSnapshotsIn"].get<std::string>( );
    std::string groupname_dual = base_params_["BackUp"]["DualSnapshotsGroupIn"].get<std::string>( );
    std::string prefix_dual = base_params_["BackUp"]["DualSnapshotsPrefixIn"].get<std::string>( );
    int prefix_dual_with_timestep = base_params_["BackUp"]["DualSnapshotsTimeStepIn"].get<int>( 1 );
    filename_dual = this->root_ + "/" + filename_dual;

    //this->prepare_higher_order_space();
    this->prepare_assembler_est ( );

    // -------------------------------------------
    // prepare error estimators
    std::vector< std::vector< std::vector<double> > > est_cell_residual ( num_time_steps );
    std::vector< std::vector< std::vector<double> > > est_cell_timejump ( num_time_steps );
    std::vector< std::vector< std::vector<double> > > est_interface ( num_time_steps );

    // [time_step][entity_index] { {primal_eq_1_res_h, primal_eq_1_weight_h, primal_eq_1_res_tau, primal_eq_1_weight_tau;
    //                              primal_eq_n_res_h, primal_eq_n_weight_h, primal_eq_n_res_tau, primal_eq_n_weight_tau}
    //                                  ........
    //                             {dual_eq_1_res_h, dual_eq_1_weight_h, dual_eq_1_res_tau, dual_eq_1_weight_tau;
    //                              dual_eq_n_res_h, dual_eq_n_weight_h, dual_eq_n_res_tau, dual_eq_n_weight_tau}}

    std::vector< std::string > est_type ( 4 );
    est_type[0] = "r_h";
    est_type[1] = "w_h";
    est_type[2] = "r_tau";
    est_type[3] = "w_tau";

    int time_step = 0;

    // setup correct mesh and fe space
    this->dmh_->update ( time_step, 0, true );
    this->prepare_higher_order_space ( );

    // first solutions
    double time;
    double cur_time = 0.;
    std::stringstream ss_primal_next_0;
    ss_primal_next_0 << prefix_primal << 0;
    if ( this->min_step_primal_ == 0 )
    {
        this->read_file ( *solP_next_, 0, filename_primal, groupname_primal, ss_primal_next_0.str ( ), 1 );
    }
    else
    {
        solP_next_->Zeros ( );
    }
    std::stringstream ss_dual_next_0;
    ss_dual_next_0 << prefix_dual << 0;
    if ( this->min_step_dual_ == 0 )
    {
        this->read_file ( *solD_next_, 0, filename_dual, groupname_dual, ss_dual_next_0.str ( ), -1 );
    }
    else
    {
        solD_next_->Zeros ( );
    }

    this->indicator_mesh_indices_.clear ( );
    while ( time_step < this->get_num_intervals ( ) )
    {
        // get time step size
        double dT_pc = this->get_delta_t ( time_step + 1 );
        double dT_cn = this->get_delta_t ( time_step + 2 );

        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << std::endl << "======================================================" << std::endl;
            std::cout << "> Time interval     " << time_step + 1 << " and " << time_step + 2 << " of " << this->get_num_intervals ( ) << std::endl;
            std::cout << "> Time interval      (" << cur_time << ", " << cur_time + dT_pc << ") and ("
                    << cur_time + dT_pc << ", " << cur_time + dT_pc + dT_cn << ")" << std::endl;
            std::cout << "> Mesh index        " << this->dmh_->mesh_index ( time_step + 1 ) << " and " << this->dmh_->mesh_index ( time_step + 2 ) << std::endl;
            std::cout << "> Active mesh index " << this->active_mesh_index_ << std::endl;
        }
        // select correct mesh and fe space
        bool updated = this->dmh_->update ( time_step + 1, 0, false );
        if ( updated )
        {
            this->prepare_higher_order_space ( );
        }

        // ---------------------------------------
        // preparation
        solP_prev_->CopyFrom ( *solP_next_ );
        solD_prev_->CopyFrom ( *solD_next_ );

        // read current solutions
        if ( time_step + 1 <= this->max_step_primal_ )
        {
            std::stringstream ss_primal;
            ss_primal << prefix_primal << time_step + 1;
            this->read_file ( *solP_, time_step + 1, filename_primal, groupname_primal, ss_primal.str ( ), 1 );
        }
        else
        {
            solP_->Zeros ( );
        }

        if ( time_step + 1 <= this->max_step_dual_ )
        {
            std::stringstream ss_dual;
            ss_dual << prefix_dual << time_step + 1;
            this->read_file ( *this->solD_, time_step + 1, filename_dual, groupname_dual, ss_dual.str ( ), -1 );
        }
        else
        {
            solD_->Zeros ( );
        }

        if ( this->is_stationary_ )
        {
            solP_->CloneFrom ( *solP_next_ );
            solD_->CloneFrom ( *solD_next_ );
        }

        // read next solutions
        if ( time_step + 2 <= this->max_step_primal_ )
        {
            std::stringstream ss_primal_next;
            ss_primal_next << prefix_primal << time_step + 2;

            this->read_file ( *solP_next_, time_step + 2, filename_primal, groupname_primal, ss_primal_next.str ( ), 1 );
        }
        else
        {
            solP_next_->Zeros ( );
        }
        if ( time_step + 2 <= this->max_step_dual_ )
        {
            std::stringstream ss_dual_next;
            ss_dual_next << prefix_dual << time_step + 2;

            std::stringstream ss_dual;
            ss_dual << prefix_dual << time_step + 1;

            this->read_file ( *this->solD_next_, time_step + 2, filename_dual, groupname_dual, ss_dual_next.str ( ), -1 );
        }
        else
        {
            solD_next_->Zeros ( );
        }

        // higher order interpolation
        patch_interpolation_.interpolate ( *solP_, fineP_ );
        patch_interpolation_.interpolate ( *solP_prev_, fineP_prev_ );
        patch_interpolation_.interpolate ( *solP_next_, fineP_next_ );

        patch_interpolation_dual_.interpolate ( *solD_, fineD_ );
        patch_interpolation_dual_.interpolate ( *solD_prev_, fineD_prev_ );
        patch_interpolation_dual_.interpolate ( *solD_next_, fineD_next_ );

        solP_ ->Update ( );
        solP_prev_ ->Update ( );
        solP_next_ ->Update ( );
        solD_ ->Update ( );
        solD_prev_ ->Update ( );
        solD_next_ ->Update ( );
        fineP_ .Update ( );
        fineP_prev_.Update ( );
        fineP_next_.Update ( );
        fineD_ .Update ( );
        fineD_prev_.Update ( );
        fineD_next_.Update ( );

        this->visualize_patch_interpolation ( time_step );

        // update assembler
        local_asm_est_->set_dT_pc ( dT_pc );
        local_asm_est_->set_dT_cn ( dT_cn );
        local_asm_est_->set_time_offset ( cur_time );
        local_asm_est_->setup_grid_search ( );
        local_asm_est_->setup_fe_evaluators ( );

        // CAUTION: modify this if spatial mesh is adapted in time
        local_asm_est_->set_unique_fine_space ( true );

        // ---------------------------------------
        // ---------------------------------------
        // assemble estimators for first subinterval

        // intra time interval indicators
        this->assemble_intra_time_indicators ( 0, est_cell_residual[time_step], est_interface[time_step] );
        this->indicator_mesh_indices_.push_back ( this->active_mesh_index_ );

        // write out mesh
        this->visualize_error_indicators ( time_step );

        // ---------------------------------------
        // ---------------------------------------
        // assemble estimators for second subinterval

        if ( !this->is_stationary_ )
        {
            // intra time interval indicators
            this->assemble_intra_time_indicators ( 1, est_cell_residual[time_step + 1], est_interface[time_step + 1] );
            this->indicator_mesh_indices_.push_back ( this->active_mesh_index_ );

            // write out mesh
            this->visualize_error_indicators ( time_step + 1 );
        }

        time_step += 2;
        cur_time += dT_pc + dT_cn;
    }

    // ---------------------------------------
    // assemble estimators for time jumps in dual solution
    // TODO
    if ( !this->is_stationary_ )
    {

    }

    std::vector< std::vector<double> > element_residual_h;
    std::vector< std::vector<double> > element_residual_tau;
    std::vector< std::vector<double> > element_weight_h;
    std::vector< std::vector<double> > element_weight_tau;

    this->accumulate_residuals_and_weights ( 0, est_cell_residual, est_cell_timejump, est_interface, element_residual_h, element_weight_h );
    this->accumulate_residuals_and_weights ( 1, est_cell_residual, est_cell_timejump, est_interface, element_residual_tau, element_weight_tau );

    this->compute_estimators ( element_residual_h, element_residual_tau, element_weight_h, element_weight_tau );

    // compute global error estimator
    MPI_Allreduce ( &local_time_estimator_, &global_time_estimator_, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    MPI_Allreduce ( &local_space_estimator_, &global_space_estimator_, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    this->global_estimator_ = this->global_time_estimator_ + this->global_space_estimator_;

    if ( rank_ == master_rank_ )
    {
        std::cout << " -------------------------------" << std::endl;
        std::cout << "  Global space error estimation: " << global_space_estimator_ << std::endl;
        std::cout << "  Global time  error estimation: " << global_time_estimator_ << std::endl;
        std::cout << "  Global       error estimation: " << global_estimator_ << std::endl;
        std::cout << " -------------------------------" << std::endl;

        string path = this->root_ + "/log/estimator.txt";
        ofstream out;
        out.open ( path.c_str ( ), ios::out | ios::app );
        out.precision ( 10 );
        out << std::scientific;
        out << this->adapt_counter_ << ": " << global_space_estimator_ << " " << global_time_estimator_ << " " << global_estimator_ << "\n";
        out.close ( );
    }

    timer.stop ( );
    if ( rank_ == master_rank_ )
    {
        std::cout << "> Time for error estimation " << timer.get_duration ( ) << " sec" << std::endl;
    }
}

void MetFlowApp::compute_estimators ( const std::vector< std::vector<double> >& element_residual_h,
                                      const std::vector< std::vector<double> >& element_residual_tau,
                                      const std::vector< std::vector<double> >& element_weight_h,
                                      const std::vector< std::vector<double> >& element_weight_tau )
{
    if ( rank_ == master_rank_ )
    {
        std::cout << "> Compute error estimators " << std::endl;
    }

    // reset indicators and estimators
    int num_steps = this->get_num_intervals ( );
    int num_mesh = this->dmh_->num_mesh ( );

    this->time_indicator_.clear ( );
    this->time_indicator_.resize ( num_steps );

    this->space_indicator_.clear ( );
    this->space_indicator_.resize ( num_steps );

    for ( int t = 0; t < num_steps; ++t )
    {
        int num_cells = this->dmh_->num_cells_by_index ( this->indicator_mesh_indices_[t] );
        this->time_indicator_[t].resize ( num_cells, 0. );
        this->space_indicator_[t].resize ( num_cells, 0. );
    }

    this->local_time_estimator_ = 0.;
    this->local_space_estimator_ = 0.;
    this->global_time_estimator_ = 0.;
    this->global_space_estimator_ = 0.;

    const double scaling = 0.5;

    // compute indicators for each space-time slab -----------------------------
    // loop over all time steps
    for ( int t = 0; t < num_steps; ++t )
    {
        int num_cells = this->dmh_->num_cells_by_index ( this->indicator_mesh_indices_[t] );

        // loop over all cells
        for ( int c = 0; c < num_cells; ++c )
        {
            // cell residual contribution
            this->space_indicator_[t][c] = scaling * element_residual_h[t][c] * element_weight_h[t][c];
            this->time_indicator_[t][c] = scaling * element_residual_tau[t][c] * element_weight_tau[t][c];
        }
        this->visualize_reduced_error_indicators ( t, this->indicator_mesh_indices_[t] );
    }

    // compute total error estimation
    for ( int t = 0; t < num_steps; ++t )
    {
        int mesh_index = this->indicator_mesh_indices_[t];
        int num_cells = this->dmh_->num_cells_by_index ( mesh_index );
        for ( int c = 0; c < num_cells; ++c )
        {
            if ( this->dmh_->get_mesh_by_index ( mesh_index )->cell_is_local ( c ) )
            {
                this->local_time_estimator_ += this->time_indicator_[t][c];
                this->local_space_estimator_ += this->space_indicator_[t][c];
            }
            else
            {
                //assert (this->space_indicator_[t][c] == 0.);
                //assert (this->time_indicator_[t][c] == 0.);
            }
        }
    }
}

void MetFlowApp::compute_reduced_estimators ( )
{
    if ( rank_ == master_rank_ )
    {
        std::cout << "> Compute reduced error estimators " << std::endl;
    }

    std::string reduction_type = base_params_["Adaptivity"]["Estimator"]["ReductionType"].get<std::string>( );

    // reset indicators and estimators
    int num_steps = this->get_num_intervals ( );
    int num_mesh = this->dmh_->num_mesh ( );

    this->reduced_space_indicator_.resize ( num_mesh );
    for ( int m = 0; m < num_mesh; ++m )
    {
        int num_cells = this->dmh_->num_cells_by_index ( m );
        this->reduced_space_indicator_[m].resize ( num_cells, 0. );
    }

    this->reduced_time_indicator_.resize ( num_steps, 0. );
    std::vector<double> local_reduced_time_indicator ( num_steps, 0. );

    // compute indicators for cells and time steps ---------------------------------------
    if ( reduction_type == "SUM" )
    {
        // loop over all time steps
        for ( int t = 0; t < num_steps; ++t )
        {
            double delta_t = this->get_delta_t ( t + 1 );
            int mesh_index = this->indicator_mesh_indices_[t];
            int num_cells = this->dmh_->num_cells_by_index ( mesh_index );

            //std::cout << " time interval " << t << " belongs to mesh " << mesh_index << std::endl;

            // loop over all cells
            for ( int c = 0; c < num_cells; ++c )
            {
                if ( this->dmh_->get_mesh_by_index ( mesh_index )->cell_is_local ( c ) )
                {
                    this->reduced_space_indicator_[mesh_index][c] += this->space_indicator_[t][c] /*/ delta_t*/;
                    local_reduced_time_indicator[t] += this->time_indicator_[t][c];
                }
            }
        }
        // compute global reduced time indicator
        MPI_Allreduce ( &local_reduced_time_indicator[0], &this->reduced_time_indicator_[0], num_steps, MPI_DOUBLE, MPI_SUM, comm_ );
    }
    else if ( reduction_type == "MAX" )
    {
        assert ( false );
        /*
        // time indicator
        for (int t=0; t<num_steps; ++t)
        {
            double max_est_cell = 0.;
            int mesh_index = mesh_indices[t];
            int num_cells  = this->get_num_cells(mesh_index, 0);

            for (int c=0; c<num_cells; ++c)
            {
                if (this->get_mesh(mesh_index, 0)->cell_is_local(c))
                {
                    if (this->time_indicator_[t][c] >= max_est_cell)
                    {
                        max_est_cell = this->time_indicator_[t][c];
                    }
                }
            }
            local_reduced_time_indicator[t] = max_est_cell;
        }

        // compute global reduced time indicator
        MPI_Allreduce ( &local_reduced_time_indicator[0], &this->reduced_time_indicator_[0], num_steps, MPI_DOUBLE, MPI_MAX, comm_ );

        // TODO: evtl buggy
        // space indicator
        for (int m=0; m<num_mesh; ++m)
        {
            int first_t   = this->get_first_step_for_mesh(m);
            int last_t    = this->get_last_step_for_mesh(m);
            int num_cells = this->get_num_cells(first_t);

            std::cout << "mesh " << m << " first time interval " << first_t << " last interval " << last_t << " num_cells " << num_cells << std::endl;
            for (int c=0; c<num_cells; ++c)
            {
                if (this->get_mesh(first_t)->cell_is_local(c))
                {
                    double max_est_interval = 0.;
                    for (int t=first_t; t<=last_t; ++t)
                    {
                        assert (m == mesh_indices[t]);
                        double delta_t = this->get_delta_t(t+1);
                        if (this->space_indicator_[t][c] / delta_t >= max_est_interval)
                        {
                            max_est_interval = this->space_indicator_[t][c] / delta_t;
                        }
                    }
                }
                this->reduced_space_indicator_[m][c] = max_est_interval;
            }
        }
         * */
    }

    this->visualize_reduced_space_indicators ( );
    if ( rank_ == master_rank_ )
    {
        this->visualize_reduced_time_indicators ( );
    }
}

void MetFlowApp::assemble_intra_time_indicators ( int rel_time,
                                                  std::vector< std::vector<double> >& est_cell_residual,
                                                  std::vector< std::vector<double> >& est_interface )
{
    local_asm_est_->set_est_rel_time ( rel_time );
    local_asm_est_->set_cell_mode ( RESIDUAL );

    MeshPtr active_mesh = this->dmh_->get_mesh_by_index ( this->active_mesh_index_ );
    VectorSpace<DATATYPE>* active_space = this->dmh_->get_space_by_index ( this->active_mesh_index_, 1 );

    // ---------------------------------------
    // cell contributions
    est_cell_residual.clear ( );
    global_asm_.assemble_multiple_scalar ( *active_space, boost::ref ( *local_asm_est_ ), 2 * this->num_eq_ * 4, est_cell_residual );

    std::vector< std::string> est_type ( 4 );
    est_type[0] = "res_h";
    est_type[1] = "weight_h";
    est_type[2] = "res_tau";
    est_type[3] = "weight_tau";

    int num_cell = est_cell_residual.size ( );

    for ( int q = 0; q<this->num_eq_; ++q )
    {
        // primal indicators
        for ( int e = 0; e < 4; ++e )
        {
            std::vector<double> attr_val ( num_cell, 0. );
            for ( int c = 0; c < num_cell; ++c )
            {
                attr_val[c] = est_cell_residual[c][q * 4 + e];
            }
            AttributePtr attr_est ( new DoubleAttribute ( attr_val ) );
            active_mesh->add_attribute ( "cell_primal_" + static_cast < ostringstream* > ( &( ostringstream ( ) << q ) )->str ( ) + "_" + est_type[e], DIM, attr_est );
        }

        // dual indicators
        for ( int e = 0; e < 4; ++e )
        {
            std::vector<double> attr_val ( num_cell, 0. );
            for ( int c = 0; c < num_cell; ++c )
            {
                attr_val[c] = est_cell_residual[c][( num_eq_ + q )*4 + e];
            }
            AttributePtr attr_est ( new DoubleAttribute ( attr_val ) );
            active_mesh->add_attribute ( "cell_dual_" + static_cast < ostringstream* > ( &( ostringstream ( ) << q ) )->str ( ) + "_" + est_type[e], DIM, attr_est );
        }
    }

    // ---------------------------------------
    // interface contributions
    est_interface.clear ( );
    jump_term_asm_.assemble_interface_multiple_scalar_cells ( *active_space, boost::ref ( *local_asm_est_ ), 2 * this->num_eq_ * 4, est_interface );
    num_cell = est_interface.size ( );

    for ( int q = 0; q<this->num_eq_; ++q )
    {
        // primal indicators
        for ( int e = 0; e < 4; ++e )
        {
            std::vector<double> attr_val ( num_cell, 0. );
            for ( int c = 0; c < num_cell; ++c )
            {
                attr_val[c] = est_interface[c][q * 4 + e];
            }
            AttributePtr attr_est ( new DoubleAttribute ( attr_val ) );
            active_mesh->add_attribute ( "jump_primal_" + static_cast < ostringstream* > ( &( ostringstream ( ) << q ) )->str ( ) + "_" + est_type[e], DIM, attr_est );
        }
        // dual indicators
        for ( int e = 0; e < 4; ++e )
        {
            std::vector<double> attr_val ( num_cell, 0. );
            for ( int c = 0; c < num_cell; ++c )
            {
                attr_val[c] = est_interface[c][( num_eq_ + q )*4 + e];
            }
            AttributePtr attr_est ( new DoubleAttribute ( attr_val ) );
            active_mesh->add_attribute ( "jump_dual_" + static_cast < ostringstream* > ( &( ostringstream ( ) << q ) )->str ( ) + "_" + est_type[e], DIM, attr_est );
        }
    }
}

void MetFlowApp::accumulate_residuals_and_weights ( int est_type,
                                                    const std::vector< std::vector< std::vector<double> > >& est_cell_residual,
                                                    const std::vector< std::vector< std::vector<double> > >& est_cell_timejump,
                                                    const std::vector< std::vector< std::vector<double> > >& est_interface,
                                                    std::vector< std::vector<double> >& element_residual,
                                                    std::vector< std::vector<double> >& element_weight )
{
    if ( rank_ == master_rank_ )
    {
        std::cout << "> Accumulate elementwise residuals and weights " << std::endl;
    }

    // Data structure of input vectors (est_cell_residual, est_cell_timejump, est_interface)
    // [time_step][entity_index] { {primal_eq_1_res_h, primal_eq_1_weight_h, primal_eq_1_res_tau, primal_eq_1_weight_tau;
    //                              primal_eq_n_res_h, primal_eq_n_weight_h, primal_eq_n_res_tau, primal_eq_n_weight_tau}
    //                                  ........
    //                             {dual_eq_1_res_h, dual_eq_1_weight_h, dual_eq_1_res_tau, dual_eq_1_weight_tau;
    //                              dual_eq_n_res_h, dual_eq_n_weight_h, dual_eq_n_res_tau, dual_eq_n_weight_tau}}

    int num_steps = this->get_num_intervals ( );
    int num_mesh = this->dmh_->num_mesh ( );

    bool csi = base_params_["Adaptivity"]["Estimator"]["UseCauchySchwarz"].get<int>( 0 );

    std::vector< std::vector<double> > hK ( num_mesh );
    for ( int m = 0; m < num_mesh; ++m )
    {
        int num_cells = this->dmh_->get_mesh_by_index ( m )->num_entities ( DIM );
        hK[m].resize ( num_cells, 1. );
        if ( csi )
        {
            this->compute_cell_diameters ( hK[m] );
        }
    }

    // est_type = 0: space error estimation
    // est_type = 1: time error estimation
    int off = 0;
    if ( est_type == 1 )
    {
        off = 2;
    }

    // allocate arrays
    element_residual.clear ( );
    element_weight.clear ( );

    element_residual.resize ( num_steps );
    element_weight.resize ( num_steps );

    for ( int t = 0; t < num_steps; ++t )
    {
        int num_cells = est_cell_residual[t].size ( );
        element_residual[t].resize ( num_cells, 0. );
        element_weight[t].resize ( num_cells, 0. );
    }

    // combine indicators
    // loop over all time steps
    for ( int t = 0; t < num_steps; ++t )
    {
        int num_cells = est_cell_residual[t].size ( );
        int mesh_index = this->indicator_mesh_indices_[t];

        // loop over all cells
        for ( int c = 0; c < num_cells; ++c )
        {
            const double h = hK[mesh_index][c];
            const double inv_h = 1. / h;

            // loop over all equations
            for ( int e = 0; e<this->num_eq_; ++e )
            {
                // element residual
                element_residual[t][c] += est_cell_residual[t][c][4 * e + off] // primal cell residual ||res||_K^{2} (csi) or (res_K, weight)_K (!csi)
                        + est_cell_residual[t][c][4 * ( e + num_eq_ ) + off] // dual cell residual ||res||_K^{2} (csi) or (res_K, weight)_K (!csi)
                        + inv_h * est_interface[t][c][4 * e + off] // primal interface residual ||res||_dK^{2} (csi) or (res_K, weight)_K (!csi)
                        + inv_h * est_interface[t][c][4 * ( e + num_eq_ ) + off]; // dual interface residual   ||res||_dK^{2} (csi) or (res_K, weight)_K (!csi)
                // + est_cell_timejump[t][c][4* e         +off]      // primal cell timejump ||res||_K^{2} (csi) or (res_K, weight)_K (!csi)
                // + est_cell_timejump[t][c][4*(e+num_eq_)+off]      // dual cell timejump ||res||_K^{2} (csi) or (res_K, weight)_K (!csi)

                // element weight
                if ( csi )
                {
                    element_weight[t][c] += est_cell_residual[t][c][4 * e + off + 1] // primal cell residual ||weight||_K^{2} (csi) or 1 (!csi)
                            + est_cell_residual[t][c][4 * ( e + num_eq_ ) + off + 1] // dual cell residual ||weight||_K^{2} (csi) or 1 (!csi)
                            + h * est_interface[t][c][4 * e + off + 1] // primal interface residual ||weight||_dK^{2} (csi) or 1 (!csi)
                            + h * est_interface[t][c][4 * ( e + num_eq_ ) + off + 1]; // dual interface residual   ||weight||_dK^{2} (csi) or 1 (!csi)
                    // + est_cell_timejump[t][c][4* e         +off+1]      // primal cell timejump |weight||_K^{2} (csi) or 1 (!csi)
                    // + est_cell_timejump[t][c][4*(e+num_eq_)+off+1]      // dual cell timejump ||weight||_K^{2} (csi) or 1 (!csi)
                }
            }
            if ( csi )
            {
                element_residual[t][c] = std::sqrt ( element_residual[t][c] );
                element_weight[t][c] = std::sqrt ( element_weight[t][c] );
            }
            else
            {
                element_residual[t][c] = std::abs ( element_residual[t][c] );
                element_weight[t][c] = 1.;
            }
        }
    }
}

void MetFlowApp::adapt_temporal_mesh ( )
{
    std::string type = base_params_["Adaptivity"]["TemporalAdaption"]["Type"].get<std::string>( );
    int coarsen_marker = base_params_["Adaptivity"]["TemporalAdaption"]["CoarsenMarker"].get<int>( );

    int num_intervals = this->t_mesh_.num_intervals ( this->adapt_counter_ );
    assert ( num_intervals % 2 == 0 );

    std::vector<int> adapt_markers ( num_intervals, 0 );

    // determine which intervals should be refined
    if ( type == "None" )
    {
        adapt_markers.clear ( );
        adapt_markers.resize ( num_intervals, 0 );
    }
    else if ( type == "Uniform" )
    {
        adapt_markers.clear ( );
        adapt_markers.resize ( num_intervals, 1 );
    }
    else if ( type == "FixedFraction" )
    {
        double refine_frac = base_params_["Adaptivity"]["TemporalAdaption"]["FixedFraction"]["FractionToRefine"].get<double>( 0.2 );
        double coarsen_frac = base_params_["Adaptivity"]["TemporalAdaption"]["FixedFraction"]["FractionToCoarsen"].get<double>( 0.1 );
        int threshold = base_params_["Adaptivity"]["TemporalAdaption"]["FixedFraction"]["CoarsenThreshold"].get<int>( 100 );

        local_fixed_fraction_strategy ( refine_frac, coarsen_frac, threshold, coarsen_marker, this->reduced_time_indicator_, adapt_markers );
    }
    else if ( type == "FixedError" )
    {
        double tol = base_params_["Adaptivity"]["TemporalAdaption"]["FixedError"]["Tolerance"].get<double>( 1e-3 );
        double conv_order = base_params_["Adaptivity"]["TemporalAdaption"]["FixedError"]["ConvergenceOrder"].get<double>( 1. );
        int threshold = base_params_["Adaptivity"]["TemporalAdaption"]["FixedError"]["CoarsenThreshold"].get<int>( 100 );

        fixed_error_strategy ( tol, num_intervals, conv_order, threshold, coarsen_marker, this->reduced_time_indicator_, adapt_markers );
    }
    assert ( adapt_markers.size ( ) == num_intervals );

    // refine
    this->t_mesh_.refine ( adapt_markers, this->adapt_counter_ );
    int index = this->adapt_counter_ + 1;

    // write out mesh
    if ( rank_ == master_rank_ )
    {
        std::string prefix = "time_mesh";

        std::stringstream pre;
        pre << this->root_ << "/mesh/" << prefix;
        std::string visu_filename = pre.str ( );
        std::cout << "> Save time mesh as " << visu_filename << std::endl;

        this->t_mesh_.save ( visu_filename, index );

        std::stringstream log_pre;
        log_pre << this->root_ << "/log/time_steps.csv";
        std::string log_filename = log_pre.str ( );

        ofstream logfile;
        logfile.open ( log_filename.c_str ( ), std::ios_base::app );

        for ( int t = 0; t<this->t_mesh_.num_intervals ( index ) + 1; ++t )
        {
            logfile << this->t_mesh_.time ( t, index ) << ", ";
        }
        logfile << "\n";
        for ( int t = 0; t<this->t_mesh_.num_intervals ( index ); ++t )
        {
            logfile << this->t_mesh_.time ( t + 1, index ) - this->t_mesh_.time ( t, index ) << ", ";
        }
        logfile << "\n";
        logfile.close ( );
    }
}

void MetFlowApp::adapt_spatial_mesh ( )
{
    std::string type = base_params_["Adaptivity"]["SpatialAdaption"]["Type"].get<std::string>( );
    int coarsen_marker = base_params_["Adaptivity"]["SpatialAdaption"]["CoarsenMarker"].get<int>( );
    int num_mesh = this->dmh_->num_mesh ( );

    for ( int m = 0; m < num_mesh; ++m )
    {
        if ( rank ( ) == master_rank ( ) )
        {
            std::cout << "> Adapt mesh " << m << std::endl;
        }

        int num_cells = this->dmh_->num_cells_by_index ( m );
        MeshPtr active_mesh = this->dmh_->get_mesh_by_index ( m );

        std::vector<int> adapt_markers ( num_cells, 0 );

        // determine which intervals should be refined
        if ( type == "None" )
        {
            adapt_markers.clear ( );
            adapt_markers.resize ( num_cells, 0 );
        }
        else if ( type == "Uniform" )
        {
            adapt_markers.clear ( );
            adapt_markers.resize ( num_cells, 1 );
        }
        else if ( type == "FixedFraction" )
        {
            double refine_frac = base_params_["Adaptivity"]["SpatialAdaption"]["FixedFraction"]["FractionToRefine"].get<double>( 0.2 );
            double coarsen_frac = base_params_["Adaptivity"]["SpatialAdaption"]["FixedFraction"]["FractionToCoarsen"].get<double>( 0.1 );
            int threshold = base_params_["Adaptivity"]["SpatialAdaption"]["FixedFraction"]["CoarsenThreshold"].get<int>( 100 );

            local_fixed_fraction_strategy ( refine_frac, coarsen_frac, threshold, coarsen_marker, this->reduced_space_indicator_[m], adapt_markers );
        }
        else if ( type == "FixedError" )
        {
            double tol = base_params_["Adaptivity"]["SpatialAdaption"]["FixedError"]["Tolerance"].get<double>( 1e-3 );
            double conv_order = base_params_["Adaptivity"]["SpatialAdaption"]["FixedError"]["ConvergenceOrder"].get<double>( 1. );
            int threshold = base_params_["Adaptivity"]["SpatialAdaption"]["FixedError"]["CoarsenThreshold"].get<int>( 100 );
            int num_global_cells = active_mesh->num_global_cells ( this->comm_ );

            fixed_error_strategy ( tol, num_global_cells, conv_order, threshold, coarsen_marker, this->reduced_space_indicator_[m], adapt_markers );
        }
        assert ( adapt_markers.size ( ) == num_cells );

        // refine mesh
        this->dmh_->refine ( m, adapt_markers );
    }
    std::string prefix = this->root_ + "/mesh/space_mesh";
    this->dmh_->save_mesh_to_file ( this->adapt_counter_ + 1, prefix );
    this->dmh_->visualize_mesh ( this->adapt_counter_ + 1, prefix );
}

void MetFlowApp::adapt_mesh_change_list ( )
{
    if ( rank_ == master_rank_ )
    {
        std::cout << "> Adapt mesh list " << std::endl;
    }
    std::string adapt_type = base_params_["Adaptivity"]["DynamicMesh"]["Type"].get<std::string>( );
    int min_steps_for_mesh = base_params_["Adaptivity"]["DynamicMesh"]["MinStepsBetweenChanges"].get<int>( );
    int max_mesh_number = base_params_["Adaptivity"]["DynamicMesh"]["MaxMeshNumber"].get<int>( );
    int rate = base_params_["Adaptivity"]["DynamicMesh"]["FixedFraction"]["GrowthRate"].get<int>( );
    double tol = base_params_["Adaptivity"]["DynamicMesh"]["FixedError"]["RelChangeThreshold"].get<double>( );

    DynamicMeshProblem<LAD, MESHIMPL, DIM>::adapt_mesh_change_list ( adapt_type, this->adapt_counter_, min_steps_for_mesh, max_mesh_number, rate, tol,
                                                                     this->space_indicator_, this->indicator_mesh_indices_ );

    if ( rank ( ) == master_rank ( ) )
    {
        std::vector<double> mesh_change_times = this->dmh_->get_mesh_change_times ( );
        this->write_mesh_change_list ( this->adapt_counter_ + 1, mesh_change_times );
    }
}

void MetFlowApp::prepare_higher_order_space ( )
{
    Timer timer;
    timer.start ( );

    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "> Higher order interpolation for active mesh " << this->active_mesh_index_ << std::endl;
    }

    assert ( this->mesh_->is_uniformly_coarsenable ( ) );

    // init interpolating space and interpolation map
    patch_interpolation_.init ( this->space_ );
    patch_interpolation_dual_.init ( this->space_dual_ );

    // get interpolating space
    fine_space_ = patch_interpolation_.get_space ( );
    fine_space_dual_ = patch_interpolation_dual_.get_space ( );
    const Couplings<double>& inter_couplings = patch_interpolation_.get_couplings ( );
    const Couplings<double>& inter_couplings_dual = patch_interpolation_dual_.get_couplings ( );

    // interpolate solution
#ifdef WITH_HYPRE
    fineP_.Init ( comm_, inter_couplings );
    fineP_prev_.Init ( comm_, inter_couplings );
    fineP_next_.Init ( comm_, inter_couplings );

    fineD_.Init ( comm_, inter_couplings_dual );
    fineD_prev_.Init ( comm_, inter_couplings_dual );
    fineD_next_.Init ( comm_, inter_couplings_dual );

#else
    fineP_.Init ( comm_, inter_couplings, CPU, NAIVE );
    fineP_.InitStructure ( );
    fineP_prev_.Init ( comm_, inter_couplings, CPU, NAIVE );
    fineP_prev_.InitStructure ( );
    fineP_next_.Init ( comm_, inter_couplings, CPU, NAIVE );
    fineP_next_.InitStructure ( );

    fineD_.Init ( comm_, inter_couplings_dual, CPU, NAIVE );
    fineD_.InitStructure ( );
    fineD_prev_.Init ( comm_, inter_couplings_dual, CPU, NAIVE );
    fineD_prev_.InitStructure ( );
    fineD_next_.Init ( comm_, inter_couplings_dual, CPU, NAIVE );
    fineD_next_.InitStructure ( );
#endif

    fineP_.Zeros ( );
    fineP_prev_.Zeros ( );
    fineP_next_.Zeros ( );

    fineD_.Zeros ( );
    fineD_prev_.Zeros ( );
    fineD_next_.Zeros ( );

    local_asm_est_->set_fineP ( fineP_ );
    local_asm_est_->set_fineP_prev ( fineP_prev_ );
    local_asm_est_->set_fineP_next ( fineP_next_ );
    local_asm_est_->set_fineD ( fineD_ );
    local_asm_est_->set_fineD_next ( fineD_next_ );
    local_asm_est_->set_fineD_prev ( fineD_prev_ );

    local_asm_est_->set_fineP_space ( *fine_space_ );
    local_asm_est_->set_fineP_space_next ( *fine_space_ );
    local_asm_est_->set_fineP_space_prev ( *fine_space_ );
    local_asm_est_->set_fineD_space ( *fine_space_dual_ );
    local_asm_est_->set_fineD_space_next ( *fine_space_dual_ );
    local_asm_est_->set_fineD_space_prev ( *fine_space_dual_ );

    timer.stop ( );

    if ( rank ( ) == master_rank ( ) )
    {
        std::cout << "  Patch interpolation time " << timer.get_duration ( ) << std::endl;
    }
}
