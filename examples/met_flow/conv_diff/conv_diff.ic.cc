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

/// ***************************************************************************
/// INITIAL CONDITION
/// ***************************************************************************
// Create file name for HDF5 of start vector

void MetFlowConvDiffApp::create_ic_hdf5_filename ( )
{
    const int temperature_deg = params_["FiniteElements"]["TemperatureDegree"].get<int>( );

    std::string prefix = this->base_params_["InitialCondition"]["FilePrefix"].get<std::string>( );
    std::string mesh_name = this->base_params_["Mesh"]["Filename"].get<std::string>( );
    mesh_name.erase ( 0, 3 );
    int ref_level = this->base_params_["Mesh"]["InitialRefLevel"].get<int>( );
    this->filename_start_ = this->root_ + "/" + prefix + "-" + mesh_name + "-" + static_cast < ostringstream* > ( &( ostringstream ( ) << ref_level ) )->str ( ) + "-"
            + static_cast < ostringstream* > ( &( ostringstream ( ) << this->num_partitions_ ) )->str ( ) + "-"
            + static_cast < ostringstream* > ( &( ostringstream ( ) << temperature_deg ) )->str ( );

}

// write initial condition into HDF5 file

void MetFlowConvDiffApp::write_ic ( )
{
    this->create_ic_hdf5_filename ( );

    std::stringstream ss;
    std::string groupname_start = "start";
    std::string prefix_start = "start";
    ss << prefix_start;
    double time = 0.;
    MetFlowApp::write_file ( *this->solP_prev_, this->filename_start_, groupname_start, ss.str ( ) );
    std::stringstream ss_base;
    ss_base << "base";
    MetFlowApp::write_file ( *this->base_, this->filename_start_, groupname_start, ss_base.str ( ) );
}

// compute initial codnition

void MetFlowConvDiffApp::prepare_ic ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare initial conditions for TEHD problem " << std::endl;

    int perturb_ic = this->base_params_["Perturbation"]["PerturbIC"].get<int>( 0 );

    std::string starting_type = this->base_params_["InitialCondition"]["Type"].get<std::string>( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  StartingType: " << starting_type.c_str ( ) << std::endl;

    this->update_convection ( 0., false );
    // enable the counting of jacobian assemblies for efficient preconditioning
    this->ilupp_time_step_ = 0;

    // **********************************************
    // Load IC
    // **********************************************
    if ( starting_type == "Load" )
    {
        // Load start solution and return
        this->create_ic_hdf5_filename ( );

        std::stringstream ss;
        std::stringstream ss_base;
        std::string groupname_start = "start";
        std::string prefix_start = "start";
        ss << prefix_start;
        ss_base << "base";
        double cur_time;
        this->read_file ( *this->solP_, this->filename_start_, groupname_start, ss.str ( ) );
        this->read_file ( *this->base_, this->filename_start_, groupname_start, ss_base.str ( ) );

        // perturb initial solution
        if ( perturb_ic > 0 )
        {
            MetFlowConvDiffApp::perturb_ic ( );
        }
        return;
    }

    // **********************************************
    // set initial (Zero) fields and boundary conditions
    // **********************************************
    if ( rank ( ) == master_rank ( ) ) std::cout << " > Set up velocity and temperature field" << std::endl;

    for ( int var = 0; var < this->num_vars_; ++var )
    {
        const int tdim = this->space_->mesh ( ).tdim ( );
        for ( mesh::EntityIterator it = this->space_->mesh ( ).begin ( tdim ), end_it = this->space_->mesh ( ).end ( tdim ); it != end_it; ++it )
        {
            std::vector<int> global_dof_ids;
            this->space_->GetDofIndices ( var, *it, &global_dof_ids );
            int num_dofs = global_dof_ids.size ( );
            std::vector<double> values;
            values.resize ( num_dofs, 0. );

            if ( var == this->t_var_ )
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

    MetFlowConvDiffApp::prepare_bc ( );
    MetFlowApp::set_bc ( 1 );

    if ( starting_type == "Constant" )
    {
        this->base_->CloneFrom ( *this->solP_ );

        // perturb initial solution
        if ( perturb_ic > 0 )
        {
            MetFlowConvDiffApp::perturb_ic ( );
        }
        return;
    }
}

// perturb initial condition TODO scale (+-1) checken

void MetFlowConvDiffApp::perturb_ic ( )
{
}

// compute z(T) = solution of Mz=j^2_u
// TODO evtl ändern

void MetFlowConvDiffApp::prepare_ic_dual ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare initial conditions for dual TEHD problem " << std::endl;

    this->update_convection ( this->cur_time_, true );

    // Assmeble mass matrix for complete T dual space
    if ( this->computed_T_mass_dual_ == false )
    {
        SparsityStructure sparsity;
        this->global_asm_.compute_sparsity_structure ( *this->space_dual_, sparsity );

#ifdef USE_HYPRE
        T_mass_matrix_dual_.Init ( this->comm_, *this->dmh_->get_couplings_by_index ( this->dmh_->num_mesh ( ) - 1, -1 ) );
#else
        T_mass_matrix_dual_.Init ( this->comm_, *this->dmh_->get_couplings_by_index ( this->dmh_->num_mesh ( ) - 1, -1 ), this->la_sys_.Platform, MetFlowApp::APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
#endif
        T_mass_matrix_dual_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                                            vec2ptr ( sparsity.diagonal_cols ),
                                            sparsity.diagonal_rows.size ( ),
                                            vec2ptr ( sparsity.off_diagonal_rows ),
                                            vec2ptr ( sparsity.off_diagonal_cols ),
                                            sparsity.off_diagonal_rows.size ( ) );

        this->global_asm_.assemble_matrix ( *this->space_dual_, boost::ref ( T_local_asm_mass_ ), T_mass_matrix_dual_ );

        if ( !this->dirichlet_dofs_dual_.empty ( ) )
        {
            T_mass_matrix_dual_.diagonalize_rows ( vec2ptr ( this->dirichlet_dofs_dual_ ), this->dirichlet_dofs_dual_.size ( ), 1.0 );
        }
        this->computed_T_mass_dual_ = true;
    }

    // Assemble right hand side vector containing contributions of j(T)
    this->rhs_dual_->Zeros ( );

    T_local_asm_.set_solP ( *this->solP_ );
    T_local_asm_.set_solP_prev ( *this->solP_prev_ );
    T_local_asm_.set_solP_next ( *this->solP_next_ );
    T_local_asm_.set_vector_mode_to_goal ( );

    this->global_asm_.assemble_vector ( *this->space_dual_, boost::ref ( T_local_asm_ ), *this->rhs_dual_ );
    T_local_asm_.set_vector_mode_to_std ( );

    this->rhs_dual_->SetValues ( vec2ptr ( this->dirichlet_dofs_dual_ ), this->dirichlet_dofs_dual_.size ( ), vec2ptr ( this->dirichlet_values_dual_ ) );
    this->solD_next_->SetValues ( vec2ptr ( this->dirichlet_dofs_dual_ ), this->dirichlet_dofs_dual_.size ( ), vec2ptr ( this->dirichlet_values_dual_ ) );

    this->rhs_dual_->Update ( );
    this->solD_next_->Update ( );

    /*
        std::string prefix = params_["Visualization"]["DualFilePrefix"].get<std::string>();
        std::stringstream pre;
        pre << this->root_ << "/" << prefix << "_RHS." << this->adapt_counter_;
        this->visualize_solution_dual(T_rhs_dual_, pre.str(), 99);
     */

#ifdef USE_HYPRE
    T_mass_precond_.Clear ( );
    T_mass_precond_.SetPreconditioningParameters ( );
    T_mass_precond_.SetNumFunctions ( 1 );
    T_mass_precond_.SetCycleType ( params_["AMG"]["CycleType"].get<int>( ) );
    T_mass_precond_.InitControl ( params_["AMG"]["MaxIterations"].get<int>( ), 0, 0 );
    T_mass_precond_.SetRelaxType ( params_["AMG"]["RelaxType"].get<int>( ) );
    T_mass_precond_.SetRelaxWt ( params_["AMG"]["RelaxWeight"].get<DATATYPE>( ) );
    T_mass_precond_.SetInterpType ( params_["AMG"]["InterpolationType"].get<int>( ) );
    T_mass_precond_.SetStrongThreshold ( params_["AMG"]["StrongThreshold"].get<DATATYPE>( ) );
    T_mass_precond_.SetAggNumLevels ( params_["AMG"]["AggNumLevels"].get<int>( ) );
    T_mass_precond_.SetCoarsenType ( params_["AMG"]["CoarsenType"].get<int>( ) );
    T_mass_precond_.SetCycleNumSweeps ( params_["AMG"]["NumDownSweeps"].get<int>( ), 1 );
    T_mass_precond_.SetCycleNumSweeps ( params_["AMG"]["NumUpSweeps"].get<int>( ), 2 );
    T_mass_precond_.SetSmoothType ( params_["AMG"]["SmoothType"].get<int>( ) );
    T_mass_precond_.SetSmoothNumLevels ( params_["AMG"]["SmoothNumLevels"].get<int>( ) );
    T_mass_precond_.SetMaxCoarseSize ( params_["AMG"]["MaxCoarseSize"].get<int>( ) );
    T_mass_precond_.SetMaxLevels ( params_["AMG"]["MaxLevels"].get<int>( ) );
    //T_mass_precond_.SetNodal(1);
    //T_mass_precond_.SetNodalDiag(1);
    T_mass_precond_.Init ( );
    T_mass_solver_.InitControl ( 100, 1e-12, 1e-6, 1e6 );

    T_mass_solver_.SetupPreconditioner ( T_mass_precond_ );
    T_mass_solver_.InitParameter ( "RightPreconditioning" );
    T_mass_solver_.SetupOperator ( T_mass_matrix_dual_ );
#else
#    ifdef WITH_ILUPP
    const int prepro_type = this->base_params_["ILUPP"]["PreprocessingType"].get<int>( );
    const int precond_no = this->base_params_["ILUPP"]["PreconditionerNumber"].get<int>( );
    const int max_levels = this->base_params_["ILUPP"]["MaxMultilevels"].get<int>( );
    const double mem_factor = this->base_params_["ILUPP"]["MemFactor"].get<double>( );
    const double threshold = this->base_params_["ILUPP"]["PivotThreshold"].get<double>( );
    const double min_pivot = this->base_params_["ILUPP"]["MinPivot"].get<double>( );
    T_mass_precond_.InitParameter ( prepro_type, precond_no, max_levels, mem_factor, threshold, min_pivot );
    T_mass_precond_.SetupOperator ( T_mass_matrix_dual_ );
    T_mass_solver_.InitControl ( 100, 1e-12, 1e-6, 1e6 );

    T_mass_solver_.InitControl ( 100, 1e-20, 1e-10, 1e6 );
    T_mass_solver_.InitParameter ( 100, "RightPreconditioning" );
    T_mass_solver_.SetupPreconditioner ( T_mass_precond_ );
    T_mass_precond_.SetupOperator ( T_mass_matrix_dual_ );
#    else
    T_mass_solver_.InitParameter ( 100, "RightPreconditioning" );
#    endif
    T_mass_solver_.SetupOperator ( T_mass_matrix_dual_ );

#endif

    double norm = this->rhs_dual_->Norm2 ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  norm of rhs " << norm << std::endl;
    LinearSolverState state = T_mass_solver_.Solve ( *this->rhs_dual_, this->solD_next_ );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  done after " << T_mass_solver_.iter ( ) << " iterations with residual norm " << T_mass_solver_.res ( ) << std::endl;

    this->solD_next_->Update ( );

    //    this->visualize_solution_dual(T_solD_next_, pre.str(), 100);

}
