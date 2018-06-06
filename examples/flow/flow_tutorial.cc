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

/// \author Staffan Ronnas

#include <algorithm>
#include <vector>

#include "flow_tutorial.h"

#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif

// if you work with extra makefile BASEDIR is defined
// static const char* DATADIR = BASEDIR "/meshes/";
static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;
static const char* PARAM_FILENAME = "flow_tutorial.xml";

// program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    // set default parameter file
    std::string param_filename ( PARAM_FILENAME );
    std::string path_mesh;
    // if set take parameter file specified on console
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }
    // if set take mesh following path specified on console
    if ( argc > 2 )
    {
        path_mesh = std::string ( argv[2] );
    }
    try
    {
        FlowTutorial app ( param_filename, path_mesh );
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

FlowTutorial::FlowTutorial ( const std::string& param_filename, const std::string& path_mesh )
: path_mesh ( path_mesh ),
comm_ ( MPI_COMM_WORLD ),
params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
use_pressure_filter_ ( false ),
refinement_level_ ( 0 ),
matrix_ ( 0 ), sol_ ( 0 ), prev_sol_ ( 0 ), res_ ( 0 ), rhs_ ( 0 ),
is_done_ ( false )
{
}

FlowTutorial::~FlowTutorial ( )
{
    delete sol_;
    delete prev_sol_;
    delete res_;
    delete rhs_;
    delete matrix_;
}

void FlowTutorial::run ( )
{
    simul_name_ = params_["OutputPrefix"].get<std::string>( );

    std::ofstream info_log ( ( simul_name_ + "_info_log" ).c_str ( ) );
    LogKeeper::get_log ( "info" ).set_target ( &info_log );
    std::ofstream debug_log ( ( simul_name_ + "_debug_log" ).c_str ( ) );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

    // output parameters for debugging
    LOG_INFO ( "parameters", params_ );

    flow_model_ = params_["FlowModel"]["Type"].get<std::string>( );

    if ( flow_model_ != std::string ( "Channel" ) &&
         flow_model_ != std::string ( "Cavity" ) )
    {
        throw UnexpectedParameterValue ( "FlowModel.Type", flow_model_ );
    }

    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    read_mesh ( );

    // The simulation has two modes: stationary and
    // instationary. Which one is used, depends on the parameter
    // Instationary.SolveInstationary. In stationary mode, solve()
    // is only called once, whereas in instationary mode, it is
    // called several times via the time-stepping method implemented in
    // run_time_loop() .
    solve_instationary_ =
            params_["Instationary"]["SolveInstationary"].get<bool>( );
    while ( !is_done_ )
    {
        prepare ( );
        if ( solve_instationary_ )
        {
            LOG_INFO ( "simulation", "Solving instationary problem" );
            run_time_loop ( );
        }
        else
        {
            LOG_INFO ( "simulation", "Solving stationary problem" );
            solve ( );
            visualize ( );
        }
        adapt ( );
    }
}

void FlowTutorial::EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F )
{
    // compute the residual vector
    compute_residual ( u, F );
}

void FlowTutorial::EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
    // assemble the matrix for the linearized system
    compute_jacobian ( u, DF );
}

void FlowTutorial::read_mesh ( )
{
    if ( rank ( ) == MASTER_RANK )
    {
        const std::string mesh_name =
                params_["Mesh"]["Filename"].get<std::string>( );
        std::string mesh_filename;
        if ( path_mesh.empty ( ) )
        {
            mesh_filename = std::string ( DATADIR ) + mesh_name;
        }
        else
        {
            mesh_filename = path_mesh + mesh_name;
        }

        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0 );

        refinement_level_ = 0;
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );

        for ( int r = 0; r < initial_ref_lvl; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
            ++refinement_level_;
        }
        LOG_INFO ( "mesh", "Initial refinement level = " << refinement_level_ );
    }

    MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK,
                                                    comm_ );
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

void FlowTutorial::prepare ( )
{
    // prepare timestep
    ts_ = 0;
    dt_ = params_["Instationary"]["Timestep"].get<double>( );

    // Set the alpha coefficients correctly for the
    // Crank-Nicolson method.
    alpha1_ = 0.5 * dt_;
    alpha2_ = dt_;
    alpha3_ = 0.5 * dt_;

    // prepare problem parameters
    rho_ = params_["FlowModel"]["Density"].get<double>( );
    nu_ = params_["FlowModel"]["Viscosity"].get<double>( );

    if ( flow_model_ == std::string ( "Channel" ) )
    {
        Um_ = params_["FlowModel"]["InflowSpeed"].get<double>( );
        H_ = params_["FlowModel"]["InflowHeight"].get<double>( );
#if DIMENSION == 3
        W_ = params_["FlowModel"]["InflowWidth"].get<double>( );
#endif
    }
    else if ( flow_model_ == std::string ( "Cavity" ) )
    {
        Um_ = params_["FlowModel"]["LidSpeed"].get<double>( );
        std::cout << "Reynold's number = " << 1. / nu_ << "\n";
    }

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

    // pressure filter
    use_pressure_filter_ = params_["UsePressureFilter"].get<bool>( );

    // prepare linear algebra structures
    couplings_.Clear ( );
    couplings_.Init ( communicator ( ), space_.dof ( ) );

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

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    matrix_->Init ( communicator ( ), couplings_ );
    CoupledVectorFactory<Scalar> CoupVecFact;
    sol_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    sol_->Init ( communicator ( ), couplings_ );
    prev_sol_ = CoupVecFact.Get (
                                  params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    prev_sol_->Init ( communicator ( ), couplings_ );
    res_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    res_->Init ( communicator ( ), couplings_ );

    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );
    matrix_->Zeros ( );

    sol_->InitStructure ( );
    sol_->Zeros ( );

    prev_sol_->InitStructure ( );
    prev_sol_->Zeros ( );

    rhs_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    rhs_->Init ( communicator ( ), couplings_ );
    rhs_->InitStructure ( );
    rhs_->Zeros ( );

    res_->InitStructure ( );
    res_->Zeros ( );

    // prepare dirichlet BC
    prepare_bc ( );
}

void FlowTutorial::prepare_bc ( )
{
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    if ( flow_model_ == std::string ( "Channel" ) )
    {
        const int inflow_bdy = params_["Boundary"]["InflowMaterial"].get<int>( );
        const int outflow_bdy = params_["Boundary"]["OutflowMaterial"].get<int>( );
#if DIMENSION == 2
        ChannelFlowBC2d bc[2] = { ChannelFlowBC2d ( 0, H_, Um_, inflow_bdy, outflow_bdy ),
                                 ChannelFlowBC2d ( 1, H_, Um_, inflow_bdy, outflow_bdy ) };
#endif

#if DIMENSION == 3
        ChannelFlowBC3d bc[3] = { ChannelFlowBC3d ( 0, W_, H_, Um_, inflow_bdy, outflow_bdy ),
                                 ChannelFlowBC3d ( 1, W_, H_, Um_, inflow_bdy, outflow_bdy ),
                                 ChannelFlowBC3d ( 2, W_, H_, Um_, inflow_bdy, outflow_bdy ) };
#endif
        for ( int var = 0; var < DIMENSION; ++var )
        {
            compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                                dirichlet_dofs_, dirichlet_values_ );
        }
    }
    else if ( flow_model_ == std::string ( "Cavity" ) )
    {
        const int lid_bdy = params_["Boundary"]["LidMaterial"].get<int>( );

        LidCavityBC bc[2] = { LidCavityBC ( 0, Um_, lid_bdy ),
                             LidCavityBC ( 1, Um_, lid_bdy ) };

        for ( int var = 0; var < DIMENSION; ++var )
        {
            compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                                dirichlet_dofs_, dirichlet_values_ );
        }
    }
    else
    {
        assert ( false );
    }

    // apply BC to initial solution
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
}

void FlowTutorial::solve ( )
{
    // setup linear solver
    LinearSolver<LAD>* linsolver_;
    LinearSolverFactory<LAD> LinSolFact;
    linsolver_ = LinSolFact.Get (
                                  params_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( params_["LinearSolver"] );
    linsolver_->SetupOperator ( *matrix_ );

    // setup nonlinear solver parameters
    NonlinearSolverFactory<LAD> NLSFact;
    nls_ = NLSFact.Get ( params_["NonlinearSolver"]["Name"].get<std::string>( ) )->
            params ( res_, matrix_, params_["NonlinearSolver"] );

    // we use our own initial solution -- this needs to be indicated
    // to the Newton-solver
    if ( params_["NonlinearSolver"]["Name"].get<std::string>( ) == "Newton" )
        ( ( Newton<LAD>* )nls_ )->InitParameter ( Newton<LAD>::NewtonInitialSolutionOwn );

    nls_->SetOperator ( *this );
    nls_->SetLinearSolver ( *linsolver_ );

    NonlinearSolverState state = nls_->Solve ( *rhs_, sol_ );

    std::cout << "Nonlinear solver ended with state " << state
            << " and residual norm " << nls_->GetResidual ( )
            << " after " << nls_->iter ( ) << " iterations\n";

    LOG_INFO ( "Nonlinear solver residual", nls_->GetResidual ( ) );
    LOG_INFO ( "Nonlinear solver steps", nls_->iter ( ) );

    delete nls_;
    delete linsolver_;
}

void FlowTutorial::run_time_loop ( )
{
    // Visualize initial solution.
    visualize ( );
    const double end_time = params_["Instationary"]["Endtime"].get<double>( );

    LOG_INFO ( "timestep", "End time = " << end_time );
    LOG_INFO ( "timestep", "Step length = " << dt_ );

    // Crank-Nicolson time-stepping method. At the beginning of each
    // time-step, the solution from the previous time-step is stored
    // in prev_sol_, which is used in InstationaryFlowAssembler. The
    // variable ts_ is used to keep track of the current
    // time-step. The solution is visualized at the end of each
    // time-step, in order to be able to animate it in Paraview.

    while ( ts_ * dt_ <= end_time )
    {
        LOG_INFO ( "timestep", "Solving time step " << ts_ );

        solve ( );

        prev_sol_->CloneFrom ( *sol_ );

        ++ts_;

        LOG_INFO ( "timestep", "Visualizing solution at time " << ts_ * dt_
                   << " (time-step " << ts_ << ")" );
        visualize ( );
    }
}

void FlowTutorial::visualize ( )
{

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

    // Create visualization from post-processing vector.
    sol_->UpdateCouplings ( );
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 0 ), "u" );
#if DIMENSION >= 2
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 1 ), "v" );
#endif
#if DIMENION == 3
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 2 ), "w" );
#endif
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), DIMENSION ), "p" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );
}

void FlowTutorial::adapt ( )
{
    is_done_ = true;
}

void FlowTutorial::compute_residual ( const LAD::VectorType& u,
                                      LAD::VectorType* F )
{
    if ( solve_instationary_ )
    {
        compute_instationary_residual ( u, F );
    }
    else
    {
        compute_stationary_residual ( u, F );
    }

    // correct BC -- set Dirichlet dofs to 0
    if ( !dirichlet_dofs_.empty ( ) )
    {
        std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
        F->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                       vec2ptr ( zeros ) );
    }
}

void FlowTutorial::compute_stationary_residual ( const LAD::VectorType& u,
                                                 LAD::VectorType* F )
{
    StationaryFlowAssembler local_asm ( u, nu_, rho_ );
    global_asm_.assemble_vector ( space_, local_asm, *F );
}

void FlowTutorial::compute_instationary_residual ( const LAD::VectorType& u,
                                                   LAD::VectorType* F )
{
    InstationaryFlowAssembler local_asm ( nu_, rho_ );
    local_asm.set_newton_solution ( &u );
    local_asm.set_time_solution ( prev_sol_ );
    local_asm.set_time_stepping_weights ( alpha1_, alpha2_, alpha3_ );

    global_asm_.assemble_vector ( space_, local_asm, *F );
}

void FlowTutorial::compute_jacobian ( const LAD::VectorType& u,
                                      LAD::MatrixType* DF )
{
    if ( solve_instationary_ )
    {
        compute_instationary_matrix ( u, DF );
    }
    else
    {
        compute_stationary_matrix ( u, DF );
    }

    // correct BC -- set Dirichlet rows to identity
    if ( !dirichlet_dofs_.empty ( ) )
    {
        DF->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
    }
}

void FlowTutorial::compute_stationary_matrix ( const LAD::VectorType& u,
                                               LAD::MatrixType* DF )
{
    StationaryFlowAssembler local_asm ( u, nu_, rho_ );
    global_asm_.assemble_matrix ( space_, local_asm, *DF );
}

void FlowTutorial::compute_instationary_matrix ( const LAD::VectorType& u,
                                                 LAD::MatrixType* DF )
{
    InstationaryFlowAssembler local_asm ( nu_, rho_ );

    local_asm.set_newton_solution ( &u );
    local_asm.set_time_solution ( prev_sol_ );
    local_asm.set_time_stepping_weights ( alpha1_, alpha2_, alpha3_ );

    global_asm_.assemble_matrix ( space_, local_asm, *DF );
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

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            double& vol )
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

void FlowTutorial::ApplyFilter ( LAD::VectorType& u )
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

    LOG_INFO ( "pressure_filter", "Average pressure before filter = "
               << average_pressure );

    pressure_correction_.CloneFromWithoutContent ( u );
    pressure_correction_.Zeros ( );

    // set value for pressure dofs to average pressure
    std::vector<int> cell_p_dofs;
    std::vector<int> local_p_dofs;
    for ( EntityIterator it = mesh_->begin ( DIMENSION ), end = mesh_->end ( DIMENSION );
          it != end; ++it )
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
    std::fill ( p_correction_values.begin ( ), p_correction_values.end ( ),
                average_pressure );

    pressure_correction_.SetValues ( vec2ptr ( local_p_dofs ), local_p_dofs.size ( ),
                                     vec2ptr ( p_correction_values ) );

    pressure_correction_.WriteHDF5 ( "debug", "vectors", "pressure_correction" );
    u.WriteHDF5 ( "debug", "vectors", "sol_before" );
    u.Axpy ( pressure_correction_, -1. );

    u.UpdateCouplings ( );
    u.WriteHDF5 ( "debug", "vectors", "sol_after" );

    PressureIntegral int_p_check ( u );
    global_asm_.integrate_scalar ( space_, int_p_check, total_pressure );
    MPI_Allreduce ( &total_pressure, &recv, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    total_pressure = recv;
    LOG_INFO ( "pressure_filter", "Average pressure after filter = "
               << total_pressure / integrated_vol );
}
