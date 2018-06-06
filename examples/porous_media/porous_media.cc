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

/// \author Ulrich Blunck, Tobias Hahn, Julian Kraemer

#include "porous_media.h"

// if you work with extra makefile BASEDIR is defined
// static const char* DATADIR = BASEDIR "/meshes/";
static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;
const char* PARAM_FILENAME;

// Main Application class
// GPME = Generalized Porous Media Equation

class GPME
{
  public:

    GPME ( )
    : comm_ ( MPI_COMM_WORLD ),
    num_partitions_ ( -1 ),
    params_ ( PARAM_FILENAME, MASTER_RANK, MPI_COMM_WORLD ),
    refinement_level_ ( 0 )
    {
    }

    // Main algorithm, all other methods of the class GMPE are executed here

    virtual void run ( )
    {
        // string simul_name_ used as prefix to all output data
        simul_name_ = params_["Output"]["SolutionFilename"].get<std::string>( );

        // set the name of the log files
        std::ofstream info_log ( ( simul_name_ + "_info_log" ).c_str ( ) );
        LogKeeper::get_log ( "info" ).set_target ( &info_log );

        // output parameters for logging
        LOG_INFO ( "parameters", params_ );

        // The type of mesh is given in the xml-files
        geometry_ = params_["FlowModel"]["Type"].get<std::string>( );

        // Test wether a suitable type of mesh was stated.
        if ( geometry_ != std::string ( "Channel" ) )
        {
            if ( geometry_ != std::string ( "Column" ) )
            {
                throw UnexpectedParameterValue ( "Geometry.Type", geometry_ );
            }
        }

        // store rank of processor
        MPI_Comm_rank ( comm_, &rank_ );
        // store number of processes/partitions
        MPI_Comm_size ( comm_, &num_partitions_ );

        // the time that is needed to prepare the mesh is measured
        setbuf ( stdout, NULL );
        start_timer ( "Reading Mesh ...", start_t, ptm );

        // Read, refine and partition mesh.
        prepare_mesh ( );
        end_timer ( start_t, end_t, diff_t, ptm );
        LOG_INFO ( "Reading Mesh ...      ", diff_t << "s" );

        // Set up datastructures and read in some parameters.
        prepare ( );
        LOG_INFO ( "simulation", "Solving stationary problem" );
        std::cout << "Solving problem" << std::endl;
        info_log.flush ( );

        // Solve nonlinear problem
        solve ( );

        // Write visualization data of the solution in a file.
        visualize ( );
    }

    // member functions
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

    // Read, refine and partition mesh.
    void prepare_mesh ( );

    // Set up datastructures and read in some parameters.
    void prepare ( );

    // Set up boundary conditions
    void prepare_bc ( );

    // Solve nonlinear problem
    void solve ( );

    // Write visualization data of the solution in a file.
    void visualize ( );

    // functions to start and end time measurement
    const void start_timer ( const char* inp_, time_t &start_t, struct tm * ptm );
    const void end_timer ( time_t start_t, time_t end_t, time_t &diff_t,
                           struct tm * ptm );

    // Helper functions for nonlinear solver
    // updates the residual
    void compute_residual ( );

    // updates matrix_ with the jacobian matrix
    void compute_jacobian ( );

    // Functions used to evaluate the flow-profile
    // at x = 4
    void compute_flowprofile ( );
    // at x = 3
    void compute_flowprofile2 ( );
    // Finds the entity "inside" which includes the point "coord"
    const bool find_entity ( std::vector<Coordinate>& coord, Entity& inside );
    // Evaluate the solution-function at the point "coord" within entity "cell"
    const double eval_sol ( int var, Entity& cell,
                            std::vector<Coordinate>& coord, CVector& sol );

    // member variables
    // Variable used to weight the continuum equation against the GPME
    double conti_weight;

    // MPI variables
    const MPI_Comm comm_;
    int rank_, num_partitions_;

    // Parameter data read in from file.
    PropertyTree params_;

    // parameter 'OutputPrefix': prefix for output files
    std::string simul_name_;
    // parameter 'Geometry.Type': "Channel"
    std::string geometry_;

    // Geometry variables
    double Um_, D_, Href_, rho_, mu_, U_max_;

    // Variables defining the porous media
    double eps_free, eps_por, kap_free, kap_por;

    // Iterational variables
    int iter_;

    // Mesh variables
    MeshPtr master_mesh_, mesh_;
    int refinement_level_;

    // vectorspace variables
    VectorSpace<double> space_;

    // Linear algebra variables
    Couplings<double> couplings_;
    CoupledMatrix<Scalar>* matrix_;
    CoupledVector<Scalar> *sol_, *prev_sol_, *cor_, *res_, *pressure_correction_;

    // assembling variables
    StandardGlobalAssembler<double> global_asm_;

    // boundary condition variables
    std::vector<int> dirichlet_dofs_;
    std::vector<Scalar> dirichlet_values_;

    // Variables for the evaluation of the time required for each computational
    // step
    struct tm * ptm;
    time_t start_t, end_t, diff_t;
};

////////////////////////////////MAIN/////////////////////////////////
// The main function starts the program by creating an instance of GMPE and
// running the method run on it.

int main ( int argc, char** argv )
{
    // initialize MPI
    MPI_Init ( &argc, &argv );

    // choose parameterfile depending on dimension
    if ( DIMENSION == 2 )
        PARAM_FILENAME = "porous_media.xml";
    else
        PARAM_FILENAME = "porous_media3d.xml";
    try
    {
        // run application
        GPME app;
        app.run ( );
    }
    catch ( std::exception& e )
    {
        std::cerr << "\nProgram ended with uncaught exception.\n";
        std::cerr << e.what ( ) << "\n";
        return -1;
    }

    // finalize MPI
    MPI_Finalize ( );

    return 0;
}
//////////////////////////////////////////////////////////////////////

// read, refine and distribute mesh

void GPME::prepare_mesh ( )
{
    // The following needs to be done only by one process, therefore it is only
    // done by the process whose rank = MASTER_RANK
    if ( rank ( ) == MASTER_RANK )
    {
        // The name of the mesh is given in the xml-file. The program reads it in
        // and uses this mesh
        const std::string mesh_name = params_["Mesh"][geometry_].get<std::string>( );
        std::string mesh_filename = std::string ( DATADIR ) + mesh_name;
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0 );

        // The initial mesh can be refined if it is too coarse in the beginning. If
        // this is desired, the variable RefLevel can be set to the desired level
        // in the xml-file
        refinement_level_ = 0;
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );

        // refine mesh globally until the desired refinement level RefLevel is
        // reached
        for ( int r = 0; r < initial_ref_lvl; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
            ++refinement_level_;
        }
        LOG_INFO ( "mesh", "Refinement level = " << refinement_level_ );
    }

    // If the computation is executed parallel, each process computes the solution
    // on one part of the mesh. The distribution is done via the following
    // function call
    MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK,
                                                    comm_ );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;

    // Ghost cells occur on the boundaries between two cell partitions. Since in
    // the finite element method cells depend on their neighbor-cells, some need
    // to be computed by two processes, but belong only to one. The other process
    // computes on a so called ghost cell.
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );

    // output on console
    std::ostringstream rank_str;
    rank_str << rank ( );

    // write out mesh data
    PVtkWriter writer ( comm_ );
    std::string output_file = std::string ( simul_name_ + "_mesh_local.pvtu" );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );

    // create boundary mesh and write out boundary mesh data
    MeshPtr bdy_mesh = MeshPtr ( mesh_->extract_boundary_mesh ( ) );
    VtkWriter bdy_writer;
    bdy_writer.add_attribute ( "_mesh_facet_index_", DIMENSION - 1 );
    std::stringstream bdy_sstr;
    bdy_sstr << simul_name_ << "_bdy_sd_" << rank_ << ".vtu";
    bdy_writer.write ( bdy_sstr.str ( ).c_str ( ), *bdy_mesh );
}

// Set up datastructures and read in some parameters.

void GPME::prepare ( )
{
    // prepare model parameters
    rho_ = params_["FlowModel"]["Density"].get<double>( );
    mu_ = params_["FlowModel"]["Viscosity"].get<double>( );
    conti_weight = params_["FlowModel"]["ContiWeight"].get<double>( );

    // The parameter u_max is read in and the avarage flow speed um_ is computed
    // depending on the setting. R_in is the radius of the inflow tube and R_out
    // is the radius of the outflow tube.
    // 2D-Channel: um_=umax * 2/3
    // 2D-Column:  um_=umax * 2/3 * R_in/R_out = 2/3*0.3= 1/5
    // 3D-Channel: um_=umax * 1/2
    // 3D-Column:  um_=umax * 1/2 * (R_in/R_out)^2 = 1/2*0.3*0.3 = 9./200
    if ( DIMENSION == 2 )
    {
        if ( geometry_ == std::string ( "Channel" ) )
        {
            Um_ = 2 * params_["FlowModel"]["InflowSpeed"].get<double>( ) / 3;
            Href_ = params_["FlowModel"]["HeightRef"].get<double>( );
            D_ = params_["FlowModel"]["InflowDiameter"].get<double>( );
        }
        else
        { // Geometry Column
            Um_ = params_["FlowModel"]["InflowSpeed"].get<double>( ) / 5;
            Href_ = params_["FlowModel"]["HeightRef"].get<double>( );
            // If this geometry is used, the inflow diameter is only
            // 0.3 * the diameter of the tube so this has to be corrected
            D_ = .3 * params_["FlowModel"]["InflowDiameter"].get<double>( );
        }
    }
    else
    {
        if ( geometry_ == std::string ( "Channel" ) )
        {
            Um_ = params_["FlowModel"]["InflowSpeed"].get<double>( ) / 2;
            Href_ = params_["FlowModel"]["HeightRef"].get<double>( );
            D_ = params_["FlowModel"]["InflowDiameter"].get<double>( );
        }
        else
        { // Geometry Column
            Um_ = .3 * .3 * params_["FlowModel"]["InflowSpeed"].get<double>( ) / 2;
            Href_ = params_["FlowModel"]["HeightRef"].get<double>( );
            // The inflow diameter is only 0.3 * the diameter of the tube if
            // this geometry is used so this has to be corrected
            D_ = .3 * params_["FlowModel"]["InflowDiameter"].get<double>( );
        }
    }

    start_timer ( "Preparing Space ...", start_t, ptm );

    // prepare space, therefore read in polynomial degree of finite elements of
    // different variables
    std::vector< int > degrees ( DIMENSION + 1 );
    const int u_deg = params_["FiniteElements"]["VelocityDegree"].get<int>( );
    const int p_deg = params_["FiniteElements"]["PressureDegree"].get<int>( );
    for ( int c = 0; c < DIMENSION; ++c )
    {
        degrees.at ( c ) = u_deg;
    }
    degrees.at ( DIMENSION ) = p_deg;

    // initialize finite element space
    space_.Init ( degrees, *mesh_ );

    // compute matrix graph
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    // prepare linear algebra structures
    couplings_.Clear ( );
    couplings_.Init ( communicator ( ), space_.dof ( ) );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    // initialization of needed global matrix and vectors and setting them to 0

    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
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
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    sol_->Init ( comm_, couplings_ );
    sol_->InitStructure ( );
    sol_->Zeros ( );

    prev_sol_ = CoupVecFact.Get (
                                  params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    prev_sol_->Init ( comm_, couplings_ );
    prev_sol_->InitStructure ( );
    prev_sol_->Zeros ( );

    cor_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    cor_->Init ( comm_, couplings_ );
    cor_->InitStructure ( );
    cor_->Zeros ( );

    res_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    res_->Init ( comm_, couplings_ );
    res_->InitStructure ( );
    res_->Zeros ( );

    end_timer ( start_t, end_t, diff_t, ptm );
    LOG_INFO ( "Preparing Space ...      ", diff_t << "s" );

    start_timer ( "Preparing BC ...", start_t, ptm );
    // prepare dirichlet BC
    prepare_bc ( );
    end_timer ( start_t, end_t, diff_t, ptm );
    LOG_INFO ( "Preparing BC ...         ", diff_t << "s" );

    // Assign Values for Epsilon and Kappa
    // eps_free is the porosity where no porous media is found
    // eps_por is the porosity of the porous media
    // kap_free is the permeability where no porous media is found
    // kap_por is the permeability of the porous media
    eps_free = params_["FlowModel"]["PorosityFree"].get<double>( );
    eps_por = params_["FlowModel"]["PorosityPor"].get<double>( );
    kap_por = params_["FlowModel"]["PermeabilityPor"].get<double>( );
    kap_free = params_["FlowModel"]["PermeabilityFree"].get<double>( );
}

// Set up boundary conditions

void GPME::prepare_bc ( )
{
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    // compute Dirichlet values, distinguish between different geometries and
    // dimensions. The material numbers of the facets which lie on the boundary
    // of the mesh determine wether this facets is an inflow boundary, and outflow
    // boundary or a dirichlet boundary (i.e. a wall).
    if ( geometry_ == std::string ( "Channel" ) )
    {
        // read in material numbers of inflow and outflow boundary
        const int inflow_bdy = params_["Boundary"]["InflowMaterial"].get<int>( );
        const int outflow_bdy = params_["Boundary"]["OutflowMaterial"].get<int>( );
        U_max_ = params_["FlowModel"]["InflowSpeed"].get<double>( );

        if ( DIMENSION == 2 )
        {
            ChannelFlowBC bc[2] = {
                                   ChannelFlowBC ( 0, D_, U_max_, inflow_bdy, outflow_bdy ),
                                   ChannelFlowBC ( 1, D_, U_max_, inflow_bdy, outflow_bdy )
            };

            for ( int var = 0; var < DIMENSION; ++var )
            {
                compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                                    dirichlet_dofs_, dirichlet_values_ );
            }
        }
        else
        {
            assert ( DIMENSION == 3 );
            ChannelFlowBC3d bc[3] = {
                                     ChannelFlowBC3d ( 0, D_, U_max_, inflow_bdy, outflow_bdy ),
                                     ChannelFlowBC3d ( 1, D_, U_max_, inflow_bdy, outflow_bdy ),
                                     ChannelFlowBC3d ( 2, D_, U_max_, inflow_bdy, outflow_bdy )
            };
            for ( int var = 0; var < DIMENSION; ++var )
            {
                compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                                    dirichlet_dofs_, dirichlet_values_ );
            }
        }
    }
    else if ( geometry_ == std::string ( "Column" ) )
    {
        // read in material numbers of inflow and outflow boundary
        const int inflow_bdy = params_["Boundary"]["InflowMaterial"].get<int>( );
        const int outflow_bdy = params_["Boundary"]["OutflowMaterial"].get<int>( );
        U_max_ = params_["FlowModel"]["InflowSpeed"].get<double>( );

        if ( DIMENSION == 2 )
        {
            ChannelFlowBC bc[2] = {
                                   ChannelFlowBC ( 0, D_, U_max_, inflow_bdy, outflow_bdy ),
                                   ChannelFlowBC ( 1, D_, U_max_, inflow_bdy, outflow_bdy )
            };

            for ( int var = 0; var < DIMENSION; ++var )
            {
                compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                                    dirichlet_dofs_, dirichlet_values_ );
            }
        }
        else
        {
            assert ( DIMENSION == 3 );
            ChannelFlowBC3d bc[3] = {
                                     ChannelFlowBC3d ( 0, D_, U_max_, inflow_bdy, outflow_bdy ),
                                     ChannelFlowBC3d ( 1, D_, U_max_, inflow_bdy, outflow_bdy ),
                                     ChannelFlowBC3d ( 2, D_, U_max_, inflow_bdy, outflow_bdy )
            };
            for ( int var = 0; var < DIMENSION; ++var )
            {
                compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                                    dirichlet_dofs_, dirichlet_values_ );
            }
        }
    }
    else
    {
        assert ( false );
    }

    // apply boundary conditions to initial solution
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet boundary conditions
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
}

// Solve nonlinear problem

void GPME::solve ( )
{
    // read in nonlinear solver parameters
    const int nls_max_iter =
            params_["NonlinearSolver"]["MaximumIterations"].get<int>( );
    const double nls_abs_tol =
            params_["NonlinearSolver"]["AbsoluteTolerance"].get<double>( );
    const double nls_rel_tol =
            params_["NonlinearSolver"]["RelativeTolerance"].get<double>( );
    /*const double nls_div_tol =
                   params_["NonlinearSolver"]["DivergenceLimit"].get<double>();*/

    // for clarification write out the given parameters.
    std::cout << "           Nondimensional scales:" << std::endl;
    std::cout << "           U_ref      [m/s]  :  " << Um_ << std::endl;
    std::cout << "           L_ref      [m]    :  " << Href_ << std::endl;
    std::cout << "           P_ref      [Pa]   :  " << rho_ * Um_ * Um_ << std::endl;
    std::cout << "           Re         [-]    :  " << ( Um_ * Href_ * rho_ / mu_ )
            << std::endl;
    std::cout << "           Da         [-]    :  " << params_["FlowModel"]
            ["PermeabilityPor"].get<double>( ) / ( Href_ * Href_ ) << std::endl;
    LOG_INFO ( "Reference Velocity        ", Um_ << "m/s" );
    LOG_INFO ( "Reference Length          ", Href_ << "m" );
    LOG_INFO ( "Reference Pressure        ", rho_ * Um_ * Um_ << "N/s^2" );
    LOG_INFO ( "Reynolds Number           ", ( Um_ * Href_ * rho_ / mu_ ) );
    LOG_INFO ( "Darcy Number              ", params_["FlowModel"]
               ["PermeabilityPor"].get<double>( ) / ( Href_ * Href_ ) );

    start_timer ( "Setting up GMRES ...", start_t, ptm );

    // read in linear solver parameters
    GMRES<LAD> gmres;
    const int lin_max_iter =
            params_["LinearSolver"]["MaximumIterations"].get<int>( );
    const double lin_abs_tol =
            params_["LinearSolver"]["AbsoluteTolerance"].get<double>( );
    const double lin_rel_tol =
            params_["LinearSolver"]["RelativeTolerance"].get<double>( );
    const double lin_div_tol =
            params_["LinearSolver"]["DivergenceLimit"].get<double>( );
    const int basis_size = params_["LinearSolver"]["BasisSize"].get<int>( );

    // set up linear solver
    gmres.InitControl ( lin_max_iter, lin_abs_tol, lin_rel_tol, lin_div_tol );
    gmres.InitParameter ( basis_size, "NoPreconditioning" );
    gmres.SetupOperator ( *matrix_ );

    end_timer ( start_t, end_t, diff_t, ptm );
    LOG_INFO ( "Setting up GMRES ...      ", diff_t << "s" );

    std::ofstream iteration ( ( "Data/" + simul_name_ + "_iteration" ).c_str ( ) );

    // Write visualization data of the solution in a file.
    visualize ( );

    start_timer ( "Computing Residual ...", start_t, ptm );

    // Newton's method is used to solve the nonlinear problem.  The
    // functions compute_residual() and compute_jacobian()
    // update the variables matrix_ and res_, respectively.
    // The vector cor_ is set up to be used for the correction, and
    // the solution state is stored in sol_ .
    iter_ = 0;
    compute_residual ( );

    end_timer ( start_t, end_t, diff_t, ptm );
    LOG_INFO ( "Computing Residual ...      ", diff_t << "s" );

    // counter for Newton iterations
    int iter = 0;

    // variables for storing the average time it takes to compute the residual,
    // flow and Jacobian and to visualize and solve
    double avgres ( 0. ), avgflow ( 0. ), avgvisu ( 0. ), avgsolv ( 0. ), avgjaco ( 0. );

    // compute start residual
    const double initial_res_norm = res_->Norm2 ( );
    LOG_INFO ( "nonlinear", "Nonlinear solver starts with residual norm "
               << initial_res_norm );
    std::cout << "Nonlinear solver starts with residual norm "
            << initial_res_norm << std::endl;
    double res_norm = initial_res_norm;
    double prev_res_norm = initial_res_norm;
    double norm_ratio ( 1. );
    double step_;

    // Newton method
    while ( iter < nls_max_iter
            && res_norm > nls_abs_tol
            && res_norm > nls_rel_tol * initial_res_norm )
    {
        // Solve DF * cor = res
        iter_ = iter;
        start_timer ( "Computing Jacobian ...", start_t, ptm );

        // updates matrix_ with the jacobian matrix
        compute_jacobian ( );

        end_timer ( start_t, end_t, diff_t, ptm );

        // the time it took to compute the Jacobian in this step (diff_t) is added
        // to the variable avgjaco to compute the average time
        avgjaco += diff_t;
        LOG_INFO ( "Computing Jacobian ...      ", diff_t << "s    ("
                   << avgjaco / ( iter + 1 ) << ")" );

        // Compute correction cor
        cor_->Zeros ( );
        std::cout << "[ Iteration " << iter
                << " ]:  Solving linear problem with residual: "
                << res_norm << " ." << std::endl;
        iteration << iter << ",  " << res_norm << ",  " << norm_ratio << std::endl;
        start_timer ( "Solving ...", start_t, ptm );
        gmres.Solve ( *res_, cor_ );
        end_timer ( start_t, end_t, diff_t, ptm );
        avgsolv += diff_t;
        LOG_INFO ( "Solving ...                 ", diff_t << "s    ("
                   << avgsolv / ( iter + 1 ) << ")" );

        cor_->UpdateCouplings ( );

        // update sol = sol - cor
        step_ = 1.;
        std::cout << "Ratio of = " << norm_ratio << ", resulting Steplength = "
                << step_ << std::endl;
        sol_->Axpy ( *cor_, -step_ ); /// damped Newton method
        sol_->UpdateCouplings ( );

        start_timer ( "Computing Flow Profile ...", start_t, ptm );

        // For computing the flow profile at x = 4 use compute_flowprofile()
        // For computing the flow profile at x = 3 use compute_flowprofile2()
        compute_flowprofile ( );
        // compute_flowprofile2();
        end_timer ( start_t, end_t, diff_t, ptm );

        // the time it took to compute the flowprofile in this step (diff_t) is
        // added to the variable avgflow to compute the average time
        avgflow += diff_t;
        LOG_INFO ( "Computing Flow Profile ...  ", diff_t << "s    ("
                   << avgflow / ( iter + 1 ) << ")" );

        start_timer ( "Visualizing ...", start_t, ptm );

        // visualize solution of each Newton step
        visualize ( );

        end_timer ( start_t, end_t, diff_t, ptm );

        // the time it took to viualize the solution in this step (diff_t) is added
        // to the variable avgvisu to compute the average time
        avgvisu += diff_t;
        LOG_INFO ( "Visualizing ...             ", diff_t << "s    (" << avgvisu / ( iter + 1 ) << ")" );

        // Compute new residual

        start_timer ( "Computing new residual ...", start_t, ptm );

        compute_residual ( );

        prev_res_norm = res_norm;
        res_norm = res_->Norm2 ( );
        norm_ratio = res_norm / prev_res_norm;

        end_timer ( start_t, end_t, diff_t, ptm );
        avgres += diff_t;
        LOG_INFO ( "Computing new Residual ...  ", diff_t << "s    ("
                   << avgres / ( iter + 1 ) << ")" );
        ++iter;
    }

    std::cout << "Nonlinear Solver Residual " << res_norm << std::endl;
    std::cout << "Nonlinear Solver Steps " << iter << std::endl;
    LOG_INFO ( "Nonlinear solver residual", res_norm );
    LOG_INFO ( "Nonlinear solver steps", iter );
}

// Write visualization data of the solution in a file.

void GPME::visualize ( )
{
    // initialize visualization
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    // setup output file name
    std::stringstream input;
    input << simul_name_ << "_solution_tutorial";
    input << "_stationary";

    std::vector<double> material_number ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

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
        remote_index.at ( it->index ( ) ) = temp1;
        sub_domain.at ( it->index ( ) ) = temp2;
        material_number.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
    }

    sol_->UpdateCouplings ( );

    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 0 ), "u" );
#if(DIMENSION >= 2)
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 1 ), "v" );
#endif
#if(DIMENSION == 3)
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), 2 ), "w" );
#endif
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ), DIMENSION ), "p" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );

    visu.write ( input.str ( ) );
}

void GPME::compute_residual ( )
{
    PorousMediaAssembler local_asm ( *sol_, mu_, rho_, Um_, Href_, conti_weight,
                                     eps_free, eps_por, kap_free, kap_por,
                                     geometry_ );

    global_asm_.assemble_vector ( space_, local_asm, *res_ );

    // correct boundary conditions -- set Dirichlet dofs to 0
    if ( !dirichlet_dofs_.empty ( ) )
    {
        std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
        res_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( zeros ) );
    }

    res_->UpdateCouplings ( );
}

void GPME::compute_jacobian ( )
{
    PorousMediaAssembler local_asm ( *sol_, mu_, rho_, Um_, Href_, conti_weight,
                                     eps_free, eps_por, kap_free, kap_por, geometry_ );

    // assemble system matrix
    global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );

    // correct BC -- set Dirichlet rows to identity
    if ( !dirichlet_dofs_.empty ( ) )
    {
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
    }
}

// function to start time measurement

const void GPME::start_timer ( const char* inp_, time_t &start_t,
                               struct tm * ptm )
{
#ifdef debug
    time ( &start_t );
    ptm = gmtime ( &start_t );
    std::printf ( "%02d:%02d:%02d : %-30s",
                  ( ptm->tm_hour + 1 ) % 24, ptm->tm_min, ptm->tm_sec, inp_ );
#endif
}

// function to end time measurement

const void GPME::end_timer ( time_t start_t, time_t end_t, time_t &diff_t,
                             struct tm * ptm )
{
#ifdef debug
    time ( &end_t );
    diff_t = end_t - start_t;
    ptm = gmtime ( &diff_t );
    std::printf ( ": Done   ( %02d:%02d:%02d )\n", ( ptm->tm_hour ) % 24,
                  ptm->tm_min, ptm->tm_sec );
#endif
}

// identifies the entity which contains the point with the coordinates coord.
// It is stored in inside

const bool GPME::find_entity ( std::vector<Coordinate>& coord, Entity& inside )
{
    EntityIterator cell_it = mesh_->begin ( DIMENSION );

    std::vector<Coordinate> test_coord;
    bool found = false;
    std::vector<Coordinate> coord_ref ( coord.size ( ) );

    while ( !found && cell_it != mesh_->end ( DIMENSION ) )
    {
        cell_it->get_coordinates ( test_coord );

        if ( ( std::abs ( test_coord[0] - coord[0] ) <= 0.2 ) and
             // avert bad inversions
             (std::abs ( test_coord[1] - coord[1] ) <= 0.2 ) )
        {
            if ( DIMENSION == 2 )
            {
                space_.GetCellTransformation ( *cell_it ).inverse ( coord[0], coord[1],
                                                                    coord_ref[0], coord_ref[1] );
                if ( ( 0 <= coord_ref[0] ) and (0 <= coord_ref[1] ) )
                    if ( ( coord_ref[0] + coord_ref[1] ) <= 1 )
                    { // specific to triangles -
                        // needs to be adapted for other entities
                        inside = *cell_it;
                        found = true;
                    }
            }
        }

        cell_it++;
    }
    return found;
}

// computes the flowprofile at x = 4.

void GPME::compute_flowprofile ( )
{
#if DIMENSION == 2  // Currently only for 2-dimensional cases
    std::vector<Coordinate> coord;
    Entity inside;
    int n_eval = 100; // Number of evaluations-2
    std::vector<double> velocity_values ( n_eval );
    std::vector<double> velocity_valuesy ( n_eval );
    std::vector<double> pressure_values ( n_eval );
    double y_step = 1. / ( 2 * n_eval - 2 );
    coord.push_back ( 4 );
    coord.push_back ( 0 );

    // Set the directory for the computed flow profiles
    std::ofstream profile ( ( simul_name_ + "_flow_profile" ).c_str ( ) );

    for ( int i = 0; i < ( n_eval ); ++i )
    {
        double y_step2 = y_step*i;
        coord[1] = y_step2;
        if ( find_entity ( coord, inside ) )
        {
            velocity_values[i] = eval_sol ( 0, inside, coord, *sol_ );
            velocity_valuesy[i] = eval_sol ( 1, inside, coord, *sol_ );
            pressure_values[i] = eval_sol ( 2, inside, coord, *sol_ );
            profile << y_step2 << ",  " << velocity_values[i] << ",  "
                    << velocity_valuesy[i] << ",  " << pressure_values[i] << std::endl;
        }
    }
#endif
}

// computes the flow profile at x = 3

void GPME::compute_flowprofile2 ( )
{
#if DIMENSION == 2
    std::vector<Coordinate> coord;
    Entity inside;
    int n_eval = 100; // number of substeps
    std::vector<double> velocity_values ( n_eval );
    std::vector<double> velocity_valuesy ( n_eval );
    std::vector<double> pressure_values ( n_eval );
    double y_step = 1. / ( 2 * n_eval - 2 );
    coord.push_back ( 3 );
    coord.push_back ( 0 );

    // Set the directory of the computed flow profile
    std::ofstream profile ( ( simul_name_ + "_flow_profile2" ).c_str ( ) );

    for ( int i = 0; i < ( n_eval ); ++i )
    {
        double y_step2 = y_step*i;
        coord[1] = y_step2;
        if ( find_entity ( coord, inside ) )
        {
            velocity_values[i] = eval_sol ( 0, inside, coord, *sol_ );
            velocity_valuesy[i] = eval_sol ( 1, inside, coord, *sol_ );
            pressure_values[i] = eval_sol ( 2, inside, coord, *sol_ );
            profile << y_step2 << ",  " << velocity_values[i] << ",  "
                    << velocity_valuesy[i] << ",  " << pressure_values[i] << std::endl;
        }
    }
#endif
}

// used for the computation of flowprofiles. Evaluates the solution at the given
// coordinates coord for the variable var

const double GPME::eval_sol ( int var, Entity& cell,
                              std::vector<Coordinate>& coord, CVector& sol )
{
    // summation of
    // normal_vector *  ((P(left)+P(right)/2)*l)   /(density * PI * Radius^2)
    double sum = 0.0;
    int testvar = var;
    Entity& testcell = cell;
    std::vector<Coordinate> test_coord = coord;
    std::vector<int> ind ( sol_->size_global ( ) );
    std::vector<Scalar> visu_vec ( sol_->size_global ( ) );
    for ( int i = 0; i < static_cast < int > ( ind.size ( ) ); ++i )
    {
        ind[i] = i;
    }

    sol_->GetValues ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( visu_vec ) );

    //  Global DoF Ids on the given mesh cell
    std::vector<int> global_dof_ids;
    space_.GetDofIndices ( testvar, testcell, &global_dof_ids );

    //  Determine corresponding coordinates on reference cell via transformation
    std::vector<Coordinate> test_coord_ref ( coord.size ( ) );

    if ( space_.get_dim ( ) == 2 )
        space_.GetCellTransformation ( testcell ).inverse ( test_coord[0],
                                                            test_coord[1], test_coord_ref[0], test_coord_ref[1] );
    else
        space_.GetCellTransformation ( testcell ).inverse ( test_coord[0],
                                                            test_coord[1], test_coord[2], test_coord_ref[0],
                                                            test_coord_ref[1], test_coord_ref[2] );

    std::vector<double> weights ( global_dof_ids.size ( ) );

    space_.fe_manager ( ).get_fe_on_cell ( testcell.index ( ), testvar )
            ->N ( test_coord_ref, weights );

    //  Summation over weights of shapefunctions on reference cell
    for ( int i_loc = 0; i_loc < static_cast < int > ( global_dof_ids.size ( ) ); ++i_loc )
        sum += visu_vec[global_dof_ids[i_loc]] * weights[i_loc];

    return sum;
}
