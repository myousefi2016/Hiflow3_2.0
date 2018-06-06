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

/// \author Aksel Alpay, Martin Wlotzka

#include "geometric_multigrid.h"

static const char* PARAM_FILENAME = "geometric_multigrid.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class MultrigridTutorial
{
  public:

    MultrigridTutorial ( const std::string& param_filename, const std::string& path_mesh )
    : refinement_level_ ( 0 ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    path_mesh ( path_mesh )
    {
    }

    // Main algorithm

    void run ( );

    ~MultrigridTutorial ( )
    {
    }

  private:
    // Member functions

    // Read and distribute mesh.

    std::string path_mesh;
    void run_tests ( );

    // Parameter data read in from file.
    PropertyTree params_;

    // Local mesh and mesh on master process.
    MeshPtr master_mesh_;

    // Current refinement level.
    int refinement_level_;

}; // end class

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    int rk;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rk );

    INFO = ( rk == 0 );

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

        std::ofstream info_log;
        std::ofstream debug_log;
        std::ofstream error_log;

        if ( INFO )
        {
            info_log.open ( "info.log" );
            LogKeeper::get_log ( "info" ).set_target ( &info_log );
            debug_log.open ( "debug.log" );
            LogKeeper::get_log ( "debug" ).set_target ( &debug_log );
            error_log.open ( "error.log" );
            LogKeeper::get_log ( "error" ).set_target ( &error_log );
        }

        // Create application object and run it
        MultrigridTutorial app ( param_filename, path_mesh );
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

//////////////// MultigridTutorial implementation //////////////

void MultrigridTutorial::run ( )
{
    // Read in the mesh on the master process. The mesh is chosen according to the dimension of the problem.

    int rank = 0;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    int nproc = 0;
    MPI_Comm_size ( MPI_COMM_WORLD, &nproc );

    std::string mesh_name;
    int finest_level = params_["Multigrid"]["finest_level"].get<int>( 6 );
    int num_levels = params_["Multigrid"]["num_levels"].get<int>( 3 );
    int use_full_mg = params_["Multigrid"]["use_full_mg"].get<int>( 0 );
    char ctype = params_["Multigrid"]["cycle_type"].get<char>( 'W' );
    int use_matrix_operators = params_["Multigrid"]["use_matrix_operators"].get<int>( 0 );
    int grid_transfer_on_device = params_["Multigrid"]["grid_transfer_on_device"].get<int>( 0 );
    int use_sleep = params_["Multigrid"]["use_sleep"].get<int>( 0 );
    int use_gpus = params_["Multigrid"]["use_gpus"].get<int>( 1 );
    int solve_mode = params_["Multigrid"]["solve_mode"].get<int>( 3 );
    int global_iter = params_["Multigrid"]["global_iter"].get<int>( 3 );
    int block_iter = params_["Multigrid"]["block_iter"].get<int>( 1 );
    int inner_iter = params_["Multigrid"]["inner_iter"].get<int>( 5 );
    Scalar w = params_["Multigrid"]["damping"].get<Scalar>( 0.5 );
    int async = params_["Multigrid"]["async"].get<int>( 0 );
    int cuda_block_size = params_["Multigrid"]["cuda_block_size"].get<int>( 128 );
    int gmg_max_iter = params_["Multigrid"]["max_iter"].get<int>( 50 );
    Scalar gmg_abs_tol = params_["Multigrid"]["abs_tol"].get<Scalar>( 1.0e-12 );
    Scalar gmg_rel_tol = params_["Multigrid"]["rel_tol"].get<Scalar>( 1.0e-9 );
    Scalar gmg_div_tol = params_["Multigrid"]["div_tol"].get<Scalar>( 1.0e6 );

    std::vector<int> growth_factors;
    std::vector<int> smoother_types;

    for ( int i = 0; i < num_levels - 1; ++i )
    {
        std::stringstream str;
        str << "growth_fac" << i;
        int fac = params_["Multigrid"][str.str ( )].get<int>( 2 );
        growth_factors.push_back ( fac );
    }

    for ( int i = 0; i < num_levels - 1; ++i )
    {
        std::stringstream str;
        str << "smoother_type" << i;
        int type = params_["Multigrid"][str.str ( )].get<int>( 0 );
        smoother_types.push_back ( type );
    }

    int num_variables = params_["FEM"]["num_variables"].get<int>( 1 );
    std::vector<int> degrees;
    for ( int i = 0; i < num_variables; ++i )
    {
        std::stringstream str;
        str << "deg" << i;
        int deg = params_["FEM"][str.str ( )].get<int>( 1 );
        degrees.push_back ( deg );
    }

    if ( rank == MASTER_RANK )
    {
        switch ( DIMENSION )
        {
            case 1:
            {
                mesh_name =
                        params_["Mesh"]["Filename1"].get<std::string>( "unit_line.inp" );
                break;
            }
            case 2:
            {
                mesh_name =
                        params_["Mesh"]["Filename2"].get<std::string>( "unit_square.inp" );
                break;
            }
            case 3:
            {
                mesh_name =
                        params_["Mesh"]["Filename3"].get<std::string>( "unit_cube.inp" );
                break;
            }

            default: assert ( 0 );
        }

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
        int coarsest_level = finest_level - num_levels + 1;
        for ( int i = 0; i < coarsest_level; ++i )
        {
            master_mesh_ = master_mesh_->refine ( );
        }
    }

    gmg::Settings settings;
    settings.implementation = NAIVE;
    settings.platform = CPU;
    settings.matrix_format = CSR;

    std::string method = params_["LinearSolver"]["Method"].get<std::string>( "NoPreconditioning" );
    int max_iter = params_["LinearSolver"]["max_iter"].get<int>( 1000 );
    double abs_tol = params_["LinearSolver"]["abs_tol"].get<double>( 1.0e-13 );
    double rel_tol = params_["LinearSolver"]["rel_tol"].get<double>( 1.0e-10 );
    double div_tol = params_["LinearSolver"]["div_tol"].get<double>( 1.0e6 );

    LinearSolver<LAD>* solver;
    solver = new CG<LAD>( );
    solver->InitParameter ( method );
    solver->InitControl ( max_iter, abs_tol, rel_tol, div_tol );

    // Create communicator hierarchy
    //     gmg::communicator_hierarchy_generators::ConstantGrowthFactorGenerator
    //       generator(MPI_COMM_WORLD, growth_factor);
    gmg::communicator_hierarchy_generators::IndividualGrowthFactorGenerator
    generator ( MPI_COMM_WORLD, growth_factors );

    // Create Solver
    GeometricMultiGrid<LAD> multigrid_solver (
                                               MASTER_RANK,
                                               num_levels,
                                               master_mesh_,
                                               degrees,
                                               settings,
                                               &generator );

    // Assemble System
    LocalPoissonAssembler local_asm;
    multigrid_solver.assemble_system ( local_asm, local_asm );

    // Initialize boundary conditions
    DirichletZero dirichlet_zero;
    multigrid_solver.set_dirichlet_bc ( dirichlet_zero );

    multigrid_solver.use_restriction_matrix ( use_matrix_operators > 0 );
    multigrid_solver.use_interpolation_matrix ( use_matrix_operators > 0 );

    if ( use_matrix_operators > 0 )
    {
        multigrid_solver.build_restriction_matrix ( false );
        multigrid_solver.build_interpolation_matrix ( true );
    }

    if ( grid_transfer_on_device > 0 )
    {
        multigrid_solver.prepare_grid_transfer_on_device ( );
    }

    // use sleep ?
    multigrid_solver.use_sleep ( use_sleep > 0 );

    // Setup coarse grid solver
    multigrid_solver.setup_solver ( solver );

    // Setup smoothers on each level
    //gmg::SetupJacobiSmoother<LAD> setup_smoother(solve_mode, global_iter, 1, w, (async > 0));
    //gmg::SetupAsyncIterGPUSmoother<LAD> setup_smoother;
    gmg::SetupIndividualSmoother<LAD> setup_smoother ( smoother_types );
    setup_smoother.SetJacobiParameters ( solve_mode, global_iter, 1, w, ( async > 0 ) );
    setup_smoother.SetAsyncIterGPUParameters ( use_gpus, 1, solve_mode, global_iter, block_iter, inner_iter, w, cuda_block_size );
    multigrid_solver.for_each_but_coarsest_level ( setup_smoother );

    GeometricMultiGrid<LAD>::CycleType c_type;
    if ( ctype == 'V' )
        c_type = GeometricMultiGrid<LAD>::V_CYCLE;
    else if ( ctype == 'W' )
        c_type = GeometricMultiGrid<LAD>::W_CYCLE;
    else
    {
        std::cerr << "ERROR: invalid cycle type " << ctype << "\n";
        exit ( -1 );
    }

    multigrid_solver.set_mode ( use_full_mg > 0, c_type );
    multigrid_solver.InitControl ( gmg_max_iter, gmg_abs_tol, gmg_rel_tol, gmg_div_tol );
    // Solve system
    double t_start = MPI_Wtime ( );
    MPI_Barrier ( MPI_COMM_WORLD );

    multigrid_solver.Solve ( );

    MPI_Barrier ( MPI_COMM_WORLD );
    double t_stop = MPI_Wtime ( );

    double dt_local, dt;
    dt_local = t_stop - t_start;
    MPI_Allreduce ( &dt_local, &dt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

    if ( rank == 0 ) std::cout << "solve time = " << dt << "\n";

    // Visualize current state of the solution vector on each level
    //     gmg::visualize_multilevel_solutions(multigrid_solver.get_multilevel_hierarchy(), "solution");
    delete solver;
}
