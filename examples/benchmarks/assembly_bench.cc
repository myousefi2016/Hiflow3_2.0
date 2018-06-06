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

#include "assembly_bench.h"

// This program is meant to run on one process only.

// List of tests. TEST_FUNCS contains the number of tests defined, and
// TEST_NAMES contains the names of the tests. Modify these lists to
// add a new test. Also setup() has to be modified to set up the space correctly.
const int NUM_TESTS = 3;
typedef void (*TestFunc )( );

TestFunc TEST_FUNCS[NUM_TESTS] = {
                                  &test_laplace_aa,
                                  &test_navier_stokes_aa,
                                  &test_navier_stokes_aa_alt
};

const char* TEST_NAMES[NUM_TESTS] = {
                                     "Laplace AA",
                                     "Inst. Navier-Stokes AA",
                                     "Inst. Navier-Stokes Alt. AA"
};

// Name of current test
const char* g_test_name;

// Index of current test
int g_test_num;

// Laplace test for AssemblyAssistant

void test_laplace_aa ( )
{
    StandardGlobalAssembler<double> global_asm;
    Timer global_assembly;
    for ( g_run = 0; g_run < g_num_runs; ++g_run )
    {
        if ( g_dim == 2 )
        {

            LocalLaplaceAssembler<2> local_asm;

            global_assembly.reset ( );
            global_assembly.start ( );

            global_asm.assemble_matrix ( g_space, local_asm, g_matrix );
            global_assembly.stop ( );

        }
        else if ( g_dim == 3 )
        {

            LocalLaplaceAssembler<3> local_asm;

            global_assembly.reset ( );
            global_assembly.start ( );

            global_asm.assemble_matrix ( g_space, local_asm, g_matrix );
            global_assembly.stop ( );
        }
        g_global_times.at ( g_run ) = global_assembly.get_duration ( );
    }
}

namespace
{
    const double nu = 0.0001, rho = 1.0;
    std::vector<double> alphas ( 5, 0.5 );
}

// Instationary Navier Stokes test for AssemblyAssistant.

void test_navier_stokes_aa ( )
{
    std::vector<double> vec_values ( g_vec.size_global ( ), 3.14 );
    g_vec.SetValues ( vec2ptr ( vec_values ) );

    Timer global_assembly;
    StandardGlobalAssembler<double> global_asm;
    for ( g_run = 0; g_run < g_num_runs; ++g_run )
    {
        if ( g_dim == 2 )
        {

            InstationaryFlowAssembler<2> local_asm ( nu, rho );
            local_asm.set_newton_solution ( &g_vec );
            local_asm.set_prev_solution ( &g_vec );
            local_asm.set_timestep_parameters ( alphas );

            global_assembly.reset ( );
            global_assembly.start ( );

            global_asm.assemble_matrix ( g_space, local_asm, g_matrix );
            global_assembly.stop ( );

        }
        else if ( g_dim == 3 )
        {

            InstationaryFlowAssembler<3> local_asm ( nu, rho );
            local_asm.set_newton_solution ( &g_vec );
            local_asm.set_prev_solution ( &g_vec );
            local_asm.set_timestep_parameters ( alphas );

            global_assembly.reset ( );
            global_assembly.start ( );

            global_asm.assemble_matrix ( g_space, local_asm, g_matrix );
            global_assembly.stop ( );
        }
        g_global_times.at ( g_run ) = global_assembly.get_duration ( );
    }
}

// Alternative formulation of Instationary Navier-Stokes test for AssemblyAssistant.

void test_navier_stokes_aa_alt ( )
{
    std::vector<double> vec_values ( g_vec.size_global ( ), 3.14 );
    g_vec.SetValues ( vec2ptr ( vec_values ) );

    Timer global_assembly;
    StandardGlobalAssembler<double> global_asm;
    for ( g_run = 0; g_run < g_num_runs; ++g_run )
    {
        if ( g_dim == 2 )
        {

            AltInstationaryFlowAssembler<2> local_asm ( nu, rho );
            local_asm.set_newton_solution ( &g_vec );
            local_asm.set_time_solution ( &g_vec );
            local_asm.set_time_stepping_weights ( 0.5, 0.5, 0.5 );

            global_assembly.reset ( );
            global_assembly.start ( );

            global_asm.assemble_matrix ( g_space, local_asm, g_matrix );
            global_assembly.stop ( );

        }
        else if ( g_dim == 3 )
        {

            AltInstationaryFlowAssembler<3> local_asm ( nu, rho );
            local_asm.set_newton_solution ( &g_vec );
            local_asm.set_time_solution ( &g_vec );
            local_asm.set_time_stepping_weights ( 0.5, 0.5, 0.5 );

            global_assembly.reset ( );
            global_assembly.start ( );

            global_asm.assemble_matrix ( g_space, local_asm, g_matrix );
            global_assembly.stop ( );
        }
        g_global_times.at ( g_run ) = global_assembly.get_duration ( );
    }
}

// Common setup function, called before each test is run.

void setup ( )
{

    // Read in mesh only on first test.
    if ( g_test_num == g_first_test )
    {

        g_sys.Platform = APP_PLATFORM;
        init_platform ( g_sys );

        g_mesh = read_mesh_from_file ( g_mesh_filename, g_dim, g_dim, 0 );
        for ( int r = 0; r < g_refine_level; ++r )
        {
            g_mesh = g_mesh->refine ( );
        }
    }

    std::vector < std::vector<bool> > coupling_vars;

    // Setup space according to which test case is to be run.
    switch ( g_test_num )
    {
        case 0:
        {
            g_space.Init ( g_fe_degree, *g_mesh );

            coupling_vars.clear ( );
            coupling_vars.resize ( 1 );
            coupling_vars[0].push_back ( true );
            break;
        }
        case 1:
        case 2:
        {
            std::vector<int> deg ( g_dim + 1, g_fe_degree );
            deg[g_dim] = g_fe_degree - 1;
            g_space.Init ( deg, *g_mesh );

            coupling_vars.clear ( );
            coupling_vars.resize ( g_dim + 1 );
            for ( int i = 0; i < g_dim; ++i )
            {
                for ( int j = 0; j < g_dim + 1; ++j )
                {
                    coupling_vars[i].push_back ( true );
                }
            }
            for ( int i = 0; i < g_dim; ++i )
            {
                coupling_vars[g_dim].push_back ( true );
            }
            coupling_vars[g_dim].push_back ( false );
            break;
        }
        default:
            abort ( );
    };

    g_couplings.Init ( WORLD_COMM, g_space.dof ( ) );
    g_matrix.Init ( WORLD_COMM, g_couplings, g_sys.Platform, APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    g_vec.Init ( WORLD_COMM, g_couplings, g_sys.Platform, APP_LINALG_IMPLEMENTATION );

    InitStructure ( g_space, &g_diag_rows, &g_diag_cols, &g_off_diag_rows, &g_off_diag_cols, &coupling_vars );
    g_couplings.InitializeCouplings ( g_off_diag_rows, g_off_diag_cols );

    g_matrix.InitStructure ( vec2ptr ( g_diag_rows ), vec2ptr ( g_diag_cols ), g_diag_rows.size ( ),
                             vec2ptr ( g_off_diag_rows ), vec2ptr ( g_off_diag_cols ), g_off_diag_rows.size ( ) );
    g_matrix.Zeros ( );
    g_vec.InitStructure ( );
    g_vec.Zeros ( );

    g_global_times.clear ( );
    g_global_times.resize ( g_num_runs, 0. );
    g_local_times.clear ( );
    const int num_cells = g_mesh->num_entities ( g_dim );
    g_local_times.resize ( g_num_runs, std::vector<double>( num_cells, 0. ) );
    g_run = 0;
}

// Compute mean and variance of a vector of values.

void compute_stats ( const std::vector<double>& values, double& mean, double& variance )
{
    double sum = 0.;
    double sum_squared = 0.;
    const int num_elem = values.size ( );

    for ( int i = 0; i < num_elem; ++i )
    {
        const double t = values[i];
        sum += t;
        sum_squared += t * t;
    }

    mean = sum / num_elem;
    variance = 1. / ( num_elem - 1 ) * sum_squared - 1. / ( num_elem * ( num_elem - 1 ) ) * sum * sum;
}

// Print statistics from last test. Different verbosity levels can be used:
// Level 1: Global assembly, mean and std.dev only.
// Level 2: Global and local assemblies, mean and std.dev only.
// Level 3: Global and local assemblies, output for all runs as well as mean and std.dev.

void print_stats ( )
{
    const int num_cells = g_mesh->num_entities ( g_dim );
    std::cout << "==================================================\n";
    if ( g_verbosity >= 1 )
    {
        std::cout << "Dimension = " << g_dim << "\n";
        std::cout << "Num cells = " << num_cells << "\n";
        std::cout << "FE deg = " << g_fe_degree << "\n";

        std::cout << "\nGlobal assembly times\n";

        if ( g_verbosity >= 3 )
        {
            for ( int r = 0; r < g_num_runs; ++r )
            {
                const double t = g_global_times.at ( r );
                std::cout << r + 1 << "\t" << t << "\n";
            }
        }

        double global_mean_time;
        double global_var_time;

        compute_stats ( g_global_times, global_mean_time, global_var_time );

        const double global_std_dev_time = std::sqrt ( global_var_time );

        std::cout << "\nMean time = " << global_mean_time << " s.\n"
                << "Std.dev = " << global_std_dev_time << " s.\n";

        std::cout << "\n";
        if ( g_verbosity >= 2 )
        {
            std::vector<double> local_times_mean ( g_num_runs, 0. );
            std::vector<double> local_times_var ( g_num_runs, 0. );

            for ( int r = 0; r < g_num_runs; ++r )
            {
                compute_stats ( g_local_times[r], local_times_mean[r], local_times_var[r] );
            }

            std::cout << "\nMean Local Assembly times\n";

            if ( g_verbosity >= 3 )
            {

                for ( int r = 0; r < g_num_runs; ++r )
                {
                    const double t = local_times_mean.at ( r );
                    const double t_stddev = std::sqrt ( local_times_var.at ( r ) );
                    const double total_local_time =
                            std::accumulate ( g_local_times[r].begin ( ), g_local_times[r].end ( ), 0. );

                    std::cout << r + 1 << "\t" << t
                            << "\t\t" << t_stddev
                            << "\t\t" << total_local_time / g_global_times.at ( r ) * 100.0 << " %\n";
                }
            }

            double local_mean_of_means, local_var_of_means;
            compute_stats ( local_times_mean, local_mean_of_means, local_var_of_means );

            std::cout << "\nMean time = " << local_mean_of_means << " s.\n"
                    << "Std.dev = " << local_var_of_means << " s.\n";
        }
        std::cout << "==================================================\n";
        std::cout << "\n\n";
    }
}

// Write out the matrix computed in the last test.

void write_matrix ( )
{
    std::stringstream sstr;
    sstr << "matrix_test_" << g_test_num;
    g_matrix.diagonal ( ).WriteFile ( sstr.str ( ).c_str ( ) );
}

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    int num_partitions = -1;
    MPI_Comm_size ( MPI_COMM_WORLD, &num_partitions );

    if ( num_partitions != 1 )
    {
        if ( rank == 0 )
        {
            std::cerr << "assembly_bench only works in sequential mode: exiting!\n\n\n";
        }
        return 1;
    }

    {
        // redirect info log
        std::ofstream info_log ( "assembly_bench_info_log" );
        LogKeeper::get_log ( "info" ).set_target ( &info_log );

        // default parameter values
        g_mesh_filename = "";
        g_fe_degree = 1;
        g_refine_level = 0;
        g_dim = 2;
        g_num_runs = 10;
        g_verbosity = 1;
        g_first_test = 0;
        g_last_test = NUM_TESTS;
        g_write_matrix = false;

        // parse command line
        int c;
        while ( ( c = getopt ( argc, argv, "p:l:d:r:v:b:e:m" ) ) != -1 )
        {
            switch ( c )
            {
                case 'p':
                    g_fe_degree = atoi ( optarg );
                    break;
                case 'l':
                    g_refine_level = atoi ( optarg );
                    break;
                case 'd':
                    g_dim = atoi ( optarg );
                    break;
                case 'r':
                    g_num_runs = atoi ( optarg );
                    break;
                case 'v':
                    g_verbosity = atoi ( optarg );
                    break;
                case 'b':
                    g_first_test = std::max ( 1, atoi ( optarg ) ) - 1;
                    break;
                case 'e':
                    g_last_test = std::min ( NUM_TESTS, atoi ( optarg ) );
                    break;
                case 'm':
                    g_write_matrix = true;
                    break;
                case '?':
                    if ( optopt == 'p' || optopt == 'l' || optopt == 'd'
                         || optopt == 'r' || optopt == 'v' || optopt == 'b' || optopt == 'e' )
                    {
                        std::cerr << "Option " << char(optopt ) << " requires an argument.\n";
                    }
                    else if ( isprint ( optopt ) )
                    {
                        std::cerr << "Unknown option " << char(optopt ) << "\n";
                    }
                    else
                    {
                        std::cerr << "Unknown option character" << char(optopt ) << "\n";
                    }
                    return 1;
                default:
                    std::cerr << "c = " << c << ", optopt = " << optopt << "\n";
                    abort ( );
            }
        }

        if ( optind == argc )
        {
            std::cerr << "Missing mesh filename: exiting!\n";
            std::cerr << "For overview of options, see assembly_bench.h\n";
            LogKeeper::get_log ( "info" ).flush ( );
            return 1;
        }

        // read mesh filename
        g_mesh_filename = argv[optind];

        // Run all tests
        for ( g_test_num = g_first_test; g_test_num < g_last_test; ++g_test_num )
        {
            g_test_name = TEST_NAMES[g_test_num];
            std::cout << "Running test " << g_test_name << "...\n";
            setup ( );
            std::cout << "set...\n";
            ( TEST_FUNCS[g_test_num] )( );
            std::cout << "done!    \n";
            print_stats ( );

            if ( g_write_matrix )
            {
                write_matrix ( );
            }
        }
    }
    MPI_Finalize ( );
    return 0;
}
