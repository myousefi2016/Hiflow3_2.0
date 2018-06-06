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

#include <iostream>
#include <string>

#include "conv_diff.h"

/// Program main loop: setup MPI, read parameters and start the application

int main ( int argc, char** argv )
{
    // Check incompatible compiler options    
#ifdef RUN_MODE0
#    ifdef SCHUR_SOLVER
    std::cout << " RUN_MODE0 + SCHUR_SOLVER DOES NOT WORK " << std::endl;
    exit ( -1 );
#    endif
#    ifdef USE_HYPRE
    std::cout << " RUN_MODE0 + USE_HYPRE DOES NOT WORK " << std::endl;
    exit ( -1 );
#    endif
#endif
#ifndef USE_PXEST
    std::cout << " ConvDiff requires compilation with flag USE_PXEST (requires P4EST) " << std::endl;
    exit ( -1 );
#endif

    // initialize
    MPI_Init ( &argc, &argv );

    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    if ( rank == 0 )
    {
        std::cout << "COMPILER FLAGS" << std::endl;
#ifdef NAT_CONV
        // Neglect DEP force
        std::cout << " # NAT_CONV " << std::endl;
#endif        
#ifdef RUN_MODE0
        // Compute IC for primal problem
        std::cout << " # RUN_MODE0 " << std::endl;
#endif
#ifdef RUN_MODE1
        // Compute Primal solution
        std::cout << " # RUN_MODE1 " << std::endl;
#endif
#ifdef RUN_MODE2
        // Compute Adatpive solution
        std::cout << " # RUN_MODE2 " << std::endl;
#endif        
#ifdef RUN_MODE3 
        // PostProcessing
        std::cout << " # RUN_MODE3 " << std::endl;
#endif    
#ifdef RUN_MODE4 
        // Evaluate QoI
        std::cout << " # RUN_MODE4 " << std::endl;
#endif    
#ifdef RUN_MODE5 
        // Compute derivative of QoI w.r.t perturbation
        std::cout << " # RUN_MODE5 " << std::endl;
#endif            
#ifdef DEVEL
        std::cout << " # DEVEL " << std::endl;
#else
        std::cout << " # PRODUCTION " << std::endl;
#endif
#ifdef AUGMENT_PRESS
        std::cout << " # AUGMENT_PRESS " << std::endl;
#endif
#ifdef OPT_ASM0
        std::cout << " # OPT_ASM0 " << std::endl;
#endif
#ifdef OPT_ASM1
        std::cout << " # OPT_ASM1 " << std::endl;
#endif
#ifdef OPT_ASM2
        std::cout << " # OPT_ASM2 " << std::endl;
#endif        
#ifdef SCHUR_SOLVER 
        std::cout << " # SCHUR_SOLVER " << std::endl;
#endif        
#ifdef USE_HYPRE 
        std::cout << " # USE_HYPRE  " << std::endl;
#endif
#ifdef USE_PXEST 
        std::cout << " # USE_PXEST  " << std::endl;
#endif
#ifdef PARALLEL_PARTITIONING 
        std::cout << " # PARALLEL_PARTITIONING " << std::endl;
#endif
#ifdef ROTATING_FOR         
        std::cout << " # ROTATING_FOR " << std::endl;
#endif    

        std::cout << " GHOST_LAYER_WIDTH = " << GHOST_LAYER_WIDTH << std::endl;
    }

    std::string root_path = ".";
    if ( argc > 1 )
    {
        root_path = std::string ( argv[1] );
    }
    stringstream spath;
    string path;
    spath.str ( "" );

    stringstream base_spath;
    string base_path;
    base_spath.str ( "" );
#ifdef DEVEL    
    spath << root_path << "/in/conv_diff_devel.xml";
    base_spath << root_path << "/in/met_flow_convdiff_devel.xml";
#else
    spath << root_path << "/in/conv_diff.xml";
    base_spath << root_path << "/in/met_flow_convdiff.xml";
#endif
    path = spath.str ( );
    base_path = base_spath.str ( );

    if ( rank == 0 )
    {
        std::cout << "XML path " << path << std::endl;
        std::cout << "XML base path " << base_path << std::endl;
    }

    string log_path = root_path + "/log/convdiff.log";
    std::ofstream info_log ( log_path.c_str ( ) );

    LogKeeper::get_log ( "info" ).set_target ( &info_log );
    LogKeeper::get_log ( "debug" ).set_target ( &( std::cout ) );
    LogKeeper::get_log ( "error" ).set_target ( &( std::cout ) );

    bool resume_run = false;
    if ( argc > 2 )
    {
        std::string resume_flag ( argv[2] );
        if ( resume_flag == "1" )
        {
            resume_run = true;
        }
    }

    MetFlowConvDiffApp application ( root_path, path, base_path, resume_run );
    application.run ( );

    MPI_Finalize ( );

    return 0;
}
