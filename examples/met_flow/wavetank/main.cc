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

#include "wavetank.h"
#undef SCHUR_SOLVER
#undef ROTATING_FOR
#undef AUGMENT_PRESS
#undef USE_HYPRE
const bool CREATE_DEFAULT_OUTPUT_DIRECTORIES = true;

/// Program main loop: setup MPI, read parameters and start the application

int main ( int argc, char** argv )
{

    // initialize
    MPI_Init ( &argc, &argv );

    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    if ( rank == 0 )
    {
        std::cout << "COMPILER FLAGS" << std::endl;
#ifdef AUGMENT_PRESS
        std::cout << " # AUGMENT_PRESS " << std::endl;
#endif
#ifdef SCHUR_SOLVER
        std::cout << " Schur solver does not work with wavetank example" << std::endl;
        exit ( -1 );
#endif
#ifdef USE_HYPRE
        std::cout << " Hypre does not work with wavetank example" << std::endl;
        exit ( -1 );
#endif
#ifdef ZERO_TOP 
        std::cout << " # ZERO_TOP " << std::endl;
#endif
#ifdef ROTATING_FOR   
        std::cout << " # ROTATING_FOR " << std::endl;
#endif
#ifdef APPDIM
        std::cout << " # APPDIM = " << APPDIM << std::endl;
#endif
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
    spath << root_path << "/in/wavetank_devel.xml";
    base_spath << root_path << "/in/met_flow_wavetank_devel.xml";
#else
    spath << root_path << "/in/wavetank.xml";
    base_spath << root_path << "/in/met_flow_wavetank.xml";
#endif
    path = spath.str ( );
    base_path = base_spath.str ( );

    if ( rank == 0 )
    {
        std::cout << "XML path " << path << std::endl;
        std::cout << "XML base path " << base_path << std::endl;
    }

    bool resume_run = false;
    if ( argc > 2 )
    {
        std::string resume_flag ( argv[2] );
        if ( resume_flag == "1" )
        {
            resume_run = true;
        }
    }

    std::string info_file = root_path + "/log/wavetank_info.log";
    std::ofstream info_log ( info_file.c_str ( ) );
    LogKeeper::get_log ( "info" ).set_target ( &info_log );

    std::string error_file = root_path + "/log/wavetank_error.log";
    std::ofstream error_log ( error_file.c_str ( ) );
    LogKeeper::get_log ( "error" ).set_target ( &error_log );

    if ( CREATE_DEFAULT_OUTPUT_DIRECTORIES )
    {
        //create log directory  
        std::string log_path;
        size_t last_slash_idx = info_file.rfind ( '/' );
        if ( std::string::npos != last_slash_idx )
        {
            log_path = info_file.substr ( 0, last_slash_idx + 1 );
        }
        std::string command = "mkdir -p " + log_path + ";";
        bool ret = system ( command.c_str ( ) );
        if ( ret != 0 ) LOG_ERROR ( "ERROR in line " << __LINE__ << ", error code = " << ret << std::endl; );

        //create directory for initial mesh
        command = "mkdir -p " + root_path + "/mesh;";
        ret = system ( command.c_str ( ) );
        if ( ret != 0 ) LOG_ERROR ( "ERROR in line " << __LINE__ << ", error code = " << ret << std::endl; );

        //create directory for snapshots
        command = "mkdir -p " + root_path + "/snapshots;";
        ret = system ( command.c_str ( ) );
        if ( ret != 0 ) LOG_ERROR ( "ERROR in line " << __LINE__ << ", error code = " << ret << std::endl; );
        //create directory for output
        command = "mkdir -p " + root_path + "/out;";
        ret = system ( command.c_str ( ) );
        if ( ret != 0 ) LOG_ERROR ( "ERROR in line " << __LINE__ << ", error code = " << ret << std::endl; );
        //create directory for initial condition output
        command = "mkdir -p " + root_path + "/start;";
        ret = system ( command.c_str ( ) );
        if ( ret != 0 ) LOG_ERROR ( "ERROR in line " << __LINE__ << ", error code = " << ret << std::endl; );
    }

    // create application
    MetFlowBousCyl3dApp application ( root_path, path, base_path, resume_run );
    // run standard loop
    application.run ( );

    // finalize
    MPI_Finalize ( );

    return 0;
}
