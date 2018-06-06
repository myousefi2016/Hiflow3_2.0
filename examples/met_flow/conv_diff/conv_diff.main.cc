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
/// MAIN LOOPS
/// ***************************************************************************
/// Post Process HDF5 data -> RUN_MODE3

void MetFlowConvDiffApp::pp_run ( )
{

    // Initialize problem
    this->initial_prepare ( );

    std::string mode = this->base_params_["PostProcessing"]["Mode"].get<std::string>( "Primal" );
    int compute_quant = params_["PostProcessing"]["FlowQuantities"].get<int>( 1 );
    int compute_norm = params_["PostProcessing"]["Norm"].get<int>( 0 );
    std::string filename_quant;
    std::string filename_norm;

    if ( rank ( ) == master_rank ( ) )
    {
        if ( mode == "Primal" )
        {
            filename_quant = "/log/FlowCharacteristics";
            filename_norm = "/log/Norm2";
        }
        if ( mode == "Dual" )
        {
            filename_quant = "/log/DualFlowCharacteristics";
            filename_norm = "/log/DualNorm2";
        }
    }

    this->create_log_file ( 0, filename_quant, "txt" );
    this->create_log_file ( 0, filename_norm, "txt" );

    int frequence = this->base_params_["PostProcessing"]["PPTimeStep"].get<int>( 10 );
    int initial_step = this->base_params_["PostProcessing"]["InitialStep"].get<int>( 0 );
    int final_step = this->base_params_["PostProcessing"]["FinalStep"].get<int>( 1000 );
    int visu = this->base_params_["PostProcessing"]["Visualize"].get<int>( 0 );
    int visu_freq = this->base_params_["PostProcessing"]["VisTimeStep"].get<int>( 10 );

    // For HDF5
    std::string filename = this->base_params_["PostProcessing"]["SnapshotsIn"].get<std::string>( );
    std::string groupname = this->base_params_["PostProcessing"]["SnapshotsGroup"].get<std::string>( );
    std::string prefix = this->base_params_["PostProcessing"]["SnapshotsPrefix"].get<std::string>( );
    filename = this->root_ + "/" + filename;

    for ( int k = initial_step; k <= final_step; ++k )
    {
        int l = k;
        if ( mode == "Dual" )
        {
            l = final_step - k;
        }
        this->cur_time_ = this->get_time ( l );
        if ( l % frequence == 0 )
        {
            std::stringstream ss;
            ss << prefix << l;
            std::stringstream ss_prev;
            if ( l >= 1 )
            {
                ss_prev << prefix << l - 1;
            }
            else
            {
                ss_prev << prefix << l;
            }
            MetFlowApp::read_file ( *this->solP_prev_, filename, groupname, ss_prev.str ( ) );
            MetFlowApp::read_file ( *this->solP_, filename, groupname, ss.str ( ) );

            this->solP_prev_->Update ( );
            this->solP_->Update ( );

            T_local_asm_.set_solP ( *this->solP_ );
            T_local_asm_.set_solP_prev ( *this->solP_prev_ );

            if ( compute_quant == 1 )
                MetFlowConvDiffApp::compute_flow_quantities ( filename_quant, this->solP_, this->cur_time_ );

            if ( compute_norm == 1 )
                MetFlowConvDiffApp::compute_solution_norm ( filename_norm, this->solP_, this->solP_prev_, this->cur_time_ );

            if ( visu == 1 )
            {
                if ( l % visu_freq == 0 )
                {
                    if ( mode == "Primal" )
                        MetFlowConvDiffApp::visualize_solution ( l );
                    if ( mode == "Dual" )
                    {
                        this->solD_next_->CloneFrom ( *this->solP_ );
                        this->solD_->CloneFrom ( *this->solP_prev_ );
                        MetFlowConvDiffApp::visualize_solution_dual ( l );
                    }
                }
            }
            if ( rank ( ) == master_rank ( ) )
                std::cout << "  " << l << " / " << final_step << " done " << std::endl;
        }
    }

}
