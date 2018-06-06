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
/// IO
/// ***************************************************************************

void MetFlowConvDiffApp::create_log_files ( int offset )
{
    this->create_log_file ( offset, "/log/FlowCharacteristics", "txt" );
    this->create_log_file ( offset, "/log/NonlinearSolver", "txt" );
    this->create_log_file ( offset, "/log/LinearSolver", "txt" );
    this->create_log_file ( offset, "/log/Norm2", "txt" );
}

void MetFlowConvDiffApp::create_log_files_dual ( int offset )
{
    this->create_log_file ( offset, "/log/DualFlowCharacteristics", "txt" );
    this->create_log_file ( offset, "/log/DualNorm2", "txt" );
}

/// Visualize

void MetFlowConvDiffApp::visualize_solution ( int step )
{
    Timer timer;
    timer.start ( );

    std::string prefix = this->base_params_["Visualization"]["FilePrefix"].get<std::string>( );
    std::stringstream pre;
    pre << this->root_ << "/" << prefix << "." << this->adapt_counter_;
    this->visualize_solution ( *this->solP_, this->space_, pre.str ( ), step );
    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  took " << timer.get_duration ( ) << "sec" << std::endl;
}

/// Visualize

void MetFlowConvDiffApp::visualize_solution ( LAD::VectorType& sol, const VectorSpace<double>* space, std::string const& prefix, int time_step )
{

    if ( rank ( ) == master_rank ( ) ) std::cout << "> Visualize primal ConvDiff solution at time step " << time_step << std::endl;

    std::stringstream input;
    std::stringstream input_conv;
    bool write_cart = false;

    if ( time_step < 10 )
    {
        input << prefix << ".000" << time_step;
        input_conv << prefix << "_conv" << ".000" << time_step;
    }
    else if ( time_step < 100 )
    {
        input << prefix << ".00" << time_step;
        input_conv << prefix << "_conv" << ".00" << time_step;
    }
    else if ( time_step < 1000 )
    {
        input << prefix << ".0" << time_step;
        input_conv << prefix << "_conv" << ".0" << time_step;
    }
    else
    {
        input << prefix << "." << time_step;
        input_conv << prefix << "_conv" << "." << time_step;
    }

    if ( num_partitions_ > 1 )
    {
        input << ".pvtu";
        input_conv << ".pvtu";
    }
    else
    {
        input << ".vtu";
        input_conv << ".vtu";
    }

    LAD::VectorType tmp;
    std::string visu_filename = input.str ( );

    ParallelCellVisualization<double> visu ( *space, 1, comm_, master_rank ( ) );
    //ParallelCellVisualization<double>  V_visu ( V_space_, 1, comm_, master_rank() );
    tmp.CloneFrom ( sol );

#ifdef USE_HYPRE
    tmp.Update ( );
#else
    tmp.UpdateCouplings ( );
#endif

    visu.visualize ( EvalFeFunction<LAD>( *space, tmp, t_var_ ), "T" );
    visu.visualize ( EvalDerivativeFeFunction<LAD, DIM>( *space, tmp, t_var_, 0 ), "dT_dx" );
    visu.visualize ( EvalDerivativeFeFunction<LAD, DIM>( *space, tmp, t_var_, 1 ), "dT_dy" );
    if ( DIM > 2 )
    {
        visu.visualize ( EvalDerivativeFeFunction<LAD, DIM>( *space, tmp, t_var_, 2 ), "dT_dz" );
    }
    /*
    V_visu.visualize(EvalFeFunction<LAD>(V_space_ , this->V_ , 0), "v_x");
    V_visu.visualize(EvalFeFunction<LAD>(V_space_ , this->V_ , 1), "v_y");
    if (DIM > 2)
        V_visu.visualize(EvalFeFunction<LAD>(V_space_ , this->V_ , 2), "v_z");
     */
    T_local_asm_.set_solP ( *this->solP_ );
    T_local_asm_.set_solP_prev ( *this->solP_prev_ );

    // parallel statistics
    std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    AttributePtr sub = mesh_->get_attribute ( "_sub_domain_", mesh_->tdim ( ) );
    AttributePtr remote = mesh_->get_attribute ( "_remote_index_", mesh_->tdim ( ) );
    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) ); it != mesh_->end ( mesh_->tdim ( ) ); ++it )
    {
        remote_index.at ( it->index ( ) ) = remote->get_int_value ( it->index ( ) );
        sub_domain.at ( it->index ( ) ) = sub->get_int_value ( it->index ( ) );
    }
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );

    // write file
    visu.write ( visu_filename );
    //V_visu.write(input_conv.str());
    tmp.Clear ( );
}

/// Visualize

void MetFlowConvDiffApp::visualize_solution_dual ( int step )
{
    Timer timer;
    timer.start ( );

    std::string prefix = this->base_params_["Visualization"]["DualFilePrefix"].get<std::string>( );
    std::stringstream pre;
    pre << this->root_ << "/" << prefix << "." << this->adapt_counter_;

    std::vector<std::string> var_names;
    var_names.push_back ( "dual_T" );
    this->visualize_function ( *this->solD_, this->space_dual_, pre.str ( ), step, var_names );

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) ) std::cout << "  took " << timer.get_duration ( ) << "sec" << std::endl;
}
