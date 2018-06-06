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
/// IN-TIME LOOP FUNCTIONS
/// ***************************************************************************

void MetFlowConvDiffApp::update_assembler ( )
{
    this->update_convection ( this->cur_time_, false );
    T_local_asm_.setup_fe_evaluators ( );
}

void MetFlowConvDiffApp::update_assembler_dual ( )
{
    this->update_convection ( this->cur_time_, true );
    T_local_asm_dual_.setup_fe_evaluators ( );
}

void MetFlowConvDiffApp::filter_solution ( )
{
}

void MetFlowConvDiffApp::post_processing ( )
{
    int compute_quant = params_["PostProcessing"]["FlowQuantities"].get<int>( );
    int compute_norm = params_["PostProcessing"]["Norm"].get<int>( );

    if ( compute_quant == 1 )
        MetFlowConvDiffApp::compute_flow_quantities ( "/log/FlowCharacteristics", this->solP_, this->cur_time_ );
    if ( compute_norm == 1 )
        MetFlowConvDiffApp::compute_solution_norm ( "/log/Norm2", this->solP_, this->solP_prev_, this->cur_time_ );
}

void MetFlowConvDiffApp::post_processing_dual ( )
{
    int compute_quant = params_["PostProcessing"]["FlowQuantities"].get<int>( );
    int compute_norm = params_["PostProcessing"]["Norm"].get<int>( );

    if ( compute_quant == 1 )
        MetFlowConvDiffApp::compute_flow_quantities ( "/log/DualFlowCharacteristics", this->solP_, this->cur_time_ );
    if ( compute_norm == 1 )
        MetFlowConvDiffApp::compute_solution_norm ( "/log/DualNorm2", this->solD_, this->solD_next_, this->cur_time_ );
}

void MetFlowConvDiffApp::update_convection ( double time, bool backward )
{
    std::string conv_type = params_["ExperimentSetup"]["Convection"]["Type"].get<std::string>( );
    std::vector<double> const_conv ( 3, 0. );
    std::vector<double> rot_normal ( 3, 0. );
    std::vector<double> rot_center ( 3, 0. );

    double conv_mag = params_["ExperimentSetup"]["Convection"]["Magnitude"].get<double>( );
    const_conv[0] = params_["ExperimentSetup"]["Convection"]["Conv_x"].get<double>( );
    const_conv[1] = params_["ExperimentSetup"]["Convection"]["Conv_y"].get<double>( );
    const_conv[2] = params_["ExperimentSetup"]["Convection"]["Conv_z"].get<double>( );

    rot_normal[0] = params_["ExperimentSetup"]["Convection"]["Normal_x"].get<double>( );
    rot_normal[1] = params_["ExperimentSetup"]["Convection"]["Normal_y"].get<double>( );
    rot_normal[2] = params_["ExperimentSetup"]["Convection"]["Normal_z"].get<double>( );

    rot_center[0] = params_["ExperimentSetup"]["Convection"]["Center_x"].get<double>( );
    rot_center[1] = params_["ExperimentSetup"]["Convection"]["Center_y"].get<double>( );
    rot_center[2] = params_["ExperimentSetup"]["Convection"]["Center_z"].get<double>( );

    if ( conv_type == "Zero" )
    {
        Vec<DIM, double> constant_conv;
        for ( int l = 0; l < DIM; ++l )
        {
            constant_conv[l] = 0.;
        }
        ConstantConvectionTerm<DIM, double>* convection_struct = new ConstantConvectionTerm<DIM, double> ( constant_conv );

        T_local_asm_ .set_convection_term ( convection_struct );
        T_local_asm_dual_.set_convection_term ( convection_struct );
        T_local_asm_est_ .set_convection_term ( convection_struct );
        return;
    }

    if ( conv_type == "Constant" )
    {
        Vec<DIM, double> constant_conv;
        for ( int l = 0; l < DIM; ++l )
        {
            constant_conv[l] = const_conv[l] * conv_mag;
        }
        ConstantConvectionTerm<DIM, double>* convection_struct = new ConstantConvectionTerm<DIM, double> ( constant_conv );

        T_local_asm_ .set_convection_term ( convection_struct );
        T_local_asm_dual_.set_convection_term ( convection_struct );
        T_local_asm_est_ .set_convection_term ( convection_struct );
        return;
    }

    if ( conv_type == "Rotation" )
    {
        Vec<3, double> normal;
        Vec<3, double> center;

        for ( int l = 0; l < 3; ++l )
        {
            normal[l] = rot_normal[l] * conv_mag;
            center[l] = rot_center[l];
        }
        RotationConvectionTerm<DIM, double>* convection_struct = new RotationConvectionTerm<DIM, double> ( normal, center );

        T_local_asm_ .set_convection_term ( convection_struct );
        T_local_asm_dual_.set_convection_term ( convection_struct );
        T_local_asm_est_ .set_convection_term ( convection_struct );
        return;
    }

    assert ( false );

    LAD::VectorType* cur_conv;

    if ( backward )
    {
        cur_conv = &this->V_prev_;
        this->V_next_.CloneFrom ( this->V_ );
        this->V_.CloneFrom ( this->V_prev_ );
    }
    else
    {
        cur_conv = &this->V_;
        this->V_prev_.CloneFrom ( this->V_ );
    }

    if ( conv_type != "Load" && conv_type != "Zero" )
    {
        for ( int var = 0; var < DIM; ++var )
        {
            for ( mesh::EntityIterator it = V_space_.mesh ( ).begin ( DIM ), end_it = V_space_.mesh ( ).end ( DIM ); it != end_it; ++it )
            {
                std::vector<int> global_dof_ids;
                V_space_.GetDofIndices ( var, *it, &global_dof_ids );
                int num_dofs = global_dof_ids.size ( );
                std::vector<double> values;
                values.resize ( num_dofs, 0. );

                std::vector< Coord > coords;
                V_space_.dof ( ).get_coord_on_cell ( var, it->index ( ), coords );

                for ( int i = 0; i < num_dofs; i++ )
                {
                    if ( rank_ == V_space_.dof ( ).owner_of_dof ( global_dof_ids.at ( i ) ) )
                    {
                        std::vector<double> r ( 3, 0. );
                        r[0] = coords[i][0] - rot_center[0];
                        r[1] = coords[i][1] - rot_center[1];

                        if ( DIM > 2 )
                        {
                            r[2] = coords[i][2] - rot_center[2];
                        }
                        double val = 0.;

                        if ( conv_type == "Rotation" )
                        {
                            switch ( var )
                            {
                                case 0:
                                    val = rot_normal[1] * r[2] - rot_normal[2] * r[1];
                                    break;
                                case 1:
                                    val = rot_normal[2] * r[0] - rot_normal[0] * r[2];
                                    break;
                                default:
                                    val = rot_normal[0] * r[1] - rot_normal[1] * r[0];
                                    break;
                            }

                        }
                        else if ( conv_type == "Constant" )
                        {
                            val = const_conv[var];
                        }
                        val *= conv_mag;
                        cur_conv->SetValues ( &global_dof_ids.at ( i ), 1, &val );
                    }
                }
            }
        }
    }
    V_.Update ( );
    V_prev_.Update ( );
    V_next_.Update ( );
}
