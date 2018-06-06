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

#include "met_flow_convdiff_assembler.h"

template<int DIM, class DataType>
MetFlowConvDiffAssembler<DIM, DataType>::MetFlowConvDiffAssembler ( )
: MetFlowAssembler<DIM, DataType>( )
{
    this->kappa_ = -999999999.;
    this->gamma_ = -999999999.;
    this->lambda_ = -999999999.;
    this->vector_asm_mode_ = VECTOR_STD;
    this->galerkin_mode_ = 0;
    this->set_time_discretization_to_theta ( 0.5 );

    this->L2_perturb_.resize ( 1, false );
    this->H1_perturb_.resize ( 1, false );
    this->set_time_discretization_to_zero ( );
    this->num_var_ = 1;
    MetFlowAssembler<DIM, DataType>::allocate_function_values ( this->num_var_ );
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::clear ( )
{
    MetFlowAssembler<DIM, DataType>::clear ( );
    this->kappa_ = -999999999.;
    this->gamma_ = -999999999.;
    this->lambda_ = -999999999.;
    this->vector_asm_mode_ = VECTOR_STD;
    this->galerkin_mode_ = 0;
    this->set_time_discretization_to_theta ( 0.5 );

    this->L2_perturb_.resize ( 1, false );
    this->H1_perturb_.resize ( 1, false );
    this->set_time_discretization_to_zero ( );
    this->num_var_ = 1;
    MetFlowAssembler<DIM, DataType>::allocate_function_values ( this->num_var_ );
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::set_time_discretization_to_zero ( )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_zero ( );
    this->theta_diff_c_ = 0.;
    this->theta_diff_p_ = 0.;
    this->theta_adv_cc_ = 0.;
    this->theta_adv_pc_ = 0.;
    this->theta_adv_cp_ = 0.;
    this->theta_adv_pp_ = 0.;

    this->theta_rea_c_ = 0.;
    this->theta_rea_p_ = 0.;
    this->delta_diff_c_ = 0.;
    this->delta_diff_n_ = 0.;
    this->delta_adv_cc_ = 0.;
    this->delta_adv_nc_ = 0.;
    this->delta_adv_cn_ = 0.;
    this->delta_adv_nn_ = 0.;
    this->delta_adv_pc_ = 0.;
    this->delta_rea_c_ = 0.;
    this->delta_rea_n_ = 0.;
    this->theta_sou_c_ = 0.;
    this->theta_sou_p_ = 0.;
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::set_time_discretization_to_scaling ( DataType beta1, DataType beta2, DataType beta3, DataType beta4 )
{
    this->set_time_discretization_to_zero ( );

    this->dT_pc_ = 1.;
    this->dT_cn_ = 1.;
    this->theta_d_dt_u_ = beta1;
    this->theta_diff_c_ = beta2;
    this->kappa_ = 1.;
    this->theta_adv_cc_ = beta3;
    this->theta_adv_pc_ = beta4;
    this->gamma_ = 1.;
    this->lambda_ = 0.;
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::set_time_discretization_to_simple ( )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_simple ( );
    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::set_time_discretization_to_theta ( DataType theta )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_theta ( theta );

    // Primal Mode
    if ( this->mode_ == PRIMAL )
    {
        this->theta_diff_c_ = theta;
        this->theta_diff_p_ = 1. - theta;

        this->theta_adv_cc_ = theta;
        this->theta_adv_pc_ = 0.;
        this->theta_adv_cp_ = 0.;
        this->theta_adv_pp_ = 1 - theta;

        this->theta_rea_c_ = theta;
        this->theta_rea_p_ = 1. - theta;

        this->theta_sou_c_ = theta;
        this->theta_sou_p_ = 1. - theta;
    }
    // Dual Mode 
    if ( this->mode_ == DUAL )
    {
        this->delta_diff_c_ = theta;
        this->delta_diff_n_ = 1. - theta;

        this->delta_adv_cc_ = theta;
        this->delta_adv_nc_ = 0.;
        this->delta_adv_cn_ = 0.;
        this->delta_adv_nn_ = 1. - theta;
        this->delta_adv_pc_ = 0.;

        this->delta_rea_c_ = theta;
        this->delta_rea_n_ = 1. - theta;
    }

    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::set_time_discretization_to_galerkin_cd ( bool modified )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_galerkin_cd ( modified );

    // Primal Mode
    if ( this->mode_ == PRIMAL )
    {
        this->theta_diff_c_ = 0.5;
        this->theta_diff_p_ = 0.5;
        this->theta_adv_cc_ = 1. / 3.;
        this->theta_adv_pc_ = 1. / 6.;
        this->theta_adv_cp_ = 1. / 6.;
        this->theta_adv_pp_ = 1. / 3.;
        this->theta_rea_c_ = 0.5;
        this->theta_rea_p_ = 0.5;
        this->theta_sou_c_ = 0.5;
        this->theta_sou_p_ = 0.5;
    }

    // Dual Mode
    if ( this->mode_ == DUAL )
    {
        this->delta_diff_c_ = 0.5;
        this->delta_diff_n_ = 0.5;

        this->delta_adv_cc_ = 1. / 3.;
        this->delta_adv_nc_ = 1. / 6.;
        this->delta_adv_cn_ = 1. / 6.;
        this->delta_adv_nn_ = 1. / 3.;
        this->delta_adv_pc_ = 0.;

        this->delta_rea_c_ = 0.5;
        this->delta_rea_n_ = 0.5;
    }
    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::set_time_discretization_to_galerkin_dc ( bool modified )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_galerkin_dc ( modified );

    if ( this->mode_ == DUAL )
    {
        this->delta_diff_c_ = 0.5;
        this->delta_diff_n_ = 0.5;

        this->delta_adv_cc_ = 1. / 3.;
        this->delta_adv_nc_ = 0.;
        this->delta_adv_cn_ = 1. / 3.;
        this->delta_adv_nn_ = 1. / 6.;
        this->delta_adv_pc_ = 1. / 6.;

        this->delta_rea_c_ = 0.5;
        this->delta_rea_n_ = 0.5;
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::set_time_discretization_off ( )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_off ( );

    this->theta_diff_c_ = 1.;
    this->theta_diff_p_ = 0.;
    this->theta_adv_cc_ = 1.;
    this->theta_adv_pc_ = 0.;
    this->theta_adv_cp_ = 0.;
    this->theta_adv_pp_ = 0.;
    this->theta_rea_c_ = 1.;
    this->theta_rea_p_ = 0.;
    this->theta_sou_c_ = 1.;
    this->theta_sou_p_ = 0.;

    if ( this->mode_ == DUAL )
    {
        // TODO
        this->delta_diff_c_ = 1.;
        this->delta_diff_n_ = 0.;

        this->delta_adv_cc_ = 1.;
        this->delta_adv_nc_ = 0.;
        this->delta_adv_cn_ = 0.;
        this->delta_adv_nn_ = 0.;
        this->delta_adv_pc_ = 0.;

        this->delta_rea_c_ = 1.;
        this->delta_rea_n_ = 0.;
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::print_parameters ( )
{
    MetFlowAssembler<DIM, DataType>::print_parameters ( );

    std::cout.scientific;
    std::cout.precision ( 6 );
    if ( this->mode_ == PRIMAL )
        std::cout << "> MetFlowConvDiffAssembler primal parameters" << std::endl;
    else
        std::cout << "> MetFlowConvDiffAssembler dual parameters" << std::endl;

    std::cout << "  kappa:   " << this->kappa_ << std::endl;
    std::cout << "  lambda:  " << this->lambda_ << std::endl;
    std::cout << "  gamma:   " << this->gamma_ << std::endl;
    std::cout << std::endl;

    this->print_coeff ( );
    this->print_perturb ( );
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::print_coeff ( )
{
    MetFlowAssembler<DIM, DataType>::print_coeff ( );

    if ( this->mode_ == PRIMAL )
    {
        std::cout << "  theta_diff_c:   " << this->theta_diff_c_ << std::endl;
        std::cout << "  theta_diff_n:   " << this->theta_diff_p_ << std::endl;
        std::cout << "  theta_adv_cc:   " << this->theta_adv_cc_ << std::endl;
        std::cout << "  theta_adv_pc:   " << this->theta_adv_pc_ << std::endl;
        std::cout << "  theta_adv_cp:   " << this->theta_adv_cp_ << std::endl;
        std::cout << "  theta_adv_pp:   " << this->theta_adv_pp_ << std::endl;
        std::cout << "  theta_rea_c:    " << this->theta_rea_c_ << std::endl;
        std::cout << "  theta_rea_p:    " << this->theta_rea_p_ << std::endl;
        std::cout << "  theta_sou_c:    " << this->theta_sou_c_ << std::endl;
        std::cout << "  theta_sou_p:    " << this->theta_sou_p_ << std::endl;
    }
    else
    {
        std::cout << "  delta_diff_n:   " << this->delta_diff_n_ << std::endl;
        std::cout << "  delta_diff_c:   " << this->delta_diff_c_ << std::endl;
        std::cout << "  delta_adv_nn:   " << this->delta_adv_nn_ << std::endl;
        std::cout << "  delta_adv_cn:   " << this->delta_adv_cn_ << std::endl;
        std::cout << "  delta_adv_nc:   " << this->delta_adv_nc_ << std::endl;
        std::cout << "  delta_adv_cc:   " << this->delta_adv_cc_ << std::endl;
        std::cout << "  delta_adv_pc:   " << this->delta_adv_pc_ << std::endl;
        std::cout << "  delta_rea_n:    " << this->delta_rea_n_ << std::endl;
        std::cout << "  delta_rea_c:    " << this->delta_rea_c_ << std::endl;
    }
}

template<int DIM, class DataType>
DataType MetFlowConvDiffAssembler<DIM, DataType>::get_coeff ( std::string term )
{
    DataType ret = 9999.;
    DataType dt = this->dT_pc_;

    if ( term == "p_dt" ) ret = this->theta_d_dt_u_;
    if ( term == "p_diff_c" ) ret = this->theta_diff_c_ * this->kappa_ * dt;
    if ( term == "p_adv_cc" ) ret = this->theta_adv_cc_ * this->gamma_ * dt;
    if ( term == "p_adv_pc" ) ret = this->theta_adv_pc_ * this->gamma_ * dt;
    if ( term == "p_adv_cp" ) ret = this->theta_adv_cp_ * this->gamma_ * dt;
    if ( term == "p_adv_pp" ) ret = this->theta_adv_pp_ * this->gamma_ * dt;

    if ( this->galerkin_mode_ == 0 )
        dt = this->dT_cn_;

    if ( term == "d_dt" ) ret = this->delta_d_dt_u_;
    if ( term == "d_diff_c" ) ret = this->delta_diff_c_ * this->kappa_ * dt;
    if ( term == "d_adv_cc" ) ret = -1. * this->delta_adv_cc_ * this->gamma_ * dt;
    if ( term == "d_adv_nc" ) ret = -1. * this->delta_adv_nc_ * this->gamma_ * dt;
    if ( term == "d_adv_cn" ) ret = -1. * this->delta_adv_cn_ * this->gamma_ * dt;
    if ( term == "d_adv_nn" ) ret = -1. * this->delta_adv_nn_ * this->gamma_ * dt;
    if ( term == "d_adv_pc" ) ret = -1. * this->delta_adv_pc_ * this->gamma_ * dt;

    if ( this->theta_d_dt_u_ == 0. )
        ret /= dt;

    return ret;
}

/// ********************************************************
/// Assembly routines 
/// ********************************************************

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    if ( this->mode_ == PRIMAL )
    {
        this->assemble_local_matrix_primal ( element, lm );
    }
    else if ( this->mode_ == DUAL )
    {
        this->assemble_local_matrix_dual ( element, lm );
    }
    else
    {
        interminable_assert ( 0 );
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const
{
    if ( this->vector_asm_mode_ == VECTOR_STD )
    {
        if ( this->mode_ == PRIMAL )
        {
            this->assemble_local_vector_primal ( element, lv );
        }
        else if ( this->mode_ == DUAL )
        {
            this->assemble_local_vector_dual ( element, lv );
        }
        else
        {
            interminable_assert ( 0 );
        }
    }
    else if ( this->vector_asm_mode_ == VECTOR_GOAL )
    {
        this->assemble_local_vector_goal ( element, lv );
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const
{
}

template<int DIM, class DataType>
void MetFlowConvDiffAssembler<DIM, DataType>::assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const
{
}

template class MetFlowConvDiffAssembler<2, double>;
template class MetFlowConvDiffAssembler<3, double>;
