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

#include "met_flow_incomp_assembler.h"

template<int DIM, class DataType>
MetFlowIncompAssembler<DIM, DataType>::MetFlowIncompAssembler ( )
: MetFlowAssembler<DIM, DataType>( )
{
    this->nu_ = -999999999.;
    this->rho_ = -999999999.;
    this->inv_rho_ = -999999999.;
    this->omega_ = -999999999.;

    this->skew_mode_ = SS_OFF;
    this->graddiv_mode_ = GD_OFF;
    this->graddiv_gamma_ = 0.;
    this->graddiv_C_ = 0.;
    this->vector_asm_mode_ = VECTOR_STD;
    this->galerkin_mode_ = 0;
    this->conv_mode_ = NAVIERSTOKES;

#ifdef AUGMENT_PRESS
    this->num_var_ = DIM + 2;
#else
    this->num_var_ = DIM + 1;
#endif

    this->L2_perturb_.resize ( this->num_var_, false );
    this->H1_perturb_.resize ( this->num_var_, false );

    this->set_time_discretization_to_zero ( );
    MetFlowAssembler<DIM, DataType>::allocate_function_values ( this->num_var_ );
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::clear ( )
{
    MetFlowAssembler<DIM, DataType>::clear ( );

    this->nu_ = -999999999.;
    this->rho_ = -999999999.;
    this->inv_rho_ = -999999999.;
    this->omega_ = -999999999.;

    this->skew_mode_ = SS_OFF;
    this->graddiv_mode_ = GD_OFF;
    this->graddiv_gamma_ = 0.;
    this->graddiv_C_ = 0.;
    this->vector_asm_mode_ = VECTOR_STD;
    this->galerkin_mode_ = 0;
    this->conv_mode_ = NAVIERSTOKES;

#ifdef AUGMENT_PRESS
    this->num_var_ = DIM + 2;
#else
    this->num_var_ = DIM + 1;
#endif

    this->L2_perturb_.resize ( this->num_var_, false );
    this->H1_perturb_.resize ( this->num_var_, false );

    this->set_time_discretization_to_zero ( );

    MetFlowAssembler<DIM, DataType>::allocate_function_values ( this->num_var_ );
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_zero ( )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_zero ( );
    this->theta_mom_vis_c_ = 0.;
    this->theta_mom_vis_p_ = 0.;
    this->theta_mom_adv_cc_ = 0.;
    this->theta_mom_adv_pc_ = 0.;
    this->theta_mom_adv_cp_ = 0.;
    this->theta_mom_adv_pp_ = 0.;
    this->theta_mom_pre_c_ = 0.;
    this->theta_mom_rot_c_ = 0.;
    this->theta_mom_rot_p_ = 0.;
    this->theta_mom_graddiv_c_ = 0.;
    this->theta_mom_graddiv_p_ = 0.;
    this->theta_inc_c_ = 0.;
    this->theta_inc_p_ = 0.;

    this->delta_mom_vis_c_ = 0.;
    this->delta_mom_vis_n_ = 0.;
    this->delta_mom_adv_cc_ = 0.;
    this->delta_mom_adv_nc_ = 0.;
    this->delta_mom_adv_cn_ = 0.;
    this->delta_mom_adv_nn_ = 0.;
    this->delta_mom_pre_c_ = 0.;
    this->delta_mom_rot_c_ = 0.;
    this->delta_mom_rot_n_ = 0.;
    this->delta_mom_graddiv_c_ = 0.;
    this->delta_mom_graddiv_n_ = 0.;
    this->delta_mom_rea_cc_ = 0.;
    this->delta_mom_rea_nc_ = 0.;
    this->delta_mom_rea_cn_ = 0.;
    this->delta_mom_rea_nn_ = 0.;
    this->delta_inc_c_ = 0.;
    this->delta_inc_n_ = 0.;
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_scaling ( DataType beta1, DataType beta2,
                                                                                 DataType beta3, DataType beta4,
                                                                                 DataType beta5, DataType beta6 )
{
    this->set_time_discretization_to_zero ( );
    this->dT_pc_ = 1.;
    this->dT_cn_ = 1.;
    this->theta_d_dt_u_ = beta1;
    this->theta_mom_vis_c_ = beta2;
    this->nu_ = 1.;
    this->theta_mom_adv_cc_ = beta3;
    this->theta_mom_adv_pc_ = beta4;
    this->theta_mom_pre_c_ = beta5;
    this->inv_rho_ = 1.;
    this->theta_inc_c_ = beta6;
    this->skew_mode_ = SS_OFF;
    this->graddiv_mode_ = GD_OFF;
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_simple ( )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_simple ( );

    this->delta_mom_vis_c_ = 0.5;
    this->delta_mom_vis_n_ = 0.5;
    this->delta_mom_adv_cc_ = 1. / 3.;
    this->delta_mom_adv_nc_ = 1. / 6.;
    this->delta_mom_adv_cn_ = 1. / 6.;
    this->delta_mom_adv_nn_ = 1. / 3.;
    this->delta_mom_pre_c_ = 1.;
    this->delta_mom_rot_c_ = 0.5;
    this->delta_mom_rot_n_ = 0.5;
    this->delta_mom_graddiv_c_ = 1.;
    this->delta_mom_graddiv_n_ = 0.;
    this->delta_mom_rea_cc_ = 1. / 3.;
    this->delta_mom_rea_nc_ = 1. / 6.;
    this->delta_mom_rea_cn_ = 1. / 6.;
    this->delta_mom_rea_nn_ = 1. / 3.;

    this->delta_inc_c_ = 0.5;
    this->delta_inc_n_ = 0.5;

    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_theta ( DataType theta )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_theta ( theta );

    // Primal Mode
    if ( this->mode_ == PRIMAL )
    {
        this->theta_mom_vis_c_ = theta;
        this->theta_mom_vis_p_ = 1. - theta;

        this->theta_mom_adv_cc_ = theta;
        this->theta_mom_adv_pc_ = 0.;
        this->theta_mom_adv_cp_ = 0.;
        this->theta_mom_adv_pp_ = 1 - theta;

        this->theta_mom_pre_c_ = 1.;

        this->theta_mom_rot_c_ = theta;
        this->theta_mom_rot_p_ = 1. - theta;

        this->theta_mom_graddiv_c_ = 1.;
        this->theta_mom_graddiv_p_ = 0.;

        this->theta_inc_c_ = 1. / this->dT_pc_;
        this->theta_inc_p_ = 0.;
    }
    // Dual Mode TODO
    if ( this->mode_ == DUAL )
    {
        this->delta_mom_vis_c_ = 1. - theta;
        this->delta_mom_vis_n_ = theta;

        this->delta_mom_adv_cc_ = 1. - theta;
        this->delta_mom_adv_nc_ = 0.;
        this->delta_mom_adv_cn_ = 0.;
        this->delta_mom_adv_nn_ = theta;
        this->delta_mom_adv_pc_ = 0.;

        this->delta_mom_pre_c_ = 1.;
        this->delta_mom_pre_n_ = 0.;

        this->delta_mom_rot_c_ = 1. - theta;
        this->delta_mom_rot_n_ = theta;

        this->delta_mom_graddiv_c_ = 0.;
        this->delta_mom_graddiv_n_ = 1.;

        this->delta_mom_rea_cc_ = 1. - theta;
        this->delta_mom_rea_nc_ = 0.;
        this->delta_mom_rea_cn_ = 0.;
        this->delta_mom_rea_nn_ = theta;
        this->delta_mom_rea_pc_ = 0.;

        this->delta_inc_c_ = 1. - theta;
        this->delta_inc_n_ = theta;
    }

    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_galerkin_cd ( bool modified )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_galerkin_cd ( modified );

    // Primal Mode
    if ( this->mode_ == PRIMAL )
    {
        this->theta_mom_vis_c_ = 0.5;
        this->theta_mom_vis_p_ = 0.5;
        this->theta_mom_adv_cc_ = 1. / 3.;
        this->theta_mom_adv_pc_ = 1. / 6.;
        this->theta_mom_adv_cp_ = 1. / 6.;
        this->theta_mom_adv_pp_ = 1. / 3.;
        this->theta_mom_pre_c_ = 1.;
        this->theta_mom_rot_c_ = 0.5;
        this->theta_mom_rot_p_ = 0.5;
        this->theta_mom_graddiv_c_ = 1.;
        this->theta_mom_graddiv_p_ = 0.;

        if ( modified )
        {
            this->theta_inc_c_ = 1.;
            this->theta_inc_p_ = 0.;
        }
        else
        {
            this->theta_inc_c_ = 0.5;
            this->theta_inc_p_ = 0.5;
        }
    }

    // Dual Mode
    if ( this->mode_ == DUAL )
    {
        this->delta_mom_vis_c_ = 0.5;
        this->delta_mom_vis_n_ = 0.5;

        this->delta_mom_adv_cc_ = 1. / 3.;
        this->delta_mom_adv_nc_ = 1. / 6.;
        this->delta_mom_adv_cn_ = 1. / 6.;
        this->delta_mom_adv_nn_ = 1. / 3.;
        this->delta_mom_adv_pc_ = 0.;

        this->delta_mom_pre_c_ = 1.;
        this->delta_mom_pre_n_ = 0.;

        this->delta_mom_rot_c_ = 0.5;
        this->delta_mom_rot_n_ = 0.5;

        this->delta_mom_graddiv_c_ = 1.;
        this->delta_mom_graddiv_n_ = 0.;

        this->delta_mom_rea_cc_ = 1. / 3.;
        this->delta_mom_rea_nc_ = 1. / 6.;
        this->delta_mom_rea_cn_ = 1. / 6.;
        this->delta_mom_rea_nn_ = 1. / 3.;
        this->delta_mom_rea_pc_ = 0.;

        this->delta_inc_c_ = 0.5;
        this->delta_inc_n_ = 0.5;
    }
    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_galerkin_dc ( bool modified )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_to_galerkin_dc ( modified );

    // Dual Mode
    if ( this->mode_ == DUAL )
    {
        this->delta_mom_vis_n_ = 0.5;
        this->delta_mom_vis_c_ = 0.5;

        this->delta_mom_adv_cc_ = 1. / 3.;
        this->delta_mom_adv_nc_ = 0.;
        this->delta_mom_adv_cn_ = 1. / 3.;
        this->delta_mom_adv_nn_ = 1. / 6.;
        this->delta_mom_adv_pc_ = 1. / 6.;

        this->delta_mom_pre_c_ = .5;
        this->delta_mom_pre_n_ = .5;

        this->delta_mom_rot_c_ = 0.5;
        this->delta_mom_rot_n_ = 0.5;

        this->delta_mom_graddiv_c_ = 1.;
        this->delta_mom_graddiv_n_ = 0.;

        this->delta_mom_rea_cc_ = 1. / 3.;
        this->delta_mom_rea_nc_ = 0.;
        this->delta_mom_rea_cn_ = 1. / 3.;
        this->delta_mom_rea_nn_ = 1. / 6.;
        this->delta_mom_rea_pc_ = 1. / 6.;

        this->delta_inc_c_ = 1.;
        this->delta_inc_n_ = 0.;
    }
    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_off ( )
{
    MetFlowAssembler<DIM, DataType>::set_time_discretization_off ( );

    this->theta_mom_vis_c_ = 1.;
    this->theta_mom_vis_p_ = 0.;
    this->theta_mom_adv_cc_ = 1.;
    this->theta_mom_adv_pc_ = 0.;
    this->theta_mom_adv_cp_ = 0.;
    this->theta_mom_adv_pp_ = 0.;
    this->theta_mom_pre_c_ = 1.;
    this->theta_mom_rot_c_ = 1.;
    this->theta_mom_rot_p_ = 0.;
    this->theta_mom_graddiv_c_ = 1.;
    this->theta_mom_graddiv_p_ = 0.;
    this->theta_inc_c_ = 1.;
    this->theta_inc_p_ = 0.;
    this->is_matrix_constant_ = false;

    // TODO
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::print_parameters ( )
{
    MetFlowAssembler<DIM, DataType>::print_parameters ( );

    std::cout.scientific;
    std::cout.precision ( 6 );
    if ( this->mode_ == PRIMAL )
        std::cout << "> MetFlowIncompAssembler primal parameters" << std::endl;
    else
        std::cout << "> MetFlowIncompAssembler dual parameters" << std::endl;

    std::cout << "  nu:       " << this->nu_ << std::endl;
    std::cout << "  inv_rho:  " << this->inv_rho_ << std::endl;
    std::cout << "  omega:    " << this->omega_ << std::endl;
    std::cout << "  Skew symmetric:    " << this->skew_mode_ << std::endl;
    std::cout << "  Grad-Div mode:     " << this->graddiv_mode_ << std::endl;
    std::cout << "  Grad-Div gamma:    " << this->graddiv_gamma_ << std::endl;
    std::cout << "  Grad-Div C:        " << this->graddiv_C_ << std::endl;
    std::cout << std::endl;

    MetFlowIncompAssembler<DIM, DataType>::print_coeff ( );
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::print_coeff ( )
{
    if ( this->mode_ == PRIMAL )
    {
        std::cout << "  theta_mom_vis_c:    " << this->theta_mom_vis_c_ << std::endl;
        std::cout << "  theta_mom_vis_n:    " << this->theta_mom_vis_p_ << std::endl;
        std::cout << "  theta_mom_adv_cc:   " << this->theta_mom_adv_cc_ << std::endl;
        std::cout << "  theta_mom_adv_pc:   " << this->theta_mom_adv_pc_ << std::endl;
        std::cout << "  theta_mom_adv_cp:   " << this->theta_mom_adv_cp_ << std::endl;
        std::cout << "  theta_mom_adv_pp:   " << this->theta_mom_adv_pp_ << std::endl;
        std::cout << "  theta_mom_rot_c:    " << this->theta_mom_rot_c_ << std::endl;
        std::cout << "  theta_mom_rot_p:    " << this->theta_mom_rot_p_ << std::endl;
        std::cout << "  theta_mom_pre_c :   " << this->theta_mom_pre_c_ << std::endl;
        std::cout << "  theta_mom_gradiv_c: " << this->theta_mom_graddiv_c_ << std::endl;
        std::cout << "  theta_mom_gradiv_p: " << this->theta_mom_graddiv_p_ << std::endl;
        std::cout << "  theta_inc_c:        " << this->theta_inc_c_ << std::endl;
        std::cout << "  theta_inc_p:        " << this->theta_inc_p_ << std::endl;
    }
    else
    {
        std::cout << "  delta_mom_vis_n:    " << this->delta_mom_vis_n_ << std::endl;
        std::cout << "  delta_mom_vis_c:    " << this->delta_mom_vis_c_ << std::endl;
        std::cout << "  delta_mom_adv_nn:   " << this->delta_mom_adv_nn_ << std::endl;
        std::cout << "  delta_mom_adv_cn:   " << this->delta_mom_adv_cn_ << std::endl;
        std::cout << "  delta_mom_adv_nc:   " << this->delta_mom_adv_nc_ << std::endl;
        std::cout << "  delta_mom_adv_cc:   " << this->delta_mom_adv_cc_ << std::endl;
        std::cout << "  delta_mom_adv_pc:   " << this->delta_mom_adv_pc_ << std::endl;
        std::cout << "  delta_mom_rea_nn:   " << this->delta_mom_rea_nn_ << std::endl;
        std::cout << "  delta_mom_rea_cn:   " << this->delta_mom_rea_cn_ << std::endl;
        std::cout << "  delta_mom_rea_nc:   " << this->delta_mom_rea_nc_ << std::endl;
        std::cout << "  delta_mom_rea_cc:   " << this->delta_mom_rea_cc_ << std::endl;
        std::cout << "  delta_mom_rea_pc:   " << this->delta_mom_rea_pc_ << std::endl;
        std::cout << "  delta_mom_rot_n:    " << this->delta_mom_rot_n_ << std::endl;
        std::cout << "  delta_mom_rot_c:    " << this->delta_mom_rot_c_ << std::endl;
        std::cout << "  delta_mom_pre_c :   " << this->delta_mom_pre_c_ << std::endl;
        std::cout << "  delta_mom_pre_n :   " << this->delta_mom_pre_n_ << std::endl;
        std::cout << "  delta_mom_gradiv_n: " << this->delta_mom_graddiv_n_ << std::endl;
        std::cout << "  delta_mom_gradiv_c: " << this->delta_mom_graddiv_c_ << std::endl;
        std::cout << "  delta_inc_n:        " << this->delta_inc_n_ << std::endl;
        std::cout << "  delta_inc_c:        " << this->delta_inc_c_ << std::endl;
    }
}

template<int DIM, class DataType>
DataType MetFlowIncompAssembler<DIM, DataType>::get_coeff ( std::string term )
{
    DataType ret = 9999.;
    DataType skew_fac = 1.;
    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;
    DataType dt = this->dT_pc_;

    if ( term == "p_dt" ) ret = this->theta_d_dt_u_;
    if ( term == "p_mom_vis_c" ) ret = this->theta_mom_vis_c_ * this->nu_ * dt;
    if ( term == "p_mom_adv_cc" ) ret = skew_fac * this->theta_mom_adv_cc_ * dt;
    if ( term == "p_mom_adv_pc" ) ret = skew_fac * this->theta_mom_adv_pc_ * dt;
    if ( term == "p_mom_adv_cp" ) ret = skew_fac * this->theta_mom_adv_cp_ * dt;
    if ( term == "p_mom_adv_pp" ) ret = skew_fac * this->theta_mom_adv_pp_ * dt;
    if ( term == "p_mom_pre" ) ret = -this->theta_mom_pre_c_ * this->inv_rho_ * dt;
    if ( term == "p_inc_c" ) ret = this->theta_inc_c_ * dt;

    if ( this->galerkin_mode_ == 0 )
        dt = this->dT_cn_;

    if ( term == "d_dt" ) ret = this->delta_d_dt_u_;
    if ( term == "d_mom_vis_c" ) ret = this->delta_mom_vis_c_ * this->nu_ * dt;
    if ( term == "d_mom_adv_cc" ) ret = -skew_fac * this->delta_mom_adv_cc_ * dt;
    if ( term == "d_mom_adv_nc" ) ret = -skew_fac * this->delta_mom_adv_nc_ * dt;
    if ( term == "d_mom_adv_cn" ) ret = -skew_fac * this->delta_mom_adv_cn_ * dt;
    if ( term == "d_mom_adv_nn" ) ret = -skew_fac * this->delta_mom_adv_nn_ * dt;
    if ( term == "d_mom_adv_pc" ) ret = -skew_fac * this->delta_mom_adv_pc_ * dt;
    if ( term == "d_mom_pre_c" ) ret = this->delta_mom_pre_c_ * dt;
    if ( term == "d_inc_c" ) ret = -this->delta_inc_c_ * dt * this->inv_rho_;

    if ( this->theta_d_dt_u_ == 0. )
        ret /= dt;

    return ret;
}

template<int DIM, class DataType>
void MetFlowIncompAssembler<DIM, DataType>::compute_tau_graddiv ( const Element<DataType>& element, DataType& tau ) const
{
    if ( this->graddiv_mode_ == CONST )
    {
        tau = this->graddiv_gamma_;
        return;
    }
    DataType dt = this->dT_pc_;
    if ( this->galerkin_mode_ == 0 && this->mode_ == DUAL )
        dt = this->dT_cn_;

    const int num_q = this->num_quadrature_points ( );
    DataType sign = 1.0;

    DataType u_square = 0.0;
    DataType u_max = 0.0;
    DataType vol = 0.0;

    // loop over quadrature points to compute vel^2
    for ( int q = 0; q < num_q; ++q )
    {

        const DataType r = this->x ( q )[1];
        const DataType phi = this->x ( q )[0];
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::abs ( this->detJ ( q ) );

        DataType sin = std::sin ( phi );
        DataType cos = std::cos ( phi );

        DataType u2 = ( this->solP_prev_[0][q] * sin + this->solP_prev_[1][q] * cos ) * ( this->solP_prev_[0][q] * sin + this->solP_prev_[1][q] * cos )
                + ( this->solP_prev_[0][q] * cos + this->solP_prev_[1][q] * sin ) * ( this->solP_prev_[0][q] * cos + this->solP_prev_[1][q] * sin );
        if ( DIM == 3 ) u2 += this->solP_prev_[2][q] * this->solP_prev_[2][q];

        u_square += wq * u2 * dJ;
        vol += wq * dJ;

        if ( u2 > u_max ) u_max = u2;
    }

    u_max = std::sqrt ( u_max );
    DataType h = pow ( vol, 1. / ( DataType ) DIM );
    DataType h2 = h*h;

    if ( this->graddiv_mode_ == VMS )
    {
        DataType tau_M = 4.0 / ( dt * dt ) + u_square / h2 + this->graddiv_C_ * this->nu_ * this->nu_ / ( h2 * h2 );
        tau = this->graddiv_gamma_ * h2 * std::sqrt ( tau_M );
    }
    if ( this->graddiv_mode_ == SUPG )
    {
        DataType tau_M = 1.0 / dt + this->nu_ / h2 + u_max / h;
        tau_M = this->graddiv_gamma_ / tau_M;
        tau = u_max * u_max * tau_M;
    }
}

template class MetFlowIncompAssembler<2, double>;
template class MetFlowIncompAssembler<3, double>;
