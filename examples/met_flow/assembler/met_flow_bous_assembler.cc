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

#include "met_flow_bous_assembler.h"

// convert cartesian to cylindrical coordinates

std::vector<double> cartesian2cyl ( std::vector<double> pt )
{
    std::vector<double> pos_cyl;
    pos_cyl.resize ( pt.size ( ), 0.0 );

    double pi = 3.14159265;

    // r
    pos_cyl[1] = std::sqrt ( pt[0] * pt[0] + pt[1] * pt[1] );

    // phi
    if ( pt[0] == 0.0 )
    {
        if ( pt[1] > 0.0 ) pos_cyl[0] = 0.5 * pi;
        if ( pt[1] < 0.0 ) pos_cyl[0] = 1.5 * pi;
    }
    if ( pt[0] > 0.0 )
    {
        if ( pt[1] >= 0.0 ) pos_cyl[0] = std::atan ( std::abs ( pt[1] ) / std::abs ( pt[0] ) );
        if ( pt[1] < 0.0 ) pos_cyl[0] = 2.0 * pi - std::atan ( std::abs ( pt[1] ) / std::abs ( pt[0] ) );
    }
    if ( pt[0] < 0.0 )
    {
        if ( pt[1] >= 0.0 ) pos_cyl[0] = pi - std::atan ( std::abs ( pt[1] ) / std::abs ( pt[0] ) );
        if ( pt[1] < 0.0 ) pos_cyl[0] = pi + std::atan ( std::abs ( pt[1] ) / std::abs ( pt[0] ) );
    }
    // z
    if ( DIM == 3 ) pos_cyl[2] = pt[2];

    return pos_cyl;
}

template<int DIM, class DataType>
MetFlowBousAssembler<DIM, DataType>::MetFlowBousAssembler ( )
: MetFlowIncompAssembler<DIM, DataType>( )
{
    this->kappa_ = -999999999.;
    this->alpha_g_ = -999999999.;
    this->ref_T_ = -999999999.;
    this->temp_supg_mode_ = TS_OFF;
    this->temp_supg_gamma_ = 0.;
    this->temp_supg_comp_.resize ( DIM, 0. );
    this->g_.resize ( DIM, 0.0 );

#ifdef AUGMENT_PRESS
    this->num_var_ = DIM + 3;
#else
    this->num_var_ = DIM + 2;
#endif

    this->L2_perturb_.resize ( this->num_var_, false );
    this->H1_perturb_.resize ( this->num_var_, false );

    this->set_time_discretization_to_zero ( );
    MetFlowAssembler<DIM, DataType>::allocate_function_values ( this->num_var_ );
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::clear ( )
{
    MetFlowIncompAssembler<DIM, DataType>::clear ( );

    this->kappa_ = -999999999.;
    this->alpha_g_ = -999999999.;
    this->ref_T_ = -999999999.;
    this->temp_supg_mode_ = TS_OFF;
    this->temp_supg_gamma_ = 0.;
    this->temp_supg_comp_ .resize ( DIM, 0. );
    this->g_ .resize ( DIM, 0.0 );

#ifdef AUGMENT_PRESS
    this->num_var_ = DIM + 3;
#else
    this->num_var_ = DIM + 2;
#endif

    this->L2_perturb_.resize ( this->num_var_, false );
    this->H1_perturb_.resize ( this->num_var_, false );

    this->set_time_discretization_to_zero ( );
    MetFlowAssembler<DIM, DataType>::allocate_function_values ( this->num_var_ );
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::set_time_discretization_to_zero ( )
{
    MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_zero ( );
    this->theta_mom_buo_c_ = 0.;
    this->theta_mom_buo_p_ = 0.;
    this->theta_heat_diff_c_ = 0.;
    this->theta_heat_diff_p_ = 0.;
    this->theta_heat_adv_cc_ = 0.;
    this->theta_heat_adv_pc_ = 0.;
    this->theta_heat_adv_cp_ = 0.;
    this->theta_heat_adv_pp_ = 0.;
    this->theta_heat_supg_c_ = 0.;
    this->theta_heat_supg_p_ = 0.;

    this->delta_mom_sou_cc_ = 0.;
    this->delta_mom_sou_nc_ = 0.;
    this->delta_mom_sou_cn_ = 0.;
    this->delta_mom_sou_nn_ = 0.;
    this->delta_mom_sou_pc_ = 0.;
    this->delta_heat_diff_c_ = 0.;
    this->delta_heat_diff_n_ = 0.;
    this->delta_heat_adv_cc_ = 0.;
    this->delta_heat_adv_nc_ = 0.;
    this->delta_heat_adv_cn_ = 0.;
    this->delta_heat_adv_nn_ = 0.;
    this->delta_heat_adv_pc_ = 0.;
    this->delta_heat_rea_cc_ = 0.;
    this->delta_heat_rea_nc_ = 0.;
    this->delta_heat_rea_cn_ = 0.;
    this->delta_heat_rea_nn_ = 0.;
    this->delta_heat_rea_pc_ = 0.;
    this->delta_heat_sou_c_ = 0.;
    this->delta_heat_sou_n_ = 0.;
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::set_time_discretization_to_simple ( )
{
    MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_simple ( );

    this->delta_mom_sou_cc_ = 1. / 3.;
    this->delta_mom_sou_nc_ = 1. / 6.;
    this->delta_mom_sou_cn_ = 1. / 6.;
    this->delta_mom_sou_nn_ = 1. / 3.;

    this->delta_heat_diff_c_ = 0.5;
    this->delta_heat_diff_n_ = 0.5;
    this->delta_heat_adv_cc_ = 1. / 3.;
    this->delta_heat_adv_nc_ = 1. / 6.;
    this->delta_heat_adv_cn_ = 1. / 6.;
    this->delta_heat_adv_nn_ = 1. / 3.;
    this->delta_heat_rea_cc_ = 1. / 3.;
    this->delta_heat_rea_nc_ = 1. / 6.;
    this->delta_heat_rea_cn_ = 1. / 6.;
    this->delta_heat_rea_nn_ = 1. / 3.;
    this->delta_heat_sou_c_ = 0.5;
    this->delta_heat_sou_n_ = 0.5;
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::set_time_discretization_to_theta ( DataType theta )
{
    MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_theta ( theta );

    if ( this->mode_ == PRIMAL )
    {
        this->theta_mom_buo_c_ = theta;
        this->theta_mom_buo_p_ = 1. - theta;
        this->theta_heat_diff_c_ = theta;
        this->theta_heat_diff_p_ = 1 - theta;
        this->theta_heat_adv_cc_ = theta;
        this->theta_heat_adv_pc_ = 0.;
        this->theta_heat_adv_cp_ = 0.;
        this->theta_heat_adv_pp_ = 1 - theta;
        this->theta_heat_supg_c_ = theta;
        this->theta_heat_supg_p_ = 1 - theta;
    }
    if ( this->mode_ == DUAL )
    {
        this->delta_mom_sou_cc_ = 1. - theta;
        //      this->delta_mom_sou_cc_ = 0.;
        this->delta_mom_sou_nc_ = 0.;
        this->delta_mom_sou_cn_ = 0.;
        this->delta_mom_sou_nn_ = theta;
        //      this->delta_mom_sou_nn_ = 0.0;
        this->delta_mom_sou_pc_ = 0.;

        this->delta_heat_diff_c_ = 1. - theta;
        this->delta_heat_diff_n_ = theta;

        this->delta_heat_adv_cc_ = 1. - theta;
        this->delta_heat_adv_nc_ = 0.;
        this->delta_heat_adv_cn_ = 0.;
        this->delta_heat_adv_nn_ = theta;
        this->delta_heat_adv_pc_ = 0.;

        this->delta_heat_rea_cc_ = 1. - theta;
        this->delta_heat_rea_nc_ = 0.;
        this->delta_heat_rea_cn_ = 0.;
        this->delta_heat_rea_nn_ = theta;
        this->delta_heat_rea_pc_ = 0.;

        this->delta_heat_sou_c_ = 1. - theta;
        this->delta_heat_sou_n_ = theta;
    }
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::set_time_discretization_to_galerkin_cd ( bool modified )
{
    MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_galerkin_cd ( modified );

    // Primal Mode
    if ( this->mode_ == PRIMAL )
    {
        this->theta_mom_buo_c_ = 0.5;
        this->theta_mom_buo_p_ = 0.5;
        this->theta_heat_diff_c_ = 0.5;
        this->theta_heat_diff_p_ = 0.5;
        this->theta_heat_adv_cc_ = 1. / 3.;
        this->theta_heat_adv_pc_ = 1. / 6.;
        this->theta_heat_adv_cp_ = 1. / 6.;
        this->theta_heat_adv_pp_ = 1. / 3.;
        this->theta_heat_supg_c_ = 0.5;
        this->theta_heat_supg_p_ = 0.5;
    }

    // Dual Mode
    if ( this->mode_ == DUAL )
    {
        this->delta_mom_sou_cc_ = 1. / 3.;
        this->delta_mom_sou_nc_ = 1. / 6.;
        this->delta_mom_sou_cn_ = 1. / 6.;
        this->delta_mom_sou_nn_ = 1. / 3.;
        this->delta_mom_sou_pc_ = 0.;

        this->delta_heat_diff_c_ = 0.5;
        this->delta_heat_diff_n_ = 0.5;

        this->delta_heat_adv_cc_ = 1. / 3.;
        this->delta_heat_adv_nc_ = 1. / 6.;
        this->delta_heat_adv_cn_ = 1. / 6.;
        this->delta_heat_adv_nn_ = 1. / 3.;
        this->delta_heat_adv_pc_ = 0.;

        this->delta_heat_rea_cc_ = 1. / 3.;
        this->delta_heat_rea_nc_ = 1. / 6.;
        this->delta_heat_rea_cn_ = 1. / 6.;
        this->delta_heat_rea_nn_ = 1. / 3.;
        this->delta_heat_rea_pc_ = 0.;

        this->delta_heat_sou_c_ = 0.5;
        this->delta_heat_sou_n_ = 0.5;
    }
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::set_time_discretization_to_galerkin_dc ( bool modified )
{
    MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_to_galerkin_dc ( modified );

    if ( this->mode_ == DUAL )
    {
        this->delta_mom_sou_cc_ = 1. / 3.;
        this->delta_mom_sou_nc_ = 0.;
        this->delta_mom_sou_cn_ = 1. / 3.;
        this->delta_mom_sou_nn_ = 1. / 6.;
        this->delta_mom_sou_pc_ = 1. / 6.;

        this->delta_heat_diff_c_ = 0.5;
        this->delta_heat_diff_n_ = 0.5;

        this->delta_heat_adv_cc_ = 1. / 3.;
        this->delta_heat_adv_nc_ = 0.;
        this->delta_heat_adv_cn_ = 1. / 3.;
        this->delta_heat_adv_nn_ = 1. / 6.;
        this->delta_heat_adv_pc_ = 1. / 6.;

        this->delta_heat_rea_cc_ = 1. / 3.;
        this->delta_heat_rea_nc_ = 0.;
        this->delta_heat_rea_cn_ = 1. / 3.;
        this->delta_heat_rea_nn_ = 1. / 6.;
        this->delta_heat_rea_pc_ = 1. / 6.;

        this->delta_heat_sou_c_ = 0.5;
        this->delta_heat_sou_n_ = 0.5;
    }
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::set_time_discretization_off ( )
{
    MetFlowIncompAssembler<DIM, DataType>::set_time_discretization_off ( );
    this->theta_mom_buo_c_ = 1.;
    this->theta_mom_buo_p_ = 0.;
    this->theta_heat_diff_c_ = 1.;
    this->theta_heat_diff_p_ = 0.;
    this->theta_heat_adv_cc_ = 1.;
    this->theta_heat_adv_pc_ = 0.;
    this->theta_heat_adv_cp_ = 0.;
    this->theta_heat_adv_pp_ = 0.;
    this->theta_heat_supg_c_ = 1.;
    this->theta_heat_supg_p_ = 0.;
    this->theta_inc_c_ = 1.;
    this->theta_inc_p_ = 0.;
    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::set_gravity ( std::vector<DataType> g_cart )
{
    this->g_.resize ( DIM, 0.0 );
    this->g_ = cartesian2cyl ( g_cart );
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::print_parameters ( )
{
    MetFlowIncompAssembler<DIM, DataType>::print_parameters ( );

    std::cout.scientific;
    std::cout.precision ( 6 );
    if ( this->mode_ == PRIMAL )
        std::cout << "> MetFlowBousAssembler primal parameters" << std::endl;
    else
        std::cout << "> MetFlowBousAssembler dual parameters" << std::endl;

    std::cout << "  alpha_g:  " << this->alpha_g_ << std::endl;
    std::cout << "  kappa:    " << this->kappa_ << std::endl;
    std::cout << "  ref_T:    " << this->ref_T_ << std::endl;
    std::cout << "  omega:    " << this->omega_ << std::endl;
    std::cout << "  g:        " << this->g_[0] << " " << this->g_[1] << " " << this->g_[2] << std::endl;
    std::cout << "  Temp SUPG mode:    " << this->temp_supg_mode_ << std::endl;
    std::cout << "  Temp SUPG gamma:   " << this->temp_supg_gamma_ << std::endl;
    std::cout << "  Temp SUPG azi:     " << this->temp_supg_comp_[0] << std::endl;
    std::cout << "  Temp SUPG rad:     " << this->temp_supg_comp_[1] << std::endl;
    std::cout << "  Temp SUPG axi:     " << this->temp_supg_comp_[2] << std::endl;
    std::cout << std::endl;

    MetFlowBousAssembler<DIM, DataType>::print_coeff ( );
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::print_coeff ( )
{
    if ( this->mode_ == PRIMAL )
    {
        std::cout << "  theta_mom_buo_c :   " << this->theta_mom_buo_c_ << std::endl;
        std::cout << "  theta_mom_buo_p :   " << this->theta_mom_buo_p_ << std::endl;
        std::cout << "  theta_hea_diff_c:   " << this->theta_heat_diff_c_ << std::endl;
        std::cout << "  theta_hea_diff_p:   " << this->theta_heat_diff_p_ << std::endl;
        std::cout << "  theta_hea_adv_cc:   " << this->theta_heat_adv_cc_ << std::endl;
        std::cout << "  theta_hea_adv_pc:   " << this->theta_heat_adv_pc_ << std::endl;
        std::cout << "  theta_hea_adv_cp:   " << this->theta_heat_adv_cp_ << std::endl;
        std::cout << "  theta_hea_adv_pp:   " << this->theta_heat_adv_pp_ << std::endl;
    }
    else
    {
        std::cout << "  delta_mom_sou_nn:   " << this->delta_mom_sou_nn_ << std::endl;
        std::cout << "  delta_mom_sou_cn:   " << this->delta_mom_sou_cn_ << std::endl;
        std::cout << "  delta_mom_sou_nc:   " << this->delta_mom_sou_nc_ << std::endl;
        std::cout << "  delta_mom_sou_cc:   " << this->delta_mom_sou_cc_ << std::endl;
        std::cout << "  delta_mom_sou_pc:   " << this->delta_mom_sou_pc_ << std::endl;
        std::cout << "  delta_hea_diff_n:   " << this->delta_heat_diff_n_ << std::endl;
        std::cout << "  delta_hea_diff_c:   " << this->delta_heat_diff_c_ << std::endl;
        std::cout << "  delta_hea_adv_nn:   " << this->delta_heat_adv_nn_ << std::endl;
        std::cout << "  delta_hea_adv_cn:   " << this->delta_heat_adv_cn_ << std::endl;
        std::cout << "  delta_hea_adv_nc:   " << this->delta_heat_adv_nc_ << std::endl;
        std::cout << "  delta_hea_adv_cc:   " << this->delta_heat_adv_cc_ << std::endl;
        std::cout << "  delta_hea_adv_pc:   " << this->delta_heat_adv_pc_ << std::endl;
        std::cout << "  delta_hea_rea_nn:   " << this->delta_heat_rea_nn_ << std::endl;
        std::cout << "  delta_hea_rea_cn:   " << this->delta_heat_rea_cn_ << std::endl;
        std::cout << "  delta_hea_rea_nc:   " << this->delta_heat_rea_nc_ << std::endl;
        std::cout << "  delta_hea_rea_cc:   " << this->delta_heat_rea_cc_ << std::endl;
        std::cout << "  delta_hea_rea_pc:   " << this->delta_heat_rea_pc_ << std::endl;
        std::cout << "  delta_hea_sou_n:    " << this->delta_heat_sou_n_ << std::endl;
        std::cout << "  delta_hea_sou_c:    " << this->delta_heat_sou_c_ << std::endl;
    }
}

template<int DIM, class DataType>
DataType MetFlowBousAssembler<DIM, DataType>::get_coeff ( std::string term )
{
    DataType ret = MetFlowIncompAssembler<DIM, DataType>::get_coeff ( term );
    return ret;
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::initialize_for_facet ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, int facet_number )
{
    MetFlowAssembler<DIM, DataType>::initialize_for_facet ( element, quadrature, facet_number );

#ifdef AUGMENT_PRESS
    const int t_var = DIM + 2;
#else
    const int t_var = DIM + 1;
#endif

    // temperature
    this->grad_solP_prev_[t_var].clear ( );
    this->evaluate_fe_function_gradients ( this->vector_solP_prev ( ), t_var, this->grad_solP_prev_[t_var] );
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::initialize_for_element ( const Element<DataType>& element, const Quadrature<DataType>& quadrature )
{
    MetFlowAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );

#ifdef AUGMENT_PRESS
    const int t_var = DIM + 2;
#else
    const int t_var = DIM + 1;
#endif

    if ( this->temp_supg_mode_ > 0 )
    {
        this->evaluate_fe_function_hessians ( this->vector_solP ( ), t_var, this->hess_solP_[t_var] );
        this->evaluate_fe_function_hessians ( this->vector_solP_prev ( ), t_var, this->hess_solP_prev_[t_var] );
    }
}

template<int DIM, class DataType>
void MetFlowBousAssembler<DIM, DataType>::compute_tau_temp_supg ( const Element<DataType>& element, DataType& tau ) const
{
    const int num_q = this->num_quadrature_points ( );

    DataType u_inf = 0.0;
    DataType vol = 0.0;

    // loop over quadrature points to compute ||vel||_inf and volume
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::abs ( this->detJ ( q ) );

        vol += wq * dJ;

        DataType max_val = 0;
        for ( int s = 0; s < DIM; ++s ) max_val += this->solP_prev_[s][q] * this->solP_prev_[s][q];
        max_val = sqrt ( max_val );
        if ( max_val > u_inf ) u_inf = max_val;
    }
    DataType h = pow ( vol, 1. / ( DataType ) DIM );
    DataType pec = h * u_inf / this->kappa_;

    if ( pec > 1. )
    {
        tau = temp_supg_gamma_ * 0.5 * h / u_inf;
    }
    else
    { // no stab required
        tau = 0.;
    }
}

template class MetFlowBousAssembler<2, double>;
template class MetFlowBousAssembler<3, double>;
