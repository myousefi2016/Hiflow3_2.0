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

#include "met_flow_incomp_cart_assembler.h"

template<int DIM, class DataType>
MetFlowIncompCartAssembler<DIM, DataType>::MetFlowIncompCartAssembler ( )
: MetFlowIncompAssembler<DIM, DataType>( )
{
}

/// ********************************************************
/// General Assembly routines 
/// ********************************************************

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    if ( this->mode_ == PRIMAL )
    {
        MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_matrix_primal ( element, lm );
    }
    else
    {
        interminable_assert ( 0 );
    }
}

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    if ( this->vector_asm_mode_ == VECTOR_STD )
    {
        if ( this->mode_ == PRIMAL )
        {
            MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_primal ( element, lv );
        }
        else
        {
            interminable_assert ( 0 );
        }
    }
    else if ( this->vector_asm_mode_ == VECTOR_QUANT )
    {
        lv.clear ( );
        lv.resize ( this->num_scalars_, 0. );
        MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_quant ( element, lv );
    }
    else if ( this->vector_asm_mode_ == VECTOR_VORT )
    {
        MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_vort ( element, lv );
    }
    else if ( this->vector_asm_mode_ == VECTOR_GOAL )
    {
        MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_goal ( element, lv );
    }
}

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const
{
    ls = 0.;
    if ( this->sca_mode_ == DIV_MEAN )
    {
        MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_div_mean ( element, ls );
    }
    if ( this->sca_mode_ == DIV_MAX )
    {
        MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_div_max ( element, ls );
    }
    if ( this->sca_mode_ == GOAL_INT )
    {
        MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_goal_int ( element, ls );
    }
    if ( this->sca_mode_ == GOAL_FIN )
    {
        MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_goal_fin ( element, ls );
    }
}

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const
{
    lv[facet_number] = 0.;
}

/// ********************************************************
/// Assembly routines for primal problem
/// ********************************************************

/// Jacobian of primal problem 
#ifdef OPT_ASM0

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int p_var = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
#    endif

    const DataType dt = this->dT_pc_;
    DataType sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    DataType tau_div = 0.;
    DataType skew_fac = 1.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCartAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::fabs ( this->detJ ( q ) );

        // ***********************************************************************
        // TIME DIFFERENCE

        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq
                        * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                        * this->phi ( j, q, u_var )
                    * this->phi ( i, q, u_var )
                    * dJ;

            }
        }

        // ***********************************************************************
        // LAPLACE: theta_mom_vis_c * dT * nu * int{grad{u}:grad{v}}
        // symmetric terms
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    DataType tmp = 0.;
                    for ( int s = 0; s < DIM; ++s )
                    { // grad(u) : grad(v)
                        tmp += this->grad_phi ( j, q, u_var )[s]
                                * this->grad_phi ( i, q, u_var )[s];
                    }
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq
                            * dt
                            * this->theta_mom_vis_c_
                            * this->nu_
                            * tmp
                            * dJ;
                }
            }
        }

        // ***********************************************************************
        // CONVECTION TERM
        // nonlinear l2k part one and two: (0.5 *) theta_mom_adv_cc * dT * int{ (u_c*\grad{u})*v } + (0.5 *) theta_mom_adv_pc * dT * int{ (u_p*\grad{u})*v }
        if ( this->conv_mode_ == OSEEN || this->conv_mode_ == NAVIERSTOKES )
        {
            std::vector<DataType> convection ( DIM, 0. );
            std::vector<DataType> convection_prev ( DIM, 0. );
            if ( this->conv_mode_ == OSEEN )
            {
                for ( int s = 0; s < DIM; ++s )
                {
                    convection[s] = this->conv_[s][q];
                    convection_prev[s] = this->conv_prev_[s][q];
                }
            }
            if ( this->conv_mode_ == NAVIERSTOKES )
            {
                for ( int s = 0; s < DIM; ++s )
                {
                    convection[s] = this->solP_[s][q];
                    convection_prev[s] = this->solP_prev_[s][q];
                }
            }

            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    DataType tmp = 0;
                    for ( int s = 0; s < DIM; ++s ) // scalar product
                        tmp += ( this->theta_mom_adv_cc_ * convection[s] + this->theta_mom_adv_pc_ * convection_prev[s] ) * this->grad_phi ( j, q, u_var )[s];

                    for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    {
                        lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq
                                * skew_fac
                                * dt
                                * tmp
                                * this->phi ( i, q, u_var )
                                * dJ;
                    }
                }
            }
        } // if OSEEN || NAVIERSTOKES

        if ( this->conv_mode_ == NAVIERSTOKES )
        {
            // nonlinear l2k part three and four : (0.5 *) theta_mom_adv_cc_ * dT * int{ (u\grad{u_c}*v } + (0.5 *) theta_mom_adv_cp_ * dT * int{ (u\grad{u_p}*v }
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType tmp = 0.0;
                            tmp = this->phi ( j, q, u_var ) * this->phi ( i, q, v_var )
                                    * ( this->theta_mom_adv_cc_ * this->grad_solP_[v_var][q][u_var] + this->theta_mom_adv_cp_ * this->grad_solP_prev_[v_var][q][u_var] );

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq
                                    * skew_fac
                                    * dt
                                    * tmp
                                    * dJ;
                        }
                    }
                }
            }

            // Skew symmetric form    
            if ( this->skew_mode_ > 0 )
            {
                // nonlinear l2k part one and two: -0.5 * theta_mom_adv_cc * dT * int{ (u*\grad{v})*u_c } - 0.5 * theta_mom_adv_cp * dT * int{ (u*\grad{v})*u_p }

                for ( int v_var = 0; v_var < DIM; ++v_var )
                {
                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                        {
                            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                            {
                                DataType tmp = this->phi ( j, q, u_var ) * this->grad_phi ( i, q, v_var )[u_var];

                                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += -wq
                                        * skew_fac
                                        * dt
                                        * tmp
                                        * ( this->theta_mom_adv_cc_ * this->solP_[v_var][q] + this->theta_mom_adv_cp_ * this->solP_prev_[v_var][q] )
                                        * dJ;
                            }
                        }
                    }
                }

                // nonlinear l2k part three and four: -0.5 * theta_mom_adv_cc * dT * int{ (u_c*\grad{v})*u } - 0.5 * theta_mom_adv_pc * dT * int{ (u_p*\grad{v})*u }
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType tmp = 0.;
                            for ( int s = 0; s < DIM; s++ )
                                tmp += ( this->theta_mom_adv_cc_ * this->solP_[s][q] + this->theta_mom_adv_pc_ * this->solP_prev_[s][q] )
                                * this->grad_phi ( i, q, u_var )[s];

                            tmp *= this->phi ( j, q, u_var );
                            lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq
                                    * skew_fac
                                    * dt
                                    * tmp
                                    * dJ;
                        }
                    }
                }
            }
        } // if NAVIERSTOKES
#    ifdef ROTATING_FOR
        // Coriolis force:  -dT * theta1 * 2 * omega (u[0] * v[1] - u[1] * v[0])
        // TODO
        /*
        int u_var = 0;
        int v_var = 1;
        for (int j = 0; j < this->num_dofs(u_var); ++j) 
        {
            for (int i = 0; i < this->num_dofs(v_var); ++i) 
            {
                lm(this->dof_index(i, v_var), this->dof_index(j, u_var)) += -wq
         * this->theta_mom_rot_c_
         * dt 
         * 2.0
         * this->omega_
         * this->phi(j,q,u_var)
         * this->phi(i,q,v_var)
         * dJ;
            }
        }

        u_var = 1;
        v_var = 0;
        for (int j = 0; j < this->num_dofs(u_var); ++j) 
        {
            for (int i = 0; i < this->num_dofs(v_var); ++i) 
            {
                lm(this->dof_index(i, v_var), this->dof_index(j, u_var)) += wq
         * this->theta_mom_rot_c_
         * dt 
         * 2.0
         * this->omega_
         * this->phi(j,q,u_var)
         * this->phi(i,q,v_var)
         * dJ;
            }
        }
         * */
#    endif  
        // ***********************************************************************
        // PRESSURE: -inv_rho * theta_mpm_pre_ * dT * \int{p div{v}}

        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                DataType div_v = 0.0;

                if ( v_var == 0 ) div_v = this->grad_phi ( i, q, 0 )[0];
                if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1];
                if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                int ind_i = this->dof_index ( i, v_var );

                for ( int j = 0; j<this->num_dofs ( p_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p_var ) ) += -wq
                            * dt
                            * this->theta_mom_pre_c_
                            * this->inv_rho_
                            * this->phi ( j, q, p_var )
                            * div_v
                            * dJ;
                }
#    ifdef AUGMENT_PRESS
                for ( int j = 0; j<this->num_dofs ( p0_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p0_var ) ) += -wq
                            * dt
                            * this->theta_mom_pre_c_
                            * this->inv_rho_
                            * this->phi ( j, q, p0_var )
                            * div_v
                            * dJ;
                }
#    endif            
            }
        }

        // ***********************************************************************
        // Grad-div stabilization: theta1 * dT * gamma_div * (div{u}, div{v})
        if ( this->graddiv_mode_ > 0 )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    DataType div_v = 0.0;
                    if ( v_var == 0 ) div_v = this->grad_phi ( i, q, 0 )[0];
                    if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1];
                    if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType div_u = 0.0;
                            if ( u_var == 0 ) div_u = this->grad_phi ( j, q, 0 )[0];
                            if ( u_var == 1 ) div_u = this->grad_phi ( j, q, 1 )[1];
                            if ( u_var == 2 ) div_u = this->grad_phi ( j, q, 2 )[2];

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq
                                    * tau_div
                                    * this->theta_mom_graddiv_c_
                                    * dt
                                    * div_u
                                    * div_v
                                    * dJ;
                        }
                    }
                }
            }
        }

        // ***********************************************************************
        // constraint: theta_inc_c * dT * \int{q div(u)}
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
            {
                DataType div_u = 0.0;

                if ( u_var == 0 ) div_u = this->grad_phi ( j, q, 0 )[0];
                if ( u_var == 1 ) div_u = this->grad_phi ( j, q, 1 )[1];
                if ( u_var == 2 ) div_u = this->grad_phi ( j, q, 2 )[2];

                int ind_j = this->dof_index ( j, u_var );

                for ( int i = 0; i<this->num_dofs ( p_var ); ++i )
                {
                    lm ( this->dof_index ( i, p_var ), ind_j ) += wq
                            * dt
                            * this->theta_inc_c_
                            * this->phi ( i, q, p_var )
                            * div_u
                            * dJ;
                }
#    ifdef AUGMENT_PRESS
                for ( int i = 0; i<this->num_dofs ( p0_var ); ++i )
                {
                    lm ( this->dof_index ( i, p0_var ), ind_j ) += wq
                            * dt
                            * this->theta_inc_c_
                            * this->phi ( i, q, p0_var )
                            * div_u
                            * dJ;
                }
#    endif
            }
        }
    }
}
#endif
#ifdef OPT_ASM1

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int p_var = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
#    endif    
    const DataType dt = this->dT_pc_;

    DataType sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    DataType tau_div = 0.;
    DataType skew_fac = 1.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCartAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    const DataType mom_vis_c = dt * this->theta_mom_vis_c_ * this->nu_;
    const DataType mom_adv_cc = this->theta_mom_adv_cc_ * skew_fac * dt;
    const DataType mom_adv_pc = this->theta_mom_adv_pc_ * skew_fac * dt;
    const DataType mom_adv_cp = this->theta_mom_adv_cp_ * skew_fac * dt;
    const DataType mom_rot_c = this->theta_mom_rot_c_ * dt * 2.0 * this->omega_;
    const DataType mom_pre = dt * this->theta_mom_pre_c_ * this->inv_rho_;
    const DataType mom_graddiv = tau_div * this->theta_mom_graddiv_c_ * dt;
    const DataType inc_c = dt * this->theta_inc_c_;

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq_dJ = this->w ( q ) * std::fabs ( this->detJ ( q ) );

        // ***********************************************************************
        // TIME DIFFERENCE

        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                            * this->phi ( j, q, u_var )
                            * this->phi ( i, q, u_var );

                }
            }
        }

        // ***********************************************************************
        // LAPLACE: theta_mom_vis_c * dT * nu * int{grad{u}:grad{v}}
        // symmetric terms
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    DataType tmp = 0.;
                    for ( int s = 0; s < DIM; ++s )
                    { // grad(u) : grad(v)
                        tmp += this->grad_phi ( j, q, u_var )[s]
                                * this->grad_phi ( i, q, u_var )[s];
                    }
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * mom_vis_c
                            * tmp;
                }
            }
        }

        // ***********************************************************************
        // CONVECTION TERM
        // nonlinear l2k part one and two: (0.5 *) theta_mom_adv_cc * dT * int{ (u_c*\grad{u})*v } + (0.5 *) theta_mom_adv_pc * dT * int{ (u_p*\grad{u})*v }
        if ( this->conv_mode_ == OSEEN || this->conv_mode_ == NAVIERSTOKES )
        {
            std::vector<DataType> convection ( DIM, 0. );
            std::vector<DataType> convection_prev ( DIM, 0. );
            if ( this->conv_mode_ == OSEEN )
            {
                for ( int s = 0; s < DIM; ++s )
                {
                    convection[s] = this->conv_[s][q];
                    convection_prev[s] = this->conv_prev_[s][q];
                }
            }
            if ( this->conv_mode_ == NAVIERSTOKES )
            {
                for ( int s = 0; s < DIM; ++s )
                {
                    convection[s] = this->solP_[s][q];
                    convection_prev[s] = this->solP_prev_[s][q];
                }
            }

            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    DataType tmp = 0;
                    for ( int s = 0; s < DIM; ++s ) // scalar product
                        tmp += ( mom_adv_cc * convection[s] + mom_adv_pc * convection_prev[s] ) * this->grad_phi ( j, q, u_var )[s];

                    for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    {
                        lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                * tmp
                                * this->phi ( i, q, u_var );
                    }
                }
            }
        }

        if ( this->conv_mode_ == NAVIERSTOKES )
        {
            // nonlinear l2k part three and four : (0.5 *) theta_mom_adv_cc_ * dT * int{ (u\grad{u_c}*v } + (0.5 *) theta_mom_adv_cp_ * dT * int{ (u\grad{u_p}*v }
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType tmp = 0.0;
                            tmp = this->phi ( j, q, u_var ) * this->phi ( i, q, v_var )
                                    * ( mom_adv_cc * this->grad_solP_[v_var][q][u_var] + mom_adv_cp * this->grad_solP_prev_[v_var][q][u_var] );
                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                    * tmp;
                        }
                    }
                }
            }

            // Skew symmetric form    
            if ( this->skew_mode_ > 0 )
            {
                // nonlinear l2k part one and two: -0.5 * theta_mom_adv_cc * dT * int{ (u*\grad{v})*u_c } - 0.5 * theta_mom_adv_cp * dT * int{ (u*\grad{v})*u_p }

                for ( int v_var = 0; v_var < DIM; ++v_var )
                {
                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                        {
                            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                            {
                                DataType tmp = this->phi ( j, q, u_var ) * this->grad_phi ( i, q, v_var )[u_var];

                                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                                        * tmp
                                        * ( mom_adv_cc * this->solP_[v_var][q] + mom_adv_cp * this->solP_prev_[v_var][q] );
                            }
                        }
                    }
                }

                // nonlinear l2k part three and four: -0.5 * theta_mom_adv_cc * dT * int{ (u_c*\grad{v})*u } - 0.5 * theta_mom_adv_pc * dT * int{ (u_p*\grad{v})*u }
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType tmp = 0.;
                            for ( int s = 0; s < DIM; s++ )
                                tmp += ( mom_adv_cc * this->solP_[s][q] + mom_adv_pc * this->solP_prev_[s][q] )
                                * this->grad_phi ( i, q, u_var )[s];

                            tmp *= this->phi ( j, q, u_var );
                            lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                                    * tmp;
                        }
                    }
                }
            }
        }
#    ifdef ROTATING_FOR
        // Coriolis force:  -dT * theta1 * 2 * omega (u[0] * v[1] - u[1] * v[0])
        // TODO
        /*
            int u_var = 0;
            int v_var = 1;
            for (int j = 0; j < this->num_dofs(u_var); ++j) 
            {
                for (int i = 0; i < this->num_dofs(v_var); ++i) 
                {
                    lm(this->dof_index(i, v_var), this->dof_index(j, u_var)) += -wq_dJ
         * mom_rot_c
         * this->phi(j,q,u_var)
         * this->phi(i,q,v_var);
                }
            }

            u_var = 1;
            v_var = 0;
            for (int j = 0; j < this->num_dofs(u_var); ++j) 
            {
                for (int i = 0; i < this->num_dofs(v_var); ++i) 
                {
                    lm(this->dof_index(i, v_var), this->dof_index(j, u_var)) += wq_dJ
         * mom_rot_c
         * this->phi(j,q,u_var)
         * this->phi(i,q,v_var);
                }
            }
         */
#    endif  
        // ***********************************************************************
        // PRESSURE: -inv_rho * theta_mpm_pre_ * dT * \int{p div{v}} 
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                DataType div_v = 0.0;

                if ( v_var == 0 ) div_v = this->grad_phi ( i, q, 0 )[0];
                if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1];
                if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                int ind_i = this->dof_index ( i, v_var );

                for ( int j = 0; j<this->num_dofs ( p_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p_var ) ) += -wq_dJ
                            * mom_pre
                            * this->phi ( j, q, p_var )
                            * div_v;
                }
#    ifdef AUGMENT_PRESS
                for ( int j = 0; j<this->num_dofs ( p0_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p0_var ) ) += -wq_dJ
                            * mom_pre
                            * this->phi ( j, q, p0_var )
                            * div_v;
                }
#    endif    
            }
        }

        // ***********************************************************************
        // Grad-div stabilization: theta1 * dT * gamma_div * (div{u}, div{v})
        if ( this->graddiv_mode_ > 0 )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    DataType div_v = 0.0;
                    if ( v_var == 0 ) div_v = this->grad_phi ( i, q, 0 )[0];
                    if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1];
                    if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType div_u = 0.0;
                            if ( u_var == 0 ) div_u = this->grad_phi ( j, q, 0 )[0];
                            if ( u_var == 1 ) div_u = this->grad_phi ( j, q, 1 )[1];
                            if ( u_var == 2 ) div_u = this->grad_phi ( j, q, 2 )[2];

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                    * mom_graddiv
                                    * div_u
                                    * div_v;
                        }
                    }
                }
            }
        }

        // ***********************************************************************
        // constraint: theta_inc_c * dT * \int{q div(u)}
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
            {
                DataType div_u = 0.0;

                if ( u_var == 0 ) div_u = this->grad_phi ( j, q, 0 )[0];
                if ( u_var == 1 ) div_u = this->grad_phi ( j, q, 1 )[1];
                if ( u_var == 2 ) div_u = this->grad_phi ( j, q, 2 )[2];

                int ind_j = this->dof_index ( j, u_var );

                for ( int i = 0; i<this->num_dofs ( p_var ); ++i )
                {
                    lm ( this->dof_index ( i, p_var ), ind_j ) += wq_dJ
                            * inc_c
                            * this->phi ( i, q, p_var )
                            * div_u;
                }
#    ifdef AUGMENT_PRESS
                for ( int i = 0; i<this->num_dofs ( p0_var ); ++i )
                {
                    lm ( this->dof_index ( i, p0_var ), ind_j ) += wq_dJ
                            * inc_c
                            * this->phi ( i, q, p0_var )
                            * div_u;
                }
#    endif
            }
        }
    }
}
#endif
#ifdef OPT_ASM2

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int p_var = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#    else
    const int num_vars = DIM + 1;
#    endif    
    const int num_std_vars = DIM + 1;

    const DataType dt = this->dT_pc_;

    DataType sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    DataType tau_div = 0.;
    DataType skew_fac = 1.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCartAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    const DataType mom_vis_c = dt * this->theta_mom_vis_c_ * this->nu_;
    const DataType mom_adv_cc = this->theta_mom_adv_cc_ * skew_fac * dt;
    const DataType mom_adv_pc = this->theta_mom_adv_pc_ * skew_fac * dt;
    const DataType mom_adv_cp = this->theta_mom_adv_cp_ * skew_fac * dt;
    const DataType mom_rot_c = this->theta_mom_rot_c_ * dt * 2.0 * this->omega_;
    const DataType mom_pre = dt * this->theta_mom_pre_c_ * this->inv_rho_;
    const DataType mom_graddiv = tau_div * this->theta_mom_graddiv_c_ * dt;
    const DataType inc_c = dt * this->theta_inc_c_;

    std::vector< std::vector<DataType> > phi;
    std::vector< std::vector< std::vector<DataType> > > grad_phi;
    std::vector< std::vector< Mat<DIM, DIM, DataType> > > H_phi;
    phi.resize ( num_vars );
    grad_phi.resize ( num_std_vars );

#    ifdef AUGMENT_PRESS    
    phi[p0_var].resize ( this->num_dofs ( p0_var ), 0. );
#    endif

    for ( int k = 0; k < num_std_vars; ++k )
    {
        phi[k].resize ( this->num_dofs ( k ), 0. );
        grad_phi[k].resize ( this->num_dofs ( k ) );
        for ( int j = 0; j<this->num_dofs ( k ); ++j )
        {
            grad_phi[k][j].resize ( DIM, 0. );
        }
    }

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq_dJ = this->w ( q ) * std::fabs ( this->detJ ( q ) );

        for ( int k = 0; k < num_std_vars; ++k )
        {
            for ( int j = 0; j<this->num_dofs ( k ); ++j )
            {
                phi[k][j] = this->phi ( j, q, k );
                for ( int d = 0; d < DIM; ++d )
                {
                    grad_phi[k][j][d] = this->grad_phi ( j, q, k )[d];
                }
            }
        }

#    ifdef AUGMENT_PRESS
        for ( int j = 0; j<this->num_dofs ( p0_var ); ++j )
        {
            phi[p0_var][j] = this->phi ( j, q, p0_var );
        }
#    endif

        // ***********************************************************************
        // TIME DIFFERENCE
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                        * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                        * phi[u_var][j]
                        * phi[u_var][i];

            }
        }

        // ***********************************************************************
        // LAPLACE: theta_mom_vis_c * dT * nu * int{grad{u}:grad{v}}
        // symmetric terms
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    DataType tmp = 0.;
                    for ( int s = 0; s < DIM; ++s )
                    { // grad(u) : grad(v)
                        tmp += grad_phi[u_var][j][s]
                                * grad_phi[u_var][i][s];
                    }
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * mom_vis_c
                            * tmp;
                }
            }
        }

        // ***********************************************************************
        // CONVECTION TERM
        // nonlinear l2k part one and two: (0.5 *) theta_mom_adv_cc * dT * int{ (u_c*\grad{u})*v } + (0.5 *) theta_mom_adv_pc * dT * int{ (u_p*\grad{u})*v }
        if ( this->conv_mode_ == OSEEN || this->conv_mode_ == NAVIERSTOKES )
        {
            std::vector<DataType> convection ( DIM, 0. );
            std::vector<DataType> convection_prev ( DIM, 0. );
            if ( this->conv_mode_ == OSEEN )
            {
                for ( int s = 0; s < DIM; ++s )
                {
                    convection[s] = this->conv_[s][q];
                    convection_prev[s] = this->conv_prev_[s][q];
                }
            }
            if ( this->conv_mode_ == NAVIERSTOKES )
            {
                for ( int s = 0; s < DIM; ++s )
                {
                    convection[s] = this->solP_[s][q];
                    convection_prev[s] = this->solP_prev_[s][q];
                }
            }

            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    DataType tmp = 0;
                    for ( int s = 0; s < DIM; ++s ) // scalar product
                        tmp += ( mom_adv_cc * convection[s] + mom_adv_pc * convection_prev[s] ) * grad_phi[u_var][j][s];

                    for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    {
                        lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                * tmp
                                * phi[u_var][i];
                    }
                }
            }
        }
        if ( this->conv_mode_ == NAVIERSTOKES )
        {
            // nonlinear l2k part three and four : (0.5 *) theta_mom_adv_cc_ * dT * int{ (u\grad{u_c}*v } + (0.5 *) theta_mom_adv_cp_ * dT * int{ (u\grad{u_p}*v }
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType tmp = 0.0;
                            tmp = phi[u_var][j] * phi[v_var][i]
                                    * ( mom_adv_cc * this->grad_solP_[v_var][q][u_var] + mom_adv_cp * this->grad_solP_prev_[v_var][q][u_var] );

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                    * tmp;
                        }
                    }
                }
            }

            // Skew symmetric form    
            if ( this->skew_mode_ > 0 )
            {
                // nonlinear l2k part one and two: -0.5 * theta_mom_adv_cc * dT * int{ (u*\grad{v})*u_c } - 0.5 * theta_mom_adv_cp * dT * int{ (u*\grad{v})*u_p }

                for ( int v_var = 0; v_var < DIM; ++v_var )
                {
                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                        {
                            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                            {
                                DataType tmp = phi[u_var][j] * grad_phi[v_var][i][u_var];

                                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                                        * tmp
                                        * ( mom_adv_cc * this->solP_[v_var][q] + mom_adv_cp * this->solP_prev_[v_var][q] );
                            }
                        }
                    }
                }

                // nonlinear l2k part three and four: -0.5 * theta_mom_adv_cc * dT * int{ (u_c*\grad{v})*u } - 0.5 * theta_mom_adv_pc * dT * int{ (u_p*\grad{v})*u }
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType tmp = 0.;
                            for ( int s = 0; s < DIM; s++ )
                                tmp += ( mom_adv_cc * this->solP_[s][q] + mom_adv_pc * this->solP_prev_[s][q] )
                                * grad_phi[u_var][i][s];

                            tmp *= phi[u_var][j];
                            lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                                    * tmp;
                        }
                    }
                }
            }
        }
#    ifdef ROTATING_FOR
        // Coriolis force:  -dT * theta1 * 2 * omega (u[0] * v[1] - u[1] * v[0])
        // TODO
        /*
            int u_var = 0;
            int v_var = 1;
            for (int j = 0; j < this->num_dofs(u_var); ++j) 
            {
                for (int i = 0; i < this->num_dofs(v_var); ++i) 
                {
                    lm(this->dof_index(i, v_var), this->dof_index(j, u_var)) += -wq_dJ
         * mom_rot_c
         * phi[u_var][j]
         * phi[v_var][i];
                }
            }

            u_var = 1;
            v_var = 0;
            for (int j = 0; j < this->num_dofs(u_var); ++j)
            {
                for (int i = 0; i < this->num_dofs(v_var); ++i) 
                {
                    lm(this->dof_index(i, v_var), this->dof_index(j, u_var)) += wq_dJ
         * mom_rot_c
         * phi[u_var][j]
         * phi[v_var][i];
                }
            }
         */
#    endif  
        // ***********************************************************************
        // PRESSURE: -inv_rho * theta_mpm_pre_ * dT * \int{p div{v}}
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                DataType div_v = 0.0;

                if ( v_var == 0 ) div_v = grad_phi[0][i][0];
                if ( v_var == 1 ) div_v = grad_phi[1][i][1];
                if ( v_var == 2 ) div_v = grad_phi[2][i][2];

                int ind_i = this->dof_index ( i, v_var );

                for ( int j = 0; j<this->num_dofs ( p_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p_var ) ) += -wq_dJ
                            * mom_pre
                            * phi[p_var][j]
                            * div_v;
                }
#    ifdef AUGMENT_PRESS
                for ( int j = 0; j<this->num_dofs ( p0_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p0_var ) ) += -wq_dJ
                            * mom_pre
                            * phi[p0_var][j]
                            * div_v;
                }
#    endif
            }
        }

        // ***********************************************************************
        // Grad-div stabilization: theta1 * dT * gamma_div * (div{u}, div{v})
        if ( this->graddiv_mode_ > 0 )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    DataType div_v = 0.0;
                    if ( v_var == 0 ) div_v = grad_phi[0][i][0];
                    if ( v_var == 1 ) div_v = grad_phi[1][i][1];
                    if ( v_var == 2 ) div_v = grad_phi[2][i][2];

                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            DataType div_u = 0.0;
                            if ( u_var == 0 ) div_u = grad_phi[0][j][0];
                            if ( u_var == 1 ) div_u = grad_phi[1][j][1];
                            if ( u_var == 2 ) div_u = grad_phi[2][j][2];

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                    * mom_graddiv
                                    * div_u
                                    * div_v;
                        }
                    }
                }
            }
        }

        // ***********************************************************************
        // constraint: theta_inc_c * dT * \int{q div(u)}
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
            {
                DataType div_u = 0.0;

                switch ( u_var )
                {
                    case 0:
                        div_u = grad_phi[0][j][0];
                        break;
                    case 1:
                        div_u = grad_phi[1][j][1];
                        break;
                    case 2:
                        div_u = grad_phi[2][j][2];
                        break;
                }

                int ind_j = this->dof_index ( j, u_var );

                for ( int i = 0; i<this->num_dofs ( p_var ); ++i )
                {
                    lm ( this->dof_index ( i, p_var ), ind_j ) += wq_dJ
                            * inc_c
                            * phi[p_var][i]
                            * div_u;
                }
#    ifdef AUGMENT_PRESS
                for ( int i = 0; i<this->num_dofs ( p0_var ); ++i )
                {
                    lm ( this->dof_index ( i, p0_var ), ind_j ) += wq_dJ
                            * inc_c
                            * phi[p0_var][i]
                            * div_u;
                }
#    endif
            }
        }
    }
}
#endif

/// Residual of primal problem (lhs - rhs) 
#ifdef OPT_ASM0

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );

    const int p_var = DIM;
    const int grad_vars = DIM;

#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#    else
    const int p0_var = -1;
    const int num_vars = DIM + 1;
#    endif        

    DataType sign = 1.0;
    const DataType dt = this->dT_pc_;
    DataType div_c, div_p;

    DataType tau_div = 0.;
    DataType skew_fac = 1.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCartAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::abs ( this->detJ ( q ) );

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
#    ifdef AUGMENT_PRESS
        Vec < DIM + 2, DataType> sol_c;
        Vec < DIM + 2, DataType> sol_p;
        Vec < DIM + 2, DataType> perturb_c;
        Vec < DIM + 2, DataType> perturb_p;
#    else
        Vec < DIM + 1, DataType> sol_c;
        Vec < DIM + 1, DataType> sol_p;
        Vec < DIM + 1, DataType> perturb_c;
        Vec < DIM + 1, DataType> perturb_p;
#    endif
        Vec<DIM, DataType> convection;
        Vec<DIM, DataType> convection_prev;

        // velocity    
        for ( int var = 0; var < DIM; ++var )
        {
            sol_c[var] = this->solP_[var][q];
            sol_p[var] = this->solP_prev_[var][q];
            if ( this->L2_perturb_[var] || this->H1_perturb_[var] )
            {
                perturb_c[var] = this->perturb_vel_[var][q];
                perturb_p[var] = this->perturb_vel_prev_[var][q];
            }
            if ( this->conv_mode_ == OSEEN )
            {
                convection[var] = this->conv_[var][q];
                convection_prev[var] = this->conv_prev_[var][q];
            }
            else if ( this->conv_mode_ == NAVIERSTOKES )
            {
                convection[var] = this->solP_[var][q];
                convection_prev[var] = this->solP_prev_[var][q];
            }
        }

        // Q1 pressure     
        sol_c[p_var] = this->solP_pressure_[q];
        sol_p[p_var] = this->solP_pressure_prev_[q];
        if ( this->L2_perturb_[p_var] )
        {
            perturb_c[p_var] = this->perturb_pressure_[q];
            perturb_p[p_var] = this->perturb_pressure_prev_[q];
        }

        // Q0 pressure     
#    ifdef AUGMENT_PRESS 
        sol_c[p0_var] = this->solP_pressure0_[q];
        sol_p[p0_var] = this->solP_pressure0_prev_[q];
        if ( this->L2_perturb_[p0_var] )
        {
            perturb_c[p0_var] = this->perturb_pressure0_[q];
            perturb_p[p0_var] = this->perturb_pressure0_prev_[q];
        }
#    endif

        std::vector< Vec<DIM, DataType> > grad_sol_c ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_sol_p ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_perturb_c ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_perturb_p ( grad_vars );

        for ( int var = 0; var < grad_vars; var++ )
        {
            for ( int d = 0; d < grad_vars; d++ )
            {
                grad_sol_c[var][d] = this->grad_solP_[var][q][d];
                grad_sol_p[var][d] = this->grad_solP_prev_[var][q][d];

                if ( this->H1_perturb_[var] )
                {
                    grad_perturb_c[var][d] = this->grad_perturb_vel_[var][q][d];
                    grad_perturb_p[var][d] = this->grad_perturb_vel_prev_[var][q][d];
                }
            }
        }

        // ********************************************************************** 
        // time-diff: l0(v) = \int( dot(u_n - u_k, v))
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                lv[this->dof_index ( i, v_var )] += wq
                    * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                    * ( sol_c[v_var] - sol_p[v_var] )
                * this->phi ( i, q, v_var )
                * dJ;
        }

        // ********************************************************************** 
        // LAPLACE  
        // explicit linear: l1n(v) = dT * theta_mom_vis_p * \nu * \int( \grad{u_p} : \grad{v}
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType tmp = 0;
                for ( int s = 0; s < DIM; ++s )
                { // scalar product
                    tmp += grad_sol_p[v_var][s]
                            * this->grad_phi ( i, q, v_var )[s];
                }

                lv[this->dof_index ( i, v_var )] += wq
                        * this->theta_mom_vis_p_
                        * dt
                        * this->nu_
                        * tmp
                        * dJ;
            }
        }

        // l1k(v) = dT * theta_mom_vis_c * \nu * \int( \grad{u_c} : \grad{v}
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType tmp = 0;
                for ( int s = 0; s < DIM; ++s )
                { // scalar product
                    tmp += grad_sol_c[v_var][s]
                            * this->grad_phi ( i, q, v_var )[s];
                }

                lv[this->dof_index ( i, v_var )] += wq
                        * this->theta_mom_vis_c_
                        * dt
                        * this->nu_
                        * tmp
                        * dJ;
            }
        }

        // ********************************************************************** 
        // CONVECTIVE TERM    
        // l2n(v) = (0.5) * theta_mom_adv_cc * dT * \int(u_c*\grad{u_c}*v) 
        //            +    (0.5) * theta_mom_adv_pc * dT * \int(u_p*\grad{u_c}*v) 
        //            +    (0.5) *  theta_mom_adv_cp * dT * \int(u_c*\grad{u_p}*v) 
        //            +   (0.5) *  theta_mom_conv_4 * dT * \int(u_p*\grad{u_p}*v)

        if ( this->conv_mode_ == OSEEN || this->conv_mode_ == NAVIERSTOKES )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                DataType tmp_1 = 0.0;
                DataType tmp_2 = 0.0;
                DataType tmp_3 = 0.0;
                DataType tmp_4 = 0.0;
                for ( int s = 0; s < DIM; s++ )
                {
                    tmp_1 += convection[s] * grad_sol_c[v_var][s];
                    tmp_2 += convection_prev[s] * grad_sol_c[v_var][s];
                    tmp_3 += convection[s] * grad_sol_p[v_var][s];
                    tmp_4 += convection_prev[s] * grad_sol_p[v_var][s];
                }

                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    lv[this->dof_index ( i, v_var )] += wq
                            * skew_fac
                            * dt
                            * ( this->theta_mom_adv_cc_ * tmp_1 +
                            this->theta_mom_adv_pc_ * tmp_2 +
                            this->theta_mom_adv_cp_ * tmp_3 +
                            this->theta_mom_adv_pp_ * tmp_4 )
                            * this->phi ( i, q, v_var )
                            * dJ;
                }
            }
        }

        if ( this->conv_mode_ == NAVIERSTOKES )
        {
            if ( this->skew_mode_ > 0 )
            {
                for ( int v_var = 0; v_var < DIM; ++v_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        DataType tmp_1 = 0.0;
                        DataType tmp_2 = 0.0;
                        DataType tmp_3 = 0.0;
                        DataType tmp_4 = 0.0;

                        for ( int s = 0; s < DIM; s++ )
                        {
                            tmp_1 += this->theta_mom_adv_cc_ * sol_c[s] * this->grad_phi ( i, q, v_var )[s] * sol_c[v_var];
                            tmp_2 += this->theta_mom_adv_pc_ * sol_p[s] * this->grad_phi ( i, q, v_var )[s] * sol_c[v_var];
                            tmp_3 += this->theta_mom_adv_cp_ * sol_c[s] * this->grad_phi ( i, q, v_var )[s] * sol_p[v_var];
                            tmp_4 += this->theta_mom_adv_pp_ * sol_p[s] * this->grad_phi ( i, q, v_var )[s] * sol_p[v_var];
                        }
                        lv[this->dof_index ( i, v_var )] += -wq
                                * skew_fac
                                * dt
                                * ( tmp_1 +
                                tmp_2 +
                                tmp_3 +
                                tmp_4 )
                                * dJ;
                    }
                }
            }
        }

#    ifdef ROTATING_FOR
        // Coriolis force:  -dT * theta_mom_rot_c * 2 * omega (u_c[0] * v[1] - u_c[1] * v[0]) - dT * theta_mom_rot_p * 2 * omega (u_p[0] * v[1] - u_p[1] * v[0])
        // TODO
        /*
            int u_var = 0;
            int v_var = 1;
            for (int i = 0; i < this->num_dofs(v_var); ++i) 
            {
                lv[this->dof_index(i, v_var)] += -wq
         * dt 
         * 2.0
         * this->omega_
         * (this->theta_mom_rot_c_ * sol_c[u_var] + this->theta_mom_rot_p_ * sol_p[u_var])
         * this->phi(i,q,v_var)
         * dJ;
            }

            u_var = 1;
            v_var = 0;
            for (int i = 0; i < this->num_dofs(v_var); ++i)
            {
                lv[this->dof_index(i, v_var)] += wq
         * dt 
         * 2.0
         * this->omega_
         * (this->theta_mom_rot_c_ * sol_c[u_var] + this->theta_mom_rot_p_ * sol_p[u_var])
         * this->phi(i,q,v_var)
         * dJ;
         * */
    }
#    endif 

    // ********************************************************************** 
    // PRESSURE: - theta_mom_pre * 1/rho * dT * \int(p_c*div(v))
    for ( int v_var = 0; v_var < DIM; v_var++ )
    {
        for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
        {
            DataType div_v = 0.0;

            switch ( v_var )
            {
                case 0:
                    div_v = this->grad_phi ( i, q, 0 )[0];
                    break;
                case 1:
                    div_v = this->grad_phi ( i, q, 1 )[1];
                    break;
                case 2:
                    div_v = this->grad_phi ( i, q, 2 )[2];
                    break;
            }
            const int ind_i = this->dof_index ( i, v_var );

            lv[ind_i] += -wq
                    * dt
                    * this->theta_mom_pre_c_
                    * this->inv_rho_
#    ifdef AUGMENT_PRESS
                    * ( sol_c[p_var] + sol_c[p0_var] )
#    else
                    * sol_c[p_var]
#    endif
                    * div_v
                    * dJ;
        }
    }

    div_c = 0.;
    for ( int s = 0; s < DIM; s++ )
        div_c += grad_sol_c[s][s];

    div_p = 0.;
    for ( int s = 0; s < DIM; s++ )
        div_p += grad_sol_p[s][s];

    // ***********************************************************************
    // Grad-div stabilization: theta1 * dT * gamma_div * (div{sol_c}, div{v}) + theta2 * dT * gamma_div * (div{sol_p}, div{v})
    if ( this->graddiv_mode_ > 0 )
    {
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType div_v = 0.0;
                if ( v_var == 0 ) div_v = this->grad_phi ( i, q, 0 )[0];
                if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1];
                if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                lv[this->dof_index ( i, v_var )] += wq
                        * tau_div
                        * dt
                        * ( this->theta_mom_graddiv_c_ * div_c + this->theta_mom_graddiv_p_ * div_p )
                        * div_v
                        * dJ;

            }
        }
    }

    // Perturbation    
    for ( int v_var = 0; v_var < DIM; ++v_var )
    {
        for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
        {
            DataType tmp_H1_p = 0.;
            DataType tmp_L2_p = 0.;
            DataType tmp_H1_c = 0.;
            DataType tmp_L2_c = 0.;

            if ( this->H1_perturb_[v_var] )
            {
                for ( int s = 0; s < DIM; ++s )
                {
                    tmp_H1_p += grad_perturb_p[v_var][s]
                            * this->grad_phi ( i, q, v_var )[s];

                    tmp_H1_c += grad_perturb_c[v_var][s]
                            * this->grad_phi ( i, q, v_var )[s];

                }
            }
            if ( this->L2_perturb_[v_var] )
            {
                tmp_L2_p = perturb_p[v_var] * this->phi ( i, q, v_var );
                tmp_L2_c = perturb_c[v_var] * this->phi ( i, q, v_var );
            }
            lv[this->dof_index ( i, v_var )] += -wq
                    * dt
                    * this->perturb_scale_
                    * 0.5
                    * ( tmp_H1_p + tmp_H1_c + tmp_L2_p + tmp_L2_c )
                    * dJ;
        }
    }
    // ********************************************************************** 
    // CONSTRAINT: theta_inc_c * dt * \int(q * div(u_c)) + theta_inc_p * dt * \int(q * div(u_p))
    for ( int i = 0; i < this->num_dofs ( p_var ); ++i )
    {
        lv[this->dof_index ( i, p_var )] += wq
                * dt
                * ( this->theta_inc_c_ * div_c + this->theta_inc_p_ * div_p )
                * this->phi ( i, q, p_var )
                * dJ;
    }
#    ifdef AUGMENT_PRESS
    for ( int i = 0; i < this->num_dofs ( p0_var ); ++i )
    {
        lv[this->dof_index ( i, p0_var )] += wq
                * dt
                * ( this->theta_inc_c_ * div_c + this->theta_inc_p_ * div_p )
                * this->phi ( i, q, p0_var )
                * dJ;
    }
#    endif

    // Perturbation
    if ( this->L2_perturb_[p_var] )
    {
        for ( int i = 0; i < this->num_dofs ( p_var ); ++i )
        {
            lv[this->dof_index ( i, p_var )] += wq
                    * dt
                    * this->perturb_scale_
                    * 0.5
                    * ( perturb_p[p_var] + perturb_c[p_var] )
                    * this->phi ( i, q, p_var )
                    * dJ;
        }
#    ifdef AUGMENT_PRESS
        for ( int i = 0; i < this->num_dofs ( p0_var ); ++i )
        {
            lv[this->dof_index ( i, p0_var )] += wq
                    * dt
                    * this->perturb_scale_
                    * 0.5
                    * ( perturb_p[p0_var] + perturb_c[p0_var] )
                    * this->phi ( i, q, p0_var )
                    * dJ;
        }
#    endif
    }
}
}
#endif
#ifdef OPT_ASM1

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
    const int grad_vars = DIM;

#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#    else
    const int p0_var = -1;
    const int num_vars = DIM + 1;
#    endif        

    DataType sign = 1.0;
    DataType div_c, div_p;
    const DataType dt = this->dT_pc_;

    DataType tau_div = 0.;
    DataType skew_fac = 1.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCartAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    const DataType mom_vis_c = dt * this->theta_mom_vis_c_ * this->nu_;
    const DataType mom_vis_p = dt * this->theta_mom_vis_p_ * this->nu_;
    const DataType mom_adv_cc = this->theta_mom_adv_cc_ * skew_fac * dt;
    const DataType mom_adv_pc = this->theta_mom_adv_pc_ * skew_fac * dt;
    const DataType mom_adv_cp = this->theta_mom_adv_cp_ * skew_fac * dt;
    const DataType mom_adv_pp = this->theta_mom_adv_pp_ * skew_fac * dt;
    const DataType mom_rot_c = this->theta_mom_rot_c_ * dt * 2.0 * this->omega_;
    const DataType mom_rot_p = this->theta_mom_rot_p_ * dt * 2.0 * this->omega_;
    const DataType mom_pre = dt * this->theta_mom_pre_c_ * this->inv_rho_;
    const DataType mom_graddiv_c = tau_div * this->theta_mom_graddiv_c_ * dt;
    const DataType mom_graddiv_p = tau_div * this->theta_mom_graddiv_p_ * dt;
    const DataType inc_c = dt * this->theta_inc_c_;
    const DataType inc_p = dt * this->theta_inc_p_;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq_dJ = this->w ( q ) * std::abs ( this->detJ ( q ) );

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
#    ifdef AUGMENT_PRESS
        Vec < DIM + 2, DataType> sol_c;
        Vec < DIM + 2, DataType> sol_p;
#    else
        Vec < DIM + 1, DataType> sol_c;
        Vec < DIM + 1, DataType> sol_p;
#    endif

        Vec<DIM, DataType> convection;
        Vec<DIM, DataType> convection_prev;

        for ( int var = 0; var < DIM; ++var )
        {
            sol_c[var] = this->solP_[var][q];
        }
        sol_c[DIM] = this->solP_pressure_[q];

        for ( int var = 0; var < DIM; ++var )
        {
            sol_p[var] = this->solP_prev_[var][q];
        }
        sol_p[DIM] = this->solP_pressure_prev_[q];

        for ( int var = 0; var < DIM; ++var )
        {
            if ( this->conv_mode_ == OSEEN )
            {
                convection[var] = this->conv_[var][q];
                convection_prev[var] = this->conv_prev_[var][q];
            }
            else if ( this->conv_mode_ == NAVIERSTOKES )
            {
                convection[var] = this->solP_[var][q];
                convection_prev[var] = this->solP_prev_[var][q];
            }
        }

        // Q0 pressure     
#    ifdef AUGMENT_PRESS 
        sol_c[p0_var] = this->solP_pressure0_[q];
        sol_p[p0_var] = this->solP_pressure0_prev_[q];
#    endif

        std::vector< Vec<DIM, DataType> > grad_sol_c ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_sol_p ( grad_vars );

        for ( int var = 0; var < grad_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                if ( var < p_var )
                {
                    grad_sol_c[var][d] = this->grad_solP_[var][q][d];
                    grad_sol_p[var][d] = this->grad_solP_prev_[var][q][d];
                }
            }
        }

        // ********************************************************************** 
        // time-diff: l0(v) = \int( dot(u_n - u_k, v))
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                lv[this->dof_index ( i, v_var )] += wq_dJ
                    * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                    * ( sol_c[v_var] - sol_p[v_var] )
                * this->phi ( i, q, v_var );
        }

        // ********************************************************************** 
        // LAPLACE  
        // explicit linear: l1n(v) = dT * theta_mom_vis_p * \nu * \int( \grad{u_p} : \grad{v}
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType tmp = 0;
                for ( int s = 0; s < DIM; ++s )
                { // scalar product
                    tmp += grad_sol_p[v_var][s]
                            * this->grad_phi ( i, q, v_var )[s];
                }

                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * mom_vis_p
                        * tmp;
            }
        }

        // l1k(v) = dT * theta_mom_vis_c * \nu * \int( \grad{u_c} : \grad{v}
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType tmp = 0;
                for ( int s = 0; s < DIM; ++s )
                { // scalar product
                    tmp += grad_sol_c[v_var][s]
                            * this->grad_phi ( i, q, v_var )[s];
                }
                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * mom_vis_c
                        * tmp;
            }
        }

        // ********************************************************************** 
        // CONVECTIVE TERM    
        // l2n(v) = (0.5) * theta_mom_adv_cc * dT * \int(u_c*\grad{u_c}*v) 
        //            +    (0.5) * theta_mom_adv_pc * dT * \int(u_p*\grad{u_c}*v) 
        //            +    (0.5) *  theta_mom_adv_cp * dT * \int(u_c*\grad{u_p}*v) 
        //            +   (0.5) *  theta_mom_conv_4 * dT * \int(u_p*\grad{u_p}*v)

        if ( this->conv_mode_ == OSEEN || this->conv_mode_ == NAVIERSTOKES )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                DataType tmp_1 = 0.0;
                DataType tmp_2 = 0.0;
                DataType tmp_3 = 0.0;
                DataType tmp_4 = 0.0;
                for ( int s = 0; s < DIM; s++ )
                {
                    tmp_1 += sol_c[s] * grad_sol_c[v_var][s];
                    tmp_2 += sol_p[s] * grad_sol_c[v_var][s];
                    tmp_3 += sol_c[s] * grad_sol_p[v_var][s];
                    tmp_4 += sol_p[s] * grad_sol_p[v_var][s];
                }

                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    lv[this->dof_index ( i, v_var )] += wq_dJ
                            * ( mom_adv_cc * tmp_1 +
                            mom_adv_pc * tmp_2 +
                            mom_adv_cp * tmp_3 +
                            mom_adv_pp * tmp_4 )
                            * this->phi ( i, q, v_var );
                }
            }
        }

        if ( this->conv_mode_ == NAVIERSTOKES )
        {
            if ( this->skew_mode_ > 0 )
            {
                for ( int v_var = 0; v_var < DIM; ++v_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        DataType tmp_1 = 0.0;
                        DataType tmp_2 = 0.0;
                        DataType tmp_3 = 0.0;
                        DataType tmp_4 = 0.0;

                        for ( int s = 0; s < DIM; s++ )
                        {
                            tmp_1 += sol_c[s] * this->grad_phi ( i, q, v_var )[s] * sol_c[v_var];
                            tmp_2 += sol_p[s] * this->grad_phi ( i, q, v_var )[s] * sol_c[v_var];
                            tmp_3 += sol_c[s] * this->grad_phi ( i, q, v_var )[s] * sol_p[v_var];
                            tmp_4 += sol_p[s] * this->grad_phi ( i, q, v_var )[s] * sol_p[v_var];
                        }

                        lv[this->dof_index ( i, v_var )] += -wq_dJ
                                * ( mom_adv_cc * tmp_1 +
                                mom_adv_pc * tmp_2 +
                                mom_adv_cp * tmp_3 +
                                mom_adv_pp * tmp_4 );
                    }
                }
            }
        }

#    ifdef ROTATING_FOR
        // Coriolis force:  -dT * theta_mom_rot_c * 2 * omega (u_c[0] * v[1] - u_c[1] * v[0]) - dT * theta_mom_rot_p * 2 * omega (u_p[0] * v[1] - u_p[1] * v[0])
        // TODO
        /*
        int u_var = 0;
        int v_var = 1;
        for (int i = 0; i < this->num_dofs(v_var); ++i) 
        {
            lv[this->dof_index(i, v_var)] += -wq_dJ
         * (mom_rot_c * sol_c[u_var] + mom_rot_p * sol_p[u_var])
         * this->phi(i,q,v_var);
        }

        u_var = 1;
        v_var = 0;
        for (int i = 0; i < this->num_dofs(v_var); ++i) 
        {
            lv[this->dof_index(i, v_var)] += wq_dJ
         * (mom_rot_c * sol_c[u_var] + mom_rot_p * sol_p[u_var])
         * this->phi(i,q,v_var);
        }
         * */
#    endif 

        // ********************************************************************** 
        // PRESSURE: - theta_mom_pre * 1/rho * dT * \int(p_c*div(v))
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType div_v = 0.0;

                switch ( v_var )
                {
                    case 0:
                        div_v = this->grad_phi ( i, q, 0 )[0];
                        break;
                    case 1:
                        div_v = this->grad_phi ( i, q, 1 )[1];
                        break;
                    case 2:
                        div_v = this->grad_phi ( i, q, 2 )[2];
                        break;
                }
                const int ind_i = this->dof_index ( i, v_var );

                lv[ind_i] += -wq_dJ
                        * mom_pre
#    ifdef AUGMENT_PRESS
                        * ( sol_c[p_var] + sol_c[p0_var] )
#    else
                        * sol_c[p_var]
#    endif
                        * div_v;
            }
        }

        div_c = 0.;
        for ( int s = 0; s < DIM; s++ )
            div_c += grad_sol_c[s][s];

        div_p = 0.;
        for ( int s = 0; s < DIM; s++ )
            div_p += grad_sol_p[s][s];

        // ***********************************************************************
        // Grad-div stabilization: theta1 * dT * gamma_div * (div{sol_c}, div{v}) + theta2 * dT * gamma_div * (div{sol_p}, div{v})
        if ( this->graddiv_mode_ > 0 )
        {

            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    DataType div_v = 0.0;
                    if ( v_var == 0 ) div_v = this->grad_phi ( i, q, 0 )[0];
                    if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1];
                    if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                    lv[this->dof_index ( i, v_var )] += wq_dJ
                            * ( mom_graddiv_c * div_c + mom_graddiv_p * div_p )
                            * div_v;

                }
            }
        }

        // ********************************************************************** 
        // CONSTRAINT: theta_inc_c * dt * \int(q * div(u_c)) + theta_inc_p * dt * \int(q * div(u_p))

        for ( int i = 0; i < this->num_dofs ( p_var ); ++i )
        {
            lv[this->dof_index ( i, p_var )] += wq_dJ
                    * ( inc_c * div_c + inc_p * div_p )
                    * this->phi ( i, q, p_var );
        }
#    ifdef AUGMENT_PRESS
        for ( int i = 0; i < this->num_dofs ( p0_var ); ++i )
        {
            lv[this->dof_index ( i, p0_var )] += wq_dJ
                    * ( inc_c * div_c + inc_p * div_p )
                    * this->phi ( i, q, p0_var );
        }
#    endif
    }
}
#endif
#ifdef OPT_ASM2

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
    const int grad_vars = DIM;

#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#    else
    const int p0_var = -1;
    const int num_vars = DIM + 1;
#    endif

    DataType sign = 1.0;

    DataType div_c, div_p;
    const DataType dt = this->dT_pc_;

    DataType tau_div = 0.;
    DataType skew_fac = 1.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCartAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    const DataType mom_vis_c = dt * this->theta_mom_vis_c_ * this->nu_;
    const DataType mom_vis_p = dt * this->theta_mom_vis_p_ * this->nu_;
    const DataType mom_adv_cc = this->theta_mom_adv_cc_ * skew_fac * dt;
    const DataType mom_adv_pc = this->theta_mom_adv_pc_ * skew_fac * dt;
    const DataType mom_adv_cp = this->theta_mom_adv_cp_ * skew_fac * dt;
    const DataType mom_adv_pp = this->theta_mom_adv_pp_ * skew_fac * dt;
    const DataType mom_rot_c = this->theta_mom_rot_c_ * dt * 2.0 * this->omega_;
    const DataType mom_rot_p = this->theta_mom_rot_p_ * dt * 2.0 * this->omega_;
    const DataType mom_pre = dt * this->theta_mom_pre_c_ * this->inv_rho_;
    const DataType mom_graddiv_c = tau_div * this->theta_mom_graddiv_c_ * dt;
    const DataType mom_graddiv_p = tau_div * this->theta_mom_graddiv_p_ * dt;
    const DataType inc_c = dt * this->theta_inc_c_;
    const DataType inc_p = dt * this->theta_inc_p_;

    std::vector< std::vector<DataType> > phi;
    std::vector< std::vector< std::vector<DataType> > > grad_phi;
    std::vector< std::vector< Mat<DIM, DIM, DataType> > > H_phi;
#    ifdef AUGMENT_PRESS    
    phi.resize ( DIM + 2 );
    phi[p0_var].resize ( this->num_dofs ( p0_var ), 0. );
#    else 
    phi.resize ( DIM + 1 );
#    endif

    grad_phi.resize ( grad_vars );
    for ( int k = 0; k < num_vars; ++k )
    {
        phi[k].resize ( this->num_dofs ( k ), 0. );
        if ( k < grad_vars )
        {
            grad_phi[k].resize ( this->num_dofs ( k ) );
            for ( int j = 0; j<this->num_dofs ( k ); ++j )
            {
                grad_phi[k][j].resize ( DIM, 0. );
            }
        }
    }

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq_dJ = this->w ( q ) * std::abs ( this->detJ ( q ) );

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
#    ifdef AUGMENT_PRESS    
        Vec < DIM + 2, DataType> sol_c;
        Vec < DIM + 2, DataType> sol_p;
#    else
        Vec < DIM + 1, DataType> sol_c;
        Vec < DIM + 1, DataType> sol_p;
#    endif    

        Vec<DIM, DataType> convection;
        Vec<DIM, DataType> convection_prev;

        for ( int var = 0; var < DIM; ++var )
            sol_c[var] = this->solP_[var][q];

        sol_c[p_var] = this->solP_pressure_[q];

        for ( int var = 0; var < DIM; ++var )
            sol_p[var] = this->solP_prev_[var][q];

        sol_p[p_var] = this->solP_pressure_prev_[q];

        for ( int var = 0; var < DIM; ++var )
        {
            if ( this->conv_mode_ == OSEEN )
            {
                convection[var] = this->conv_[var][q];
                convection_prev[var] = this->conv_prev_[var][q];
            }
            else if ( this->conv_mode_ == NAVIERSTOKES )
            {
                convection[var] = this->solP_[var][q];
                convection_prev[var] = this->solP_prev_[var][q];
            }
        }

        // Q0 pressure     
#    ifdef AUGMENT_PRESS 
        sol_c[p0_var] = this->solP_pressure0_[q];
        sol_p[p0_var] = this->solP_pressure0_prev_[q];
#    endif

        std::vector< Vec<DIM, DataType> > grad_sol_c ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_sol_p ( grad_vars );

        for ( int var = 0; var < grad_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_sol_c[var][d] = this->grad_solP_[var][q][d];
                grad_sol_p[var][d] = this->grad_solP_prev_[var][q][d];
            }
        }

        for ( int k = 0; k < num_vars; ++k )
        {
            for ( int j = 0; j<this->num_dofs ( k ); ++j )
            {
                phi[k][j] = this->phi ( j, q, k );
                if ( k < grad_vars )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_phi[k][j][d] = this->grad_phi ( j, q, k )[d];
                    }
                }
            }
        }

        // ********************************************************************** 
        // time-diff: l0(v) = \int( dot(u_n - u_k, v))
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                lv[this->dof_index ( i, v_var )] += wq_dJ
                    * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                    * ( sol_c[v_var] - sol_p[v_var] )
                * phi[v_var][i];
        }

        // ********************************************************************** 
        // LAPLACE  
        // explicit linear: l1n(v) = dT * theta_mom_vis_p * \nu * \int( \grad{u_p} : \grad{v}
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType tmp = 0;
                for ( int s = 0; s < DIM; ++s )
                { // scalar product
                    tmp += grad_sol_p[v_var][s]
                            * grad_phi[v_var][i][s];
                }

                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * mom_vis_p
                        * tmp;
            }
        }

        // l1k(v) = dT * theta_mom_vis_c * \nu * \int( \grad{u_c} : \grad{v}
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType tmp = 0;
                for ( int s = 0; s < DIM; ++s )
                { // scalar product
                    tmp += grad_sol_c[v_var][s]
                            * grad_phi[v_var][i][s];
                }

                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * mom_vis_c
                        * tmp;
            }
        }

        // ********************************************************************** 
        // CONVECTIVE TERM    
        // l2n(v) = (0.5) * theta_mom_adv_cc * dT * \int(u_c*\grad{u_c}*v) 
        //            +    (0.5) * theta_mom_adv_pc * dT * \int(u_p*\grad{u_c}*v) 
        //            +    (0.5) *  theta_mom_adv_cp * dT * \int(u_c*\grad{u_p}*v) 
        //            +   (0.5) *  theta_mom_conv_4 * dT * \int(u_p*\grad{u_p}*v)

        if ( this->conv_mode_ == OSEEN || this->conv_mode_ == NAVIERSTOKES )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                DataType tmp_1 = 0.0;
                DataType tmp_2 = 0.0;
                DataType tmp_3 = 0.0;
                DataType tmp_4 = 0.0;
                for ( int s = 0; s < DIM; s++ )
                {
                    tmp_1 += sol_c[s] * grad_sol_c[v_var][s];
                    tmp_2 += sol_p[s] * grad_sol_c[v_var][s];
                    tmp_3 += sol_c[s] * grad_sol_p[v_var][s];
                    tmp_4 += sol_p[s] * grad_sol_p[v_var][s];
                }

                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    lv[this->dof_index ( i, v_var )] += wq_dJ
                            * ( mom_adv_cc * tmp_1 +
                            mom_adv_pc * tmp_2 +
                            mom_adv_cp * tmp_3 +
                            mom_adv_pp * tmp_4 )
                            * phi[v_var][i];
                }
            }
        }

        if ( this->conv_mode_ == NAVIERSTOKES )
        {
            if ( this->skew_mode_ > 0 )
            {
                for ( int v_var = 0; v_var < DIM; ++v_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        DataType tmp_1 = 0.0;
                        DataType tmp_2 = 0.0;
                        DataType tmp_3 = 0.0;
                        DataType tmp_4 = 0.0;

                        for ( int s = 0; s < DIM; s++ )
                        {
                            tmp_1 += sol_c[s] * grad_phi[v_var][i][s] * sol_c[v_var];
                            tmp_2 += sol_p[s] * grad_phi[v_var][i][s] * sol_c[v_var];
                            tmp_3 += sol_c[s] * grad_phi[v_var][i][s] * sol_p[v_var];
                            tmp_4 += sol_p[s] * grad_phi[v_var][i][s] * sol_p[v_var];
                        }
                        lv[this->dof_index ( i, v_var )] += -wq_dJ
                                * ( mom_adv_cc * tmp_1 +
                                mom_adv_pc * tmp_2 +
                                mom_adv_cp * tmp_3 +
                                mom_adv_pp * tmp_4 );
                    }
                }
            }
        }

#    ifdef ROTATING_FOR
        // Coriolis force:  -dT * theta_mom_rot_c * 2 * omega (u_c[0] * v[1] - u_c[1] * v[0]) - dT * theta_mom_rot_p * 2 * omega (u_p[0] * v[1] - u_p[1] * v[0])
        // TODO
        /*
        int u_var = 0;
        int v_var = 1;
        for (int i = 0; i < this->num_dofs(v_var); ++i) 
        {
            lv[this->dof_index(i, v_var)] += -wq_dJ
         * (mom_rot_c * sol_c[u_var] + mom_rot_p * sol_p[u_var])
         * phi[v_var][i];
        }

        u_var = 1;
        v_var = 0;
        for (int i = 0; i < this->num_dofs(v_var); ++i) 
        {
            lv[this->dof_index(i, v_var)] += wq_dJ
         * (mom_rot_c * sol_c[u_var] + mom_rot_p * sol_p[u_var])
         * phi[v_var][i];
        }
         * */
#    endif 

        // ********************************************************************** 
        // PRESSURE: - theta_mom_pre * 1/rho * dT * \int(p_c*div(v))
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                DataType div_v = 0.0;

                if ( v_var == 0 ) div_v = this->grad_phi ( i, q, 0 )[0];
                if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1];
                if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                int ind_i = this->dof_index ( i, v_var );

                lv[ind_i] += -wq_dJ
                        * mom_pre
#    ifdef AUGMENT_PRESS
                        * ( sol_c[p_var] + sol_c[p0_var] )
#    else
                        * sol_c[p_var]
#    endif
                        * div_v;
            }
        }

        div_c = 0.;
        for ( int s = 0; s < DIM; s++ )
            div_c += grad_sol_c[s][s];

        div_p = 0.;
        for ( int s = 0; s < DIM; s++ )
            div_p += grad_sol_p[s][s];

        // ***********************************************************************
        // Grad-div stabilization: theta1 * dT * gamma_div * (div{sol_c}, div{v}) + theta2 * dT * gamma_div * (div{sol_p}, div{v})
        if ( this->graddiv_mode_ > 0 )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    DataType div_v = 0.0;
                    if ( v_var == 0 ) div_v = this->grad_phi ( i, q, 0 )[0];
                    if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1];
                    if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                    lv[this->dof_index ( i, v_var )] += wq_dJ
                            * ( mom_graddiv_c * div_c + mom_graddiv_p * div_p )
                            * div_v;

                }
            }
        }

        // ********************************************************************** 
        // CONSTRAINT: theta_inc_c * dt * \int(q * div(u_c)) + theta_inc_p * dt * \int(q * div(u_p))
        for ( int i = 0; i < this->num_dofs ( p_var ); ++i )
        {
            lv[this->dof_index ( i, p_var )] += wq_dJ
                    * ( inc_c * div_c + inc_p * div_p )
                    * phi[p_var][i];
        }
#    ifdef AUGMENT_PRESS
        for ( int i = 0; i < this->num_dofs ( p0_var ); ++i )
        {
            lv[this->dof_index ( i, p0_var )] += wq_dJ
                    * ( inc_c * div_c + inc_p * div_p )
                    * phi[p0_var][i];
        }
#    endif    

    }
}
#endif

/// vorticity

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_vort ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#else
    const int p0_var = -1;
    const int num_vars = DIM + 1;
#endif    

    const int num_q = this->num_quadrature_points ( );

    // add electro-dynamic parts    
    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::fabs ( this->detJ ( q ) );

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index

        Vec<DIM, DataType> sol_c;
        for ( int var = 0; var < DIM; ++var ) sol_c[var] = this->solP_[var][q];

        std::vector< Vec<DIM, DataType> > grad_sol_c ( DIM + 3 );
        for ( int var = 0; var < DIM; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_sol_c[var][d] = this->grad_solP_[var][q][d];
            }
        }

        // ***************************************************************
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            DataType vort = 0.;
            switch ( v_var )
            {
                case 0:
                    vort = grad_sol_c[1][2] - grad_sol_c[2][1];
                    break;
                case 1:
                    vort = grad_sol_c[2][0] - grad_sol_c[0][2];
                    break;
                case 2:
                    vort = grad_sol_c[0][1] - grad_sol_c[1][0];
                    break;
            }
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += wq
                        * vort
                        * this->phi ( i, q, v_var )
                        * dJ;
            }
        }
    }
}

/// Initial condition for dual problem

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_goal ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#else
    const int p0_var = -1;
    const int num_vars = DIM + 1;
#endif    

    const int num_q = this->num_quadrature_points ( );

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::fabs ( this->detJ ( q ) );

        ParametersFinalType<DIM, DataType> p;
        p.solP.resize ( num_vars );
        p.grad_solP.resize ( num_vars );

        for ( int var = 0; var < DIM; ++var )
            p.solP[var] = this->solP_[var][q];

        p.solP[p_var] = this->solP_pressure_[q];
#ifdef AUGMENT_PRESS
        p.solP[p0_var] = this->solP_pressure0_[q];
#endif    
        for ( int var = 0; var < grad_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                p.grad_solP[var][d] = this->grad_solP_[var][q][d];
            }
        }
        p.x = this->x ( q );

        for ( int var = 0; var < num_vars; ++var )
        {
            p.var = var;
            for ( int i = 0; i<this->num_dofs ( var ); ++i )
            {
                p.phi = this->phi ( i, q, var );
                p.grad_phi = this->grad_phi ( i, q, var );
                lv[this->dof_index ( i, var )] += -wq * this->goal_functional_->j_final_type ( p ) * dJ;
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_vector_quant ( const Element<DataType>& element, LocalVector& lv ) const
{

    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#else
    const int p0_var = -1;
    const int num_vars = DIM + 1;
#endif
    DataType sign = 1.0;

    DataType div_c, div_p;

    DataType tau = 0.;
    if ( this->graddiv_mode_ > 0 )
        this->compute_tau_graddiv ( element, tau );

    const int num_q = this->num_quadrature_points ( );
    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::abs ( this->detJ ( q ) );

        for ( int l = 0; l<this->incomp_scalars_; ++l )
        {
            DataType kin = 0.;
            DataType pot = 0.;
            DataType grad_2 = 0.0;
            DataType div = 0.;
            DataType vort = 0.;
            DataType diff = 0.;
            switch ( l )
            {
                case 0:
                    // kinetic energy                    
                    for ( int var = 0; var < DIM; ++var )
                        kin += this->solP_[var][q] * this->solP_[var][q];

                    lv[l] += wq * 0.5 * kin * dJ;
                    break;
                case 1:
                    // azimuthal kinetic energy
                    lv[l] += wq * 0.5 * this->solP_[0][q] * this->solP_[0][q] * dJ;
                    break;
                case 2:
                    // radial kinetic energy 
                    lv[l] += wq * 0.5 * this->solP_[1][q] * this->solP_[1][q] * dJ;
                    break;
                case 3:
                    // axial kinetic energy 
                    lv[l] += wq * 0.5 * this->solP_[DIM - 1][q] * this->solP_[DIM - 1][q] * dJ;
                    break;
                case 4:
                    // energy dissipation                     
                    for ( int var = 0; var < DIM; var++ )
                    {
                        for ( int d = 0; d < DIM; d++ )
                        {
                            pot += this->grad_solP_[var][q][d] * this->grad_solP_[var][q][d];
                        }
                    }

                    lv[l] += wq * pot * dJ;

                    break;
                case 5:
                    // azimuthal energy dissipation                        
                    for ( int s = 0; s < DIM; s++ )
                        grad_2 += this->grad_solP_[0][q][s] * this->grad_solP_[0][q][s];

                    lv[l] += wq * this->nu_ * ( -1. * grad_2 ) * dJ;
                    break;
                case 6:
                    // radial energy dissipation                     
                    for ( int s = 0; s < DIM; s++ )
                        grad_2 += this->grad_solP_[1][q][s] * this->grad_solP_[1][q][s];

                    lv[l] += wq * this->nu_ * ( -1. * grad_2 ) * dJ;
                    break;
                case 7:
                    // axial energy dissipation                    
                    for ( int s = 0; s < DIM; s++ )
                        grad_2 += this->grad_solP_[2][q][s] * this->grad_solP_[2][q][s];

                    lv[l] += wq * this->nu_ * ( -1. * grad_2 ) * dJ;
                    break;
                case 8:
                    // azimuthal pressure 
                    lv[l] += wq * this->inv_rho_ * this->grad_solP_[0][q][0] * this->solP_pressure_[q] * dJ;
#ifdef AUGMENT_PRESS
                    lv[l] += wq * this->inv_rho_ * this->grad_solP_[0][q][0] * this->solP_pressure0_[q] * dJ;
#endif                    
                    break;
                case 9:
                    // radial pressure 
                    lv[l] += wq * this->inv_rho_ * ( this->grad_solP_[1][q][1] ) * this->solP_pressure_[q] * dJ;
#ifdef AUGMENT_PRESS
                    lv[l] += wq * this->inv_rho_ * ( this->grad_solP_[1][q][1] ) * this->solP_pressure0_[q] * dJ;
#endif                    
                    break;
                case 10:
                    // axial pressure 
                    lv[l] += wq * this->inv_rho_ * this->grad_solP_[2][q][2] * this->solP_pressure_[q] * dJ;
#ifdef AUGMENT_PRESS
                    lv[l] += wq * this->inv_rho_ * this->grad_solP_[2][q][2] * this->solP_pressure0_[q] * dJ;
#endif                    
                    break;
                case 11:
                    // azimuthal grad div contribution                    
                    if ( this->graddiv_mode_ > 0 )
                    {
                        div = this->grad_solP_[0][q][0] + this->grad_solP_[1][q][1];
                        if ( DIM == 3 ) div += this->grad_solP_[2][q][2];
                        lv[l] -= wq * tau * div * this->grad_solP_[0][q][0] * dJ;
                    }
                    else
                        lv[l] = 0.;
                    break;
                case 12:
                    // radial grad div contribution                    
                    if ( this->graddiv_mode_ > 0 )
                    {
                        div = this->grad_solP_[0][q][0] + this->grad_solP_[1][q][1];
                        if ( DIM == 3 ) div += this->grad_solP_[2][q][2];
                        lv[l] -= wq * tau * div * this->grad_solP_[0][q][0] * dJ;
                    }
                    else
                        lv[l] = 0.;
                    break;
                case 13:
                    // axial grad div contribution            
                    if ( this->graddiv_mode_ > 0 )
                    {
                        div = this->grad_solP_[0][q][0] + this->grad_solP_[1][q][1];
                        if ( DIM == 3 ) div += this->grad_solP_[2][q][2];

                        lv[l] -= wq * div * this->grad_solP_[2][q][2] * dJ;
                    }
                    else
                        lv[l] = 0.;
                    break;
                case 14:
                    // azimuthal to radial convection 
                    lv[l] -= wq * this->solP_[1][q] * this->solP_[0][q] * this->solP_[0][q] * dJ;
                    break;
                case 15:
                    // radial to azimuthal convection 
                    lv[l] += wq * this->solP_[1][q] * this->solP_[0][q] * this->solP_[0][q] * dJ;
                    break;
                case 16:
                    // divergence 
                    div = this->grad_solP_[0][q][0] + this->grad_solP_[1][q][1];
                    if ( DIM == 3 ) div += this->grad_solP_[2][q][2];

                    lv[l] += wq * div * div * dJ;
                    break;
                case 17:
                    // azimuthal vorticity
                    vort = this->grad_solP_[1][q][2] - this->grad_solP_[2][q][1];
                    lv[l] += wq * vort * vort * dJ;
                    break;
                case 18:
                    // radial vorticity 
                    vort = this->grad_solP_[2][q][0] - this->grad_solP_[0][q][2];
                    lv[l] += wq * vort * vort * dJ;
                    break;
                case 19:
                    // axial vorticity
                    vort = this->grad_solP_[0][q][1] - this->grad_solP_[1][q][0];
                    lv[l] += wq * vort * vort * dJ;
                    break;
                case 20:
                    // azimuthal velocity difference from base state
                    if ( this->vector_base_ != NULL )
                    {
                        diff = ( this->base_vel_[0][q] - this->solP_[0][q] );
                        lv[l] += wq * diff * diff * dJ;
                    }
                    break;
                case 21:
                    // radial velocity difference from base state
                    if ( this->vector_base_ != NULL )
                    {
                        diff = ( this->base_vel_[1][q] - this->solP_[1][q] );
                        lv[l] += wq * diff * diff * dJ;
                    }
                    break;
                case 22:
                    // axial velocity difference from base state
                    if ( this->vector_base_ != NULL )
                    {
                        diff = ( this->base_vel_[2][q] - this->solP_[2][q] );
                        lv[l] += wq * diff * diff * dJ;
                    }
                    break;
                case 23:
                    // pressure difference from base state
                    if ( this->vector_base_ != NULL )
                    {
                        diff = ( this->base_pressure_[q] - this->solP_pressure_[q] );
#ifdef AUGMENT_PRESS
                        diff += ( this->base_pressure0_[q] - this->solP_pressure0_[q] );
#endif                        
                        lv[l] += wq * diff * diff * dJ;
                    }
                    break;
            }
        }
    }
}

/// squared L2- W1,2 norm of complete solution vector 

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_goal_int ( const Element<DataType>& element, DataType& ls ) const
{

    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#else
    const int p0_var = -1;
    const int num_vars = DIM + 1;
#endif

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::abs ( this->detJ ( q ) );

        std::vector<DataType> solP_c ( num_vars );

        for ( int var = 0; var < DIM; ++var )
        {
            solP_c[var] = this->solP_[var][q];
        }

        solP_c[p_var] = this->solP_pressure_[q];
#ifdef AUGMENT_PRESS
        solP_c[p0_var] = this->solP_pressure0_[q];
#endif

        std::vector< Vec<DIM, DataType> > grad_solP_c ( num_vars );

        for ( int var = 0; var < p_var; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_solP_c[var][d] = this->grad_solP_[var][q][d];
            }
        }

        ParametersEvalType<DIM, DataType> gp;
        gp.x = this->x ( q );
        gp.absolute_time = this->t_;
        gp.solP = solP_c;
        gp.grad_solP = grad_solP_c;

        ls += wq * this->goal_functional_->j_force_eval ( gp ) * dJ;
    }
}

/// squared L2- W1,2 norm of complete solution vector 

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_goal_fin ( const Element<DataType>& element, DataType& ls ) const
{

    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#else
    const int p0_var = -1;
    const int num_vars = DIM + 1;
#endif

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::abs ( this->detJ ( q ) );

        std::vector<DataType> solP_c ( num_vars );

        for ( int var = 0; var < DIM; ++var )
        {
            solP_c[var] = this->solP_[var][q];
        }

        solP_c[p_var] = this->solP_pressure_[q];
#ifdef AUGMENT_PRESS
        solP_c[p0_var] = this->solP_pressure0_[q];
#endif

        std::vector< Vec<DIM, DataType> > grad_solP_c ( num_vars );

        for ( int var = 0; var < p_var; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_solP_c[var][d] = this->grad_solP_[var][q][d];
            }
        }

        ParametersEvalType<DIM, DataType> gp;
        gp.x = this->x ( q );
        gp.solP = solP_c;
        gp.grad_solP = grad_solP_c;

        ls += wq * this->goal_functional_->j_final_eval ( gp ) * dJ;
    }
}

/// mean value of divergence of velocity

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_div_mean ( const Element<DataType>& element, DataType& ls ) const
{

    const int num_q = this->num_quadrature_points ( );
    const int total_dofs = this->num_dofs_total ( );
    DataType sign = 1.0;

    DataType int_div = 0.0;
    DataType int_1 = 0.0;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::abs ( this->detJ ( q ) );

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
        Vec<DIM, DataType> sol_p;

        for ( int var = 0; var < DIM; ++var ) sol_p[var] = this->solP_[var][q];

        std::vector< Vec<DIM, DataType> > grad_sol_p ( DIM );

        for ( int var = 0; var < DIM; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_sol_p[var][d] = this->grad_solP_[var][q][d];
            }
        }

        DataType div = grad_sol_p[0][0] + grad_sol_p[1][1];
        if ( DIM == 3 ) div += grad_sol_p[2][2];

        int_div += wq * div * dJ;
        int_1 += wq * dJ;
    }
    ls = int_div / int_1;
}

/// mean value of divergence of velocity

template<int DIM, class DataType>
void MetFlowIncompCartAssembler<DIM, DataType>::assemble_local_scalar_div_max ( const Element<DataType>& element, DataType& ls ) const
{

    const int num_q = this->num_quadrature_points ( );
    const int total_dofs = this->num_dofs_total ( );
    DataType sign = 1.0;

    DataType int_div = 0.0;
    DataType int_1 = 0.0;
    ls = 0.;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::abs ( this->detJ ( q ) );

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
        Vec<DIM, DataType> sol_p;

        for ( int var = 0; var < DIM; ++var ) sol_p[var] = this->solP_[var][q];

        std::vector< Vec<DIM, DataType> > grad_sol_p ( DIM );

        for ( int var = 0; var < DIM; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_sol_p[var][d] = this->grad_solP_[var][q][d];
            }
        }

        DataType div = grad_sol_p[0][0] + grad_sol_p[1][1];
        if ( DIM == 3 ) div += grad_sol_p[2][2];
        div = std::abs ( div );

        if ( div > ls )
            ls = div;
    }
}

template class MetFlowIncompCartAssembler<2, double>;
template class MetFlowIncompCartAssembler<3, double>;

