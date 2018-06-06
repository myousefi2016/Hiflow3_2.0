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

#include "met_flow_bous_cyl_assembler.h"

template<int DIM, class DataType>
MetFlowBousCylAssembler<DIM, DataType>::MetFlowBousCylAssembler ( )
: MetFlowIncompCylAssembler<DIM, DataType>( ),
MetFlowBousAssembler<DIM, DataType>( )
{
}

/// ********************************************************
/// General Assembly routines 
/// ********************************************************

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    if ( this->mode_ == PRIMAL )
    {
        MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_matrix_primal ( element, lm );
        MetFlowBousCylAssembler <DIM, DataType>::assemble_local_matrix_primal ( element, lm );
    }
    else
    {
        interminable_assert ( 0 );
    }
}

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    if ( this->vector_asm_mode_ == VECTOR_STD )
    {
        if ( this->mode_ == PRIMAL )
        {
            MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_vector_primal ( element, lv );
            MetFlowBousCylAssembler <DIM, DataType>::assemble_local_vector_primal ( element, lv );
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
        MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_vector_quant ( element, lv );
        MetFlowBousCylAssembler <DIM, DataType>::assemble_local_vector_quant ( element, lv );
    }
    else if ( this->vector_asm_mode_ == VECTOR_VORT )
    {
        MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_vector_vort ( element, lv );
    }
    else if ( this->vector_asm_mode_ == VECTOR_GOAL )
    {
        MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_vector_goal ( element, lv );
        MetFlowBousCylAssembler <DIM, DataType>::assemble_local_vector_goal ( element, lv );
    }
    else if ( this->vector_asm_mode_ == VECTOR_GRADTEMP )
    {
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_grad_temp ( element, lv );
    }
    else if ( this->vector_asm_mode_ == VECTOR_FORCE )
    {
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_buoyancy ( element, lv );
    }
}

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const
{
    ls = 0.;
    MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_scalar ( element, ls );

    if ( this->sca_mode_ == TEMP_MEAN )
    {
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_temp_mean ( element, ls );
    }
    if ( this->sca_mode_ == RAD_HEAT_VOLUME )
    {
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_rad_heat_volume ( element, ls );
    }
    if ( this->sca_mode_ == GOAL_INT )
    {
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_goal_int ( element, ls );
    }
    if ( this->sca_mode_ == GOAL_FIN )
    {
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_goal_fin ( element, ls );
    }
}

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const
{
    MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_scalar_boundary ( element, facet_number, lv );
    if ( this->sca_mode_ == RAD_HEAT_SURFACE )
    {
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_rad_heat_surface ( element, facet_number, lv );
    }
}

/// ********************************************************
/// Assembly routines for primal problem
/// ********************************************************

/// Jacobian of primal problem 
#ifdef OPT_ASM0

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const
{
#    ifdef AUGMENT_PRESS
    const int t_var = DIM + 2;
#    else
    const int t_var = DIM + 1;
#    endif

    const DataType dt = this->dT_pc_;

    std::vector<DataType> grav;
    grav.resize ( 3, 0.0 );
    DataType sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    DataType skew_fac = 1.;
    DataType tau_temp_supg = 0.;

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    if ( this->temp_supg_mode_ > 0 )
        this->compute_tau_temp_supg ( element, tau_temp_supg );

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1 / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::fabs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_rr, 1., 1. };

#    ifdef ROTATING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;
#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif  

        // ***********************************************************************
        // Boussinesq forcing quadrant: theta_mom_bou_1 * dT * alpha * int{(t*g,v)} 
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, t_var ) ) += wq
                            * dt
                            * this->theta_mom_buo_c_
                            * this->alpha_g_
                            * grav[u_var]
                            * this->phi ( i, q, u_var )
                            * this->phi ( j, q, t_var )
                            * dJ;
                }
            }
        }

        // ***********************************************************************
        // **********************************************************************
        // NOW THE TEMPERATURE PARTS ...

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {

                // ***********************************************************************
                // TIME-DERIVATIVE: theta_d_dt_u * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * this->theta_d_dt_u_
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var )
                        * dJ;

                // ***********************************************************************
                // Thermal diffusion: theta1 * kappa * dT * int{grad{u}:grad{v})
                DataType tmp = inv_rr * this->grad_phi ( i, q, t_var )[0] * this->grad_phi ( j, q, t_var )[0]
                        + this->grad_phi ( i, q, t_var )[1] * this->grad_phi ( j, q, t_var )[1];
                if ( DIM == 3 ) tmp += this->grad_phi ( i, q, t_var )[2] * this->grad_phi ( j, q, t_var )[2];

                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * dt
                        * this->theta_heat_diff_c_
                        * this->kappa_
                        * tmp
                        * dJ;
            }
        }

        // Nonlinear advection part one: linearized velocity: (0.5 *) theta_heat_adv1 * dT * int{(u_c * grad{T}, v)} 
        //                                                     + (0.5) * theta_heat_adv1 * dT * int{(u_p * grad{T}, v)}
        for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
        {
            DataType tmp = 0.0;
            for ( int s = 0; s < DIM; ++s )
            {
                tmp += ( this->theta_heat_adv_cc_ * this->solP_[s][q] + this->theta_heat_adv_pc_ * this->solP_prev_[s][q] )
                        * inv_r_comp[s] * this->grad_phi ( j, q, t_var )[s];
            }
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * skew_fac
                        * dt
                        * tmp
                        * this->phi ( i, q, t_var )
                        * dJ;
            }
        }

        // Nonlinear advection part three: linearized temp: (0.5 *) theta_heat_adv1 * dT * int{(u * grad{T_c}, v)} 
        //                                                   + (0.5 *) theta_heat_adv3 * dT * int{(u * grad{T_p}, v)}    
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {
                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += wq
                            * skew_fac
                            * dt
                            * ( this->theta_heat_adv_cc_ * this->grad_solP_[t_var][q][u_var]
                            + this->theta_heat_adv_cp_ * this->grad_solP_prev_[t_var][q][u_var] )
                            * this->phi ( j, q, u_var )
                            * inv_r_comp[u_var]
                            * this->phi ( i, q, t_var )
                            * dJ;
                }
            }
        }
        if ( this->skew_mode_ > 0 )
        {
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                DataType tmp = 0.;
                for ( int s = 0; s < DIM; ++s )
                {
                    tmp += inv_r_comp[s] * this->grad_phi ( i, q, t_var )[s]
                            * ( this->theta_heat_adv_cc_ * this->solP_[s][q] + this->theta_heat_adv_pc_ * this->solP_prev_[s][q] );
                }
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq
                            * skew_fac
                            * dt
                            * tmp
                            * this->phi ( j, q, t_var )
                            * dJ;
                }
            }

            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {
                    for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
                    {
                        lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += -wq
                                * skew_fac
                                * dt
                                * inv_r_comp[u_var]
                                * this->phi ( j, q, u_var )
                                * this->grad_phi ( i, q, t_var )[u_var]
                                * ( this->theta_heat_adv_cc_ * this->solP_[t_var][q] + this->theta_heat_adv_cp_ * this->solP_prev_[t_var][q] )
                                * dJ;
                    }
                }
            }
        }

        // SUPG stabilization 
        if ( this->temp_supg_mode_ > 0 && this->theta_d_dt_u_ > 0 ) // don't use SUPG in stationary case
        {
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                // tmp_i = u_p * grad(v)
                DataType tmp_i = 0.;
                for ( int var = 0; var < DIM; ++var )
                {
                    tmp_i += this->solP_prev_[var][q] * inv_r_comp[var] * grad_phi ( i, q, t_var )[var];
                }
                //Time Derivative part 
                // tau * theta_d_dt_u * int{(T, u_p * grad(v))}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += wq
                            * tau_temp_supg
                            * this->theta_d_dt_u_
                            * phi ( j, q, t_var )
                            * tmp_i
                            * dJ;
                }

                //THERMAL DIFFUSION: - tau * theta_heat_supg_c * kappa * dT * int{Laplace{T},u_p*grad{v}}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    DataType laplace_T = 0;
                    for ( int d = 0; d < DIM; ++d )
                    {
                        laplace_T += inv_rr_comp[d] * this->H_phi ( j, q, t_var )( d, d );
                    }
                    laplace_T += inv_r * grad_phi ( j, q, t_var )[1];

                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += -wq
                            * tau_temp_supg
                            * this->kappa_
                            * dt
                            * this->theta_heat_supg_c_
                            * laplace_T
                            * tmp_i
                            * dJ;
                }

                //LINEARIZED PART
                // Nonlinear advection part one: linearized velocity: tau * theta_heat_supg_c * dT * int{(u_c * grad{T},u_p * grad{v})}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    DataType tmp_j = 0;
                    for ( int var = 0; var < DIM; ++var )
                    {
                        tmp_j += this->solP_[var][q] * inv_r_comp[var] * grad_phi ( j, q, t_var )[var];
                    }
                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += wq
                            * tau_temp_supg
                            * dt
                            * this->theta_heat_supg_c_
                            * tmp_j
                            * tmp_i
                            * dJ;
                }

                // Nonlinear advection part three: linearized temp: tau * theta_heat_supg_c_ * dT * int{(u * grad{T_c},u_p * grad{v})}
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, t_var ), dof_index ( j, u_var ) ) += wq
                                * tau_temp_supg
                                * dt
                                * this->theta_heat_supg_c_
                                * phi ( j, q, u_var )
                                * inv_r_comp[u_var]
                                * this->grad_solP_[t_var][q][u_var]
                                * tmp_i
                                * dJ;
                    }
                }
            }
        }
    }
}
#endif
#ifdef OPT_ASM1

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const
{
#    ifdef AUGMENT_PRESS
    const int t_var = DIM + 2;
#    else
    const int t_var = DIM + 1;
#    endif

    const DataType dt = this->dT_pc_;

    std::vector<DataType> grav;
    grav.resize ( 3, 0.0 );
    DataType sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    DataType skew_fac = 1.;
    DataType tau_temp_supg = 0.;

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    if ( this->temp_supg_mode_ > 0 )
        this->compute_tau_temp_supg ( element, tau_temp_supg );

    const DataType mom_bou_1 = dt * this->theta_mom_buo_c_ * this->alpha_g_;
    const DataType heat_diff_c = dt * this->theta_heat_diff_c_ * this->kappa_;
    const DataType heat_adv_cc = this->theta_heat_adv_cc_ * skew_fac * dt;
    const DataType heat_adv_pc = this->theta_heat_adv_pc_ * skew_fac * dt;
    const DataType heat_adv_cp = this->theta_heat_adv_cp_ * skew_fac * dt;
    const DataType heat_supg_c = dt * this->theta_heat_supg_c_ * tau_temp_supg;

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1 / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq_dJ = this->w ( q ) * r * std::fabs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_rr, 1., 1. };

#    ifdef ROTATING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;
#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif  

        // ***********************************************************************
        // Boussinesq forcing quadrant: theta_mom_bou_1 * dT * alpha * int{(t*g,v)} 
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                            * mom_bou_1
                            * grav[u_var]
                            * this->phi ( i, q, u_var )
                            * this->phi ( j, q, t_var );
                }
            }
        }

        // ***********************************************************************
        // **********************************************************************
        // NOW THE TEMPERATURE PARTS ...   
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {
                // ***********************************************************************
                // TIME-DERIVATIVE: theta_d_dt_u * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * this->theta_d_dt_u_
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var );

                // ***********************************************************************
                // Thermal diffusion: theta1 * kappa * dT * int{grad{u}:grad{v})
                DataType tmp = inv_rr * this->grad_phi ( i, q, t_var )[0] * this->grad_phi ( j, q, t_var )[0]
                        + this->grad_phi ( i, q, t_var )[1] * this->grad_phi ( j, q, t_var )[1];
                if ( DIM == 3 ) tmp += this->grad_phi ( i, q, t_var )[2] * this->grad_phi ( j, q, t_var )[2];

                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * heat_diff_c
                        * tmp;
            }
        }

        // Nonlinear advection part one: linearized velocity: (0.5 *) theta_heat_adv1 * dT * int{(u_c * grad{T}, v)} + (0.5) * theta_heat_adv1 * dT * int{(u_p * grad{T}, v)}
        for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
        {
            DataType tmp = 0.0;
            for ( int s = 0; s < DIM; ++s )
                tmp += ( heat_adv_cc * this->solP_[s][q] + heat_adv_pc * this->solP_prev_[s][q] ) * inv_r_comp[s] * this->grad_phi ( j, q, t_var )[s];

            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * tmp
                        * this->phi ( i, q, t_var );
            }
        }

        // Nonlinear advection part three: linearized temp: (0.5 *) theta_heat_adv1 * dT * int{(u * grad{T_c}, v)}  + (0.5 *) theta_heat_adv3 * dT * int{(u * grad{T_p}, v)}    
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {
                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * ( heat_adv_cc * this->grad_solP_[t_var][q][u_var]
                            + heat_adv_cp * this->grad_solP_prev_[t_var][q][u_var] )
                            * this->phi ( j, q, u_var )
                            * inv_r_comp[u_var]
                            * this->phi ( i, q, t_var );
                }
            }
        }
        if ( this->skew_mode_ > 0 )
        {
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                DataType tmp = 0.;
                for ( int s = 0; s < DIM; ++s )
                    tmp += inv_r_comp[s] * this->grad_phi ( i, q, t_var )[s] * ( heat_adv_cc * this->solP_[s][q] + heat_adv_pc * this->solP_prev_[s][q] );

                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq_dJ
                            * tmp
                            * this->phi ( j, q, t_var );
                }
            }

            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {
                    for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
                    {
                        lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                                * inv_r_comp[u_var]
                                * this->phi ( j, q, u_var )
                                * this->grad_phi ( i, q, t_var )[u_var]
                                * ( heat_adv_cc * this->solP_[t_var][q] + heat_adv_cp * this->solP_prev_[t_var][q] );
                    }
                }
            }
        }

        // SUPG stabilization 
        if ( this->temp_supg_mode_ > 0 && this->theta_d_dt_u_ > 0 ) // don't use SUPG in stationary case
        {
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                // tmp_i = u_p * grad(v)
                DataType tmp_i = 0.;
                for ( int var = 0; var < DIM; ++var )
                {
                    tmp_i += this->solP_prev_[var][q] * inv_r_comp[var] * grad_phi ( i, q, t_var )[var];
                }
                //Time Derivative part 
                // tau * theta_d_dt_u * int{(T, u_p * grad(v))}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += wq_dJ
                            * tau_temp_supg
                            * this->theta_d_dt_u_
                            * phi ( j, q, t_var )
                            * tmp_i;
                }

                //THERMAL DIFFUSION: - tau * theta_heat_supg_c * kappa * dT * int{Laplace{T},u_p*grad{v}}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    DataType laplace_T = 0;
                    for ( int d = 0; d < DIM; ++d )
                    {
                        laplace_T += inv_rr_comp[d] * this->H_phi ( j, q, t_var )( d, d );
                    }
                    laplace_T += inv_r * grad_phi ( j, q, t_var )[1];

                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += -wq_dJ
                            * this->kappa_
                            * heat_supg_c
                            * laplace_T
                            * tmp_i;
                }

                //LINEARIZED PART
                // Nonlinear advection part one: linearized velocity: tau * theta_heat_supg_c * dT * int{(u_c * grad{T},u_p * grad{v})}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    DataType tmp_j = 0;
                    for ( int var = 0; var < DIM; ++var )
                    {
                        tmp_j += this->solP_[var][q] * inv_r_comp[var] * grad_phi ( j, q, t_var )[var];
                    }
                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += wq_dJ
                            * heat_supg_c
                            * tmp_j
                            * tmp_i;
                }

                // Nonlinear advection part three: linearized temp: tau * theta_heat_supg_c_ * dT * int{(u * grad{T_c},u_p * grad{v})}
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, t_var ), dof_index ( j, u_var ) ) += wq_dJ
                                * heat_supg_c
                                * phi ( j, q, u_var )
                                * inv_r_comp[u_var]
                                * this->grad_solP_[t_var][q][u_var]
                                * tmp_i;
                    }
                }
            }
        }
    }
}
#endif
#ifdef OPT_ASM2

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const
{
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#    else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#    endif

    const DataType dt = this->dT_pc_;

    std::vector<DataType> grav;
    grav.resize ( 3, 0.0 );
    DataType sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    DataType skew_fac = 1.;
    DataType tau_temp_supg = 0.;

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    if ( this->temp_supg_mode_ > 0 )
        this->compute_tau_temp_supg ( element, tau_temp_supg );

    const DataType mom_bou_1 = dt * this->theta_mom_buo_c_ * this->alpha_g_;
    const DataType heat_diff_c = dt * this->theta_heat_diff_c_ * this->kappa_;
    const DataType heat_adv_cc = this->theta_heat_adv_cc_ * skew_fac * dt;
    const DataType heat_adv_pc = this->theta_heat_adv_pc_ * skew_fac * dt;
    const DataType heat_adv_cp = this->theta_heat_adv_cp_ * skew_fac * dt;
    const DataType heat_supg_c = dt * this->theta_heat_supg_c_ * tau_temp_supg;

    std::vector< std::vector<DataType> > phi;
    std::vector< std::vector< std::vector<DataType> > > grad_phi;
    std::vector< std::vector< Mat<DIM, DIM, DataType> > > H_phi;

    phi.resize ( num_vars );
    grad_phi.resize ( num_vars );

    for ( int k = 0; k < num_vars; ++k )
    {
        phi[k].resize ( this->num_dofs ( k ), 0. );
        grad_phi[k].resize ( this->num_dofs ( k ) );
        if ( k != p0_var )
        {
            for ( int j = 0; j<this->num_dofs ( k ); ++j )
            {
                grad_phi[k][j].resize ( DIM, 0. );
            }
        }
    }
    if ( this->temp_supg_mode_ > 0 )
    {
        H_phi.resize ( num_vars );
        H_phi[t_var].resize ( this->num_dofs ( t_var ) );
    }

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1 / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq_dJ = this->w ( q ) * r * std::fabs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_rr, 1., 1. };

        // PERF_TODO: check which phi's are necessary
        for ( int k = 0; k < num_vars; ++k )
        {
            for ( int j = 0; j<this->num_dofs ( k ); ++j )
            {
                phi[k][j] = this->phi ( j, q, k );
                if ( k != p0_var )
                {
                    for ( int d = 0; d < DIM; ++d )
                    {
                        grad_phi[k][j][d] = this->grad_phi ( j, q, k )[d];
                    }
                }
            }
        }
        if ( this->temp_supg_mode_ > 0 )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {
                H_phi[t_var][j] = this->H_phi ( j, q, t_var );
            }
        }

#    ifdef ROTATING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;
#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif  

        // ***********************************************************************
        // Boussinesq forcing quadrant: theta_mom_bou_1 * dT * alpha * int{(t*g,v)} 
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                            * mom_bou_1
                            * grav[u_var]
                            * phi[u_var][i]
                            * phi[t_var][j];
                }
            }
        }

        // ***********************************************************************
        // **********************************************************************
        // NOW THE TEMPERATURE PARTS ...    
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {

                // ***********************************************************************
                // TIME-DERIVATIVE: theta_d_dt_u * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * this->theta_d_dt_u_
                        * phi[t_var][j]
                        * phi[t_var][i];

                // ***********************************************************************
                // Thermal diffusion: theta1 * kappa * dT * int{grad{u}:grad{v})
                DataType tmp = inv_rr * grad_phi[t_var][i][0] * grad_phi[t_var][j][0]
                        + grad_phi[t_var][i][1] * grad_phi[t_var][j][1];
                if ( DIM == 3 ) tmp += grad_phi[t_var][i][2] * grad_phi[t_var][j][2];

                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * heat_diff_c
                        * tmp;
            }
        }

        // Nonlinear advection part one: linearized velocity: (0.5 *) theta_heat_adv1 * dT * int{(u_c * grad{T}, v)} + (0.5) * theta_heat_adv1 * dT * int{(u_p * grad{T}, v)}
        for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
        {
            DataType tmp = 0.0;
            for ( int s = 0; s < DIM; ++s )
                tmp += ( heat_adv_cc * this->solP_[s][q] + heat_adv_pc * this->solP_prev_[s][q] ) * inv_r_comp[s] * grad_phi[t_var][j][s];

            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * tmp
                        * phi[t_var][i];
            }
        }

        // Nonlinear advection part three: linearized temp: (0.5 *) theta_heat_adv1 * dT * int{(u * grad{T_c}, v)}  + (0.5 *) theta_heat_adv3 * dT * int{(u * grad{T_p}, v)}    
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {
                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * ( heat_adv_cc * this->grad_solP_[t_var][q][u_var]
                            + heat_adv_cp * this->grad_solP_prev_[t_var][q][u_var] )
                            * phi[u_var][j]
                            * inv_r_comp[u_var]
                            * phi[t_var][i];
                }
            }
        }
        if ( this->skew_mode_ > 0 )
        {
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                DataType tmp = 0.;
                for ( int s = 0; s < DIM; ++s )
                    tmp += inv_r_comp[s] * grad_phi[t_var][i][s] * ( heat_adv_cc * this->solP_[s][q] + heat_adv_pc * this->solP_prev_[s][q] );

                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq_dJ
                            * tmp
                            * phi[t_var][j];
                }
            }

            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {
                    for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
                    {
                        lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                                * inv_r_comp[u_var]
                                * phi[u_var][j]
                                * grad_phi[t_var][i][u_var]
                                * ( heat_adv_cc * this->solP_[t_var][q] + heat_adv_cp * this->solP_prev_[t_var][q] );
                    }
                }
            }
        }

        // SUPG stabilization 
        if ( this->temp_supg_mode_ > 0 && this->theta_d_dt_u_ > 0 ) // don't use SUPG in stationary case
        {
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                // tmp_i = u_p * grad(v)
                DataType tmp_i = 0.;
                for ( int var = 0; var < DIM; ++var )
                {
                    tmp_i += this->solP_prev_[var][q] * inv_r_comp[var] * grad_phi[t_var][i][var];
                }
                //Time Derivative part 
                // tau * theta_d_dt_u * int{(T, u_p * grad(v))}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += wq_dJ
                            * tau_temp_supg
                            * this->theta_d_dt_u_
                            * phi[t_var][j]
                            * tmp_i;
                }

                //THERMAL DIFFUSION: - tau * theta_heat_supg_c * kappa * dT * int{Laplace{T},u_p*grad{v}}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    DataType laplace_T = 0;
                    for ( int d = 0; d < DIM; ++d )
                    {
                        laplace_T += inv_rr_comp[d] * H_phi[t_var][j]( d, d );
                    }
                    laplace_T += inv_r * grad_phi[t_var][j][1];

                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += -wq_dJ
                            * this->kappa_
                            * heat_supg_c
                            * laplace_T
                            * tmp_i;
                }

                //LINEARIZED PART
                // Nonlinear advection part one: linearized velocity: tau * theta_heat_supg_c * dT * int{(u_c * grad{T},u_p * grad{v})}
                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    DataType tmp_j = 0;
                    for ( int var = 0; var < DIM; ++var )
                    {
                        tmp_j += this->solP_[var][q] * inv_r_comp[var] * grad_phi[t_var][j][var];
                    }
                    lm ( dof_index ( i, t_var ), dof_index ( j, t_var ) ) += wq_dJ
                            * heat_supg_c
                            * tmp_j
                            * tmp_i;
                }

                // Nonlinear advection part three: linearized temp: tau * theta_heat_supg_c_ * dT * int{(u * grad{T_c},u_p * grad{v})}
                for ( int u_var = 0; u_var < DIM; ++u_var )
                {
                    for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, t_var ), dof_index ( j, u_var ) ) += wq_dJ
                                * heat_supg_c
                                * phi[u_var][j]
                                * inv_r_comp[u_var]
                                * this->grad_solP_[t_var][q][u_var]
                                * tmp_i;
                    }
                }
            }
        }
    }
}
#endif

/// Residual of primal problem (lhs - rhs) 
#ifdef OPT_ASM0

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#    else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#    endif

    DataType sign = 1.0;
    const DataType dt = this->dT_pc_;

    std::vector<DataType> grav;
    grav.resize ( 3, 0.0 );
    DataType div_c, div_p;

    DataType tau_div = 0.;
    DataType skew_fac = 1.;
    DataType tau_temp_supg = 0.;

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    if ( this->temp_supg_mode_ > 0 )
        this->compute_tau_temp_supg ( element, tau_temp_supg );

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1. / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::abs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_r*inv_r, 1., 1. };

#    ifdef ROTAING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;

#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif 

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
#    ifdef AUGMENT_PRESS    
        Vec < DIM + 3, DataType> sol_c;
        Vec < DIM + 3, DataType> sol_p;
        Vec < DIM + 3, DataType> perturb_c;
        Vec < DIM + 3, DataType> perturb_p;
#    else
        Vec < DIM + 2, DataType> sol_c;
        Vec < DIM + 2, DataType> sol_p;
        Vec < DIM + 2, DataType> perturb_c;
        Vec < DIM + 2, DataType> perturb_p;
#    endif

        // velocity    
        for ( int var = 0; var < DIM; ++var )
        {
            sol_c[var] = this->solP_[var][q];
            sol_p[var] = this->solP_prev_[var][q];
            if ( this->L2_perturb_[var] || this->H1_perturb_[var] )
            {
                perturb_c[var] = this->perturb_[var][q];
                perturb_p[var] = this->perturb_prev_[var][q];
            }
        }

        // temperature    
        sol_c[t_var] = this->solP_[t_var][q];
        sol_p[t_var] = this->solP_prev_[t_var][q];
        if ( this->L2_perturb_[t_var] )
        {
            perturb_c[t_var] = this->perturb_[t_var][q];
            perturb_p[t_var] = this->perturb_prev_[t_var][q];
        }
        Mat<DIM, DIM, DataType> hess_temp_c;
        Mat<DIM, DIM, DataType> hess_temp_p;
        if ( this->temp_supg_mode_ > 0 )
        {
            hess_temp_c = this->hess_solP_[t_var][q];
            hess_temp_p = this->hess_solP_prev_[t_var][q];
        }

        std::vector< Vec<DIM, DataType> > grad_sol_c ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_sol_p ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_perturb_c ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_perturb_p ( num_vars );

        for ( int var = 0; var < num_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                if ( var < p_var )
                {
                    grad_sol_c[var][d] = this->grad_solP_[var][q][d];
                    grad_sol_p[var][d] = this->grad_solP_prev_[var][q][d];

                    if ( this->H1_perturb_[var] )
                    {
                        grad_perturb_c[var][d] = this->grad_perturb_[var][q][d];
                        grad_perturb_p[var][d] = this->grad_perturb_prev_[var][q][d];
                    }
                }
                if ( var == t_var )
                {
                    grad_sol_c[var][d] = this->grad_solP_[t_var][q][d];
                    grad_sol_p[var][d] = this->grad_solP_prev_[t_var][q][d];

                    if ( this->H1_perturb_[var] )
                    {
                        grad_perturb_c[var][d] = this->grad_perturb_[t_var][q][d];
                        grad_perturb_p[var][d] = this->grad_perturb_prev_[t_var][q][d];
                    }
                }
            }
        }

        // ********************************************************************** 
        // BOUSSINESQ-FORCING: theta_mom_bou_1 * dT * alpha_g * int{((T_c-T-ref)*g,v)} + theta_mom_bou_2 * dT * alpha_g * int{((T_p-T-ref)*g,v)}
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += wq
                        * dt
                        * this->alpha_g_
                        * grav[v_var]
                        * this->phi ( i, q, v_var )
                        * ( ( sol_c[t_var] - this->ref_T_ ) * this->theta_mom_buo_c_ // current  time
                        + ( sol_p[t_var] - this->ref_T_ ) * this->theta_mom_buo_p_ ) // previous time
                        * dJ;
            }
        }

        // **********************************************************************
        // **********************************************************************
        // NOW THE TEMPERATURE PARTS ...

        // **********************************************************************
        // TIME DERIVATIVE: theta_d_dt_u * ((T_c - T_p),v) 

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq
                    * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                    * ( sol_c[t_var] - sol_p[t_var] ) // current - previous
                    * this->phi ( i, q, t_var )
                    * dJ;
        }

        // **********************************************************************
        // LAPLACE: theta_heat_diff_c * dT * diff * int{grad{T_c} : grad{v}} + theta_heat_diff_p * dT * diff * int{grad{T_p} : grad{v}}

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            DataType laplace_c = 0.0;
            for ( int s = 0; s < DIM; s++ ) laplace_c += inv_rr_comp[s] * this->grad_phi ( i, q, t_var )[s] * grad_sol_c[t_var][s];

            DataType laplace_p = 0.0;
            for ( int s = 0; s < DIM; s++ ) laplace_p += inv_rr_comp[s] * this->grad_phi ( i, q, t_var )[s] * grad_sol_p[t_var][s];

            lv[this->dof_index ( i, t_var )] += wq
                    * dt
                    * this->kappa_
                    * ( this->theta_heat_diff_c_ * laplace_c
                    + this->theta_heat_diff_p_ * laplace_p )
                    * dJ;
        }

        // **********************************************************************
        // CONVECTIVE TERM 

        //    (0.5)* theta_heat_adv_cc * dT * int{(u_c * grad{T_c}, v)} + theta_heat_adv_pc * dT * int{(u_p * grad{T_c}, v)}
        //  + (0.5)* theta_heat_adv_cp * dT * int{(u_c * grad{T_p}, v)} + theta_heat_adv_pp * dT * int{(u_p * grad{T_p}, v)}
        DataType convection = 0.0;
        for ( int s = 0; s < DIM; ++s )
        {
            convection += ( sol_c[s] * grad_sol_c[t_var][s] * this->theta_heat_adv_cc_
                    + sol_c[s] * grad_sol_p[t_var][s] * this->theta_heat_adv_cp_
                    + sol_p[s] * grad_sol_c[t_var][s] * this->theta_heat_adv_pc_
                    + sol_p[s] * grad_sol_p[t_var][s] * this->theta_heat_adv_pp_ )
                    * inv_r_comp[s];
        }

        for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq
                    * skew_fac
                    * dt
                    * convection
                    * this->phi ( i, q, t_var )
                    * dJ;

        }

        if ( this->skew_mode_ > 0 )
        {
            for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
            {
                DataType tmp_1 = 0.;
                DataType tmp_2 = 0.;
                DataType tmp_3 = 0.;
                DataType tmp_4 = 0.;

                for ( int s = 0; s < DIM; ++s )
                {
                    tmp_1 += inv_r_comp[s] * this->theta_heat_adv_cc_ * sol_c[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_2 += inv_r_comp[s] * this->theta_heat_adv_pc_ * sol_p[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_3 += inv_r_comp[s] * this->theta_heat_adv_cp_ * sol_c[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_4 += inv_r_comp[s] * this->theta_heat_adv_pp_ * sol_p[s] * this->grad_phi ( i, q, t_var )[s];
                }
                lv[this->dof_index ( i, t_var )] += -wq
                        * skew_fac
                        * dt
                        * ( tmp_1 * sol_c[t_var] +
                        tmp_2 * sol_c[t_var] +
                        tmp_3 * sol_p[t_var] +
                        tmp_4 * sol_p[t_var] )
                        * dJ;
            }
        }

        // SUPG stabilization
        if ( this->temp_supg_mode_ > 0 && this->theta_d_dt_u_ > 0 )
        {
            // TIME DERIVATIVE: tau * theta_d_dt_u * ((T_c - T_p),u_p*grad{v}) 
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                // tmp_i = u_p * grad{v}
                DataType tmp_i = 0.;
                for ( int var = 0; var < DIM; ++var )
                {
                    tmp_i += sol_p[var] * inv_r_comp[var] * grad_phi ( i, q, t_var )[var];
                }
                lv[this->dof_index ( i, t_var )] += wq
                        * tau_temp_supg
                        * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                        * ( sol_c[t_var] - sol_p[t_var] ) // current - previous
                        * tmp_i
                        * dJ;

                // **********************************************************************
                // LAPLACE: - tau * dT * kappa * (theta_heat_supg_c * int{Laplace{T_c}, u_p*grad{v}} - theta_heat_supg_p * int{Laplace{T_p}, u_p*grad{v}})
                DataType laplace_c = 0.0;
                for ( int d = 0; d < DIM; d++ )
                    laplace_c += inv_rr_comp[d] * hess_temp_c ( d, d );

                laplace_c += inv_r * grad_sol_c[t_var][1];

                DataType laplace_p = 0.0;
                for ( int d = 0; d < DIM; d++ )
                    laplace_p += inv_rr_comp[d] * hess_temp_p ( d, d );

                laplace_p += inv_r * grad_sol_p[t_var][1];

                lv[this->dof_index ( i, t_var )] += -wq
                        * tau_temp_supg
                        * this->kappa_
                        * dt
                        * ( this->theta_heat_supg_c_ * laplace_c
                        + this->theta_heat_supg_p_ * laplace_p )
                        * tmp_i
                        * dJ;

                //  tau * dT * ( theta_heat_supg_c * int{(u_c * grad{T_c}, u_p * grad{v})} + theta_heat_supg_p  * int{(u_p * grad{T_p}, u_p * grad{v})})
                DataType convection_tmp = 0.0;
                for ( int d = 0; d < DIM; ++d )
                {
                    convection_tmp += ( sol_c[d] * grad_sol_c[t_var][d] * this->theta_heat_supg_c_ // current
                            + sol_p[d] * grad_sol_p[t_var][d] * this->theta_heat_supg_p_ )
                            * inv_r_comp[d];
                }

                lv[this->dof_index ( i, t_var )] += wq
                        * tau_temp_supg
                        * dt
                        * convection_tmp
                        * tmp_i
                        * dJ;

            }
        }

        // perturbation 
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            DataType tmp_L2_p = 0.;
            DataType tmp_L2_c = 0.;
            DataType tmp_H1_p = 0.;
            DataType tmp_H1_c = 0.;

            if ( this->H1_perturb_[t_var] )
            {
                for ( int s = 0; s < DIM; s++ )
                {
                    tmp_H1_c += inv_rr_comp[s] * this->grad_phi ( i, q, t_var )[s] * grad_perturb_c[t_var][s];
                    tmp_H1_p += inv_rr_comp[s] * this->grad_phi ( i, q, t_var )[s] * grad_perturb_p[t_var][s];
                }
            }
            if ( this->L2_perturb_[t_var] )
            {
                tmp_L2_c += this->phi ( i, q, t_var ) * perturb_c[t_var];
                tmp_L2_p += this->phi ( i, q, t_var ) * perturb_p[t_var];
            }

            lv[this->dof_index ( i, t_var )] += wq
                    * dt
                    * 0.5
                    * this->perturb_scale_
                    * ( tmp_H1_c + tmp_H1_p + tmp_L2_c + tmp_L2_p )
                    * dJ;
        }
    }
}
#endif
#ifdef OPT_ASM1

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#    else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#    endif

    DataType sign = 1.0;

    std::vector<DataType> grav;
    grav.resize ( 3, 0.0 );
    DataType div_c, div_p;
    const DataType dt = this->dT_pc_;

    DataType skew_fac = 1.;
    DataType tau_temp_supg = 0.;

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    if ( this->temp_supg_mode_ > 0 )
        this->compute_tau_temp_supg ( element, tau_temp_supg );

    const DataType mom_buo_c = dt * this->theta_mom_buo_c_ * this->alpha_g_;
    const DataType mom_buo_p = dt * this->theta_mom_buo_p_ * this->alpha_g_;
    const DataType heat_diff_c = dt * this->theta_heat_diff_c_ * this->kappa_;
    const DataType heat_diff_p = dt * this->theta_heat_diff_p_ * this->kappa_;
    const DataType heat_adv_cc = this->theta_heat_adv_cc_ * skew_fac * dt;
    const DataType heat_adv_pc = this->theta_heat_adv_pc_ * skew_fac * dt;
    const DataType heat_adv_cp = this->theta_heat_adv_cp_ * skew_fac * dt;
    const DataType heat_adv_pp = this->theta_heat_adv_pp_ * skew_fac * dt;
    const DataType heat_supg_c = dt * this->theta_heat_supg_c_ * tau_temp_supg;
    const DataType heat_supg_p = dt * this->theta_heat_supg_p_ * tau_temp_supg;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1. / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq_dJ = this->w ( q ) * r * std::abs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_r*inv_r, 1., 1. };

#    ifdef ROTAING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;

#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif 

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
#    ifdef AUGMENT_PRESS    
        Vec < DIM + 3, DataType> sol_c;
        Vec < DIM + 3, DataType> sol_p;
        Vec < DIM + 3, DataType> perturb_c;
        Vec < DIM + 3, DataType> perturb_p;
#    else
        Vec < DIM + 2, DataType> sol_c;
        Vec < DIM + 2, DataType> sol_p;
        Vec < DIM + 2, DataType> perturb_c;
        Vec < DIM + 2, DataType> perturb_p;
#    endif

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
        }

        // temperature    
        sol_c[t_var] = this->solP_[t_var][q];
        sol_p[t_var] = this->solP_prev_[t_var][q];
        if ( this->L2_perturb_[t_var] )
        {
            perturb_c[t_var] = this->perturb_[t_var][q];
            perturb_p[t_var] = this->perturb_prev_[t_var][q];
        }
        Mat<DIM, DIM, DataType> hess_temp_c;
        Mat<DIM, DIM, DataType> hess_temp_p;
        if ( this->temp_supg_mode_ > 0 )
        {
            hess_temp_c = this->hess_solP_[t_var][q];
            hess_temp_p = this->hess_solP_prev_[t_var][q];
        }

        std::vector< Vec<DIM, DataType> > grad_sol_c ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_sol_p ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_perturb_c ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_perturb_p ( num_vars );

        for ( int var = 0; var < num_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                if ( var < p_var )
                {
                    grad_sol_c[var][d] = this->grad_solP_[var][q][d];
                    grad_sol_p[var][d] = this->grad_solP_prev_[var][q][d];

                    if ( this->H1_perturb_[var] )
                    {
                        grad_perturb_c[var][d] = this->grad_perturb_vel_[var][q][d];
                        grad_perturb_p[var][d] = this->grad_perturb_vel_prev_[var][q][d];
                    }
                }
                if ( var == t_var )
                {
                    grad_sol_c[var][d] = this->grad_solP_[t_var][q][d];
                    grad_sol_p[var][d] = this->grad_solP_prev_[t_var][q][d];

                    if ( this->H1_perturb_[var] )
                    {
                        grad_perturb_c[var][d] = this->grad_perturb_[t_var][q][d];
                        grad_perturb_p[var][d] = this->grad_perturb_prev_[t_var][q][d];
                    }
                }
            }
        }

        // ********************************************************************** 
        // BOUSSINESQ-FORCING: theta_mom_bou_1 * dT * alpha_g * int{((T_c-T-ref)*g,v)} + theta_mom_bou_2 * dT * alpha_g * int{((T_p-T-ref)*g,v)}

        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * grav[v_var]
                        * this->phi ( i, q, v_var )
                        * ( ( sol_c[t_var] - this->ref_T_ ) * mom_buo_c // current  time
                        + ( sol_p[t_var] - this->ref_T_ ) * mom_buo_p ); // previous time

            }
        }

        // **********************************************************************
        // **********************************************************************
        // NOW THE TEMPERATURE PARTS ...

        // **********************************************************************
        // TIME DERIVATIVE: theta_d_dt_u * ((T_c - T_p),v) 

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                    * ( sol_c[t_var] - sol_p[t_var] ) // current - previous
                    * this->phi ( i, q, t_var );
        }

        // **********************************************************************
        // LAPLACE: theta_heat_diff_c * dT * diff * int{grad{T_c} : grad{v}} + theta_heat_diff_p * dT * diff * int{grad{T_p} : grad{v}}

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            DataType laplace_c = 0.0;
            for ( int s = 0; s < DIM; s++ )
                laplace_c += inv_rr_comp[s] * this->grad_phi ( i, q, t_var )[s] * grad_sol_c[t_var][s];

            DataType laplace_p = 0.0;
            for ( int s = 0; s < DIM; s++ )
                laplace_p += inv_rr_comp[s] * this->grad_phi ( i, q, t_var )[s] * grad_sol_p[t_var][s];

            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * ( heat_diff_c * laplace_c
                    + heat_diff_p * laplace_p );
        }

        // **********************************************************************
        // CONVECTIVE TERM 

        //    (0.5)* theta_heat_adv_cc * dT * int{(u_c * grad{T_c}, v)} + theta_heat_adv_pc * dT * int{(u_p * grad{T_c}, v)}
        //  + (0.5)* theta_heat_adv_cp * dT * int{(u_c * grad{T_p}, v)} + theta_heat_adv_pp * dT * int{(u_p * grad{T_p}, v)}
        DataType convection = 0.0;
        for ( int s = 0; s < DIM; ++s )
        {
            convection += ( sol_c[s] * grad_sol_c[t_var][s] * heat_adv_cc
                    + sol_c[s] * grad_sol_p[t_var][s] * heat_adv_cp
                    + sol_p[s] * grad_sol_c[t_var][s] * heat_adv_pc
                    + sol_p[s] * grad_sol_p[t_var][s] * heat_adv_pp )
                    * inv_r_comp[s];
        }

        for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * convection
                    * this->phi ( i, q, t_var );

        }

        if ( this->skew_mode_ > 0 )
        {
            for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
            {
                DataType tmp_1 = 0.;
                DataType tmp_2 = 0.;
                DataType tmp_3 = 0.;
                DataType tmp_4 = 0.;

                for ( int s = 0; s < DIM; ++s )
                {
                    tmp_1 += inv_r_comp[s] * sol_c[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_2 += inv_r_comp[s] * sol_p[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_3 += inv_r_comp[s] * sol_c[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_4 += inv_r_comp[s] * sol_p[s] * this->grad_phi ( i, q, t_var )[s];
                }
                lv[this->dof_index ( i, t_var )] += -wq_dJ
                        * ( heat_adv_cc * tmp_1 * sol_c[t_var] +
                        heat_adv_pc * tmp_2 * sol_c[t_var] +
                        heat_adv_cp * tmp_3 * sol_p[t_var] +
                        heat_adv_pp * tmp_4 * sol_p[t_var] );
            }
        }
        // SUPG stabilization
        if ( this->temp_supg_mode_ > 0 && this->theta_d_dt_u_ > 0 )
        {

            // TIME DERIVATIVE: tau * theta_d_dt_u * ((T_c - T_p),u_p*grad{v}) 
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                // tmp_i = u_p * grad{v}
                DataType tmp_i = 0.;
                for ( int var = 0; var < DIM; ++var )
                {
                    tmp_i += sol_p[var] * inv_r_comp[var] * grad_phi ( i, q, t_var )[var];
                }
                lv[this->dof_index ( i, t_var )] += wq_dJ
                        * tau_temp_supg
                        * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                        * ( sol_c[t_var] - sol_p[t_var] ) // current - previous
                        * tmp_i;

                // **********************************************************************
                // LAPLACE: - tau * dT * kappa * (theta_heat_supg_c * int{Laplace{T_c}, u_p*grad{v}} - theta_heat_supg_p * int{Laplace{T_p}, u_p*grad{v}})
                DataType laplace_c = 0.0;
                for ( int d = 0; d < DIM; d++ )
                    laplace_c += inv_rr_comp[d] * hess_temp_c ( d, d );

                laplace_c += inv_r * grad_sol_c[t_var][1];

                DataType laplace_p = 0.0;
                for ( int d = 0; d < DIM; d++ )
                    laplace_p += inv_rr_comp[d] * hess_temp_p ( d, d );

                laplace_p += inv_r * grad_sol_p[t_var][1];

                lv[this->dof_index ( i, t_var )] += -wq_dJ
                        * this->kappa_
                        * ( heat_supg_c * laplace_c
                        + heat_supg_p * laplace_p )
                        * tmp_i;

                //  tau * dT * ( theta_heat_supg_c * int{(u_c * grad{T_c}, u_p * grad{v})} + theta_heat_supg_p  * int{(u_p * grad{T_p}, u_p * grad{v})})
                DataType convection_tmp = 0.0;
                for ( int d = 0; d < DIM; ++d )
                {
                    convection_tmp += ( sol_c[d] * grad_sol_c[t_var][d] * heat_supg_c // current
                            + sol_p[d] * grad_sol_p[t_var][d] * heat_supg_p )
                            * inv_r_comp[d];
                }

                lv[this->dof_index ( i, t_var )] += wq_dJ
                        * convection_tmp
                        * tmp_i;

            }
        }
    }
}
#endif
#ifdef OPT_ASM2

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#    else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#    endif
    DataType sign = 1.0;

    std::vector<DataType> grav;
    grav.resize ( 3, 0.0 );
    DataType div_c, div_p;
    const DataType dt = this->dT_pc_;

    DataType tau_div = 0.;
    DataType skew_fac = 1.;
    DataType tau_temp_supg = 0.;

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    if ( this->temp_supg_mode_ > 0 )
        this->compute_tau_temp_supg ( element, tau_temp_supg );

    const DataType mom_buo_c = dt * this->theta_mom_buo_c_ * this->alpha_g_;
    const DataType mom_buo_p = dt * this->theta_mom_buo_p_ * this->alpha_g_;
    const DataType heat_diff_c = dt * this->theta_heat_diff_c_ * this->kappa_;
    const DataType heat_diff_p = dt * this->theta_heat_diff_p_ * this->kappa_;
    const DataType heat_adv_cc = this->theta_heat_adv_cc_ * skew_fac * dt;
    const DataType heat_adv_pc = this->theta_heat_adv_pc_ * skew_fac * dt;
    const DataType heat_adv_cp = this->theta_heat_adv_cp_ * skew_fac * dt;
    const DataType heat_adv_pp = this->theta_heat_adv_pp_ * skew_fac * dt;
    const DataType heat_supg_c = dt * this->theta_heat_supg_c_ * tau_temp_supg;
    const DataType heat_supg_p = dt * this->theta_heat_supg_p_ * tau_temp_supg;

    std::vector< std::vector<DataType> > phi;
    std::vector< std::vector< std::vector<DataType> > > grad_phi;
    std::vector< std::vector< Mat<DIM, DIM, DataType> > > H_phi;

    phi.resize ( num_vars );
    grad_phi.resize ( num_vars );

    for ( int k = 0; k < num_vars; ++k )
    {
        phi[k].resize ( this->num_dofs ( k ), 0. );
        grad_phi[k].resize ( this->num_dofs ( k ) );
        for ( int j = 0; j<this->num_dofs ( k ); ++j )
        {
            grad_phi[k][j].resize ( DIM, 0. );
        }
    }
    if ( this->temp_supg_mode_ > 0 )
    {
        H_phi.resize ( num_vars );
        H_phi[t_var].resize ( this->num_dofs ( t_var ) );
    }

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1. / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq_dJ = this->w ( q ) * r * std::abs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_r*inv_r, 1., 1. };

#    ifdef ROTAING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;

#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif 

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
#    ifdef AUGMENT_PRESS    
        Vec < DIM + 3, DataType> sol_c;
        Vec < DIM + 3, DataType> sol_p;
        Vec < DIM + 3, DataType> perturb_c;
        Vec < DIM + 3, DataType> perturb_p;
#    else
        Vec < DIM + 2, DataType> sol_c;
        Vec < DIM + 2, DataType> sol_p;
        Vec < DIM + 2, DataType> perturb_c;
        Vec < DIM + 2, DataType> perturb_p;
#    endif

        // velocity    
        for ( int var = 0; var < DIM; ++var )
        {
            sol_c[var] = this->solP_[var][q];
            sol_p[var] = this->solP_prev_[var][q];
            if ( this->L2_perturb_[var] || this->H1_perturb_[var] )
            {
                perturb_c[var] = this->perturb_[var][q];
                perturb_p[var] = this->perturb_prev_[var][q];
            }
        }

        // temperature    
        sol_c[t_var] = this->solP_[t_var][q];
        sol_p[t_var] = this->solP_prev_[t_var][q];
        if ( this->L2_perturb_[t_var] )
        {
            perturb_c[t_var] = this->perturb_[t_var][q];
            perturb_p[t_var] = this->perturb_prev_[t_var][q];
        }
        Mat<DIM, DIM, DataType> hess_temp_c;
        Mat<DIM, DIM, DataType> hess_temp_p;
        if ( this->temp_supg_mode_ > 0 )
        {
            hess_temp_c = this->hess_solP_[t_var][q];
            hess_temp_p = this->hess_solP_prev_[t_var][q];
        }

        std::vector< Vec<DIM, DataType> > grad_sol_c ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_sol_p ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_perturb_c ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_perturb_p ( num_vars );

        for ( int var = 0; var < num_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                if ( var < p_var )
                {
                    grad_sol_c[var][d] = this->grad_solP_[var][q][d];
                    grad_sol_p[var][d] = this->grad_solP_prev_[var][q][d];

                    if ( this->H1_perturb_[var] )
                    {
                        grad_perturb_c[var][d] = this->grad_perturb_[var][q][d];
                        grad_perturb_p[var][d] = this->grad_perturb_prev_[var][q][d];
                    }
                }
                if ( var == t_var )
                {
                    grad_sol_c[var][d] = this->grad_solP_[t_var][q][d];
                    grad_sol_p[var][d] = this->grad_solP_prev_[t_var][q][d];

                    if ( this->H1_perturb_[var] )
                    {
                        grad_perturb_c[var][d] = this->grad_perturb_[t_var][q][d];
                        grad_perturb_p[var][d] = this->grad_perturb_prev_[t_var][q][d];
                    }
                }
            }
        }

        // shape functions
        for ( int k = 0; k < num_vars; ++k )
        {
            if ( k < p_var || k == t_var )
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
        }

        if ( this->temp_supg_mode_ > 0 )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {
                H_phi[t_var][j] = this->H_phi ( j, q, t_var );
            }
        }

        // ********************************************************************** 
        // BOUSSINESQ-FORCING: theta_mom_bou_1 * dT * alpha_g * int{((T_c-T-ref)*g,v)} + theta_mom_bou_2 * dT * alpha_g * int{((T_p-T-ref)*g,v)}      
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * grav[v_var]
                        * phi[v_var][i]
                        * ( ( sol_c[t_var] - this->ref_T_ ) * mom_buo_c // current  time
                        + ( sol_p[t_var] - this->ref_T_ ) * mom_buo_p ); // previous time

            }
        }

        // **********************************************************************
        // **********************************************************************
        // NOW THE TEMPERATURE PARTS ...

        // **********************************************************************
        // TIME DERIVATIVE: theta_d_dt_u * ((T_c - T_p),v) 

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                    * ( sol_c[t_var] - sol_p[t_var] ) // current - previous
                    * phi[t_var][i];
        }

        // **********************************************************************
        // LAPLACE: theta_heat_diff_c * dT * diff * int{grad{T_c} : grad{v}} + theta_heat_diff_p * dT * diff * int{grad{T_p} : grad{v}}

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            DataType laplace_c = 0.0;
            for ( int s = 0; s < DIM; s++ ) laplace_c += inv_rr_comp[s] * grad_phi[t_var][i][s] * grad_sol_c[t_var][s];

            DataType laplace_p = 0.0;
            for ( int s = 0; s < DIM; s++ ) laplace_p += inv_rr_comp[s] * grad_phi[t_var][i][s] * grad_sol_p[t_var][s];

            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * ( heat_diff_c * laplace_c
                    + heat_diff_p * laplace_p );
        }

        // **********************************************************************
        // CONVECTIVE TERM 

        //    (0.5)* theta_heat_adv_cc * dT * int{(u_c * grad{T_c}, v)} + theta_heat_adv_pc * dT * int{(u_p * grad{T_c}, v)}
        //  + (0.5)* theta_heat_adv_cp * dT * int{(u_c * grad{T_p}, v)} + theta_heat_adv_pp * dT * int{(u_p * grad{T_p}, v)}
        DataType convection = 0.0;
        for ( int s = 0; s < DIM; ++s )
        {
            convection += ( sol_c[s] * grad_sol_c[t_var][s] * heat_adv_cc
                    + sol_c[s] * grad_sol_p[t_var][s] * heat_adv_cp
                    + sol_p[s] * grad_sol_c[t_var][s] * heat_adv_pc
                    + sol_p[s] * grad_sol_p[t_var][s] * heat_adv_pp )
                    * inv_r_comp[s];
        }

        for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * convection
                    * phi[t_var][i];

        }

        if ( this->skew_mode_ > 0 )
        {
            for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
            {
                DataType tmp_1 = 0.;
                DataType tmp_2 = 0.;
                DataType tmp_3 = 0.;
                DataType tmp_4 = 0.;

                for ( int s = 0; s < DIM; ++s )
                {
                    tmp_1 += inv_r_comp[s] * sol_c[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_2 += inv_r_comp[s] * sol_p[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_3 += inv_r_comp[s] * sol_c[s] * this->grad_phi ( i, q, t_var )[s];
                    tmp_4 += inv_r_comp[s] * sol_p[s] * this->grad_phi ( i, q, t_var )[s];
                }
                lv[this->dof_index ( i, t_var )] += -wq_dJ
                        * ( heat_adv_cc * tmp_1 * sol_c[t_var] +
                        heat_adv_pc * tmp_2 * sol_c[t_var] +
                        heat_adv_cp * tmp_3 * sol_p[t_var] +
                        heat_adv_pp * tmp_4 * sol_p[t_var] );
            }
        }
        // SUPG stabilization
        if ( this->temp_supg_mode_ > 0 && this->theta_d_dt_u_ > 0 )
        {

            // TIME DERIVATIVE: tau * theta_d_dt_u * ((T_c - T_p),u_p*grad{v}) 
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                // tmp_i = u_p * grad{v}
                DataType tmp_i = 0.;
                for ( int var = 0; var < DIM; ++var )
                {
                    tmp_i += sol_p[var] * inv_r_comp[var] * grad_phi[t_var][i][var];
                }
                lv[this->dof_index ( i, t_var )] += wq_dJ
                        * tau_temp_supg
                        * this->theta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                        * ( sol_c[t_var] - sol_p[t_var] ) // current - previous
                        * tmp_i;

                // **********************************************************************
                // LAPLACE: - tau * dT * kappa * (theta_heat_supg_c * int{Laplace{T_c}, u_p*grad{v}} - theta_heat_supg_p * int{Laplace{T_p}, u_p*grad{v}})
                DataType laplace_c = 0.0;
                for ( int d = 0; d < DIM; d++ )
                    laplace_c += inv_rr_comp[d] * hess_temp_c ( d, d );

                laplace_c += inv_r * grad_sol_c[t_var][1];

                DataType laplace_p = 0.0;
                for ( int d = 0; d < DIM; d++ )
                    laplace_p += inv_rr_comp[d] * hess_temp_p ( d, d );

                laplace_p += inv_r * grad_sol_p[t_var][1];

                lv[this->dof_index ( i, t_var )] += -wq_dJ
                        * this->kappa_
                        * ( heat_supg_c * laplace_c
                        + heat_supg_p * laplace_p )
                        * tmp_i;

                //  tau * dT * ( theta_heat_supg_c * int{(u_c * grad{T_c}, u_p * grad{v})} + theta_heat_supg_p  * int{(u_p * grad{T_p}, u_p * grad{v})})
                DataType convection_tmp = 0.0;
                for ( int d = 0; d < DIM; ++d )
                {
                    convection_tmp += ( sol_c[d] * grad_sol_c[t_var][d] * heat_supg_c // current
                            + sol_p[d] * grad_sol_p[t_var][d] * heat_supg_p )
                            * inv_r_comp[d];
                }

                lv[this->dof_index ( i, t_var )] += wq_dJ
                        * convection_tmp
                        * tmp_i;

            }
        }
    }
}
#endif

/// gradient of temperature

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_grad_temp ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;

#else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#endif    

    const int num_q = this->num_quadrature_points ( );

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1 / r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::fabs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };

        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            DataType grad = inv_r_comp[v_var] * this->grad_solP_prev_[t_var][q][v_var];

            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += wq
                        * grad
                        * this->phi ( i, q, v_var )
                        * dJ;
            }
        }
    }
}

/// Initial condition for dual problem

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_goal ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#endif    

    const int num_q = this->num_quadrature_points ( );

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1 / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::fabs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_rr, 1., 1. };

        ParametersFinalType<DIM, DataType> p;
        p.solP.resize ( num_vars );
        p.grad_solP.resize ( num_vars );

        for ( int var = 0; var < DIM; ++var )
            p.solP[var] = this->solP_[var][q];

        p.solP[p_var] = this->solP_[p_var][q];
#ifdef AUGMENT_PRESS
        p.solP[p0_var] = this->solP_pressure0_[q];
#endif    
        p.solP[t_var] = this->solP_[t_var][q];

        for ( int var = 0; var < grad_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                p.grad_solP[var][d] = this->grad_solP_[var][q][d];
            }
        }
        for ( int d = 0; d < DIM; d++ )
        {
            p.grad_solP[t_var][d] = this->grad_solP_[t_var][q][d];
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
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_quant ( const Element<DataType>& element, LocalVector& lv ) const
{

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#endif    
    DataType sign = 1.0;

    std::vector<DataType> grav;
    grav.resize ( 3, 0.0 );

    const int num_q = this->num_quadrature_points ( );
    const int offset = this->incomp_scalars_;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1. / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::abs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_r*inv_r, 1., 1. };

        for ( int l = 0; l<this->bous_scalars_; ++l )
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
                    // axial buoyancy
                    lv[offset + l] -= wq * this->alpha_g_ * this->g_[2] * ( this->solP_[t_var][q] - this->ref_T_ ) * this->solP_[2][q] * dJ;
                    break;
                case 1:
                    // temperature
                    lv[offset + l] += wq * this->solP_[t_var][q] * this->solP_[t_var][q] * dJ;
                    break;
                case 2:
                    // azimuthal gradient of temperature 
                    for ( int d = 0; d < DIM; d++ )
                        pot += inv_rr_comp[d] * this->grad_solP_[t_var][q][d] * this->grad_solP_[t_var][q][d];

                    lv[offset + l] += wq * pot * dJ;
                    break;
                case 3:
                    // radial heat transfer in given volume 
                    if ( r < this->nusselt_r_min_ || r > this->nusselt_r_max_ )
                        lv[offset + l] += 0.;
                    else
                        lv[offset + l] += wq * this->grad_solP_[t_var][q][1] * dJ;
                    break;
                case 4:
                    // temp difference from base state
                    if ( this->base_vector_set_ )
                    {
                        diff = ( this->base_[t_var][q] - this->solP_[t_var][q] );
                        lv[offset + l] += wq * diff * diff * dJ;
                    }
                    break;
            }
        }
    }
}

/// DEP + buoyancy force

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector_buoyancy ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#endif
    const int num_q = this->num_quadrature_points ( );
    const int total_dofs = this->num_dofs_total ( );

    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1 / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::fabs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_rr, 1., 1. };

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
        Vec<num_vars, DataType> sol_c;

        sol_c[t_var] = this->solP_[t_var][q];

        std::vector<DataType> grav;
        grav.resize ( 3, 0.0 );

#ifdef ROTATING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += this->omega_ * this->omega_ * r;

#else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#endif     

        // ***************************************************************
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += -wq
                        * this->alpha_g_
                        * grav[v_var]
                        * this->phi ( i, q, v_var )
                        * ( sol_c[t_var] - this->ref_T_ )
                        * dJ;
            }
        }
    }
}

/// squared L2- W1,2 norm of complete solution vector 

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_goal_int ( const Element<DataType>& element, DataType& ls ) const
{

    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#endif    

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1. / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::abs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_r*inv_r, 1., 1. };

        std::vector<DataType> solP_c ( num_vars );

        for ( int var = 0; var < DIM; ++var )
        {
            solP_c[var] = this->solP_[var][q];
        }

        solP_c[p_var] = this->solP_[p_var][q];
#ifdef AUGMENT_PRESS
        solP_c[p0_var] = this->solP_pressure0_[q];
#endif
        solP_c[t_var] = this->solP_[t_var][q];

        std::vector< Vec<DIM, DataType> > grad_solP_c ( num_vars );

        for ( int var = 0; var < p_var; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_solP_c[var][d] = this->grad_solP_[var][q][d];
            }
        }
        for ( int d = 0; d < DIM; d++ )
        {
            grad_solP_c[t_var][d] = this->grad_solP_[t_var][q][d];
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
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_goal_fin ( const Element<DataType>& element, DataType& ls ) const
{

    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#endif    

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1. / r;
        const DataType inv_rr = inv_r * inv_r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::abs ( this->detJ ( q ) );
        const DataType inv_r_comp[3] = { inv_r, 1., 1. };
        const DataType inv_rr_comp[3] = { inv_r*inv_r, 1., 1. };

        std::vector<DataType> solP_c ( num_vars );

        for ( int var = 0; var < DIM; ++var )
        {
            solP_c[var] = this->solP_[var][q];
        }

        solP_c[p_var] = this->solP_[p_var][q];
#ifdef AUGMENT_PRESS
        solP_c[p0_var] = this->solP_pressure0_[q];
#endif
        solP_c[t_var] = this->solP_[t_var][q];

        std::vector< Vec<DIM, DataType> > grad_solP_c ( num_vars );

        for ( int var = 0; var < p_var; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_solP_c[var][d] = this->grad_solP_[var][q][d];
            }
        }
        for ( int d = 0; d < DIM; d++ )
        {
            grad_solP_c[t_var][d] = this->grad_solP_[t_var][q][d];
        }

        ParametersEvalType<DIM, DataType> gp;
        gp.x = this->x ( q );
        gp.solP = solP_c;
        gp.grad_solP = grad_solP_c;

        ls += wq * this->goal_functional_->j_final_eval ( gp ) * dJ;
    }
}

/// mean value of temperature

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_temp_mean ( const Element<DataType>& element, DataType& ls ) const
{
#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#endif    

    const int num_q = this->num_quadrature_points ( );
    DataType sign = 1.0;

    DataType int_T = 0.0;
    DataType int_1 = 0.0;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {

        const DataType r = this->x ( q )[1];
        const DataType inv_r = 1.0 / r;
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::abs ( this->detJ ( q ) );

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index

        DataType temp = this->solP_[t_var][q];

        int_T += wq * temp * dJ;
        int_1 += wq * dJ;
    }
    ls = int_T / int_1;
}

/// heat transfer 

template<int DIM, class DataType>
void MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_rad_heat_surface ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const
{
#ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
    const int num_vars = DIM + 3;
#else
    const int p0_var = -1;
    const int t_var = DIM + 1;
    const int num_vars = DIM + 2;
#endif    

    // get material number of current facet 
    mesh::TDim tdim = DIM;
    mesh::IncidentEntityIterator iter = element.get_cell ( ).begin_incident ( tdim - 1 );
    mesh::IncidentEntityIterator end = element.get_cell ( ).end_incident ( tdim - 1 );

    int counter = 0;
    int mat_number = -1;

    for (; iter != end; iter++ )
    {
        if ( counter == facet_number )
        {
            mat_number = iter->get_material_number ( );
            break;
        }
        counter++;
    }
    if ( mat_number != this->nusselt_surface_id_ )
    {
        return;
    }
    const int num_q = this->num_quadrature_points ( );

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType r = this->x ( q )[1];
        const DataType wq = this->w ( q );
        const DataType dJ = r * std::abs ( this->ds ( q ) );

        lv[facet_number] += wq * this->grad_solP_[t_var][q][1] * dJ;
    }
}

template class MetFlowBousCylAssembler<2, double>;
template class MetFlowBousCylAssembler<3, double>;
