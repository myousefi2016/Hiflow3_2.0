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

#include "met_flow_incomp_cyl_dual_assembler.h"

template<int DIM, class DataType>
MetFlowIncompCylDualAssembler<DIM, DataType>::MetFlowIncompCylDualAssembler ( )
: MetFlowIncompAssembler<DIM, DataType>( )
{
    this->mode_ = DUAL;
}

/// ********************************************************
/// General Assembly routines 
/// ********************************************************

template<int DIM, class DataType>
void MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_matrix ( const Element<double>& element, LocalMatrix& lm ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( element, lm );
}

template<int DIM, class DataType>
void MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_vector ( const Element<double>& element, LocalVector& lv ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_vector_dual ( element, lv );
}

template<int DIM, class DataType>
void MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_scalar ( const Element<double>& element, double& ls ) const
{
    ls = 0.;
}
/// ********************************************************
/// Assembly routines for dual problem
/// ********************************************************
/// Jacobian of dual problem
#ifdef OPT_ASM0

template<int DIM, class DataType>
void MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( const Element<double>& element, LocalMatrix& lm ) const
{
    const int p_var = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
#    endif

    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    double dt = this->dT_cn_;
    if ( this->galerkin_mode_ == 1 )
        dt = this->dT_pc_;

    double sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    double tau_div = 0.;
    double skew_fac = 1.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCylDualAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const double r = this->x ( q )[1];
        const double inv_r = 1 / r;
        const double inv_rr = inv_r * inv_r;
        const double wq = this->w ( q );
        const double dJ = r * std::fabs ( this->detJ ( q ) );
        const double inv_r_comp[3] = { inv_r, 1., 1. };
        const double inv_rr_comp[3] = { inv_rr, 1., 1. };

        // ***********************************************************************
        // TIME DIFFERENCE 
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq
                        * this->delta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                        * this->phi ( j, q, u_var )
                    * this->phi ( i, q, u_var )
                    * dJ;

            }
        }

        // ***********************************************************************
        // LAPLACE: delta_mom_vis_p * dT * nu * int{grad{u}:grad{v}}
        // symmetric terms
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    double tmp = 0;
                    for ( int s = 0; s < DIM; ++s )
                    { // grad(u) : grad(v)
                        tmp += inv_rr_comp[s]
                                * this->grad_phi ( j, q, u_var )[s]
                                * this->grad_phi ( i, q, u_var )[s];
                    }
                    if ( u_var < DIM - 1 )
                    { // (u_phi * v_phi + u_r * v_r) / r²
                        tmp += inv_rr
                                * this->phi ( j, q, u_var )
                                * this->phi ( i, q, u_var );
                    }
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq
                            * this->delta_mom_vis_c_
                            * dt
                            * this->nu_
                            * tmp
                            * dJ;
                }
            }
        }

        // (mixed terms: involving d_phi, u_r, u_phi, v_r, v_phi    
        for ( int v_var = 0; v_var < DIM - 1; ++v_var )
        {
            for ( int u_var = 0; u_var < DIM - 1; ++u_var )
            {
                if ( v_var != u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            if ( u_var == 0 ) sign = 1.0;
                            else sign = -1.0;

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq
                                    * this->delta_mom_vis_c_
                                    * dt
                                    * this->nu_
                                    * sign
                                    * inv_rr
                                    * ( this->phi ( i, q, v_var ) * this->grad_phi ( j, q, u_var )[0]
                                    - this->grad_phi ( i, q, v_var )[0] * this->phi ( j, q, u_var ) )
                                    * dJ;
                        }
                    }
                }
            }
        }
        // ***********************************************************************
        // CONVECTION TERM
        // ------------------------------------
        // PART 1: (u_P * grad(u_D), v)
        // -(0.5 *) delta_mom_adv_cp * dT * int{ (u_c*\grad{u})*v } - (0.5 *) delta_mom_adv_pp * dT * int{ (u_p*\grad{u})*v }
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
            {
                double tmp = 0;
                for ( int s = 0; s < DIM; ++s ) // scalar product
                    tmp += inv_r_comp[s] * ( this->delta_mom_adv_nc_ * this->solP_next_[s][q] + this->delta_mom_adv_cc_ * this->solP_[s][q] + this->delta_mom_adv_pc_ * this->solP_prev_[s][q] ) * this->grad_phi ( j, q, u_var )[s];

                for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                {
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq
                            * skew_fac
                            * dt
                            * tmp
                            * this->phi ( i, q, u_var )
                            * dJ;
                }
            }
        }
        for ( int i = 0; i < this->num_dofs ( 0 ); ++i )
        {
            for ( int j = 0; j < this->num_dofs ( 1 ); ++j )
            {
                lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 1 ) ) += -wq
                        * skew_fac
                        * dt
                        * ( this->delta_mom_adv_nc_ * ( inv_r * this->solP_next_[0][q] * this->phi ( j, q, 1 ) )
                        + this->delta_mom_adv_cc_ * ( inv_r * this->solP_[0][q] * this->phi ( j, q, 1 ) )
                        + this->delta_mom_adv_pc_ * ( inv_r * this->solP_prev_[0][q] * this->phi ( j, q, 1 ) )
                        )
                        * this->phi ( i, q, 0 )
                        * dJ;
            }
        }
        for ( int i = 0; i < this->num_dofs ( 1 ); ++i )
        {
            for ( int j = 0; j < this->num_dofs ( 0 ); ++j )
            {
                lm ( this->dof_index ( i, 1 ), this->dof_index ( j, 0 ) ) += wq
                        * skew_fac
                        * dt
                        * ( this->delta_mom_adv_nc_ * ( inv_r * this->solP_next_[0][q] * this->phi ( j, q, 0 ) )
                        + this->delta_mom_adv_cc_ * ( inv_r * this->solP_[0][q] * this->phi ( j, q, 0 ) )
                        + this->delta_mom_adv_pc_ * ( inv_r * this->solP_prev_[0][q] * this->phi ( j, q, 0 ) ) )
                        * this->phi ( i, q, 1 )
                        * dJ;
            }
        }

        // ------------------------------------
        // Part 2: (u_D * grad^T (u_P), v)
        // nonlinear l2k part three and four : (0.5 *) delta_mom_adv_cp_ * dT * int{ (u\grad^T {u_c}*v } + (0.5 *) delta_mom_adv_pp_ * dT * int{ (u \grad^{T} {u_p}*v }
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                    {
                        double tmp = 0.0;
                        tmp = inv_r_comp[v_var] * this->phi ( j, q, u_var ) * this->phi ( i, q, v_var )
                                * ( this->delta_mom_adv_nc_ * this->grad_solP_next_[u_var][q][v_var] + this->delta_mom_adv_cc_ * this->grad_solP_[u_var][q][v_var] + this->delta_mom_adv_pc_ * this->grad_solP_prev_[u_var][q][v_var] );

                        if ( v_var == 0 && u_var == 0 ) tmp += inv_r * this->phi ( j, q, 0 )
                            * ( this->delta_mom_adv_nc_ * this->solP_next_[1][q]
                                + this->delta_mom_adv_cc_ * this->solP_[1][q]
                                + this->delta_mom_adv_pc_ * this->solP_prev_[1][q] )
                            * this->phi ( i, q, 0 );
                        if ( v_var == 0 && u_var == 1 ) tmp -= inv_r * this->phi ( j, q, 1 )
                            * ( this->delta_mom_adv_nc_ * this->solP_next_[0][q]
                                + this->delta_mom_adv_cc_ * this->solP_[0][q]
                                + this->delta_mom_adv_pc_ * this->solP_prev_[0][q] )
                            * this->phi ( i, q, 0 );

                        lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq
                                * skew_fac
                                * dt
                                * tmp
                                * dJ;
                    }
                }
            }
        }

        // ************************************************************************
        // REACTION TERM: (div(u_P) u_D, v)
        // (0.5 *) delta_mom_rea_3_ * dT * int{ (div(u_c) u*v } + (0.5 *) delta_mom_rea_4_ * dT * int{ (div(u_p) u *v }
        double div_n = inv_r * this->solP_next_[1][q];
        double div_c = inv_r * this->solP_[1][q];
        double div_p = inv_r * this->solP_prev_[1][q];
        for ( int s = 0; s < DIM; s++ )
        {
            div_n += inv_r_comp[s] * this->grad_solP_next_[s][q][s];
            div_c += inv_r_comp[s] * this->grad_solP_[s][q][s];
            div_p += inv_r_comp[s] * this->grad_solP_prev_[s][q][s];
        }

        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
            {
                for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq
                        * skew_fac
                        * dt
                        * ( div_n * this->delta_mom_rea_nc_ + div_c * this->delta_mom_rea_cc_ + div_p * this->delta_mom_rea_pc_ )
                    * this->phi ( i, q, u_var )
                    * this->phi ( j, q, u_var )
                    * dJ;
            }
        }

#    ifdef ROTATING_FOR
        // CORIOLIS FORCE:  dT * delta_mom_rot_p * 2 * omega (u[0] * v[1] - u[1] * v[0])
        int u_var = 0;
        int v_var = 1;
        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq
                        * this->delta_mom_rot_c_
                        * dt
                        * 2.0
                        * this->omega_
                        * this->phi ( j, q, u_var )
                        * this->phi ( i, q, v_var )
                        * dJ;
            }
        }

        u_var = 1;
        v_var = 0;
        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += -wq
                        * this->delta_mom_rot_c_
                        * dt
                        * 2.0
                        * this->omega_
                        * this->phi ( j, q, u_var )
                        * this->phi ( i, q, v_var )
                        * dJ;
            }
        }
#    endif  

        // ***********************************************************************
        // PRESSURE: + dT * delta_mom_pre * \int{p, div{v}}
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                double div_v = 0.0;

                switch ( v_var )
                {
                    case 1:
                        div_v = inv_r * this->grad_phi ( i, q, 0 )[0];
                        break;
                    case 2:
                        div_v = this->grad_phi ( i, q, 1 )[1] + inv_r * this->phi ( i, q, 1 );
                        break;
                    case 3:
                        div_v = this->grad_phi ( i, q, 2 )[2];
                        break;
                }

                const int ind_i = this->dof_index ( i, v_var );

                for ( int j = 0; j<this->num_dofs ( p_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p_var ) ) += wq
                            * dt
                            * this->delta_mom_pre_c_
                            * this->phi ( j, q, p_var )
                            * div_v
                            * dJ;
                }
#    ifdef AUGMENT_PRESS
                for ( int j = 0; j<this->num_dofs ( p0_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p0_var ) ) += wq
                            * dt
                            * this->delta_mom_pre_c_
                            * this->phi ( j, q, p0_var )
                            * div_v
                            * dJ;
                }
#    endif
            }
        }

        // ***********************************************************************
        // Grad-div stabilization: delta_mom_graddiv2 * dT * gamma_div * (div{u}, div{v})
        if ( this->graddiv_mode_ > 0 )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    double div_v = 0.0;
                    switch ( v_var )
                    {
                        case 1:
                            div_v = inv_r * this->grad_phi ( i, q, 0 )[0];
                            break;
                        case 2:
                            div_v = this->grad_phi ( i, q, 1 )[1] + inv_r * this->phi ( i, q, 1 );
                            break;
                        case 3:
                            div_v = this->grad_phi ( i, q, 2 )[2];
                            break;
                    }

                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            double div_u = 0.0;
                            if ( u_var == 0 ) div_u = inv_r * this->grad_phi ( j, q, 0 )[0];
                            if ( u_var == 1 ) div_u = inv_r * this->phi ( j, q, 1 ) + this->grad_phi ( j, q, 1 )[1];
                            if ( u_var == 2 ) div_u = this->grad_phi ( j, q, 2 )[2];

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq
                                    * tau_div
                                    * this->delta_mom_graddiv_c_
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
        // ***********************************************************************
        // CONSTRAINT: - inv_rho * dT * delta_inc_p * \int{q div(u)}
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
            {
                double div_u = 0.0;

                switch ( u_var )
                {
                    case 0:
                        div_u = inv_r * this->grad_phi ( j, q, 0 )[0];
                        break;
                    case 1:
                        div_u = this->grad_phi ( j, q, 1 )[1] + inv_r * this->phi ( j, q, 1 );
                        break;
                    case 2:
                        div_u = this->grad_phi ( j, q, 2 )[2];
                        break;
                }

                const int ind_j = this->dof_index ( j, u_var );

                for ( int i = 0; i<this->num_dofs ( p_var ); ++i )
                {
                    lm ( this->dof_index ( i, p_var ), ind_j ) += -wq
                            * this->inv_rho_
                            * dt
                            * this->delta_inc_c_
                            * this->phi ( i, q, p_var )
                            * div_u
                            * dJ;
                }
#    ifdef AUGMENT_PRESS
                for ( int i = 0; i<this->num_dofs ( p0_var ); ++i )
                {
                    lm ( this->dof_index ( i, p0_var ), ind_j ) += -wq
                            * this->inv_rho_
                            * dt
                            * this->delta_inc_c_
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
void MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( const Element<double>& element, LocalMatrix& lm ) const
{
    const int p_var = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
#    endif

    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    double dt = this->dT_cn_;
    if ( this->galerkin_mode_ == 1 )
        dt = this->dT_pc_;

    double sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    double tau_div = 0.;
    double skew_fac = 1.;
    double temp_supg = 0.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCylDualAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    const double mom_vis_c = this->delta_mom_vis_c_ * dt * this->nu_;
    const double mom_adv_pc = skew_fac * dt * this->delta_mom_adv_pc_;
    const double mom_adv_nc = skew_fac * dt * this->delta_mom_adv_nc_;
    const double mom_adv_cc = skew_fac * dt * this->delta_mom_adv_cc_;
    const double mom_rea_pc = skew_fac * dt * this->delta_mom_rea_pc_;
    const double mom_rea_nc = skew_fac * dt * this->delta_mom_rea_nc_;
    const double mom_rea_cc = skew_fac * dt * this->delta_mom_rea_cc_;
    const double mom_rot_c = this->delta_mom_rot_c_ * dt * 2.0 * this->omega_;
    const double mom_graddiv_c = tau_div * this->delta_mom_graddiv_c_ * dt;
    const double mom_pre_c = dt * this->delta_mom_pre_c_;
    const double inc_c = this->inv_rho_ * dt * this->delta_inc_c_;

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const double r = this->x ( q )[1];
        const double inv_r = 1 / r;
        const double inv_rr = inv_r * inv_r;
        const double wq_dJ = this->w ( q ) * r * std::fabs ( this->detJ ( q ) );
        const double inv_r_comp[3] = { inv_r, 1., 1. };
        const double inv_rr_comp[3] = { inv_rr, 1., 1. };

        // ***********************************************************************
        // TIME DIFFERENCE 
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                        * this->delta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                        * this->phi ( j, q, u_var )
                    * this->phi ( i, q, u_var );

            }
        }

        // ***********************************************************************
        // LAPLACE: delta_mom_vis_p * dT * nu * int{grad{u}:grad{v}}
        // symmetric terms
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    double tmp = 0;
                    for ( int s = 0; s < DIM; ++s )
                    { // grad(u) : grad(v)
                        tmp += inv_rr_comp[s]
                                * this->grad_phi ( j, q, u_var )[s]
                                * this->grad_phi ( i, q, u_var )[s];
                    }
                    if ( u_var < DIM - 1 )
                    { // (u_phi * v_phi + u_r * v_r) / r²
                        tmp += inv_rr
                                * this->phi ( j, q, u_var )
                                * this->phi ( i, q, u_var );
                    }
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * mom_vis_c
                            * tmp;
                }
            }
        }

        // (mixed terms: involving d_phi, u_r, u_phi, v_r, v_phi    
        for ( int v_var = 0; v_var < DIM - 1; ++v_var )
        {
            for ( int u_var = 0; u_var < DIM - 1; ++u_var )
            {
                if ( v_var != u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            if ( u_var == 0 ) sign = 1.0;
                            else sign = -1.0;

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                    * mom_vis_c
                                    * sign
                                    * inv_rr
                                    * ( this->phi ( i, q, v_var ) * this->grad_phi ( j, q, u_var )[0]
                                    - this->grad_phi ( i, q, v_var )[0] * this->phi ( j, q, u_var ) );
                        }
                    }
                }
            }
        }
        // ***********************************************************************
        // CONVECTION TERM
        // ------------------------------------
        // PART 1: (u_P * grad(u_D), v)
        // -(0.5 *) delta_mom_adv_cp * dT * int{ (u_c*\grad{u})*v } - (0.5 *) delta_mom_adv_pp * dT * int{ (u_p*\grad{u})*v }
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
            {
                double tmp = 0;
                for ( int s = 0; s < DIM; ++s ) // scalar product
                    tmp += inv_r_comp[s] * ( mom_adv_nc * this->solP_next_[s][q] + mom_adv_cc * this->solP_[s][q] + mom_adv_pc * this->solP_prev_[s][q] ) * this->grad_phi ( j, q, u_var )[s];

                for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                {
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                            * tmp
                            * this->phi ( i, q, u_var );
                }
            }
        }
        for ( int i = 0; i < this->num_dofs ( 0 ); ++i )
        {
            for ( int j = 0; j < this->num_dofs ( 1 ); ++j )
            {
                lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 1 ) ) += -wq_dJ
                        * ( mom_adv_nc * ( inv_r * this->solP_next_[0][q] * this->phi ( j, q, 1 ) )
                        + mom_adv_cc * ( inv_r * this->solP_[0][q] * this->phi ( j, q, 1 ) )
                        + mom_adv_pc * ( inv_r * this->solP_prev_[0][q] * this->phi ( j, q, 1 ) )
                        )
                        * this->phi ( i, q, 0 );
            }
        }
        for ( int i = 0; i < this->num_dofs ( 1 ); ++i )
        {
            for ( int j = 0; j < this->num_dofs ( 0 ); ++j )
            {
                lm ( this->dof_index ( i, 1 ), this->dof_index ( j, 0 ) ) += wq_dJ
                        * ( mom_adv_nc * ( inv_r * this->solP_next_[0][q] * this->phi ( j, q, 0 ) )
                        + mom_adv_cc * ( inv_r * this->solP_[0][q] * this->phi ( j, q, 0 ) )
                        + mom_adv_pc * ( inv_r * this->solP_prev_[0][q] * this->phi ( j, q, 0 ) )
                        )
                        * this->phi ( i, q, 1 );
            }
        }

        // ------------------------------------
        // Part 2: (u_D * grad^T (u_P), v)
        // nonlinear l2k part three and four : (0.5 *) delta_mom_adv_cp_ * dT * int{ (u\grad^T {u_c}*v } + (0.5 *) delta_mom_adv_pp_ * dT * int{ (u \grad^{T} {u_p}*v }
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                    {
                        double tmp = 0.0;
                        tmp = inv_r_comp[v_var] * this->phi ( j, q, u_var ) * this->phi ( i, q, v_var )
                                * ( mom_adv_nc * this->grad_solP_next_[u_var][q][v_var] + mom_adv_cc * this->grad_solP_[u_var][q][v_var] + mom_adv_pc * this->grad_solP_prev_[u_var][q][v_var] );

                        if ( v_var == 0 && u_var == 0 ) tmp += inv_r * this->phi ( j, q, 0 ) * ( mom_adv_nc * this->solP_next_[1][q] + mom_adv_cc * this->solP_[1][q] + mom_adv_pc * this->solP_prev_[1][q] ) * this->phi ( i, q, 0 );
                        if ( v_var == 0 && u_var == 1 ) tmp -= inv_r * this->phi ( j, q, 1 ) * ( mom_adv_nc * this->solP_next_[0][q] + mom_adv_cc * this->solP_[0][q] + mom_adv_pc * this->solP_prev_[0][q] ) * this->phi ( i, q, 0 );

                        lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ * tmp;
                    }
                }
            }
        }

        // ************************************************************************
        // REACTION TERM: (div(u_P) u_D, v)
        // (0.5 *) delta_mom_rea_3_ * dT * int{ (div(u_c) u*v } + (0.5 *) delta_mom_rea_4_ * dT * int{ (div(u_p) u *v }
        double div_n = inv_r * this->solP_next_[1][q];
        double div_c = inv_r * this->solP_[1][q];
        double div_p = inv_r * this->solP_prev_[1][q];
        for ( int s = 0; s < DIM; s++ )
        {
            div_n += inv_r_comp[s] * this->grad_solP_next_[s][q][s];
            div_c += inv_r_comp[s] * this->grad_solP_[s][q][s];
            div_p += inv_r_comp[s] * this->grad_solP_prev_[s][q][s];
        }

        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
            {
                for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                        * ( div_n * mom_rea_nc + div_c * mom_rea_cc + div_p * mom_rea_pc )
                    * this->phi ( i, q, u_var )
                    * this->phi ( j, q, u_var );
            }
        }

#    ifdef ROTATING_FOR
        // CORIOLIS FORCE:  dT * delta_mom_rot_p * 2 * omega (u[0] * v[1] - u[1] * v[0])
        int u_var = 0;
        int v_var = 1;
        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                        * mom_rot_c
                        * this->phi ( j, q, u_var )
                        * this->phi ( i, q, v_var );
            }
        }

        u_var = 1;
        v_var = 0;
        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                        * mom_rot_c
                        * this->phi ( j, q, u_var )
                        * this->phi ( i, q, v_var );
            }
        }
#    endif  

        // ***********************************************************************
        // PRESSURE: + dT * delta_mom_pre * \int{p, div{v}}
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                double div_v = 0.0;

                switch ( v_var )
                {
                    case 1:
                        div_v = inv_r * this->grad_phi ( i, q, 0 )[0];
                        break;
                    case 2:
                        div_v = this->grad_phi ( i, q, 1 )[1] + inv_r * this->phi ( i, q, 1 );
                        break;
                    case 3:
                        div_v = this->grad_phi ( i, q, 2 )[2];
                        break;
                }
                const int ind_i = this->dof_index ( i, v_var );

                for ( int j = 0; j<this->num_dofs ( p_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p_var ) ) += wq_dJ
                            * dt
                            * this->delta_mom_pre_c_
                            * this->phi ( j, q, p_var )
                            * div_v;
                }
#    ifdef AUGMENT_PRESS
                for ( int j = 0; j<this->num_dofs ( p0_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p0_var ) ) += wq_dJ
                            * dt
                            * this->delta_mom_pre_c_
                            * this->phi ( j, q, p0_var )
                            * div_v;
                }
#    endif
            }
        }

        // ***********************************************************************
        // Grad-div stabilization: delta_mom_graddiv2 * dT * gamma_div * (div{u}, div{v})
        if ( this->graddiv_mode_ > 0 )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    double div_v = 0.0;
                    switch ( v_var )
                    {
                        case 1:
                            div_v = inv_r * this->grad_phi ( i, q, 0 )[0];
                            break;
                        case 2:
                            div_v = this->grad_phi ( i, q, 1 )[1] + inv_r * this->phi ( i, q, 1 );
                            break;
                        case 3:
                            div_v = this->grad_phi ( i, q, 2 )[2];
                            break;
                    }

                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            double div_u = 0.0;
                            switch ( u_var )
                            {
                                case 1:
                                    div_u = inv_r * this->grad_phi ( j, q, 0 )[0];
                                    break;
                                case 2:
                                    div_u = this->grad_phi ( j, q, 1 )[1] + inv_r * this->phi ( j, q, 1 );
                                    break;
                                case 3:
                                    div_u = this->grad_phi ( j, q, 2 )[2];
                                    break;
                            }

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                    * mom_graddiv_c
                                    * div_u
                                    * div_v;
                        }
                    }
                }
            }
        }

        // ***********************************************************************       
        // ***********************************************************************
        // CONSTRAINT: - inv_rho * dT * delta_inc_p * \int{q div(u)}
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
            {
                double div_u = 0.0;
                switch ( u_var )
                {
                    case 1:
                        div_u = inv_r * this->grad_phi ( j, q, 0 )[0];
                        break;
                    case 2:
                        div_u = this->grad_phi ( j, q, 1 )[1] + inv_r * this->phi ( j, q, 1 );
                        break;
                    case 3:
                        div_u = this->grad_phi ( j, q, 2 )[2];
                        break;
                }

                const int ind_j = this->dof_index ( j, u_var );

                for ( int i = 0; i<this->num_dofs ( p_var ); ++i )
                {
                    lm ( this->dof_index ( i, p_var ), ind_j ) += -wq_dJ
                            * this->inv_rho_
                            * dt
                            * this->delta_inc_c_
                            * this->phi ( i, q, p_var )
                            * div_u;
                }
#    ifdef AUGMENT_PRESS
                for ( int i = 0; i<this->num_dofs ( p0_var ); ++i )
                {
                    lm ( this->dof_index ( i, p0_var ), ind_j ) += -wq_dJ
                            * this->inv_rho_
                            * dt
                            * this->delta_inc_c_
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
void MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( const Element<double>& element, LocalMatrix& lm ) const
{
    const int p_var = DIM;
    const int grad_vars = DIM;
#    ifdef AUGMENT_PRESS
    const int p0_var = DIM + 1;
    const int num_vars = DIM + 2;
#    else
    const int num_vars = DIM + 1;
#    endif

    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    double dt = this->dT_cn_;
    if ( this->galerkin_mode_ == 1 )
        dt = this->dT_pc_;

    double sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    double tau_div = 0.;
    double skew_fac = 1.;
    double temp_supg = 0.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCylDualAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    const double mom_vis_c = this->delta_mom_vis_c_ * dt * this->nu_;
    const double mom_adv_pc = skew_fac * dt * this->delta_mom_adv_pc_;
    const double mom_adv_nc = skew_fac * dt * this->delta_mom_adv_nc_;
    const double mom_adv_cc = skew_fac * dt * this->delta_mom_adv_cc_;
    const double mom_rea_pc = skew_fac * dt * this->delta_mom_rea_pc_;
    const double mom_rea_nc = skew_fac * dt * this->delta_mom_rea_nc_;
    const double mom_rea_cc = skew_fac * dt * this->delta_mom_rea_cc_;
    const double mom_rot_c = this->delta_mom_rot_c_ * dt * 2.0 * this->omega_;
    const double mom_graddiv_c = tau_div * this->delta_mom_graddiv_c_ * dt;
    const double mom_pre_c = dt * this->delta_mom_pre_c_;
    const double inc_c = this->inv_rho_ * dt * this->delta_inc_c_;

    std::vector< std::vector<double> > phi ( num_vars );
    std::vector< std::vector< std::vector<double> > > grad_phi ( grad_vars );

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
        const double r = this->x ( q )[1];
        const double inv_r = 1 / r;
        const double inv_rr = inv_r * inv_r;
        const double wq_dJ = this->w ( q ) * r * std::fabs ( this->detJ ( q ) );
        const double inv_r_comp[3] = { inv_r, 1., 1. };
        const double inv_rr_comp[3] = { inv_rr, 1., 1. };

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

        // ***********************************************************************
        // TIME DIFFERENCE 
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i<this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                        * this->delta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                        * phi[u_var][j]
                        * phi[u_var][i];

            }
        }

        // ***********************************************************************
        // LAPLACE: delta_mom_vis_p * dT * nu * int{grad{u}:grad{v}}
        // symmetric terms
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
            {
                for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                {
                    double tmp = 0;
                    for ( int s = 0; s < DIM; ++s )
                    { // grad(u) : grad(v)
                        tmp += inv_rr_comp[s]
                                * grad_phi[u_var][j][s]
                                * grad_phi[u_var][i][s];
                    }
                    if ( u_var < DIM - 1 )
                    { // (u_phi * v_phi + u_r * v_r) / r²
                        tmp += inv_rr
                                * phi[u_var][j]
                                * phi[u_var][i];
                    }
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * mom_vis_c
                            * tmp;
                }
            }
        }

        // (mixed terms: involving d_phi, u_r, u_phi, v_r, v_phi    
        for ( int v_var = 0; v_var < DIM - 1; ++v_var )
        {
            for ( int u_var = 0; u_var < DIM - 1; ++u_var )
            {
                if ( v_var != u_var )
                {
                    for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            if ( u_var == 0 ) sign = 1.0;
                            else sign = -1.0;

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                    * mom_vis_c
                                    * sign
                                    * inv_rr
                                    * ( phi[v_var][i] * grad_phi[u_var][j][0]
                                    - grad_phi[v_var][i][0] * phi[u_var][j] );
                        }
                    }
                }
            }
        }
        // ***********************************************************************
        // CONVECTION TERM
        // ------------------------------------
        // PART 1: (u_P * grad(u_D), v)

        // -(0.5 *) delta_mom_adv_cp * dT * int{ (u_c*\grad{u})*v } - (0.5 *) delta_mom_adv_pp * dT * int{ (u_p*\grad{u})*v }
        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
            {
                double tmp = 0;
                for ( int s = 0; s < DIM; ++s ) // scalar product
                    tmp += inv_r_comp[s] * ( mom_adv_nc * this->solP_next_[s][q] + mom_adv_cc * this->solP_[s][q] + mom_adv_pc * this->solP_prev_[s][q] ) * grad_phi[u_var][j][s];

                for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                {
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                            * tmp
                            * phi[u_var][i];
                }
            }
        }
        for ( int i = 0; i < this->num_dofs ( 0 ); ++i )
        {
            for ( int j = 0; j < this->num_dofs ( 1 ); ++j )
            {
                lm ( this->dof_index ( i, 0 ), this->dof_index ( j, 1 ) ) += -wq_dJ
                        * ( mom_adv_nc * ( inv_r * this->solP_next_[0][q] * phi[1][j] )
                        + mom_adv_cc * ( inv_r * this->solP_[0][q] * phi[1][j] )
                        + mom_adv_pc * ( inv_r * this->solP_prev_[0][q] * phi[1][j] )
                        )
                        * phi[0][i];
            }
        }
        for ( int i = 0; i < this->num_dofs ( 1 ); ++i )
        {
            for ( int j = 0; j < this->num_dofs ( 0 ); ++j )
            {
                lm ( this->dof_index ( i, 1 ), this->dof_index ( j, 0 ) ) += wq_dJ
                        * ( mom_adv_nc * ( inv_r * this->solP_next_[0][q] * phi[0][j] )
                        + mom_adv_cc * ( inv_r * this->solP_[0][q] * phi[0][j] )
                        + mom_adv_pc * ( inv_r * this->solP_prev_[0][q] * phi[0][j] )
                        )
                        * phi[1][i];
            }
        }

        // ------------------------------------
        // Part 2: (u_D * grad^T (u_P), v)
        // nonlinear l2k part three and four : (0.5 *) delta_mom_adv_cp_ * dT * int{ (u\grad^T {u_c}*v } + (0.5 *) delta_mom_adv_pp_ * dT * int{ (u \grad^{T} {u_p}*v }
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int u_var = 0; u_var < DIM; ++u_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                    {
                        double tmp = 0.0;
                        tmp = inv_r_comp[v_var] * phi[u_var][j] * phi[v_var][i]
                                * ( mom_adv_nc * this->grad_solP_next_[u_var][q][v_var] + mom_adv_cc * this->grad_solP_[u_var][q][v_var] + mom_adv_pc * this->grad_solP_prev_[u_var][q][v_var] );

                        if ( v_var == 0 && u_var == 0 ) tmp += inv_r * phi[0][j] * ( mom_adv_nc * this->solP_next_[1][q] + mom_adv_cc * this->solP_[1][q] + mom_adv_pc * this->solP_prev_[1][q] ) * phi[0][i];
                        if ( v_var == 0 && u_var == 1 ) tmp -= inv_r * phi[1][j] * ( mom_adv_nc * this->solP_next_[0][q] + mom_adv_cc * this->solP_[0][q] + mom_adv_pc * this->solP_prev_[0][q] ) * phi[0][i];

                        lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                * tmp;
                    }
                }
            }
        }

        // ************************************************************************
        // REACTION TERM: (div(u_P) u_D, v)
        // (0.5 *) delta_mom_rea_3_ * dT * int{ (div(u_c) u*v } + (0.5 *) delta_mom_rea_4_ * dT * int{ (div(u_p) u *v }
        double div_n = inv_r * this->solP_next_[1][q];
        double div_c = inv_r * this->solP_[1][q];
        double div_p = inv_r * this->solP_prev_[1][q];
        for ( int s = 0; s < DIM; s++ )
        {
            div_n += inv_r_comp[s] * this->grad_solP_next_[s][q][s];
            div_c += inv_r_comp[s] * this->grad_solP_[s][q][s];
            div_p += inv_r_comp[s] * this->grad_solP_prev_[s][q][s];
        }

        for ( int u_var = 0; u_var < DIM; ++u_var )
        {
            for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
            {
                for ( int i = 0; i < this->num_dofs ( u_var ); ++i )
                    lm ( this->dof_index ( i, u_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                        * ( div_n * mom_rea_nc + div_c * mom_rea_cc + div_p * mom_rea_pc )
                    * phi[u_var][i]
                        * phi[u_var][j];
            }
        }

#    ifdef ROTATING_FOR
        // CORIOLIS FORCE:  dT * delta_mom_rot_p * 2 * omega (u[0] * v[1] - u[1] * v[0])
        int u_var = 0;
        int v_var = 1;
        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                        * mom_rot_c
                        * phi[u_var][j]
                        * phi[v_var][i];
            }
        }

        u_var = 1;
        v_var = 0;
        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += -wq_dJ
                        * mom_rot_c
                        * phi[u_var][j]
                        * phi[v_var][i];
            }
        }
#    endif  

        // ***********************************************************************
        // PRESSURE: + dT * delta_mom_pre * \int{p, div{v}}
        for ( int v_var = 0; v_var < DIM; v_var++ )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                double div_v = 0.0;

                if ( v_var == 0 ) div_v = inv_r * grad_phi[0][i][0];
                if ( v_var == 1 ) div_v = grad_phi[1][i][1] + inv_r * phi[1][i];
                if ( v_var == 2 ) div_v = grad_phi[2][i][2];

                int ind_i = this->dof_index ( i, v_var );

                for ( int j = 0; j<this->num_dofs ( p_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p_var ) ) += wq_dJ
                            * dt
                            * this->delta_mom_pre_c_
                            * phi[p_var][j]
                            * div_v;
                }
#    ifdef AUGMENT_PRESS
                for ( int j = 0; j<this->num_dofs ( p0_var ); ++j )
                {
                    lm ( ind_i, this->dof_index ( j, p0_var ) ) += wq_dJ
                            * dt
                            * this->delta_mom_pre_c_
                            * phi[p0_var][j]
                            * div_v;
                }
#    endif            
            }
        }

        // ***********************************************************************
        // Grad-div stabilization: delta_mom_graddiv2 * dT * gamma_div * (div{u}, div{v})
        if ( this->graddiv_mode_ > 0 )
        {
            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    double div_v = 0.0;
                    if ( v_var == 0 ) div_v = inv_r * grad_phi[0][i][0];
                    if ( v_var == 1 ) div_v = inv_r * phi[1][i] + grad_phi[1][i][1];
                    if ( v_var == 2 ) div_v = grad_phi[2][i][2];

                    for ( int u_var = 0; u_var < DIM; ++u_var )
                    {
                        for ( int j = 0; j < this->num_dofs ( u_var ); ++j )
                        {
                            double div_u = 0.0;
                            if ( u_var == 0 ) div_u = inv_r * grad_phi[0][j][0];
                            if ( u_var == 1 ) div_u = inv_r * phi[1][j] + grad_phi[1][j][1];
                            if ( u_var == 2 ) div_u = grad_phi[2][j][2];

                            lm ( this->dof_index ( i, v_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                                    * mom_graddiv_c
                                    * div_u
                                    * div_v;
                        }
                    }
                }
            }
        }

        // ***********************************************************************       
        // ***********************************************************************
        // CONSTRAINT: - inv_rho * dT * delta_inc_p * \int{q div(u)}
        for ( int u_var = 0; u_var < DIM; u_var++ )
        {
            for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
            {
                double div_u = 0.0;
                switch ( u_var )
                {
                    case 0:
                        div_u = inv_r * grad_phi[0][j][0];
                        break;
                    case 1:
                        div_u = grad_phi[1][j][1] + inv_r * phi[1][j];
                        break;
                    case 2:
                        div_u = grad_phi[2][j][2];
                        break;
                }

                const int ind_j = this->dof_index ( j, u_var );

                for ( int i = 0; i<this->num_dofs ( p_var ); ++i )
                {
                    lm ( this->dof_index ( i, p_var ), ind_j ) += -wq_dJ
                            * this->inv_rho_
                            * dt
                            * this->delta_inc_c_
                            * phi[p_var][i]
                            * div_u;
                }
#    ifdef AUGMENT_PRESS
                for ( int i = 0; i<this->num_dofs ( p0_var ); ++i )
                {
                    lm ( this->dof_index ( i, p0_var ), ind_j ) += -wq_dJ
                            * this->inv_rho_
                            * dt
                            * this->delta_inc_c_
                            * phi[p0_var][i]
                            * div_u;
                }
#    endif
            }
        }
    }
}
#endif

/// Residual of dual problem (rhs - lhs)

template<int DIM, class DataType>
void MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_vector_dual ( const Element<double>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );

    const int p_var = DIM;
    const int grad_vars = DIM;

#ifdef AUGMENT_PRESS
    const int num_vars = DIM + 2;
    const int p0_var = DIM + 1;
#else
    const int num_vars = DIM + 1;
    const int p0_var = -1;
#endif

    double sign = 1.0;

    double div_n, div_c;

    double tau_div = 0.;
    double skew_fac = 1.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowIncompCylDualAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    double dt_n;
    double dt_c;
    double dt_p;

    if ( this->galerkin_mode_ == 0 && this->delta_d_dt_u_ == 1. )
    {
        dt_n = this->dT_cn_;
        dt_c = this->dT_cn_;
        dt_p = 0.;
    }
    else
    {
        dt_n = this->dT_cn_;
        dt_c = this->dT_pc_;
        dt_p = this->dT_pc_;
    }

    const double mom_vis_n = this->delta_mom_vis_n_ * dt_n * this->nu_;
    const double mom_vis_c = this->delta_mom_vis_c_ * dt_c * this->nu_;
    const double mom_adv_nn = skew_fac * dt_n * this->delta_mom_adv_nn_;
    const double mom_adv_cn = skew_fac * dt_n * this->delta_mom_adv_cn_;
    const double mom_adv_nc = skew_fac * dt_n * this->delta_mom_adv_nc_;
    const double mom_adv_cc = skew_fac * dt_c * this->delta_mom_adv_cc_;
    const double mom_adv_pc = skew_fac * dt_p * this->delta_mom_adv_pc_;

    const double mom_rea_nn = skew_fac * dt_n * this->delta_mom_rea_nn_;
    const double mom_rea_cn = skew_fac * dt_n * this->delta_mom_rea_cn_;
    const double mom_rea_nc = skew_fac * dt_n * this->delta_mom_rea_nc_;
    const double mom_rea_cc = skew_fac * dt_c * this->delta_mom_rea_cc_;
    const double mom_rea_pc = skew_fac * dt_p * this->delta_mom_rea_pc_;

    const double mom_rot_n = this->delta_mom_rot_n_ * dt_n * 2.0 * this->omega_;
    const double mom_rot_c = this->delta_mom_rot_c_ * dt_c * 2.0 * this->omega_;
    const double mom_graddiv_n = tau_div * this->delta_mom_graddiv_n_ * dt_n;
    const double mom_graddiv_c = tau_div * this->delta_mom_graddiv_c_ * dt_c;

    const double mom_pre_c = dt_c * this->delta_mom_pre_c_;
    const double mom_pre_n = dt_n * this->delta_mom_pre_n_;
    const double inc_n = this->inv_rho_ * dt_n * this->delta_inc_n_;
    const double inc_c = this->inv_rho_ * dt_c * this->delta_inc_c_;

    const double j_c_c = dt_p * this->delta_j_c_c_;
    const double j_c_pc = dt_p * this->delta_j_c_pc_;
    const double j_n_c = dt_n * this->delta_j_n_c_;
    const double j_n_nc = dt_n * this->delta_j_n_nc_;
    const double j_n_n = dt_n * this->delta_j_n_n_;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const double r = this->x ( q )[1];
        const double inv_r = 1. / r;
        const double inv_rr = inv_r * inv_r;
        const double wq_dJ = this->w ( q ) * r * std::abs ( this->detJ ( q ) );
        const double inv_r_comp[3] = { inv_r, 1., 1. };
        const double inv_rr_comp[3] = { inv_r*inv_r, 1., 1. };

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
        std::vector<double> solP_c ( num_vars );
        std::vector<double> solP_n ( num_vars );
        std::vector<double> solP_p ( num_vars );
        std::vector<double> solD_c ( num_vars );
        std::vector<double> solD_n ( num_vars );

        for ( int var = 0; var < DIM; ++var )
        {
            solP_n[var] = this->solP_next_[var][q];
            solD_n[var] = this->solD_next_[var][q];
            solP_c[var] = this->solP_[var][q];
            solD_c[var] = this->solD_[var][q];
            solP_p[var] = this->solP_prev_[var][q];
        }
        solP_n[p_var] = this->solP_next_[p_var][q];
        solD_n[p_var] = this->solD_next_[p_var][q];

        solP_c[p_var] = this->solP_[p_var][q];
        solD_c[p_var] = this->solD_[p_var][q];

        solP_p[p_var] = this->solP_prev_[p_var][q];

#ifdef AUGMENT_PRESS
        solP_n[p0_var] = this->solP_next_[p0_var][q];
        solP_p[p0_var] = this->solP_prev_[p0_var][q];
        solP_c[p0_var] = this->solP_[p0_var][q];
        solD_n[p0_var] = this->solD_next_[p0_var][q];
        solD_c[p0_var] = this->solD_[p0_var][q];
#endif    

        std::vector< Vec<DIM, DataType> > grad_solP_n ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_solP_c ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_solP_p ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_solD_n ( grad_vars );
        std::vector< Vec<DIM, DataType> > grad_solD_c ( grad_vars );

        for ( int var = 0; var < grad_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                grad_solP_n[var][d] = this->grad_solP_next_[var][q][d];
                grad_solP_c[var][d] = this->grad_solP_[var][q][d];
                grad_solP_p[var][d] = this->grad_solP_prev_[var][q][d];
                grad_solD_n[var][d] = this->grad_solD_next_[var][q][d];
                grad_solD_c[var][d] = this->grad_solD_[var][q][d];
            }
        }

        // setup parameter struct for goal functional
        ParametersForceType<DIM, DataType> gp;
        gp.x = this->x ( q );
        gp.absolute_time = this->t_;
        gp.solP = solP_c;
        gp.solP_next = solP_n;
        gp.solP_prev = solP_p;
        gp.grad_solP = grad_solP_c;
        gp.grad_solP_next = grad_solP_n;
        gp.grad_solP_prev = grad_solP_p;

        // ********************************************************************** 
        // time-diff: l0(v) = -\int( dot(u_n - u_k, v))
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                lv[this->dof_index ( i, v_var )] += wq_dJ
                    * this->delta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                    * ( solD_c[v_var] - solD_n[v_var] )
                * this->phi ( i, q, v_var );
        }

        // ********************************************************************** 
        // LAPLACE  
        // explicit linear: l1n(v) = delta _mom_vis_p * dT * \nu * \int( \grad{uD_p} : \grad{v} )
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                double tmp = 0;
                for ( int s = 0; s < DIM; ++s )
                { // scalar product
                    tmp += inv_rr_comp[s]
                            * ( mom_vis_c * grad_solD_c[v_var][s] + mom_vis_n * grad_solD_n[v_var][s] )
                            * this->grad_phi ( i, q, v_var )[s];
                }
                if ( v_var == 0 ) tmp += -inv_rr * ( mom_vis_c * grad_solD_c[1][0] + mom_vis_n * grad_solD_n[1][0] ) * this->phi ( i, q, 0 )
                    + inv_rr * ( mom_vis_c * solD_c[1] + mom_vis_n * solD_n[1] ) * this->grad_phi ( i, q, 0 )[0]
                    + inv_rr * ( mom_vis_c * solD_c[0] + mom_vis_n * solD_n[0] ) * this->phi ( i, q, 0 );
                if ( v_var == 1 ) tmp += inv_rr * ( mom_vis_c * grad_solD_c[0][0] + mom_vis_n * grad_solD_n[0][0] ) * this->phi ( i, q, 1 )
                    - inv_rr * ( mom_vis_c * solD_c[0] + mom_vis_n * solD_n[0] ) * this->grad_phi ( i, q, 1 )[0]
                    + inv_rr * ( mom_vis_c * solD_c[1] + mom_vis_n * solD_n[1] ) * this->phi ( i, q, 1 );

                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * tmp;
            }
        }

        // ********************************************************************** 
        // CONVECTIVE TERM    I: (u_P * grad(u_D), v)
        // -(0.5 *) delta_mom_adv_cc * dT * int{ (uP_c*\grad{uD_c})*v } - (0.5 *) delta_mom_adv_pc * dT * int{ (uP_p*\grad{uD_c})*v }
        // -(0.5 *) delta_mom_adv_cp * dT * int{ (uP_c*\grad{uD_p})*v } - (0.5 *) delta_mom_adv_pp * dT * int{ (uP_p*\grad{uD_p})*v }

        double factor_nn[3] = { inv_r * solP_n[0] * solD_n[1], -inv_r * solP_n[0] * solD_n[0], 0. };
        double factor_cn[3] = { inv_r * solP_c[0] * solD_n[1], -inv_r * solP_c[0] * solD_n[0], 0. };
        double factor_nc[3] = { inv_r * solP_n[0] * solD_c[1], -inv_r * solP_n[0] * solD_c[0], 0. };
        double factor_cc[3] = { inv_r * solP_c[0] * solD_c[1], -inv_r * solP_c[0] * solD_c[0], 0. };
        double factor_pc[3] = { inv_r * solP_p[0] * solD_c[1], -inv_r * solP_p[0] * solD_c[0], 0. };

        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            double tmp_nn = 0.0;
            double tmp_cn = 0.0;
            double tmp_nc = 0.0;
            double tmp_cc = 0.0;
            double tmp_pc = 0.0;
            for ( int s = 0; s < DIM; s++ )
            {
                tmp_nn += inv_r_comp[s] * solP_n[s] * grad_solD_n[v_var][s];
                tmp_cn += inv_r_comp[s] * solP_c[s] * grad_solD_n[v_var][s];
                tmp_nc += inv_r_comp[s] * solP_n[s] * grad_solD_c[v_var][s];
                tmp_cc += inv_r_comp[s] * solP_c[s] * grad_solD_c[v_var][s];
                tmp_pc += inv_r_comp[s] * solP_p[s] * grad_solD_c[v_var][s];
            }

            // three additional guys of which only two remain
            tmp_nn += factor_nn[v_var];
            tmp_cn += factor_cn[v_var];
            tmp_nc += factor_nc[v_var];
            tmp_cc += factor_cc[v_var];
            tmp_pc += factor_pc[v_var];

            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += -wq_dJ
                        * ( mom_adv_nn * tmp_nn +
                        mom_adv_cn * tmp_cn +
                        mom_adv_nc * tmp_nc +
                        mom_adv_cc * tmp_cc +
                        mom_adv_pc * tmp_pc )
                        * this->phi ( i, q, v_var );
            }
        }

        // CONVECTIVE TERM    II: (u_P * grad^T (u_D), v)
        // (0.5 *) delta_mom_adv_cc * dT * int{ (uP_c*\grad^T{uD_c})*v } + (0.5 *) delta_mom_adv_pc * dT * int{ (uP_p*\grad^T{uD_c})*v }
        // (0.5 *) delta_mom_adv_cp * dT * int{ (uP_c*\grad^T{uD_p})*v } + (0.5 *) delta_mom_adv_pp * dT * int{ (uP_p*\grad^T{uD_p})*v }
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            double tmp_nn = 0.0;
            double tmp_nc = 0.0;
            double tmp_cn = 0.0;
            double tmp_cc = 0.0;
            double tmp_pc = 0.0;

            for ( int s = 0; s < DIM; s++ )
            {
                tmp_nn += inv_r_comp[v_var] * solD_n[s] * grad_solP_n[s][v_var];
                tmp_cn += inv_r_comp[v_var] * solD_n[s] * grad_solP_c[s][v_var];
                tmp_nc += inv_r_comp[v_var] * solD_c[s] * grad_solP_n[s][v_var];
                tmp_cc += inv_r_comp[v_var] * solD_c[s] * grad_solP_c[s][v_var];
                tmp_pc += inv_r_comp[v_var] * solD_c[s] * grad_solP_p[s][v_var];
            }

            if ( v_var == 0 )
            {
                tmp_nn += inv_r * solD_n[0] * solP_n[1] - inv_r * solD_n[1] * solP_n[0];
                tmp_cn += inv_r * solD_n[0] * solP_c[1] - inv_r * solD_n[1] * solP_c[0];
                tmp_nc += inv_r * solD_c[0] * solP_n[1] - inv_r * solD_c[1] * solP_n[0];
                tmp_cc += inv_r * solD_c[0] * solP_c[1] - inv_r * solD_c[1] * solP_c[0];
                tmp_pc += inv_r * solD_c[0] * solP_p[1] - inv_r * solD_c[1] * solP_p[0];
            }

            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * ( mom_adv_nn * tmp_nn +
                        mom_adv_cn * tmp_cn +
                        mom_adv_nc * tmp_nc +
                        mom_adv_cc * tmp_cc +
                        mom_adv_pc * tmp_pc )
                        * this->phi ( i, q, v_var );
            }
        }

        // **************************************************************
        // REACTIVE TERM: (div(uP) * uD, v)
        // -(0.5 *) delta_mom_rea_1 * dT * int{ div(uP_c) * uD_c * v } - (0.5 *) delta_mom_rea_2 * dT * int{ div(uP_p)* uD_c * v }
        // -(0.5 *) delta_mom_rea_3 * dT * int{ div(uP_c) * uD_p * v } - (0.5 *) delta_mom_rea_4 * dT * int{ div(uP_p)* uD_p * v }

        double divP_n = inv_r * solP_n[1];
        double divP_c = inv_r * solP_c[1];
        double divP_p = inv_r * solP_p[1];
        for ( int s = 0; s < DIM; s++ )
        {
            divP_n += inv_r_comp[s] * grad_solP_n[s][s];
            divP_c += inv_r_comp[s] * grad_solP_c[s][s];
            divP_p += inv_r_comp[s] * grad_solP_p[s][s];
        }
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += -wq_dJ
                        * ( mom_rea_nn * divP_n * solD_n[v_var] +
                        mom_rea_cn * divP_c * solD_n[v_var] +
                        mom_rea_nc * divP_n * solD_c[v_var] +
                        mom_rea_cc * divP_c * solD_c[v_var] +
                        mom_rea_pc * divP_p * solD_c[v_var] )
                        * this->phi ( i, q, v_var );
            }
        }

        // **************************************************************
        // PRESSURE TERM: (pD, div(v))
        // delta_mom_pre_ * dT * (pD, div(v))
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                double div_v = 0.0;

                if ( v_var == 0 ) div_v = inv_r * this->grad_phi ( i, q, 0 )[0];
                if ( v_var == 1 ) div_v = this->grad_phi ( i, q, 1 )[1] + inv_r * this->phi ( i, q, 1 );
                if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                lv[this->dof_index ( i, v_var )] += wq_dJ
#ifdef AUGMENT_PRESS
                        * ( mom_pre_c * solD_c[p_var] + mom_pre_n * solD_n[p_var] + mom_pre_c * solD_c[p0_var] + mom_pre_n * solD_n[p0_var] )
#else
                        * ( mom_pre_c * solD_c[p_var] + mom_pre_n * solD_n[p_var] )
#endif
                        * div_v;
            }
        }

#ifdef ROTATING_FOR
        // TODO überprüfen
        // Coriolis force:  dT * delta_mom_rot_c_ * 2 * omega (u_c[0] * v[1] - u_c[1] * v[0]) + dT * delta_mom_rot_p_ * 2 * omega (u_p[0] * v[1] - u_p[1] * v[0])
        int u_var = 0;
        int v_var = 1;
        for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
        {
            lv[this->dof_index ( i, v_var )] += wq_dJ
                    * ( mom_rot_n * solD_n[u_var] + mom_rot_c * solD_c[u_var] )
                    * this->phi ( i, q, v_var );
        }

        u_var = 1;
        v_var = 0;
        for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
        {
            lv[this->dof_index ( i, v_var )] += -wq_dJ
                    * ( mom_rot_n * solD_n[u_var] + mom_rot_c * solD_c[u_var] )
                    * this->phi ( i, q, v_var );
        }
#endif 

        double divD_n = inv_r * solD_n[1];
        for ( int s = 0; s < DIM; s++ )
            divD_n += inv_r_comp[s] * grad_solD_n[s][s];

        double divD_c = inv_r * solD_c[1];
        for ( int s = 0; s < DIM; s++ )
            divD_c += inv_r_comp[s] * grad_solD_c[s][s];

        // ***********************************************************************
        // GRADDIV: delta_mom_graddiv_c * dT * gamma_div * (div{solD_c}, div{v}) + delta_mom_graddiv_p_ * dT * gamma_div * (div{solD_p}, div{v})
#ifdef GRADDIV

        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                double div_v = 0.0;
                if ( v_var == 0 ) div_v = inv_r * this->grad_phi ( i, q, 0 )[0];
                if ( v_var == 1 ) div_v = inv_r * this->phi ( i, q, 1 ) + this->grad_phi ( i, q, 1 )[1];
                if ( v_var == 2 ) div_v = this->grad_phi ( i, q, 2 )[2];

                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * ( mom_graddiv_n * divD_n + mom_graddiv_c * divD_c )
                        * div_v;

            }
        }
#endif

        // **********************************************************************
        // GOAL FUNCTIONAL: delta_j_c * dT * j_u(solP, phi) + delta_j_c * dT * j_u(solP_prev, phi)
        if ( this->goal_functional_->GetIntegralType ( ) == 1 )
        {
            for ( int v_var = 0; v_var < DIM; v_var++ )
            {
                gp.var = v_var;
                for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
                {
                    gp.phi = this->phi ( i, q, v_var );
                    gp.grad_phi = this->grad_phi ( i, q, v_var );

                    double J_pc = 0.;
                    double J_nc = 0.;
                    double J_c = 0.;
                    double J_n = 0.;
                    if ( j_c_pc != 0. )
                    {
                        gp.relative_time = -0.5;
                        J_pc = this->goal_functional_->j_force_type ( gp );
                    }
                    if ( j_n_nc != 0. )
                    {
                        gp.relative_time = 0.5;
                        J_nc = this->goal_functional_->j_force_type ( gp );
                    }
                    if ( j_c_c != 0. || j_n_c != 0. )
                    {
                        gp.relative_time = 0.;
                        J_c = this->goal_functional_->j_force_type ( gp );
                    }
                    if ( j_n_n != 0. )
                    {
                        gp.relative_time = 1.;
                        J_n = this->goal_functional_->j_force_type ( gp );
                    }

                    lv[this->dof_index ( i, v_var )] += wq_dJ * ( j_c_c * J_c + j_c_pc * J_pc + j_n_c * J_c + j_n_nc * J_nc + j_n_n * J_n );
                }
            }
        }

        // add contribution of final condition
        if ( this->final_goal_contrib_ )
        {
            ParametersFinalType<DIM, DataType> gp;
            gp.solP.resize ( num_vars );
            gp.grad_solP.resize ( num_vars );

            for ( int var = 0; var < num_vars; ++var )
                gp.solP[var] = solP_c[var];

            for ( int var = 0; var < grad_vars; var++ )
            {
                for ( int d = 0; d < DIM; d++ )
                {
                    gp.grad_solP[var][d] = grad_solP_c[var][d];
                }
            }
            gp.x = this->x ( q );

            for ( int v_var = 0; v_var < DIM; ++v_var )
            {
                gp.var = v_var;
                for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
                {
                    gp.phi = this->phi ( i, q, v_var );
                    gp.grad_phi = this->grad_phi ( i, q, v_var );
                    lv[this->dof_index ( i, v_var )] += /*-*/wq_dJ * this->goal_functional_->j_final_type ( gp );
                }
            }
        }

        // **********************************************************************
        // ********************************************************************** 
        // CONSTRAINT: -delta_inc_c * eta * \int(q * div(uD_c)) - delta_inc_p * eta * \int(q * div(uD_p)) 

        for ( int i = 0; i < this->num_dofs ( p_var ); ++i )
        {
            lv[this->dof_index ( i, p_var )] += -wq_dJ
                    * ( inc_n * divD_n + inc_c * divD_c )
                    * this->phi ( i, q, p_var );
        }
#ifdef AUGMENT_PRESS
        for ( int i = 0; i < this->num_dofs ( p0_var ); ++i )
        {
            lv[this->dof_index ( i, p0_var )] += -wq_dJ
                    * ( inc_n * divD_n + inc_c * divD_c )
                    * this->phi ( i, q, p0_var );
        }
#endif

        // **********************************************************************
        // GOAL FUNCTIONAL: delta_j_c * dT * j_u(solP, phi) + delta_j_c * dT * j_u(solP_prev, phi)
        if ( this->goal_functional_->GetIntegralType ( ) == 1 )
        {
            gp.var = p_var;
            for ( int i = 0; i < this->num_dofs ( p_var ); ++i )
            {
                gp.phi = this->phi ( i, q, p_var );
                gp.grad_phi = this->grad_phi ( i, q, p_var );

                double J_pc = 0.;
                double J_nc = 0.;
                double J_c = 0.;
                double J_n = 0.;
                if ( j_c_pc != 0. )
                {
                    gp.relative_time = -0.5;
                    J_pc = this->goal_functional_->j_force_type ( gp );
                }
                if ( j_n_nc != 0. )
                {
                    gp.relative_time = 0.5;
                    J_nc = this->goal_functional_->j_force_type ( gp );
                }
                if ( j_c_c != 0. || j_n_c != 0. )
                {
                    gp.relative_time = 0.;
                    J_c = this->goal_functional_->j_force_type ( gp );
                }
                if ( j_n_n != 0. )
                {
                    gp.relative_time = 1.;
                    J_n = this->goal_functional_->j_force_type ( gp );
                }
                lv[this->dof_index ( i, p_var )] += wq_dJ * ( j_c_c * J_c + j_c_pc * J_pc + j_n_c * J_c + j_n_nc * J_nc + j_n_n * J_n );
            }
#ifdef AUGMENT_PRESS
            gp.var = p0_var;
            for ( int i = 0; i < this->num_dofs ( p0_var ); ++i )
            {
                gp.phi = this->phi ( i, q, p0_var );
                gp.grad_phi = this->grad_phi ( i, q, p0_var );

                double J_pc = 0.;
                double J_nc = 0.;
                double J_c = 0.;
                double J_n = 0.;
                if ( j_c_pc != 0. )
                {
                    gp.relative_time = -0.5;
                    J_pc = this->goal_functional_->j_force_type ( gp );
                }
                if ( j_n_nc != 0. )
                {
                    gp.relative_time = 0.5;
                    J_nc = this->goal_functional_->j_force_type ( gp );
                }
                if ( j_c_c != 0. || j_n_c != 0. )
                {
                    gp.relative_time = 0.;
                    J_c = this->goal_functional_->j_force_type ( gp );
                }
                if ( j_n_n != 0. )
                {
                    gp.relative_time = 1.;
                    J_n = this->goal_functional_->j_force_type ( gp );
                }
                lv[this->dof_index ( i, p0_var )] += wq_dJ * ( j_c_c * J_c + j_c_pc * J_pc + j_n_c * J_c + j_n_nc * J_nc + j_n_n * J_n );
            }
#endif
        }
    }
}

template class MetFlowIncompCylDualAssembler<2, double>;
template class MetFlowIncompCylDualAssembler<3, double>;

