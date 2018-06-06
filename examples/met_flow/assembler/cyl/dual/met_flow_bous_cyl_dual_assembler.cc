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

#include "met_flow_bous_cyl_dual_assembler.h"

template<int DIM, class DataType>
MetFlowBousCylDualAssembler<DIM, DataType>::MetFlowBousCylDualAssembler ( )
: MetFlowIncompCylDualAssembler<DIM, DataType>( ),
MetFlowBousAssembler<DIM, DataType>( )
{
    this->mode_ = DUAL;
}

/// ********************************************************
/// General Assembly routines 
/// ********************************************************

template<int DIM, class DataType>
void MetFlowBousCylDualAssembler<DIM, DataType>::assemble_local_matrix ( const Element<double>& element, LocalMatrix& lm ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( element, lm );
    MetFlowBousCylDualAssembler <DIM, DataType>::assemble_local_matrix_dual ( element, lm );
}

template<int DIM, class DataType>
void MetFlowBousCylDualAssembler<DIM, DataType>::assemble_local_vector ( const Element<double>& element, LocalVector& lv ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_vector_dual ( element, lv );
    MetFlowBousCylDualAssembler <DIM, DataType>::assemble_local_vector_dual ( element, lv );
}

template<int DIM, class DataType>
void MetFlowBousCylDualAssembler<DIM, DataType>::assemble_local_scalar ( const Element<double>& element, double& ls ) const
{
    ls = 0.;
    MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_scalar ( element, ls );
}

/// ********************************************************
/// Assembly routines for dual problem
/// ********************************************************
/// Jacobian of dual problem
#ifdef OPT_ASM0

template<int DIM, class DataType>
void MetFlowBousCylDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( const Element<double>& element, LocalMatrix& lm ) const
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

    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    double dt = this->dT_cn_;
    if ( this->galerkin_mode_ == 1 )
        dt = this->dT_pc_;

    std::vector<double> grav;
    grav.resize ( 3, 0.0 );
    double sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    double skew_fac = 1.;
    double temp_supg = 0.;

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

#    ifdef ROTATING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;
#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif  

        double div_n = inv_r * this->solP_next_[1][q];
        double div_c = inv_r * this->solP_[1][q];
        double div_p = inv_r * this->solP_prev_[1][q];
        for ( int s = 0; s < DIM; s++ )
        {
            div_n += inv_r_comp[s] * this->grad_solP_next_[s][q][s];
            div_c += inv_r_comp[s] * this->grad_solP_[s][q][s];
            div_p += inv_r_comp[s] * this->grad_solP_prev_[s][q][s];
        }

        // ***********************************************************************
        // SOURCE TERM: -(-grad(T_P) T_D, v) 
        //  delta_mom_sou_3_ * dT * (grad(T_c) * T,v) + delta_mom_sou_4_ * dT * (grad(T_p) * T,v)  
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                double tmp = ( this->delta_mom_sou_nc_ * this->grad_solP_next_[t_var][q][v_var] + this->delta_mom_sou_cc_ * this->grad_solP_[t_var][q][v_var] + this->delta_mom_sou_pc_ * this->grad_solP_prev_[t_var][q][v_var] ) * inv_r_comp[v_var] * this->phi ( i, q, v_var );

                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, v_var ), this->dof_index ( j, t_var ) ) += wq
                            * dt
                            * tmp
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
                // TIME-DERIVATIVE: -theta_d_dt_u * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * this->delta_d_dt_u_
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var )
                        * dJ;

                // ***********************************************************************
                // Thermal diffusion:  delta_heat_diff_p_ * kappa * dT * int{grad{u}:grad{v})
                double tmp = inv_rr * this->grad_phi ( i, q, t_var )[0] * this->grad_phi ( j, q, t_var )[0]
                        + this->grad_phi ( i, q, t_var )[1] * this->grad_phi ( j, q, t_var )[1];
                if ( DIM == 3 ) tmp += this->grad_phi ( i, q, t_var )[2] * this->grad_phi ( j, q, t_var )[2];

                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * this->kappa_
                        * dt
                        * this->delta_heat_diff_c_
                        * tmp
                        * dJ;
            }
        }

        // *****************************************************************************
        // Nonlinear advection : -(0.5 *) delta_heat_adv3 * dT * int{(u_c * grad{T}, v)} - (0.5) * delta_heat_adv4 * dT * int{(u_p * grad{T}, v)}
        for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
        {
            double tmp = 0.0;
            for ( int s = 0; s < DIM; ++s )
                tmp += ( this->delta_heat_adv_nc_ * this->solP_next_[s][q]
                    + this->delta_heat_adv_cc_ * this->solP_[s][q]
                    + this->delta_heat_adv_pc_ * this->solP_prev_[s][q] )
                * inv_r_comp[s] * this->grad_phi ( j, q, t_var )[s];

            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq
                        * skew_fac
                        * dt
                        * tmp
                        * this->phi ( i, q, t_var )
                        * dJ;
            }
        }

        // *****************************************************************************
        // reaction term: - delta_heat_rea3 * dT (div(u_c) T, v) - delta_heat_rea4 * dT (div(u_p) T, v)
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {

                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq
                        * skew_fac
                        * dt
                        * ( this->delta_heat_rea_nc_ * div_n + this->delta_heat_rea_cc_ * div_c + this->delta_heat_rea_pc_ * div_p )
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var )
                        * dJ;
            }
        }
        // *******************************************************************************
        // source term I: delta_heat_sou_2_ * dT * alpha_g * (grav * u, v)  
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int u_var = 0; u_var < DIM; u_var++ )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {

                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += wq
                            * dt
                            * this->delta_heat_sou_c_
                            * this->alpha_g_
                            * this->phi ( j, q, u_var )
                            * grav[u_var]
                            * this->phi ( i, q, t_var )
                            * dJ;
                }
            }
        }
    }
}
#endif
#ifdef OPT_ASM1

template<int DIM, class DataType>
void MetFlowBousCylDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( const Element<double>& element, LocalMatrix& lm ) const
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
    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    double dt = this->dT_cn_;
    if ( this->galerkin_mode_ == 1 )
        dt = this->dT_pc_;

    std::vector<double> grav;
    grav.resize ( 3, 0.0 );
    double sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    double skew_fac = 1.;
    double temp_supg = 0.;

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    const double mom_sou_pc = this->delta_mom_sou_pc_ * dt;
    const double mom_sou_nc = this->delta_mom_sou_nc_ * dt;
    const double mom_sou_cc = this->delta_mom_sou_cc_ * dt;
    const double heat_diff_c = this->kappa_ * dt * this->delta_heat_diff_c_;
    const double heat_adv_pc = skew_fac * dt * this->delta_heat_adv_pc_;
    const double heat_adv_nc = skew_fac * dt * this->delta_heat_adv_nc_;
    const double heat_adv_cc = skew_fac * dt * this->delta_heat_adv_cc_;
    const double heat_rea_pc = skew_fac * dt * this->delta_heat_rea_pc_;
    const double heat_rea_nc = skew_fac * dt * this->delta_heat_rea_nc_;
    const double heat_rea_cc = skew_fac * dt * this->delta_heat_rea_cc_;
    const double heat_sou_c = dt * this->delta_heat_sou_c_ * this->alpha_g_;

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const double r = this->x ( q )[1];
        const double inv_r = 1 / r;
        const double inv_rr = inv_r * inv_r;
        const double wq_dJ = this->w ( q ) * r * std::fabs ( this->detJ ( q ) );
        const double inv_r_comp[3] = { inv_r, 1., 1. };
        const double inv_rr_comp[3] = { inv_rr, 1., 1. };

#    ifdef ROTATING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;

#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif  

        double div_n = inv_r * this->solP_next_[1][q];
        double div_c = inv_r * this->solP_[1][q];
        double div_p = inv_r * this->solP_prev_[1][q];
        for ( int s = 0; s < DIM; s++ )
        {
            div_n += inv_r_comp[s] * this->grad_solP_next_[s][q][s];
            div_c += inv_r_comp[s] * this->grad_solP_[s][q][s];
            div_p += inv_r_comp[s] * this->grad_solP_prev_[s][q][s];
        }

        // ***********************************************************************
        // SOURCE TERM: -(-grad(T_P) T_D, v) 
        //  delta_mom_sou_3_ * dT * (grad(T_c) * T,v) + delta_mom_sou_4_ * dT * (grad(T_p) * T,v)  
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                double tmp = ( mom_sou_nc * this->grad_solP_[t_var]next_[q][v_var] + mom_sou_cc * this->grad_solP_[t_var][q][v_var] + mom_sou_pc * this->grad_solP_prev_[t_var][q][v_var] ) * inv_r_comp[v_var] * this->phi ( i, q, v_var );

                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, v_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                            * tmp
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
                // TIME-DERIVATIVE: -theta_d_dt_u * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * this->delta_d_dt_u_
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var );

                // ***********************************************************************
                // Thermal diffusion:  delta_heat_diff_p_ * kappa * dT * int{grad{u}:grad{v})
                double tmp = inv_rr * this->grad_phi ( i, q, t_var )[0] * this->grad_phi ( j, q, t_var )[0]
                        + this->grad_phi ( i, q, t_var )[1] * this->grad_phi ( j, q, t_var )[1];
                if ( DIM == 3 ) tmp += this->grad_phi ( i, q, t_var )[2] * this->grad_phi ( j, q, t_var )[2];

                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * heat_diff_c
                        * tmp;
            }
        }

        // *****************************************************************************
        // Nonlinear advection : -(0.5 *) delta_heat_adv3 * dT * int{(u_c * grad{T}, v)} - (0.5) * delta_heat_adv4 * dT * int{(u_p * grad{T}, v)}
        for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
        {
            double tmp = 0.0;
            for ( int s = 0; s < DIM; ++s )
                tmp += ( heat_adv_nc * this->solP_next_[s][q] + heat_adv_cc * this->solP_[s][q] + heat_adv_pc * this->solP_prev_[s][q] ) * inv_r_comp[s] * this->grad_phi ( j, q, t_var )[s];

            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq_dJ
                        * tmp
                        * this->phi ( i, q, t_var );
            }
        }

        // *****************************************************************************
        // reaction term: - delta_heat_rea3 * dT (div(u_c) T, v) - delta_heat_rea4 * dT (div(u_p) T, v)
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq_dJ
                        * ( heat_rea_nc * div_n + heat_rea_cc * div_c + heat_rea_pc * div_p )
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var );
            }
        }
        // *******************************************************************************
        // source term I: delta_heat_sou_2_ * dT * alpha_g * (grav * u, v)  
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int u_var = 0; u_var < DIM; u_var++ )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {

                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * heat_sou_c
                            * this->phi ( j, q, u_var )
                            * grav[u_var]
                            * this->phi ( i, q, t_var );
                }
            }
        }
    }
}
#endif
#ifdef OPT_ASM2

template<int DIM, class DataType>
void MetFlowBousCylDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( const Element<double>& element, LocalMatrix& lm ) const
{
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
    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    double dt = this->dT_cn_;
    if ( this->galerkin_mode_ == 1 )
        dt = this->dT_pc_;

    std::vector<double> grav;
    grav.resize ( 3, 0.0 );
    double sign = 1.0;
    const int num_q = this->num_quadrature_points ( );

    double skew_fac = 1.;
    double temp_supg = 0.;

    if ( this->skew_mode_ > 0 )
        skew_fac = 0.5;

    const double mom_sou_pc = this->delta_mom_sou_pc_ * dt;
    const double mom_sou_nc = this->delta_mom_sou_nc_ * dt;
    const double mom_sou_cc = this->delta_mom_sou_cc_ * dt;
    const double heat_diff_c = this->kappa_ * dt * this->delta_heat_diff_c_;
    const double heat_adv_pc = skew_fac * dt * this->delta_heat_adv_pc_;
    const double heat_adv_nc = skew_fac * dt * this->delta_heat_adv_nc_;
    const double heat_adv_cc = skew_fac * dt * this->delta_heat_adv_cc_;
    const double heat_rea_pc = skew_fac * dt * this->delta_heat_rea_pc_;
    const double heat_rea_nc = skew_fac * dt * this->delta_heat_rea_nc_;
    const double heat_rea_cc = skew_fac * dt * this->delta_heat_rea_cc_;
    const double heat_sou_c = dt * this->delta_heat_sou_c_ * this->alpha_g_;

    std::vector< std::vector<double> > phi ( num_vars );
    std::vector< std::vector< std::vector<double> > > grad_phi ( num_vars );

    for ( int k = 0; k < num_vars; ++k )
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
        const double r = this->x ( q )[1];
        const double inv_r = 1 / r;
        const double inv_rr = inv_r * inv_r;
        const double wq_dJ = this->w ( q ) * r * std::fabs ( this->detJ ( q ) );
        const double inv_r_comp[3] = { inv_r, 1., 1. };
        const double inv_rr_comp[3] = { inv_rr, 1., 1. };

        for ( int k = 0; k < num_vars; ++k )
        {
            if ( k != p_var && k != p0_var )
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

#    ifdef ROTATING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;
#    else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#    endif  

        double div_n = inv_r * this->solP_next_[1][q];
        double div_c = inv_r * this->solP_[1][q];
        double div_p = inv_r * this->solP_prev_[1][q];
        for ( int s = 0; s < DIM; s++ )
        {
            div_n += inv_r_comp[s] * this->grad_solP_next_[s][q][s];
            div_c += inv_r_comp[s] * this->grad_solP_[s][q][s];
            div_p += inv_r_comp[s] * this->grad_solP_prev_[s][q][s];
        }

        // ***********************************************************************
        // SOURCE TERM: -(-grad(T_P) T_D, v) 
        //  delta_mom_sou_3_ * dT * (grad(T_c) * T,v) + delta_mom_sou_4_ * dT * (grad(T_p) * T,v)  
        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i<this->num_dofs ( v_var ); ++i )
            {
                double tmp = ( mom_sou_nc * this->grad_solP_next_[t_var][q][v_var] + mom_sou_cc * this->grad_solP_[t_var][q][v_var] + mom_sou_pc * this->grad_solP_prev_[t_var][q][v_var] ) * inv_r_comp[v_var] * phi[v_var][i];

                for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
                {
                    lm ( this->dof_index ( i, v_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                            * tmp
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
                // TIME-DERIVATIVE: -theta_d_dt_u * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * this->delta_d_dt_u_
                        * phi[t_var][j]
                        * phi[t_var][i];

                // ***********************************************************************
                // Thermal diffusion:  delta_heat_diff_p_ * kappa * dT * int{grad{u}:grad{v})
                double tmp = inv_rr * grad_phi[t_var][i][0] * grad_phi[t_var][j][0]
                        + grad_phi[t_var][i][1] * grad_phi[t_var][j][1];
                if ( DIM == 3 ) tmp += grad_phi[t_var][i][2] * grad_phi[t_var][j][2];

                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq_dJ
                        * heat_diff_c
                        * tmp;
            }
        }

        // *****************************************************************************
        // Nonlinear advection : -(0.5 *) delta_heat_adv3 * dT * int{(u_c * grad{T}, v)} - (0.5) * delta_heat_adv4 * dT * int{(u_p * grad{T}, v)}
        for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
        {
            double tmp = 0.0;
            for ( int s = 0; s < DIM; ++s )
                tmp += ( heat_adv_nc * this->solP_next_[s][q] + heat_adv_cc * this->solP_[s][q] + heat_adv_pc * this->solP_prev_[s][q] ) * inv_r_comp[s] * grad_phi[t_var][j][s];

            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq_dJ
                        * tmp
                        * phi[t_var][i];
            }
        }

        // *****************************************************************************
        // reaction term: - delta_heat_rea3 * dT (div(u_c) T, v) - delta_heat_rea4 * dT (div(u_p) T, v)
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq_dJ
                        * ( heat_rea_nc * div_n + heat_rea_cc * div_c + heat_rea_pc * div_p )
                        * phi[t_var][j]
                        * phi[t_var][i];
            }
        }
        // *******************************************************************************
        // source term I: delta_heat_sou_2_ * dT * alpha_g * (grav * u, v)  
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int u_var = 0; u_var < DIM; u_var++ )
            {
                for ( int j = 0; j<this->num_dofs ( u_var ); ++j )
                {

                    lm ( this->dof_index ( i, t_var ), this->dof_index ( j, u_var ) ) += wq_dJ
                            * heat_sou_c
                            * phi[u_var][j]
                            * grav[u_var]
                            * phi[t_var][i];
                }
            }
        }
    }
}
#endif

/// Residual of dual problem (rhs - lhs)

template<int DIM, class DataType>
void MetFlowBousCylDualAssembler<DIM, DataType>::assemble_local_vector_dual ( const Element<double>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );
    const int p_var = DIM;

#ifdef AUGMENT_PRESS
    const int num_vars = DIM + 3;
    const int p0_var = DIM + 1;
    const int t_var = DIM + 2;
#else
    const int num_vars = DIM + 2;
    const int p0_var = -1;
    const int t_var = DIM + 1;
#endif

    double sign = 1.0;

    std::vector<double> grav;
    grav.resize ( 3, 0.0 );
    double div_n, div_c;

    double tau_div = 0.;
    double skew_fac = 1.;
    double temp_supg = 0.;

    if ( this->graddiv_mode_ > 0 )
        MetFlowBousCylDualAssembler<DIM, DataType>::compute_tau_graddiv ( element, tau_div );

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

    const double mom_sou_nn = this->delta_mom_sou_nn_ * dt_n;
    const double mom_sou_cn = this->delta_mom_sou_cn_ * dt_n;
    const double mom_sou_nc = this->delta_mom_sou_nc_ * dt_n;
    const double mom_sou_cc = this->delta_mom_sou_cc_ * dt_c;
    const double mom_sou_pc = this->delta_mom_sou_pc_ * dt_p;

    const double heat_diff_n = this->kappa_ * dt_n * this->delta_heat_diff_n_;
    const double heat_diff_c = this->kappa_ * dt_c * this->delta_heat_diff_c_;
    const double heat_adv_nn = skew_fac * dt_n * this->delta_heat_adv_nn_;
    const double heat_adv_cn = skew_fac * dt_n * this->delta_heat_adv_cn_;
    const double heat_adv_nc = skew_fac * dt_n * this->delta_heat_adv_nc_;
    const double heat_adv_cc = skew_fac * dt_c * this->delta_heat_adv_cc_;
    const double heat_adv_pc = skew_fac * dt_p * this->delta_heat_adv_pc_;

    const double heat_rea_nn = skew_fac * dt_n * this->delta_heat_rea_nn_;
    const double heat_rea_cn = skew_fac * dt_n * this->delta_heat_rea_cn_;
    const double heat_rea_nc = skew_fac * dt_n * this->delta_heat_rea_nc_;
    const double heat_rea_cc = skew_fac * dt_c * this->delta_heat_rea_cc_;
    const double heat_rea_pc = skew_fac * dt_p * this->delta_heat_rea_pc_;

    const double heat_sou_n = dt_n * this->delta_heat_sou_n_ * this->alpha_g_;
    const double heat_sou_c = dt_c * this->delta_heat_sou_c_ * this->alpha_g_;
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

        solP_c[p_var] = this->solP_[p_var][q];
        solP_c[t_var] = this->solP_[t_var][q];
        solD_c[p_var] = this->solD_[p_var][q];
        solD_c[t_var] = this->solD_[t_var][q];

        solP_p[p_var] = this->solP_prev_[p_var][q];
        solP_p[t_var] = this->solP_prev_[t_var][q];

        solP_n[p_var] = this->solP_next_[p_var][q];
        solP_n[t_var] = this->solP_next_[t_var][q];
        solD_n[p_var] = this->solD_next_[p_var][q];
        solD_n[t_var] = this->solD_next_[t_var][q];

#ifdef AUGMENT_PRESS
        solP_n[p0_var] = this->solP_next_[p0_var][q];
        solP_p[p0_var] = this->solP_prev_[p0_var][q];
        solP_c[p0_var] = this->solP_[p0_var][q];
        solD_n[p0_var] = this->solD_next_[p0_var][q];
        solD_c[p0_var] = this->solD_[p0_var][q];
#endif    

        std::vector< Vec<DIM, DataType> > grad_solP_n ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_solP_c ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_solP_p ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_solD_n ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_solD_c ( num_vars );

        for ( int var = 0; var < num_vars; var++ )
        {
            for ( int d = 0; d < DIM; d++ )
            {
                if ( var < p_var )
                {
                    grad_solP_n[var][d] = this->grad_solP_next_[var][q][d];
                    grad_solP_c[var][d] = this->grad_solP_[var][q][d];
                    grad_solP_p[var][d] = this->grad_solP_prev_[var][q][d];
                    grad_solD_n[var][d] = this->grad_solD_next_[var][q][d];
                    grad_solD_c[var][d] = this->grad_solD_[var][q][d];
                }
                if ( var == t_var )
                {
                    grad_solP_n[var][d] = this->grad_solP_next_[t_var][q][d];
                    grad_solP_c[var][d] = this->grad_solP_[t_var][q][d];
                    grad_solP_p[var][d] = this->grad_solP_prev_[t_var][q][d];
                    grad_solD_n[var][d] = this->grad_solD_next_[t_var][q][d];
                    grad_solD_c[var][d] = this->grad_solD_[t_var][q][d];
                }
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

#ifdef ROTATING_FOR
        // Centrifugal force -> modify gravity 
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
        grav[1] += pow ( this->omega_, 2.0 ) * r;
#else
        for ( int i = 0; i < DIM; i++ ) grav[i] = this->g_[i];
#endif 

        // **************************************************************
        // SOURCE TERM: (grad(TP) * TD, v)
        // (0.5 *) delta_mom_sou_1 * dT * int{ grad(TP_c) * TD_c * v } + (0.5 *) delta_mom_sou_2 * dT * int{ grad(TP_p)* TD_c * v }
        // (0.5 *) delta_mom_sou_3 * dT * int{ grad(TP_c) * TD_p * v } + (0.5 *) delta_mom_sou_4 * dT * int{ grad(TP_p)* TD_p * v }

        for ( int v_var = 0; v_var < DIM; ++v_var )
        {
            for ( int i = 0; i < this->num_dofs ( v_var ); ++i )
            {
                lv[this->dof_index ( i, v_var )] += wq_dJ
                        * inv_r_comp[v_var]
                        * ( mom_sou_nn * grad_solP_n[t_var][v_var] * solD_n[t_var] +
                        mom_sou_cn * grad_solP_c[t_var][v_var] * solD_n[t_var] +
                        mom_sou_nc * grad_solP_n[t_var][v_var] * solD_c[t_var] +
                        mom_sou_cc * grad_solP_c[t_var][v_var] * solD_c[t_var] +
                        mom_sou_pc * grad_solP_p[t_var][v_var] * solD_c[t_var] )
                        * this->phi ( i, q, v_var );
            }
        }

        double divD_n = inv_r * solD_n[1];
        for ( int s = 0; s < DIM; s++ )
            divD_n += inv_r_comp[s] * grad_solD_n[s][s];

        double divD_c = inv_r * solD_c[1];
        for ( int s = 0; s < DIM; s++ )
            divD_c += inv_r_comp[s] * grad_solD_c[s][s];

        // **********************************************************************
        // **********************************************************************
        // NOW THE TEMPERATURE PARTS ...

        // **********************************************************************
        // TIME DERIVATIVE: -theta_d_dt_u * ((TD_c - TD_p),v) 

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * this->delta_d_dt_u_ // 0. or 1. -> stationary or instationary configuration
                    * ( solD_c[t_var] - solD_n[t_var] )
                    * this->phi ( i, q, t_var );
        }

        // **********************************************************************
        // LAPLACE: delta_heat_diff_c_ * dT * kappa * int{grad{TD_c} : grad{v}} + delta_heat_diff_p * dT * kappa * int{grad{TD_p} : grad{v}}

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            double laplaceD_n = 0.0;
            for ( int s = 0; s < DIM; s++ ) laplaceD_n += inv_rr_comp[s] * this->grad_phi ( i, q, t_var )[s] * grad_solD_n[t_var][s];

            double laplaceD_c = 0.0;
            for ( int s = 0; s < DIM; s++ ) laplaceD_c += inv_rr_comp[s] * this->grad_phi ( i, q, t_var )[s] * grad_solD_c[t_var][s];

            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * ( heat_diff_n * laplaceD_n
                    + heat_diff_c * laplaceD_c );
        }

        // **********************************************************************
        // CONVECTIVE TERM (uP * grad(TD), v)

        // - (0.5) * delta_heat_adv_cc * dT * int{(uP_c * grad{TD_c}, v)} - (0.5) * delta_heat_adv_pc * dT * int{(uP_p * grad{TD_c}, v)} 
        // - (0.5) * delta_heat_adv_cp * dT * int{(uP_c * grad{TD_p}, v)} - (0.5) * delta_heat_adv_pp * dT * int{(uP_p * grad{TD_p}, v)} 
        double convection = 0.0;
        for ( int s = 0; s < DIM; ++s )
        {
            convection += ( solP_n[s] * grad_solD_n[t_var][s] * heat_adv_nn
                    + solP_c[s] * grad_solD_n[t_var][s] * heat_adv_cn
                    + solP_n[s] * grad_solD_c[t_var][s] * heat_adv_nc
                    + solP_c[s] * grad_solD_c[t_var][s] * heat_adv_cc
                    + solP_p[s] * grad_solD_c[t_var][s] * heat_adv_pc )
                    * inv_r_comp[s];
        }

        for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += -wq_dJ
                    * convection
                    * this->phi ( i, q, t_var );

        }

        double divP_n = inv_r * this->solP_next_[1][q];
        double divP_c = inv_r * this->solP_[1][q];
        double divP_p = inv_r * this->solP_prev_[1][q];
        for ( int s = 0; s < DIM; s++ )
        {
            divP_n += inv_r_comp[s] * this->grad_solP_next_[s][q][s];
            divP_c += inv_r_comp[s] * this->grad_solP_[s][q][s];
            divP_p += inv_r_comp[s] * this->grad_solP_prev_[s][q][s];
        }

        // **************************************************************
        // REACTION TERM (div(uP) * TD, v)

        for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += -wq_dJ
                    * ( heat_rea_nn * divP_n * solD_n[t_var] +
                    heat_rea_cn * divP_c * solD_n[t_var] +
                    heat_rea_nc * divP_n * solD_c[t_var] +
                    heat_rea_cc * divP_c * solD_c[t_var] +
                    heat_rea_pc * divP_p * solD_c[t_var] )
                    * this->phi ( i, q, t_var );
        }

        // **************************************************************
        // SOURCE TERM alpha_g (uD * grav, v)
        double gravD_n = 0.0;
        double gravD_c = 0.0;
        for ( int s = 0; s < DIM; s++ )
        {
            gravD_n += solD_n[s] * grav[s];
            gravD_c += solD_c[s] * grav[s];
        }
        for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * ( heat_sou_n * gravD_n + heat_sou_c * gravD_c )
                    * this->phi ( i, q, t_var );
        }

        // **********************************************************************
        // GOAL FUNCTIONAL: delta_j_c * dT * j_u(solP, phi) + delta_j_c * dT * j_u(solP_prev, phi)
        if ( this->goal_functional_->GetIntegralType ( ) == 1 )
        {
            gp.var = t_var;
            for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
            {
                gp.phi = this->phi ( i, q, t_var );
                gp.grad_phi = this->grad_phi ( i, q, t_var );

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
                lv[this->dof_index ( i, t_var )] += wq_dJ * ( j_c_c * J_c + j_c_pc * J_pc + j_n_c * J_c + j_n_nc * J_nc + j_n_n * J_n );
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

            for ( int var = 0; var < num_vars; var++ )
            {
                for ( int d = 0; d < DIM; d++ )
                {
                    gp.grad_solP[var][d] = grad_solP_c[var][d];
                }
            }
            gp.x = this->x ( q );
            gp.var = t_var;
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                gp.phi = this->phi ( i, q, t_var );
                gp.grad_phi = this->grad_phi ( i, q, t_var );
                lv[this->dof_index ( i, t_var )] += /*-*/wq_dJ * this->goal_functional_->j_final_type ( gp );
            }
        }
    }
}

template class MetFlowBousCylDualAssembler<2, double>;
template class MetFlowBousCylDualAssembler<3, double>;
