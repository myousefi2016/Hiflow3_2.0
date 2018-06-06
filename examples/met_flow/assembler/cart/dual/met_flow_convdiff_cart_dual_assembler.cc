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

#include "met_flow_convdiff_cart_dual_assembler.h"

template<int DIM, class DataType>
MetFlowConvDiffCartDualAssembler<DIM, DataType>::MetFlowConvDiffCartDualAssembler ( )
: MetFlowConvDiffAssembler<DIM, DataType>( )
{
}

/// ********************************************************
/// General Assembly routines 
/// ********************************************************

template<int DIM, class DataType>
void MetFlowConvDiffCartDualAssembler<DIM, DataType>::assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    if ( this->mode_ == DUAL )
    {
        MetFlowConvDiffCartDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( element, lm );
    }
    else
    {
        interminable_assert ( 0 );
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCartDualAssembler<DIM, DataType>::assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    if ( this->vector_asm_mode_ == VECTOR_STD )
    {
        if ( this->mode_ == DUAL )
        {
            MetFlowConvDiffCartDualAssembler<DIM, DataType>::assemble_local_vector_dual ( element, lv );
        }
        else
        {
            interminable_assert ( 0 );
        }
    }
    else
    {
        interminable_assert ( 0 );
    }
}

/// ********************************************************
/// Assembly routines for primal problem
/// ********************************************************

template<int DIM, class DataType>
void MetFlowConvDiffCartDualAssembler<DIM, DataType>::assemble_local_matrix_dual ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int t_var = 0;
    const DataType dt = this->dT_pc_;
    const int num_q = this->num_quadrature_points ( );

    // loop over quadrature points  
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::fabs ( this->detJ ( q ) );

        double div_n = 0.;
        double div_c = 0.;
        double div_p = 0.;
        for ( int s = 0; s < DIM; s++ )
        {
            div_n += this->grad_conv_next_[s][q][s];
            div_c += this->grad_conv_[s][q][s];
            div_p += this->grad_conv_prev_[s][q][s];
        }

        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
            {

                // ***********************************************************************
                // TIME-DERIVATIVE: theta_d_dt_u * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * this->delta_d_dt_u_
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var )
                        * dJ;

                // ***********************************************************************
                // Thermal diffusion: theta1 * kappa * dT * int{grad{u}:grad{v})
                DataType tmp = this->grad_phi ( i, q, t_var )[0] * this->grad_phi ( j, q, t_var )[0]
                        + this->grad_phi ( i, q, t_var )[1] * this->grad_phi ( j, q, t_var )[1];
                if ( DIM == 3 ) tmp += this->grad_phi ( i, q, t_var )[2] * this->grad_phi ( j, q, t_var )[2];

                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * dt
                        * this->delta_diff_c_
                        * this->kappa_
                        * tmp
                        * dJ;

                // ***********************************************************************
                // reaction: lambda * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * dt
                        * this->lambda_
                        * this->delta_rea_c_
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var )
                        * dJ;
            }
        }

        // *****************************************************************************
        // Nonlinear advection : -(0.5 *) delta_heat_adv3 * dT * int{(u_c * grad{T}, v)} - (0.5) * delta_heat_adv4 * dT * int{(u_p * grad{T}, v)}
        for ( int j = 0; j<this->num_dofs ( t_var ); ++j )
        {
            double tmp = 0.0;
            for ( int s = 0; s < DIM; ++s )
                tmp += ( this->delta_adv_nc_ * this->conv_next_[s][q]
                    + this->delta_adv_cc_ * this->conv_[s][q]
                    + this->delta_adv_pc_ * this->conv_prev_[s][q] )
                * this->grad_phi ( j, q, t_var )[s];

            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += -wq
                        * this->gamma_
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
                        * this->gamma_
                        * dt
                        * ( this->delta_adv_nc_ * div_n + this->delta_adv_cc_ * div_c + this->delta_adv_pc_ * div_p )
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var )
                        * dJ;
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCartDualAssembler<DIM, DataType>::assemble_local_vector_dual ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int num_q = this->num_quadrature_points ( );
    const int t_var = 0;
    const int num_vars = 1;

    DataType sign = 1.0;

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

    const double diff_n = dt_n * this->delta_diff_n_;
    const double diff_c = dt_c * this->delta_diff_c_;
    const double adv_nn = dt_n * this->delta_adv_nn_;
    const double adv_cn = dt_n * this->delta_adv_cn_;
    const double adv_nc = dt_n * this->delta_adv_nc_;
    const double adv_cc = dt_c * this->delta_adv_cc_;
    const double adv_pc = dt_p * this->delta_adv_pc_;
    const double rea_n = dt_n * this->delta_rea_n_;
    const double rea_c = dt_c * this->delta_rea_c_;

    const double j_c_c = dt_p * this->delta_j_c_c_;
    const double j_c_pc = dt_p * this->delta_j_c_pc_;
    const double j_n_c = dt_n * this->delta_j_n_c_;
    const double j_n_nc = dt_n * this->delta_j_n_nc_;
    const double j_n_n = dt_n * this->delta_j_n_n_;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq_dJ = this->w ( q ) * std::abs ( this->detJ ( q ) );

        // get previous newton step solution in vector form
        // ns: Newton index, ts: time stepping index
        std::vector<double> solP_c ( num_vars );
        std::vector<double> solP_n ( num_vars );
        std::vector<double> solP_p ( num_vars );
        std::vector<double> solD_c ( num_vars );
        std::vector<double> solD_n ( num_vars );

        solP_c[t_var] = this->solP_[t_var][q];
        solD_c[t_var] = this->solD_[t_var][q];
        solP_p[t_var] = this->solP_prev_[t_var][q];
        solP_n[t_var] = this->solP_next_[t_var][q];
        solD_n[t_var] = this->solD_next_[t_var][q];

        std::vector< Vec<DIM, DataType> > grad_solP_n ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_solP_c ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_solP_p ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_solD_n ( num_vars );
        std::vector< Vec<DIM, DataType> > grad_solD_c ( num_vars );

        for ( int d = 0; d < DIM; d++ )
        {
            grad_solP_n[t_var][d] = this->grad_solP_next_[t_var][q][d];
            grad_solP_c[t_var][d] = this->grad_solP_[t_var][q][d];
            grad_solP_p[t_var][d] = this->grad_solP_prev_[t_var][q][d];
            grad_solD_n[t_var][d] = this->grad_solD_next_[t_var][q][d];
            grad_solD_c[t_var][d] = this->grad_solD_[t_var][q][d];
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

        double div_n = 0.;
        double div_c = 0.;
        double div_p = 0.;
        for ( int s = 0; s < DIM; s++ )
        {
            div_n += this->grad_conv_next_[s][q][s];
            div_c += this->grad_conv_[s][q][s];
            div_p += this->grad_conv_prev_[s][q][s];
        }

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
            for ( int s = 0; s < DIM; s++ ) laplaceD_n += this->grad_phi ( i, q, t_var )[s] * grad_solD_n[t_var][s];

            double laplaceD_c = 0.0;
            for ( int s = 0; s < DIM; s++ ) laplaceD_c += this->grad_phi ( i, q, t_var )[s] * grad_solD_c[t_var][s];

            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * this->kappa_
                    * ( diff_n * laplaceD_n
                    + diff_c * laplaceD_c );
        }

        // **********************************************************************
        // Reaction:        
        for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += wq_dJ
                    * this->lambda_
                    * ( rea_n * solD_n[t_var]
                    + rea_c * solD_c[t_var] )
                    * this->phi ( i, q, t_var );
        }

        // **********************************************************************
        // CONVECTIVE TERM (uP * grad(TD), v)

        // - (0.5) * delta_heat_adv_cc * dT * int{(uP_c * grad{TD_c}, v)} - (0.5) * delta_heat_adv_pc * dT * int{(uP_p * grad{TD_c}, v)} 
        // - (0.5) * delta_heat_adv_cp * dT * int{(uP_c * grad{TD_p}, v)} - (0.5) * delta_heat_adv_pp * dT * int{(uP_p * grad{TD_p}, v)} 
        double convection = 0.0;
        for ( int s = 0; s < DIM; ++s )
        {
            convection += ( this->conv_next_[s][q] * grad_solD_n[t_var][s] * adv_nn
                    + this->conv_[s][q] * grad_solD_n[t_var][s] * adv_cn
                    + this->conv_next_[s][q] * grad_solD_c[t_var][s] * adv_nc
                    + this->conv_[s][q] * grad_solD_c[t_var][s] * adv_cc
                    + this->conv_prev_[s][q] * grad_solD_c[t_var][s] * adv_pc );
        }

        for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += -wq_dJ
                    * convection
                    * this->gamma_
                    * this->phi ( i, q, t_var );

        }

        // **************************************************************
        // REACTION TERM (div(uP) * TD, v)

        for ( int i = 0; i < this->num_dofs ( t_var ); ++i )
        {
            lv[this->dof_index ( i, t_var )] += -wq_dJ
                    * ( adv_nn * div_n * solD_n[t_var] +
                    adv_cn * div_c * solD_n[t_var] +
                    adv_nc * div_n * solD_c[t_var] +
                    adv_cc * div_c * solD_c[t_var] +
                    adv_pc * div_p * solD_c[t_var] )
                    * this->gamma_
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
                /*
                DataType tmp = (j_c_c * J_c + j_c_pc * J_pc + j_n_c * J_c + j_n_nc * J_nc + j_n_n * J_n);
                if (tmp != 0.)
                {
                   std::cout << this->x(q)[0] << ", " << this->x(q)[1] << " :: " << tmp << std::endl;
                }*/
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
            }/*
            if (this->goal_functional_->j_final_type(gp) != 0.)
            { 
                std::cout << "FINAL " << this->goal_functional_->j_final_type(gp) << std::endl;
            }*/
        }
    }
}

template class MetFlowConvDiffCartDualAssembler<2, double>;
template class MetFlowConvDiffCartDualAssembler<3, double>;

