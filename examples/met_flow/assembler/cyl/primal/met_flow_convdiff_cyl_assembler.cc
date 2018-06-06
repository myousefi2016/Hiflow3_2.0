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

#include "met_flow_convdiff_cyl_assembler.h"

template<int DIM, class DataType>
MetFlowConvDiffCylAssembler<DIM, DataType>::MetFlowConvDiffCylAssembler ( )
: MetFlowConvDiffAssembler<DIM, DataType>( )
{
}

/// ********************************************************
/// General Assembly routines 
/// ********************************************************

template<int DIM, class DataType>
void MetFlowConvDiffCylAssembler<DIM, DataType>::assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    if ( this->mode_ == PRIMAL )
    {
        MetFlowConvDiffCylAssembler<DIM, DataType>::assemble_local_matrix_primal ( element, lm );
    }
    else
    {
        interminable_assert ( 0 );
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCylAssembler<DIM, DataType>::assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const
{
    const int total_dofs = this->num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    if ( this->vector_asm_mode_ == VECTOR_STD )
    {
        if ( this->mode_ == PRIMAL )
        {
            MetFlowConvDiffCylAssembler<DIM, DataType>::assemble_local_vector_primal ( element, lv );
        }
        else
        {
            interminable_assert ( 0 );
        }
    }
    else if ( this->vector_asm_mode_ == VECTOR_GOAL )
    {
        MetFlowConvDiffCylAssembler<DIM, DataType>::assemble_local_vector_goal ( element, lv );
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCylAssembler<DIM, DataType>::assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const
{
    if ( this->sca_mode_ == GOAL_INT )
    {
        MetFlowConvDiffAssembler<DIM, DataType>::assemble_local_scalar_goal_int ( element, ls );
    }
    if ( this->sca_mode_ == GOAL_FIN )
    {
        MetFlowConvDiffAssembler<DIM, DataType>::assemble_local_scalar_goal_fin ( element, ls );
    }
}

/// ********************************************************
/// Assembly routines for primal problem
/// ********************************************************

template<int DIM, class DataType>
void MetFlowConvDiffCylAssembler<DIM, DataType>::assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int t_var = 0;
    const DataType dt = this->dT_pc_;
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
                        * this->theta_diff_c_
                        * this->kappa_
                        * tmp
                        * dJ;

                // ***********************************************************************
                // reaction: lambda * int{(u,v)} 
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * dt
                        * this->lambda_
                        * this->theta_rea_c_
                        * this->phi ( j, q, t_var )
                        * this->phi ( i, q, t_var )
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
                tmp += ( this->theta_adv_cc_ * this->conv_[s][q] + this->theta_adv_pc_ * this->conv_prev_[s][q] )
                        * inv_r_comp[s] * this->grad_phi ( j, q, t_var )[s];
            }
            for ( int i = 0; i<this->num_dofs ( t_var ); ++i )
            {
                lm ( this->dof_index ( i, t_var ), this->dof_index ( j, t_var ) ) += wq
                        * dt
                        * tmp
                        * this->gamma_
                        * this->phi ( i, q, t_var )
                        * dJ;
            }
        }

    }
}

template class MetFlowConvDiffCylAssembler<2, double>;
template class MetFlowConvDiffCylAssembler<3, double>;

