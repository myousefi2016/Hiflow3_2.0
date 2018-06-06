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

#include "met_flow_convdiff_cart_est_assembler.h"

template<int DIM, class DataType>
MetFlowConvDiffCartEstimatorAssembler<DIM, DataType>::MetFlowConvDiffCartEstimatorAssembler ( )
: MetFlowConvDiffEstimatorAssembler<DIM, DataType>( )
{

}

/// ********************************************************
/// General Assembly routines 
/// ********************************************************

template<int DIM, class DataType>
void MetFlowConvDiffCartEstimatorAssembler<DIM, DataType>::assemble_primal_cell_residual ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const
{
    const int num_q = AssemblyAssistant<DIM, DataType>::num_quadrature_points ( );
    const int t_var = 0;

    assert ( num_q == this->hess_solP_[0].size ( ) );

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq_dJ = AssemblyAssistant<DIM, DataType>::w ( q ) * std::abs ( AssemblyAssistant<DIM, DataType>::detJ ( q ) );

        // Residual: d_t(T) - kappa * Laplace(T) + gamma * conv * grad(T) + lambda * source * T 
        for ( int i = 0; i<this->quad_c_.size ( ); ++i )
        {
            const DataType c_i = this->c_i ( i );
            const DataType w_i = this->w_i ( i );

            //std::cout << c_i << " " << w_i << std::endl;

            double laplace = 0.;
            double conv = 0.;
            for ( int d = 0; d < DIM; ++d )
            {
                laplace += this->hess_trialP ( c_i, q, t_var )( d, d );
                conv += this->convection ( c_i, q, d ) * this->grad_trialP ( c_i, q, t_var )[d];
            }

            const DataType residual = ( this->theta_d_dt_u_ * this->dt_trialP ( c_i, q, t_var )
                    - this->kappa_ * laplace
                    + this->gamma_ * conv
                    + this->lambda_ * this->trialP ( c_i, q, t_var )
                    - this->source ( c_i, q, t_var ) );

            const DataType weight = this->weightD ( c_i, q, t_var, est_mode );

            if ( !this->use_csi_ )
            {
                r += residual // primal residual
                        * weight // dual weight function
                        * w_i
                        * wq_dJ;

                w = 1.;
            }
            else
            {
                r += residual * residual * w_i * wq_dJ;
                w += weight * weight * w_i * wq_dJ;
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCartEstimatorAssembler<DIM, DataType>::assemble_dual_cell_residual ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const
{
    const int num_q = AssemblyAssistant<DIM, DataType>::num_quadrature_points ( );
    const int t_var = 0;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq_dJ = AssemblyAssistant<DIM, DataType>::w ( q ) * std::abs ( AssemblyAssistant<DIM, DataType>::detJ ( q ) );

        // Residual: d_t(T) - kappa * Laplace(T) + gamma * conv * grad(T) + lambda * source * T 
        for ( int i = 0; i<this->quad_c_.size ( ); ++i )
        {
            const DataType c_i = this->c_i ( i );
            const DataType w_i = this->w_i ( i );
            const DataType t_i = this->get_absolute_time ( c_i );

            // setup parameter struct for goal functional
            ParametersEvalType<DIM, DataType> gp;
            gp.var = t_var;
            gp.x = AssemblyAssistant<DIM, DataType>::x ( q );
            gp.absolute_time = t_i;
            gp.solP = std::vector<DataType> ( 1, this->trialP ( c_i, q, t_var ) );
            gp.grad_solP = std::vector<Vec<DIM, DataType> >( 1, this->grad_trialP ( c_i, q, t_var ) );
            gp.d_solP = std::vector<DataType> ( 1, this->weightP ( c_i, q, t_var, est_mode ) );
            //          gp.grad_d_solP   = this->grad_weightD(c_i, q, t_var, est_mode);

            DataType laplace = 0.;
            DataType conv = 0.;
            DataType div = 0.;
            for ( int d = 0; d < DIM; ++d )
            {
                laplace += this->hess_trialD ( c_i, q, t_var )( d, d );
                conv += this->convection ( c_i, q, d ) * this->grad_trialD ( c_i, q, t_var )[d];
                div += this->grad_convection ( c_i, q, d )[d];
            }

            const DataType residual = ( -this->delta_d_dt_u_ * this->dt_trialD ( c_i, q, t_var )
                    - this->kappa_ * laplace
                    - this->gamma_ * conv
                    + ( this->lambda_ - div ) * this->trialD ( c_i, q, t_var ) );
            const DataType goal = this->goal_functional_->j_force_eval_deriv ( gp );
            const DataType weight = this->weightP ( c_i, q, t_var, est_mode );
            if ( !this->use_csi_ )
            {
                r += ( residual * weight + goal )
                        * w_i
                        * wq_dJ;

                w = 1.;
            }
            else
            {
                r += ( residual * residual + std::abs ( goal ) ) * w_i * wq_dJ;
                w += ( weight * weight + std::abs ( goal ) ) * w_i * wq_dJ;
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCartEstimatorAssembler<DIM, DataType>::assemble_dual_cell_time_jump ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const
{
    const int num_q = AssemblyAssistant<DIM, DataType>::num_quadrature_points ( );
    const int t_var = 0;
    const DataType scaling_factor = 0.5;

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq_dJ = AssemblyAssistant<DIM, DataType>::w ( q ) * std::abs ( AssemblyAssistant<DIM, DataType>::detJ ( q ) );

        // setup parameter struct for goal functional
        if ( this->final_time_interval_ )
        {
            ParametersEvalType<DIM, DataType> gp;
            gp.var = t_var;
            gp.x = AssemblyAssistant<DIM, DataType>::x ( q );
            gp.absolute_time = 0.;
            gp.solP = std::vector<DataType> ( 1, this->trialP ( 1., q, t_var ) );
            gp.grad_solP = std::vector<Vec<DIM, DataType> >( 1, this->grad_trialP ( 1., q, t_var ) );
            gp.d_solP = std::vector<DataType> ( 1, this->weightP ( 1., q, t_var, est_mode ) );
            //          gp.grad_d_solP      = this->grad_weightD(c_i, q, t_var, est_mode);
            DataType j_val = this->goal_functional_->j_final_eval_deriv ( gp );

            r += scaling_factor
                    * ( j_val - this->solD_[t_var][q] * this->weightP ( 0., q, t_var, est_mode ) ) // TODO: c_i checken
                    * wq_dJ;
        }
        else
        {
            const DataType residual = this->solD_next_[t_var][q] - this->solD_[t_var][q];
            const DataType weight = this->weightP ( 0., q, t_var, est_mode ); // TODO: c_i checken

            if ( !this->use_csi_ )
            {
                r += scaling_factor
                        * residual
                        * weight
                        * wq_dJ;

                w = 1.;
            }
            else
            {
                r += scaling_factor * residual * residual * wq_dJ;
                w += scaling_factor * weight * weight * wq_dJ;
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCartEstimatorAssembler<DIM, DataType>::assemble_primal_interface_jump ( const Element<DataType>& left_elem,
                                                                                            const Element<DataType>& right_elem,
                                                                                            const Quadrature<DataType>& left_quad,
                                                                                            const Quadrature<DataType>& right_quad,
                                                                                            int left_facet_number,
                                                                                            int right_facet_number,
                                                                                            InterfaceSide left_if_side,
                                                                                            InterfaceSide right_if_side,
                                                                                            int eq,
                                                                                            int est_mode,
                                                                                            DataType& r,
                                                                                            DataType& w )
{
    const int num_q = DGAssemblyAssistant<DIM, DataType>::num_quadrature_points ( );
    assert ( left_if_side != right_if_side );

    const int t_var = 0;
    const DataType scaling_factor = 0.5;

    // Loop over quadrature points on each edge
    for ( int q = 0.; q < num_q; ++q )
    {
        const DataType wq = DGAssemblyAssistant<DIM, DataType>::w ( q );
        const DataType dS = std::abs ( DGAssemblyAssistant<DIM, DataType>::ds ( q ) );

        for ( int i = 0; i<this->quad_c_.size ( ); ++i )
        {
            const DataType c_i = this->c_i ( i );
            const DataType w_i = this->w_i ( i );
            const DataType t_i = this->get_absolute_time ( c_i );

            Vec<DIM, DataType> n_e;
            Vec<DIM, DataType> grad_jump;

            if ( left_if_side == DGGlobalAssembler<DataType>::INTERFACE_MASTER )
            {
                n_e = this->trial ( ).n ( q );
                grad_jump = this->left_grad_trialP ( c_i, q, t_var ) - this->right_grad_trialP ( c_i, q, t_var ); // master - slave
            }
            else
            {
                n_e = this->test ( ).n ( q );
                grad_jump = this->right_grad_trialP ( c_i, q, t_var ) - this->left_grad_trialP ( c_i, q, t_var );
            }

            const DataType residual = this->kappa_ * dot ( n_e, grad_jump );
            const DataType weight = this->weightD ( c_i, q, t_var, est_mode );

            if ( !this->use_csi_ )
            {
                r += scaling_factor /*  * this->h_E_ */ * residual * weight * wq * w_i * dS;
                w = scaling_factor;
            }
            else
            {
                r += scaling_factor /*  * this->h_E_ */ * residual * residual * wq * w_i * dS;
                w += scaling_factor /*  * this->h_E_ */ * weight * weight * wq * w_i * dS;
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCartEstimatorAssembler<DIM, DataType>::assemble_dual_interface_jump ( const Element<DataType>& left_elem,
                                                                                          const Element<DataType>& right_elem,
                                                                                          const Quadrature<DataType>& left_quad,
                                                                                          const Quadrature<DataType>& right_quad,
                                                                                          int left_facet_number,
                                                                                          int right_facet_number,
                                                                                          InterfaceSide left_if_side,
                                                                                          InterfaceSide right_if_side,
                                                                                          int eq,
                                                                                          int est_mode,
                                                                                          DataType& r,
                                                                                          DataType& w )
{
    assert ( left_if_side != right_if_side );

    const int num_q = DGAssemblyAssistant<DIM, DataType>::num_quadrature_points ( );
    const int t_var = 0;
    const DataType scaling_factor = 0.5;

    // Loop over quadrature points on each edge
    for ( int q = 0.; q < num_q; ++q )
    {
        const DataType wq = DGAssemblyAssistant<DIM, DataType>::w ( q );
        const DataType dS = std::abs ( DGAssemblyAssistant<DIM, DataType>::ds ( q ) );

        for ( int i = 0; i<this->quad_c_.size ( ); ++i )
        {
            const DataType c_i = this->c_i ( i );
            const DataType w_i = this->w_i ( i );
            const DataType t_i = this->get_absolute_time ( c_i );

            Vec<DIM, DataType> n_e;
            Vec<DIM, DataType> grad_jump;

            if ( left_if_side == DGGlobalAssembler<DataType>::INTERFACE_MASTER )
            {
                n_e = this->trial ( ).n ( q );
                grad_jump = this->left_grad_trialD ( c_i, q, t_var ) - this->right_grad_trialD ( c_i, q, t_var ); // master - slave
            }
            else
            {
                n_e = this->test ( ).n ( q );
                grad_jump = this->right_grad_trialD ( c_i, q, t_var ) - this->left_grad_trialD ( c_i, q, t_var );
            }

            const DataType residual = this->kappa_ * dot ( n_e, grad_jump );
            const DataType weight = this->weightP ( c_i, q, t_var, est_mode );

            if ( !this->use_csi_ )
            {
                r += scaling_factor /*  * this->h_E_ */ * residual * weight * wq * w_i * dS;
                w = scaling_factor;
            }
            else
            {
                r += scaling_factor /*  * this->h_E_ */ * residual * residual * wq * w_i * dS;
                w += scaling_factor /*  * this->h_E_ */ * weight * weight * wq * w_i * dS;
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCartEstimatorAssembler<DIM, DataType>::assemble_primal_interface_boundary ( const Element<DataType>& left_elem,
                                                                                                const Element<DataType>& right_elem,
                                                                                                const Quadrature<DataType>& left_quad,
                                                                                                const Quadrature<DataType>& right_quad,
                                                                                                int left_facet_number,
                                                                                                int right_facet_number,
                                                                                                InterfaceSide left_if_side,
                                                                                                InterfaceSide right_if_side,
                                                                                                int eq,
                                                                                                int est_mode,
                                                                                                DataType& r,
                                                                                                DataType& w )
{
    const int num_q = DGAssemblyAssistant<DIM, DataType>::num_quadrature_points ( );
    const int t_var = 0;
    const DataType scaling_factor = 1.;

    // Loop over quadrature points on each edge
    for ( int q = 0.; q < num_q; ++q )
    {
        const DataType wq = DGAssemblyAssistant<DIM, DataType>::w ( q );
        const DataType dS = std::abs ( DGAssemblyAssistant<DIM, DataType>::ds ( q ) );

        for ( int i = 0; i<this->quad_c_.size ( ); ++i )
        {
            const DataType c_i = this->c_i ( i );
            const DataType w_i = this->w_i ( i );
            const DataType t_i = this->get_absolute_time ( c_i );

            const DataType residual = this->kappa_ * dot ( this->test ( ).n ( q ), this->right_grad_trialP ( c_i, q, t_var ) );
            const DataType weight = this->weightD ( c_i, q, t_var, est_mode );

            if ( !this->use_csi_ )
            {
                r += scaling_factor /*  * this->h_E_ */ * wq * w_i * residual * weight * dS;
                w = scaling_factor;
            }
            else
            {
                r += scaling_factor /*  * this->h_E_ */ * wq * w_i * residual * residual * dS;
                w += scaling_factor /*  * this->h_E_ */ * wq * w_i * weight * weight * dS;
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowConvDiffCartEstimatorAssembler<DIM, DataType>::assemble_dual_interface_boundary ( const Element<DataType>& left_elem,
                                                                                              const Element<DataType>& right_elem,
                                                                                              const Quadrature<DataType>& left_quad,
                                                                                              const Quadrature<DataType>& right_quad,
                                                                                              int left_facet_number,
                                                                                              int right_facet_number,
                                                                                              InterfaceSide left_if_side,
                                                                                              InterfaceSide right_if_side,
                                                                                              int eq,
                                                                                              int est_mode,
                                                                                              DataType& r,
                                                                                              DataType& w )
{
    const int num_q = DGAssemblyAssistant<DIM, DataType>::num_quadrature_points ( );
    const int t_var = 0;
    const DataType scaling_factor = 1.;

    // Loop over quadrature points on each edge
    for ( int q = 0.; q < num_q; ++q )
    {
        const DataType wq = DGAssemblyAssistant<DIM, DataType>::w ( q );
        const DataType dS = std::abs ( DGAssemblyAssistant<DIM, DataType>::ds ( q ) );

        for ( int i = 0; i<this->quad_c_.size ( ); ++i )
        {
            const DataType c_i = this->c_i ( i );
            const DataType w_i = this->w_i ( i );
            const DataType t_i = this->get_absolute_time ( c_i );

            const DataType residual = this->kappa_ * dot ( this->test ( ).n ( q ), this->right_grad_trialD ( c_i, q, t_var ) );
            const DataType weight = this->weightP ( c_i, q, t_var, est_mode );

            if ( !this->use_csi_ )
            {
                r += scaling_factor /*  * this->h_E_ */ * wq * w_i * residual * weight * dS;
                w = scaling_factor;
            }
            else
            {
                r += scaling_factor /*  * this->h_E_ */ * wq * w_i * residual * residual * dS;
                w += scaling_factor /*  * this->h_E_ */ * wq * w_i * weight * weight * dS;
            }
        }
    }
}

template class MetFlowConvDiffCartEstimatorAssembler<2, double>;
template class MetFlowConvDiffCartEstimatorAssembler<3, double>;

