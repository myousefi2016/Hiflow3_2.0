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

#include "met_flow_est_assembler.h"
#include "met_flow_assembler.h"

template<int DIM, class DataType>
MetFlowEstimatorAssembler<DIM, DataType>::MetFlowEstimatorAssembler ( )
: DGAssemblyAssistant<DIM, DataType>( ),
space_fineP_ ( NULL ),
space_fineP_prev_ ( NULL ),
space_fineP_next_ ( NULL ),
space_fineD_ ( NULL ),
space_fineD_next_ ( NULL ),
space_fineD_prev_ ( NULL ),
vector_fineP_ ( NULL ),
vector_fineP_prev_ ( NULL ),
vector_fineP_next_ ( NULL ),
vector_fineD_ ( NULL ),
vector_fineD_next_ ( NULL ),
vector_fineD_prev_ ( NULL ),
search_fineP_ ( NULL ),
search_fineP_prev_ ( NULL ),
search_fineP_next_ ( NULL ),
search_fineD_ ( NULL ),
search_fineD_prev_ ( NULL ),
search_fineD_next_ ( NULL ),
assemble_primal_indicator_ ( false ),
assemble_dual_indicator_ ( false ),
assemble_temporal_indicator_ ( false ),
assemble_spatial_indicator_ ( false ),
final_time_interval_ ( false ),
unique_fine_space_ ( false ),
use_csi_ ( false ),
quad_order_ ( 0 ),
num_eq_ ( 0 ),
time_offset_ ( 0. ),
use_dwr_ ( true )
{
    this->num_var_ = 0;
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::clear ( )
{
    MetFlowAssembler<DIM, DataType>::clear ( );

    this->quad_c_.clear ( );
    this->quad_w_.clear ( );
    this->quad_order_ = 0;

    this->est_rel_time_ = 0;
    this->t_ = 0.;
    this->final_time_interval_ = false;
    this->time_offset_ = 0.;

    this->assemble_primal_indicator_ = false;
    this->assemble_dual_indicator_ = false;
    this->assemble_temporal_indicator_ = false;
    this->assemble_spatial_indicator_ = false;
    this->use_csi_ = false;
    this->use_dwr_ = true;

    this->vector_fineP_ = NULL;
    this->vector_fineP_prev_ = NULL;
    this->vector_fineP_next_ = NULL;
    this->vector_fineD_ = NULL;
    this->vector_fineD_prev_ = NULL;
    this->vector_fineD_next_ = NULL;

    this->space_fineP_ = NULL;
    this->space_fineP_next_ = NULL;
    this->space_fineP_prev_ = NULL;

    this->space_fineD_ = NULL;
    this->space_fineD_next_ = NULL;
    this->space_fineD_prev_ = NULL;

    if ( this->search_fineP_ != NULL )
    {
        delete this->search_fineP_;
        this->search_fineP_ = NULL;
    }
    if ( this->search_fineP_prev_ != NULL )
    {
        delete this->search_fineP_prev_;
        this->search_fineP_prev_ = NULL;
    }
    if ( this->search_fineP_next_ != NULL )
    {
        delete this->search_fineP_next_;
        this->search_fineP_next_ = NULL;
    }
    if ( this->search_fineD_ != NULL )
    {
        delete this->search_fineD_;
        this->search_fineD_ = NULL;
    }
    if ( this->search_fineD_prev_ != NULL )
    {
        delete this->search_fineD_prev_;
        this->search_fineD_prev_ = NULL;
    }
    if ( this->search_fineD_next_ != NULL )
    {
        delete this->search_fineD_next_;
        this->search_fineD_next_ = NULL;
    }

    // TODO evtl delete verwenden
    this->funP_.clear ( );
    this->funP_prev_.clear ( );
    this->funP_next_.clear ( );
    this->funD_.clear ( );
    this->funD_prev_.clear ( );
    this->funD_next_.clear ( );

    this->num_eq_ = 0;
    this->num_var_ = 0;

    this->h_E_ = 0.;
    this->h_K_ = 0.;

    this->fineP_.clear ( );
    this->fineP_next_.clear ( );
    this->fineP_prev_.clear ( );

    this->fineD_.clear ( );
    this->fineD_prev_.clear ( );
    this->fineD_next_.clear ( );

    this->left_grad_solP_.clear ( );
    this->left_grad_solP_prev_.clear ( );
    this->left_grad_solP_next_.clear ( );
    this->right_grad_solP_.clear ( );
    this->right_grad_solP_prev_.clear ( );
    this->right_grad_solP_next_.clear ( );

    this->left_grad_solD_.clear ( );
    this->left_grad_solD_next_.clear ( );
    this->left_grad_solD_prev_.clear ( );
    this->right_grad_solD_.clear ( );
    this->right_grad_solD_next_.clear ( );
    this->right_grad_solD_prev_.clear ( );

    this->time_order_P_.clear ( );
    this->time_order_D_.clear ( );
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_values ( int num_var )
{
    MetFlowAssembler<DIM, DataType>::allocate_function_values ( num_var );

    this->fineP_.resize ( num_var );
    this->fineP_next_.resize ( num_var );
    this->fineP_prev_.resize ( num_var );

    this->fineD_.resize ( num_var );
    this->fineD_next_.resize ( num_var );
    this->fineD_prev_.resize ( num_var );

    this->left_grad_solP_.resize ( num_var );
    this->left_grad_solP_prev_.resize ( num_var );
    this->left_grad_solP_next_.resize ( num_var );
    this->right_grad_solP_.resize ( num_var );
    this->right_grad_solP_prev_.resize ( num_var );
    this->right_grad_solP_next_.resize ( num_var );

    this->left_grad_solD_.resize ( num_var );
    this->left_grad_solD_next_.resize ( num_var );
    this->left_grad_solD_prev_.resize ( num_var );
    this->right_grad_solD_.resize ( num_var );
    this->right_grad_solD_next_.resize ( num_var );
    this->right_grad_solD_prev_.resize ( num_var );
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_evaluators ( int num_var )
{
    this->funP_.resize ( num_var, NULL );
    this->funP_prev_.resize ( num_var, NULL );
    this->funP_next_.resize ( num_var, NULL );
    this->funD_.resize ( num_var, NULL );
    this->funD_next_.resize ( num_var, NULL );
    this->funD_prev_.resize ( num_var, NULL );
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::set_indicators ( bool primal, bool dual, bool temporal, bool spatial, bool csi )
{
    this->assemble_primal_indicator_ = primal;
    this->assemble_dual_indicator_ = dual;
    this->assemble_temporal_indicator_ = temporal;
    this->assemble_spatial_indicator_ = spatial;
    this->use_csi_ = csi;
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::setup_grid_search ( )
{
    MetFlowAssembler<DIM, DataType>::setup_grid_search ( );
    if ( this->space_fineP_ != NULL )
    {
        if ( this->search_fineP_ != NULL )
        {
            delete this->search_fineP_;
            this->search_fineP_ = NULL;
        }
        this->search_fineP_ = new GridGeometricSearch ( this->space_fineP_->meshPtr ( ) );
    }
    if ( this->space_fineP_prev_ != NULL )
    {
        if ( this->search_fineP_prev_ != NULL )
        {
            delete this->search_fineP_prev_;
            this->search_fineP_prev_ = NULL;
        }
        this->search_fineP_prev_ = new GridGeometricSearch ( this->space_fineP_prev_->meshPtr ( ) );
    }
    if ( this->space_fineP_next_ != NULL )
    {
        if ( this->search_fineP_next_ != NULL )
        {
            delete this->search_fineP_next_;
            this->search_fineP_next_ = NULL;
        }
        this->search_fineP_next_ = new GridGeometricSearch ( this->space_fineP_next_->meshPtr ( ) );
    }
    if ( this->space_fineD_ != NULL )
    {
        if ( this->search_fineD_ != NULL )
        {
            delete this->search_fineD_;
            this->search_fineD_ = NULL;
        }
        this->search_fineD_ = new GridGeometricSearch ( this->space_fineD_->meshPtr ( ) );
    }
    if ( this->space_fineD_next_ != NULL )
    {
        if ( this->search_fineD_next_ != NULL )
        {
            delete this->search_fineD_next_;
            this->search_fineD_next_ = NULL;
        }
        this->search_fineD_next_ = new GridGeometricSearch ( this->space_fineD_next_->meshPtr ( ) );
    }
    if ( this->space_fineD_prev_ != NULL )
    {
        if ( this->search_fineD_prev_ != NULL )
        {
            delete this->search_fineD_prev_;
            this->search_fineD_prev_ = NULL;
        }
        this->search_fineD_prev_ = new GridGeometricSearch ( this->space_fineD_prev_->meshPtr ( ) );
    }
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::setup_time_quadrature ( int order )
{
    this->quad_order_ = order;
    this->quad_c_.clear ( );
    this->quad_w_.clear ( );

    switch ( order )
    {
        case 0:
            this->quad_c_.resize ( 1, 0. );
            this->quad_w_.resize ( 1, 0. );

            this->quad_c_[0] = 1.;
            this->quad_w_[0] = 1.;
            break;

        case 1:
            this->quad_c_.resize ( 1, 0. );
            this->quad_w_.resize ( 1, 0. );

            this->quad_c_[0] = 0.5;
            this->quad_w_[0] = 1.;
            break;
        case 2:
            this->quad_c_.resize ( 2, 0. );
            this->quad_w_.resize ( 2, 0. );

            this->quad_c_[0] = ( -1. / std::sqrt ( 3. ) + 1. ) / 2.;
            this->quad_c_[1] = ( 1. / std::sqrt ( 3. ) + 1. ) / 2.;
            this->quad_w_[0] = 0.5;
            this->quad_w_[1] = 0.5;
            break;
        case 3:
            this->quad_c_.resize ( 3, 0. );
            this->quad_w_.resize ( 3, 0. );

            this->quad_c_[0] = ( -std::sqrt ( 3. / 5. ) + 1. ) / 2.;
            this->quad_c_[1] = 0.5;
            this->quad_c_[2] = ( std::sqrt ( 3. / 5. ) + 1. ) / 2.;

            this->quad_w_[0] = 5. / 18.;
            this->quad_w_[1] = 4. / 9.;
            this->quad_w_[2] = 5. / 18.;
            break;
    }
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::w_i ( int index ) const
{
    assert ( index < this->quad_w_.size ( ) );
    switch ( this->est_rel_time_ )
    {
        case 0:
            return this->quad_w_[index] * this->dT_pc_;
            break;
        case 1:
            return this->quad_w_[index] * this->dT_cn_;
            break;
    }
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::c_i ( int index ) const
{
    assert ( index < this->quad_c_.size ( ) );
    return this->quad_c_[index];
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::get_absolute_time ( DataType c ) const
{
    switch ( this->rel_time_ )
    {
        case 0:
            return this->time_offset_ + this->dT_pc_ * c;
            break;
        case 1:
            return this->time_offset_ + this->dT_pc_ + c * this->dT_cn_;
            break;
    }
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::setup_fe_evaluators ( )
{
    MetFlowAssembler<DIM, DataType>::setup_fe_evaluators ( );
    for ( int v = 0; v<this->num_eq_; ++v )
    {
        if ( this->funP_[v] != NULL )
        {
            delete this->funP_[v];
        }
        if ( this->funP_prev_[v] != NULL )
        {
            delete this->funP_prev_[v];
        }
        if ( this->funP_next_[v] != NULL )
        {
            delete this->funP_next_[v];
        }
        if ( this->funD_[v] != NULL )
        {
            delete this->funD_[v];
        }
        if ( this->funD_prev_[v] != NULL )
        {
            delete this->funD_prev_[v];
        }
        if ( this->funD_next_[v] != NULL )
        {
            delete this->funD_next_[v];
        }

        if ( this->vector_fineP_ != NULL && this->space_fineP_ != NULL )
        {
            this->funP_[v] = new EvalFeFunction<LAD> ( *this->space_fineP_, *this->vector_fineP_, v );
        }
        if ( this->vector_fineP_prev_ != NULL && this->space_fineP_prev_ != NULL )
        {
            this->funP_prev_[v] = new EvalFeFunction<LAD> ( *this->space_fineP_prev_, *this->vector_fineP_prev_, v );
        }
        if ( this->vector_fineP_next_ != NULL && this->space_fineP_next_ != NULL )
        {
            this->funP_next_[v] = new EvalFeFunction<LAD> ( *this->space_fineP_next_, *this->vector_fineP_next_, v );
        }
        if ( this->vector_fineD_ != NULL && this->space_fineD_ != NULL )
        {
            this->funD_[v] = new EvalFeFunction<LAD> ( *this->space_fineD_, *this->vector_fineD_, v );
        }
        if ( this->vector_fineD_next_ != NULL && this->space_fineD_next_ != NULL )
        {
            this->funD_next_[v] = new EvalFeFunction<LAD> ( *this->space_fineD_next_, *this->vector_fineD_next_, v );
        }
        if ( this->vector_fineD_prev_ != NULL && this->space_fineD_prev_ != NULL )
        {
            this->funD_prev_[v] = new EvalFeFunction<LAD> ( *this->space_fineD_prev_, *this->vector_fineD_prev_, v );
        }
    }
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::assemble_local_est_cell ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& vals )
{
    int num_mode = 2;
    int size = this->num_eq_ * num_mode;
    vals.clear ( );
    vals.resize ( size * 4, 0. );

    this->initialize_for_element_estimator ( element, quadrature );
    if ( this->assemble_primal_indicator_ )
    {
        for ( int je = 0; je<this->num_eq_; ++je )
        {
            if ( this->assemble_spatial_indicator_ )
            {
                if ( this->cell_mode_ == RESIDUAL )
                {
                    this->assemble_primal_cell_residual ( element, je, 0, vals[4 * je + 0], vals[4 * je + 1] );
                }
            }
            if ( this->assemble_temporal_indicator_ )
            {
                if ( this->cell_mode_ == RESIDUAL )
                {
                    this->assemble_primal_cell_residual ( element, je, 1, vals[4 * je + 2], vals[4 * je + 3] );
                }
            }
        }

    }
    int offset = this->num_eq_ * 4;
    if ( this->assemble_dual_indicator_ )
    {
        for ( int je = 0; je<this->num_eq_; ++je )
        {
            if ( this->assemble_spatial_indicator_ )
            {
                if ( this->cell_mode_ == RESIDUAL )
                {
                    this->assemble_dual_cell_residual ( element, je, 0, vals[4 * je + offset], vals[4 * je + offset + 1] );
                }
                else if ( this->cell_mode_ == TIME_JUMP )
                {
                    this->assemble_dual_cell_time_jump ( element, je, 0, vals[4 * je + offset], vals[4 * je + offset + 1] );
                }
            }
            if ( this->assemble_temporal_indicator_ )
            {
                if ( this->cell_mode_ == RESIDUAL )
                {
                    this->assemble_dual_cell_residual ( element, je, 1, vals[4 * je + offset + 2], vals[4 * je + offset + 3] );
                }
                else if ( this->cell_mode_ == TIME_JUMP )
                {
                    this->assemble_dual_cell_time_jump ( element, je, 1, vals[4 * je + offset + 2], vals[4 * je + offset + 3] );
                }
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::assemble_local_est_interface ( const Element<DataType>& left_elem,
                                                                              const Element<DataType>& right_elem,
                                                                              const Quadrature<DataType>& left_quad,
                                                                              const Quadrature<DataType>& right_quad,
                                                                              int left_facet_number,
                                                                              int right_facet_number,
                                                                              InterfaceSide left_if_side,
                                                                              InterfaceSide right_if_side,
                                                                              LocalVector& vals )
{
    const bool is_boundary = ( right_if_side == DGGlobalAssembler<DataType>::INTERFACE_BOUNDARY );
    int num_mode = 2;

    int size = this->num_eq_ * num_mode;
    vals.clear ( );
    vals.resize ( size * 4, 0. );

    if ( right_if_side == left_if_side && ( !is_boundary ) )
    {
        return;
    }
    this->initialize_for_interface_estimator ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side );

    if ( this->assemble_primal_indicator_ )
    {
        for ( int je = 0; je<this->num_eq_; ++je )
        {
            if ( this->assemble_spatial_indicator_ )
            {
                if ( !is_boundary )
                {
                    this->assemble_primal_interface_jump ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side,
                                                           je, 0, vals[4 * je + 0], vals[4 * je + 1] );
                }
                else
                {
                    this->assemble_primal_interface_boundary ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side,
                                                               je, 0, vals[4 * je + 0], vals[4 * je + 1] );
                }
            }
            if ( this->assemble_temporal_indicator_ )
            {
                if ( !is_boundary )
                {
                    this->assemble_primal_interface_jump ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side,
                                                           je, 1, vals[4 * je + 2], vals[4 * je + 3] );
                }
                else
                {
                    this->assemble_primal_interface_boundary ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side,
                                                               je, 1, vals[4 * je + 2], vals[4 * je + 3] );
                }
            }
        }

    }
    int offset = this->num_eq_ * 4;
    if ( this->assemble_dual_indicator_ )
    {
        for ( int je = 0; je<this->num_eq_; ++je )
        {
            if ( this->assemble_spatial_indicator_ )
            {
                if ( !is_boundary )
                {
                    this->assemble_dual_interface_jump ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side,
                                                         je, 0, vals[4 * je + 0 + offset], vals[4 * je + 1 + offset] );
                }
                else
                {
                    this->assemble_dual_interface_boundary ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side,
                                                             je, 0, vals[4 * je + 0 + offset], vals[4 * je + 1 + offset] );
                }
            }
            if ( this->assemble_temporal_indicator_ )
            {
                if ( !is_boundary )
                {
                    this->assemble_dual_interface_jump ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side,
                                                         je, 1, vals[4 * je + 2 + offset], vals[4 * je + 3 + offset] );
                }
                else
                {
                    this->assemble_dual_interface_boundary ( left_elem, right_elem, left_quad, right_quad, left_facet_number, right_facet_number, left_if_side, right_if_side,
                                                             je, 1, vals[4 * je + 2 + offset], vals[4 * je + 3 + offset] );
                }
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::initialize_for_element_diameter ( const Element<DataType>& element, const Quadrature<DataType>& quadrature )
{
    AssemblyAssistant<DIM, DataType>::initialize_for_element ( element, quadrature );
    const int num_q = MetFlowAssembler<DIM, DataType>::num_quadrature_points ( );
    DataType vol = 0.;
    for ( int q = 0.; q < num_q; ++q )
    {
        vol += AssemblyAssistant<DIM, DataType>::w ( q ) * std::abs ( AssemblyAssistant<DIM, DataType>::detJ ( q ) );
    }
    this->h_K_ = std::pow ( vol, 1. / DIM );
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::initialize_for_element_estimator ( const Element<DataType>& element, const Quadrature<DataType>& quadrature )
{
    AssemblyAssistant<DIM, DataType>::initialize_for_element ( element, quadrature );
    MetFlowAssembler<DIM, DataType> ::setup_time_fem_evaluator ( );
    MetFlowAssembler<DIM, DataType> ::initialize_for_element_convection ( );
    MetFlowAssembler<DIM, DataType> ::initialize_for_element_source ( );

    const int rel_time = this->est_rel_time_;
    const int num_q = MetFlowAssembler<DIM, DataType>::num_quadrature_points ( );
    DataType vol = 0.;
    for ( int q = 0.; q < num_q; ++q )
    {
        vol += AssemblyAssistant<DIM, DataType>::w ( q ) * std::abs ( AssemblyAssistant<DIM, DataType>::detJ ( q ) );
    }
    this->h_K_ = std::pow ( vol, 1. / DIM );

    // velocity ********************************************
    for ( int v = 0; v<this->num_var_; ++v )
    {
        // primal ------------------------------------------
        // current
        this->solP_ [v].clear ( );
        this->grad_solP_ [v].clear ( );
        this->hess_solP_ [v].clear ( );
        if ( this->vector_solP_ != NULL )
        {
            this->evaluate_fe_function ( this->vector_solP ( ), v, this->solP_ [v] );
            this->evaluate_fe_function_gradients ( this->vector_solP ( ), v, this->grad_solP_ [v] );
            this->evaluate_fe_function_hessians ( this->vector_solP ( ), v, this->hess_solP_ [v] );
        }

        // previuos
        this->solP_prev_ [v].clear ( );
        this->grad_solP_prev_[v].clear ( );
        this->hess_solP_prev_[v].clear ( );
        if ( this->vector_solP_prev_ != NULL )
        {
            this->evaluate_fe_function ( this->vector_solP_prev ( ), v, this->solP_prev_ [v] );
            if ( rel_time == 0 )
            {
                this->evaluate_fe_function_gradients ( this->vector_solP_prev ( ), v, this->grad_solP_prev_[v] );
                this->evaluate_fe_function_hessians ( this->vector_solP_prev ( ), v, this->hess_solP_prev_[v] );
            }
            else
            {
                this->grad_solP_prev_[v].zeros ( num_q );
                this->hess_solP_prev_[v].zeros ( num_q );
            }
        }

        // next
        this->solP_next_ [v].clear ( );
        this->grad_solP_next_[v].clear ( );
        this->hess_solP_next_[v].clear ( );
        if ( this->vector_solP_next_ != NULL )
        {
            this->evaluate_fe_function ( this->vector_solP_next ( ), v, this->solP_next_ [v] );
            if ( rel_time == 1 )
            {
                this->evaluate_fe_function_gradients ( this->vector_solP_next ( ), v, this->grad_solP_next_[v] );
                this->evaluate_fe_function_hessians ( this->vector_solP_next ( ), v, this->hess_solP_next_[v] );
            }
            else
            {
                this->grad_solP_next_[v].zeros ( num_q );
                this->hess_solP_next_[v].zeros ( num_q );
            }
        }

        // dual --------------------------------------------
        // current 
        this->solD_ [v].clear ( );
        this->grad_solD_ [v].clear ( );
        this->hess_solD_ [v].clear ( );
        if ( this->vector_solD_ != NULL )
        {
            this->evaluate_fe_function ( this->vector_solD ( ), v, this->solD_ [v] );
            this->evaluate_fe_function_gradients ( this->vector_solD ( ), v, this->grad_solD_ [v] );
            this->evaluate_fe_function_hessians ( this->vector_solD ( ), v, this->hess_solD_ [v] );
        }

        // previuous
        this->solD_prev_ [v].clear ( );
        this->grad_solD_prev_[v].clear ( );
        this->hess_solD_prev_[v].clear ( );
        if ( this->vector_solD_prev_ != NULL )
        {
            this->evaluate_fe_function ( this->vector_solD_prev ( ), v, this->solD_prev_ [v] );
            if ( rel_time == 0 )
            {
                this->evaluate_fe_function_gradients ( this->vector_solD_prev ( ), v, this->grad_solD_prev_[v] );
                this->evaluate_fe_function_hessians ( this->vector_solD_prev ( ), v, this->hess_solD_prev_[v] );
            }
            else
            {
                this->grad_solD_prev_[v].zeros ( num_q );
                this->hess_solD_prev_[v].zeros ( num_q );
            }
        }

        // next
        this->solD_next_ [v].clear ( );
        this->grad_solD_next_[v].clear ( );
        this->hess_solD_next_[v].clear ( );
        if ( this->vector_solD_next_ != NULL )
        {
            this->evaluate_fe_function ( this->vector_solD_next ( ), v, this->solD_next_ [v] );
            if ( rel_time == 1 )
            {
                this->evaluate_fe_function_gradients ( this->vector_solD_next ( ), v, this->grad_solD_next_[v] );
                this->evaluate_fe_function_hessians ( this->vector_solD_next ( ), v, this->hess_solD_next_[v] );
            }
            else
            {
                this->grad_solD_next_[v].zeros ( num_q );
                this->hess_solD_next_[v].zeros ( num_q );
            }
        }
    }

    // higher order interpolation ---------------------------
    std::vector< Vec<DIM, DataType> > quad_points ( num_q );

    for ( int q = 0; q < num_q; ++q )
    {
        quad_points[q] = MetFlowAssembler<DIM, DataType>::x ( q );
    }

    if ( this->unique_fine_space_ )
    {
        std::vector< SortedArray<int> > trial_cells ( 1 );
        for ( int v = 0; v<this->num_var_; ++v )
        {
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_, *this->funP_[v], trial_cells[0], this->fineP_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_prev_, *this->funP_prev_[v], trial_cells[0], this->fineP_prev_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_next_, *this->funP_next_[v], trial_cells[0], this->fineP_next_[v] );

            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_, *this->funD_[v], trial_cells[0], this->fineD_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_prev_, *this->funD_prev_[v], trial_cells[0], this->fineD_prev_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_next_, *this->funD_next_[v], trial_cells[0], this->fineD_next_[v] );
        }
    }
    else
    {
        std::vector< SortedArray<int> > trial_cells ( 6 );
        for ( int v = 0; v<this->num_var_; ++v )
        {
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_, *this->funP_[v], trial_cells[0], this->fineP_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_prev_, *this->funP_prev_[v], trial_cells[1], this->fineP_prev_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_next_, *this->funP_next_[v], trial_cells[2], this->fineP_next_[v] );

            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_, *this->funD_[v], trial_cells[3], this->fineD_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_prev_, *this->funD_prev_[v], trial_cells[4], this->fineD_prev_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_next_, *this->funD_next_[v], trial_cells[5], this->fineD_next_[v] );
        }
    }
}

template<int DIM, class DataType>
void MetFlowEstimatorAssembler<DIM, DataType>::initialize_for_interface_estimator ( const Element<DataType>& left_elem,
                                                                                    const Element<DataType>& right_elem,
                                                                                    const Quadrature<DataType>& left_quad,
                                                                                    const Quadrature<DataType>& right_quad,
                                                                                    int left_facet_number,
                                                                                    int right_facet_number,
                                                                                    InterfaceSide left_if_side,
                                                                                    InterfaceSide right_if_side )
{
    const bool is_boundary = ( right_if_side == DGGlobalAssembler<DataType>::INTERFACE_BOUNDARY );

    DGAssemblyAssistant<DIM, DataType>::initialize_for_interface ( left_elem, right_elem,
                                                                   left_quad, right_quad,
                                                                   left_facet_number, right_facet_number,
                                                                   left_if_side, right_if_side );

    MetFlowAssembler<DIM, DataType> ::setup_time_fem_evaluator ( );

    const int num_q = DGAssemblyAssistant<DIM, DataType>::num_quadrature_points ( );
    DataType vol = 0.;
    for ( int q = 0.; q < num_q; ++q )
    {
        vol += DGAssemblyAssistant<DIM, DataType>::w ( q ) * std::abs ( DGAssemblyAssistant<DIM, DataType>::ds ( q ) );
    }
    this->h_E_ = std::pow ( vol, 1. / ( DIM - 1 ) );

    // primal stuff
    for ( int v = 0; v<this->num_var_; ++v )
    {
        this->left_grad_solP_[v].clear ( );
        this->left_grad_solP_prev_[v].clear ( );
        this->left_grad_solP_next_[v].clear ( );

        this->right_grad_solP_[v].clear ( );
        this->right_grad_solP_prev_[v].clear ( );
        this->right_grad_solP_next_[v].clear ( );

        this->solP_ [v].clear ( );
        this->solP_prev_[v].clear ( );
        this->solP_next_[v].clear ( );

        if ( !is_boundary )
        {
            if ( this->vector_solP_ != NULL )
            {
                this->trial ( ).evaluate_fe_function_gradients ( this->vector_solP ( ), v, this->left_grad_solP_[v] );
            }
            if ( this->vector_solP_prev_ != NULL )
            {
                this->trial ( ).evaluate_fe_function_gradients ( this->vector_solP_prev ( ), v, this->left_grad_solP_prev_[v] );
            }
            if ( this->vector_solP_next_ != NULL )
            {
                this->trial ( ).evaluate_fe_function_gradients ( this->vector_solP_next ( ), v, this->left_grad_solP_next_[v] );
            }
        }

        if ( this->vector_solP_ != NULL )
        {
            this->test ( ).evaluate_fe_function_gradients ( this->vector_solP ( ), v, this->right_grad_solP_[v] );
            this->test ( ).evaluate_fe_function ( this->vector_solP ( ), v, this->solP_[v] );
        }
        if ( this->vector_solP_prev_ != NULL )
        {
            this->test ( ).evaluate_fe_function_gradients ( this->vector_solP_prev ( ), v, this->right_grad_solP_prev_[v] );
            this->test ( ).evaluate_fe_function ( this->vector_solP_prev ( ), v, this->solP_prev_[v] );
        }
        if ( this->vector_solP_next_ != NULL )
        {
            this->test ( ).evaluate_fe_function_gradients ( this->vector_solP_next ( ), v, this->right_grad_solP_next_[v] );
            this->test ( ).evaluate_fe_function ( this->vector_solP_next ( ), v, this->solP_next_[v] );
        }

        int rank;
        MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
        /*
        if (!is_boundary)
        {
            for (int q=0; q<num_q; ++q)
            {
                DataType co_x  = DGAssemblyAssistant<DIM,DataType>::x(q)[0];
                DataType co_y  = DGAssemblyAssistant<DIM,DataType>::x(q)[1];
                DataType left  = this->left_grad_solP_[v][q][0];
                DataType right = this->right_grad_solP_[v][q][0];

                DataType rel_diff = std::abs(left - right) / std::abs(right);
                if (co_x == 0.5 && co_y <= 0.5 && (rank == 1 || rank == 0))
                {
                    std::cout << " ------------------------ " << std::endl;
                    std::cout << rank << "|| " << co_x << ", " << co_y << ": " << left << " <-> " << right << " " << rel_diff << std::endl;
                    std::cout << " ------------------------ " << std::endl;
                }
                else
                {
                   // std::cout << co_x << ", " << co_y << ": " << left << " <-> " << right << " " << rel_diff << std::endl;
                }
            }
        }
         */
        // dual stuff
        this->left_grad_solD_[v].clear ( );
        this->left_grad_solD_next_[v].clear ( );
        this->left_grad_solD_prev_[v].clear ( );

        this->right_grad_solD_[v].clear ( );
        this->right_grad_solD_next_[v].clear ( );
        this->right_grad_solD_prev_[v].clear ( );

        this->solD_ [v].clear ( );
        this->solD_prev_[v].clear ( );
        this->solD_next_[v].clear ( );

        if ( !is_boundary )
        {
            if ( this->vector_solD_ != NULL )
            {
                this->trial ( ).evaluate_fe_function_gradients ( this->vector_solD ( ), v, this->left_grad_solD_[v] );
            }
            if ( this->vector_solD_next_ != NULL )
            {
                this->trial ( ).evaluate_fe_function_gradients ( this->vector_solD_next ( ), v, this->left_grad_solD_next_[v] );
            }
            if ( this->vector_solD_prev_ != NULL )
            {
                this->trial ( ).evaluate_fe_function_gradients ( this->vector_solD_prev ( ), v, this->left_grad_solD_prev_[v] );
            }
        }

        if ( this->vector_solD_ != NULL )
        {
            this->test ( ).evaluate_fe_function_gradients ( this->vector_solD ( ), v, this->right_grad_solD_[v] );
            this->test ( ).evaluate_fe_function ( this->vector_solD ( ), v, this->solD_[v] );
        }
        if ( this->vector_solD_next_ != NULL )
        {
            this->test ( ).evaluate_fe_function_gradients ( this->vector_solD_next ( ), v, this->right_grad_solD_next_[v] );
            this->test ( ).evaluate_fe_function ( this->vector_solD_next ( ), v, this->solD_next_[v] );
        }
        if ( this->vector_solD_prev_ != NULL )
        {
            this->test ( ).evaluate_fe_function_gradients ( this->vector_solD_prev ( ), v, this->right_grad_solD_prev_[v] );
            this->test ( ).evaluate_fe_function ( this->vector_solD_prev ( ), v, this->solD_prev_[v] );
        }
    }

    // higher order interpolation ---------------------------
    std::vector< Vec<DIM, DataType> > quad_points ( num_q );

    for ( int q = 0; q < num_q; ++q )
    {
        quad_points[q] = DGAssemblyAssistant<DIM, DataType>::x ( q );
    }

    if ( this->unique_fine_space_ )
    {
        std::vector< SortedArray<int> > trial_cells ( 1 );
        for ( int v = 0; v<this->num_var_; ++v )
        {
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_, *this->funP_[v], trial_cells[0], this->fineP_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_prev_, *this->funP_prev_[v], trial_cells[0], this->fineP_prev_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_next_, *this->funP_next_[v], trial_cells[0], this->fineP_next_[v] );

            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_, *this->funD_[v], trial_cells[0], this->fineD_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_prev_, *this->funD_prev_[v], trial_cells[0], this->fineD_prev_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_next_, *this->funD_next_[v], trial_cells[0], this->fineD_next_[v] );
        }
    }
    else
    {
        std::vector< SortedArray<int> > trial_cells ( 6 );
        for ( int v = 0; v<this->num_var_; ++v )
        {
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_, *this->funP_[v], trial_cells[0], this->fineP_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_prev_, *this->funP_prev_[v], trial_cells[1], this->fineP_prev_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineP_next_, *this->funP_next_[v], trial_cells[2], this->fineP_next_[v] );

            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_, *this->funD_[v], trial_cells[3], this->fineD_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_prev_, *this->funD_prev_[v], trial_cells[4], this->fineD_prev_[v] );
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_fineD_next_, *this->funD_next_[v], trial_cells[5], this->fineD_next_[v] );
        }
    }
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::weightP ( DataType c, int q, int var, int est_mode ) const
{
    if ( !use_dwr_ )
    {
        return 1.;
    }
    switch ( est_mode )
    {
        case 0:
            return this->weightP_h ( c, q, var );
            break;
        case 1:
            return this->weightP_tau ( c, q, var );
            break;
    }
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::weightD ( DataType c, int q, int var, int est_mode ) const
{
    if ( !use_dwr_ )
    {
        return 1.;
    }
    switch ( est_mode )
    {
        case 0:
            return this->weightD_h ( c, q, var );
            break;
        case 1:
            return this->weightD_tau ( c, q, var );
            break;
    }
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::weightP_tau ( DataType c, int q, int var ) const
{
    switch ( this->time_order_P_[var] )
    {
        case 0:
            return this->scalar_eval_.patch_linear ( c, 0., this->solP_[var][q], this->solP_next_[var][q] )
                    - this->scalar_eval_.constant ( c, 0., this->solP_[var][q], this->solP_next_[var][q] );
            break;
        default:
            return this->scalar_eval_.patch_quadratic ( c, this->solP_prev_[var][q], this->solP_[var][q], this->solP_next_[var][q] )
                    - this->scalar_eval_.linear ( c, this->solP_prev_[var][q], this->solP_[var][q], this->solP_next_[var][q] );
            break;
    }
    return 9e6;
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::weightD_tau ( DataType c, int q, int var ) const
{
    switch ( this->time_order_D_[var] )
    {
        case 0:
            return this->scalar_eval_.patch_linear ( c, 0., this->solD_[var][q], this->solD_next_[var][q] )
                    - this->scalar_eval_.constant ( c, 0., this->solD_[var][q], this->solD_next_[var][q] );
            break;
        default:
            return this->scalar_eval_.patch_quadratic ( c, this->solD_prev_[var][q], this->solD_[var][q], this->solD_next_[var][q] )
                    - this->scalar_eval_.linear ( c, this->solD_prev_[var][q], this->solD_[var][q], this->solD_next_[var][q] );
            break;
    }
    return 9e6;
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::weightP_h ( DataType c, int q, int var ) const
{
    switch ( this->time_order_P_[var] )
    {
        case 0:
            return this->scalar_eval_.constant ( c, 0., this->solP_[var][q], this->solP_next_[var][q] )
                    - this->scalar_eval_.constant ( c, 0., this->fineP_[var][q], this->fineP_next_[var][q] );
            break;
        default:
            return this->scalar_eval_.linear ( c, this->solP_prev_[var][q], this->solP_[var][q], this->solP_next_[var][q] )
                    - this->scalar_eval_.linear ( c, this->fineP_prev_[var][q], this->fineP_[var][q], this->fineP_next_[var][q] );
            break;
    }
    return 9e6;
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::weightD_h ( DataType c, int q, int var ) const
{
    switch ( this->time_order_D_[var] )
    {
        case 0:
            return this->scalar_eval_.constant ( c, 0., this->solD_[var][q], this->solD_next_[var][q] )
                    - this->scalar_eval_.constant ( c, 0., this->fineD_[var][q], this->fineD_next_[var][q] );
            break;
        default:
            return this->scalar_eval_.linear ( c, this->solD_prev_[var][q], this->solD_[var][q], this->solD_next_[var][q] )
                    - this->scalar_eval_.linear ( c, this->fineD_prev_[var][q], this->fineD_[var][q], this->fineD_next_[var][q] );
            break;
    }
    return 9e6;
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::trialP ( DataType c, int q, int var ) const
{
    switch ( this->time_order_P_[var] )
    {
        case 0:
            return this->scalar_eval_.constant ( c, 0., this->solP_[var][q], this->solP_next_[var][q] );
            break;
        default:
            return this->scalar_eval_.linear ( c, this->solP_prev_[var][q], this->solP_[var][q], this->solP_next_[var][q] );
            break;
    }
    return 9e6;
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::dt_trialP ( DataType c, int q, int var ) const
{
    switch ( this->time_order_P_[var] )
    {
        case 0:
            return 0.;
            break;
        default:
            return this->scalar_eval_.constant ( c, 0., this->solP_[var][q], this->solP_next_[var][q] );
            break;
    }
    return 9e6;
}

template<int DIM, class DataType>
Vec<DIM, DataType> MetFlowEstimatorAssembler<DIM, DataType>::grad_trialP ( DataType c, int q, int var ) const
{
    switch ( this->time_order_P_[var] )
    {
        case 0:
            return this->vec_eval_.constant ( c, Vec<DIM, DataType> ( ), this->grad_solP_[var][q], this->grad_solP_next_[var][q] );
            break;
        default:
            return this->vec_eval_.linear ( c, this->grad_solP_prev_[var][q], this->grad_solP_[var][q], this->grad_solP_next_[var][q] );
            break;
    }
    return Vec<DIM, DataType> ( );
}

template<int DIM, class DataType>
Vec<DIM, DataType> MetFlowEstimatorAssembler<DIM, DataType>::left_grad_trialP ( DataType c, int q, int var ) const
{
    switch ( this->time_order_P_[var] )
    {
        case 0:
            return this->vec_eval_.constant ( c, Vec<DIM, DataType> ( ), this->left_grad_solP_[var][q], this->left_grad_solP_next_[var][q] );
            break;
        default:
            return this->vec_eval_.linear ( c, this->left_grad_solP_prev_[var][q], this->left_grad_solP_[var][q], this->left_grad_solP_next_[var][q] );
            break;
    }
    return Vec<DIM, DataType> ( );
}

template<int DIM, class DataType>
Vec<DIM, DataType> MetFlowEstimatorAssembler<DIM, DataType>::right_grad_trialP ( DataType c, int q, int var ) const
{
    switch ( this->time_order_P_[var] )
    {
        case 0:
            return this->vec_eval_.constant ( c, Vec<DIM, DataType> ( ), this->right_grad_solP_[var][q], this->right_grad_solP_next_[var][q] );
            break;
        default:
            return this->vec_eval_.linear ( c, this->right_grad_solP_prev_[var][q], this->right_grad_solP_[var][q], this->right_grad_solP_next_[var][q] );
            break;
    }
    return Vec<DIM, DataType> ( );
}

template<int DIM, class DataType>
Mat<DIM, DIM, DataType> MetFlowEstimatorAssembler<DIM, DataType>::hess_trialP ( DataType c, int q, int var ) const
{
    switch ( this->time_order_P_[var] )
    {
        case 0:
            return this->mat_eval_.constant ( c, Mat<DIM, DIM, DataType> ( ), this->hess_solP_[var][q], this->hess_solP_next_[var][q] );
            break;
        default:
            return this->mat_eval_.linear ( c, this->hess_solP_prev_[var][q], this->hess_solP_[var][q], this->hess_solP_next_[var][q] );
            break;
    }
    return Mat<DIM, DIM, DataType> ( );
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::trialD ( DataType c, int q, int var ) const
{
    switch ( this->time_order_D_[var] )
    {
        case 0:
            return this->scalar_eval_.constant ( c, 0., this->solD_[var][q], this->solD_next_[var][q] );
            break;
        default:
            return this->scalar_eval_.linear ( c, this->solD_prev_[var][q], this->solD_[var][q], this->solD_next_[var][q] );
            break;
    }
    return 9e6;
}

template<int DIM, class DataType>
DataType MetFlowEstimatorAssembler<DIM, DataType>::dt_trialD ( DataType c, int q, int var ) const
{
    switch ( this->time_order_D_[var] )
    {
        case 0:
            return 0.;
            break;
        default:
            return this->scalar_eval_.constant ( c, 0., this->solD_[var][q], this->solD_next_[var][q] );
            break;
    }
    return 9e6;
}

template<int DIM, class DataType>
Vec<DIM, DataType> MetFlowEstimatorAssembler<DIM, DataType>::grad_trialD ( DataType c, int q, int var ) const
{
    switch ( this->time_order_D_[var] )
    {
        case 0:
            return this->vec_eval_.constant ( c, Vec<DIM, DataType> ( ), this->grad_solD_[var][q], this->grad_solD_next_[var][q] );
            break;
        default:
            return this->vec_eval_.linear ( c, this->grad_solD_prev_[var][q], this->grad_solD_[var][q], this->grad_solD_next_[var][q] );
            break;
    }
    return Vec<DIM, DataType>( );
}

template<int DIM, class DataType>
Vec<DIM, DataType> MetFlowEstimatorAssembler<DIM, DataType>::left_grad_trialD ( DataType c, int q, int var ) const
{
    switch ( this->time_order_D_[var] )
    {
        case 0:
            return this->vec_eval_.constant ( c, Vec<DIM, DataType> ( ), this->left_grad_solD_[var][q], this->left_grad_solD_next_[var][q] );
            break;
        default:
            return this->vec_eval_.linear ( c, this->left_grad_solD_prev_[var][q], this->left_grad_solD_[var][q], this->left_grad_solD_next_[var][q] );
            break;
    }
    return Vec<DIM, DataType>( );
}

template<int DIM, class DataType>
Vec<DIM, DataType> MetFlowEstimatorAssembler<DIM, DataType>::right_grad_trialD ( DataType c, int q, int var ) const
{
    switch ( this->time_order_D_[var] )
    {
        case 0:
            return this->vec_eval_.constant ( c, Vec<DIM, DataType> ( ), this->right_grad_solD_[var][q], this->right_grad_solD_next_[var][q] );
            break;
        default:
            return this->vec_eval_.linear ( c, this->right_grad_solD_prev_[var][q], this->right_grad_solD_[var][q], this->right_grad_solD_next_[var][q] );
            break;
    }
    return Vec<DIM, DataType>( );
}

template<int DIM, class DataType>
Mat<DIM, DIM, DataType> MetFlowEstimatorAssembler<DIM, DataType>::hess_trialD ( DataType c, int q, int var ) const
{
    switch ( this->time_order_D_[var] )
    {
        case 0:
            return this->mat_eval_.constant ( c, Mat<DIM, DIM, DataType> ( ), this->hess_solD_[var][q], this->hess_solD_next_[var][q] );
            break;
        default:
            return this->mat_eval_.linear ( c, this->hess_solD_prev_[var][q], this->hess_solD_[var][q], this->hess_solD_next_[var][q] );
            break;
    }
    return Mat<DIM, DIM, DataType>( );
}

template class MetFlowEstimatorAssembler<2, double>;
template class MetFlowEstimatorAssembler<3, double>;

