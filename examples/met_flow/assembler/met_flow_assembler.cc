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

#include "met_flow_assembler.h"

template<int DIM, class DataType>
MetFlowAssembler<DIM, DataType>::MetFlowAssembler ( )
: vector_base_ ( NULL ),
vector_perturb_ ( NULL ),
vector_perturb_prev_ ( NULL ),
vector_solD_ ( NULL ),
vector_solD_next_ ( NULL ),
vector_solD_prev_ ( NULL ),
vector_solP_ ( NULL ),
vector_solP_next_ ( NULL ),
vector_solP_prev_ ( NULL ),
vector_conv_ ( NULL ),
vector_conv_prev_ ( NULL ),
vector_conv_next_ ( NULL ),
source_term_ ( NULL ),
convection_term_ ( NULL ),
search_conv_ ( NULL ),
space_conv_ ( NULL ),
final_goal_contrib_ ( false ),
num_var_ ( 0 )
{
    this->is_matrix_constant_ = false;
    this->perturb_scale_ = 1.;
    this->t_ = -999999999.;
    this->dT_pc_ = -999999999.;
    this->dT_cn_ = -999999999.;
    this->mode_ = PRIMAL;
    this->L2_perturb_.resize ( DIM + 1, false );
    this->H1_perturb_.resize ( DIM + 1, false );
    this->base_vector_set_ = false;
    this->set_time_discretization_to_zero ( );

    this->funConv_.resize ( DIM, NULL );
    this->funConv_prev_.resize ( DIM, NULL );
    this->funConv_next_.resize ( DIM, NULL );

    this->funGradConv_.resize ( DIM );
    this->funGradConv_prev_.resize ( DIM );
    this->funGradConv_next_.resize ( DIM );

    for ( int v = 0; v < DIM; ++v )
    {
        this->funGradConv_[v].resize ( DIM, NULL );
        this->funGradConv_prev_[v].resize ( DIM, NULL );
        this->funGradConv_next_[v].resize ( DIM, NULL );
    }
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::clear ( )
{
    this->is_matrix_constant_ = false;
    this->perturb_scale_ = 1.;
    this->t_ = -999999999.;
    this->dT_pc_ = -999999999.;
    this->dT_cn_ = -999999999.;
    this->mode_ = PRIMAL;
    this->L2_perturb_.clear ( );
    this->L2_perturb_.resize ( DIM + 1, false );
    this->H1_perturb_.clear ( );
    this->H1_perturb_.resize ( DIM + 1, false );
    this->base_vector_set_ = false;
    this->set_time_discretization_to_zero ( );

    // TODO delete verwenden
    this->funConv_.clear ( );
    this->funConv_.resize ( DIM, NULL );
    this->funConv_prev_.clear ( );
    this->funConv_prev_.resize ( DIM, NULL );
    this->funConv_next_.clear ( );
    this->funConv_next_.resize ( DIM, NULL );

    this->funGradConv_.clear ( );
    this->funGradConv_prev_.clear ( );
    this->funGradConv_next_.clear ( );
    this->funGradConv_.resize ( DIM );
    this->funGradConv_prev_.resize ( DIM );
    this->funGradConv_next_.resize ( DIM );

    for ( int v = 0; v < DIM; ++v )
    {
        this->funGradConv_[v].resize ( DIM, NULL );
        this->funGradConv_prev_[v].resize ( DIM, NULL );
        this->funGradConv_next_[v].resize ( DIM, NULL );
    }

    this->galerkin_mode_ = 0;
    this->rel_time_ = 0;

    this->set_time_discretization_to_zero ( );

    this->final_goal_contrib_ = false;

    this->scalar_eval_.clear ( );
    this->vec_eval_.clear ( );
    this->mat_eval_.clear ( );

    this->goal_functional_ = NULL;

    /*
    if (search_conv_ != NULL)
    {
        delete search_conv_;
    }
     * */
    this->search_conv_ = NULL;
    this->space_conv_ = NULL;
    this->source_term_ = NULL;
    this->convection_term_ = NULL;

    this->vector_solP_ = NULL;
    this->vector_solP_prev_ = NULL;
    this->vector_solP_next_ = NULL;
    this->vector_solD_ = NULL;
    this->vector_solD_next_ = NULL;
    this->vector_solD_prev_ = NULL;
    this->vector_perturb_ = NULL;
    this->vector_perturb_prev_ = NULL;
    this->vector_base_ = NULL;
    this->vector_conv_next_ = NULL;
    this->vector_conv_ = NULL;
    this->vector_conv_prev_ = NULL;

    this->solP_.clear ( );
    this->solP_prev_.clear ( );
    this->solP_next_.clear ( );
    this->grad_solP_.clear ( );
    this->grad_solP_prev_.clear ( );
    this->grad_solP_next_.clear ( );
    this->hess_solP_.clear ( );
    this->hess_solP_prev_.clear ( );
    this->hess_solP_next_.clear ( );
    this->solD_.clear ( );
    this->solD_prev_.clear ( );
    this->solD_next_.clear ( );
    this->grad_solD_.clear ( );
    this->grad_solD_prev_.clear ( );
    this->grad_solD_next_.clear ( );
    this->hess_solD_.clear ( );
    this->hess_solD_prev_.clear ( );
    this->hess_solD_next_.clear ( );
    this->perturb_.clear ( );
    this->perturb_prev_.clear ( );
    this->grad_perturb_.clear ( );
    this->grad_perturb_prev_.clear ( );
    this->base_.clear ( );

    // Function values
    for ( int d = 0; d < DIM; ++d )
    {
        this->conv_[d].clear ( );
        this->conv_prev_[d].clear ( );
        this->conv_next_[d].clear ( );
        this->grad_conv_[d].clear ( );
        this->grad_conv_prev_[d].clear ( );
        this->grad_conv_next_[d].clear ( );
    }
    this->source_.clear ( );
    this->source_prev_.clear ( );
    this->source_next_.clear ( );
}

template<int DIM, class DataType>
MetFlowAssembler<DIM, DataType>::~MetFlowAssembler ( )
{
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::allocate_function_values ( int num_var )
{
    solP_.resize ( num_var );
    solP_prev_.resize ( num_var );
    solP_next_.resize ( num_var );

    grad_solP_.resize ( num_var );
    grad_solP_prev_.resize ( num_var );
    grad_solP_next_.resize ( num_var );

    hess_solP_.resize ( num_var );
    hess_solP_prev_.resize ( num_var );
    hess_solP_next_.resize ( num_var );

    solD_.resize ( num_var );
    solD_prev_.resize ( num_var );
    solD_next_.resize ( num_var );

    grad_solD_.resize ( num_var );
    grad_solD_prev_.resize ( num_var );
    grad_solD_next_.resize ( num_var );

    hess_solD_.resize ( num_var );
    hess_solD_prev_.resize ( num_var );
    hess_solD_next_.resize ( num_var );

    perturb_.resize ( num_var );
    perturb_prev_.resize ( num_var );

    grad_perturb_.resize ( num_var );
    grad_perturb_prev_.resize ( num_var );
    base_.resize ( num_var );
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::set_perturb_type ( int var, bool L2, bool H1 )
{
    this->L2_perturb_[var] = L2;
    this->H1_perturb_[var] = H1;
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::print_perturb ( )
{
    std::cout << "  Perturbation" << std::endl;
    for ( int l = 0; l<this->L2_perturb_.size ( ); ++l )
    {
        std::cout << "  L2 (" << l << ") : " << this->L2_perturb_[l] << std::endl;
        std::cout << "  H1 (" << l << ") : " << this->H1_perturb_[l] << std::endl;
    }
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::set_time_discretization_to_zero ( )
{
    this->theta_d_dt_u_ = 0.;
    this->delta_d_dt_u_ = 0.;
    this->delta_j_c_c_ = 0.;
    this->delta_j_c_pc_ = 0.;
    this->delta_j_n_c_ = 0.;
    this->delta_j_n_nc_ = 0.;
    this->delta_j_n_n_ = 0.;
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::set_time_discretization_to_simple ( )
{
    this->delta_d_dt_u_ = 1.;

    this->delta_j_c_c_ = 0.;
    this->delta_j_c_pc_ = 0.;
    this->delta_j_n_nc_ = 0.;
    this->delta_j_n_c_ = 0.5;
    this->delta_j_n_n_ = 0.5;
    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::set_time_discretization_to_theta ( DataType theta )
{
    if ( this->mode_ == PRIMAL )
    {
        this->theta_d_dt_u_ = 1.;
    }
    if ( this->mode_ == DUAL )
    {
        this->delta_d_dt_u_ = 1.;
        this->delta_j_c_c_ = 0.;
        this->delta_j_c_pc_ = 0.;
        this->delta_j_n_nc_ = 0.;
        this->delta_j_n_c_ = 1. - theta;
        this->delta_j_n_n_ = theta;
    }

    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::set_time_discretization_to_galerkin_cd ( bool modified )
{
    this->galerkin_mode_ = 0;

    // Primal Mode
    if ( this->mode_ == PRIMAL )
    {
        this->theta_d_dt_u_ = 1.;
    }

    // Dual Mode
    if ( this->mode_ == DUAL )
    {
        this->delta_d_dt_u_ = 1.;
        this->delta_j_c_c_ = 0.;
        this->delta_j_c_pc_ = 0.;
        this->delta_j_n_c_ = 1. / 6.;
        this->delta_j_n_nc_ = 4. / 6.;
        this->delta_j_n_n_ = 1. / 6.;
    }
    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::set_time_discretization_to_galerkin_dc ( bool modified )
{
    this->galerkin_mode_ = 1;

    // Dual Mode
    if ( this->mode_ == DUAL )
    {
        this->delta_d_dt_u_ = 1.;
        this->delta_j_c_c_ = 1. / 6.;
        this->delta_j_c_pc_ = 1. / 3.;
        this->delta_j_n_c_ = 1. / 6.;
        this->delta_j_n_nc_ = 1. / 3.;
        this->delta_j_n_n_ = 0.;
    }
    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::set_time_discretization_off ( )
{
    this->dT_pc_ = 1.;
    this->dT_cn_ = 1.;

    if ( this->mode_ == PRIMAL )
    {
        this->theta_d_dt_u_ = 0.;
    }
    if ( this->mode_ == DUAL )
    {
        this->delta_d_dt_u_ = 0.;
        this->delta_j_c_c_ = 1.;
        this->delta_j_c_pc_ = 0.;
        this->delta_j_n_nc_ = 0.;
        this->delta_j_n_c_ = 0.;
        this->delta_j_n_n_ = 0.;
    }

    this->is_matrix_constant_ = false;
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::print_parameters ( )
{
    std::cout.scientific;
    std::cout.precision ( 6 );
    if ( this->mode_ == PRIMAL )
    {
        std::cout << "> MetFlowAssembler primal parameters" << std::endl;
    }
    else
    {
        std::cout << "> MetFlowAssembler dual parameters" << std::endl;
    }

    std::cout << "  dT_pc:    " << this->dT_pc_ << std::endl;
    std::cout << "  dT_cn:    " << this->dT_cn_ << std::endl;
    std::cout << "  Coord:    cylindrical" << std::endl;
#ifdef ROTATING_FOR
    std::cout << "  FoR:      rotating" << std::endl;
#else
    std::cout << "  FoR:      inertial" << std::endl;
#endif
    std::cout << std::endl;

    MetFlowAssembler<DIM, DataType>::print_coeff ( );
    this->print_perturb ( );
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::print_coeff ( )
{
    if ( this->mode_ == PRIMAL )
    {
        std::cout << "  theta_d_dt_u:       " << this->theta_d_dt_u_ << std::endl;
    }
    else
    {
        std::cout << "  delta_d_dt_u:       " << this->delta_d_dt_u_ << std::endl;
        std::cout << "  delta_j_c_c:        " << this->delta_j_c_c_ << std::endl;
        std::cout << "  delta_j_c_pc:       " << this->delta_j_c_pc_ << std::endl;
        std::cout << "  delta_j_n_c:        " << this->delta_j_n_c_ << std::endl;
        std::cout << "  delta_j_n_nc:       " << this->delta_j_n_nc_ << std::endl;
        std::cout << "  delta_j_n_n:        " << this->delta_j_n_n_ << std::endl;
    }
}

template<int DIM, class DataType>
DataType MetFlowAssembler<DIM, DataType>::get_coeff ( std::string term )
{
    return 0.;
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::setup_grid_search ( )
{
    if ( this->space_conv_ != NULL )
    {
        if ( this->search_conv_ != NULL )
        {
            delete this->search_conv_;
            this->search_conv_ = NULL;
        }
        this->search_conv_ = new GridGeometricSearch ( this->space_conv_->meshPtr ( ) );
    }
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::setup_fe_evaluators ( )
{
    if ( this->vector_conv_ == NULL )
    {
        return;
    }

    for ( int v = 0; v < DIM; ++v )
    {
        if ( this->funConv_[v] != NULL )
        {
            delete this->funConv_[v];
        }
        if ( this->funConv_prev_[v] != NULL )
        {
            delete this->funConv_prev_[v];
        }
        if ( this->funConv_next_[v] != NULL )
        {
            delete this->funConv_next_[v];
        }

        if ( this->vector_conv_ != NULL && this->space_conv_ != NULL )
        {
            this->funConv_[v] = new EvalFeFunction<LAD> ( *this->space_conv_, *this->vector_conv_, v );
        }
        if ( this->vector_conv_prev_ != NULL && this->space_conv_ != NULL )
        {
            this->funConv_prev_[v] = new EvalFeFunction<LAD> ( *this->space_conv_, *this->vector_conv_prev_, v );
        }
        if ( this->vector_conv_next_ != NULL && this->space_conv_ != NULL )
        {
            this->funConv_next_[v] = new EvalFeFunction<LAD> ( *this->space_conv_, *this->vector_conv_next_, v );
        }
    }
    for ( int v = 0; v < DIM; ++v )
    {
        for ( int d = 0; d < DIM; ++d )
        {
            if ( this->funGradConv_[v][d] != NULL )
            {
                delete this->funGradConv_[v][d];
            }
            if ( this->funGradConv_prev_[v][d] != NULL )
            {
                delete this->funGradConv_prev_[v][d];
            }
            if ( this->funGradConv_next_[v][d] != NULL )
            {
                delete this->funGradConv_next_[v][d];
            }

            if ( this->vector_conv_ != NULL && this->space_conv_ != NULL )
            {
                this->funGradConv_[v][d] = new EvalDerivativeFeFunction<LAD, DIM> ( *this->space_conv_, *this->vector_conv_, v, d );
            }
            if ( this->vector_conv_prev_ != NULL && this->space_conv_ != NULL )
            {
                this->funGradConv_prev_[v][d] = new EvalDerivativeFeFunction<LAD, DIM> ( *this->space_conv_, *this->vector_conv_prev_, v, d );
            }
            if ( this->vector_conv_next_ != NULL && this->space_conv_ != NULL )
            {
                this->funGradConv_next_[v][d] = new EvalDerivativeFeFunction<LAD, DIM> ( *this->space_conv_, *this->vector_conv_next_, v, d );
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::setup_time_fem_evaluator ( )
{
    this->scalar_eval_.set_time_steps ( this->dT_pc_, this->dT_cn_, this->rel_time_ );
    this->vec_eval_ .set_time_steps ( this->dT_pc_, this->dT_cn_, this->rel_time_ );
    this->mat_eval_ .set_time_steps ( this->dT_pc_, this->dT_cn_, this->rel_time_ );

    this->scalar_eval_.compute_weight_tau_coeff ( );
    this->vec_eval_ .compute_weight_tau_coeff ( );
    this->mat_eval_ .compute_weight_tau_coeff ( );
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::set_newton_solution ( const LAD::VectorType& newton_sol )
{
    if ( this->mode_ == PRIMAL )
    {
        this->vector_solP_ = &newton_sol;
    }
    else
    {
        this->vector_solD_ = &newton_sol;
    }
}

template<int DIM, class DataType>
DataType MetFlowAssembler<DIM, DataType>::source ( DataType c, int q, int var ) const
{
    return this->scalar_eval_.linear ( c, this->source_prev_[var][q], this->source_[var][q], this->source_next_[var][q] );
}

template<int DIM, class DataType>
DataType MetFlowAssembler<DIM, DataType>::convection ( DataType c, int q, int var ) const
{
    return this->scalar_eval_.linear ( c, this->conv_prev_[var][q], this->conv_[var][q], this->conv_next_[var][q] );
}

template<int DIM, class DataType>
Vec<DIM, DataType> MetFlowAssembler<DIM, DataType>::grad_convection ( DataType c, int q, int var ) const
{
    return this->vec_eval_.linear ( c, this->grad_conv_prev_[var][q], this->grad_conv_[var][q], this->grad_conv_next_[var][q] );
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::initialize_for_facet ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, int facet_number )
{
    AssemblyAssistant<DIM, DataType>::initialize_for_facet ( element, quadrature, facet_number );
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::initialize_for_element ( const Element<DataType>& element, const Quadrature<DataType>& quadrature )
{
    AssemblyAssistant<DIM, DataType>::initialize_for_element ( element, quadrature );
    MetFlowAssembler<DIM, DataType> ::initialize_for_element_convection ( );
    MetFlowAssembler<DIM, DataType> ::initialize_for_element_source ( );
    MetFlowAssembler<DIM, DataType> ::setup_time_fem_evaluator ( );

    // *************************************************************************
    // PRIMAL VARIABLES 
    assert ( this->vector_solP_ != NULL );
    assert ( this->vector_solP_prev_ != NULL );

    // velocity
    for ( int v = 0; v<this->num_var_; ++v )
    {
        solP_ [v].clear ( );
        solP_prev_ [v].clear ( );
        grad_solP_ [v].clear ( );
        grad_solP_prev_[v].clear ( );
        if ( this->vector_solP_ != NULL )
        {
            this->evaluate_fe_function ( this->vector_solP ( ), v, solP_ [v] );
            this->evaluate_fe_function_gradients ( this->vector_solP ( ), v, grad_solP_ [v] );
        }
        if ( this->vector_solP_prev_ != NULL )
        {
            this->evaluate_fe_function ( this->vector_solP_prev ( ), v, solP_prev_ [v] );
            this->evaluate_fe_function_gradients ( this->vector_solP_prev ( ), v, grad_solP_prev_[v] );
        }
    }

    if ( this->mode_ == PRIMAL )
    {
        for ( int v = 0; v<this->num_var_; ++v )
        {
            if ( this->L2_perturb_[v] || this->H1_perturb_[v] )
            {
                perturb_ [v].clear ( );
                perturb_prev_[v].clear ( );
                grad_perturb_[v].clear ( );
                grad_perturb_prev_[v].clear ( );
                this->evaluate_fe_function ( this->vector_perturb ( ), v, perturb_[v] );
                this->evaluate_fe_function ( this->vector_perturb_prev ( ), v, perturb_prev_[v] );
                this->evaluate_fe_function_gradients ( this->vector_perturb ( ), v, grad_perturb_ [v] );
                this->evaluate_fe_function_gradients ( this->vector_perturb_prev ( ), v, grad_perturb_prev_[v] );
            }
        }

        if ( this->vector_base_ != NULL )
        {
            for ( int v = 0; v<this->num_var_; ++v )
            {
                base_[v].clear ( );
                this->evaluate_fe_function ( this->vector_base ( ), v, base_[v] );
            }
        }
    }

    // *************************************************************************
    // DUAL VARIABLES 

    if ( this->mode_ == DUAL )
    {
        assert ( this->vector_solP_next_ != NULL );
        assert ( this->vector_solD_ != NULL );
        assert ( this->vector_solD_next_ != NULL );

        // velocity
        for ( int v = 0; v<this->num_var_; ++v )
        {
            solD_[v].clear ( );
            solD_prev_[v].clear ( );
            solD_next_[v].clear ( );
            grad_solD_[v].clear ( );
            grad_solD_prev_[v].clear ( );
            grad_solD_next_[v].clear ( );
            if ( this->vector_solD_ != NULL )
            {
                this->evaluate_fe_function ( this->vector_solD ( ), v, solD_[v] );
                this->evaluate_fe_function_gradients ( this->vector_solD ( ), v, grad_solD_[v] );
            }
            if ( this->vector_solD_next_ != NULL )
            {
                this->evaluate_fe_function ( this->vector_solD_next ( ), v, solD_next_[v] );
                this->evaluate_fe_function_gradients ( this->vector_solD_next ( ), v, grad_solD_next_[v] );
            }
            if ( this->vector_solD_prev_ != NULL )
            {
                this->evaluate_fe_function ( this->vector_solD_prev ( ), v, solD_prev_[v] );
                this->evaluate_fe_function_gradients ( this->vector_solD_prev ( ), v, grad_solD_prev_[v] );
            }

            solP_next_[v].clear ( );
            grad_solP_next_[v].clear ( );
            if ( this->vector_solP_next_ != NULL )
            {
                this->evaluate_fe_function ( this->vector_solP_next ( ), v, solP_next_[v] );
                this->evaluate_fe_function_gradients ( this->vector_solP_next ( ), v, grad_solP_next_[v] );
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::initialize_for_element_convection ( )
{
    const int num_q = this->num_quadrature_points ( );
    for ( int v = 0; v < DIM; ++v )
    {
        this->conv_[v].clear ( );
        this->conv_[v].resize ( num_q, 0. );

        this->conv_prev_[v].clear ( );
        this->conv_prev_[v].resize ( num_q, 0. );

        this->conv_next_[v].clear ( );
        this->conv_next_[v].resize ( num_q, 0. );

        this->grad_conv_[v].clear ( );
        this->grad_conv_[v].resize ( num_q );

        this->grad_conv_prev_[v].clear ( );
        this->grad_conv_prev_[v].resize ( num_q );

        this->grad_conv_next_[v].clear ( );
        this->grad_conv_next_[v].resize ( num_q );
    }

    // evaluate struct
    if ( this->convection_term_ != NULL )
    {
        for ( int d = 0; d < DIM; ++d )
        {
            for ( int q = 0; q < num_q; ++q )
            {
                this->convection_term_->evaluate ( this->t_, this->x ( q ), d, this->conv_[d][q] );
                this->convection_term_->evaluate ( this->t_ + this->dT_cn_, this->x ( q ), d, this->conv_next_[d][q] );
                this->convection_term_->evaluate ( this->t_ - this->dT_pc_, this->x ( q ), d, this->conv_prev_[d][q] );

                this->convection_term_->evaluate_grad ( this->t_, this->x ( q ), d, this->grad_conv_[d][q] );
                this->convection_term_->evaluate_grad ( this->t_ + this->dT_cn_, this->x ( q ), d, this->grad_conv_next_[d][q] );
                this->convection_term_->evaluate_grad ( this->t_ - this->dT_pc_, this->x ( q ), d, this->grad_conv_prev_[d][q] );
            }
        }
        return;
    }

    if ( this->space_conv_ == NULL )
    {
        return;
    }

    // evaluate fe function
    std::vector< Vec<DIM, DataType> > quad_points ( num_q );
    for ( int q = 0; q < num_q; ++q )
    {
        quad_points[q] = this->x ( q );
    }
    std::vector< SortedArray<int> > trial_cells ( 3 );

    for ( int v = 0; v < DIM; ++v )
    {
        if ( this->funConv_[v] != NULL )
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_conv_, *this->funConv_[v], trial_cells[0], this->conv_[v] );
        if ( this->funConv_prev_[v] != NULL )
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_conv_, *this->funConv_prev_[v], trial_cells[1], this->conv_prev_[v] );
        if ( this->funConv_next_[v] != NULL )
            evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_conv_, *this->funConv_next_[v], trial_cells[2], this->conv_next_[v] );

        for ( int d = 0; d < DIM; ++d )
        {
            std::vector<DataType> vals;
            std::vector<DataType> vals_prev;
            std::vector<DataType> vals_next;

            if ( this->funGradConv_[v][d] != NULL )
            {
                evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_conv_, *this->funGradConv_[v][d], trial_cells[0], vals );
                for ( int q = 0; q < num_q; ++q )
                {
                    this->grad_conv_[v][q][d] = vals[q];
                }
            }
            if ( this->funGradConv_prev_[v][d] != NULL )
            {
                evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_conv_, *this->funGradConv_prev_[v][d], trial_cells[1], vals_prev );
                for ( int q = 0; q < num_q; ++q )
                {
                    this->grad_conv_prev_[v][q][d] = vals_prev[q];
                }
            }

            if ( this->funGradConv_next_[v][d] != NULL )
            {
                evaluate_general_fe_function<DIM, DataType> ( quad_points, *this->space_conv_, *this->funGradConv_next_[v][d], trial_cells[2], vals_next );
                for ( int q = 0; q < num_q; ++q )
                {
                    this->grad_conv_next_[v][q][d] = vals_next[q];
                }
            }
        }
    }
}

template<int DIM, class DataType>
void MetFlowAssembler<DIM, DataType>::initialize_for_element_source ( )
{
    const int num_q = this->num_quadrature_points ( );
    this->source_.clear ( );
    this->source_.resize ( this->num_var_ );
    this->source_prev_.clear ( );
    this->source_prev_.resize ( this->num_var_ );
    this->source_next_.clear ( );
    this->source_next_.resize ( this->num_var_ );

    for ( int var = 0; var < this->num_var_; ++var )
    {
        this->source_[var].resize ( num_q, 0. );
        this->source_prev_[var].resize ( num_q, 0. );
        this->source_next_[var].resize ( num_q, 0. );
    }

    if ( this->source_term_ == NULL )
    {
        return;
    }

    for ( int var = 0; var < this->num_var_; ++var )
    {
        for ( int q = 0; q < num_q; ++q )
        {
            this->source_term_->evaluate ( var, this->t_, this->x ( q ), this->source_[var][q] );
            this->source_term_->evaluate ( var, this->t_ + this->dT_cn_, this->x ( q ), this->source_next_[var][q] );
            this->source_term_->evaluate ( var, this->t_ - this->dT_pc_, this->x ( q ), this->source_prev_[var][q] );
        }
    }
}

/// Template instanciation.
template class MetFlowAssembler<2, double>;
template class MetFlowAssembler<3, double>;

