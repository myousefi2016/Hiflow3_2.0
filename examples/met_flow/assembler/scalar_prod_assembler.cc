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

#include "scalar_prod_assembler.h"

template<int DIM, int VARDIM, class DataType>
ScalarProdAssembler<DIM, VARDIM, DataType>::ScalarProdAssembler ( )
{
    this->L2_.resize ( VARDIM, true );
    this->H1_.resize ( VARDIM, false );
    this->H2_.resize ( VARDIM, false );
    this->belongs_to_field_.resize ( VARDIM, false );
}

template<int DIM, int VARDIM, class DataType>
void ScalarProdAssembler<DIM, VARDIM, DataType>::initialize_for_element ( const Element<DataType>& element, const Quadrature<DataType>& quadrature )
{
    AssemblyAssistant<DIM, DataType>::initialize_for_element ( element, quadrature );

    for ( int v = 0; v < VARDIM; v++ )
    {
        this->left_[v].clear ( );
        this->right_[v].clear ( );
        this->grad_left_[v].clear ( );
        this->grad_right_[v].clear ( );
        this->hess_left_[v].clear ( );
        this->hess_right_[v].clear ( );

        this->evaluate_fe_function ( this->vector_left ( ), v, this->left_[v] );
        this->evaluate_fe_function ( this->vector_right ( ), v, this->right_[v] );

        if ( this->H1_[v] || this->H2_[v] )
        {
            this->evaluate_fe_function_gradients ( this->vector_left ( ), v, this->grad_left_[v] );
            this->evaluate_fe_function_gradients ( this->vector_right ( ), v, this->grad_right_[v] );
        }
        if ( this->H2_[v] )
        {
            this->evaluate_fe_function_hessians ( this->vector_left ( ), v, this->hess_left_[v] );
            this->evaluate_fe_function_hessians ( this->vector_right ( ), v, this->hess_right_[v] );
        }
    }
}

/// Template instanciation.
template class ScalarProdAssembler<2, 1, double>;
template class ScalarProdAssembler<3, 1, double>;
template class ScalarProdAssembler<2, 2, double>;
template class ScalarProdAssembler<3, 2, double>;
template class ScalarProdAssembler<2, 3, double>;
template class ScalarProdAssembler<3, 3, double>;
template class ScalarProdAssembler<2, 4, double>;
template class ScalarProdAssembler<3, 4, double>;
template class ScalarProdAssembler<2, 5, double>;
template class ScalarProdAssembler<3, 5, double>;
template class ScalarProdAssembler<2, 6, double>;
template class ScalarProdAssembler<3, 6, double>;
template class ScalarProdAssembler<2, 7, double>;
template class ScalarProdAssembler<3, 7, double>;
