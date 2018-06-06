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

#include "scalar_prod_cart_assembler.h"

template<int DIM, int VARDIM, class DataType>
ScalarProdCartAssembler<DIM, VARDIM, DataType>::ScalarProdCartAssembler ( )
: ScalarProdAssembler<DIM, VARDIM, DataType>( )
{
}

template<int DIM, int VARDIM, class DataType>
void ScalarProdCartAssembler<DIM, VARDIM, DataType>::assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const
{
    const int num_q = this->num_quadrature_points ( );

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::abs ( this->detJ ( q ) );

        // loop over scalar variables
        for ( int var = 0; var < VARDIM; ++var )
        {
            if ( !this->belongs_to_field_[var] )
            {
                if ( this->L2_[var] )
                {
                    ls += wq * ( this->left_[var][q] * this->right_[var][q] ) * dJ;
                }
                if ( this->H1_[var] )
                {
                    for ( int s = 0; s < DIM; ++s )
                    {
                        ls += wq * this->grad_left_[var][q][s] * this->grad_right_[var][q][s] * dJ;
                    }
                }
                if ( this->H2_[var] )
                {

                }
            }
        }
        // loop over vector fields
        for ( int l = 0; l<this->vector_field_ind_.size ( ); ++l )
        {
            const int phi_var = this->vector_field_ind_[l][0];
            const int rad_var = this->vector_field_ind_[l][1];
            for ( int v = 0; v < DIM; ++v )
            {
                const int var = this->vector_field_ind_[l][v];

                if ( this->L2_[var] )
                {
                    ls += wq * ( this->left_[var][q] * this->right_[var][q] ) * dJ;
                }
                if ( this->H1_[var] )
                {
                    for ( int s = 0; s < DIM; ++s )
                    {
                        ls += wq * this->grad_left_[var][q][s] * this->grad_right_[var][q][s] * dJ;
                    }
                }
                if ( this->H2_[var] )
                {

                }
            }
        }
    }
}

/// Template instanciation.
template class ScalarProdCartAssembler<2, 1, double>;
template class ScalarProdCartAssembler<3, 1, double>;
template class ScalarProdCartAssembler<2, 2, double>;
template class ScalarProdCartAssembler<3, 2, double>;
template class ScalarProdCartAssembler<2, 3, double>;
template class ScalarProdCartAssembler<3, 3, double>;
template class ScalarProdCartAssembler<2, 4, double>;
template class ScalarProdCartAssembler<3, 4, double>;
template class ScalarProdCartAssembler<2, 5, double>;
template class ScalarProdCartAssembler<3, 5, double>;
template class ScalarProdCartAssembler<2, 6, double>;
template class ScalarProdCartAssembler<3, 6, double>;
template class ScalarProdCartAssembler<2, 7, double>;
template class ScalarProdCartAssembler<3, 7, double>;
