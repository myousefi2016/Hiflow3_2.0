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

#include "mass_matrix_cart_assembler.h"

/// Global mass matrix

template<int DIM, int VARDIM, class DataType>
void MassMatrixCartAssembler<DIM, VARDIM, DataType>::assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const
{
    const int num_q = this->num_quadrature_points ( );
    const int total_dofs = this->num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const DataType wq = this->w ( q );
        const DataType dJ = std::fabs ( this->detJ ( q ) );

        for ( int var = 0; var < VARDIM; var++ )
        {
            for ( int j = 0; j<this->num_dofs ( var ); j++ )
            {
                for ( int i = 0; i<this->num_dofs ( var ); i++ )
                {

                    lm ( this->dof_index ( i, var ), this->dof_index ( j, var ) ) += wq
                            * this->phi ( j, q, var )
                            * this->phi ( i, q, var )
                            * dJ;
                }
            }
        }
    }
}

/// Template instanciation.
template class MassMatrixCartAssembler<2, 1, double>;
template class MassMatrixCartAssembler<3, 1, double>;
template class MassMatrixCartAssembler<2, 2, double>;
template class MassMatrixCartAssembler<3, 2, double>;
template class MassMatrixCartAssembler<2, 3, double>;
template class MassMatrixCartAssembler<3, 3, double>;
template class MassMatrixCartAssembler<2, 4, double>;
template class MassMatrixCartAssembler<3, 4, double>;
template class MassMatrixCartAssembler<2, 5, double>;
template class MassMatrixCartAssembler<3, 5, double>;
template class MassMatrixCartAssembler<2, 6, double>;
template class MassMatrixCartAssembler<3, 6, double>;
template class MassMatrixCartAssembler<2, 7, double>;
template class MassMatrixCartAssembler<3, 7, double>;
