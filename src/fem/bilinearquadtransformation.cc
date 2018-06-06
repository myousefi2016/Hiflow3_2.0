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

#include "bilinearquadtransformation.h"
#include <cmath>
#include <iomanip>

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        BiLinearQuadTransformation<DataType>::BiLinearQuadTransformation ( int gdim ) : CellTransformation<DataType>( gdim )
        {
        }

        template<class DataType>
        void BiLinearQuadTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy,
                                                             DataType& x_ref, DataType& y_ref ) const
        {
            this->inverse_newton_2d ( x_phy, y_phy, x_ref, y_ref );
        }

        template<class DataType>
        void BiLinearQuadTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy, DataType z_phy,
                                                             DataType& x_ref, DataType& y_ref, DataType& z_ref ) const
        {
            throw "This cell transformation does not support 3d inversion!\n";
        }

        template<class DataType>
        void BiLinearQuadTransformation<DataType>::inverse ( DataType x_phy,
                                                             DataType& x_ref ) const
        {
            throw "This cell transformation does not support 1d inversion!\n";
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_1 = coord_ref[1];

            return this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( 1. - coord_0 ) * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 1, 0 )] * coord_0 * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_0 * coord_1
                    + this->coord_vtx_[this->ij2ind ( 3, 0 )] * ( 1. - coord_0 ) * coord_1;
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::x_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );
            const DataType coord_1 = coord_ref[1];

            return this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( -1. + coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 1, 0 )] * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_1
                    - this->coord_vtx_[this->ij2ind ( 3, 0 )] * coord_1;
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::x_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );
            const DataType coord_0 = coord_ref[0];

            return this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( -1. + coord_0 )
                    - this->coord_vtx_[this->ij2ind ( 1, 0 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 3, 0 )] * ( 1. - coord_0 );
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::x_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::x_xy ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );
            return this->coord_vtx_[this->ij2ind ( 0, 0 )]
                    - this->coord_vtx_[this->ij2ind ( 1, 0 )]
                    + this->coord_vtx_[this->ij2ind ( 2, 0 )]
                    - this->coord_vtx_[this->ij2ind ( 3, 0 )];
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::x_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_1 = coord_ref[1];

            return this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( 1. - coord_0 ) * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 1, 1 )] * coord_0 * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_0 * coord_1
                    + this->coord_vtx_[this->ij2ind ( 3, 1 )] * ( 1. - coord_0 ) * coord_1;
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::y_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );
            const DataType coord_1 = coord_ref[1];

            return this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( -1. + coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 1, 1 )] * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_1
                    - this->coord_vtx_[this->ij2ind ( 3, 1 )] * coord_1;
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::y_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );
            const DataType coord_0 = coord_ref[0];

            return this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( -1. + coord_0 )
                    - this->coord_vtx_[this->ij2ind ( 1, 1 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 3, 1 )] * ( 1. - coord_0 );
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::y_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::y_xy ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );

            return this->coord_vtx_[this->ij2ind ( 0, 1 )]
                    - this->coord_vtx_[this->ij2ind ( 1, 1 )]
                    + this->coord_vtx_[this->ij2ind ( 2, 1 )]
                    - this->coord_vtx_[this->ij2ind ( 3, 1 )];
        }

        template<class DataType>
        DataType BiLinearQuadTransformation<DataType>::y_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        // TODO avoid trunc

        template<class DataType>
        bool BiLinearQuadTransformation<DataType>::contains_reference_point ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == this->gdim_ );
            return coord_ref[0] >= 0. && coord_ref[0] <= 1.
                    && coord_ref[1] >= 0. && coord_ref[1] <= 1.;
        }

        template class BiLinearQuadTransformation<double>;
        template class BiLinearQuadTransformation<float>;

    } // namespace doffem
} // namespace hiflow

