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

#include "lineartriangletransformation.h"
#include <cmath>
#include <iostream>
#include <iomanip>

/// \author Michael Schick<br>Martin Baumann

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        LinearTriangleTransformation<DataType>::LinearTriangleTransformation ( int gdim ) : CellTransformation<DataType>( gdim )
        {
        }

        template<class DataType>
        void LinearTriangleTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy,
                                                               DataType& x_ref, DataType& y_ref ) const
        {
            DataType a11 = this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )];
            DataType a12 = this->coord_vtx_[this->ij2ind ( 2, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )];
            DataType a21 = this->coord_vtx_[this->ij2ind ( 1, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )];
            DataType a22 = this->coord_vtx_[this->ij2ind ( 2, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )];

            DataType det = a11 * a22 - a21*a12;

            assert ( det != 0.0 );

            x_ref = ( 1.0 / det ) * ( a22 * ( x_phy - this->coord_vtx_[this->ij2ind ( 0, 0 )] )
                    - a12 * ( y_phy - this->coord_vtx_[this->ij2ind ( 0, 1 )] ) );
            y_ref = ( 1.0 / det ) * ( -a21 * ( x_phy - this->coord_vtx_[this->ij2ind ( 0, 0 )] )
                    + a11 * ( y_phy - this->coord_vtx_[this->ij2ind ( 0, 1 )] ) );
        }

        template<class DataType>
        void LinearTriangleTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy, DataType z_phy,
                                                               DataType& x_ref, DataType& y_ref, DataType& z_ref ) const
        {
            throw "This cell transformation does not support 3d inversion!\n";
        }

        template<class DataType>
        void LinearTriangleTransformation<DataType>::inverse ( DataType x_phy,
                                                               DataType& x_ref ) const
        {
            throw "This cell transformation does not support 1d inversion!\n";
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );

            return (this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] ) * coord_ref[0]
                    + ( this->coord_vtx_[this->ij2ind ( 2, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] ) * coord_ref[1]
                    + this->coord_vtx_[this->ij2ind ( 0, 0 )];
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::x_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );

            return (this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] );
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::x_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );

            return (this->coord_vtx_[this->ij2ind ( 2, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] );
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::x_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::x_xy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::x_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );

            return (this->coord_vtx_[this->ij2ind ( 1, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )] ) * coord_ref[0]
                    + ( this->coord_vtx_[this->ij2ind ( 2, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )] ) * coord_ref[1]
                    + this->coord_vtx_[this->ij2ind ( 0, 1 )];
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::y_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );

            return (this->coord_vtx_[this->ij2ind ( 1, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )] );
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::y_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 2 );

            return (this->coord_vtx_[this->ij2ind ( 2, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )] );
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::y_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::y_xy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearTriangleTransformation<DataType>::y_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        bool LinearTriangleTransformation<DataType>::contains_reference_point ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == this->gdim_ );
            return coord_ref[0] >= 0. && coord_ref[0] <= 1.
                    && coord_ref[1] >= 0. && coord_ref[1] <= 1. - coord_ref[0];
        }

        template class LinearTriangleTransformation<double>;
        template class LinearTriangleTransformation<float>;

    } // namespace doffem
} // namespace hiflow
