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

#include "linearlinetransformation.h"
#include <cmath>
#include <iostream>
#include <iomanip>

/// \author Michael Schick<br>Martin Baumann<br>Julian Kraemer

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        LinearLineTransformation<DataType>::LinearLineTransformation ( int gdim ) : CellTransformation<DataType>( gdim )
        {
        }

        template<class DataType>
        void LinearLineTransformation<DataType>::inverse ( DataType x_phy,
                                                           DataType& x_ref ) const
        {

            x_ref = ( x_phy - this->coord_vtx_[this->ij2ind ( 0, 0 )] ) / ( this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] );

        }

        template<class DataType>
        void LinearLineTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy,
                                                           DataType& x_ref, DataType& y_ref ) const
        {
            throw "This cell transformation does not support 2d inversion!\n";
        }

        template<class DataType>
        void LinearLineTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy, DataType z_phy,
                                                           DataType& x_ref, DataType& y_ref, DataType& z_ref ) const
        {
            throw "This cell transformation does not support 3d inversion!\n";
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) >= 1 );

            return coord_ref[0] * ( this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] ) + this->coord_vtx_[this->ij2ind ( 0, 0 )];
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::x_x ( const Coord& coord_ref ) const
        {
            return this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )];
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::x_y ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::x_xy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::x_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::x_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::y ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::y_x ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::y_y ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::y_xy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::y_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearLineTransformation<DataType>::y_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        bool LinearLineTransformation<DataType>::contains_reference_point ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == this->gdim_ );
            return coord_ref[0] >= 0. && coord_ref[0] <= 1.;
        }

        template class LinearLineTransformation<double>;
        template class LinearLineTransformation<float>;

    } // namespace doffem
} // namespace hiflow
