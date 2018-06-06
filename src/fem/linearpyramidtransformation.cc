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

#include "linearpyramidtransformation.h"
#include <cmath>
#include <iomanip>
#include <cassert>

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        LinearPyramidTransformation<DataType>::LinearPyramidTransformation ( int gdim ) : CellTransformation<DataType>( gdim )
        {
        }

        template<class DataType>
        void LinearPyramidTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy,
                                                              DataType& x_ref, DataType& y_ref ) const
        {
            throw "This cell transformation does not support 2d inversion!\n";
        }

        template<class DataType>
        void LinearPyramidTransformation<DataType>::inverse ( DataType x_phy,
                                                              DataType& x_ref ) const
        {
            throw "This cell transformation does not support 1d inversion!\n";
        }

        template<class DataType>
        void LinearPyramidTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy, DataType z_phy,
                                                              DataType& x_ref, DataType& y_ref, DataType& z_ref ) const
        {

            //x = x0 + (x1 - x0) \xi + (x3 - x0) \eta + (x4 - 0.5*x1 - 0.5*x3) \zeta
            //y = y0 + (y1 - y0) \xi + (y3 - y0) \eta + (y4 - 0.5*y1 - 0.5*y3) \zeta
            //z = z0 + (z1 - z0) \xi + (z3 - z0) \eta + (z4 - 0.5*z1 - 0.5*z3) \zeta
            //only works for crystalline element

            DataType a11 = this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )];
            DataType a12 = this->coord_vtx_[this->ij2ind ( 3, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )];
            DataType a13 = this->coord_vtx_[this->ij2ind ( 4, 0 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 0 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 0 )];

            DataType a21 = this->coord_vtx_[this->ij2ind ( 1, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )];
            DataType a22 = this->coord_vtx_[this->ij2ind ( 3, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )];
            DataType a23 = this->coord_vtx_[this->ij2ind ( 4, 1 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 1 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 1 )];

            DataType a31 = this->coord_vtx_[this->ij2ind ( 1, 2 )] - this->coord_vtx_[this->ij2ind ( 0, 2 )];
            DataType a32 = this->coord_vtx_[this->ij2ind ( 3, 2 )] - this->coord_vtx_[this->ij2ind ( 0, 2 )];
            DataType a33 = this->coord_vtx_[this->ij2ind ( 4, 2 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 2 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 2 )];

            DataType det = a11 * ( a33 * a22 - a32 * a23 )
                    - a21 * ( a33 * a12 - a32 * a13 )
                    + a31 * ( a23 * a12 - a22 * a13 );

            assert ( det != 0.0 );

            x_ref = ( 1.0 / det ) * ( ( a33 * a22 - a32 * a23 ) * ( x_phy - this->coord_vtx_[this->ij2ind ( 0, 0 )] )
                    - ( a33 * a12 - a32 * a13 ) * ( y_phy - this->coord_vtx_[this->ij2ind ( 0, 1 )] )
                    + ( a23 * a12 - a22 * a13 ) * ( z_phy - this->coord_vtx_[this->ij2ind ( 0, 2 )] ) );

            y_ref = ( 1.0 / det ) * ( -( a33 * a21 - a31 * a23 ) * ( x_phy - this->coord_vtx_[this->ij2ind ( 0, 0 )] )
                    + ( a33 * a11 - a31 * a13 ) * ( y_phy - this->coord_vtx_[this->ij2ind ( 0, 1 )] )
                    - ( a23 * a11 - a21 * a13 ) * ( z_phy - this->coord_vtx_[this->ij2ind ( 0, 2 )] ) );

            z_ref = ( 1.0 / det ) * ( ( a32 * a21 - a31 * a22 ) * ( x_phy - this->coord_vtx_[this->ij2ind ( 0, 0 )] )
                    - ( a32 * a11 - a31 * a12 ) * ( y_phy - this->coord_vtx_[this->ij2ind ( 0, 1 )] )
                    + ( a22 * a11 - a21 * a12 ) * ( z_phy - this->coord_vtx_[this->ij2ind ( 0, 2 )] ) );

        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return +this->coord_vtx_[this->ij2ind ( 0, 0 )]
                    + ( this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] ) * coord_ref[0]
                    + ( this->coord_vtx_[this->ij2ind ( 3, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] ) * coord_ref[1]
                    + ( this->coord_vtx_[this->ij2ind ( 4, 0 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 0 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 0 )] ) * coord_ref[2];

        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 1, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 3, 0 )] - this->coord_vtx_[this->ij2ind ( 0, 0 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_z ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 4, 0 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 0 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 0 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_xy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_xz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_yz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::x_zz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return +this->coord_vtx_[this->ij2ind ( 0, 1 )]
                    + ( this->coord_vtx_[this->ij2ind ( 1, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )] ) * coord_ref[0]
                    + ( this->coord_vtx_[this->ij2ind ( 3, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )] ) * coord_ref[1]
                    + ( this->coord_vtx_[this->ij2ind ( 4, 1 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 1 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 1 )] ) * coord_ref[2];
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 1, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 3, 1 )] - this->coord_vtx_[this->ij2ind ( 0, 1 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_z ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 4, 1 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 1 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 1 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_xy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_xz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_yz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::y_zz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return +this->coord_vtx_[this->ij2ind ( 0, 2 )]
                    + ( this->coord_vtx_[this->ij2ind ( 1, 2 )] - this->coord_vtx_[this->ij2ind ( 0, 2 )] ) * coord_ref[0]
                    + ( this->coord_vtx_[this->ij2ind ( 3, 2 )] - this->coord_vtx_[this->ij2ind ( 0, 2 )] ) * coord_ref[1]
                    + ( this->coord_vtx_[this->ij2ind ( 4, 2 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 2 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 2 )] ) * coord_ref[2];
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 1, 2 )] - this->coord_vtx_[this->ij2ind ( 0, 2 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 3, 2 )] - this->coord_vtx_[this->ij2ind ( 0, 2 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_z ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );

            return (this->coord_vtx_[this->ij2ind ( 4, 2 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 1, 2 )] - 0.5 * this->coord_vtx_[this->ij2ind ( 3, 2 )] );
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_xy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_xz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_yz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType LinearPyramidTransformation<DataType>::z_zz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        bool LinearPyramidTransformation<DataType>::contains_reference_point ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == this->gdim_ );

            std::setprecision ( 2 * sizeof (DataType ) );
            return trunc ( coord_ref[0] ) >= 0. + 0.5 * trunc ( coord_ref[2] ) && trunc ( coord_ref[0] ) <= 1. - 0.5 * trunc ( coord_ref[2] )
                    && trunc ( coord_ref[1] ) >= 0. + 0.5 * trunc ( coord_ref[2] ) && trunc ( coord_ref[1] ) <= 1. - 0.5 * trunc ( coord_ref[2] )
                    && trunc ( coord_ref[2] ) >= 0. && trunc ( coord_ref[2] ) <= 1.;
        }

        template class LinearPyramidTransformation<double>;
        template class LinearPyramidTransformation<float>;

    } // namespace doffem
} // namespace hiflow
