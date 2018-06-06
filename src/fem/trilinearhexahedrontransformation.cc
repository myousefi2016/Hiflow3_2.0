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

#include "trilinearhexahedrontransformation.h"
#include <cmath>
#include <iomanip>

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

namespace hiflow
{
    namespace doffem
    {

        // Reordering of vertices to make transformation coorespond to mesh
        // ordering, with (0,0,0) mapped to vertex 0, and (1,1,1) mapped to vertex 7.

        template<class DataType>
        TriLinearHexahedronTransformation<DataType>::TriLinearHexahedronTransformation ( int gdim ) : CellTransformation<DataType>( gdim )
        {
        }

        template<class DataType>
        void TriLinearHexahedronTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy,
                                                                    DataType& x_ref, DataType& y_ref ) const
        {
            throw "This cell transformation does not support 2d inversion!\n";
        }

        template<class DataType>
        void TriLinearHexahedronTransformation<DataType>::inverse ( DataType x_phy,
                                                                    DataType& x_ref ) const
        {
            throw "This cell transformation does not support 1d inversion!\n";
        }

        template<class DataType>
        void TriLinearHexahedronTransformation<DataType>::inverse ( DataType x_phy, DataType y_phy, DataType z_phy,
                                                                    DataType& x_ref, DataType& y_ref, DataType& z_ref ) const
        {
            this->inverse_newton_3d ( x_phy, y_phy, z_phy, x_ref, y_ref, z_ref );
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_1 = coord_ref[1];
            const DataType coord_2 = coord_ref[2];

            return +this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( 1. - coord_0 ) * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 1, 0 )] * coord_0 * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_0 * coord_1 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 3, 0 )] * ( 1. - coord_0 ) * coord_1 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 4, 0 )] * ( 1. - coord_0 ) * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 5, 0 )] * coord_0 * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 0 )] * coord_0 * coord_1 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 7, 0 )] * ( 1. - coord_0 ) * coord_1 * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_1 = coord_ref[1];
            const DataType coord_2 = coord_ref[2];

            return -this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 1, 0 )] * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_1 * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 3, 0 )] * coord_1 * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 4, 0 )] * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 5, 0 )] * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 0 )] * coord_1 * coord_2
                    - this->coord_vtx_[this->ij2ind ( 7, 0 )] * coord_1 * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_2 = coord_ref[2];

            return -this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( 1. - coord_0 ) * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 1, 0 )] * coord_0 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_0 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 3, 0 )] * ( 1. - coord_0 ) * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 4, 0 )] * ( 1. - coord_0 ) * coord_2
                    - this->coord_vtx_[this->ij2ind ( 5, 0 )] * coord_0 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 0 )] * coord_0 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 7, 0 )] * ( 1. - coord_0 ) * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_z ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_1 = coord_ref[1];

            return -this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( 1. - coord_0 ) * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 1, 0 )] * coord_0 * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_0 * coord_1
                    - this->coord_vtx_[this->ij2ind ( 3, 0 )] * ( 1. - coord_0 ) * coord_1
                    + this->coord_vtx_[this->ij2ind ( 4, 0 )] * ( 1. - coord_0 ) * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 5, 0 )] * coord_0 * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 6, 0 )] * coord_0 * coord_1
                    + this->coord_vtx_[this->ij2ind ( 7, 0 )] * ( 1. - coord_0 ) * coord_1;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_xy ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_2 = coord_ref[2];

            return +this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 1, 0 )] * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 0 )] * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 3, 0 )] * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 4, 0 )] * coord_2
                    - this->coord_vtx_[this->ij2ind ( 5, 0 )] * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 0 )] * coord_2
                    - this->coord_vtx_[this->ij2ind ( 7, 0 )] * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_xz ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_1 = coord_ref[1];

            return +this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 1, 0 )] * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_1
                    + this->coord_vtx_[this->ij2ind ( 3, 0 )] * coord_1
                    - this->coord_vtx_[this->ij2ind ( 4, 0 )] * ( 1. - coord_1 ) *
                    + this->coord_vtx_[this->ij2ind ( 5, 0 )] * ( 1. - coord_1 ) *
                    + this->coord_vtx_[this->ij2ind ( 6, 0 )] * coord_1
                    - this->coord_vtx_[this->ij2ind ( 7, 0 )] * coord_1;

        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_yz ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];

            return +this->coord_vtx_[this->ij2ind ( 0, 0 )] * ( 1. - coord_0 )
                    + this->coord_vtx_[this->ij2ind ( 1, 0 )] * coord_0
                    - this->coord_vtx_[this->ij2ind ( 2, 0 )] * coord_0
                    - this->coord_vtx_[this->ij2ind ( 3, 0 )] * ( 1. - coord_0 )
                    - this->coord_vtx_[this->ij2ind ( 4, 0 )] * ( 1. - coord_0 )
                    - this->coord_vtx_[this->ij2ind ( 5, 0 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 6, 0 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 7, 0 )] * ( 1. - coord_0 );
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::x_zz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_1 = coord_ref[1];
            const DataType coord_2 = coord_ref[2];

            return +this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( 1. - coord_0 ) * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 1, 1 )] * coord_0 * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_0 * coord_1 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 3, 1 )] * ( 1. - coord_0 ) * coord_1 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 4, 1 )] * ( 1. - coord_0 ) * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 5, 1 )] * coord_0 * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 1 )] * coord_0 * coord_1 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 7, 1 )] * ( 1. - coord_0 ) * coord_1 * coord_2;

        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_1 = coord_ref[1];
            const DataType coord_2 = coord_ref[2];

            return -this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 1, 1 )] * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_1 * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 3, 1 )] * coord_1 * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 4, 1 )] * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 5, 1 )] * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 1 )] * coord_1 * coord_2
                    - this->coord_vtx_[this->ij2ind ( 7, 1 )] * coord_1 * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_2 = coord_ref[2];

            return -this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( 1. - coord_0 ) * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 1, 1 )] * coord_0 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_0 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 3, 1 )] * ( 1. - coord_0 ) * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 4, 1 )] * ( 1. - coord_0 ) * coord_2
                    - this->coord_vtx_[this->ij2ind ( 5, 1 )] * coord_0 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 1 )] * coord_0 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 7, 1 )] * ( 1. - coord_0 ) * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_z ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_1 = coord_ref[1];

            return -this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( 1. - coord_0 ) * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 1, 1 )] * coord_0 * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_0 * coord_1
                    - this->coord_vtx_[this->ij2ind ( 3, 1 )] * ( 1. - coord_0 ) * coord_1
                    + this->coord_vtx_[this->ij2ind ( 4, 1 )] * ( 1. - coord_0 ) * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 5, 1 )] * coord_0 * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 6, 1 )] * coord_0 * coord_1
                    + this->coord_vtx_[this->ij2ind ( 7, 1 )] * ( 1. - coord_0 ) * coord_1;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_xy ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_2 = coord_ref[2];

            return +this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 1, 1 )] * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 1 )] * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 3, 1 )] * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 4, 1 )] * coord_2
                    - this->coord_vtx_[this->ij2ind ( 5, 1 )] * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 1 )] * coord_2
                    - this->coord_vtx_[this->ij2ind ( 7, 1 )] * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_xz ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_1 = coord_ref[1];

            return +this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 1, 1 )] * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_1
                    + this->coord_vtx_[this->ij2ind ( 3, 1 )] * coord_1
                    - this->coord_vtx_[this->ij2ind ( 4, 1 )] * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 5, 1 )] * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 6, 1 )] * coord_1
                    - this->coord_vtx_[this->ij2ind ( 7, 1 )] * coord_1;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_yz ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];

            return +this->coord_vtx_[this->ij2ind ( 0, 1 )] * ( 1. - coord_0 )
                    + this->coord_vtx_[this->ij2ind ( 1, 1 )] * coord_0
                    - this->coord_vtx_[this->ij2ind ( 2, 1 )] * coord_0
                    - this->coord_vtx_[this->ij2ind ( 3, 1 )] * ( 1. - coord_0 )
                    - this->coord_vtx_[this->ij2ind ( 4, 1 )] * ( 1. - coord_0 )
                    - this->coord_vtx_[this->ij2ind ( 5, 1 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 6, 1 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 7, 1 )] * ( 1. - coord_0 );
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::y_zz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_1 = coord_ref[1];
            const DataType coord_2 = coord_ref[2];

            return +this->coord_vtx_[this->ij2ind ( 0, 2 )] * ( 1. - coord_0 ) * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 1, 2 )] * coord_0 * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 2 )] * coord_0 * coord_1 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 3, 2 )] * ( 1. - coord_0 ) * coord_1 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 4, 2 )] * ( 1. - coord_0 ) * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 5, 2 )] * coord_0 * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 2 )] * coord_0 * coord_1 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 7, 2 )] * ( 1. - coord_0 ) * coord_1 * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_x ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_1 = coord_ref[1];
            const DataType coord_2 = coord_ref[2];

            return -this->coord_vtx_[this->ij2ind ( 0, 2 )] * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 1, 2 )] * ( 1. - coord_1 ) * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 2 )] * coord_1 * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 3, 2 )] * coord_1 * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 4, 2 )] * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 5, 2 )] * ( 1. - coord_1 ) * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 2 )] * coord_1 * coord_2
                    - this->coord_vtx_[this->ij2ind ( 7, 2 )] * coord_1 * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_y ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_2 = coord_ref[2];

            return -this->coord_vtx_[this->ij2ind ( 0, 2 )] * ( 1. - coord_0 ) * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 1, 2 )] * coord_0 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 2 )] * coord_0 * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 3, 2 )] * ( 1. - coord_0 ) * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 4, 2 )] * ( 1. - coord_0 ) * coord_2
                    - this->coord_vtx_[this->ij2ind ( 5, 2 )] * coord_0 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 2 )] * coord_0 * coord_2
                    + this->coord_vtx_[this->ij2ind ( 7, 2 )] * ( 1. - coord_0 ) * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_z ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];
            const DataType coord_1 = coord_ref[1];

            return -this->coord_vtx_[this->ij2ind ( 0, 2 )] * ( 1. - coord_0 ) * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 1, 2 )] * coord_0 * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 2, 2 )] * coord_0 * coord_1
                    - this->coord_vtx_[this->ij2ind ( 3, 2 )] * ( 1. - coord_0 ) * coord_1
                    + this->coord_vtx_[this->ij2ind ( 4, 2 )] * ( 1. - coord_0 ) * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 5, 2 )] * coord_0 * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 6, 2 )] * coord_0 * coord_1
                    + this->coord_vtx_[this->ij2ind ( 7, 2 )] * ( 1. - coord_0 ) * coord_1;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_xx ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_xy ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_2 = coord_ref[2];

            return +this->coord_vtx_[this->ij2ind ( 0, 2 )] * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 1, 2 )] * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 2, 2 )] * ( 1. - coord_2 )
                    - this->coord_vtx_[this->ij2ind ( 3, 2 )] * ( 1. - coord_2 )
                    + this->coord_vtx_[this->ij2ind ( 4, 2 )] * coord_2
                    - this->coord_vtx_[this->ij2ind ( 5, 2 )] * coord_2
                    + this->coord_vtx_[this->ij2ind ( 6, 2 )] * coord_2
                    - this->coord_vtx_[this->ij2ind ( 7, 2 )] * coord_2;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_xz ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_1 = coord_ref[1];

            return +this->coord_vtx_[this->ij2ind ( 0, 2 )] * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 1, 2 )] * ( 1. - coord_1 )
                    - this->coord_vtx_[this->ij2ind ( 2, 2 )] * coord_1
                    + this->coord_vtx_[this->ij2ind ( 3, 2 )] * coord_1
                    - this->coord_vtx_[this->ij2ind ( 4, 2 )] * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 5, 2 )] * ( 1. - coord_1 )
                    + this->coord_vtx_[this->ij2ind ( 6, 2 )] * coord_1
                    - this->coord_vtx_[this->ij2ind ( 7, 2 )] * coord_1;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_yy ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_yz ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == 3 );
            const DataType coord_0 = coord_ref[0];

            return +this->coord_vtx_[this->ij2ind ( 0, 2 )] * ( 1. - coord_0 )
                    + this->coord_vtx_[this->ij2ind ( 1, 2 )] * coord_0
                    - this->coord_vtx_[this->ij2ind ( 2, 2 )] * coord_0
                    - this->coord_vtx_[this->ij2ind ( 3, 2 )] * ( 1. - coord_0 )
                    - this->coord_vtx_[this->ij2ind ( 4, 2 )] * ( 1. - coord_0 )
                    - this->coord_vtx_[this->ij2ind ( 5, 2 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 6, 2 )] * coord_0
                    + this->coord_vtx_[this->ij2ind ( 7, 2 )] * ( 1. - coord_0 );
        }

        template<class DataType>
        DataType TriLinearHexahedronTransformation<DataType>::z_zz ( const Coord& coord_ref ) const
        {
            return 0.;
        }

        template<class DataType>
        bool TriLinearHexahedronTransformation<DataType>::contains_reference_point ( const Coord& coord_ref ) const
        {
            assert ( coord_ref.size ( ) == this->gdim_ );
            return coord_ref[0] >= 0. && coord_ref[0] <= 1.
                    && coord_ref[1] >= 0. && coord_ref[1] <= 1.
                    && coord_ref[2] >= 0. && coord_ref[2] <= 1.;
        }

        template class TriLinearHexahedronTransformation<double>;
        template class TriLinearHexahedronTransformation<float>;

    } // namespace doffem
} // namespace hiflow
