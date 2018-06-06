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

#include "qgausstriangle.h"

namespace hiflow
{

    template<class DataType>
    QuadratureGaussTriangle<DataType>::QuadratureGaussTriangle ( ) : QuadratureType<DataType>::QuadratureType ( )
    {
        //////////////////////////////////////////////////////////////
        // Implemented sizes: (please fill up if a new size was added)
        // 1, 3, 4
        //////////////////////////////////////////////////////////////
        int total_number_of_quadratures = 3;

        int cntr = 0;
        int size = 0;

        // First resize vector fields to the number of implemented quadratures
        x_.resize ( total_number_of_quadratures );
        y_.resize ( total_number_of_quadratures );
        z_.resize ( total_number_of_quadratures );
        w_.resize ( total_number_of_quadratures );

        // Next fill vector fields with quadrature data

        // ---------------------------------------------
        size = 1;

        DataType x1[] = {
                         .333333333333333333333333333333
        };
        DataType y1[] = {
                         .333333333333333333333333333333
        };
        DataType w1[] = {
                         .5
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x1, x1 + size );
        y_[cntr].reserve ( size );
        y_[cntr].insert ( y_[cntr].begin ( ), y1, y1 + size );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w1, w1 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 3;

        DataType x3[] = {
                         .166666666666666666666666666667,
                         .166666666666666666666666666667,
                         .666666666666666666666666666667
        };
        DataType y3[] = {
                         .166666666666666666666666666667,
                         .666666666666666666666666666667,
                         .166666666666666666666666666667
        };
        DataType w3[] = {
                         .166666666666666666666666666667,
                         .166666666666666666666666666667,
                         .166666666666666666666666666667
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x3, x3 + size );
        y_[cntr].reserve ( size );
        y_[cntr].insert ( y_[cntr].begin ( ), y3, y3 + size );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w3, w3 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 4;

        DataType x4[] = {
                         .2,
                         .333333333333333333333333333334,
                         .6,
                         .2
        };
        DataType y4[] = {
                         .2,
                         .333333333333333333333333333334,
                         .2,
                         .6
        };
        DataType w4[] = {
                         .260416666666666666666666666667,
                         -0.28125,
                         .260416666666666666666666666667,
                         .260416666666666666666666666667
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x4, x4 + size );
        y_[cntr].reserve ( size );
        y_[cntr].insert ( y_[cntr].begin ( ), y4, y4 + size );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w4, w4 + size );

        index_field_[size] = cntr;
        cntr++;

    }

    template<class DataType>
    QuadratureGaussTriangle<DataType>* QuadratureGaussTriangle<DataType>::clone ( ) const
    {
        return new QuadratureGaussTriangle<DataType>( *this );
    }

    // template instanciation
    template class QuadratureGaussTriangle<double>;
    template class QuadratureGaussTriangle<float>;

} // namespace hiflow
