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

#include "qgausstetrahedron.h"

namespace hiflow
{

    template<class DataType>
    QuadratureGaussTetrahedron<DataType>::QuadratureGaussTetrahedron ( ) : QuadratureType<DataType>::QuadratureType ( )
    {
        //////////////////////////////////////////////////////////////
        // Implemented sizes: (please fill up if a new size was added)
        // 1, 5, 15
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
                         .25
        };
        DataType y1[] = {
                         .25
        };
        DataType z1[] = {
                         .25
        };
        DataType w1[] = {
                         .166666666666666666666666666667
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x1, x1 + size );
        y_[cntr].reserve ( size );
        y_[cntr].insert ( y_[cntr].begin ( ), y1, y1 + size );
        z_[cntr].reserve ( size );
        z_[cntr].insert ( z_[cntr].begin ( ), z1, z1 + size );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w1, w1 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 5;
        DataType x5[] = {
                         .166666666666666666666666666667,
                         .166666666666666666666666666667,
                         .166666666666666666666666666667,
                         .25,
                         .5
        };
        DataType y5[] = {
                         .166666666666666666666666666667,
                         .166666666666666666666666666667,
                         .5,
                         .25,
                         .166666666666666666666666666667
        };
        DataType z5[] = {
                         .166666666666666666666666666667,
                         .5,
                         .166666666666666666666666666667,
                         .25,
                         .166666666666666666666666666667
        };
        DataType w5[] = {
                         .075,
                         .075,
                         .075,
                         -.133333333333333333333333333337,
                         .075
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x5, x5 + size );
        y_[cntr].reserve ( size );
        y_[cntr].insert ( y_[cntr].begin ( ), y5, y5 + size );
        z_[cntr].reserve ( size );
        z_[cntr].insert ( z_[cntr].begin ( ), z5, z5 + size );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w5, w5 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 15;
        DataType x15[] = {
                          .25,
                          .091971078052723032788845135301,
                          .091971078052723032788845135301,
                          .091971078052723032788845135301,
                          .816057843894553934422309729399,
                          .319793627829629908387625452935,
                          .319793627829629908387625452935,
                          .319793627829629908387625452935,
                          .360412744340740183224749094130,
                          .138196601125010515179541316563,
                          .138196601125010515179541316563,
                          .361803398874989484820458683437,
                          .361803398874989484820458683437,
                          .361803398874989484820458683437,
                          .138196601125010515179541316563
        };
        DataType y15[] = {
                          .25,
                          .091971078052723032788845135301,
                          .091971078052723032788845135301,
                          .816057843894553934422309729399,
                          .091971078052723032788845135301,
                          .319793627829629908387625452935,
                          .319793627829629908387625452935,
                          .360412744340740183224749094130,
                          .319793627829629908387625452935,
                          .138196601125010515179541316563,
                          .361803398874989484820458683437,
                          .138196601125010515179541316563,
                          .361803398874989484820458683437,
                          .138196601125010515179541316563,
                          .361803398874989484820458683437
        };
        DataType z15[] = {
                          .25,
                          .091971078052723032788845135301,
                          .816057843894553934422309729399,
                          .091971078052723032788845135301,
                          .091971078052723032788845135301,
                          .319793627829629908387625452935,
                          .360412744340740183224749094130,
                          .319793627829629908387625452935,
                          .319793627829629908387625452935,
                          .361803398874989484820458683437,
                          .138196601125010515179541316563,
                          .138196601125010515179541316563,
                          .138196601125010515179541316563,
                          .361803398874989484820458683437,
                          .361803398874989484820458683437
        };
        DataType w15[] = {
                          .019753086419753086419753086420,
                          .011989513963169770001730642485,
                          .011989513963169770001730642485,
                          .011989513963169770001730642485,
                          .011989513963169770001730642485,
                          .011511367871045397546770239349,
                          .011511367871045397546770239349,
                          .011511367871045397546770239349,
                          .011511367871045397546770239349,
                          .008818342151675485008818342152,
                          .008818342151675485008818342152,
                          .008818342151675485008818342152,
                          .008818342151675485008818342152,
                          .008818342151675485008818342152,
                          .008818342151675485008818342152
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x15, x15 + size );
        y_[cntr].reserve ( size );
        y_[cntr].insert ( y_[cntr].begin ( ), y15, y15 + size );
        z_[cntr].reserve ( size );
        z_[cntr].insert ( z_[cntr].begin ( ), z15, z15 + size );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w15, w15 + size );

        index_field_[size] = cntr;
        cntr++;

    }

    template<class DataType>
    QuadratureGaussTetrahedron<DataType>* QuadratureGaussTetrahedron<DataType>::clone ( ) const
    {
        return new QuadratureGaussTetrahedron<DataType>( *this );
    }

    // template instanciation
    template class QuadratureGaussTetrahedron<double>;
    template class QuadratureGaussTetrahedron<float>;

} // namespace hiflow
