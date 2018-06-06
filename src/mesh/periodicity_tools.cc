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

#include "periodicity_tools.h"
#include <cmath>

namespace hiflow
{
    namespace mesh
    {

        std::vector<Coordinate> periodify ( const std::vector<Coordinate>& coordinates, const GDim gdim, const std::vector<MasterSlave>& period )
        {
            int num_coords = coordinates.size ( ) / gdim;
            std::vector<Coordinate> temp ( coordinates );
            const Coordinate eps = 1.0e-13;
            int size = period.size ( );
            for ( int k = 0; k < size; k++ )
            {
                for ( int i = 0; i < num_coords; ++i )
                {
                    if ( std::abs ( temp[i * gdim + period[k].index ( )] - period[k].slave ( ) ) < eps )
                    {
                        temp[i * gdim + period[k].index ( )] = period[k].master ( );
                    }
                }
            }

            return temp;
        }

        std::vector<Coordinate> unperiodify ( const std::vector<Coordinate>& coordinates, const GDim gdim, const std::vector<MasterSlave>& period )
        {
            int num_coords = coordinates.size ( ) / gdim;
            std::vector<Coordinate> temp ( coordinates );

            int size = period.size ( );
            for ( int k = 0; k < size; k++ )
            {
                bool slave_boundary = false;
                Coordinate boundary_band = period[k].h ( );

                for ( int i = 0; i < num_coords; ++i )
                {
                    if ( std::abs ( temp[i * gdim + period[k].index ( )] - period[k].slave ( ) ) < boundary_band )
                    {
                        slave_boundary = true;
                        break;
                    }
                }

                // modify coordinates if necessary
                const Coordinate eps = 1.0e-13;
                if ( slave_boundary )
                {
                    for ( int i = 0; i < num_coords; ++i )
                    {
                        if ( std::abs ( temp[i * gdim + period[k].index ( )] - period[k].master ( ) ) < eps )
                        {
                            temp[i * gdim + period[k].index ( )] = period[k].slave ( );
                        }
                    }
                }
            }

            return temp;
        }

    } //namespace mesh
} //namespace hiflow
