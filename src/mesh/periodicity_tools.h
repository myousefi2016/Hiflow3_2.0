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

#ifndef PERIODICITY_TOOLS_H
#    define PERIODICITY_TOOLS_H

/// \file periodicity_tools.h
/// \brief Periodicity tools.
///
/// \author Teresa Beck, Jonathan Schwegler

#    include <boost/bind.hpp>
#    include <mpi.h>
#    include "types.h"

/// \brief Class MasterSlave contains relevant information for periodicity.
/// \param master Master bound
/// \param slave Slave bound
/// \param h cell width in vicinity of periodic bound
/// \param index direction of periodicity: 0 = periodicity in x-direction,
/// 1 = periodicity in y-direction, 2 = periodicity in z-direction.

namespace hiflow
{
    namespace mesh
    {

        class MasterSlave
        {
          public:

            MasterSlave ( ) : master_ ( 0 ), slave_ ( 0 ), h_ ( 0 ), index_ ( -1 )
            {
                ;
            }

            MasterSlave ( double master, double slave, double h, int index ) :
            master_ ( master ), slave_ ( slave ), h_ ( h ), index_ ( index )
            {
                ;
            }

            double master ( ) const
            {
                return master_;
            }

            double slave ( ) const
            {
                return slave_;
            }

            double h ( ) const
            {
                return h_;
            }

            int index ( ) const
            {
                return index_;
            }

            void set_h ( double new_h )
            {
                h_ = new_h;
            }

          private:
            double master_;
            double slave_;
            double h_;
            int index_;
        };

        /// \brief Map points from slave boundary to master boundary. See also MasterSlave
        std::vector<Coordinate> periodify ( const std::vector<Coordinate>& coordinates,
                                            const GDim gdim, const std::vector<MasterSlave>& period );

        /// \brief Reverse the process of periodify. This function will only work
        /// correctly if h of MasterSlave is chosen "good". Also the mesh should be
        /// refined enough. ATTENTION: The coordinates have to form an entity.
        /// To ensure this only use coordinates obtained via entity.get_coordinates(std::vector<T>& coords).
        std::vector<Coordinate> unperiodify ( const std::vector<Coordinate>& coordinates,
                                              const GDim gdim, const std::vector<MasterSlave>& period );

    } //namespace mesh
} //namespace hiflow
#endif
