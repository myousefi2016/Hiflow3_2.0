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

/// \author Philipp Gerstner

#ifndef P4EST_TYPES_H
#    define P4EST_TYPES_H

#    include <vector>
#    include <map>
#    include <iostream>
#    include <string>
#    include "common/sorted_array.h"
#    include "config.h"
#    include "mesh/types.h"
#    include <mpi.h>
#    include "communication.h"

#    ifdef WITH_P4EST
#        include "p4est.h"
#        include "p8est.h"
#    endif

namespace hiflow
{
    namespace mesh
    {

        class pXestGhost
        {
          public:

            pXestGhost ( p4est_ghost_t* ghost4 ) : tdim_ ( 2 ), ghost4_ ( ghost4 )
            {
            }

            pXestGhost ( p8est_ghost_t* ghost8 ) : tdim_ ( 3 ), ghost8_ ( ghost8 )
            {
            }

            int num_elem ( )
            {

            }

            QuadData* get_ghost_data_ptr ( int j )
            {

            }

            QuadData* get_mirror_data_ptr ( int j )
            {

            }

          private:
            int tdim_;
            p4est_ghost_t* ghost4_;
            p8est_ghost_t* ghost8_;

        };

        class pXestQuad
        {
          public:

            pXestQuad ( p4est_quadrant_t* quad4 ) : tdim_ ( 2 ), quad4_ ( quad4 )
            {
            }

            pXestQuad ( p8est_quadrant_t* quad8 ) : tdim_ ( 3 ), quad8_ ( quad8 )
            {
            }

            int level ( )
            {
                switch ( tdim_ )
                {
                    case 2:
                        return quad4_->level;
                        break;
                    case 3:
                        return quad8_->level;
                        break;
                }
            }

          private:
            int tdim_;
            p4est_quadrant_t* quad4_;
            p8est_quadrant_t* quad8_;

        }

    }
}
#endif
