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

/// @author Chandramowli Subramanian

#ifndef HIFLOW_LINEARSOLVER_MUMPS_STRUCTURE_H_
#    define HIFLOW_LINEARSOLVER_MUMPS_STRUCTURE_H_

#    include <cstdio>

#    include "config.h"

#    ifdef WITH_MUMPS
#        include "dmumps_c.h"
#        include "smumps_c.h"
#    endif

namespace hiflow
{
    namespace la
    {

        // macros such that indices match MUMPS documentation
#    define ICNTL(I) icntl[(I)-1]
#    define CNTL(I) cntl[(I)-1]
#    define INFO(I) info[(I)-1]
#    define INFOG(I) infog[(I)-1]

        /// @brief Mumps structure for real data types in double precision.

        struct MumpsStructureD
        {
#    ifdef WITH_MUMPS
            DMUMPS_STRUC_C id_;
#    endif
            inline void apply ( );
        };

        inline void MumpsStructureD::apply ( )
        {
#    ifdef WITH_MUMPS
            dmumps_c ( &( this->id_ ) );
            assert ( this->id_.INFOG ( 1 ) == 0 );
#    else
            printf ( "MumpsStructureD::apply: No MUMPS support.\n" );
            exit ( -1 );
#    endif
        }

        /// @brief Mumps structure for real data types in single precision.
        /// @author Chandramowli Subramanian

        struct MumpsStructureS
        {
#    ifdef WITH_MUMPS
            SMUMPS_STRUC_C id_;
#    endif
            inline void apply ( );
        };

        inline void MumpsStructureS::apply ( )
        {
#    ifdef WITH_MUMPS
            smumps_c ( &( this->id_ ) );
            assert ( this->id_.INFOG ( 1 ) >= 0 );
#    else
            printf ( "MumpsStructureS::apply: No MUMPS support.\n" );
            exit ( -1 );
#    endif
        }

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_MUMPS_STRUCTURE_H_
