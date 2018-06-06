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

/// @author Philipp Gerstner
/// @date 2016-09-26

#include "slepc_environment.h"
#include "slepc.h"

namespace hiflow
{
    namespace la
    {

        SLEPcEnvironment::~SLEPcEnvironment ( )
        {
            finalize ( );
        }

        void SLEPcEnvironment::initialize ( int argc, char **argv )
        {
            if ( initialized_ ) return;
            PetscErrorCode ierr = SlepcInitialize ( &argc, &argv, NULL, NULL );
            initialized_ = true;
        }

        void SLEPcEnvironment::initialize ( )
        {
            if ( initialized_ ) return;
            PetscErrorCode ierr = SlepcInitializeNoArguments ( );
            initialized_ = true;
        }

        void SLEPcEnvironment::finalize ( )
        {
            if ( !initialized_ ) return;
            PetscErrorCode ierr = SlepcFinalize ( );
            initialized_ = false;
        }

        bool SLEPcEnvironment::initialized_ = false;

    } // namespace la
} // namespace hiflow
