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

/// @author Nico Trost, Benedikt Galler, Dimitar Lukarski

#ifndef __PLATFORM_MANAGMENT_OPENCL_H
#    define __PLATFORM_MANAGMENT_OPENCL_H

#    include "config.h"

namespace hiflow
{
    namespace la
    {

        void print_platform_opencl_info ( void );

        void init_platform_opencl ( struct SYSTEM &my_system );

        void init_platform_opencl ( struct SYSTEM &my_system, int opencl_plat, int opencl_dev );

        void stop_platform_opencl ( struct SYSTEM &my_system );

    }
}
#endif
