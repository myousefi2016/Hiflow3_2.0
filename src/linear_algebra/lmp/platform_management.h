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

/// @author Dimitar Lukarski

#ifndef __PLATFORM_MANAGEMENT_H
#    define __PLATFORM_MANAGEMENT_H

#    include "la_global.h"

/// @brief Management of the different Platforms
/// @author Dimitar Lukarski
///
/// Initialize and Stop a specific platform

namespace hiflow
{
    namespace la
    {

        // init the platform (CPU,GPU,...)
        void init_platform ( struct SYSTEM &my_system );

        // stop the platform
        void stop_platform ( struct SYSTEM &my_system );

        // simple platform selector
        void select_platform ( struct SYSTEM &my_system,
                               enum PLATFORM &v_platform, enum IMPLEMENTATION &v_implementation,
                               enum PLATFORM &m_platform, enum IMPLEMENTATION &m_implementation,
                               enum MATRIX_FORMAT &m_format, enum MATRIX_FREE_PRECOND &m_precond );

        template <typename ValueType>
        void select_platform ( struct SYSTEM &my_system,
                               enum PLATFORM &v_platform, enum IMPLEMENTATION &v_implementation,
                               enum PLATFORM &m_platform, enum IMPLEMENTATION &m_implementation,
                               enum MATRIX_FORMAT &m_format, enum MATRIX_FREE_PRECOND &m_precond );

    } // namespace la
} // namespace hiflow

#endif
