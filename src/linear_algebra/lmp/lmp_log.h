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

// To use the global HiFlow3 logging utility
// please uncomment the following line
#include "common/log.h"

#ifndef HIFLOW_LOG_H

// Use local log funictions
#    ifndef _LOG_H_
#        define _LOG_H_

#        ifndef NDEBUG

// Levels

#            define LOG_DEBUG(lvl, stream) { if (lvl <= DEBUG_LEVEL) \
      std::cout << "lmpLAtoolbox debug level:" << lvl << " - "<< stream << std::endl; }

#        else

#            define LOG_DEBUG(lvl, stream) ;

#        endif

#        define LOG_ERROR(stream) { std::cout << "lmpLAtoolbox error:" << stream << std::endl; }

#        define LOG_INFO(name, stream) { std::cout << "lmpLAtoolbox info:" << stream << std::endl; }

#    endif // Use local log funictions

#endif
