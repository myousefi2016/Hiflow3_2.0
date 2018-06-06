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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-10-06

#ifndef SRC_COMMON_SMART_POINTERS_H_
#    define SRC_COMMON_SMART_POINTERS_H_

/// TODO: Must be defined within compiler vendor flags
#    define HIFLOW_USE_BOOST_SMART_POINTERS

namespace hiflow
{

    /// Compiler switch between boost or standard smart pointers
#    ifdef HIFLOW_USE_BOOST_SMART_POINTERS

#        include <boost/scoped_array.hpp>
#        include <boost/scoped_ptr.hpp>
#        include <boost/shared_ptr.hpp>

    using boost::scoped_array;
    using boost::scoped_ptr;
    using boost::shared_ptr;

#    else

#        include <memory>

    using std::scoped_array;
    using std::scoped_ptr;
    using std::shared_ptr;

#    endif

} // hiflow

#endif /* SRC_COMMON_SMART_POINTERS_H_ */
