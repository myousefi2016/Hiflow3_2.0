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

#ifndef HIFLOW_COMMON_POINTERS_H
#    define HIFLOW_COMMON_POINTERS_H

#    include <boost/scoped_array.hpp>
#    include <boost/scoped_ptr.hpp>
#    include <boost/shared_ptr.hpp>
#    include <vector>

/// \file This file contains some useful types for dealing with smart
/// pointers, and other pointer-related operations.

/// \author Staffan Ronnas

namespace hiflow
{
    /// Include scoped_ptr template from boost

    template <typename T>
    struct ScopedPtr
    {
        typedef boost::scoped_ptr<T> Type;
    };

    /// Include scoped_array template from boost

    template <typename T>
    struct ScopedArray
    {
        typedef boost::scoped_array<T> Type;
    };

    /// Include shared_ptr from std::tr1

    template <typename T>
    struct SharedPtr
    {
        typedef boost::shared_ptr<T> Type;
    };

    /// Conversion from vector to pointer to first element.

    template <class T>
    T* vec2ptr ( std::vector<T>& vec )
    {
        return vec.empty ( ) ? 0 : &( vec[0] );
    }

    /// Conversion from const vector to pointer to first element.

    template <class T>
    const T* vec2ptr ( const std::vector<T>& vec )
    {
        return vec.empty ( ) ? 0 : &( vec[0] );
    }

}

#endif
