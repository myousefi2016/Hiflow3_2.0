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

/// \author Staffan Ronnas

#ifndef __HORNER_H_
#    define __HORNER_H_

#    include <cassert>

#    include "polycoeff.h"

namespace hiflow
{

    template <class T>
    T _horner ( int deg, const PolyCoeff<T>& p, const T& x )
    {
        // Check Input

        assert ( deg >= 0 );
        assert ( p.size ( ) > deg );

        // Compute

        T result = p[deg];
        for ( int k = deg; k > 0; --k )
        {
            result = result * x + p[k - 1];
        }

        return result;
    }

} // namespace hiflow

#endif
