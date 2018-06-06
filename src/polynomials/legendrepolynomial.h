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

/// \author Martin Baumann

#ifndef __POLYNOMIALS_LEGENDRE_POLYNOMIAL_H_
#    define __POLYNOMIALS_LEGENDRE_POLYNOMIAL_H_

#    include <cassert>
#    include <vector>

namespace hiflow
{

    template<class T>
    T _legendre_polynomial ( int i, const T& x )
    {
        if ( i == 0 )
            return 1.;
        else if ( i == 1 )
            return x;
        else return ( ( 2 * i - 1 ) * x * _legendre_polynomial<T>( i - 1, x )
                - ( i - 1 ) * _legendre_polynomial<T>( i - 2, x ) ) / i;
    }

    template<class T>
    T _legendre_polynomial_x ( int i, const T& x )
    {
        if ( i == 0 )
            return 0.;
        else if ( i == 1 )
            return 1.;
        else return (2 * i - 1 ) * _legendre_polynomial<T>( i - 1, x )
            + _legendre_polynomial_x<T>( i - 2, x );
    }

} // namespace hiflow

#endif
