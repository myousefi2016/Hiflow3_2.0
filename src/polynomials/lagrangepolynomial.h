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

#ifndef __POLYNOMIALS_LAGRANGE_POLYNOMIAL_H_
#    define __POLYNOMIALS_LAGRANGE_POLYNOMIAL_H_

#    include <cassert>
#    include <vector>

#    include "polynomials/lagrangecoeff.h"

namespace hiflow
{

    //-------------------------------------------------
    //
    // 1D, equidistant  [0,1]
    //

    /// _lagrange_polynomial

    template<class T>
    T _lagrange_polynomial ( int deg, int i, const T& x )
    {
        assert ( ( i >= 0 )&&( i <= deg ) );

        T r = 1.;
        T p = ( T ) 1;

        T dx = deg*x;

        for ( int j = 0; j <= deg; ++j )
            if ( i != j )
            {
                r *= ( dx - j );
                p *= ( i - j );
            }

        return r / p;
    }

    /// _x

    template<class T>
    T _lagrange_polynomial_x ( int deg, int i, const T& x )
    {
        assert ( ( i >= 0 )&&( i <= deg ) );

        T r = 0.;
        T p = ( T ) 1;
        T dx = deg*x;

        for ( int j = 0; j <= deg; ++j )

            if ( i != j )
            {
                T s = 1.;

                for ( int k = 0; k <= deg; ++k )
                    if ( ( k != j ) && ( k != i ) )
                        s *= ( dx - k );

                r += s;
                p *= ( i - j );
            }

        return deg * r / p;
    }

    /// _xx

    template<class T>
    T _lagrange_polynomial_xx ( int deg, int i, const T& x )
    {
        assert ( ( i >= 0 )&&( i <= deg ) );

        T r = 0.;
        T p = ( T ) 1;
        T dx = deg*x;

        for ( int j = 0; j <= deg; ++j )

            if ( i != j )
            {
                T s = 0.;

                for ( int k = 0; k <= deg; ++k )
                    if ( ( k != i ) && ( k != j ) )
                    {
                        T t = 1.;

                        for ( int l = 0; l <= deg; ++l )
                            if ( ( l != i ) && ( l != j )&& ( l != k ) )
                                t *= ( dx - k );

                        s += t;
                    }

                r += s;
                p *= ( i - j );
            }

        return deg * deg * r / p;
    }

    template <class T>
    class LagrangePolynomial
    {
      private:

        std::vector<const LagrangeCoeff*> _lag_coeff;

      public:

        LagrangePolynomial ( );

        void reinit ( );

        T poly ( int, int, const T& ) const;
        T poly_x ( int, int, const T& ) const;
        T poly_xx ( int, int, const T& ) const;
    };

    /* --------------------------------------------------------------
    /                                                               /
    /  INLINE FUNCTIONS for LagrangePolynomial                      /
    /                                                               /
    /  ----------------------------------------------------------- */

} // namespace hiflow

#endif
