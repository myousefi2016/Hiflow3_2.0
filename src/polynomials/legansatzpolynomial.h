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

#ifndef __POLYNOMIALS_LEG_ANSATZ_POLYNOMIAL_H_
#    define __POLYNOMIALS_LEG_ANSATZ_POLYNOMIAL_H_

#    include <cassert>
#    include <cmath>
#    include <vector>

#    include "legendrepolynomial.h"

namespace hiflow
{

    template <class T>
    class LegAnsatzPolynomial
    {
      private:

        int _degree;
        std::vector<T> _sqrt4ip2;

      public:

        LegAnsatzPolynomial ( );

        void reinit ( int );

        inline T poly ( int, const T& ) const;
        inline T poly_x ( int, const T& ) const;
        inline T poly_xx ( int, const T& ) const;

        inline int deg ( ) const;

    };

    /* --------------------------------------------------------------
    /                                                               /
    /  INLINE FUNCTIONS for LegAnsatzPolynomial                     /
    /                                                               /
    /  ----------------------------------------------------------- */

    template<class T>
    LegAnsatzPolynomial<T>::LegAnsatzPolynomial ( )
    {
        reinit ( 0 );
    }

    template<class T>
    void LegAnsatzPolynomial<T>::reinit ( int degr )
    {
        _degree = degr;

        if ( !_sqrt4ip2.empty ( ) )
            _sqrt4ip2.erase ( _sqrt4ip2.begin ( ), _sqrt4ip2.end ( ) );

        _sqrt4ip2.reserve ( deg ( ) );

        for ( int i = 0; i < deg ( ); ++i ) _sqrt4ip2.push_back ( sqrt ( ( T ) 4 * i + 2 ) );
    }

    template<class T>
    T LegAnsatzPolynomial<T>::poly ( int i, const T& x ) const
    {
        assert ( ( i >= 0 ) && ( i <= deg ( ) ) );

        if ( i == 0 )
            return (1. - x );
        else if ( i == deg ( ) )
            return x;
        else
            return (_legendre_polynomial<T>( i + 1, 2 * x - 1 ) -
                _legendre_polynomial<T>( i - 1, 2 * x - 1 ) ) / ( _sqrt4ip2[i] );
    }

    template<class T>
    T LegAnsatzPolynomial<T>::poly_x ( int i, const T& x ) const
    {
        assert ( ( i >= 0 ) && ( i <= deg ( ) ) );

        if ( i == 0 )
            return -1.;
        else if ( i == deg ( ) )
            return 1.;
        else
            return (_legendre_polynomial<T>( i, 2 * x - 1 ) * _sqrt4ip2[i] );
    }

    template<class T>
    T LegAnsatzPolynomial<T>::poly_xx ( int i, const T& x ) const
    {
        assert ( ( i >= 0 ) && ( i <= deg ( ) ) );

        if ( i == 0 )
            return 0.;
        else if ( i == deg ( ) )
            return 0.;
        else
            return (_legendre_polynomial_x<T>( i, 2 * x - 1 ) * 2 * _sqrt4ip2[i] );
    }

    template<class T>
    int LegAnsatzPolynomial<T>::deg ( ) const
    {
        return _degree;
    }

} // namespace hiflow

#endif
