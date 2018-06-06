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

#include "horner.h"

#include "lagrangecoeff1.h"
#include "lagrangecoeff2.h"
#include "lagrangecoeff3.h"
#include "lagrangecoeff4.h"
#include "lagrangecoeff5.h"
#include "lagrangecoeff6.h"
#include "lagrangecoeff7.h"
#include "lagrangecoeff8.h"
#include "lagrangecoeff9.h"
#include "lagrangecoeff10.h"
#include "lagrangecoeff11.h"
#include "lagrangecoeff12.h"
#include "lagrangecoeff13.h"
#include "lagrangecoeff14.h"
#include "lagrangecoeff15.h"
#include "lagrangecoeff16.h"
#include "lagrangecoeff17.h"
#include "lagrangecoeff18.h"
#include "lagrangecoeff19.h"
#include "lagrangecoeff20.h"

#include "lagrangepolynomial.h"

using namespace std;

namespace hiflow
{

    /// Constructor

    template<class T>
    LagrangePolynomial<T>::LagrangePolynomial ( )
    {
        reinit ( );
    }

    /// reinit

    template<class T>
    void LagrangePolynomial<T>::reinit ( )
    {
        // Initialize pre computed polynomials

        int nlag = 20;

        _lag_coeff.reserve ( nlag );

        _lag_coeff.push_back ( &_lagrange_coeff1 );
        _lag_coeff.push_back ( &_lagrange_coeff2 );
        _lag_coeff.push_back ( &_lagrange_coeff3 );
        _lag_coeff.push_back ( &_lagrange_coeff4 );
        _lag_coeff.push_back ( &_lagrange_coeff5 );
        _lag_coeff.push_back ( &_lagrange_coeff6 );
        _lag_coeff.push_back ( &_lagrange_coeff7 );
        _lag_coeff.push_back ( &_lagrange_coeff8 );
        _lag_coeff.push_back ( &_lagrange_coeff9 );
        _lag_coeff.push_back ( &_lagrange_coeff10 );
        _lag_coeff.push_back ( &_lagrange_coeff11 );
        _lag_coeff.push_back ( &_lagrange_coeff12 );
        _lag_coeff.push_back ( &_lagrange_coeff13 );
        _lag_coeff.push_back ( &_lagrange_coeff14 );
        _lag_coeff.push_back ( &_lagrange_coeff15 );
        _lag_coeff.push_back ( &_lagrange_coeff16 );
        _lag_coeff.push_back ( &_lagrange_coeff17 );
        _lag_coeff.push_back ( &_lagrange_coeff18 );
        _lag_coeff.push_back ( &_lagrange_coeff19 );
        _lag_coeff.push_back ( &_lagrange_coeff20 );

    }

    /// poly ...

    template<class T>
    T LagrangePolynomial<T>::poly ( int deg, int i, const T& x ) const
    {
        assert ( deg > 0 );
        assert ( ( i >= 0 ) && ( i <= deg ) );
        // assert((x>=0.) && (x<=1. ));

        if ( deg <= _lag_coeff.size ( ) )
            return _horner ( deg, _lag_coeff[deg - 1]->poly ( i ), x );
        else
            return _lagrange_polynomial ( deg, i, x );
    }

    template<class T>
    T LagrangePolynomial<T>::poly_x ( int deg, int i, const T& x ) const
    {
        assert ( deg > 0 );
        assert ( ( i >= 0 ) && ( i <= deg ) );
        // assert((x>=0.) && (x<=1. ));

        if ( deg <= _lag_coeff.size ( ) )
            return _horner ( deg, _lag_coeff[deg - 1]->poly_x ( i ), x );
        else
            return _lagrange_polynomial_x ( deg, i, x );
    }

    template<class T>
    T LagrangePolynomial<T>::poly_xx ( int deg, int i, const T& x ) const
    {
        assert ( deg > 0 );
        assert ( ( i >= 0 ) && ( i <= deg ) );
        // assert((x>=0.) && (x<=1. ));

        if ( deg <= _lag_coeff.size ( ) )
            return _horner ( deg, _lag_coeff[deg - 1]->poly_xx ( i ), x );
        else
            return _lagrange_polynomial_xx ( deg, i, x );
    }

    /// template instanciation

    template class LagrangePolynomial<double>;
} // namespace hiflow
