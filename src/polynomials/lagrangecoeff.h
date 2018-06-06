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

/// \author Martin Baumann, Michael Schick

#ifndef __POLYNOMIALS_LAGRANGECOEFF_H_
#    define __POLYNOMIALS_LAGRANGECOEFF_H_

#    include <cassert>
#    include <vector>

#    include "polynomials/polycoeff.h"

namespace hiflow
{

    class LagrangeCoeff
    {
      protected:

        std::vector<PolyCoeff<double> > _poly;
        std::vector<PolyCoeff<double> > _poly_x;
        std::vector<PolyCoeff<double> > _poly_xx;

      public:

        LagrangeCoeff ( int, const double*, const double*, const double* );

        inline const PolyCoeff<double>& poly ( int ) const;
        inline PolyCoeff<double>& poly ( int );

        inline const PolyCoeff<double>& poly_x ( int ) const;
        inline PolyCoeff<double>& poly_x ( int );

        inline const PolyCoeff<double>& poly_xx ( int ) const;
        inline PolyCoeff<double>& poly_xx ( int );

    };

    /* --------------------------------------------------------------
    /                                                               /
    /  INLINE FUNCTIONS for LagrangeCoeff                           /
    /                                                               /
    /  ----------------------------------------------------------- */

    /// poly

    inline
    const PolyCoeff<double>& LagrangeCoeff::poly ( int k ) const
    {
        assert ( ( k >= 0 ) && ( k<static_cast < int > ( _poly.size ( ) ) ) );
        return _poly[k];
    }

    inline
    PolyCoeff<double>& LagrangeCoeff::poly ( int k )
    {
        assert ( ( k >= 0 ) && ( k<static_cast < int > ( _poly.size ( ) ) ) );
        return _poly[k];
    }

    /// poly_x

    inline
    const PolyCoeff<double>& LagrangeCoeff::poly_x ( int k ) const
    {
        assert ( ( k >= 0 ) && ( k<static_cast < int > ( _poly_x.size ( ) ) ) );
        return _poly_x[k];
    }

    inline
    PolyCoeff<double>& LagrangeCoeff::poly_x ( int k )
    {
        assert ( ( k >= 0 ) && ( k<static_cast < int > ( _poly_x.size ( ) ) ) );
        return _poly_x[k];
    }

    /// poly_xx

    inline
    const PolyCoeff<double>& LagrangeCoeff::poly_xx ( int k ) const
    {
        assert ( ( k >= 0 ) && ( k<static_cast < int > ( _poly_xx.size ( ) ) ) );
        return _poly_xx[k];
    }

    inline
    PolyCoeff<double>& LagrangeCoeff::poly_xx ( int k )
    {
        assert ( ( k >= 0 ) && ( k<static_cast < int > ( _poly_xx.size ( ) ) ) );
        return _poly_xx[k];
    }

} // namespace hiflow

#endif
