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

#include <algorithm>
#include <cassert>

#include "lagrangecoeff.h"
#include "polycoeff.h"

using namespace std;

namespace hiflow
{

    /// Constructor

    LagrangeCoeff::LagrangeCoeff ( int deg,
                                   const double* p,
                                   const double* p_x,
                                   const double* p_xx )
    {
        assert ( deg > 0 );

        // Allocate

        PolyCoeff<double> dummy ( deg );

        _poly .reserve ( deg + 1 );
        _poly .insert ( _poly .begin ( ), deg + 1, dummy );
        _poly_x .reserve ( deg + 1 );
        _poly_x .insert ( _poly_x .begin ( ), deg + 1, dummy );
        _poly_xx.reserve ( deg + 1 );
        _poly_xx.insert ( _poly_xx.begin ( ), deg + 1, dummy );

        assert ( static_cast < int > ( _poly .size ( ) ) == deg + 1 );
        assert ( static_cast < int > ( _poly_x .size ( ) ) == deg + 1 );
        assert ( static_cast < int > ( _poly_xx.size ( ) ) == deg + 1 );

        // Initialize

        const double* it_p = p;
        const double* it_p_x = p_x;
        const double* it_p_xx = p_xx;

        const int deg_1 = deg + 1;

        for ( int k = 0; k <= deg; ++k )
        {
            copy ( it_p, it_p + deg_1, poly ( k ).begin ( ) );
            it_p += deg_1;
            copy ( it_p_x, it_p_x + deg_1, poly_x ( k ).begin ( ) );
            it_p_x += deg_1;
            copy ( it_p_xx, it_p_xx + deg_1, poly_xx ( k ).begin ( ) );
            it_p_xx += deg_1;
        }

    }

} // namespace hiflow
