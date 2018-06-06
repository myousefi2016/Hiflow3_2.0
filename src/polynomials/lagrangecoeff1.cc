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

#include "lagrangecoeff.h"
#include "lagrangecoeff1.h"

namespace hiflow
{

    static double _poly_lagrange_coeff1[] = {
                                             1.000000000000000000000000000000e+00,
                                             -1.000000000000000000000000000000e+00,

                                             0.000000000000000000000000000000e+00,
                                             1.000000000000000000000000000000e+00
    };

    static double _poly_x_lagrange_coeff1[] = {
                                               -1.000000000000000000000000000000e+00,
                                               0.000000000000000000000000000000e+00,

                                               1.000000000000000000000000000000e+00,
                                               0.000000000000000000000000000000e+00
    };

    static double _poly_xx_lagrange_coeff1[] = {
                                                0.000000000000000000000000000000e+00,
                                                0.000000000000000000000000000000e+00,

                                                0.000000000000000000000000000000e+00,
                                                0.000000000000000000000000000000e+00
    };

    LagrangeCoeff _lagrange_coeff1 ( 1,
                                     _poly_lagrange_coeff1,
                                     _poly_x_lagrange_coeff1,
                                     _poly_xx_lagrange_coeff1 );

} // namespace hiflow
