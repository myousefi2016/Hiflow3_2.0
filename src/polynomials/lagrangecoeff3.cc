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

// $Id: lagrangecoeff3.cc,v 1.1 2000/03/11 14:13:41 heuvelin Exp $

/// \author Martin Baumann, Michael Schick

#include "lagrangecoeff.h"
#include "lagrangecoeff3.h"

namespace hiflow
{

    static double _poly_lagrange_coeff3[] = {
                                             1.000000000000000000000000000000e+00,
                                             -5.500000000000000000000000000000e+00,
                                             9.000000000000000000000000000000e+00,
                                             -4.500000000000000000000000000000e+00,

                                             0.000000000000000000000000000000e-01,
                                             9.000000000000000000000000000000e+00,
                                             -2.250000000000000000000000000000e+01,
                                             1.350000000000000000000000000000e+01,

                                             0.000000000000000000000000000000e-01,
                                             -4.500000000000000000000000000000e+00,
                                             1.800000000000000000000000000000e+01,
                                             -1.350000000000000000000000000000e+01,

                                             0.000000000000000000000000000000e-01,
                                             1.000000000000000000000000000000e+00,
                                             -4.500000000000000000000000000000e+00,
                                             4.500000000000000000000000000000e+00

    };

    static double _poly_x_lagrange_coeff3[] = {
                                               -5.500000000000000000000000000000e+00,
                                               1.800000000000000000000000000000e+01,
                                               -1.350000000000000000000000000000e+01,
                                               0.000000000000000000000000000000e-01,

                                               9.000000000000000000000000000000e+00,
                                               -4.500000000000000000000000000000e+01,
                                               4.050000000000000000000000000000e+01,
                                               0.000000000000000000000000000000e-01,

                                               -4.500000000000000000000000000000e+00,
                                               3.600000000000000000000000000000e+01,
                                               -4.050000000000000000000000000000e+01,
                                               0.000000000000000000000000000000e-01,

                                               1.000000000000000000000000000000e+00,
                                               -9.000000000000000000000000000000e+00,
                                               1.350000000000000000000000000000e+01,
                                               0.000000000000000000000000000000e-01

    };

    static double _poly_xx_lagrange_coeff3[] = {
                                                1.800000000000000000000000000000e+01,
                                                -2.700000000000000000000000000000e+01,
                                                0.000000000000000000000000000000e-01,
                                                0.000000000000000000000000000000e-01,

                                                -4.500000000000000000000000000000e+01,
                                                8.100000000000000000000000000000e+01,
                                                0.000000000000000000000000000000e-01,
                                                0.000000000000000000000000000000e-01,

                                                3.600000000000000000000000000000e+01,
                                                -8.100000000000000000000000000000e+01,
                                                0.000000000000000000000000000000e-01,
                                                0.000000000000000000000000000000e-01,

                                                -9.000000000000000000000000000000e+00,
                                                2.700000000000000000000000000000e+01,
                                                0.000000000000000000000000000000e-01,
                                                0.000000000000000000000000000000e-01

    };

    LagrangeCoeff _lagrange_coeff3 ( 3,
                                     _poly_lagrange_coeff3,
                                     _poly_x_lagrange_coeff3,
                                     _poly_xx_lagrange_coeff3 );

} // namespace hiflow
