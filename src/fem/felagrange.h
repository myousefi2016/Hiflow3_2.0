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

#ifndef _FEM_FELAGRANGE_H_
#    define _FEM_FELAGRANGE_H_

#    include "fem/fetype.h"
#    include "polynomials/lagrangepolynomial.h"

namespace hiflow
{
    namespace doffem
    {

        ///
        /// \class FELagrange felagrange.h
        /// \brief Lagrangian Finite Element
        /// \author Michael Schick<br>Martin Baumann
        ///

        template<class DataType>
        class FELagrange : public FEType<DataType>
        {
          public:

            /// Default Constructor
            FELagrange ( );

          protected:

            /// Lagrange polynomials which are used for evaluating shapefunctions
            LagrangePolynomial<double> lp_;
        };

    }
} // namespace hiflow
#endif
