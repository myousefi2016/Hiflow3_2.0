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

#ifndef __QUADRATURE_QGAUSSTETRAHEDRON_H_
#    define __QUADRATURE_QGAUSSTETRAHEDRON_H_

#    include "quadraturetype.h"

namespace hiflow
{

    ///
    /// \class QuadratureGaussTetrahedron qgausstetrahedron.h
    /// \brief Derived class of all gaussian quadratures on tetrahedrons
    ///        which is templated by a DataType (e.g. double)
    /// \author Michael Schick
    ///

    template<class DataType>
    class QuadratureGaussTetrahedron : public QuadratureType<DataType>
    {
      public:

        /// Constructor which sets up all necessary information about the
        /// Gaussian quadrature on a tetrahedron for all implemented sizes
        QuadratureGaussTetrahedron ( );

        QuadratureGaussTetrahedron<DataType>* clone ( ) const;

        using QuadratureType<DataType>::x_;
        using QuadratureType<DataType>::y_;
        using QuadratureType<DataType>::z_;
        using QuadratureType<DataType>::w_;
        using QuadratureType<DataType>::index_field_;
    };

} // namespace hiflow

#endif
