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

/// \author Code generated with program quadrature_solin_generator.

// This code was generated with quadrature_solin_generator by Staffan Ronnas.
// The quadrature data comes from the book
// P. Solin, et. al., "Higher-Order Finite Element Methods"
// Michael Schick developed the original class on which the generated code is based.

#ifndef HIFLOW_QUADRATURE_QUADRATUREGAUSSTRIANGLE_H
#    define HIFLOW_QUADRATURE_QUADRATUREGAUSSTRIANGLE_H

#    include "quadraturetype.h"

namespace hiflow
{

    template<class DataType>
    class QuadratureGaussTriangle : public QuadratureType<DataType>
    {
      public:
        QuadratureGaussTriangle ( );
        QuadratureGaussTriangle<DataType>* clone ( ) const;

        using QuadratureType<DataType>::x_;
        using QuadratureType<DataType>::y_;
        using QuadratureType<DataType>::z_;
        using QuadratureType<DataType>::w_;
        using QuadratureType<DataType>::index_field_;
        using QuadratureType<DataType>::order_size_map_;
    };
} // namespace hiflow
#endif
