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

/// @author Chandramowli Subramanian

#ifndef HIFLOW_LINEARALGEBRA_LA_DESCRIPTOR_H_
#    define HIFLOW_LINEARALGEBRA_LA_DESCRIPTOR_H_

#    include "linear_algebra/coupled_matrix.h"
#    include "linear_algebra/coupled_vector.h"
#    include "polynomial_chaos/pc_galerkin_matrix.h"
#    include "polynomial_chaos/pc_galerkin_vector.h"
#    include "linear_algebra/hypre_matrix.h"
#    include "linear_algebra/hypre_vector.h"
#    include "linear_algebra/pce_vector.h"
#    include "linear_algebra/pce_matrix.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Linear algebra descriptors to define matrix/vector pairs.
        ///
        /// Here you can add your own descriptor

        /// Generic Descriptor

        template<class Mat, class Vec, class Dat>
        class LADescriptor
        {
          public:
            typedef Mat MatrixType;
            typedef Vec VectorType;
            typedef Dat DataType;
        };

        /// Coupled matrices and vectors in double precision.

        class LADescriptorCoupledD
        {
          public:
            typedef CoupledMatrix<double> MatrixType;
            typedef CoupledVector<double> VectorType;
            typedef double DataType;
        };

        /// Coupled matrices and vectors in single precision.

        class LADescriptorCoupledS
        {
          public:
            typedef CoupledMatrix<float> MatrixType;
            typedef CoupledVector<float> VectorType;
            typedef float DataType;
        };

        /// Polynomial Chaos Galerkin matrices and vectors in double precision.

        class LADescriptorPolynomialChaosD
        {
          public:
            typedef polynomialchaos::PCGalerkinMatrix<double> MatrixType;
            typedef polynomialchaos::PCGalerkinVector<double> VectorType;
            typedef double DataType;
        };

        /// HYPRE matrices and vectors in double precision.

        class LADescriptorHypreD
        {
          public:
            typedef HypreMatrix<double> MatrixType;
            typedef HypreVector<double> VectorType;
            typedef double DataType;
        };

        /// Polynomial Chaos Expansion matrices and vector in double precision

        class LADescriptorPolynomialChaosExpansionD
        {
          public:
            typedef PCEMatrix<double> MatrixType;
            typedef PCEVector<double> VectorType;
            typedef double DataType;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_LA_DESCRIPTOR_H_
