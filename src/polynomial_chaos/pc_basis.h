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

#ifndef HIFLOW_POLYNOMIALCHAOS_PC_BASIS_H_
#    define HIFLOW_POLYNOMIALCHAOS_PC_BASIS_H_

/// \file pc_basis.h
/// \brief Basic class for working with Polynomial Chaos.
///
/// \author Michael Schick

#    include <vector>
#    include <map>
#    include "config.h"

#    ifdef WITH_GAUSSQ
extern "C"
{
    void gaussq_ ( int* kind, int* n, double* alpha, double* beta, int* kpts, double* endpts, double* b, double* t, double* w );
}
#    endif

namespace hiflow
{
    namespace polynomialchaos
    {

        class PCBasis
        {
          public:
            /// Implemented probability distributions

            enum Distribution
            {
                NOT_SET = -1,
                UNIFORM,
                GAUSSIAN
            };

            /// Constructor
            PCBasis ( );
            /// Destructor

            ~PCBasis ( )
            {
            }

            // Initialize with random space dimension N, polynomial degree No,
            /// vector containing probability distributions of input (size N)
            void Init ( int N, int No, std::vector<Distribution> const& distributions );

            /// Initialize with random space dimension N, polynomial degree No,
            /// vector containing probability distributions of input (size N), exponent for hyperbolic index selection and maximal number of nonzero
            /// entries in multi index
            void Init ( int N, int No, std::vector<Distribution> const& distributions, double q, int t );

            /// Compute integral of chaos polynomials according to their global indices
            double ComputeIntegral ( std::vector<int> const& global_indices ) const;

            /// Compute third order integral (faster than the general approach)
            double ComputeIntegral3rdOrder ( int i, int j, int k ) const;

            /// Returns value of a chaos polynomial associated with a certain distribution
            double poly ( Distribution dist, int degree, double x ) const;

            /// Returns value of tensor-polynomial with global index at position x
            double poly ( int index, std::vector<double> const& x ) const;

            /// Returns total dimension of basis

            int Dim ( ) const
            {
                return Dim_;
            }

            /// Returns total dimension of basis w.r.t. a certain total polynomial degree

            int SubDim ( int p ) const
            {
                return SubDim_[p];
            }

            /// Calculate the Dimension of PC space for given polynomial degree and number of rvs
            int CalculateDimension ( int No, int N ) const;

            /// Returns dimension of probability space

            int N ( ) const
            {
                return N_;
            }

            /// Returns Polynomial Degree

            int No ( ) const
            {
                return No_;
            }

            /// Global 2 Local (for multi-index)

            std::vector<int> G2L ( int glo_idx ) const
            {
                return alpha_[glo_idx];
            }

            /// Local 2 Global (for multi-index)

            int L2G ( std::vector<int>& idx ) const
            {
                return l2g_.at ( idx );
            }

            /// Returns vector of used probability distributions

            Distribution const& Distributions ( int n ) const
            {
                return distributions_[n];
            }

          protected:
            /// Precompute one-dimensional third order integrals
            void Precomputed3rdOrderIntegrals ( );
            /// Dimension of probability space
            int N_;
            /// Polynomial degree
            int No_;
            /// Total dimension of PC expansion
            int Dim_;
            /// Exponent for hyperbolic index selection
            double q_;
            /// Maximal number of nonzeros multi index entries
            int t_;

            /// Total dimension of PC expansion w.r.t. all degrees less or equal to max. dimension
            std::vector<int> SubDim_;
            /// Number of implemented distributions
            static const int ndist_ = 2;

            /// Compute Multi-Indices for multidimensional random input
            void ComputeMultiIndices ( );
            /// Vector containing Multi-Indices
            std::vector<std::vector<int> > alpha_;
            /// Map local to global for multi-indices
            std::map<std::vector<int>, int> l2g_;

            /// Stored probability distributions
            std::vector<Distribution> distributions_;

            /// Check variable for computation of quadratures
            std::vector<bool> quadrature_computed_;

            /// Compute necessary quadrature points and weights
            void ComputeQuadrature ( );
            /// Compute GaussLegendre quadrature
            void ComputeGaussLegendreQuadrature ( );
            /// Compute GaussHermite quadrature
            void ComputeGaussHermiteQuadrature ( );

            /// Gauss quadrature nodes
            std::vector<std::vector<double> > x_qpts_;
            /// Gauss quadrature weights
            std::vector<std::vector<double> > w_qpts_;

            /// Precomputed quadrature nodes (if gaussq.f is not available)
            std::vector<std::vector<double> > x_qpts_precomp_;
            /// Precomputed quadrature weights (if gaussq.f is not available)
            std::vector<std::vector<double> > w_qpts_precomp_;

            /// Precomputed one-dimensional third order integrals
            std::vector<std::vector<std::vector<std::vector<double> > > > third_order_integrals_;

            /// Legendre polynomials
            double LegendreP ( int degree, double x ) const;
            /// Hermite polynomials
            double HermiteP ( int degree, double x ) const;
        };

    }
}

#endif
