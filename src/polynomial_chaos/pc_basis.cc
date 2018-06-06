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

/// \author Michael Schick

#include "polynomial_chaos/pc_basis.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace hiflow
{
    namespace polynomialchaos
    {

        PCBasis::PCBasis ( )
        {
            N_ = -1;
            No_ = -1;
            Dim_ = -1;
            q_ = -1.0;
            t_ = -1;

            // Precomputed quadrature is exact up to polynomial order 41
            x_qpts_precomp_.resize ( ndist_ );
            w_qpts_precomp_.resize ( ndist_ );

            double x_qpts_tmp[] = {
                                   -9.937521706203894523e-01,
                                   -9.672268385663064238e-01,
                                   -9.200993341504007939e-01,
                                   -8.533633645833176296e-01,
                                   -7.684399634756780006e-01,
                                   -6.671388041974124494e-01,
                                   -5.516188358872193831e-01,
                                   -4.243421202074389997e-01,
                                   -2.880213168024010617e-01,
                                   -1.455618541608955929e-01,
                                   1.544987264365841335e-17,
                                   1.455618541608950933e-01,
                                   2.880213168024010617e-01,
                                   4.243421202074387222e-01,
                                   5.516188358872200492e-01,
                                   6.671388041974123384e-01,
                                   7.684399634756776676e-01,
                                   8.533633645833172965e-01,
                                   9.200993341504010159e-01,
                                   9.672268385663062018e-01,
                                   9.937521706203894523e-01
            };

            double x_qpts_tmp_g[] = {
                                     -7.849382895113821590e+00,
                                     -6.751444718717457327e+00,
                                     -5.829382007304467095e+00,
                                     -4.994963944782024434e+00,
                                     -4.214343981688423391e+00,
                                     -3.469846690475375084e+00,
                                     -2.750592981052374153e+00,
                                     -2.049102468257162801e+00,
                                     -1.359765823211230407e+00,
                                     -6.780456924406442765e-01,
                                     -1.474106067249749284e-16,
                                     6.780456924406442765e-01,
                                     1.359765823211230629e+00,
                                     2.049102468257163689e+00,
                                     2.750592981052373709e+00,
                                     3.469846690475375972e+00,
                                     4.214343981688421614e+00,
                                     4.994963944782024434e+00,
                                     5.829382007304468871e+00,
                                     6.751444718717459992e+00,
                                     7.849382895113824254e+00
            };

            double w_qpts_tmp[] = {
                                   8.008614128887271352e-03,
                                   1.847689488542610808e-02,
                                   2.856721271342867532e-02,
                                   3.805005681418983238e-02,
                                   4.672221172801710454e-02,
                                   5.439864958357445990e-02,
                                   6.091570802686405162e-02,
                                   6.613446931666877582e-02,
                                   6.994369739553660259e-02,
                                   7.226220199498513408e-02,
                                   7.304056682484481866e-02,
                                   7.226220199498505081e-02,
                                   6.994369739553646381e-02,
                                   6.613446931666837336e-02,
                                   6.091570802686474551e-02,
                                   5.439864958357445990e-02,
                                   4.672221172801706290e-02,
                                   3.805005681418925645e-02,
                                   2.856721271342847063e-02,
                                   1.847689488542662156e-02,
                                   8.008614128887031960e-03
            };

            double w_qpts_tmp_g[] = {
                                     2.098991219565668087e-14,
                                     4.975368604121643035e-11,
                                     1.450661284493110865e-08,
                                     1.225354836148250053e-06,
                                     4.219234742551674646e-05,
                                     7.080477954815380141e-04,
                                     6.439697051408762966e-03,
                                     3.395272978654294976e-02,
                                     1.083922856264193380e-01,
                                     2.153337156950596021e-01,
                                     2.702601835728762336e-01,
                                     2.153337156950596021e-01,
                                     1.083922856264195184e-01,
                                     3.395272978654288731e-02,
                                     6.439697051408775109e-03,
                                     7.080477954815368215e-04,
                                     4.219234742551678712e-05,
                                     1.225354836148263605e-06,
                                     1.450661284493089523e-08,
                                     4.975368604121606200e-11,
                                     2.098991219565668087e-14
            };

            x_qpts_precomp_[UNIFORM].reserve ( 21 );
            x_qpts_precomp_[UNIFORM].insert ( x_qpts_precomp_[UNIFORM].begin ( ), x_qpts_tmp, x_qpts_tmp + 21 );
            w_qpts_precomp_[UNIFORM].reserve ( 21 );
            w_qpts_precomp_[UNIFORM].insert ( w_qpts_precomp_[UNIFORM].begin ( ), w_qpts_tmp, w_qpts_tmp + 21 );

            x_qpts_precomp_[GAUSSIAN].reserve ( 21 );
            x_qpts_precomp_[GAUSSIAN].insert ( x_qpts_precomp_[GAUSSIAN].begin ( ), x_qpts_tmp_g, x_qpts_tmp_g + 21 );
            w_qpts_precomp_[GAUSSIAN].reserve ( 21 );
            w_qpts_precomp_[GAUSSIAN].insert ( w_qpts_precomp_[GAUSSIAN].begin ( ), w_qpts_tmp_g, w_qpts_tmp_g + 21 );
        }

        int PCBasis::CalculateDimension ( int No, int N ) const
        {
            // Calculate total dimension
            long int fak_No = 1;

            for ( int i = 1; i <= No; ++i )
                fak_No *= i;

            long int tmp = 1;
            for ( int i = N + 1; i <= No + N; ++i )
                tmp *= i;

            return tmp / fak_No;
        }

        void PCBasis::Init ( int N, int No, std::vector<Distribution> const& distributions )
        {
            assert ( static_cast < int > ( distributions.size ( ) ) == N );

            N_ = N;
            No_ = No;

            q_ = 1.0;
            t_ = N_;

            SubDim_.resize ( No_ + 1 );
            SubDim_[0] = 1;

            for ( int p = 1; p<static_cast < int > ( SubDim_.size ( ) ); ++p )
                SubDim_[p] = CalculateDimension ( p, N_ );

            // Calculate total dimension
            int fak_N = 1;
            int fak_No = 1;
            int fak_NplusNo = 1;

            for ( int i = 1; i <= N_; ++i )
                fak_N *= i;
            for ( int i = 1; i <= No_; ++i )
                fak_No *= i;
            for ( int i = 1; i <= N_ + No_; ++i )
                fak_NplusNo *= i;

            int tmp = 1;
            for ( int i = N_ + 1; i <= N_ + No_; ++i )
                tmp *= i;

            Dim_ = tmp / fak_No;

            distributions_ = distributions;

            quadrature_computed_.resize ( ndist_, false );

            x_qpts_.resize ( ndist_ );
            w_qpts_.resize ( ndist_ );

            ComputeQuadrature ( );
            ComputeMultiIndices ( );

            Precomputed3rdOrderIntegrals ( );
        }

        void PCBasis::Init ( int N, int No, std::vector<Distribution> const& distributions, double q, int t )
        {
            assert ( static_cast < int > ( distributions.size ( ) ) == N );

            N_ = N;
            No_ = No;
            q_ = q;
            t_ = t;

            SubDim_.resize ( No_ + 1 );
            SubDim_[0] = 1;

            for ( int p = 1; p<static_cast < int > ( SubDim_.size ( ) ); ++p )
                SubDim_[p] = CalculateDimension ( p, N_ );

            distributions_ = distributions;

            quadrature_computed_.resize ( ndist_, false );

            x_qpts_.resize ( ndist_ );
            w_qpts_.resize ( ndist_ );

            ComputeQuadrature ( );
            ComputeMultiIndices ( );

            Dim_ = alpha_.size ( );

            Precomputed3rdOrderIntegrals ( );
        }

        void PCBasis::ComputeQuadrature ( )
        {
#ifdef WITH_GAUSSQ
            for ( int i = 0; i < N_; ++i )
            {
                if ( distributions_[i] == UNIFORM ) ComputeGaussLegendreQuadrature ( );
                if ( distributions_[i] == GAUSSIAN ) ComputeGaussHermiteQuadrature ( );
                if ( distributions_[i] == NOT_SET ) assert ( 0 );
            }
#else
            for ( int i = 0; i < N_; ++i )
            {
                w_qpts_[distributions_[i]] = w_qpts_precomp_[distributions_[i]];
                x_qpts_[distributions_[i]] = x_qpts_precomp_[distributions_[i]];
                quadrature_computed_[distributions_[i]] = true;
            }
#endif
        }

        void PCBasis::ComputeGaussLegendreQuadrature ( )
        {
            if ( !quadrature_computed_[UNIFORM] )
            {
                int n;
                // Exact up to 4th order tensor
                if ( ( 4 * No_ + 1 ) % 2 == 0 ) n = ( 4 * No_ + 1 ) / 2;
                else n = ( 4 * No_ + 1 ) / 2 + 1;

                assert ( n > 0 );

                double* t = new double [n];
                double* w = new double [n];

#ifdef WITH_GAUSSQ
                double alpha = 0.0;
                double beta = 0.0;
                int kpts = 0;
                double* endpts = new double [2];
                double* b = new double [n];

                int kind = 1;

                gaussq_ ( &kind, &n, &alpha, &beta, &kpts, endpts, b, t, w );
#endif
                w_qpts_[UNIFORM].resize ( n );
                x_qpts_[UNIFORM].resize ( n );

                for ( int i = 0; i < n; ++i )
                {
                    w_qpts_[UNIFORM][i] = 0.5 * w[i];
                    x_qpts_[UNIFORM][i] = t[i];
                }

#ifdef WITH_GAUSSQ
                delete[] endpts;
                delete[] b;
#endif
                delete[] t;
                delete[] w;

                quadrature_computed_[UNIFORM] = true;
            }
        }

        void PCBasis::ComputeGaussHermiteQuadrature ( )
        {
            if ( !quadrature_computed_[GAUSSIAN] )
            {
                int n;
                // Exact up to 4th order tensor
                if ( ( 4 * No_ + 1 ) % 2 == 0 ) n = ( 4 * No_ + 1 ) / 2;
                else n = ( 4 * No_ + 1 ) / 2 + 1;

                assert ( n > 0 );

                double* t = new double [n];
                double* w = new double [n];

#ifdef WITH_GAUSSQ
                double alpha = 0.0;
                double beta = 0.0;
                int kpts = 0;
                double* endpts = new double [2];
                double* b = new double [n];

                int kind = 4;

                gaussq_ ( &kind, &n, &alpha, &beta, &kpts, endpts, b, t, w );
#endif

                w_qpts_[GAUSSIAN].resize ( n );
                x_qpts_[GAUSSIAN].resize ( n );
                for ( int i = 0; i < n; ++i )
                {
                    w_qpts_[GAUSSIAN][i] = ( 1.0 / std::sqrt ( M_PI ) ) * w[i];
                    x_qpts_[GAUSSIAN][i] = std::sqrt ( 2.0 ) * t[i];
                }

#ifdef WITH_GAUSSQ
                delete[] endpts;
                delete[] b;
#endif
                delete[] t;
                delete[] w;

                quadrature_computed_[GAUSSIAN] = true;
            }
        }

        double PCBasis::LegendreP ( int degree, double x ) const
        {
            if ( degree == 0 )
            {
                return 1.0;
            }
            else if ( degree == 1 )
            {
                return x;
            }
            else
            {
                return (( 2.0 * degree - 1.0 ) * x * LegendreP ( degree - 1, x ) / degree - ( degree - 1.0 ) * LegendreP ( degree - 2, x ) / degree );
            }
        }

        double PCBasis::HermiteP ( int degree, double x ) const
        {
            if ( degree == 0 )
            {
                return 1.0;
            }
            else if ( degree == 1 )
            {
                return x;
            }
            else
            {
                return (x * HermiteP ( degree - 1, x )-( degree - 1.0 ) * HermiteP ( degree - 2, x ) );
            }
        }

        double PCBasis::poly ( Distribution dist, int degree, double x ) const
        {
            switch ( dist )
            {
                case UNIFORM:
                    return LegendreP ( degree, x );
                case GAUSSIAN:
                    return HermiteP ( degree, x );
                default:
                {
                    assert ( 0 );
                    return 0.0;
                }
            };
        }

        void PCBasis::Precomputed3rdOrderIntegrals ( )
        {
            third_order_integrals_.resize ( ndist_ );
            for ( int n = 0; n < ndist_; ++n )
            {
                third_order_integrals_[n].resize ( No_ + 1 );
                for ( int i = 0; i <= No_; ++i )
                {
                    third_order_integrals_[n][i].resize ( No_ + 1 );
                    for ( int j = 0; j <= No_; ++j )
                    {
                        third_order_integrals_[n][i][j].resize ( No_ + 1, 0.0 );
                    }
                }
            }

            for ( int n = 0; n < ndist_; ++n )
            {
                for ( int i = 0; i <= No_; ++i )
                {
                    for ( int j = 0; j <= No_; ++j )
                    {
                        for ( int k = 0; k <= No_; ++k )
                        {
                            if ( k <= i + j && i <= j + k && j <= i + k )
                            {
                                for ( int qpt = 0; qpt<static_cast < int > ( x_qpts_[n].size ( ) ); ++qpt )
                                {
                                    third_order_integrals_[n][i][j][k] += w_qpts_[n][qpt]
                                            * poly ( ( Distribution ) n, i, x_qpts_[n][qpt] )
                                            * poly ( ( Distribution ) n, j, x_qpts_[n][qpt] )
                                            * poly ( ( Distribution ) n, k, x_qpts_[n][qpt] );
                                }
                            }
                        }
                    }
                }
            }
        }

        double PCBasis::ComputeIntegral3rdOrder ( int i, int j, int k ) const
        {
            double integral = 1.0;
            for ( int n = 0; n < N_; ++n )
                integral *= third_order_integrals_[distributions_[n]][alpha_[i][n]][alpha_[j][n]][alpha_[k][n]];
            return integral;
        }

        double PCBasis::ComputeIntegral ( std::vector<int> const& global_indices ) const
        {
            int size = global_indices.size ( );

            // 1st filter polynomial degrees according to distributions
            std::vector<std::vector<int> > degrees ( N_ );
            for ( int k = 0; k < size; ++k )
            {
                for ( int i = 0; i < N_; ++i )
                {
                    degrees[i].push_back ( alpha_[global_indices[k]][i] );
                }
            }

            // 2nd compute integral with appropriate quadrature and distribution
            std::vector<double> loc_integral ( N_, 0.0 );
            for ( int j = 0; j < N_; ++j )
            {
                for ( int qpt = 0; qpt<static_cast < int > ( x_qpts_[distributions_[j]].size ( ) ); ++qpt )
                {
                    double f = 1.0;
                    for ( int i = 0; i<static_cast < int > ( degrees[j].size ( ) ); ++i )
                        f *= poly ( distributions_[j], degrees[j][i], x_qpts_[distributions_[j]][qpt] );
                    loc_integral[j] += w_qpts_[distributions_[j]][qpt] * f;
                }
            }

            // 3rd compute total integral as product of local ones
            double integral = 1.0;
            for ( int i = 0; i < N_; ++i )
            {
                integral *= loc_integral[i];
            }

            return integral;
        }

        void PCBasis::ComputeMultiIndices ( )
        {
            if ( q_ == 0 && t_ == 0 )
            {
                alpha_.resize ( 1 );
                alpha_[0].resize ( N_, 0 );
                return;
            }

            int minDim_ = CalculateDimension ( No_, N_ );
            bool notfinished = true;
            std::vector<std::vector<int> > beta_;

            beta_.resize ( minDim_ );
            if ( q_ != 0 )
            {
                alpha_.resize ( minDim_ );
            }
            if ( q_ == 0 )
            {
                alpha_.resize ( N_ + 1 );
            }

            for ( int i = 0; i < minDim_; ++i )
            {
                beta_[i].resize ( N_, 0 );
                if ( q_ != 0 )
                {
                    alpha_[i].resize ( N_, 0 );
                }
                if ( q_ == 0 && i < N_ + 1 )
                {
                    alpha_[i].resize ( N_, 0 );
                }
            }

            if ( No_ > 0 )
            {
                for ( int i = 1; i <= N_; ++i )
                {
                    beta_[i][i - 1] = 1;
                    alpha_[i][i - 1] = 1;
                }
            }

            if ( No_ > 1 )
            {
                int P = N_;
                std::vector<std::vector<int> > p ( N_ );

                for ( int i = 0; i < N_; ++i ) p[i].resize ( No_, 0 );

                for ( int i = 0; i < N_; ++i ) p[i][0] = 1;

                for ( int k = 1; notfinished; ++k )
                {

                    if ( q_ == 0 && k == No_ - 1 )
                    {
                        notfinished = false;
                    }

                    if ( k >= No_ )
                    {
                        notfinished = false;
                        int lastDim_ = beta_.size ( );
                        int newDim_ = CalculateDimension ( k + 1, N_ );
                        beta_.resize ( newDim_ );
                        for ( int i = lastDim_; i < newDim_; ++i )
                        {
                            beta_[i].resize ( N_, 0 );
                        }
                        for ( int i = 0; i < N_; ++i )
                        {
                            p[i].resize ( k + 1, 0 );
                        }
                    }
                    int L = P;

                    for ( int i = 0; i < N_; ++i )
                    {
                        int sum = 0;
                        for ( int m = i; m < N_; ++m ) sum += p[m][k - 1];
                        p[i][k] = sum;
                    }

                    for ( int j = 0; j < N_; ++j )
                    {
                        for ( int m = L - p[j][k] + 1; m <= L; ++m )
                        {
                            P++;
                            for ( int i = 0; i < N_; ++i )
                            {
                                beta_[P][i] = beta_[m][i];
                            }
                            beta_[P][j] = beta_[P][j] + 1;

                            if ( q_ == 0 )
                            {
                                int max = *std::max_element ( beta_[P].begin ( ), beta_[P].end ( ) );
                                if ( max <= 1 )
                                {
                                    int vectorsum = 0;
                                    for ( int i = 0; i < N_; ++i )
                                    {
                                        vectorsum += beta_[P][i];
                                    }
                                    if ( vectorsum <= t_ )
                                    {
                                        int size = alpha_.size ( );
                                        alpha_.resize ( size + 1 );
                                        alpha_[size].resize ( N_, 0 );
                                        for ( int i = 0; i < N_; ++i )
                                        {
                                            alpha_[size][i] = beta_[P][i];
                                        }
                                    }
                                }
                            }

                            if ( k < No_ && q_ != 0 )
                            {
                                for ( int i = 0; i < N_; ++i )
                                {
                                    alpha_[P][i] = beta_[P][i];
                                }
                            }

                            if ( k >= No_ )
                            {
                                double vectorsum = 0;
                                for ( int i = 0; i < N_; ++i )
                                {
                                    vectorsum += pow ( beta_[P][i], q_ );
                                }
                                if ( vectorsum <= No_ )
                                {
                                    notfinished = true;
                                    int size = alpha_.size ( );
                                    alpha_.resize ( size + 1 );
                                    alpha_[size].resize ( N_, 0 );
                                    for ( int i = 0; i < N_; ++i )
                                    {
                                        alpha_[size][i] = beta_[P][i];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        double PCBasis::poly ( int index, std::vector<double> const& x ) const
        {
            assert ( static_cast < int > ( x.size ( ) ) == N_ );

            double res = 1.0;
            std::vector<int> loc_index = alpha_[index];

            for ( int j = 0; j < N_; ++j )
            {
                res *= poly ( distributions_[j], loc_index[j], x[j] );
            }

            return res;
        }

    }
}
