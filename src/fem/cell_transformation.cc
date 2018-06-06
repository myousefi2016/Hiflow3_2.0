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

#include "common/vector_algebra.h"
#include "fem/cell_transformation.h"
#include <cmath>

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

namespace hiflow
{
    namespace doffem
    {

        /// \details Setting up the geometrical dimension

        template<class DataType>
        CellTransformation<DataType>::CellTransformation ( int gdim ) : gdim_ ( gdim )
        {
        }

        /// \details Given vector of coordinates on physical cell, the are
        ///          stored the protected member variable coord_vtx_

        template<class DataType>
        void CellTransformation<DataType>::reinit ( const Coord& coord_vtx )
        {
            coord_vtx_ = coord_vtx;
        }

        /// \details The Inverse of the mapping (reference cell to physical cell)
        ///          is computed in 2D via a Newton-scheme with Armijo updates. The
        ///          Problem, which is to be solved writes<br>
        ///          G([x_ref, y_ref]) := [x_phy - x(coord_ref) , y_phy - y(coord_ref) ]<br>
        ///          solve G([x_ref, y_ref]) = 0 via Newton-scheme
        /// \param[in] x_phy x coordinate on physical cell
        /// \param[in] y_phy y coordinate on physical cell
        /// \param[out] x_ref x coordinate on reference cell
        /// \param[out] y_ref y coordinate on reference cell

        template<class DataType>
        void CellTransformation<DataType>::inverse_newton_2d ( DataType x_phy, DataType y_phy,
                                                               DataType& x_ref, DataType& y_ref ) const
        {
            // Initialisation

            Vec<2, DataType> pt_phy;
            pt_phy[0] = x_phy;
            pt_phy[1] = y_phy;

            Vec<2, DataType> ref_k1, ref_k;

            ref_k1[0] = 0.55;
            ref_k1[1] = 0.55;

            ref_k[0] = 0.5;
            ref_k[1] = 0.5;

            const DataType tol_eps = 1.e3 * std::numeric_limits<DataType>::epsilon ( );

            // Some general parameters

            const int iter_max = 1000;
            int iter = 0;

            // Residual

            Coord coord_ref_k ( 2 );
            coord_ref_k[0] = ref_k[0];
            coord_ref_k[1] = ref_k[1];

            Vec<2, DataType> pt_k;
            pt_k[0] = this->x ( coord_ref_k );
            pt_k[1] = this->y ( coord_ref_k );

            DataType residual = norm ( pt_k - pt_phy );
            DataType abserr = norm ( ref_k1 - ref_k );

            // Newton

            // Jacobian Matrix (grad G)
            Mat<2, 2, DataType> G;
            // Inverse of the jacobian
            Mat<2, 2, DataType> B;

            while ( ( residual > 5. * std::numeric_limits<DataType>::epsilon ( ) ) && ( abserr > 5. * std::numeric_limits<DataType>::epsilon ( ) ) )
            {

                G ( 0, 0 ) = this->x_x ( coord_ref_k );
                G ( 0, 1 ) = this->x_y ( coord_ref_k );

                G ( 1, 0 ) = this->y_x ( coord_ref_k );
                G ( 1, 1 ) = this->y_y ( coord_ref_k );

#ifndef NDEBUG
                const DataType detG = det ( G );
                assert ( detG != 0. );
#endif
                inv ( G, B );

                // Armijo parameter

                const int iter_max_armijo = 500;
                int iter_armijo = 0;

                DataType residual_armijo = 2. * residual;

                DataType omega = 1.;

                const Vec<2, DataType> update_vec = B * ( pt_phy - pt_k );

                // Start Armijo

                while ( ( iter_armijo <= iter_max_armijo )
                        && ( residual_armijo > residual ) )
                {

                    ref_k1 = ref_k + omega * update_vec;

                    if ( ( ref_k1[0] >= -tol_eps ) && ( ref_k1[0] <= tol_eps ) ) ref_k1[0] = 0.;
                    if ( ( ref_k1[1] >= -tol_eps ) && ( ref_k1[1] <= tol_eps ) ) ref_k1[1] = 0.;

                    if ( ( ref_k1[0] - 1. >= -tol_eps ) && ( ref_k1[0] - 1. <= tol_eps ) ) ref_k1[0] = 1.;
                    if ( ( ref_k1[1] - 1. >= -tol_eps ) && ( ref_k1[1] - 1. <= tol_eps ) ) ref_k1[1] = 1.;

                    Coord ref_check ( 2 );
                    ref_check[0] = ref_k1[0];
                    ref_check[1] = ref_k1[1];
                    while ( !( this->contains_reference_point ( ref_check ) ) )
                    {
                        omega /= 2.0;
                        ref_k1 = ref_k + omega * update_vec;

                        if ( ( ref_k1[0] >= -tol_eps ) && ( ref_k1[0] <= tol_eps ) ) ref_k1[0] = 0.;
                        if ( ( ref_k1[1] >= -tol_eps ) && ( ref_k1[1] <= tol_eps ) ) ref_k1[1] = 0.;

                        if ( ( ref_k1[0] - 1. >= -tol_eps ) && ( ref_k1[0] - 1. <= tol_eps ) ) ref_k1[0] = 1.;
                        if ( ( ref_k1[1] - 1. >= -tol_eps ) && ( ref_k1[1] - 1. <= tol_eps ) ) ref_k1[1] = 1.;

                        ref_check[0] = ref_k1[0];
                        ref_check[1] = ref_k1[1];
                    }

                    Coord coord_ref_k1 ( 2 );
                    coord_ref_k1[0] = ref_k1[0];
                    coord_ref_k1[1] = ref_k1[1];

                    Vec<2, DataType> F_k1;
                    F_k1[0] = this->x ( coord_ref_k1 );
                    F_k1[1] = this->y ( coord_ref_k1 );

                    residual_armijo = norm ( F_k1 - pt_phy );

                    ++iter_armijo;
                    omega /= 2.;

                }
                abserr = norm ( ref_k1 - ref_k );

                ref_k = ref_k1;

                coord_ref_k[0] = ref_k[0];
                coord_ref_k[1] = ref_k[1];

                pt_k[0] = this->x ( coord_ref_k );
                pt_k[1] = this->y ( coord_ref_k );

                residual = norm ( pt_phy - pt_k );

                ++iter;

                if ( iter > iter_max ) break;

            } // end newton

            // Set values ...

            x_ref = ref_k[0];
            y_ref = ref_k[1];
        }

        /// \details The Inverse of the mapping (reference cell to physical cell)
        ///          is computed in 3D via a Newton-scheme with Armijo updates. The
        ///          Problem, which is to be solved writes<br>
        ///          G([x_ref, y_ref, z_ref]) := [x_phy - x(coord_ref),
        ///                                       y_phy - y(coord_ref),
        ///                                       z_phy - z(coord_ref)]<br>
        ///          solve G([x_ref, y_ref, z_ref]) = 0 via Newton-scheme
        /// \param[in] x_phy x coordinate on physical cell
        /// \param[in] y_phy y coordinate on physical cell
        /// \param[in] z_phy z coordinate on physical cell
        /// \param[out] x_ref x coordinate on reference cell
        /// \param[out] y_ref y coordinate on reference cell
        /// \param[out] z_ref z coordinate on reference cell

        template<class DataType>
        void CellTransformation<DataType>::inverse_newton_3d ( DataType x_phy, DataType y_phy, DataType z_phy,
                                                               DataType& x_ref, DataType& y_ref, DataType& z_ref ) const
        {
            // Initialisation

            Vec<3, DataType> pt_phy;
            pt_phy[0] = x_phy;
            pt_phy[1] = y_phy;
            pt_phy[2] = z_phy;

            Vec<3, DataType> ref_k1, ref_k;

            ref_k1[0] = 0.;
            ref_k1[1] = 0.;
            ref_k1[2] = 0.;

            ref_k[0] = 0.1154;
            ref_k[1] = 0.1832;
            ref_k[2] = 0.1385;

            const DataType tol_eps = 1.e3 * std::numeric_limits<DataType>::epsilon ( );

            // Some general parameters

            const int iter_max = 1000;
            int iter = 0;

            // Residual

            Coord coord_ref_k ( 3 );
            coord_ref_k[0] = ref_k[0];
            coord_ref_k[1] = ref_k[1];
            coord_ref_k[2] = ref_k[2];

            Vec<3, DataType> pt_k;
            pt_k[0] = this->x ( coord_ref_k );
            pt_k[1] = this->y ( coord_ref_k );
            pt_k[2] = this->z ( coord_ref_k );

            DataType residual = norm ( pt_phy - pt_k );
            DataType abserr = norm ( ref_k1 - ref_k );

            // Newton

            // Jacobian Matrix (grad G)
            Mat<3, 3, DataType> G;
            // Inverse of the jacobian
            Mat<3, 3, DataType> B;

            while ( ( residual > 5. * std::numeric_limits<DataType>::epsilon ( ) ) && ( abserr > 5. * std::numeric_limits<DataType>::epsilon ( ) ) )
            {

                G ( 0, 0 ) = this->x_x ( coord_ref_k );
                G ( 0, 1 ) = this->x_y ( coord_ref_k );
                G ( 0, 2 ) = this->x_z ( coord_ref_k );

                G ( 1, 0 ) = this->y_x ( coord_ref_k );
                G ( 1, 1 ) = this->y_y ( coord_ref_k );
                G ( 1, 2 ) = this->y_z ( coord_ref_k );

                G ( 2, 0 ) = this->z_x ( coord_ref_k );
                G ( 2, 1 ) = this->z_y ( coord_ref_k );
                G ( 2, 2 ) = this->z_z ( coord_ref_k );

#ifndef NDEBUG
                const DataType detG = det ( G );

                assert ( detG != 0. );
#endif

                inv ( G, B );

                // Armijo parameter

                const int iter_max_armijo = 500;
                int iter_armijo = 0;

                DataType residual_armijo = 2. * residual;

                DataType omega = 1.;

                const Vec<3, DataType> update_vec = B * ( pt_phy - pt_k );

                // Start Armijo

                while ( ( iter_armijo <= iter_max_armijo )
                        && ( residual_armijo > residual ) )
                {

                    ref_k1 = ref_k + omega * update_vec;

                    if ( ( ref_k1[0] >= -tol_eps ) && ( ref_k1[0] <= tol_eps ) ) ref_k1[0] = 0.;
                    if ( ( ref_k1[1] >= -tol_eps ) && ( ref_k1[1] <= tol_eps ) ) ref_k1[1] = 0.;
                    if ( ( ref_k1[2] >= -tol_eps ) && ( ref_k1[2] <= tol_eps ) ) ref_k1[2] = 0.;

                    if ( ( ref_k1[0] - 1. >= -tol_eps ) && ( ref_k1[0] - 1. <= tol_eps ) ) ref_k1[0] = 1.;
                    if ( ( ref_k1[1] - 1. >= -tol_eps ) && ( ref_k1[1] - 1. <= tol_eps ) ) ref_k1[1] = 1.;
                    if ( ( ref_k1[2] - 1. >= -tol_eps ) && ( ref_k1[2] - 1. <= tol_eps ) ) ref_k1[2] = 1.;

                    Coord ref_check ( 3 );
                    ref_check[0] = ref_k1[0];
                    ref_check[1] = ref_k1[1];
                    ref_check[2] = ref_k1[2];
                    while ( !( this->contains_reference_point ( ref_check ) ) )
                    {
                        omega /= 2.0;
                        ref_k1 = ref_k + omega * update_vec;

                        if ( ( ref_k1[0] >= -tol_eps ) && ( ref_k1[0] <= tol_eps ) ) ref_k1[0] = 0.;
                        if ( ( ref_k1[1] >= -tol_eps ) && ( ref_k1[1] <= tol_eps ) ) ref_k1[1] = 0.;
                        if ( ( ref_k1[2] >= -tol_eps ) && ( ref_k1[2] <= tol_eps ) ) ref_k1[2] = 0.;

                        if ( ( ref_k1[0] - 1. >= -tol_eps ) && ( ref_k1[0] - 1. <= tol_eps ) ) ref_k1[0] = 1.;
                        if ( ( ref_k1[1] - 1. >= -tol_eps ) && ( ref_k1[1] - 1. <= tol_eps ) ) ref_k1[1] = 1.;
                        if ( ( ref_k1[2] - 1. >= -tol_eps ) && ( ref_k1[2] - 1. <= tol_eps ) ) ref_k1[2] = 1.;

                        ref_check[0] = ref_k1[0];
                        ref_check[1] = ref_k1[1];
                        ref_check[2] = ref_k1[2];
                    }

                    Coord coord_ref_k1 ( 3 );
                    coord_ref_k1[0] = ref_k1[0];
                    coord_ref_k1[1] = ref_k1[1];
                    coord_ref_k1[2] = ref_k1[2];

                    Vec<3, DataType> F_k1;
                    F_k1[0] = this->x ( coord_ref_k1 );
                    F_k1[1] = this->y ( coord_ref_k1 );
                    F_k1[2] = this->z ( coord_ref_k1 );

                    residual_armijo = norm ( pt_phy - F_k1 );

                    ++iter_armijo;
                    omega /= 2.;
                }

                abserr = norm ( ref_k1 - ref_k );
                ref_k = ref_k1;

                coord_ref_k[0] = ref_k[0];
                coord_ref_k[1] = ref_k[1];
                coord_ref_k[2] = ref_k[2];

                pt_k[0] = x ( coord_ref_k );
                pt_k[1] = y ( coord_ref_k );
                pt_k[2] = z ( coord_ref_k );

                residual = norm ( pt_phy - pt_k );

                ++iter;
                if ( iter > iter_max ) break;

            } // end newton

            // Set values ...

            x_ref = ref_k[0];
            y_ref = ref_k[1];
            z_ref = ref_k[2];
        }

        template<class DataType>
        bool CellTransformation<DataType>::contains_physical_point ( const Coord& coord_phys,
                                                                     Coord* coord_ref ) const
        {
            if ( coord_phys.size ( ) != this->gdim_ )
            {
                throw "Incorrect geometrical dimension of coord_phys in CellTransformation::contains_physical_point!";
            }

            std::vector<DataType> cr ( this->gdim_ );

            // Compute reference coordinates
            switch ( this->gdim_ )
            {
                case 1:
                {
                    this->inverse ( coord_phys[0], cr[0] );
                    break;
                }
                case 2:
                {
                    this->inverse ( coord_phys[0], coord_phys[1], cr[0], cr[1] );
                    break;
                }
                case 3:
                {
                    this->inverse ( coord_phys[0], coord_phys[1], coord_phys[2], cr[0], cr[1], cr[2] );
                    break;
                }
                default:
                {
                    throw "Invalid dimension!";
                    break;
                }
            }

            const bool contains_pt = this->contains_reference_point ( cr );

            // Output parameter
            if ( contains_pt && coord_ref )
            {
                cr.swap ( *coord_ref );
            }

            return contains_pt;
        }

        // template instantiation
        template class CellTransformation<double>;
        template class CellTransformation<float>;

    } // namespace doffem
} // namespace hiflow
