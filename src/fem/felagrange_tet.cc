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

#include "felagrange_tet.h"
#include <cassert>

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        FELagrangeTet<DataType>::FELagrangeTet ( )
        {
            this->my_id_ = FEType<DataType>::LAGRANGE_TET;

            // initialize reference cell

            assert ( this->ref_cell_ == NULL );
            this->ref_cell_ = &( mesh::CellType::get_instance ( mesh::CellType::TETRAHEDRON ) );
        }

        template<class DataType>
        FELagrangeTet<DataType>::~FELagrangeTet ( )
        {
        }

        template<class DataType>
        void FELagrangeTet<DataType>::init_coord ( )
        {
            // set topological degree

            this->tdim_ = 3;

            // Lexicographical ordering

            if ( this->fe_deg_ == 0 )
            {
                this->coord_.clear ( );

                // Coordinates of the middle-point of tet
                Coord coord;
                coord.resize ( 3 );

                // Centroid
                coord[0] = 0.25;
                coord[1] = 0.25;
                coord[2] = 0.25;

                this->coord_.push_back ( coord );
            }
            else
            {
                assert ( this->fe_deg_ > 0 );

                const DataType offset = ( 1.0 / this->fe_deg_ );

                this->coord_.clear ( );
                const int nb_dof_on_cell = ( this->fe_deg_ + 1 )*( this->fe_deg_ + 2 )*( this->fe_deg_ + 3 ) / 6;
                this->coord_.resize ( nb_dof_on_cell );

                // Filling coord vector for full tet by lexicographical strategy
                // and filtering the line coordinates by given mesh ordering strategy
                // which was computed in first step

                const int nb_dof_line = this->fe_deg_ + 1;

                for ( int k = 0; k < nb_dof_line; ++k )
                { // z axis
                    const DataType k_offset = k*offset;
                    for ( int j = 0; j < nb_dof_line - k; ++j )
                    { // y axis
                        const DataType j_offset = j*offset;
                        for ( int i = 0; i < nb_dof_line - k - j; ++i ) // x axis
                        {
                            Coord coord;
                            coord.resize ( 3 );

                            coord[0] = i*offset;
                            coord[1] = j_offset;
                            coord[2] = k_offset;

                            this->coord_[ijk2ind ( i, j, k )] = coord;
                        }
                    }
                }
            }
        }

        /// \details The counting of the dofs on the reference cell is done by the
        ///          lexicographical numbering strategy. Therefore, beginning on the
        ///          x coordinate, then continuing on the y coordinate and last continuing
        ///          on the z coordinate this is achieved by computing the corresponding
        ///          offsets to consider the restriction given by the tetrahedron which
        ///          reads z in [0,1], y < 1 - z, x < 1 - y - z

        template<class DataType>
        int FELagrangeTet<DataType>::ijk2ind ( int i, int j, int k ) const
        {
            // x component = i, y component = j, z component = k

            int offset = 0;
            const int nb_dof_line = this->fe_deg_ + 1;

            // First: offset z axis

            for ( int m = 0; m < k; ++m )
            {
                const int help = nb_dof_line - m;
                for ( int dof = 0; dof < nb_dof_line - m; ++dof )
                {
                    offset += help - dof;
                }
            }

            // Second: increasing offset by y axis on current z axis

            for ( int n = 0; n < j; ++n )
                offset += nb_dof_line - n - k;

            return (i + offset );
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree fe_deg_". Since fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$

        template<class DataType>
        void FELagrangeTet<DataType>::N ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = this->lp_.poly ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 1.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = this->lp_.poly ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( this->fe_deg_, 0, 0 )] = this->lp_.poly ( this->fe_deg_, this->fe_deg_, pt[0] );

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, this->fe_deg_, 0 )] = this->lp_.poly ( this->fe_deg_, this->fe_deg_, pt[1] );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                * this->lp_.poly ( j, j, deg * pt[1] / j );

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j );
                }

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, this->fe_deg_ )] = this->lp_.poly ( this->fe_deg_, this->fe_deg_, pt[2] );

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, deg * pt[2] / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, dp_0 / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = this->lp_.poly ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, dp_1 / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = this->lp_.poly ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = this->lp_.poly ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree fe_deg_". Since fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for x are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_x ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = -this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = -( deg / ( deg - i ) ) * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i )
                + this->lp_.poly ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i );

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( this->fe_deg_, 0, 0 )] = this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, pt[0] );

            weight[ijk2ind ( 0, this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = -( deg / ( deg - j ) ) * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = ( deg / ( deg - j ) ) * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = -( deg / ( deg - i - j ) ) * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j );
                }

            weight[ijk2ind ( 0, 0, this->fe_deg_ )] = 0.0;

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = -( deg / ( deg - k ) ) * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = ( deg / ( deg - k ) ) * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, dp_0 / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = -( deg / ( deg - i - k ) ) * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    + this->lp_.poly ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = 0.0;
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = -( deg / ( deg - j - k ) ) * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = ( deg / ( deg - k - j ) ) * this->lp_.poly_x ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = -( deg / ( deg - i - j - k ) ) * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for y are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_y ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = -this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = -( deg / ( deg - i ) ) * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            weight[ijk2ind ( this->fe_deg_, 0, 0 )] = 0.0;

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, this->fe_deg_, 0 )] = this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, pt[1] );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = -( deg / ( deg - j ) ) * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j )
                + this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = -( deg / ( deg - i - j ) ) * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );
                }

            weight[ijk2ind ( 0, 0, this->fe_deg_ )] = 0.0;

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = -( deg / ( deg - k ) ) * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = 0.0;

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = -( deg / ( deg - i - k ) ) * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = ( deg / ( deg - k ) ) * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, deg * pt[1] / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = -( deg / ( deg - j - k ) ) * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            + this->lp_.poly ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = this->lp_.poly ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = -( deg / ( deg - i - j - k ) ) * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for z are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_z ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = -this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = -( deg / ( deg - i ) ) * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            weight[ijk2ind ( this->fe_deg_, 0, 0 )] = 0.0;

            weight[ijk2ind ( 0, this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = -( deg / ( deg - j ) ) * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = 0.0;

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = -( deg / ( deg - i - j ) ) * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j );
                }

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, this->fe_deg_ )] = this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, pt[2] );

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = -( deg / ( deg - k ) ) * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k )
                        + this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, dp_0 / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = -( deg / ( deg - i - k ) ) * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    + this->lp_.poly ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, dp_1 / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = -( deg / ( deg - j - k ) ) * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            + this->lp_.poly ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = this->lp_.poly ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = -( deg / ( deg - i - j - k ) ) * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for xx are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_xx ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = ( deg / ( deg - i ) )
                *( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i )
                -( deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                *( deg / i )
                * this->lp_.poly_x ( i, i, dp_0 / i )
                -( deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                + this->lp_.poly ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                *( deg / i )*( deg / i )
                * this->lp_.poly_xx ( i, i, dp_0 / i );

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( this->fe_deg_, 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, pt[0] );

            weight[ijk2ind ( 0, this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = ( deg / ( deg - i - j ) )*( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            *( deg / i )*( deg / i )
                            * this->lp_.poly_xx ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j );
                }

            weight[ijk2ind ( 0, 0, this->fe_deg_ )] = 0.0;

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = ( deg / ( deg - k ) )*( deg / ( deg - k ) )
                        * this->lp_.poly_xx ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = ( deg / ( deg - k ) )*( deg / ( deg - k ) )
                        * this->lp_.poly_xx ( this->fe_deg_ - k, this->fe_deg_ - k, dp_0 / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = ( deg / ( deg - i - k ) )*( deg / ( deg - i - k ) )
                    * this->lp_.poly_xx ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    -( deg / ( deg - i - k ) )
                    * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    *( deg / i )
                    * this->lp_.poly_x ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    -( deg / ( deg - i - k ) )
                    * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    + this->lp_.poly ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    *( deg / i )*( deg / i )
                    * this->lp_.poly_xx ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = 0.0;
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = ( deg / ( deg - j - k ) )*( deg / ( deg - j - k ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = ( deg / ( deg - k - j ) )*( deg / ( deg - k - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = ( deg / ( deg - i - j - k ) )*( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_xx ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                *( deg / i )*( deg / i )
                                * this->lp_.poly_xx ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for xy are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_xy ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = ( deg / ( deg - i ) )*( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i )
                -( deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i );

            weight[ijk2ind ( this->fe_deg_, 0, 0 )] = 0.0;

            weight[ijk2ind ( 0, this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j )
                -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = ( deg / ( deg - j ) ) * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = ( deg / ( deg - i - j ) )*( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );
                }

            weight[ijk2ind ( 0, 0, this->fe_deg_ )] = 0.0;

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = ( deg / ( deg - k ) )*( deg / ( deg - k ) )
                        * this->lp_.poly_xx ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = 0.0;

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = ( deg / ( deg - i - k ) )*( deg / ( deg - i - k ) )
                    * this->lp_.poly_xx ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    -( deg / ( deg - i - k ) )
                    * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = 0.0;
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = ( deg / ( deg - j - k ) )*( deg / ( deg - j - k ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            -( deg / ( deg - j - k ) )
                            * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = ( deg / ( deg - k - j ) ) * this->lp_.poly_x ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = ( deg / ( deg - i - j - k ) )*( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_xx ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for xz are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_xz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = ( deg / ( deg - i ) )*( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i )
                -( deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i );

            weight[ijk2ind ( this->fe_deg_, 0, 0 )] = 0.0;

            weight[ijk2ind ( 0, this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = 0.0;

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = ( deg / ( deg - i - j ) )*( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j );
                }

            weight[ijk2ind ( 0, 0, this->fe_deg_ )] = 0.0;

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = ( deg / ( deg - k ) )*( deg / ( deg - k ) )
                        * this->lp_.poly_xx ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k )
                        -( deg / ( deg - k ) )
                        * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = ( deg / ( deg - k ) ) * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, dp_0 / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = ( deg / ( deg - i - k ) )*( deg / ( deg - i - k ) )
                    * this->lp_.poly_xx ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    -( deg / ( deg - i - k ) )
                    * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )

                    -( deg / ( deg - i - k ) )
                    * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    + this->lp_.poly ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                    *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = 0.0;
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = ( deg / ( deg - j - k ) )*( deg / ( deg - j - k ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            -( deg / ( deg - j - k ) )
                            * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = ( deg / ( deg - k - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = ( deg / ( deg - i - j - k ) )*( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_xx ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )

                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for yy are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_yy ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = ( deg / ( deg - i ) )*( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            weight[ijk2ind ( this->fe_deg_, 0, 0 )] = 0.0;

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, this->fe_deg_, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, pt[1] );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j )
                -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                + this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                *( deg / j )*( deg / j )
                * this->lp_.poly_xx ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                *( deg / j )*( deg / j )
                * this->lp_.poly_xx ( j, j, dp_1 / j );

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = ( deg / ( deg - i - j ) )*( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j )*( deg / j )
                            * this->lp_.poly_xx ( j, j, dp_1 / j );
                }

            weight[ijk2ind ( 0, 0, this->fe_deg_ )] = 0.0;

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = ( deg / ( deg - k ) )*( deg / ( deg - k ) )
                        * this->lp_.poly_xx ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = 0.0;

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = ( deg / ( deg - i - k ) )*( deg / ( deg - i - k ) )
                    * this->lp_.poly_xx ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = ( deg / ( deg - k ) )*( deg / ( deg - k ) )
                        * this->lp_.poly_xx ( this->fe_deg_ - k, this->fe_deg_ - k, dp_1 / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k );
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = ( deg / ( deg - j - k ) )*( deg / ( deg - j - k ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            -( deg / ( deg - j - k ) )
                            * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            -( deg / ( deg - j - k ) )
                            * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            + this->lp_.poly ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            *( deg / j )*( deg / j ) * this->lp_.poly_xx ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = this->lp_.poly ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            *( deg / j )*( deg / j ) * this->lp_.poly_xx ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = ( deg / ( deg - i - j - k ) )*( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_xx ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                *( deg / j )*( deg / j ) * this->lp_.poly_xx ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for yz are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_yz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = ( deg / ( deg - i ) )*( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            weight[ijk2ind ( this->fe_deg_, 0, 0 )] = 0.0;

            weight[ijk2ind ( 0, this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j )
                -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = 0.0;

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = ( deg / ( deg - i - j ) )*( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );
                }

            weight[ijk2ind ( 0, 0, this->fe_deg_ )] = 0.0;

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = ( deg / ( deg - k ) )*( deg / ( deg - k ) )
                        * this->lp_.poly_xx ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k )
                        -( deg / ( deg - k ) )
                        * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = 0.0;

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = ( deg / ( deg - i - k ) )*( deg / ( deg - i - k ) )
                    * this->lp_.poly_xx ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    -( deg / ( deg - i - k ) )
                    * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = ( deg / ( deg - k ) ) * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, dp_1 / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = ( deg / ( deg - j - k ) )*( deg / ( deg - j - k ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            -( deg / ( deg - j - k ) )
                            * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                            -( deg / ( deg - j - k ) )
                            * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            + this->lp_.poly ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = this->lp_.poly ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = ( deg / ( deg - i - j - k ) )*( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_xx ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                                *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k );
                    }
        }

        /// \details The restriction of lagrangian finite elements on a tetrahedron reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y,z) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y-z, x, y, z). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j-k}_{d-i-j-k} ((d/(d-i-j-k)^*(1-x-y-z))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)^*L^k_k(d/k^*z)\f$
        ///          Here, the derivatives for zz are considered via the chain rule.

        template<class DataType>
        void FELagrangeTet<DataType>::N_zz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help = 1.0 - pt[0] - pt[1] - pt[2];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];
            const DataType dp_2 = deg * pt[2];

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help );
            else
                weight[ijk2ind ( 0, 0, 0 )] = 0.0;

            for ( int i = 1; i< this->fe_deg_; ++i )
                weight[ijk2ind ( i, 0, 0 )] = ( deg / ( deg - i ) )*( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            weight[ijk2ind ( this->fe_deg_, 0, 0 )] = 0.0;

            weight[ijk2ind ( 0, this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( 0, j, 0 )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ijk2ind ( this->fe_deg_ - j, j, 0 )] = 0.0;

            for ( int j = 1; j< this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ijk2ind ( i, j, 0 )] = ( deg / ( deg - i - j ) )*( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j );
                }

            if ( this->fe_deg_ > 0 )
                weight[ijk2ind ( 0, 0, this->fe_deg_ )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, pt[2] );

            for ( int k = 1; k<this->fe_deg_; ++k )
            {
                weight[ijk2ind ( 0, 0, k )] = ( deg / ( deg - k ) )*( deg / ( deg - k ) )
                        * this->lp_.poly_xx ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        * this->lp_.poly ( k, k, dp_2 / k )
                        -( deg / ( deg - k ) )
                        * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                        -( deg / ( deg - k ) )
                        * this->lp_.poly_x ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                        + this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, deg * help / ( deg - k ) )
                        *( deg / k )*( deg / k ) * this->lp_.poly_xx ( k, k, dp_2 / k );

                weight[ijk2ind ( this->fe_deg_ - k, 0, k )] = this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, dp_0 / ( deg - k ) )
                        *( deg / k )*( deg / k ) * this->lp_.poly_xx ( k, k, dp_2 / k );

                for ( int i = 1; i< this->fe_deg_ - k; ++i )
                    weight[ijk2ind ( i, 0, k )] = ( deg / ( deg - i - k ) )*( deg / ( deg - i - k ) )
                    * this->lp_.poly_xx ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    * this->lp_.poly ( k, k, dp_2 / k )
                    -( deg / ( deg - i - k ) )
                    * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                    -( deg / ( deg - i - k ) )
                    * this->lp_.poly_x ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                    + this->lp_.poly ( this->fe_deg_ - i - k, this->fe_deg_ - i - k, deg * help / ( deg - i - k ) )
                    * this->lp_.poly ( i, i, dp_0 / i )
                    *( deg / k )*( deg / k ) * this->lp_.poly_xx ( k, k, dp_2 / k );

                weight[ijk2ind ( 0, this->fe_deg_ - k, k )] = this->lp_.poly ( this->fe_deg_ - k, this->fe_deg_ - k, dp_1 / ( deg - k ) )
                        *( deg / k )*( deg / k ) * this->lp_.poly_xx ( k, k, dp_2 / k );
            }

            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                {
                    weight[ijk2ind ( 0, j, k )] = ( deg / ( deg - j - k ) )*( deg / ( deg - j - k ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            * this->lp_.poly ( k, k, dp_2 / k )
                            -( deg / ( deg - j - k ) )
                            * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                            -( deg / ( deg - j - k ) )
                            * this->lp_.poly_x ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                            + this->lp_.poly ( this->fe_deg_ - j - k, this->fe_deg_ - j - k, deg * help / ( deg - j - k ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k )*( deg / k ) * this->lp_.poly_xx ( k, k, dp_2 / k );

                    weight[ijk2ind ( this->fe_deg_ - k - j, j, k )] = this->lp_.poly ( this->fe_deg_ - k - j, this->fe_deg_ - k - j, dp_0 / ( deg - k - j ) )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            *( deg / k )*( deg / k ) * this->lp_.poly_xx ( k, k, dp_2 / k );
                }

            // Main "for" loop
            for ( int k = 1; k<this->fe_deg_; ++k )
                for ( int j = 1; j<this->fe_deg_ - k; ++j )
                    for ( int i = 1; i<this->fe_deg_ - k - j; ++i )
                    {
                        weight[ijk2ind ( i, j, k )] = ( deg / ( deg - i - j - k ) )*( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_xx ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                * this->lp_.poly ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                                -( deg / ( deg - i - j - k ) )
                                * this->lp_.poly_x ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                *( deg / k ) * this->lp_.poly_x ( k, k, dp_2 / k )
                                + this->lp_.poly ( this->fe_deg_ - i - j - k, this->fe_deg_ - i - j - k, deg * help / ( deg - i - j - k ) )
                                * this->lp_.poly ( i, i, dp_0 / i )
                                * this->lp_.poly ( j, j, dp_1 / j )
                                *( deg / k )*( deg / k ) * this->lp_.poly_xx ( k, k, dp_2 / k );
                    }
        }

        template class FELagrangeTet<double>;
        template class FELagrangeTet<float>;

    } // namespace doffem
} // namespace hiflow
