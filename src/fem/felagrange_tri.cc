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

#include "felagrange_tri.h"
#include <cassert>
#include <cmath>

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        FELagrangeTri<DataType>::FELagrangeTri ( )
        {
            this->my_id_ = FEType<DataType>::LAGRANGE_TRI;

            // initialize reference cell

            assert ( this->ref_cell_ == NULL );
            this->ref_cell_ = &( mesh::CellType::get_instance ( mesh::CellType::TRIANGLE ) );
        }

        template<class DataType>
        FELagrangeTri<DataType>::~FELagrangeTri ( )
        {
        }

        template<class DataType>
        void FELagrangeTri<DataType>::init_coord ( )
        {
            // set topological degree

            this->tdim_ = 2;

            // Lexicographical ordering

            if ( this->fe_deg_ == 0 )
            {
                this->coord_.clear ( );

                // Coordinates of the middle-point of triangle
                Coord coord;
                coord.resize ( 2 );

                // Centroid
                coord[0] = 1.0 / 3.0;
                coord[1] = 1.0 / 3.0;

                this->coord_.push_back ( coord );
            }
            else
            {
                assert ( this->fe_deg_ > 0 );

                DataType offset = ( 1.0 / this->fe_deg_ );

                this->coord_.clear ( );
                int nb_dof_on_cell = ( ( ( this->fe_deg_ + 1 )*( this->fe_deg_ + 2 ) ) / 2 );
                this->coord_.resize ( nb_dof_on_cell );

                // Filling coord vector for full triangle by lexicographical strategy
                // and filtering the line coordinates by given mesh ordering strategy
                // which was computed in first step

                int nb_dof_line = this->fe_deg_ + 1;

                for ( int j = 0; j < nb_dof_line; ++j )
                { // y axis
                    const DataType j_offset = j*offset;
                    for ( int i = 0; i < nb_dof_line - j; ++i ) // x axis
                    {
                        Coord coord;
                        coord.resize ( 2 );

                        coord[0] = i*offset;
                        coord[1] = j_offset;

                        this->coord_[ij2ind ( i, j )] = coord;
                    }
                }
            }
        }

        /// \details The counting of the dofs on the reference cell is done by the
        ///          lexicographical numbering strategy. Therefore, beginning on the
        ///          x coordinate, then continuing on the y coordinate this is achieved
        ///          by computing the corresponding offsets to consider the restriction
        ///          given by the triangle which reads y in [0,1], x < 1 - y

        template<class DataType>
        int FELagrangeTri<DataType>::ij2ind ( int i, int j ) const
        {
            int offset = 0;
            const int nb_dof_line = this->fe_deg_ + 1;

            for ( int n = 0; n < j; ++n )
                offset += nb_dof_line - n;

            return (i + offset );
        }

        /// \details The restriction of lagrangian finite elements on a triangle reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y in cartesian sense, the
        ///          barycentric coordinates read (1-x-y, x, y). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j}_{d-i-j} ((d/(d-i-j)^*(1-x-y))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)\f$

        template<class DataType>
        void FELagrangeTri<DataType>::N ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );

            const DataType help1 = 1.0 - pt[0] - pt[1];
            const DataType dh_1 = deg * help1;

            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, 0 )] = this->lp_.poly ( this->fe_deg_, this->fe_deg_, help1 );
            else
                weight[ij2ind ( 0, 0 )] = 1.0;

            for ( int i = 1; i<this->fe_deg_; ++i )
                weight[ij2ind ( i, 0 )] = this->lp_.poly ( this->fe_deg_ - i, this->fe_deg_ - i, dh_1 / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( this->fe_deg_, 0 )] = this->lp_.poly ( this->fe_deg_, this->fe_deg_, pt[0] );
            else
                weight[ij2ind ( this->fe_deg_, 0 )] = 1.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( 0, j )] = this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, dh_1 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, this->fe_deg_ )] = this->lp_.poly ( this->fe_deg_, this->fe_deg_, pt[1] );
            else
                weight[ij2ind ( 0, this->fe_deg_ )] = 1.0;

            // Main "for" loop
            for ( int j = 1; j<this->fe_deg_; ++j )
            {
                const DataType lp_j = this->lp_.poly ( j, j, dp_1 / j );
                const int offset = this->fe_deg_ - j;
                const DataType offset_double = static_cast < DataType > ( deg - j );
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ij2ind ( i, j )] = this->lp_.poly ( offset - i, offset - i, dh_1 / ( offset_double - i ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * lp_j;
                }
            }

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( this->fe_deg_ - j, j )] = this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );
        }

        /// \details The restriction of lagrangian finite elements on a triangle reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y, x, y). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j}_{d-i-j} ((d/(d-i-j)^*(1-x-y))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)\f$
        ///          Here, the derivatives for x are considered via the chain rule.

        template<class DataType>
        void FELagrangeTri<DataType>::N_x ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );

            const DataType help1 = 1.0 - pt[0] - pt[1];
            const DataType dh_1 = deg * help1;

            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, 0 )] = -this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, help1 );
            else
                weight[ij2ind ( 0, 0 )] = 0.0;

            for ( int i = 1; i<this->fe_deg_; ++i )
                weight[ij2ind ( i, 0 )] = ( -deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, dh_1 / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i )
                + this->lp_.poly ( this->fe_deg_ - i, this->fe_deg_ - i, dh_1 / ( deg - i ) )
                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i );

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( this->fe_deg_, 0 )] = this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, pt[0] );
            else
                weight[ij2ind ( this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( 0, j )] = -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, dh_1 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            weight[ij2ind ( 0, this->fe_deg_ )] = 0.0;

            // Main "for" loop
            for ( int j = 1; j<this->fe_deg_; ++j )
            {
                const DataType lp_j = this->lp_.poly ( j, j, dp_1 / j );
                const int offset = this->fe_deg_ - j;
                const DataType offset_double = static_cast < DataType > ( deg - j );
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ij2ind ( i, j )] = -( deg / ( offset_double - i ) )
                            * this->lp_.poly_x ( offset - i, offset - i, dh_1 / ( offset_double - i ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * lp_j
                            + this->lp_.poly ( offset - i, offset - i, dh_1 / ( offset_double - i ) )
                            *( deg / i )
                            * this->lp_.poly_x ( i, i, dp_0 / i )
                            * lp_j;
                }
            }

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( this->fe_deg_ - j, j )] = ( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );
        }

        /// \details The restriction of lagrangian finite elements on a triangle reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y, x, y). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j}_{d-i-j} ((d/(d-i-j)^*(1-x-y))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)\f$
        ///          Here, the derivatives for y are considered via the chain rule.

        template<class DataType>
        void FELagrangeTri<DataType>::N_y ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help1 = 1.0 - pt[0] - pt[1];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, 0 )] = -this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, help1 );
            else
                weight[ij2ind ( 0, 0 )] = 0.0;

            for ( int i = 1; i<this->fe_deg_; ++i )
                weight[ij2ind ( i, 0 )] = -( deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help1 / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            weight[ij2ind ( this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( 0, j )] = -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j )
                + this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, this->fe_deg_ )] = this->lp_.poly_x ( this->fe_deg_, this->fe_deg_, pt[1] );
            else
                weight[ij2ind ( 0, this->fe_deg_ )] = 0.0;

            // Main "for" loop
            for ( int j = 1; j<this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ij2ind ( i, j )] = -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, deg * pt[1] / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j ) * this->lp_.poly_x ( j, j, deg * pt[1] / j );
                }

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( this->fe_deg_ - j, j )] = this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                *( deg / j )
                * this->lp_.poly_x ( j, j, dp_1 / j );
        }

        /// \details There are no z derivatives in 2D

        template<class DataType>
        void FELagrangeTri<DataType>::N_z ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details The restriction of lagrangian finite elements on a triangle reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y, x, y). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j}_{d-i-j} ((d/(d-i-j)^*(1-x-y))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)\f$
        ///          Here, the derivatives for xx are considered via the chain rule.

        template<class DataType>
        void FELagrangeTri<DataType>::N_xx ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help1 = 1.0 - pt[0] - pt[1];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help1 );
            else
                weight[ij2ind ( 0, 0 )] = 0.0;

            for ( int i = 1; i<this->fe_deg_; ++i )
                weight[ij2ind ( i, 0 )] = ( deg / ( deg - i ) )
                *( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help1 / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i )
                -( deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help1 / ( deg - i ) )
                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                -( deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help1 / ( deg - i ) )
                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                + this->lp_.poly ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help1 / ( deg - i ) )
                *( deg / i )*( deg / i )
                * this->lp_.poly_xx ( i, i, dp_0 / i );

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( this->fe_deg_, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, pt[0] );
            else
                weight[ij2ind ( this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( 0, j )] = ( deg / ( deg - j ) )
                *( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );

            weight[ij2ind ( 0, this->fe_deg_ )] = 0.0;

            // Main "for" loop
            for ( int j = 1; j<this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ij2ind ( i, j )] = ( deg / ( deg - i - j ) )
                            *( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            +( -deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            +( -deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            *( deg / i )
                            * this->lp_.poly_x ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            *( deg / i )*( deg / i )
                            * this->lp_.poly_xx ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j );
                }

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( this->fe_deg_ - j, j )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j );
        }

        /// \details The restriction of lagrangian finite elements on a triangle reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y, x, y). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j}_{d-i-j} ((d/(d-i-j)^*(1-x-y))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)\f$
        ///          Here, the derivatives for xy are considered via the chain rule.

        template<class DataType>
        void FELagrangeTri<DataType>::N_xy ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help1 = 1.0 - pt[0] - pt[1];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help1 );
            else
                weight[ij2ind ( 0, 0 )] = 0.0;

            for ( int i = 1; i<this->fe_deg_; ++i )
                weight[ij2ind ( i, 0 )] = ( deg / ( deg - i ) )
                *( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help1 / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i )
                -( deg / ( deg - i ) )
                * this->lp_.poly_x ( this->fe_deg_ - i, this->fe_deg_ - i, deg * help1 / ( deg - i ) )
                *( deg / i ) * this->lp_.poly_x ( i, i, dp_0 / i );

            weight[ij2ind ( this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( 0, j )] = ( deg / ( deg - j ) )
                *( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j )
                -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j );

            weight[ij2ind ( 0, this->fe_deg_ )] = 0.0;

            // Main "for" loop
            for ( int j = 1; j<this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ij2ind ( i, j )] = ( deg / ( deg - i - j ) )
                            *( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j )
                            * this->lp_.poly_x ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            *( deg / i )
                            * this->lp_.poly_x ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            *( deg / i )
                            * this->lp_.poly_x ( i, i, dp_0 / i )
                            *( deg / j )
                            * this->lp_.poly_x ( j, j, dp_1 / j );
                }

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( this->fe_deg_ - j, j )] = ( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                *( deg / j )
                * this->lp_.poly_x ( j, j, dp_1 / j );
        }

        /// \details There are no z derivatives in 2D

        template<class DataType>
        void FELagrangeTri<DataType>::N_xz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details The restriction of lagrangian finite elements on a triangle reads
        ///          "sum of all multiplied polynomial degrees is less or equal to the total
        ///          degree this->fe_deg_". Since this->fe_deg_ = 0 is also allowed, there are several
        ///          distinctions to be done. For performance reasons, the code becomes a
        ///          little bit trenched. But the main "for(int ...)", representing what is
        ///          really happening is found at the end of the function. The values
        ///          for the coordinates are transformed from the cartesian system to the
        ///          barycentric system. This means, given (x,y) in cartesian sense, the
        ///          barycentric coordinates read (1-x-y, x, y). Also, they need to be
        ///          scaled by the factor (this->fe_deg_ / polynomial degree). The
        ///          resulting combination of the polynomials which is computed is given by
        ///          \f$L^{d-i-j}_{d-i-j} ((d/(d-i-j)^*(1-x-y))^*L^i_i(d/i^*x)^*L^j_j(d/j^*y)\f$
        ///          Here, the derivatives for yy are considered via the chain rule.

        template<class DataType>
        void FELagrangeTri<DataType>::N_yy ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType deg = static_cast < DataType > ( this->fe_deg_ );
            const DataType help1 = 1.0 - pt[0] - pt[1];
            const DataType dp_0 = deg * pt[0];
            const DataType dp_1 = deg * pt[1];

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, 0 )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, help1 );
            else
                weight[ij2ind ( 0, 0 )] = 0.0;

            for ( int i = 1; i<this->fe_deg_; ++i )
                weight[ij2ind ( i, 0 )] = ( deg / ( deg - i ) )*( deg / ( deg - i ) )
                * this->lp_.poly_xx ( this->fe_deg_ - i, this->fe_deg_ - i, this->fe_deg_ * help1 / ( deg - i ) )
                * this->lp_.poly ( i, i, dp_0 / i );

            weight[ij2ind ( this->fe_deg_, 0 )] = 0.0;

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( 0, j )] = ( deg / ( deg - j ) )*( deg / ( deg - j ) )
                * this->lp_.poly_xx ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                * this->lp_.poly ( j, j, dp_1 / j )
                -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                *( deg / j )
                * this->lp_.poly_x ( j, j, dp_1 / j )
                -( deg / ( deg - j ) )
                * this->lp_.poly_x ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                + this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, deg * help1 / ( deg - j ) )
                *( deg / j )*( deg / j )
                * this->lp_.poly_xx ( j, j, dp_1 / j );

            if ( this->fe_deg_ > 0 )
                weight[ij2ind ( 0, this->fe_deg_ )] = this->lp_.poly_xx ( this->fe_deg_, this->fe_deg_, pt[1] );
            else
                weight[ij2ind ( 0, this->fe_deg_ )] = 0.0;

            // Main "for" loop
            for ( int j = 1; j<this->fe_deg_; ++j )
                for ( int i = 1; i<this->fe_deg_ - j; ++i )
                {
                    weight[ij2ind ( i, j )] = ( deg / ( deg - i - j ) )*( deg / ( deg - i - j ) )
                            * this->lp_.poly_xx ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            * this->lp_.poly ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j )
                            * this->lp_.poly_x ( j, j, dp_1 / j )
                            -( deg / ( deg - i - j ) )
                            * this->lp_.poly_x ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j ) * this->lp_.poly_x ( j, j, dp_1 / j )
                            + this->lp_.poly ( this->fe_deg_ - i - j, this->fe_deg_ - i - j, deg * help1 / ( deg - i - j ) )
                            * this->lp_.poly ( i, i, dp_0 / i )
                            *( deg / j )*( deg / j )
                            * this->lp_.poly_xx ( j, j, dp_1 / j );
                }

            for ( int j = 1; j<this->fe_deg_; ++j )
                weight[ij2ind ( this->fe_deg_ - j, j )] = this->lp_.poly ( this->fe_deg_ - j, this->fe_deg_ - j, dp_0 / ( deg - j ) )
                *( deg / j )*( deg / j )
                * this->lp_.poly_xx ( j, j, dp_1 / j );
        }

        /// \details There are no z derivatives in 2D

        template<class DataType>
        void FELagrangeTri<DataType>::N_yz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details There are no z derivatives in 2D

        template<class DataType>
        void FELagrangeTri<DataType>::N_zz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        template class FELagrangeTri<double>;
        template class FELagrangeTri<float>;

    } // namespace doffem
} // namespace hiflow
