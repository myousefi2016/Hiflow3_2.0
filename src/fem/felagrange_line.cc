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

#include "felagrange_line.h"
#include <cassert>

/// \author Michael Schick<br>Martin Baumann<br>Julian Kraemer<br>Simon Gawlok

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        FELagrangeLine<DataType>::FELagrangeLine ( )
        {
            this->my_id_ = FEType<DataType>::LAGRANGE_LINE;

            // initialize reference cell

            assert ( this->ref_cell_ == NULL );
            this->ref_cell_ = &( mesh::CellType::get_instance ( mesh::CellType::LINE ) );
        }

        template<class DataType>
        FELagrangeLine<DataType>::~FELagrangeLine ( )
        {
        }

        template<class DataType>
        void FELagrangeLine<DataType>::init_coord ( )
        {
            // set topological degree

            this->tdim_ = 1;

            if ( this->fe_deg_ == 0 )
            {
                this->coord_.clear ( );

                // Coordinates of the middle-point of line
                Coord coord;
                coord.resize ( 1 );

                // Centroid
                coord[0] = 0.5;

                this->coord_.push_back ( coord );
            }
            else
            {
                assert ( this->fe_deg_ > 0 );

                // Lexicographical ordering

                const DataType offset = ( 1.0 / this->fe_deg_ );

                this->coord_.clear ( );
                const int nb_dof_on_cell = ( this->fe_deg_ + 1 );
                this->coord_.resize ( nb_dof_on_cell );

                const int nb_dof_line = this->fe_deg_ + 1;

                // Filling coord vector by lexicographical strategy

                for ( int j = 0; j < nb_dof_line; ++j )
                {
                    Coord coord;
                    coord.resize ( 1 );

                    coord[0] = j*offset;

                    this->coord_[ij2ind ( 0, j )] = coord;
                }
            }
        }

        template<class DataType>
        int FELagrangeLine<DataType>::ij2ind ( int i, int j ) const
        {

            return (j );
        }

        /// \details Every degree of a used lagrangian polynomial has to satisfy the
        ///          condition degree < this->fe_deg_. All possible combinations are multiplied
        ///          and considered in the sum.

        template<class DataType>
        void FELagrangeLine<DataType>::N ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            if ( this->fe_deg_ > 0 )
            {
                for ( int j = 0; j <= this->fe_deg_; ++j )
                {
                    weight[ij2ind ( 0, j )] = this->lp_.poly ( this->fe_deg_, j, pt[0] );
                }
            }
            else
                weight[ij2ind ( 0, 0 )] = 1.0;
        }

        /// \details Every degree of a used lagrangian polynomial has to satisfy the
        ///          condition degree < this->fe_deg_. All possible combinations are multiplied
        ///          and considered in the sum, w.r.t. the derivatives for the x - variable.

        template<class DataType>
        void FELagrangeLine<DataType>::N_x ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            if ( this->fe_deg_ > 0 )
            {
                for ( int j = 0; j <= this->fe_deg_; ++j )

                {
                    weight[ij2ind ( 0, j )] = this->lp_.poly_x ( this->fe_deg_, j, pt[0] );
                }
            }
            else
                weight[ij2ind ( 0, 0 )] = 0.0;
        }

        /// \details There are no y - derivatives in a 1D case

        template<class DataType>
        void FELagrangeLine<DataType>::N_y ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details There are no z - derivatives in a 1D case

        template<class DataType>
        void FELagrangeLine<DataType>::N_z ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details Every degree of a used lagrangian polynomial has to satisfy the
        ///          condition degree < this->fe_deg_. All possible combinations are multiplied
        ///          and considered in the sum, w.r.t. the second derivatives for the xx - variable.

        template<class DataType>
        void FELagrangeLine<DataType>::N_xx ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            if ( this->fe_deg_ > 0 )
            {
                for ( int j = 0; j <= this->fe_deg_; ++j )
                {
                    weight[ij2ind ( 0, j )] = this->lp_.poly_xx ( this->fe_deg_, j, pt[0] );

                }
            }
            else
                weight[ij2ind ( 0, 0 )] = 0.0;
        }

        /// \details There are no y - derivatives in a 1D case

        template<class DataType>
        void FELagrangeLine<DataType>::N_xy ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details There are no z - derivatives in a 1D case

        template<class DataType>
        void FELagrangeLine<DataType>::N_xz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details There are no y - derivatives in a 1D case

        template<class DataType>
        void FELagrangeLine<DataType>::N_yy ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details There are no z - derivatives in a 1D case

        template<class DataType>
        void FELagrangeLine<DataType>::N_yz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        /// \details There are no z - derivatives in a 1D case

        template<class DataType>
        void FELagrangeLine<DataType>::N_zz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            weight.assign ( weight.size ( ), 0. );
        }

        template class FELagrangeLine<double>;
        template class FELagrangeLine<float>;

    } // namespace doffem
} // namespace hiflow
