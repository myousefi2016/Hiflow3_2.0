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

#include "felagrange_pyr.h"
#include <cassert>
#include <cmath>
#include <iomanip>

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        FELagrangePyr<DataType>::FELagrangePyr ( )
        {
            this->my_id_ = FEType<DataType>::LAGRANGE_PYR;

            // initialize reference cell

            assert ( this->ref_cell_ == NULL );

            if ( this->fe_deg_ > 2 )
            {
                std::cerr << "Only support up to degree = 2 !" << std::endl;
                ;
                exit ( -1 );
            }

            this->ref_cell_ = &( mesh::CellType::get_instance ( mesh::CellType::PYRAMID ) );
        }

        template<class DataType>
        FELagrangePyr<DataType>::~FELagrangePyr ( )
        {
        }

        template<class DataType>
        void FELagrangePyr<DataType>::init_coord ( )
        {
            // set topological degree

            this->tdim_ = 3;

            // Lexicographical ordering

            if ( this->fe_deg_ == 0 )
            {
                this->coord_.clear ( );

                Coord coord;
                coord.resize ( 3 );

                // Centroid
                coord[0] = 0.5;
                coord[1] = 0.5;
                coord[2] = 0.25;

                this->coord_.push_back ( coord );
            }
            else
            {
                assert ( this->fe_deg_ > 0 );

                const DataType offset = ( 1.0 / this->fe_deg_ );

                this->coord_.clear ( );

                int nb_dof_on_cell;

                if ( this->fe_deg_ == 1 )
                {
                    nb_dof_on_cell = 5;
                }
                else if ( this->fe_deg_ == 2 )
                {
                    nb_dof_on_cell = 14;
                }
                else
                {
                    std::cerr << "Unexpected finite element degree " << this->fe_deg_ << std::endl;
                    exit ( -1 );
                }

                this->coord_.resize ( nb_dof_on_cell );

                const int nb_dof_line = this->fe_deg_ + 1;

                for ( int k = 0; k < nb_dof_line; ++k )
                { // z axis
                    for ( int j = 0; j < nb_dof_line - k; ++j )
                    { // y axis
                        for ( int i = 0; i < nb_dof_line - k; ++i ) // x axis
                        {
                            Coord coord;
                            coord.resize ( 3 );

                            coord[0] = k * offset * 0.5 + i * offset;
                            coord[1] = k * offset * 0.5 + j * offset;
                            coord[2] = k * offset;

                            this->coord_[ijk2ind ( i, j, k )] = coord;
                        }
                    }
                }
            }
        }

        template<class DataType>
        int FELagrangePyr<DataType>::ijk2ind ( int i, int j, int k ) const
        {
            // x component = i, y component = j, z component = k

            int offset = 0;
            const int nb_dof_line = this->fe_deg_ + 1;

            // First: offset z axis

            for ( int m = 0; m < k; ++m )
            {
                const int help = nb_dof_line - m;
                offset += help * help;
            }

            // Second: increasing offset by y axis on current z axis

            for ( int n = 0; n < j; ++n )
                offset += nb_dof_line - k;

            return (i + offset );
        }

        template<class DataType>
        void FELagrangePyr<DataType>::N ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            std::setprecision ( 2 * sizeof (DataType ) );

            if ( this->fe_deg_ == 0 )
            {

                weight[0] = 1.0;

            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = -0.5 * z + ( 0.25 * x - 0.25 ) * ( y - z - 1.0 );
                        weight[1] = -0.5 * z + ( 0.25 * x + 0.25 ) * ( -y + z + 1.0 );
                        weight[2] = ( -0.25 * x + 0.25 ) * ( y - z + 1.0 );
                        weight[3] = ( 0.25 * x + 0.25 ) * ( y - z + 1.0 );
                        weight[4] = z;
                        break;
                    case 2:
                        weight[0] = 0.25 * z * ( x + y - 2.0 ) + ( 0.25 * x - 0.25 ) * ( y - z - 1.0 );
                        weight[1] = -0.25 * z * ( x + y + 2.0 ) + ( 0.25 * x + 0.25 ) * ( -y + z + 1.0 );
                        weight[2] = -0.25 * z * ( x + y ) + ( -0.25 * x + 0.25 ) * ( y - z + 1.0 );
                        weight[3] = 0.25 * z * ( x + y ) + ( 0.25 * x + 0.25 ) * ( y - z + 1.0 );
                        weight[4] = z;
                        break;
                    case 3:
                        weight[0] = ( 0.25 * x - 0.25 ) * ( y + z - 1.0 );
                        weight[1] = ( 0.25 * x + 0.25 ) * ( -y - z + 1.0 );
                        weight[2] = -0.5 * x * z + ( -0.25 * x + 0.25 ) * ( y - z + 1.0 );
                        weight[3] = 0.5 * x * z + ( 0.25 * x + 0.25 ) * ( y - z + 1.0 );
                        weight[4] = z;
                        break;
                    case 4:
                        weight[0] = 0.25 * z * ( x - y - 2.0 ) + ( 0.25 * x - 0.25 ) * ( y - z - 1.0 );
                        weight[1] = -0.25 * z * ( x - y + 2.0 ) + ( 0.25 * x + 0.25 ) * ( -y + z + 1.0 );
                        weight[2] = -0.25 * z * ( x - y ) + ( -0.25 * x + 0.25 ) * ( y - z + 1.0 );
                        weight[3] = 0.25 * z * ( x - y ) + ( 0.25 * x + 0.25 ) * ( y - z + 1.0 );
                        weight[4] = z;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 1.0 * weight[0];
                weight[1] = 1.0 * weight[1];
                weight[2] = 1.0 * weight[2];
                weight[3] = 1.0 * weight[3];
                weight[4] = 1.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0.5 * z * ( x + y ) - 0.5 * ( 0.25 * x - 0.25 * z ) * ( y - z ) * ( 4.0 * z - ( -x + z + 1.0 ) * ( -y + z + 1.0 ) )
                                - 0.125 * ( x + z ) * ( y - z ) * ( x + z - 1.0 ) * ( -y + z + 1.0 );
                        weight[1] = -0.5 * pow ( x, 2 ) * pow ( y, 2 ) + 0.5 * pow ( x, 2 ) * y * z + 0.5 * pow ( x, 2 ) * y - 1.0 * pow ( x, 2 ) * z
                                - 0.5 * pow ( y, 2 ) * z + 0.5 * pow ( y, 2 ) + 0.5 * y * pow ( z, 2 ) - 0.5 * y;
                        weight[2] = -0.5 * z * ( x - y ) - 0.5 * ( 0.25 * x - 0.25 * z ) * ( y - z ) * ( x - z + 1.0 ) * ( -y + z + 1.0 )
                                + 0.125 * ( x + z ) * ( y - z ) * ( 4.0 * z - ( x + z + 1.0 ) * ( -y + z + 1.0 ) );
                        weight[3] = ( 0.5 * y - 0.5 * z + 0.5 ) * ( 0.5 * x * ( y - 1.0 ) * ( -x + z + 1.0 ) - 0.5 * x * ( y - 1.0 ) * ( x + z - 1.0 )
                                - 0.5 * z * ( 2.0 * y + 1.0 ) + 0.5 * z );
                        weight[4] = ( ( z * ( x - y + z + 1.0 ) + ( x + 1.0 ) * ( y - 1.0 ) ) * ( x + z - 1.0 ) + ( z * ( x + y - z - 1.0 )
                                + ( x - 1.0 ) * ( y - 1.0 ) ) * ( x - z + 1.0 ) ) * ( 0.5 * y - 0.5 * z + 0.5 );
                        weight[5] = ( -0.25 * y + 0.25 * z - 0.25 ) * ( x * ( y - 1.0 ) * ( x - z + 1.0 ) + x * ( y - 1.0 ) * ( x + z + 1.0 ) + z * ( 2.0 * y + 1.0 ) - z );
                        weight[6] = ( -0.5 * ( x - z ) * ( 0.25 * y - 0.25 * z ) * ( -x + z + 1.0 ) + 0.125 * ( x + z ) * ( y - z ) * ( x + z - 1.0 ) ) * ( y - z + 1.0 );
                        weight[7] = -0.5 * y * ( ( x - 1.0 ) * ( x - z + 1.0 ) + ( x + 1.0 ) * ( x + z - 1.0 ) ) * ( 0.5 * y - 0.5 * z + 0.5 );
                        weight[8] = ( ( 0.25 * x - 0.25 * z ) * ( y - z ) * ( x - z + 1.0 ) + ( x + z ) * ( 0.25 * y - 0.25 * z ) * ( x + z + 1.0 ) ) * ( 0.5 * y - 0.5 * z + 0.5 );
                        weight[9] = 1.0 * z * ( x * y - x * z - x - y - z + 1 );
                        weight[10] = 1.0 * z * ( -x * y + x * z + x - y - z + 1 );
                        weight[11] = 1.0 * z * ( -x + 1 ) * ( y - z + 1.0 );
                        weight[12] = 0.5 * z * ( 2 * x + 2.0 ) * ( y - z + 1.0 );
                        weight[13] = z * ( 2.0 * z - 1.0 );
                        break;
                    case 2:
                        weight[0] = ( x + z ) * ( -0.125 * ( y - z ) * ( -y + z + 1.0 ) + 0.125 * ( y + z ) * ( y + z - 1.0 ) ) * ( x + z - 1.0 );
                        weight[1] = ( x + z - 1 ) * ( 0.25 * y * ( x + 1.0 ) * ( -y + z + 1.0 ) - 0.25 * y * ( x + 1.0 ) * ( y + z - 1.0 )
                                - 0.25 * z * ( 2.0 * x + 1.0 ) + 0.25 * z );
                        weight[2] = -0.5 * z * ( x - y ) + 0.125 * ( x + z ) * ( y - z ) * ( 4.0 * z - ( x + z + 1.0 ) * ( -y + z + 1.0 ) )
                                + 0.125 * ( x + z ) * ( y + z ) * ( x + z + 1.0 ) * ( y + z - 1.0 );
                        weight[3] = -x * ( 0.5 * ( y - 1.0 ) * ( x + z - 1.0 ) * ( 0.5 * y - 0.5 * z + 0.5 ) + 0.25 * ( y + 1.0 ) * ( x + z - 1 ) * ( y + z - 1.0 ) );
                        weight[4] = ( 0.5 * ( z * ( x - y + z + 1.0 ) + ( x + 1.0 ) * ( y - 1.0 ) ) * ( y - z + 1.0 )
                                - 0.5 * ( z * ( x + y + z + 1.0 ) - ( x + 1.0 ) * ( y + 1.0 ) ) * ( y + z - 1.0 ) ) * ( x + z - 1.0 );
                        weight[5] = -0.5 * pow ( x, 2 ) * pow ( y, 2 ) - 0.5 * pow ( x, 2 ) * z + 0.5 * pow ( x, 2 ) - 0.5 * x * pow ( y, 2 ) * z
                                - 0.5 * x * pow ( y, 2 ) - 0.5 * x * pow ( z, 2 ) + 0.5 * x - 1.0 * pow ( y, 2 ) * z;
                        weight[6] = 0.125 * ( x + z ) * ( ( y - z ) * ( y - z + 1.0 ) + ( y + z ) * ( y + z + 1.0 ) ) * ( x + z - 1.0 );
                        weight[7] = -0.5 * y * ( x + 1.0 ) * ( x + z - 1.0 ) * ( 0.5 * y - 0.5 * z + 0.5 ) - 0.25 * ( x + z - 1 ) * ( y * ( x + 1.0 ) * ( y + z + 1.0 )
                                + z * ( 2.0 * x - 1.0 ) + z );
                        weight[8] = -0.5 * z * ( x + y ) + 0.5 * ( x + z ) * ( 0.25 * y - 0.25 * z ) * ( x + z + 1.0 ) * ( y - z + 1.0 )
                                - 0.125 * ( x + z ) * ( y + z ) * ( 4.0 * z - ( x + z + 1.0 ) * ( y + z + 1.0 ) );
                        weight[9] = 1.0 * z * ( y - 1 ) * ( x + z - 1.0 );
                        weight[10] = 1.0 * z * ( -x * y + x - y * z - y - z + 1 );
                        weight[11] = -0.5 * z * ( 2 * y + 2.0 ) * ( x + z - 1.0 );
                        weight[12] = 1.0 * z * ( x * y + x + y * z + y - z + 1 );
                        weight[13] = z * ( 2.0 * z - 1.0 );
                        break;
                    case 3:
                        weight[0] = ( y + z ) * ( -0.125 * ( x - z ) * ( -x + z + 1.0 ) + 0.125 * ( x + z ) * ( x + z - 1.0 ) ) * ( y + z - 1.0 );
                        weight[1] = -y * ( 0.5 * ( x - 1.0 ) * ( 0.5 * x - 0.5 * z + 0.5 ) + 0.25 * ( x + 1.0 ) * ( x + z - 1 ) ) * ( y + z - 1.0 );
                        weight[2] = ( y + z ) * ( 0.5 * ( 0.25 * x - 0.25 * z ) * ( x - z + 1.0 ) + 0.125 * ( x + z ) * ( x + z + 1.0 ) ) * ( y + z - 1.0 );
                        weight[3] = -0.25 * x * ( y + 1.0 ) * ( x + z - 1 ) * ( y + z - 1.0 )
                                + 0.25 * ( y + z - 1 ) * ( x * ( y + 1.0 ) * ( -x + z + 1.0 ) - z * ( 2.0 * y + 1.0 ) + z );
                        weight[4] = ( 0.5 * ( z * ( -x + y + z + 1.0 ) + ( x - 1.0 ) * ( y + 1.0 ) ) * ( x - z + 1.0 ) - 0.5 * ( z * ( x + y + z + 1.0 )
                                - ( x + 1.0 ) * ( y + 1.0 ) ) * ( x + z - 1.0 ) ) * ( y + z - 1.0 );
                        weight[5] = -0.5 * x * ( y + 1.0 ) * ( 0.5 * x - 0.5 * z + 0.5 ) * ( y + z - 1.0 ) - 0.25 * ( y + z - 1 ) * ( x * ( y + 1.0 ) * ( x + z + 1.0 )
                                + z * ( 2.0 * y + 1.0 ) - z );
                        weight[6] = 0.5 * z * ( x - y ) + 0.5 * ( 0.25 * x - 0.25 * z ) * ( y + z ) * ( 4.0 * z - ( -x + z + 1.0 ) * ( y + z + 1.0 ) )
                                + 0.125 * ( x + z ) * ( y + z ) * ( x + z - 1.0 ) * ( y + z + 1.0 );
                        weight[7] = -0.5 * pow ( x, 2 ) * pow ( y, 2 ) - 0.5 * pow ( x, 2 ) * y * z - 0.5 * pow ( x, 2 ) * y - 1.0 * pow ( x, 2 ) * z
                                - 0.5 * pow ( y, 2 ) * z + 0.5 * pow ( y, 2 ) - 0.5 * y * pow ( z, 2 ) + 0.5 * y;
                        weight[8] = -0.5 * z * ( x + y ) + 0.125 * ( x - z ) * ( y + z ) * ( x - z + 1.0 ) * ( y + z + 1.0 )
                                - 0.125 * ( x + z ) * ( y + z ) * ( 4.0 * z - ( x + z + 1.0 ) * ( y + z + 1.0 ) );
                        weight[9] = 1.0 * z * ( x - 1 ) * ( y + z - 1.0 );
                        weight[10] = -0.5 * z * ( 2 * x + 2.0 ) * ( y + z - 1.0 );
                        weight[11] = 1.0 * z * ( -x * y - x * z - x + y - z + 1 );
                        weight[12] = 1.0 * z * ( x * y + x * z + x + y - z + 1 );
                        weight[13] = z * ( 2.0 * z - 1.0 );
                        break;
                    case 4:
                        weight[0] = 0.5 * z * ( x + y ) - 0.5 * ( 0.25 * x - 0.25 * z ) * ( y - z ) * ( 4.0 * z - ( -x + z + 1.0 ) * ( -y + z + 1.0 ) )
                                - 0.125 * ( x - z ) * ( y + z ) * ( -x + z + 1.0 ) * ( y + z - 1.0 );
                        weight[1] = ( -0.25 * x + 0.25 * z - 0.25 ) * ( -y * ( x - 1.0 ) * ( -y + z + 1.0 ) + y * ( x - 1.0 ) * ( y + z - 1.0 ) + z * ( 2.0 * x - 1.0 ) + z );
                        weight[2] = ( 0.25 * x - 0.25 * z ) * ( -0.5 * ( y - z ) * ( -y + z + 1.0 ) + 0.5 * ( y + z ) * ( y + z - 1.0 ) ) * ( x - z + 1.0 );
                        weight[3] = -0.5 * pow ( x, 2 ) * pow ( y, 2 ) - 0.5 * pow ( x, 2 ) * z + 0.5 * pow ( x, 2 ) + 0.5 * x * pow ( y, 2 ) * z + 0.5 * x * pow ( y, 2 )
                                + 0.5 * x * pow ( z, 2 ) - 0.5 * x - 1.0 * pow ( y, 2 ) * z;
                        weight[4] = ( ( z * ( -x + y + z + 1.0 ) + ( x - 1.0 ) * ( y + 1.0 ) ) * ( y + z - 1.0 ) + ( z * ( x + y - z - 1.0 )
                                + ( x - 1.0 ) * ( y - 1.0 ) ) * ( y - z + 1.0 ) ) * ( 0.5 * x - 0.5 * z + 0.5 );
                        weight[5] = -0.5 * x * ( ( y - 1.0 ) * ( x - z + 1.0 ) * ( 0.5 * y - 0.5 * z + 0.5 ) + ( y + 1.0 ) * ( 0.5 * x - 0.5 * z + 0.5 ) * ( y + z - 1.0 ) );
                        weight[6] = 0.5 * z * ( x - y ) + 0.5 * ( 0.25 * x - 0.25 * z ) * ( y + z ) * ( 4.0 * z - ( -x + z + 1.0 ) * ( y + z + 1.0 ) )
                                - 0.5 * ( x - z ) * ( 0.25 * y - 0.25 * z ) * ( -x + z + 1.0 ) * ( y - z + 1.0 );
                        weight[7] = -0.5 * y * ( x - 1.0 ) * ( x - z + 1.0 ) * ( 0.5 * y - 0.5 * z + 0.5 )
                                - 0.5 * ( 0.5 * x - 0.5 * z + 0.5 ) * ( y * ( x - 1.0 ) * ( y + z + 1.0 ) + z * ( 2.0 * x + 1.0 ) - z );
                        weight[8] = ( 0.5 * ( 0.25 * x - 0.25 * z ) * ( y - z ) * ( y - z + 1.0 ) + 0.125 * ( x - z ) * ( y + z ) * ( y + z + 1.0 ) ) * ( x - z + 1.0 );
                        weight[9] = 1.0 * z * ( x * y - x - y * z - y - z + 1 );
                        weight[10] = 1.0 * z * ( -y + 1 ) * ( x - z + 1.0 );
                        weight[11] = 1.0 * z * ( -x * y - x + y * z + y - z + 1 );
                        weight[12] = 0.5 * z * ( 2 * y + 2.0 ) * ( x - z + 1.0 );
                        weight[13] = z * ( 2.0 * z - 1.0 );
                        break;
                    default:
                        std::cout << x << " " << y << std::endl;
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 1.0 * weight[0];
                weight[1] = 1.0 * weight[1];
                weight[2] = 1.0 * weight[2];
                weight[3] = 1.0 * weight[3];
                weight[4] = 1.0 * weight[4];
                weight[5] = 1.0 * weight[5];
                weight[6] = 1.0 * weight[6];
                weight[7] = 1.0 * weight[7];
                weight[8] = 1.0 * weight[8];
                weight[9] = 1.0 * weight[9];
                weight[10] = 1.0 * weight[10];
                weight[11] = 1.0 * weight[11];
                weight[12] = 1.0 * weight[12];
                weight[13] = 1.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_x ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {
                weight[0] = 0.0;
            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0.25 * y - 0.25 * z - 0.25;
                        weight[1] = -0.25 * y + 0.25 * z + 0.25;
                        weight[2] = -0.25 * y + 0.25 * z - 0.25;
                        weight[3] = 0.25 * y - 0.25 * z + 0.25;
                        weight[4] = 0;
                        break;
                    case 2:
                        weight[0] = 0.25 * y - 0.25;
                        weight[1] = -0.25 * y + 0.25;
                        weight[2] = -0.25 * y - 0.25;
                        weight[3] = 0.25 * y + 0.25;
                        weight[4] = 0;
                        break;
                    case 3:
                        weight[0] = 0.25 * y + 0.25 * z - 0.25;
                        weight[1] = -0.25 * y - 0.25 * z + 0.25;
                        weight[2] = -0.25 * y - 0.25 * z - 0.25;
                        weight[3] = 0.25 * y + 0.25 * z + 0.25;
                        weight[4] = 0;
                        break;
                    case 4:
                        weight[0] = 0.25 * y - 0.25;
                        weight[1] = -0.25 * y + 0.25;
                        weight[2] = -0.25 * y - 0.25;
                        weight[3] = 0.25 * y + 0.25;
                        weight[4] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 2.0 * weight[0];
                weight[1] = 2.0 * weight[1];
                weight[2] = 2.0 * weight[2];
                weight[3] = 2.0 * weight[3];
                weight[4] = 2.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0.5 * x * pow ( y, 2 ) - 1.0 * x * y * z - 0.5 * x * y + 0.5 * x * pow ( z, 2 ) + 0.5 * x * z - 0.25 * pow ( y, 2 )
                                + 0.25 * y + 0.25 * pow ( z, 2 ) + 0.25 * z;
                        weight[1] = x * ( -1.0 * pow ( y, 2 ) + 1.0 * y * z + 1.0 * y - 2.0 * z );
                        weight[2] = 0.5 * x * pow ( y, 2 ) - 1.0 * x * y * z - 0.5 * x * y + 0.5 * x * pow ( z, 2 ) + 0.5 * x * z + 0.25 * pow ( y, 2 )
                                - 0.25 * y - 0.25 * pow ( z, 2 ) - 0.25 * z;
                        weight[3] = -0.5 * ( 4 * x - 2.0 ) * ( y - 1.0 ) * ( 0.5 * y - 0.5 * z + 0.5 );
                        weight[4] = x * ( 2.0 * pow ( y, 2 ) - 2.0 * pow ( z, 2 ) + 4.0 * z - 2.0 );
                        weight[5] = -0.5 * ( 4 * x + 2.0 ) * ( y - 1.0 ) * ( 0.5 * y - 0.5 * z + 0.5 );
                        weight[6] = ( y - z + 1.0 ) * ( 0.5 * ( x - z ) * ( 0.25 * y - 0.25 * z ) + 0.125 * ( x + z ) * ( y - z )
                                - 0.5 * ( 0.25 * y - 0.25 * z ) * ( -x + z + 1.0 ) + 0.125 * ( y - z ) * ( x + z - 1.0 ) );
                        weight[7] = 1.0 * x * y * ( -y + z - 1 );
                        weight[8] = ( y - z + 1.0 ) * ( 0.5 * ( 0.25 * x - 0.25 * z ) * ( y - z ) + 0.5 * ( x + z ) * ( 0.25 * y - 0.25 * z )
                                + 0.5 * ( 0.25 * y - 0.25 * z ) * ( x + z + 1.0 ) + 0.125 * ( y - z ) * ( x - z + 1.0 ) );
                        weight[9] = 1.0 * z * ( y - z - 1 );
                        weight[10] = 1.0 * z * ( -y + z + 1.0 );
                        weight[11] = 1.0 * z * ( -y + z - 1 );
                        weight[12] = 1.0 * z * ( y - z + 1.0 );
                        weight[13] = 0;
                        break;
                    case 2:
                        weight[0] = 0.5 * x * pow ( y, 2 ) - 0.5 * x * y + 0.5 * x * pow ( z, 2 ) + 0.5 * pow ( y, 2 ) * z - 0.25 * pow ( y, 2 )
                                - 0.5 * y * z + 0.25 * y + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 );
                        weight[1] = -1.0 * x * pow ( y, 2 ) + 1.0 * x * y - 1.0 * x * z - 0.5 * pow ( y, 2 ) * z + 0.5 * y * z - 0.5 * pow ( z, 2 ) + 0.5 * z;
                        weight[2] = 0.5 * x * pow ( y, 2 ) - 0.5 * x * y + 0.5 * x * pow ( z, 2 ) + 0.5 * pow ( y, 2 ) * z + 0.25 * pow ( y, 2 ) - 0.25 * y
                                + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 ) - 0.5 * z;
                        weight[3] = -1.0 * x * pow ( y, 2 ) - 1.0 * x * z + 1.0 * x - 0.5 * pow ( y, 2 ) * z + 0.5 * pow ( y, 2 ) - 0.5 * pow ( z, 2 ) + 1.0 * z - 0.5;
                        weight[4] = 2.0 * x * pow ( y, 2 ) - 2.0 * x * pow ( z, 2 ) + 4.0 * x * z - 2.0 * x - 2.0 * pow ( z, 3 ) + 3.0 * pow ( z, 2 ) - 1.0 * z;
                        weight[5] = -( 0.5 * ( y - 1.0 ) * ( 0.5 * y - 0.5 * z + 0.5 ) + 0.25 * ( y + 1.0 ) * ( y + z - 1 ) ) * ( 2 * x + z + 1.0 );
                        weight[6] = 0.5 * x * pow ( y, 2 ) + 0.5 * x * y + 0.5 * x * pow ( z, 2 ) + 0.5 * pow ( y, 2 ) * z - 0.25 * pow ( y, 2 ) + 0.5 * y * z
                                - 0.25 * y + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 );
                        weight[7] = -1.0 * x * pow ( y, 2 ) - 1.0 * x * y - 1.0 * x * z - 0.5 * pow ( y, 2 ) * z - 0.5 * y * z - 0.5 * pow ( z, 2 ) + 0.5 * z;
                        weight[8] = 0.5 * x * pow ( y, 2 ) + 0.5 * x * y + 0.5 * x * pow ( z, 2 ) + 0.5 * pow ( y, 2 ) * z + 0.25 * pow ( y, 2 ) + 0.25 * y
                                + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 ) - 0.5 * z;
                        weight[9] = 1.0 * z * ( y - 1 );
                        weight[10] = 1.0 * z * ( -y + 1 );
                        weight[11] = -1.0 * z * ( y + 1 );
                        weight[12] = 1.0 * z * ( y + 1 );
                        weight[13] = 0;
                        break;
                    case 3:
                        weight[0] = ( 0.5 * x - 0.25 ) * ( y + z ) * ( y + z - 1.0 );
                        weight[1] = 1.0 * x * y * ( -y - z + 1 );
                        weight[2] = ( 0.5 * x + 0.25 ) * ( y + z ) * ( y + z - 1.0 );
                        weight[3] = -0.25 * ( y + 1 ) * ( x * ( y + z - 1.0 ) - ( -2 * x + z + 1.0 ) * ( y + z - 1 ) + ( x + z - 1 ) * ( y + z - 1.0 ) );
                        weight[4] = x * ( 2.0 * pow ( y, 2 ) - 2.0 * pow ( z, 2 ) + 4.0 * z - 2.0 );
                        weight[5] = -1.0 * x * pow ( y, 2 ) - 1.0 * x * y * z - 1.0 * x * z + 1.0 * x - 0.5 * pow ( y, 2 ) - 0.5 * y * z - 0.5 * z + 0.5;
                        weight[6] = 0.5 * x * pow ( y, 2 ) + 1.0 * x * y * z + 0.5 * x * y + 0.5 * x * pow ( z, 2 ) + 0.5 * x * z - 0.25 * pow ( y, 2 )
                                - 0.25 * y + 0.25 * pow ( z, 2 ) + 0.25 * z;
                        weight[7] = -x * ( 1.0 * pow ( y, 2 ) + 1.0 * y * z + 1.0 * y + 2.0 * z );
                        weight[8] = 0.5 * x * pow ( y, 2 ) + 1.0 * x * y * z + 0.5 * x * y + 0.5 * x * pow ( z, 2 ) + 0.5 * x * z + 0.25 * pow ( y, 2 )
                                + 0.25 * y - 0.25 * pow ( z, 2 ) - 0.25 * z;
                        weight[9] = 1.0 * z * ( y + z - 1.0 );
                        weight[10] = 1.0 * z * ( -y - z + 1 );
                        weight[11] = -1.0 * z * ( y + z + 1.0 );
                        weight[12] = 1.0 * z * ( y + z + 1.0 );
                        weight[13] = 0;
                        break;
                    case 4:
                        weight[0] = 0.5 * x * pow ( y, 2 ) - 0.5 * x * y + 0.5 * x * pow ( z, 2 ) - 0.5 * pow ( y, 2 ) * z - 0.25 * pow ( y, 2 ) + 0.25 * y
                                - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 ) + 0.5 * z;
                        weight[1] = -1.0 * x * pow ( y, 2 ) + 1.0 * x * y - 1.0 * x * z + 0.5 * pow ( y, 2 ) * z - 0.5 * y * z + 0.5 * pow ( z, 2 ) - 0.5 * z;
                        weight[2] = 0.5 * x * pow ( y, 2 ) - 0.5 * x * y + 0.5 * x * pow ( z, 2 ) - 0.5 * pow ( y, 2 ) * z + 0.25 * pow ( y, 2 ) + 0.5 * y * z
                                - 0.25 * y - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 );
                        weight[3] = ( -0.5 * ( y - 1.0 ) * ( 0.5 * y - 0.5 * z + 0.5 ) - 0.25 * ( y + 1.0 ) * ( y + z - 1 ) ) * ( 2 * x - z - 1.0 );
                        weight[4] = 2.0 * x * pow ( y, 2 ) - 2.0 * x * pow ( z, 2 ) + 4.0 * x * z - 2.0 * x + 2.0 * pow ( z, 3 ) - 3.0 * pow ( z, 2 ) + 1.0 * z;
                        weight[5] = -1.0 * x * pow ( y, 2 ) - 1.0 * x * z + 1.0 * x + 0.5 * pow ( y, 2 ) * z - 0.5 * pow ( y, 2 ) + 0.5 * pow ( z, 2 ) - 1.0 * z + 0.5;
                        weight[6] = 0.5 * x * pow ( y, 2 ) + 0.5 * x * y + 0.5 * x * pow ( z, 2 ) - 0.5 * pow ( y, 2 ) * z - 0.25 * pow ( y, 2 ) - 0.25 * y
                                - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 ) + 0.5 * z;
                        weight[7] = -1.0 * x * pow ( y, 2 ) - 1.0 * x * y - 1.0 * x * z + 0.5 * pow ( y, 2 ) * z + 0.5 * y * z + 0.5 * pow ( z, 2 ) - 0.5 * z;
                        weight[8] = 0.5 * x * pow ( y, 2 ) + 0.5 * x * y + 0.5 * x * pow ( z, 2 ) - 0.5 * pow ( y, 2 ) * z + 0.25 * pow ( y, 2 ) - 0.5 * y * z
                                + 0.25 * y - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 );
                        weight[9] = 1.0 * z * ( y - 1 );
                        weight[10] = 1.0 * z * ( -y + 1 );
                        weight[11] = -1.0 * z * ( y + 1 );
                        weight[12] = 1.0 * z * ( y + 1 );
                        weight[13] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                if ( xx == 0 && yy == 0 )
                {
                    weight[2] = 0.5 * x * pow ( y, 2 ) - 0.5 * x * y + 0.5 * x * pow ( z, 2 ) - 0.5 * pow ( y, 2 ) * z + 0.25 * pow ( y, 2 ) + 0.5 * y * z
                            - 0.25 * y - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 );
                    weight[3] = -0.25 * ( y + 1 ) * ( x * ( y + z - 1.0 ) - ( -2 * x + z + 1.0 ) * ( y + z - 1 ) + ( x + z - 1 ) * ( y + z - 1.0 ) );
                    weight[5] = -1.0 * x * pow ( y, 2 ) - 1.0 * x * y * z - 1.0 * x * z + 1.0 * x - 0.5 * pow ( y, 2 ) - 0.5 * y * z - 0.5 * z + 0.5;
                    weight[8] = 0.5 * x * pow ( y, 2 ) + 0.5 * x * y + 0.5 * x * pow ( z, 2 ) - 0.5 * pow ( y, 2 ) * z + 0.25 * pow ( y, 2 ) - 0.5 * y * z
                            + 0.25 * y - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 );
                }

                //rescaling due to changing coordinates
                weight[0] = 2.0 * weight[0];
                weight[1] = 2.0 * weight[1];
                weight[2] = 2.0 * weight[2];
                weight[3] = 2.0 * weight[3];
                weight[4] = 2.0 * weight[4];
                weight[5] = 2.0 * weight[5];
                weight[6] = 2.0 * weight[6];
                weight[7] = 2.0 * weight[7];
                weight[8] = 2.0 * weight[8];
                weight[9] = 2.0 * weight[9];
                weight[10] = 2.0 * weight[10];
                weight[11] = 2.0 * weight[11];
                weight[12] = 2.0 * weight[12];
                weight[13] = 2.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_y ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {
                weight[0] = 0.0;
            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx >= yy && xx <= -yy )
                {
                    indicator = 1;
                }
                else if ( xx > yy && xx > -yy )
                {
                    indicator = 2;
                }
                else if ( xx <= yy && xx >= -yy )
                {
                    indicator = 3;
                }
                else if ( xx < yy && xx < -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0.25 * x - 0.25;
                        weight[1] = -0.25 * x - 0.25;
                        weight[2] = -0.25 * x + 0.25;
                        weight[3] = 0.25 * x + 0.25;
                        weight[4] = 0;
                        break;
                    case 2:
                        weight[0] = 0.25 * x + 0.25 * z - 0.25;
                        weight[1] = -0.25 * x - 0.25 * z - 0.25;
                        weight[2] = -0.25 * x - 0.25 * z + 0.25;
                        weight[3] = 0.25 * x + 0.25 * z + 0.25;
                        weight[4] = 0;
                        break;
                    case 3:
                        weight[0] = 0.25 * x - 0.25;
                        weight[1] = -0.25 * x - 0.25;
                        weight[2] = -0.25 * x + 0.25;
                        weight[3] = 0.25 * x + 0.25;
                        weight[4] = 0;
                        break;
                    case 4:
                        weight[0] = 0.25 * x - 0.25 * z - 0.25;
                        weight[1] = -0.25 * x + 0.25 * z - 0.25;
                        weight[2] = -0.25 * x + 0.25 * z + 0.25;
                        weight[3] = 0.25 * x - 0.25 * z + 0.25;
                        weight[4] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 2.0 * weight[0];
                weight[1] = 2.0 * weight[1];
                weight[2] = 2.0 * weight[2];
                weight[3] = 2.0 * weight[3];
                weight[4] = 2.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx >= yy && xx <= -yy )
                {
                    indicator = 1;
                }
                else if ( xx > yy && xx > -yy )
                {
                    indicator = 2;
                }
                else if ( xx <= yy && xx >= -yy )
                {
                    indicator = 3;
                }
                else if ( xx < yy && xx < -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0.5 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) - 0.5 * x * y + 0.25 * x
                                + 0.5 * y * pow ( z, 2 ) - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 ) + 0.5 * z;
                        weight[1] = ( -0.5 * ( x - 1.0 ) * ( 0.5 * x - 0.5 * z + 0.5 ) - 0.25 * ( x + 1.0 ) * ( x + z - 1 ) ) * ( 2 * y - z - 1.0 );
                        weight[2] = 0.5 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) + 0.5 * x * y - 0.25 * x
                                + 0.5 * y * pow ( z, 2 ) - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 ) + 0.5 * z;
                        weight[3] = -1.0 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z + 1.0 * x * y - 0.5 * x * z - 1.0 * y * z + 0.5 * pow ( z, 2 ) - 0.5 * z;
                        weight[4] = 2.0 * pow ( x, 2 ) * y - 2.0 * y * pow ( z, 2 ) + 4.0 * y * z - 2.0 * y + 2.0 * pow ( z, 3 ) - 3.0 * pow ( z, 2 ) + 1.0 * z;
                        weight[5] = -1.0 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 1.0 * x * y + 0.5 * x * z - 1.0 * y * z + 0.5 * pow ( z, 2 ) - 0.5 * z;
                        weight[6] = 0.5 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) * z + 0.25 * pow ( x, 2 ) - 0.5 * x * y + 0.5 * x * z
                                - 0.25 * x + 0.5 * y * pow ( z, 2 ) - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 );
                        weight[7] = -1.0 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.5 * pow ( x, 2 ) - 1.0 * y * z + 1.0 * y + 0.5 * pow ( z, 2 ) - 1.0 * z + 0.5;
                        weight[8] = 0.5 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) * z + 0.25 * pow ( x, 2 ) + 0.5 * x * y - 0.5 * x * z
                                + 0.25 * x + 0.5 * y * pow ( z, 2 ) - 0.5 * pow ( z, 3 ) + 0.25 * pow ( z, 2 );
                        weight[9] = 1.0 * z * ( x - 1 );
                        weight[10] = -1.0 * z * ( x + 1 );
                        weight[11] = 1.0 * z * ( -x + 1 );
                        weight[12] = 1.0 * z * ( x + 1 );
                        weight[13] = 0;
                        break;
                    case 2:
                        weight[0] = ( x + z ) * ( 0.5 * y - 0.25 ) * ( x + z - 1.0 );
                        weight[1] = -0.25 * ( x + 1.0 ) * ( 4 * y - 2.0 ) * ( x + z - 1 );
                        weight[2] = 0.5 * pow ( x, 2 ) * y - 0.25 * pow ( x, 2 ) + 1.0 * x * y * z + 0.5 * x * y - 0.25 * x + 0.5 * y * pow ( z, 2 )
                                + 0.5 * y * z + 0.25 * pow ( z, 2 ) + 0.25 * z;
                        weight[3] = 1.0 * x * y * ( -x - z + 1 );
                        weight[4] = y * ( 2.0 * pow ( x, 2 ) - 2.0 * pow ( z, 2 ) + 4.0 * z - 2.0 );
                        weight[5] = -y * ( 1.0 * pow ( x, 2 ) + 1.0 * x * z + 1.0 * x + 2.0 * z );
                        weight[6] = 0.125 * ( x + z ) * ( 4 * y + 2.0 ) * ( x + z - 1.0 );
                        weight[7] = -1.0 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) - 1.0 * x * y * z - 0.5 * x * z - 1.0 * y * z + 1.0 * y - 0.5 * z + 0.5;
                        weight[8] = 0.5 * pow ( x, 2 ) * y + 0.25 * pow ( x, 2 ) + 1.0 * x * y * z + 0.5 * x * y + 0.25 * x + 0.5 * y * pow ( z, 2 )
                                + 0.5 * y * z - 0.25 * pow ( z, 2 ) - 0.25 * z;
                        weight[9] = 1.0 * z * ( x + z - 1.0 );
                        weight[10] = -1.0 * z * ( x + z + 1.0 );
                        weight[11] = 1.0 * z * ( -x - z + 1 );
                        weight[12] = 1.0 * z * ( x + z + 1.0 );
                        weight[13] = 0;
                        break;
                    case 3:
                        weight[0] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) - 0.5 * x * y - 0.5 * x * z + 0.25 * x
                                + 0.5 * y * pow ( z, 2 ) + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 );
                        weight[1] = -1.0 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) * z + 0.5 * pow ( x, 2 ) - 1.0 * y * z + 1.0 * y - 0.5 * pow ( z, 2 ) + 1.0 * z - 0.5;
                        weight[2] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) + 0.5 * x * y + 0.5 * x * z - 0.25 * x
                                + 0.5 * y * pow ( z, 2 ) + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 );
                        weight[3] = -1.0 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) * z + 1.0 * x * y + 0.5 * x * z - 1.0 * y * z - 0.5 * pow ( z, 2 ) + 0.5 * z;
                        weight[4] = 2.0 * pow ( x, 2 ) * y - 2.0 * y * pow ( z, 2 ) + 4.0 * y * z - 2.0 * y - 2.0 * pow ( z, 3 ) + 3.0 * pow ( z, 2 ) - 1.0 * z;
                        weight[5] = -1.0 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) * z - 1.0 * x * y - 0.5 * x * z - 1.0 * y * z - 0.5 * pow ( z, 2 ) + 0.5 * z;
                        weight[6] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z + 0.25 * pow ( x, 2 ) - 0.5 * x * y - 0.25 * x + 0.5 * y * pow ( z, 2 )
                                + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 ) - 0.5 * z;
                        weight[7] = -( 0.5 * ( x - 1.0 ) * ( 0.5 * x - 0.5 * z + 0.5 ) + 0.25 * ( x + 1.0 ) * ( x + z - 1 ) ) * ( 2 * y + z + 1.0 );
                        weight[8] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z + 0.25 * pow ( x, 2 ) + 0.5 * x * y + 0.25 * x + 0.5 * y * pow ( z, 2 )
                                + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 ) - 0.5 * z;
                        weight[9] = 1.0 * z * ( x - 1 );
                        weight[10] = -1.0 * z * ( x + 1 );
                        weight[11] = 1.0 * z * ( -x + 1 );
                        weight[12] = 1.0 * z * ( x + 1 );
                        weight[13] = 0;
                        break;
                    case 4:
                        weight[0] = 0.5 * pow ( x, 2 ) * y - 0.25 * pow ( x, 2 ) - 1.0 * x * y * z - 0.5 * x * y + 0.25 * x + 0.5 * y * pow ( z, 2 )
                                + 0.5 * y * z + 0.25 * pow ( z, 2 ) + 0.25 * z;
                        weight[1] = -0.5 * ( x - 1.0 ) * ( 4 * y - 2.0 ) * ( 0.5 * x - 0.5 * z + 0.5 );
                        weight[2] = ( 0.25 * x - 0.25 * z ) * ( 2.0 * y - 1.0 ) * ( x - z + 1.0 );
                        weight[3] = y * ( -1.0 * pow ( x, 2 ) + 1.0 * x * z + 1.0 * x - 2.0 * z );
                        weight[4] = y * ( 2.0 * pow ( x, 2 ) - 2.0 * pow ( z, 2 ) + 4.0 * z - 2.0 );
                        weight[5] = 1.0 * x * y * ( -x + z - 1 );
                        weight[6] = 0.5 * pow ( x, 2 ) * y + 0.25 * pow ( x, 2 ) - 1.0 * x * y * z - 0.5 * x * y - 0.25 * x + 0.5 * y * pow ( z, 2 )
                                + 0.5 * y * z - 0.25 * pow ( z, 2 ) - 0.25 * z;
                        weight[7] = -1.0 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) + 1.0 * x * y * z + 0.5 * x * z - 1.0 * y * z + 1.0 * y - 0.5 * z + 0.5;
                        weight[8] = ( x - z + 1.0 ) * ( 0.5 * ( 0.25 * x - 0.25 * z ) * ( y - z ) + 0.5 * ( 0.25 * x - 0.25 * z ) * ( y - z + 1.0 )
                                + 0.125 * ( x - z ) * ( y + z ) + 0.125 * ( x - z ) * ( y + z + 1.0 ) );
                        weight[9] = 1.0 * z * ( x - z - 1 );
                        weight[10] = 1.0 * z * ( -x + z - 1 );
                        weight[11] = 1.0 * z * ( -x + z + 1.0 );
                        weight[12] = 1.0 * z * ( x - z + 1.0 );
                        weight[13] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                if ( xx == 0 && yy == 0 )
                {
                    weight[0] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) - 0.5 * x * y - 0.5 * x * z + 0.25 * x
                            + 0.5 * y * pow ( z, 2 ) + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 );
                    weight[1] = -0.5 * ( x - 1.0 ) * ( 4 * y - 2.0 ) * ( 0.5 * x - 0.5 * z + 0.5 );
                    weight[7] = -1.0 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) + 1.0 * x * y * z + 0.5 * x * z - 1.0 * y * z + 1.0 * y - 0.5 * z + 0.5;
                    weight[2] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) + 0.5 * x * y + 0.5 * x * z - 0.25 * x
                            + 0.5 * y * pow ( z, 2 ) + 0.5 * pow ( z, 3 ) - 0.25 * pow ( z, 2 );
                }

                //rescaling due to changing coordinates
                weight[0] = 2.0 * weight[0];
                weight[1] = 2.0 * weight[1];
                weight[2] = 2.0 * weight[2];
                weight[3] = 2.0 * weight[3];
                weight[4] = 2.0 * weight[4];
                weight[5] = 2.0 * weight[5];
                weight[6] = 2.0 * weight[6];
                weight[7] = 2.0 * weight[7];
                weight[8] = 2.0 * weight[8];
                weight[9] = 2.0 * weight[9];
                weight[10] = 2.0 * weight[10];
                weight[11] = 2.0 * weight[11];
                weight[12] = 2.0 * weight[12];
                weight[13] = 2.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_z ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {
                weight[0] = 0.0;
            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = -0.25 * x - 0.25;
                        weight[1] = 0.25 * x - 0.25;
                        weight[2] = 0.25 * x - 0.25;
                        weight[3] = -0.25 * x - 0.25;
                        weight[4] = 1;
                        break;
                    case 2:
                        weight[0] = 0.25 * y - 0.25;
                        weight[1] = -0.25 * y - 0.25;
                        weight[2] = -0.25 * y - 0.25;
                        weight[3] = 0.25 * y - 0.25;
                        weight[4] = 1;
                        break;
                    case 3:
                        weight[0] = 0.25 * x - 0.25;
                        weight[1] = -0.25 * x - 0.25;
                        weight[2] = -0.25 * x - 0.25;
                        weight[3] = 0.25 * x - 0.25;
                        weight[4] = 1;
                        break;
                    case 4:
                        weight[0] = -0.25 * y - 0.25;
                        weight[1] = 0.25 * y - 0.25;
                        weight[2] = 0.25 * y - 0.25;
                        weight[3] = -0.25 * y - 0.25;
                        weight[4] = 1;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );

                }

                //rescaling due to changing coordinates
                weight[0] = 1.0 * weight[0];
                weight[1] = 1.0 * weight[1];
                weight[2] = 1.0 * weight[2];
                weight[3] = 1.0 * weight[3];
                weight[4] = 1.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = -0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z + 0.25 * pow ( x, 2 ) + 0.5 * x * z + 0.25 * x + 0.5 * pow ( y, 2 ) * z
                                - 1.5 * y * pow ( z, 2 ) + 0.5 * y * z + 0.5 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[1] = 0.5 * pow ( x, 2 ) * y - 1.0 * pow ( x, 2 ) - 0.5 * pow ( y, 2 ) + 1.0 * y * z;
                        weight[2] = -0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z + 0.25 * pow ( x, 2 ) - 0.5 * x * z - 0.25 * x + 0.5 * pow ( y, 2 ) * z
                                - 1.5 * y * pow ( z, 2 ) + 0.5 * y * z + 0.5 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[3] = 0.5 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) - 0.5 * x * y + 0.5 * x - 0.5 * pow ( y, 2 ) + 1.0 * y * z - 0.5 * y;
                        weight[4] = -2.0 * pow ( x, 2 ) * z + 2.0 * pow ( x, 2 ) - 2.0 * pow ( y, 2 ) * z + 2.0 * pow ( y, 2 ) + 6.0 * y * pow ( z, 2 ) - 6.0 * y * z
                                + 1.0 * y - 4.0 * pow ( z, 3 ) + 3.0 * pow ( z, 2 ) + 4.0 * z - 3.0;
                        weight[5] = 0.5 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) + 0.5 * x * y - 0.5 * x - 0.5 * pow ( y, 2 ) + 1.0 * y * z - 0.5 * y;
                        weight[6] = -0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) + 0.5 * x * y - 0.5 * x * z + 0.25 * x
                                + 0.5 * pow ( y, 2 ) * z - 1.5 * y * pow ( z, 2 ) + 0.5 * y * z + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[7] = y * ( 0.5 * pow ( x, 2 ) - 0.5 * y + 1.0 * z - 1.0 );
                        weight[8] = -0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) - 0.5 * x * y + 0.5 * x * z - 0.25 * x
                                + 0.5 * pow ( y, 2 ) * z - 1.5 * y * pow ( z, 2 ) + 0.5 * y * z + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[9] = 1.0 * x * y - 2.0 * x * z - 1.0 * x - 1.0 * y - 2.0 * z + 1.0;
                        weight[10] = -1.0 * x * y + 2.0 * x * z + 1.0 * x - 1.0 * y - 2.0 * z + 1.0;
                        weight[11] = -1.0 * x * y + 2.0 * x * z - 1.0 * x + 1.0 * y - 2.0 * z + 1.0;
                        weight[12] = 1.0 * x * y - 2.0 * x * z + 1.0 * x + 1.0 * y - 2.0 * z + 1.0;
                        weight[13] = 4.0 * z - 1.0;
                        break;
                    case 2:
                        weight[0] = 0.5 * pow ( x, 2 ) * z + 0.5 * x * pow ( y, 2 ) - 0.5 * x * y + 1.5 * x * pow ( z, 2 ) - 0.5 * x * z + 0.5 * pow ( y, 2 ) * z
                                - 0.25 * pow ( y, 2 ) - 0.5 * y * z + 0.25 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[1] = -0.5 * pow ( x, 2 ) - 0.5 * x * pow ( y, 2 ) + 0.5 * x * y - 1.0 * x * z + 0.5 * x - 0.5 * pow ( y, 2 ) + 0.5 * y;
                        weight[2] = 0.5 * pow ( x, 2 ) * z + 0.5 * x * pow ( y, 2 ) + 1.5 * x * pow ( z, 2 ) - 0.5 * x * z - 0.5 * x + 0.5 * pow ( y, 2 ) * z
                                + 0.25 * pow ( y, 2 ) + 0.5 * y * z + 0.25 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[3] = x * ( -0.5 * x - 0.5 * pow ( y, 2 ) - 1.0 * z + 1.0 );
                        weight[4] = -2.0 * pow ( x, 2 ) * z + 2.0 * pow ( x, 2 ) - 6.0 * x * pow ( z, 2 ) + 6.0 * x * z - 1.0 * x - 2.0 * pow ( y, 2 ) * z
                                + 2.0 * pow ( y, 2 ) - 4.0 * pow ( z, 3 ) + 3.0 * pow ( z, 2 ) + 4.0 * z - 3.0;
                        weight[5] = -0.5 * pow ( x, 2 ) - 0.5 * x * pow ( y, 2 ) - 1.0 * x * z - 1.0 * pow ( y, 2 );
                        weight[6] = 0.5 * pow ( x, 2 ) * z + 0.5 * x * pow ( y, 2 ) + 0.5 * x * y + 1.5 * x * pow ( z, 2 ) - 0.5 * x * z + 0.5 * pow ( y, 2 ) * z
                                - 0.25 * pow ( y, 2 ) + 0.5 * y * z - 0.25 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[7] = -0.5 * pow ( x, 2 ) - 0.5 * x * pow ( y, 2 ) - 0.5 * x * y - 1.0 * x * z + 0.5 * x - 0.5 * pow ( y, 2 ) - 0.5 * y;
                        weight[8] = 0.5 * pow ( x, 2 ) * z + 0.5 * x * pow ( y, 2 ) + 1.5 * x * pow ( z, 2 ) - 0.5 * x * z - 0.5 * x + 0.5 * pow ( y, 2 ) * z
                                + 0.25 * pow ( y, 2 ) - 0.5 * y * z - 0.25 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[9] = 1.0 * x * y - 1.0 * x + 2.0 * y * z - 1.0 * y - 2.0 * z + 1.0;
                        weight[10] = -1.0 * x * y + 1.0 * x - 2.0 * y * z - 1.0 * y - 2.0 * z + 1.0;
                        weight[11] = -1.0 * x * y - 1.0 * x - 2.0 * y * z + 1.0 * y - 2.0 * z + 1.0;
                        weight[12] = 1.0 * x * y + 1.0 * x + 2.0 * y * z + 1.0 * y - 2.0 * z + 1.0;
                        weight[13] = 4.0 * z - 1.0;
                        break;
                    case 3:
                        weight[0] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) - 0.5 * x * y - 0.5 * x * z + 0.25 * x
                                + 0.5 * pow ( y, 2 ) * z + 1.5 * y * pow ( z, 2 ) - 0.5 * y * z + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[1] = y * ( -0.5 * pow ( x, 2 ) - 0.5 * y - 1.0 * z + 1.0 );
                        weight[2] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z - 0.25 * pow ( x, 2 ) + 0.5 * x * y + 0.5 * x * z - 0.25 * x
                                + 0.5 * pow ( y, 2 ) * z + 1.5 * y * pow ( z, 2 ) - 0.5 * y * z + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[3] = -0.5 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) + 0.5 * x * y + 0.5 * x - 0.5 * pow ( y, 2 ) - 1.0 * y * z + 0.5 * y;
                        weight[4] = -2.0 * pow ( x, 2 ) * z + 2.0 * pow ( x, 2 ) - 2.0 * pow ( y, 2 ) * z + 2.0 * pow ( y, 2 ) - 6.0 * y * pow ( z, 2 )
                                + 6.0 * y * z - 1.0 * y - 4.0 * pow ( z, 3 ) + 3.0 * pow ( z, 2 ) + 4.0 * z - 3.0;
                        weight[5] = -0.5 * pow ( x, 2 ) * y - 0.5 * pow ( x, 2 ) - 0.5 * x * y - 0.5 * x - 0.5 * pow ( y, 2 ) - 1.0 * y * z + 0.5 * y;
                        weight[6] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z + 0.25 * pow ( x, 2 ) + 0.5 * x * z + 0.25 * x + 0.5 * pow ( y, 2 ) * z
                                + 1.5 * y * pow ( z, 2 ) - 0.5 * y * z - 0.5 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[7] = -0.5 * pow ( x, 2 ) * y - 1.0 * pow ( x, 2 ) - 0.5 * pow ( y, 2 ) - 1.0 * y * z;
                        weight[8] = 0.5 * pow ( x, 2 ) * y + 0.5 * pow ( x, 2 ) * z + 0.25 * pow ( x, 2 ) - 0.5 * x * z - 0.25 * x + 0.5 * pow ( y, 2 ) * z
                                + 1.5 * y * pow ( z, 2 ) - 0.5 * y * z - 0.5 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[9] = 1.0 * x * y + 2.0 * x * z - 1.0 * x - 1.0 * y - 2.0 * z + 1.0;
                        weight[10] = -1.0 * x * y - 2.0 * x * z + 1.0 * x - 1.0 * y - 2.0 * z + 1.0;
                        weight[11] = -1.0 * x * y - 2.0 * x * z - 1.0 * x + 1.0 * y - 2.0 * z + 1.0;
                        weight[12] = 1.0 * x * y + 2.0 * x * z + 1.0 * x + 1.0 * y - 2.0 * z + 1.0;
                        weight[13] = 4.0 * z - 1.0;
                        break;
                    case 4:
                        weight[0] = 0.5 * pow ( x, 2 ) * z - 0.5 * x * pow ( y, 2 ) - 1.5 * x * pow ( z, 2 ) + 0.5 * x * z + 0.5 * x + 0.5 * pow ( y, 2 ) * z
                                + 0.25 * pow ( y, 2 ) + 0.5 * y * z + 0.25 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[1] = -0.5 * pow ( x, 2 ) + 0.5 * x * pow ( y, 2 ) - 0.5 * x * y + 1.0 * x * z - 0.5 * x - 0.5 * pow ( y, 2 ) + 0.5 * y;
                        weight[2] = 0.5 * pow ( x, 2 ) * z - 0.5 * x * pow ( y, 2 ) + 0.5 * x * y - 1.5 * x * pow ( z, 2 ) + 0.5 * x * z + 0.5 * pow ( y, 2 ) * z
                                - 0.25 * pow ( y, 2 ) - 0.5 * y * z + 0.25 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[3] = -0.5 * pow ( x, 2 ) + 0.5 * x * pow ( y, 2 ) + 1.0 * x * z - 1.0 * pow ( y, 2 );
                        weight[4] = -2.0 * pow ( x, 2 ) * z + 2.0 * pow ( x, 2 ) + 6.0 * x * pow ( z, 2 ) - 6.0 * x * z + 1.0 * x - 2.0 * pow ( y, 2 ) * z
                                + 2.0 * pow ( y, 2 ) - 4.0 * pow ( z, 3 ) + 3.0 * pow ( z, 2 ) + 4.0 * z - 3.0;
                        weight[5] = x * ( -0.5 * x + 0.5 * pow ( y, 2 ) + 1.0 * z - 1.0 );
                        weight[6] = 0.5 * pow ( x, 2 ) * z - 0.5 * x * pow ( y, 2 ) - 1.5 * x * pow ( z, 2 ) + 0.5 * x * z + 0.5 * x + 0.5 * pow ( y, 2 ) * z
                                + 0.25 * pow ( y, 2 ) - 0.5 * y * z - 0.25 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[7] = -0.5 * pow ( x, 2 ) + 0.5 * x * pow ( y, 2 ) + 0.5 * x * y + 1.0 * x * z - 0.5 * x - 0.5 * pow ( y, 2 ) - 0.5 * y;
                        weight[8] = 0.5 * pow ( x, 2 ) * z - 0.5 * x * pow ( y, 2 ) - 0.5 * x * y - 1.5 * x * pow ( z, 2 ) + 0.5 * x * z + 0.5 * pow ( y, 2 ) * z
                                - 0.25 * pow ( y, 2 ) + 0.5 * y * z - 0.25 * y + 1.0 * pow ( z, 3 ) - 0.75 * pow ( z, 2 );
                        weight[9] = 1.0 * x * y - 1.0 * x - 2.0 * y * z - 1.0 * y - 2.0 * z + 1.0;
                        weight[10] = -1.0 * x * y + 1.0 * x + 2.0 * y * z - 1.0 * y - 2.0 * z + 1.0;
                        weight[11] = -1.0 * x * y - 1.0 * x + 2.0 * y * z + 1.0 * y - 2.0 * z + 1.0;
                        weight[12] = 1.0 * x * y + 1.0 * x - 2.0 * y * z + 1.0 * y - 2.0 * z + 1.0;
                        weight[13] = 4.0 * z - 1.0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 1.0 * weight[0];
                weight[1] = 1.0 * weight[1];
                weight[2] = 1.0 * weight[2];
                weight[3] = 1.0 * weight[3];
                weight[4] = 1.0 * weight[4];
                weight[5] = 1.0 * weight[5];
                weight[6] = 1.0 * weight[6];
                weight[7] = 1.0 * weight[7];
                weight[8] = 1.0 * weight[8];
                weight[9] = 1.0 * weight[9];
                weight[10] = 1.0 * weight[10];
                weight[11] = 1.0 * weight[11];
                weight[12] = 1.0 * weight[12];
                weight[13] = 1.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_xx ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {
                weight[0] = 0.0;
            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 2:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 3:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 4:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 4.0 * weight[0];
                weight[1] = 4.0 * weight[1];
                weight[2] = 4.0 * weight[2];
                weight[3] = 4.0 * weight[3];
                weight[4] = 4.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = -0.5 * ( y - z ) * ( -y + z + 1.0 );
                        weight[1] = 1.0 * y * ( -y + z + 1.0 ) - 2.0 * z;
                        weight[2] = -0.5 * ( y - z ) * ( -y + z + 1.0 );
                        weight[3] = -1.0 * pow ( y, 2 ) + 1.0 * y * z - 1.0 * z + 1.0;
                        weight[4] = 2.0 * ( y - z + 1.0 ) * ( y + z - 1.0 );
                        weight[5] = -1.0 * pow ( y, 2 ) + 1.0 * y * z - 1.0 * z + 1.0;
                        weight[6] = ( 0.5 * y - 0.5 * z ) * ( y - z + 1.0 );
                        weight[7] = 1.0 * y * ( -y + z - 1 );
                        weight[8] = ( 0.5 * y - 0.5 * z ) * ( y - z + 1.0 );
                        weight[9] = 0;
                        weight[10] = 0;
                        weight[11] = 0;
                        weight[12] = 0;
                        weight[13] = 0;
                        break;
                    case 2:
                        weight[0] = 0.5 * pow ( y, 2 ) - 0.5 * y + 0.5 * pow ( z, 2 );
                        weight[1] = -1.0 * pow ( y, 2 ) + 1.0 * y - 1.0 * z;
                        weight[2] = 0.5 * pow ( y, 2 ) - 0.5 * y + 0.5 * pow ( z, 2 );
                        weight[3] = -1.0 * pow ( y, 2 ) - 1.0 * z + 1.0;
                        weight[4] = 2.0 * ( y - z + 1.0 ) * ( y + z - 1.0 );
                        weight[5] = -1.0 * pow ( y, 2 ) - 1.0 * z + 1.0;
                        weight[6] = 0.5 * pow ( y, 2 ) + 0.5 * y + 0.5 * pow ( z, 2 );
                        weight[7] = -1.0 * pow ( y, 2 ) - 1.0 * y - 1.0 * z;
                        weight[8] = 0.5 * pow ( y, 2 ) + 0.5 * y + 0.5 * pow ( z, 2 );
                        weight[9] = 0;
                        weight[10] = 0;
                        weight[11] = 0;
                        weight[12] = 0;
                        weight[13] = 0;
                        break;
                    case 3:
                        weight[0] = 0.5 * ( y + z ) * ( y + z - 1.0 );
                        weight[1] = 1.0 * y * ( -y - z + 1 );
                        weight[2] = 0.5 * ( y + z ) * ( y + z - 1.0 );
                        weight[3] = -1.0 * pow ( y, 2 ) - 1.0 * y * z - 1.0 * z + 1.0;
                        weight[4] = 2.0 * ( y - z + 1.0 ) * ( y + z - 1.0 );
                        weight[5] = -1.0 * pow ( y, 2 ) - 1.0 * y * z - 1.0 * z + 1.0;
                        weight[6] = 0.5 * ( y + z ) * ( y + z + 1.0 );
                        weight[7] = -1.0 * y * ( y + z + 1.0 ) - 2.0 * z;
                        weight[8] = 0.5 * ( y + z ) * ( y + z + 1.0 );
                        weight[9] = 0;
                        weight[10] = 0;
                        weight[11] = 0;
                        weight[12] = 0;
                        weight[13] = 0;
                        break;
                    case 4:
                        weight[0] = 0.5 * pow ( y, 2 ) - 0.5 * y + 0.5 * pow ( z, 2 );
                        weight[1] = -1.0 * pow ( y, 2 ) + 1.0 * y - 1.0 * z;
                        weight[2] = 0.5 * pow ( y, 2 ) - 0.5 * y + 0.5 * pow ( z, 2 );
                        weight[3] = -1.0 * pow ( y, 2 ) - 1.0 * z + 1.0;
                        weight[4] = 2.0 * ( y - z + 1.0 ) * ( y + z - 1.0 );
                        weight[5] = -1.0 * pow ( y, 2 ) - 1.0 * z + 1.0;
                        weight[6] = 0.5 * pow ( y, 2 ) + 0.5 * y + 0.5 * pow ( z, 2 );
                        weight[7] = -1.0 * pow ( y, 2 ) - 1.0 * y - 1.0 * z;
                        weight[8] = 0.5 * pow ( y, 2 ) + 0.5 * y + 0.5 * pow ( z, 2 );
                        weight[9] = 0;
                        weight[10] = 0;
                        weight[11] = 0;
                        weight[12] = 0;
                        weight[13] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 4.0 * weight[0];
                weight[1] = 4.0 * weight[1];
                weight[2] = 4.0 * weight[2];
                weight[3] = 4.0 * weight[3];
                weight[4] = 4.0 * weight[4];
                weight[5] = 4.0 * weight[5];
                weight[6] = 4.0 * weight[6];
                weight[7] = 4.0 * weight[7];
                weight[8] = 4.0 * weight[8];
                weight[9] = 4.0 * weight[9];
                weight[10] = 4.0 * weight[10];
                weight[11] = 4.0 * weight[11];
                weight[12] = 4.0 * weight[12];
                weight[13] = 4.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_xy ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {

                weight[0] = 0.0;

            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0.250000000000000;
                        weight[1] = -0.250000000000000;
                        weight[2] = -0.250000000000000;
                        weight[3] = 0.250000000000000;
                        weight[4] = 0;
                        break;
                    case 2:
                        weight[0] = 0.250000000000000;
                        weight[1] = -0.250000000000000;
                        weight[2] = -0.250000000000000;
                        weight[3] = 0.250000000000000;
                        weight[4] = 0;
                        break;
                    case 3:
                        weight[0] = 0.250000000000000;
                        weight[1] = -0.250000000000000;
                        weight[2] = -0.250000000000000;
                        weight[3] = 0.250000000000000;
                        weight[4] = 0;
                        break;
                    case 4:
                        weight[0] = 0.250000000000000;
                        weight[1] = -0.250000000000000;
                        weight[2] = -0.250000000000000;
                        weight[3] = 0.250000000000000;
                        weight[4] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 4.0 * weight[0];
                weight[1] = 4.0 * weight[1];
                weight[2] = 4.0 * weight[2];
                weight[3] = 4.0 * weight[3];
                weight[4] = 4.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 1.0 * x * y - 1.0 * x * z - 0.5 * x - 0.5 * y + 0.25;
                        weight[1] = x * ( -2.0 * y + 1.0 * z + 1.0 );
                        weight[2] = 1.0 * x * y - 1.0 * x * z - 0.5 * x + 0.5 * y - 0.25;
                        weight[3] = -2.0 * x * y + 1.0 * x * z + 1.0 * y - 0.5 * z;
                        weight[4] = 4.0 * x * y;
                        weight[5] = -2.0 * x * y + 1.0 * x * z - 1.0 * y + 0.5 * z;
                        weight[6] = 1.0 * x * y - 1.0 * x * z + 0.5 * x - 0.5 * y + 0.5 * z - 0.25;
                        weight[7] = x * ( -2.0 * y + 1.0 * z - 1.0 );
                        weight[8] = 1.0 * x * y - 1.0 * x * z + 0.5 * x + 0.5 * y - 0.5 * z + 0.25;
                        weight[9] = 1.0 * z;
                        weight[10] = -1.0 * z;
                        weight[11] = -1.0 * z;
                        weight[12] = 1.0 * z;
                        weight[13] = 0;
                        break;
                    case 2:
                        weight[0] = 1.0 * x * y - 0.5 * x + 1.0 * y * z - 0.5 * y - 0.5 * z + 0.25;
                        weight[1] = -2.0 * x * y + 1.0 * x - 1.0 * y * z + 0.5 * z;
                        weight[2] = 1.0 * x * y - 0.5 * x + 1.0 * y * z + 0.5 * y - 0.25;
                        weight[3] = y * ( -2.0 * x - 1.0 * z + 1.0 );
                        weight[4] = 4.0 * x * y;
                        weight[5] = -y * ( 2.0 * x + 1.0 * z + 1.0 );
                        weight[6] = 1.0 * x * y + 0.5 * x + 1.0 * y * z - 0.5 * y + 0.5 * z - 0.25;
                        weight[7] = -2.0 * x * y - 1.0 * x - 1.0 * y * z - 0.5 * z;
                        weight[8] = 1.0 * x * y + 0.5 * x + 1.0 * y * z + 0.5 * y + 0.25;
                        weight[9] = 1.0 * z;
                        weight[10] = -1.0 * z;
                        weight[11] = -1.0 * z;
                        weight[12] = 1.0 * z;
                        weight[13] = 0;
                        break;
                    case 3:
                        weight[0] = 1.0 * x * y + 1.0 * x * z - 0.5 * x - 0.5 * y - 0.5 * z + 0.25;
                        weight[1] = x * ( -2.0 * y - 1.0 * z + 1.0 );
                        weight[2] = 1.0 * x * y + 1.0 * x * z - 0.5 * x + 0.5 * y + 0.5 * z - 0.25;
                        weight[3] = -2.0 * x * y - 1.0 * x * z + 1.0 * y + 0.5 * z;
                        weight[4] = 4.0 * x * y;
                        weight[5] = -2.0 * x * y - 1.0 * x * z - 1.0 * y - 0.5 * z;
                        weight[6] = 1.0 * x * y + 1.0 * x * z + 0.5 * x - 0.5 * y - 0.25;
                        weight[7] = -x * ( 2.0 * y + 1.0 * z + 1.0 );
                        weight[8] = 1.0 * x * y + 1.0 * x * z + 0.5 * x + 0.5 * y + 0.25;
                        weight[9] = 1.0 * z;
                        weight[10] = -1.0 * z;
                        weight[11] = -1.0 * z;
                        weight[12] = 1.0 * z;
                        weight[13] = 0;
                        break;
                    case 4:
                        weight[0] = 1.0 * x * y - 0.5 * x - 1.0 * y * z - 0.5 * y + 0.25;
                        weight[1] = -2.0 * x * y + 1.0 * x + 1.0 * y * z - 0.5 * z;
                        weight[2] = 1.0 * x * y - 0.5 * x - 1.0 * y * z + 0.5 * y + 0.5 * z - 0.25;
                        weight[3] = y * ( -2.0 * x + 1.0 * z + 1.0 );
                        weight[4] = 4.0 * x * y;
                        weight[5] = y * ( -2.0 * x + 1.0 * z - 1.0 );
                        weight[6] = 1.0 * x * y + 0.5 * x - 1.0 * y * z - 0.5 * y - 0.25;
                        weight[7] = -2.0 * x * y - 1.0 * x + 1.0 * y * z + 0.5 * z;
                        weight[8] = 1.0 * x * y + 0.5 * x - 1.0 * y * z + 0.5 * y - 0.5 * z + 0.25;
                        weight[9] = 1.0 * z;
                        weight[10] = -1.0 * z;
                        weight[11] = -1.0 * z;
                        weight[12] = 1.0 * z;
                        weight[13] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                if ( xx == 0 && yy == 0 )
                {
                    weight[2] = 1.0 * x * y - 0.5 * x - 1.0 * y * z + 0.5 * y + 0.5 * z - 0.25;
                    weight[3] = -2.0 * x * y - 1.0 * x * z + 1.0 * y + 0.5 * z;
                    weight[5] = -2.0 * x * y - 1.0 * x * z - 1.0 * y - 0.5 * z;
                    weight[8] = 1.0 * x * y + 0.5 * x - 1.0 * y * z + 0.5 * y - 0.5 * z + 0.25;
                }

                //rescaling due to changing coordinates
                weight[0] = 4.0 * weight[0];
                weight[1] = 4.0 * weight[1];
                weight[2] = 4.0 * weight[2];
                weight[3] = 4.0 * weight[3];
                weight[4] = 4.0 * weight[4];
                weight[5] = 4.0 * weight[5];
                weight[6] = 4.0 * weight[6];
                weight[7] = 4.0 * weight[7];
                weight[8] = 4.0 * weight[8];
                weight[9] = 4.0 * weight[9];
                weight[10] = 4.0 * weight[10];
                weight[11] = 4.0 * weight[11];
                weight[12] = 4.0 * weight[12];
                weight[13] = 4.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_xz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {

                weight[0] = 0.0;

            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = -0.250000000000000;
                        weight[1] = 0.250000000000000;
                        weight[2] = 0.250000000000000;
                        weight[3] = -0.250000000000000;
                        weight[4] = 0;
                        break;
                    case 2:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 3:
                        weight[0] = 0.250000000000000;
                        weight[1] = -0.250000000000000;
                        weight[2] = -0.250000000000000;
                        weight[3] = 0.250000000000000;
                        weight[4] = 0;
                        break;
                    case 4:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 2.0 * weight[0];
                weight[1] = 2.0 * weight[1];
                weight[2] = 2.0 * weight[2];
                weight[3] = 2.0 * weight[3];
                weight[4] = 2.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = -1.0 * x * y + 1.0 * x * z + 0.5 * x + 0.5 * z + 0.25;
                        weight[1] = x * ( 1.0 * y - 2.0 );
                        weight[2] = -1.0 * x * y + 1.0 * x * z + 0.5 * x - 0.5 * z - 0.25;
                        weight[3] = ( 1.0 * x - 0.5 ) * ( y - 1.0 );
                        weight[4] = x * ( -4.0 * z + 4.0 );
                        weight[5] = ( 1.0 * x + 0.5 ) * ( y - 1.0 );
                        weight[6] = -1.0 * x * y + 1.0 * x * z - 0.5 * x + 0.5 * y - 0.5 * z + 0.25;
                        weight[7] = x * y;
                        weight[8] = -1.0 * x * y + 1.0 * x * z - 0.5 * x - 0.5 * y + 0.5 * z - 0.25;
                        weight[9] = 1.0 * y - 2.0 * z - 1.0;
                        weight[10] = -1.0 * y + 2.0 * z + 1.0;
                        weight[11] = -1.0 * y + 2.0 * z - 1.0;
                        weight[12] = 1.0 * y - 2.0 * z + 1.0;
                        weight[13] = 0;
                        break;
                    case 2:
                        weight[0] = 1.0 * x * z + 0.5 * pow ( y, 2 ) - 0.5 * y + 1.5 * pow ( z, 2 ) - 0.5 * z;
                        weight[1] = -1.0 * x - 0.5 * pow ( y, 2 ) + 0.5 * y - 1.0 * z + 0.5;
                        weight[2] = 1.0 * x * z + 0.5 * pow ( y, 2 ) + 1.5 * pow ( z, 2 ) - 0.5 * z - 0.5;
                        weight[3] = -1.0 * x - 0.5 * pow ( y, 2 ) - 1.0 * z + 1.0;
                        weight[4] = -4.0 * x * z + 4.0 * x - 6.0 * pow ( z, 2 ) + 6.0 * z - 1.0;
                        weight[5] = -1.0 * x - 0.5 * pow ( y, 2 ) - 1.0 * z;
                        weight[6] = 1.0 * x * z + 0.5 * pow ( y, 2 ) + 0.5 * y + 1.5 * pow ( z, 2 ) - 0.5 * z;
                        weight[7] = -1.0 * x - 0.5 * pow ( y, 2 ) - 0.5 * y - 1.0 * z + 0.5;
                        weight[8] = 1.0 * x * z + 0.5 * pow ( y, 2 ) + 1.5 * pow ( z, 2 ) - 0.5 * z - 0.5;
                        weight[9] = 1.0 * y - 1.0;
                        weight[10] = -1.0 * y + 1.0;
                        weight[11] = -1.0 * y - 1.0;
                        weight[12] = 1.0 * y + 1.0;
                        weight[13] = 0;
                        break;
                    case 3:
                        weight[0] = 1.0 * x * y + 1.0 * x * z - 0.5 * x - 0.5 * y - 0.5 * z + 0.25;
                        weight[1] = -1.0 * x * y;
                        weight[2] = 1.0 * x * y + 1.0 * x * z - 0.5 * x + 0.5 * y + 0.5 * z - 0.25;
                        weight[3] = ( -1.0 * x + 0.5 ) * ( y + 1.0 );
                        weight[4] = x * ( -4.0 * z + 4.0 );
                        weight[5] = -( 1.0 * x + 0.5 ) * ( y + 1.0 );
                        weight[6] = 1.0 * x * y + 1.0 * x * z + 0.5 * x + 0.5 * z + 0.25;
                        weight[7] = -x * ( 1.0 * y + 2.0 );
                        weight[8] = 1.0 * x * y + 1.0 * x * z + 0.5 * x - 0.5 * z - 0.25;
                        weight[9] = 1.0 * y + 2.0 * z - 1.0;
                        weight[10] = -1.0 * y - 2.0 * z + 1.0;
                        weight[11] = -1.0 * y - 2.0 * z - 1.0;
                        weight[12] = 1.0 * y + 2.0 * z + 1.0;
                        weight[13] = 0;
                        break;
                    case 4:
                        weight[0] = 1.0 * x * z - 0.5 * pow ( y, 2 ) - 1.5 * pow ( z, 2 ) + 0.5 * z + 0.5;
                        weight[1] = -1.0 * x + 0.5 * pow ( y, 2 ) - 0.5 * y + 1.0 * z - 0.5;
                        weight[2] = 1.0 * x * z - 0.5 * pow ( y, 2 ) + 0.5 * y - 1.5 * pow ( z, 2 ) + 0.5 * z;
                        weight[3] = -1.0 * x + 0.5 * pow ( y, 2 ) + 1.0 * z;
                        weight[4] = -4.0 * x * z + 4.0 * x + 6.0 * pow ( z, 2 ) - 6.0 * z + 1.0;
                        weight[5] = -1.0 * x + 0.5 * pow ( y, 2 ) + 1.0 * z - 1.0;
                        weight[6] = 1.0 * x * z - 0.5 * pow ( y, 2 ) - 1.5 * pow ( z, 2 ) + 0.5 * z + 0.5;
                        weight[7] = -1.0 * x + 0.5 * pow ( y, 2 ) + 0.5 * y + 1.0 * z - 0.5;
                        weight[8] = 1.0 * x * z - 0.5 * pow ( y, 2 ) - 0.5 * y - 1.5 * pow ( z, 2 ) + 0.5 * z;
                        weight[9] = 1.0 * y - 1.0;
                        weight[10] = -1.0 * y + 1.0;
                        weight[11] = -1.0 * y - 1.0;
                        weight[12] = 1.0 * y + 1.0;
                        weight[13] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                if ( xx == 0 && yy == 0 )
                {
                    weight[2] = 1.0 * x * z - 0.5 * pow ( y, 2 ) + 0.5 * y - 1.5 * pow ( z, 2 ) + 0.5 * z;
                    weight[3] = -1.0 * x + 0.5 * pow ( y, 2 ) + 1.0 * z;
                    weight[8] = 1.0 * x * z - 0.5 * pow ( y, 2 ) - 0.5 * y - 1.5 * pow ( z, 2 ) + 0.5 * z;
                }

                //rescaling due to changing coordinates
                weight[0] = 2.0 * weight[0];
                weight[1] = 2.0 * weight[1];
                weight[2] = 2.0 * weight[2];
                weight[3] = 2.0 * weight[3];
                weight[4] = 2.0 * weight[4];
                weight[5] = 2.0 * weight[5];
                weight[6] = 2.0 * weight[6];
                weight[7] = 2.0 * weight[7];
                weight[8] = 2.0 * weight[8];
                weight[9] = 2.0 * weight[9];
                weight[10] = 2.0 * weight[10];
                weight[11] = 2.0 * weight[11];
                weight[12] = 2.0 * weight[12];
                weight[13] = 2.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_yy ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {

                weight[0] = 0.0;

            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 2:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 3:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 4:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 4.0 * weight[0];
                weight[1] = 4.0 * weight[1];
                weight[2] = 4.0 * weight[2];
                weight[3] = 4.0 * weight[3];
                weight[4] = 4.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0.5 * pow ( x, 2 ) - 0.5 * x + 0.5 * pow ( z, 2 );
                        weight[1] = -1.0 * pow ( x, 2 ) - 1.0 * z + 1.0;
                        weight[2] = 0.5 * pow ( x, 2 ) + 0.5 * x + 0.5 * pow ( z, 2 );
                        weight[3] = -1.0 * pow ( x, 2 ) + 1.0 * x - 1.0 * z;
                        weight[4] = 2.0 * ( x - z + 1.0 ) * ( x + z - 1.0 );
                        weight[5] = -1.0 * pow ( x, 2 ) - 1.0 * x - 1.0 * z;
                        weight[6] = 0.5 * pow ( x, 2 ) - 0.5 * x + 0.5 * pow ( z, 2 );
                        weight[7] = -1.0 * pow ( x, 2 ) - 1.0 * z + 1.0;
                        weight[8] = 0.5 * pow ( x, 2 ) + 0.5 * x + 0.5 * pow ( z, 2 );
                        weight[9] = 0;
                        weight[10] = 0;
                        weight[11] = 0;
                        weight[12] = 0;
                        weight[13] = 0;
                        break;
                    case 2:
                        weight[0] = 0.5 * ( x + z ) * ( x + z - 1.0 );
                        weight[1] = -1.0 * ( x + 1 ) * ( x + z - 1 );
                        weight[2] = 0.5 * ( x + z ) * ( x + z + 1.0 );
                        weight[3] = 1.0 * x * ( -x - z + 1 );
                        weight[4] = 2.0 * ( x - z + 1.0 ) * ( x + z - 1.0 );
                        weight[5] = -1.0 * x * ( x + z + 1.0 ) - 2.0 * z;
                        weight[6] = 0.5 * ( x + z ) * ( x + z - 1.0 );
                        weight[7] = -1.0 * pow ( x, 2 ) - 1.0 * x * z - 1.0 * z + 1.0;
                        weight[8] = 0.5 * ( x + z ) * ( x + z + 1.0 );
                        weight[9] = 0;
                        weight[10] = 0;
                        weight[11] = 0;
                        weight[12] = 0;
                        weight[13] = 0;
                        break;
                    case 3:
                        weight[0] = 0.5 * pow ( x, 2 ) - 0.5 * x + 0.5 * pow ( z, 2 );
                        weight[1] = -1.0 * pow ( x, 2 ) - 1.0 * z + 1.0;
                        weight[2] = 0.5 * pow ( x, 2 ) + 0.5 * x + 0.5 * pow ( z, 2 );
                        weight[3] = -1.0 * pow ( x, 2 ) + 1.0 * x - 1.0 * z;
                        weight[4] = 2.0 * ( x - z + 1.0 ) * ( x + z - 1.0 );
                        weight[5] = -1.0 * pow ( x, 2 ) - 1.0 * x - 1.0 * z;
                        weight[6] = 0.5 * pow ( x, 2 ) - 0.5 * x + 0.5 * pow ( z, 2 );
                        weight[7] = -1.0 * pow ( x, 2 ) - 1.0 * z + 1.0;
                        weight[8] = 0.5 * pow ( x, 2 ) + 0.5 * x + 0.5 * pow ( z, 2 );
                        weight[9] = 0;
                        weight[10] = 0;
                        weight[11] = 0;
                        weight[12] = 0;
                        weight[13] = 0;
                        break;
                    case 4:
                        weight[0] = ( -0.5 * x + 0.5 * z ) * ( -x + z + 1.0 );
                        weight[1] = -1.0 * pow ( x, 2 ) + 1.0 * x * z - 1.0 * z + 1.0;
                        weight[2] = 2.0 * ( 0.25 * x - 0.25 * z ) * ( x - z + 1.0 );
                        weight[3] = 1.0 * x * ( -x + z + 1.0 ) - 2.0 * z;
                        weight[4] = 2.0 * ( x - z + 1.0 ) * ( x + z - 1.0 );
                        weight[5] = 1.0 * x * ( -x + z - 1 );
                        weight[6] = ( -0.5 * x + 0.5 * z ) * ( -x + z + 1.0 );
                        weight[7] = -1.0 * pow ( x, 2 ) + 1.0 * x * z - 1.0 * z + 1.0;
                        weight[8] = ( 0.5 * x - 0.5 * z ) * ( x - z + 1.0 );
                        weight[9] = 0;
                        weight[10] = 0;
                        weight[11] = 0;
                        weight[12] = 0;
                        weight[13] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                if ( xx == 0 && yy == 0 )
                {
                    weight[2] = 2.0 * ( 0.25 * x - 0.25 * z ) * ( x - z + 1.0 );
                    weight[3] = 1.0 * x * ( -x + z + 1.0 ) - 2.0 * z;
                    weight[8] = ( 0.5 * x - 0.5 * z ) * ( x - z + 1.0 );
                }

                //rescaling due to changing coordinates
                weight[0] = 4.0 * weight[0];
                weight[1] = 4.0 * weight[1];
                weight[2] = 4.0 * weight[2];
                weight[3] = 4.0 * weight[3];
                weight[4] = 4.0 * weight[4];
                weight[5] = 4.0 * weight[5];
                weight[6] = 4.0 * weight[6];
                weight[7] = 4.0 * weight[7];
                weight[8] = 4.0 * weight[8];
                weight[9] = 4.0 * weight[9];
                weight[10] = 4.0 * weight[10];
                weight[11] = 4.0 * weight[11];
                weight[12] = 4.0 * weight[12];
                weight[13] = 4.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_yz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {

                weight[0] = 0.0;

            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 2:
                        weight[0] = 0.250000000000000;
                        weight[1] = -0.250000000000000;
                        weight[2] = -0.250000000000000;
                        weight[3] = 0.250000000000000;
                        weight[4] = 0;
                        break;
                    case 3:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 4:
                        weight[0] = -0.250000000000000;
                        weight[1] = 0.250000000000000;
                        weight[2] = 0.250000000000000;
                        weight[3] = -0.250000000000000;
                        weight[4] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 2.0 * weight[0];
                weight[1] = 2.0 * weight[1];
                weight[2] = 2.0 * weight[2];
                weight[3] = 2.0 * weight[3];
                weight[4] = 2.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx >= yy && xx <= -yy )
                {
                    indicator = 1;
                }
                else if ( xx > yy && xx > -yy )
                {
                    indicator = 2;
                }
                else if ( xx <= yy && xx >= -yy )
                {
                    indicator = 3;
                }
                else if ( xx < yy && xx < -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = -0.5 * pow ( x, 2 ) + 1.0 * y * z - 1.5 * pow ( z, 2 ) + 0.5 * z + 0.5;
                        weight[1] = 0.5 * pow ( x, 2 ) - 1.0 * y + 1.0 * z;
                        weight[2] = -0.5 * pow ( x, 2 ) + 1.0 * y * z - 1.5 * pow ( z, 2 ) + 0.5 * z + 0.5;
                        weight[3] = 0.5 * pow ( x, 2 ) - 0.5 * x - 1.0 * y + 1.0 * z - 0.5;
                        weight[4] = -4.0 * y * z + 4.0 * y + 6.0 * pow ( z, 2 ) - 6.0 * z + 1.0;
                        weight[5] = 0.5 * pow ( x, 2 ) + 0.5 * x - 1.0 * y + 1.0 * z - 0.5;
                        weight[6] = -0.5 * pow ( x, 2 ) + 0.5 * x + 1.0 * y * z - 1.5 * pow ( z, 2 ) + 0.5 * z;
                        weight[7] = 0.5 * pow ( x, 2 ) - 1.0 * y + 1.0 * z - 1.0;
                        weight[8] = -0.5 * pow ( x, 2 ) - 0.5 * x + 1.0 * y * z - 1.5 * pow ( z, 2 ) + 0.5 * z;
                        weight[9] = 1.0 * x - 1.0;
                        weight[10] = -1.0 * x - 1.0;
                        weight[11] = -1.0 * x + 1.0;
                        weight[12] = 1.0 * x + 1.0;
                        weight[13] = 0;
                        break;
                    case 2:
                        weight[0] = 1.0 * x * y - 0.5 * x + 1.0 * y * z - 0.5 * y - 0.5 * z + 0.25;
                        weight[1] = ( x + 1.0 ) * ( -1.0 * y + 0.5 );
                        weight[2] = 1.0 * x * y + 1.0 * y * z + 0.5 * y + 0.5 * z + 0.25;
                        weight[3] = -1.0 * x * y;
                        weight[4] = y * ( -4.0 * z + 4.0 );
                        weight[5] = -y * ( 1.0 * x + 2.0 );
                        weight[6] = 1.0 * x * y + 0.5 * x + 1.0 * y * z - 0.5 * y + 0.5 * z - 0.25;
                        weight[7] = -( x + 1.0 ) * ( 1.0 * y + 0.5 );
                        weight[8] = 1.0 * x * y + 1.0 * y * z + 0.5 * y - 0.5 * z - 0.25;
                        weight[9] = 1.0 * x + 2.0 * z - 1.0;
                        weight[10] = -1.0 * x - 2.0 * z - 1.0;
                        weight[11] = -1.0 * x - 2.0 * z + 1.0;
                        weight[12] = 1.0 * x + 2.0 * z + 1.0;
                        weight[13] = 0;
                        break;
                    case 3:
                        weight[0] = 0.5 * pow ( x, 2 ) - 0.5 * x + 1.0 * y * z + 1.5 * pow ( z, 2 ) - 0.5 * z;
                        weight[1] = -0.5 * pow ( x, 2 ) - 1.0 * y - 1.0 * z + 1.0;
                        weight[2] = 0.5 * pow ( x, 2 ) + 0.5 * x + 1.0 * y * z + 1.5 * pow ( z, 2 ) - 0.5 * z;
                        weight[3] = -0.5 * pow ( x, 2 ) + 0.5 * x - 1.0 * y - 1.0 * z + 0.5;
                        weight[4] = -4.0 * y * z + 4.0 * y - 6.0 * pow ( z, 2 ) + 6.0 * z - 1.0;
                        weight[5] = -0.5 * pow ( x, 2 ) - 0.5 * x - 1.0 * y - 1.0 * z + 0.5;
                        weight[6] = 0.5 * pow ( x, 2 ) + 1.0 * y * z + 1.5 * pow ( z, 2 ) - 0.5 * z - 0.5;
                        weight[7] = -0.5 * pow ( x, 2 ) - 1.0 * y - 1.0 * z;
                        weight[8] = 0.5 * pow ( x, 2 ) + 1.0 * y * z + 1.5 * pow ( z, 2 ) - 0.5 * z - 0.5;
                        weight[9] = 1.0 * x - 1.0;
                        weight[10] = -1.0 * x - 1.0;
                        weight[11] = -1.0 * x + 1.0;
                        weight[12] = 1.0 * x + 1.0;
                        weight[13] = 0;
                        break;
                    case 4:
                        weight[0] = -1.0 * x * y + 1.0 * y * z + 0.5 * y + 0.5 * z + 0.25;
                        weight[1] = ( x - 1.0 ) * ( 1.0 * y - 0.5 );
                        weight[2] = -1.0 * x * y + 0.5 * x + 1.0 * y * z - 0.5 * y - 0.5 * z + 0.25;
                        weight[3] = y * ( 1.0 * x - 2.0 );
                        weight[4] = y * ( -4.0 * z + 4.0 );
                        weight[5] = 1.0 * x * y;
                        weight[6] = -1.0 * x * y + 1.0 * y * z + 0.5 * y - 0.5 * z - 0.25;
                        weight[7] = ( x - 1.0 ) * ( 1.0 * y + 0.5 );
                        weight[8] = -1.0 * x * y - 0.5 * x + 1.0 * y * z - 0.5 * y + 0.5 * z - 0.25;
                        weight[9] = 1.0 * x - 2.0 * z - 1.0;
                        weight[10] = -1.0 * x + 2.0 * z - 1.0;
                        weight[11] = -1.0 * x + 2.0 * z + 1.0;
                        weight[12] = 1.0 * x - 2.0 * z + 1.0;
                        weight[13] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                if ( xx == 0 && yy == 0 )
                {
                    weight[0] = 0.5 * pow ( x, 2 ) - 0.5 * x + 1.0 * y * z + 1.5 * pow ( z, 2 ) - 0.5 * z;
                    weight[2] = 0.5 * pow ( x, 2 ) + 0.5 * x + 1.0 * y * z + 1.5 * pow ( z, 2 ) - 0.5 * z;
                    weight[7] = -0.5 * pow ( x, 2 ) - 1.0 * y - 1.0 * z;
                }

                //rescaling due to changing coordinates
                weight[0] = 2.0 * weight[0];
                weight[1] = 2.0 * weight[1];
                weight[2] = 2.0 * weight[2];
                weight[3] = 2.0 * weight[3];
                weight[4] = 2.0 * weight[4];
                weight[5] = 2.0 * weight[5];
                weight[6] = 2.0 * weight[6];
                weight[7] = 2.0 * weight[7];
                weight[8] = 2.0 * weight[8];
                weight[9] = 2.0 * weight[9];
                weight[10] = 2.0 * weight[10];
                weight[11] = 2.0 * weight[11];
                weight[12] = 2.0 * weight[12];
                weight[13] = 2.0 * weight[13];

            }

        }

        template<class DataType>
        void FELagrangePyr<DataType>::N_zz ( const Coord& pt, std::vector<DataType>& weight ) const
        {
            assert ( weight.size ( ) == this->get_nb_dof_on_cell ( ) );

            std::setprecision ( 2 * sizeof (DataType ) );

            const DataType x = 2.0 * pt[0] - 1.0;
            const DataType y = 2.0 * pt[1] - 1.0;
            const DataType z = pt[2];

            if ( this->fe_deg_ == 0 )
            {

                weight[0] = 0.0;

            }
            else if ( this->fe_deg_ == 1 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 2:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 3:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    case 4:
                        weight[0] = 0;
                        weight[1] = 0;
                        weight[2] = 0;
                        weight[3] = 0;
                        weight[4] = 0;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 1.0 * weight[0];
                weight[1] = 1.0 * weight[1];
                weight[2] = 1.0 * weight[2];
                weight[3] = 1.0 * weight[3];
                weight[4] = 1.0 * weight[4];

            }
            else if ( this->fe_deg_ == 2 )
            {

                int indicator = -1;

                const DataType xx = trunc ( x * 1e10 );
                const DataType yy = trunc ( y * 1e10 );

                if ( xx > yy && xx < -yy )
                {
                    indicator = 1;
                }
                else if ( xx >= yy && xx >= -yy )
                {
                    indicator = 2;
                }
                else if ( xx < yy && xx > -yy )
                {
                    indicator = 3;
                }
                else if ( xx <= yy && xx <= -yy )
                {
                    indicator = 4;
                }

                switch ( indicator )
                {
                    case 1:
                        weight[0] = 0.5 * pow ( x, 2 ) + 0.5 * x + 0.5 * pow ( y, 2 ) - 3.0 * y * z + 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[1] = 1.0 * y;
                        weight[2] = 0.5 * pow ( x, 2 ) - 0.5 * x + 0.5 * pow ( y, 2 ) - 3.0 * y * z + 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[3] = 1.0 * y;
                        weight[4] = -2.0 * pow ( x, 2 ) - 2.0 * pow ( y, 2 ) + 12.0 * y * z - 6.0 * y - 12.0 * pow ( z, 2 ) + 6.0 * z + 4.0;
                        weight[5] = 1.0 * y;
                        weight[6] = 0.5 * pow ( x, 2 ) - 0.5 * x + 0.5 * pow ( y, 2 ) - 3.0 * y * z + 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[7] = 1.0 * y;
                        weight[8] = 0.5 * pow ( x, 2 ) + 0.5 * x + 0.5 * pow ( y, 2 ) - 3.0 * y * z + 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[9] = -2.0 * x - 2.0;
                        weight[10] = 2.0 * x - 2.0;
                        weight[11] = 2.0 * x - 2.0;
                        weight[12] = -2.0 * x - 2.0;
                        weight[13] = 4.00000000000000;
                        break;
                    case 2:
                        weight[0] = 0.5 * pow ( x, 2 ) + 3.0 * x * z - 0.5 * x + 0.5 * pow ( y, 2 ) - 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[1] = -1.0 * x;
                        weight[2] = 0.5 * pow ( x, 2 ) + 3.0 * x * z - 0.5 * x + 0.5 * pow ( y, 2 ) + 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[3] = -1.0 * x;
                        weight[4] = -2.0 * pow ( x, 2 ) - 12.0 * x * z + 6.0 * x - 2.0 * pow ( y, 2 ) - 12.0 * pow ( z, 2 ) + 6.0 * z + 4.0;
                        weight[5] = -1.0 * x;
                        weight[6] = 0.5 * pow ( x, 2 ) + 3.0 * x * z - 0.5 * x + 0.5 * pow ( y, 2 ) + 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[7] = -1.0 * x;
                        weight[8] = 0.5 * pow ( x, 2 ) + 3.0 * x * z - 0.5 * x + 0.5 * pow ( y, 2 ) - 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[9] = 2.0 * y - 2.0;
                        weight[10] = -2.0 * y - 2.0;
                        weight[11] = -2.0 * y - 2.0;
                        weight[12] = 2.0 * y - 2.0;
                        weight[13] = 4.00000000000000;
                        break;
                    case 3:
                        weight[0] = 0.5 * pow ( x, 2 ) - 0.5 * x + 0.5 * pow ( y, 2 ) + 3.0 * y * z - 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[1] = -1.0 * y;
                        weight[2] = 0.5 * pow ( x, 2 ) + 0.5 * x + 0.5 * pow ( y, 2 ) + 3.0 * y * z - 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[3] = -1.0 * y;
                        weight[4] = -2.0 * pow ( x, 2 ) - 2.0 * pow ( y, 2 ) - 12.0 * y * z + 6.0 * y - 12.0 * pow ( z, 2 ) + 6.0 * z + 4.0;
                        weight[5] = -1.0 * y;
                        weight[6] = 0.5 * pow ( x, 2 ) + 0.5 * x + 0.5 * pow ( y, 2 ) + 3.0 * y * z - 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[7] = -1.0 * y;
                        weight[8] = 0.5 * pow ( x, 2 ) - 0.5 * x + 0.5 * pow ( y, 2 ) + 3.0 * y * z - 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[9] = 2.0 * x - 2.0;
                        weight[10] = -2.0 * x - 2.0;
                        weight[11] = -2.0 * x - 2.0;
                        weight[12] = 2.0 * x - 2.0;
                        weight[13] = 4.00000000000000;
                        break;
                    case 4:
                        weight[0] = 0.5 * pow ( x, 2 ) - 3.0 * x * z + 0.5 * x + 0.5 * pow ( y, 2 ) + 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[1] = 1.0 * x;
                        weight[2] = 0.5 * pow ( x, 2 ) - 3.0 * x * z + 0.5 * x + 0.5 * pow ( y, 2 ) - 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[3] = 1.0 * x;
                        weight[4] = -2.0 * pow ( x, 2 ) + 12.0 * x * z - 6.0 * x - 2.0 * pow ( y, 2 ) - 12.0 * pow ( z, 2 ) + 6.0 * z + 4.0;
                        weight[5] = 1.0 * x;
                        weight[6] = 0.5 * pow ( x, 2 ) - 3.0 * x * z + 0.5 * x + 0.5 * pow ( y, 2 ) - 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[7] = 1.0 * x;
                        weight[8] = 0.5 * pow ( x, 2 ) - 3.0 * x * z + 0.5 * x + 0.5 * pow ( y, 2 ) + 0.5 * y + 3.0 * pow ( z, 2 ) - 1.5 * z;
                        weight[9] = -2.0 * y - 2.0;
                        weight[10] = 2.0 * y - 2.0;
                        weight[11] = 2.0 * y - 2.0;
                        weight[12] = -2.0 * y - 2.0;
                        weight[13] = 4.00000000000000;
                        break;
                    default:
                        std::cerr << "Unexpected coordinate " << std::endl;
                        exit ( -1 );
                }

                //rescaling due to changing coordinates
                weight[0] = 1.0 * weight[0];
                weight[1] = 1.0 * weight[1];
                weight[2] = 1.0 * weight[2];
                weight[3] = 1.0 * weight[3];
                weight[4] = 1.0 * weight[4];
                weight[5] = 1.0 * weight[5];
                weight[6] = 1.0 * weight[6];
                weight[7] = 1.0 * weight[7];
                weight[8] = 1.0 * weight[8];
                weight[9] = 1.0 * weight[9];
                weight[10] = 1.0 * weight[10];
                weight[11] = 1.0 * weight[11];
                weight[12] = 1.0 * weight[12];
                weight[13] = 1.0 * weight[13];

            }

        }

        template class FELagrangePyr<double>;
        template class FELagrangePyr<float>;

    } // namespace doffem
} // namespace hiflow

