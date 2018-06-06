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

#include "feinstance.h"

#include "felagrange_line.h"
#include "felagrange_tri.h"
#include "felagrange_quad.h"
#include "felagrange_tet.h"
#include "felagrange_hex.h"
#include "felagrange_pyr.h"

#include <iostream>

/// \author Michael Schick<br>Martin Baumann

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        FEInstance<DataType>::FEInstance ( )
        {
            id_ = 0;
        }

        template<class DataType>
        FEInstance<DataType>::~FEInstance ( )
        {
            clear ( );
        }

        /// \details Every allocated instance of a finite element is being
        ///          deleted and the vectors storing this information are cleared

        template<class DataType>
        void FEInstance<DataType>::clear ( )
        {
            for ( size_t i = 0, e_i = fe_lag_line_.size ( ); i != e_i; ++i )
                delete fe_lag_line_[i];

            for ( size_t i = 0, e_i = fe_lag_tri_.size ( ); i != e_i; ++i )
                delete fe_lag_tri_[i];

            for ( size_t i = 0, e_i = fe_lag_quad_.size ( ); i != e_i; ++i )
                delete fe_lag_quad_[i];

            for ( size_t i = 0, e_i = fe_lag_tet_.size ( ); i != e_i; ++i )
                delete fe_lag_tet_[i];

            for ( size_t i = 0, e_i = fe_lag_hex_.size ( ); i != e_i; ++i )
                delete fe_lag_hex_[i];

            for ( size_t i = 0, e_i = fe_lag_pyr_.size ( ); i != e_i; ++i )
                delete fe_lag_pyr_[i];

            fe_lag_line_.clear ( );
            fe_lag_tri_.clear ( );
            fe_lag_quad_.clear ( );
            fe_lag_tet_.clear ( );
            fe_lag_hex_.clear ( );
            fe_lag_pyr_.clear ( );

            id_ = 0;
        }

        template<class DataType>
        void FEInstance<DataType>::init_line_elements ( const std::vector<Coord>& ref_coords )
        {
            assert ( ref_coords.size ( ) == 2 );

            for ( size_t deg = 0, e_deg = fe_lag_line_.size ( ); deg != e_deg; ++deg )
            {
                if ( !fe_lag_line_[deg]->get_init_status ( ) )
                {
                    fe_lag_line_[deg]->set_ref_vtx_coords ( ref_coords );
                    fe_lag_line_[deg]->init ( );
                    fe_lag_line_[deg]->set_id ( id_ );
                    ++id_;
                }
            }
        }

        template<class DataType>
        void FEInstance<DataType>::init_triangle_elements ( const std::vector<Coord>& ref_coords )
        {
            assert ( ref_coords.size ( ) == 3 );

            for ( size_t deg = 0, e_deg = fe_lag_tri_.size ( ); deg != e_deg; ++deg )
            {
                if ( !fe_lag_tri_[deg]->get_init_status ( ) )
                {
                    fe_lag_tri_[deg]->set_ref_vtx_coords ( ref_coords );
                    fe_lag_tri_[deg]->init ( );
                    fe_lag_tri_[deg]->set_id ( id_ );
                    ++id_;
                }
            }
        }

        template<class DataType>
        void FEInstance<DataType>::init_quadrilateral_elements ( const std::vector<Coord>& ref_coords )
        {
            assert ( ref_coords.size ( ) == 4 );

            for ( size_t deg = 0, e_deg = fe_lag_quad_.size ( ); deg != e_deg; ++deg )
            {
                if ( !fe_lag_quad_[deg]->get_init_status ( ) )
                {
                    fe_lag_quad_[deg]->set_ref_vtx_coords ( ref_coords );
                    fe_lag_quad_[deg]->init ( );
                    fe_lag_quad_[deg]->set_id ( id_ );
                    ++id_;
                }
            }
        }

        template<class DataType>
        void FEInstance<DataType>::init_tetrahedron_elements ( const std::vector<Coord>& ref_coords )
        {
            assert ( ref_coords.size ( ) == 4 );

            for ( size_t deg = 0, e_deg = fe_lag_tet_.size ( ); deg != e_deg; ++deg )
            {
                if ( !fe_lag_tet_[deg]->get_init_status ( ) )
                {
                    fe_lag_tet_[deg]->set_ref_vtx_coords ( ref_coords );
                    fe_lag_tet_[deg]->init ( );
                    fe_lag_tet_[deg]->set_id ( id_ );
                    ++id_;
                }
            }
        }

        template<class DataType>
        void FEInstance<DataType>::init_hexahedron_elements ( const std::vector<Coord>& ref_coords )
        {
            assert ( ref_coords.size ( ) == 8 );

            for ( size_t deg = 0, e_deg = fe_lag_hex_.size ( ); deg != e_deg; ++deg )
            {
                if ( !fe_lag_hex_[deg]->get_init_status ( ) )
                {
                    fe_lag_hex_[deg]->set_ref_vtx_coords ( ref_coords );
                    fe_lag_hex_[deg]->init ( );
                    fe_lag_hex_[deg]->set_id ( id_ );
                    ++id_;
                }
            }
        }

        template<class DataType>
        void FEInstance<DataType>::init_pyramid_elements ( const std::vector<Coord>& ref_coords )
        {
            assert ( ref_coords.size ( ) == 5 );

            for ( size_t deg = 0, e_deg = fe_lag_pyr_.size ( ); deg != e_deg; ++deg )
            {
                if ( !fe_lag_pyr_[deg]->get_init_status ( ) )
                {
                    fe_lag_pyr_[deg]->set_ref_vtx_coords ( ref_coords );
                    fe_lag_pyr_[deg]->init ( );
                    fe_lag_pyr_[deg]->set_id ( id_ );
                    ++id_;
                }
            }
        }

        /// \details Checks wether a lagrangian finite element on a line with
        ///          given degree is already existing. If not, then all elements up to
        ///          the given degree are created.

        template<class DataType>
        FELagrangeLine<DataType>* FEInstance<DataType>::get_lagrange_line_element ( int degree )
        {
            if ( fe_lag_line_.size ( ) <= degree )
            {
                // add all FETypes up to degree
                for ( int i = fe_lag_line_.size ( ); i < degree + 1; ++i )
                {
                    FELagrangeLine<DataType>* dummy_lag = new FELagrangeLine<DataType>;
                    dummy_lag->set_fe_deg ( i );
                    fe_lag_line_.push_back ( dummy_lag );
                }
            }

            return fe_lag_line_[degree];
        }

        /// \details Checks wether a lagrangian finite element on a triangle with
        ///          given degree is already existing. If not, then all elements up to
        ///          the given degree are created.

        template<class DataType>
        FELagrangeTri<DataType>* FEInstance<DataType>::get_lagrange_tri_element ( int degree )
        {
            if ( fe_lag_tri_.size ( ) <= degree )
            {
                // add all FETypes up to degree
                for ( int i = fe_lag_tri_.size ( ); i < degree + 1; ++i )
                {
                    FELagrangeTri<DataType>* dummy_lag = new FELagrangeTri<DataType>;
                    dummy_lag->set_fe_deg ( i );
                    fe_lag_tri_.push_back ( dummy_lag );
                }
            }

            return fe_lag_tri_[degree];
        }

        /// \details Checks wether a lagrangian finite element on a quadrilateral with
        ///          given degree is already existing. If not, then all elements up to
        ///          the given degree are created.

        template<class DataType>
        FELagrangeQuad<DataType>* FEInstance<DataType>::get_lagrange_quad_element ( int degree )
        {
            if ( fe_lag_quad_.size ( ) <= degree )
            {
                // add all FETypes up to degree
                for ( int i = fe_lag_quad_.size ( ); i < degree + 1; ++i )
                {
                    FELagrangeQuad<DataType>* dummy_lag = new FELagrangeQuad<DataType>;
                    dummy_lag->set_fe_deg ( i );
                    fe_lag_quad_.push_back ( dummy_lag );
                }
            }

            return fe_lag_quad_[degree];
        }

        /// \details Checks wether a lagrangian finite element on a tetrahedron with
        ///          given degree is already existing. If not, then all elements up to
        ///          the given degree are created.

        template<class DataType>
        FELagrangeTet<DataType>* FEInstance<DataType>::get_lagrange_tet_element ( int degree )
        {
            if ( fe_lag_tet_.size ( ) <= degree )
            {
                // add all FETypes up to degree
                for ( int i = fe_lag_tet_.size ( ); i < degree + 1; ++i )
                {
                    FELagrangeTet<DataType>* dummy_lag = new FELagrangeTet<DataType>;
                    dummy_lag->set_fe_deg ( i );
                    fe_lag_tet_.push_back ( dummy_lag );
                }
            }

            return fe_lag_tet_[degree];
        }

        /// \details Checks wether a lagrangian finite element on a hexahedron with
        ///          given degree is already existing. If not, then all elements up to
        ///          the given degree are created.

        template<class DataType>
        FELagrangeHex<DataType>* FEInstance<DataType>::get_lagrange_hex_element ( int degree )
        {
            if ( fe_lag_hex_.size ( ) <= degree )
            {
                // add all FETypes up to degree
                for ( int i = fe_lag_hex_.size ( ); i < degree + 1; ++i )
                {
                    FELagrangeHex<DataType>* dummy_lag = new FELagrangeHex<DataType>;
                    dummy_lag->set_fe_deg ( i );
                    fe_lag_hex_.push_back ( dummy_lag );
                }
            }

            return fe_lag_hex_[degree];
        }

        /// \details Checks wether a lagrangian finite element on a pyramid with
        ///          given degree is already existing. If not, then all elements up to
        ///          the given degree are created.

        template<class DataType>
        FELagrangePyr<DataType>* FEInstance<DataType>::get_lagrange_pyr_element ( int degree )
        {
            if ( fe_lag_pyr_.size ( ) <= degree )
            {
                // add all FETypes up to degree
                for ( int i = fe_lag_pyr_.size ( ); i < degree + 1; ++i )
                {
                    FELagrangePyr<DataType>* dummy_lag = new FELagrangePyr<DataType>;
                    dummy_lag->set_fe_deg ( i );
                    fe_lag_pyr_.push_back ( dummy_lag );
                }
            }

            return fe_lag_pyr_[degree];
        }

        template<class DataType>
        FELagrangeLine<DataType>* FEInstance<DataType>::get_const_lagrange_line_element ( int degree ) const
        {
            return fe_lag_line_[degree];
        }

        template<class DataType>
        FELagrangeTri<DataType>* FEInstance<DataType>::get_const_lagrange_tri_element ( int degree ) const
        {
            return fe_lag_tri_[degree];
        }

        template<class DataType>
        FELagrangeQuad<DataType>* FEInstance<DataType>::get_const_lagrange_quad_element ( int degree ) const
        {
            return fe_lag_quad_[degree];
        }

        template<class DataType>
        FELagrangeTet<DataType>* FEInstance<DataType>::get_const_lagrange_tet_element ( int degree ) const
        {
            return fe_lag_tet_[degree];
        }

        template<class DataType>
        FELagrangeHex<DataType>* FEInstance<DataType>::get_const_lagrange_hex_element ( int degree ) const
        {
            return fe_lag_hex_[degree];
        }

        template<class DataType>
        FELagrangePyr<DataType>* FEInstance<DataType>::get_const_lagrange_pyr_element ( int degree ) const
        {
            return fe_lag_pyr_[degree];
        }

        template<class DataType>
        void FEInstance<DataType>::get_status ( ) const
        {
            // std::cout << "FEInstance size: " << fe_lag_hex_.size() << std::endl;

            std::cout << "  Lines" << std::endl;
            for ( int i = 0; i < fe_lag_line_.size ( ); ++i )
                std::cout << "\t" << i << "\t" << fe_lag_line_[i]->get_name ( )
                << " (" << fe_lag_line_[i]->get_fe_deg ( ) << ")"
                << "\t@ " << get_const_lagrange_line_element ( i ) << std::endl;

            std::cout << "  Triangles" << std::endl;
            for ( int i = 0; i < fe_lag_tri_.size ( ); ++i )
                std::cout << "\t" << i << "\t" << fe_lag_tri_[i]->get_name ( )
                << " (" << fe_lag_tri_[i]->get_fe_deg ( ) << ")"
                << "\t@ " << get_const_lagrange_tri_element ( i ) << std::endl;

            std::cout << "  Quadrilaterals" << std::endl;
            for ( int i = 0; i < fe_lag_quad_.size ( ); ++i )
                std::cout << "\t" << i << "\t" << fe_lag_quad_[i]->get_name ( )
                << " (" << fe_lag_quad_[i]->get_fe_deg ( ) << ")"
                << "\t@ " << get_const_lagrange_quad_element ( i ) << std::endl;

            std::cout << "  Tetrahedrons" << std::endl;
            for ( int i = 0; i < fe_lag_tet_.size ( ); ++i )
                std::cout << "\t" << i << "\t" << fe_lag_tet_[i]->get_name ( )
                << " (" << fe_lag_tet_[i]->get_fe_deg ( ) << ")"
                << "\t@ " << get_const_lagrange_tet_element ( i ) << std::endl;

            std::cout << "  Hexahedrons" << std::endl;
            for ( int i = 0; i < fe_lag_hex_.size ( ); ++i )
                std::cout << "\t" << i << "\t" << fe_lag_hex_[i]->get_name ( )
                << " (" << fe_lag_hex_[i]->get_fe_deg ( ) << ")"
                << "\t@ " << get_const_lagrange_hex_element ( i ) << std::endl;

            std::cout << "  Pyramid" << std::endl;
            for ( int i = 0; i < fe_lag_pyr_.size ( ); ++i )
                std::cout << "\t" << i << "\t" << fe_lag_pyr_[i]->get_name ( )
                << " (" << fe_lag_pyr_[i]->get_fe_deg ( ) << ")"
                << "\t@ " << get_const_lagrange_pyr_element ( i ) << std::endl;

        }

        template class FEInstance<double>;
        template class FEInstance<float>;

    }
} // namespace hiflow
