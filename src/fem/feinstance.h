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

#ifndef __FEM_FEINSTANCE_H_
#    define __FEM_FEINSTANCE_H_

#    include "dof/dof_fem_types.h"
#    include <vector>

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        class FELagrangeLine;

        template<class DataType>
        class FELagrangeTri;

        template<class DataType>
        class FELagrangeQuad;

        template<class DataType>
        class FELagrangeTet;

        template<class DataType>
        class FELagrangeHex;

        template<class DataType>
        class FELagrangePyr;

        ///
        /// \class FEInstance feinstance.h
        /// \brief Holds instances of different FE ansatz functions
        /// \author Michael Schick<br>Martin Baumann
        ///

        template<class DataType>
        class FEInstance
        {
          public:

            typedef std::vector<DataType> Coord;

            /// Default Constructor
            FEInstance ( );
            /// Default Destructor
            ~FEInstance ( );

            /// Returns pointer to an instance of a Lagrangian Line element of given degree
            FELagrangeLine<DataType>* get_lagrange_line_element ( int degree );
            /// Returns pointer to an instance of a Lagrangian Triangle element of given degree
            FELagrangeTri<DataType>* get_lagrange_tri_element ( int degree );
            /// Returns pointer to an instance of a Lagrangian Quadrilateral element of given degree
            FELagrangeQuad<DataType>* get_lagrange_quad_element ( int degree );
            /// Returns pointer to an instance of a Lagrangian Tetrahedron element of given degree
            FELagrangeTet<DataType>* get_lagrange_tet_element ( int degree );
            /// Returns pointer to an instance of a Lagrangian Hexahedron element of given degree
            FELagrangeHex<DataType>* get_lagrange_hex_element ( int degree );
            /// Returns pointer to an instance of a Lagrangian Pyramid element of given degree
            FELagrangePyr<DataType>* get_lagrange_pyr_element ( int degree );

            /// Returns pointer to an instance of a Lagrangian Line element of given degree (const)
            FELagrangeLine<DataType>* get_const_lagrange_line_element ( int degree ) const;
            /// Returns pointer to an instance of a Lagrangian Triangle element of given degree (const)
            FELagrangeTri<DataType>* get_const_lagrange_tri_element ( int degree ) const;
            /// Returns pointer to an instance of a Lagrangian Quadrilateral element of given degree (const)
            FELagrangeQuad<DataType>* get_const_lagrange_quad_element ( int degree ) const;
            /// Returns pointer to an instance of a Lagrangian Tetrahedron element of given degree (const)
            FELagrangeTet<DataType>* get_const_lagrange_tet_element ( int degree ) const;
            /// Returns pointer to an instance of a Lagrangian Hexahedron element of given degree (const)
            FELagrangeHex<DataType>* get_const_lagrange_hex_element ( int degree ) const;
            /// Returns pointer to an instance of a Lagrangian Pyramid element of given degree (const)
            FELagrangePyr<DataType>* get_const_lagrange_pyr_element ( int degree ) const;

            /// Return status information about FEInstance
            void get_status ( ) const;
            /// Return number of initialized Finite Elements (use only after complete initialization!)

            int nfe ( ) const
            {
                return id_;
            }

            /// Clearing all stored finite elements and reseting FEInstance
            void clear ( );

            /// Initialize the coordinates on the reference line for vertices
            void init_line_elements ( const std::vector<Coord>& ref_coords );
            /// Initialize the coordinates on the reference triangle for vertices
            void init_triangle_elements ( const std::vector<Coord>& ref_coords );
            /// Initialize the coordinates on the reference quadrilateral for vertices
            void init_quadrilateral_elements ( const std::vector<Coord>& ref_coords );
            /// Initialize the coordinates on the reference tetrahedron for vertices
            void init_tetrahedron_elements ( const std::vector<Coord>& ref_coords );
            /// Initialize the coordinates on the reference hexahedron for vertices
            void init_hexahedron_elements ( const std::vector<Coord>& ref_coords );
            /// Initialize the coordinates on the reference pyramid for vertices
            void init_pyramid_elements ( const std::vector<Coord>& ref_coords );

          private:

            /// Id counter
            int id_;

            /// Stored Lagrange instances on a Line
            std::vector<FELagrangeLine<DataType>* > fe_lag_line_;

            /// Stored Lagrange instances on a Triangle
            std::vector<FELagrangeTri<DataType>* > fe_lag_tri_;

            /// Stored Lagrange instances on a Quadrilateral
            std::vector<FELagrangeQuad<DataType>* > fe_lag_quad_;

            /// Stored Lagrange instances on a Tetrahedron
            std::vector<FELagrangeTet<DataType>* > fe_lag_tet_;

            /// Stored Lagrange instances on a Hexahedron
            std::vector<FELagrangeHex<DataType>* > fe_lag_hex_;

            /// Stored Lagrange instances on a Pyramid
            std::vector<FELagrangePyr<DataType>* > fe_lag_pyr_;

        };

    }
} // namespace hiflow
#endif
