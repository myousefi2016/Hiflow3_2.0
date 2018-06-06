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

#ifndef __FEM_FELAGRANGE_QUAD_H_
#    define __FEM_FELAGRANGE_QUAD_H_

#    include "felagrange.h"

namespace hiflow
{
    namespace doffem
    {

        ///
        /// \class FELagrangeQuad felagrange_quad.h
        /// \brief Lagrangian Finite Element on a Quadrilateral
        /// \author Michael Schick<br>Martin Baumann
        ///

        template<class DataType>
        class FELagrangeQuad : public FELagrange<DataType>
        {
          public:

            typedef std::vector<DataType> Coord;

            /// Default Constructor
            FELagrangeQuad ( );

            /// Default Destructor
            ~FELagrangeQuad ( );

            std::string get_name ( ) const
            {
                return "LagrangeQuadrilateral";
            }

            void N ( const Coord& pt, std::vector<DataType>& weight ) const;

            void N_x ( const Coord& pt, std::vector<DataType>& weight ) const;
            void N_y ( const Coord& pt, std::vector<DataType>& weight ) const;
            void N_z ( const Coord& pt, std::vector<DataType>& weight ) const;

            void N_xx ( const Coord& pt, std::vector<DataType>& weight ) const;
            void N_xy ( const Coord& pt, std::vector<DataType>& weight ) const;
            void N_xz ( const Coord& pt, std::vector<DataType>& weight ) const;
            void N_yy ( const Coord& pt, std::vector<DataType>& weight ) const;
            void N_yz ( const Coord& pt, std::vector<DataType>& weight ) const;
            void N_zz ( const Coord& pt, std::vector<DataType>& weight ) const;

          protected:

            void init_coord ( );

          private:

            /// Index ordering from plane (x=i, y=j) to vector index
            int ij2ind ( int i, int j ) const;

        };

    } // namespace doffem
} // namespace hiflow
#endif
