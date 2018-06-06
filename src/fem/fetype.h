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

#ifndef __FEM_FETYPE_H_
#    define __FEM_FETYPE_H_

#    include <string>
#    include <vector>
#    include <cassert>
#    include <iostream>

#    include "mesh/cell_type.h"
#    include "dof/dof_fem_types.h"

namespace hiflow
{
    namespace doffem
    {

        ///
        /// \class FEType fetype.h
        /// \brief Ancestor class of different Finite Elements
        /// \author Michael Schick<br>Martin Baumann
        ///
        /// \todo rename coord_on_subentity_ to coord_on_entity_ => also coord can be treated

        template<class DataType>
        class FEType
        {
          public:

            /// Collection of information about implemented Finite Elements

            enum FiniteElement
            {
                NOT_SET = 0,
                LAGRANGE_LINE,
                LAGRANGE_TRI,
                LAGRANGE_QUAD,
                LAGRANGE_TET,
                LAGRANGE_HEX,
                LAGRANGE_PYR
            };

            /// Collection of abstract fe ansatz classes

            enum FEAnsatz
            {
                LAGRANGE
            };

            typedef std::vector<DataType> Coord;

            /// Default Constructor
            FEType ( );

            /// Default Destructor
            virtual ~FEType ( );

            /// Setting the coordinates for the vertices on the reference cell
            void set_ref_vtx_coords ( std::vector<Coord> ref_coords );

            /// Setting Id of the Finite Element (from FEInstance)

            void set_id ( int Id )
            {
                instance_id_ = Id;
            }

            /// Set the polynomial degree of the ansatz and starting init_coord()
            void set_fe_deg ( int degree );

            /// Get the polynomial degree of the ansatz.
            inline int get_fe_deg ( ) const;

            /// Get information (name) about the used FE Type
            virtual std::string get_name ( ) const = 0;

            /// Get ID of Finite Element, which is stored in enum FiniteElement
            inline virtual FiniteElement get_my_id ( ) const;
            /// Get better Id of the Finite Element
            inline int get_global_id ( ) const;

            /// Get Ansatz of Finite Element, which is stored in enum FEAnsatz
            inline virtual FEAnsatz get_my_ansatz ( ) const;

            /// Get information if this finite element was initialized (true) or not (false)
            inline bool get_init_status ( ) const;

            /// Total number of dofs on the cell
            int get_nb_dof_on_cell ( ) const;

            /// Number of subentities
            int get_nb_subentity ( int tdim ) const;

            /// Number of dofs on a subentity
            int get_nb_dof_on_subentity ( int tdim, int index ) const;

            /// Get information about the coordinates on a specific subentity
            std::vector<Coord> const& get_coord_on_subentity ( int tdim, int index ) const;

            /// Get information about the DofIDs on a specific subentity
            std::vector<DofID> const& get_dof_on_subentity ( int tdim, int index ) const;

            /// Get information about the coordinates on reference cell
            std::vector<Coord> const& get_coord_on_cell ( ) const;

            /// Initialize the Finite Element
            void init ( );

            /// For given point, get values of all shapefunctions on reference cell
            virtual void N ( const Coord& pt, std::vector<DataType>& weight ) const = 0;

            /// For given coordinates, get x - derivative of all shapefunctions on reference cell
            virtual void N_x ( const Coord& pt, std::vector<DataType>& weight ) const = 0;
            /// For given coordinates, get y - derivative of all shapefunctions on reference cell
            virtual void N_y ( const Coord& pt, std::vector<DataType>& weight ) const = 0;
            /// For given coordinates, get z - derivative of all shapefunctions on reference cell
            virtual void N_z ( const Coord& pt, std::vector<DataType>& weight ) const = 0;

            /// For given coordinates, get xx - derivative of all shapefunctions on reference cell
            virtual void N_xx ( const Coord& pt, std::vector<DataType>& weight ) const = 0;
            /// For given coordinates, get xy - derivative of all shapefunctions on reference cell
            virtual void N_xy ( const Coord& pt, std::vector<DataType>& weight ) const = 0;
            /// For given coordinates, get xz - derivative of all shapefunctions on reference cell
            virtual void N_xz ( const Coord& pt, std::vector<DataType>& weight ) const = 0;
            /// For given coordinates, get yy - derivative of all shapefunctions on reference cell
            virtual void N_yy ( const Coord& pt, std::vector<DataType>& weight ) const = 0;
            /// For given coordinates, get yz - derivative of all shapefunctions on reference cell
            virtual void N_yz ( const Coord& pt, std::vector<DataType>& weight ) const = 0;
            /// For given coordinates, get zz - derivative of all shapefunctions on reference cell
            virtual void N_zz ( const Coord& pt, std::vector<DataType>& weight ) const = 0;

            /// Operators needed to be able to create maps where FEType is
            /// used as key. \see FEInterfacePattern::operator < (const FEInterfacePattern& test) const
            /// Comparison by protected variable my_id_ and fe_deg_ .
            virtual bool operator== ( const FEType& fe_slave ) const;

            /// Operators needed to be able to create maps where FEType is
            /// used as key. \see FEInterfacePattern::operator < (const FEInterfacePattern& test) const
            /// Comparison by protected variable my_id_ and fe_deg_ .
            virtual bool operator< ( const FEType& fe_slave ) const;

          protected:

            /// Better Id
            int instance_id_;

            /// Status information if fetype was initialized
            bool init_status_;

            /// Compute coordinates on reference cell
            virtual void init_coord ( ) = 0;

            /// filter the DoF points for the cell's faces, edges or vertices
            void init_dofs_on_subentities ( );

            /// check whether a point 'point' is located within a subentity (point, line, face), defined by 'points'
            bool is_point_on_subentity ( Coord point, const std::vector<Coord>& points ) const;

            /// Storing coordinates in vector by lexicographical numbering strategy
            std::vector<Coord> coord_;

            /// Storing coordinates of reference cell's vertices
            std::vector<Coord> coord_vertices_;

            /// Storing an instance of the reference cell
            const mesh::CellType* ref_cell_;

            /// coordinates stored for specific subentities (point, line, face)
            std::vector<std::vector<std::vector<Coord> > > coord_on_subentity_;

            /// DoF ID stored for specific subentities (point, line, face)
            std::vector<std::vector<std::vector<DofID> > > dof_on_subentity_;

            /// Id of a Finite Element \see FEType::FiniteElement
            FiniteElement my_id_;
            /// Ansatz of a Finite Element \see FEType::FEAnsatz
            FEAnsatz my_ansatz_;

            /// Topological dimension
            int tdim_;

            /// Finite Element degree
            int fe_deg_;

        };

        //-------------- INLINE FUNCTIONS FOR FETYPE----------------------

        template<class DataType>
        typename FEType<DataType>::FiniteElement FEType<DataType>::get_my_id ( ) const
        {
            return my_id_;
        }

        template<class DataType>
        int FEType<DataType>::get_global_id ( ) const
        {
            return instance_id_;
        }

        template<class DataType>
        typename FEType<DataType>::FEAnsatz FEType<DataType>::get_my_ansatz ( ) const
        {
            return my_ansatz_;
        }

        template<class DataType>
        bool FEType<DataType>::get_init_status ( ) const
        {
            return init_status_;
        }

        template<class DataType>
        int FEType<DataType>::get_fe_deg ( ) const
        {
            return fe_deg_;
        }

    } // namespace doffem
} // namespace hiflow
#endif
