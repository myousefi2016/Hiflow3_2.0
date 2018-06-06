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

#ifndef _DOF_NUMBERING_STRATEGY_H_
#    define _DOF_NUMBERING_STRATEGY_H_

#    include <vector>
#    include <map>

#    include "mesh/mesh.h"

#    include "dof/fe_interface_pattern.h"
#    include "dof/dof_interpolation_pattern.h"
#    include "dof/dof_interpolation.h"
#    include "dof/dof_fem_types.h"

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        class FEManager;

        template<class DataType>
        class FEType;

        template<class DataType>
        class DegreeOfFreedom;

        /// Enumeration of different DoF ordering strategies. HIFLOW_CLASSIC refers to the 
        /// DoF numbering as always done in HiFlow3. The other two options allow permutations 
        /// of the classic numbering by means of the Cuthill-McKee and the King method, 
        /// respectively.

        enum DOF_ORDERING
        {
            HIFLOW_CLASSIC, CUTHILL_MCKEE, KING
        };

        /// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

        template<class DataType>
        class NumberingStrategy
        {
          public:

            typedef std::vector<DataType> Coord;

            /// Constructor

            NumberingStrategy ( )
            {
            }
            /// Destructor

            virtual ~NumberingStrategy ( )
            {
            }

            /// Initialization. Here, the critical member variables of DegreeOfFreedom are being
            /// given to NumberingStrategy, such that the implementation class of the function
            /// void number() can calculate all neccessary information and store it in these variables
            void initialize ( DegreeOfFreedom<DataType>& dof, std::vector<DofID>& numer,
                              std::vector<std::vector<int> >& numer_offsets_cell_varloc,
                              DofInterpolation& dof_interpolation );

            /// Kernel of numbering strategy. Here the user can specify in an inherited class his wishes
            /// for some numbering procedure, when dealing for example with varying finite element types
            /// @param[in] order Ordering strategy for DoFs.
            virtual void number ( DOF_ORDERING order = HIFLOW_CLASSIC ) = 0;

            /// Helper function which permutes data within the interpolation container and numer_ field
            void apply_permutation ( const std::vector<DofID>& permutation, const std::string& = "" );

          protected:

            /// Topological dimension
            int tdim_;
            /// Total number of variables
            int nvar_;

            /// The DegreeOfFreedom class used throughout the numeration procedure
            DegreeOfFreedom<DataType>* dof_;

            /// Const pointer to mesh
            const mesh::Mesh* mesh_;

            /// FEManager on the given mesh
            FEManager<DataType> const* fe_manager_;

            /// Holds DoF IDs, needed for mapl2g
            std::vector<DofID>* numer_;

            /// Offset container for numer_, needed for mapl2g
            std::vector<std::vector<int> >* numer_offsets_cell_varloc_;

            /// Interpolation Container, which stores the interpolation weights
            DofInterpolation* dof_interpolation_;

            /// Update number_of_dofs_total_ and number_of_dofs_for_var_
            void update_number_of_dofs ( const std::string& = "" );

            /// Total number of dofs for all variables
            int number_of_dofs_total_;

            /// Total number of dofs per variable
            std::vector<int> number_of_dofs_for_var_;

            /// Get vertex points of mesh entity.
            void get_points ( const mesh::Entity& entity, std::vector<Coord>& points );

        };

    }
}

#endif
