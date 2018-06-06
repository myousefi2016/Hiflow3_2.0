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

#ifndef _DOF_DEGREE_OF_FREEDOM_H_
#    define _DOF_DEGREE_OF_FREEDOM_H_

#    include <vector>
#    include <map>

#    include "mesh/mesh.h"
#    include "dof/dof_fem_types.h"

#    include "numbering_strategy.h"

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        class FEManager;

        template<class DataType>
        class FEType;

        /// \author Michael Schick<br>Martin Baumann

        template<class DataType>
        class DegreeOfFreedom
        {
          public:

            typedef std::vector<DataType> Coord;
            /// Constructor
            DegreeOfFreedom ( );

            /// Destructor
            virtual ~DegreeOfFreedom ( );

            /// Set given mesh to a constant pointer
            void set_mesh ( const mesh::Mesh* mesh );

            /// Numbering of the given mesh with given finite element ansatz
            /// @param[in] order Ordering strategy for DoFs.
            virtual void number ( DOF_ORDERING order = HIFLOW_CLASSIC );

            /// \brief Mapping local 2 global. Local DofId is local to a specific
            ///        cell and a specific variable. Global DofId is the global Id
            ///        on the whole given mesh
            DofID mapl2g ( int var, int cell_index, DofID local_id ) const;

            /// Setting the FEManager for the given mesh
            void set_fe_manager ( FEManager<DataType> const* manager );

            /// Get the FEManager for the given mesh
            FEManager<DataType> const& get_fe_manager ( ) const;

            /// Initialize numbering strategy for the domain given by mesh
            void init_numbering_strategy ( const std::string& number_strategy = "" );

            /// Get the mesh

            const mesh::Mesh& get_mesh ( ) const
            {
                return *mesh_;
            }

            /// Get the DoF Interpolation
            DofInterpolation& dof_interpolation ( );
            const DofInterpolation& dof_interpolation ( ) const;

            /// Get the coordinates of dofs on a mesh cell for a specific variable and cell index
            void get_coord_on_cell ( int var, int cell_index, std::vector<Coord>& coords ) const;

            /// Get the global DofIds on a specific mesh cell and variable
            void get_dofs_on_cell ( int var, int cell_index, std::vector<DofID>& ids ) const;

            /// Get the number of dofs for a specific variable on a specific cell
            int get_nb_dofs_on_cell ( int var, int cell_index ) const;

            /// Get the total number of dofs on a specific cell (including all variables)
            int get_nb_dofs_on_cell ( int cell_index ) const;

            /// Get the number of dofs on the boundary for a specific variable on a specific cell
            int get_nb_dofs_on_subentity ( int var, int tdim, int cell_index ) const;

            /// Get the total number of dofs on the boundary on a specific cell (including all variables)
            int get_nb_dofs_on_subentity ( int tdim, int cell_index ) const;

            /// Get the coordinates of dofs on a subentity (point, edge, fase) of a cell
            void get_coord_on_subentity ( int var, int cell_index, int tdim,
                                          int sindex, std::vector<Coord>& coords ) const;

            /// Get the global DofIds on a subentity (point, edge, fase) of a cell
            void get_dofs_on_subentity ( int var, int cell_index, int tdim,
                                         int sindex, std::vector<DofID>& ids ) const;

            /// Get the number of dofs on the whole mesh for a specific variable
            int get_nb_dofs ( int var ) const;

            /// Get the number of dofs on the whole mesh for all variables
            int get_nb_dofs ( ) const;

            /// Apply permutation of DoF IDs
            void apply_permutation ( const std::vector<DofID>& permutation, const std::string& = "" );

            /// Printing information about the numbering field
            void print_numer ( ) const;

          protected:
            /// Numbering strategy for (sub)domain numbering
            NumberingStrategy<DataType>* number_strategy_;

            /// Topological dimension
            int tdim_;
            /// Total number of variables
            int nvar_;

            /// Check that ordering of vertices of mesh entity is valid.
            bool check_mesh ( );

            /// Holds DoF IDs, needed for mapl2g
            std::vector<DofID> numer_;

            /// Offset container for numer_, needed for mapl2g
            std::vector<std::vector<int> > numer_offsets_cell_varloc_;

            /// Const pointer to mesh
            const mesh::Mesh* mesh_;

            /// FEManager on the given mesh
            FEManager<DataType> const* fe_manager_;

            /// Interpolation Container, which stores the interpolation weights
            DofInterpolation dof_interpolation_;

            /// Get the geometrical dimension stored in FEManager
            int get_dim ( ) const;

            /// Update number_of_dofs_total_ and number_of_dofs_for_var_
            void update_number_of_dofs ( const std::string& = "" );

            /// Total number of dofs for all variables
            int number_of_dofs_total_;

            /// Total number of dofs per variable
            std::vector<int> number_of_dofs_for_var_;

          private:
            /// Check if mesh is set
            bool mesh_flag_;
            /// Check if fe_manager is set
            bool fe_manager_flag_;
            /// Check if numbering strategy is set
            bool strategy_;
        };

        template<class DataType>
        void permute_constrained_dofs_to_end ( DegreeOfFreedom<DataType>& dof );

    } // namespace doffem
} // namespace hiflow
#endif
