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

/// \author Thomas Gengenbach and Staffan Ronnas

#ifndef HIFLOW_MESH_REFINED_MESH_DB_VIEW_H
#    define HIFLOW_MESH_REFINED_MESH_DB_VIEW_H

#    include <vector>

#    include "mesh_db_view.h"
#    include "types.h"

namespace hiflow
{
    namespace mesh
    {

        /// \brief Class specializing some MeshDbView behavior after refinement.

        class RefinedMeshDbView : public MeshDbView
        {
          public:

            RefinedMeshDbView ( TDim tdim, GDim gdim, MeshDatabasePtr db,
                                const std::vector<Id>& cells,
                                const std::vector<ConstMeshPtr>& history,
                                const std::vector<int>& parent_history_indices,
                                const std::vector<EntityNumber>& parent_cell_indices,
                                const std::vector<int>& sub_cell_numbers,
                                std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            RefinedMeshDbView ( TDim tdim, GDim gdim, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Get parent cell id of a cell
            Id get_parent_cell_id ( EntityNumber cell_index ) const;

            Entity get_parent_cell ( EntityNumber cell_index ) const;

            bool cell_has_parent ( EntityNumber cell_index ) const;

            /// \brief Get an ancestor mesh
            ConstMeshPtr get_ancestor_mesh ( int history_index ) const;

            std::vector<ConstMeshPtr> get_all_ancestors ( ) const;

            void set_history ( std::vector<ConstMeshPtr> history );

            std::pair<ConstMeshPtr, EntityNumber> get_cell_parent ( EntityNumber cell_index ) const;

            using MeshDbView::refine;

            virtual MeshPtr refine ( std::vector<int>& refinements ) const;

            /// \brief Determines if facet is on the boundary
            virtual bool is_boundary_facet ( EntityNumber facet_index ) const;

            /// \brief copy mesh object from given reference
            // TODO

            virtual void copy_from ( const MeshPtr mesh )
            {
            }

          private:

            void mark_forbidden_ancestor_cells ( EntityNumber cell_index,
                                                 SortedArray<Id>& ancestors ) const;

            bool find_coarsenable_parent ( EntityNumber cell_index,
                                           const SortedArray<Id>& forbidden_ancestors,
                                           int* parent_history_index,
                                           EntityNumber* parent_cell_index ) const;

            AttributePtr parent_history_index_attr_;
            AttributePtr parent_cell_index_attr_;
            AttributePtr sub_cell_number_attr_;

            // STRUCT_TODO move to database, see MeshPXest
            std::vector<ConstMeshPtr> history_;
        };

    }
} // namespace hiflow

#endif /* _REFINED_MESH_DB_VIEW_H_ */
