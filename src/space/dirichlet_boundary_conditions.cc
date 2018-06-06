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

#include "dirichlet_boundary_conditions.h"

#include "mesh/entity.h"
#include "dof/degree_of_freedom.h"

/// @author Eva Ketelaer, Staffan Ronnas

namespace hiflow
{

    template<class DataType>
    void find_dofs_and_coords_on_face ( const mesh::EntityIterator& boundary_face,
                                        const VectorSpace<DataType>& space, int var,
                                        std::vector<doffem::DofID>& dof_ids,
                                        std::vector<std::vector<DataType> >& dof_points )
    {
        using namespace mesh;
        using namespace doffem;

        const Mesh& mesh = space.mesh ( );
        const TDim tdim = mesh.tdim ( );

        // get id of boundary face
        const Id boundary_id = boundary_face->id ( );

        // check if the boundary face exists and get the location where the entity
        // number should be stored
        int face_number;
        const bool check = mesh.find_entity ( tdim - 1, boundary_id, &face_number );
        assert ( check );

        // Get the face to be able to access to the data associated with the face
        Entity face = mesh.get_entity ( tdim - 1, face_number );

#ifndef NDEBUG
        // Get the cell associated with the cell and check if only one cell was
        // found
        IncidentEntityIterator dbg_cell = mesh.begin_incident ( face, tdim );
        assert ( dbg_cell != mesh.end_incident ( face, tdim ) );
        assert ( ++dbg_cell == mesh.end_incident ( face, tdim ) );
#endif

        // reset the cell because it was changed in the assert above
        IncidentEntityIterator cell = face.begin_incident ( tdim );

        // loop over all faces of the cell to get the local face index for
        // identifying the dofs
        int local_face_number = 0;
        for ( IncidentEntityIterator global_face = cell->begin_incident ( tdim - 1 );
              global_face != cell->end_incident ( tdim - 1 );
              ++global_face )
        {
            // if the global face id equals the boundary id the local face index is
            // found
            if ( global_face->id ( ) == boundary_id )
            {
                break;
            }
            else
            {
                local_face_number++;
            }
        }

        assert ( local_face_number >= 0 &&
                 local_face_number < cell->cell_type ( ).num_regular_entities ( tdim - 1 ) );

        // Get dofs and coords
        const DegreeOfFreedom<DataType>& dof = space.dof ( );

        dof.get_dofs_on_subentity ( var, cell->index ( ), tdim - 1, local_face_number,
                                    dof_ids );
        dof.get_coord_on_subentity ( var, cell->index ( ), tdim - 1, local_face_number,
                                     dof_points );
    }

    template void find_dofs_and_coords_on_face<double>( const mesh::EntityIterator& boundary_face,
            const VectorSpace<double>& space, int var,
            std::vector<doffem::DofID>& dof_ids,
            std::vector<std::vector<double> >& dof_points );
    template void find_dofs_and_coords_on_face<float>( const mesh::EntityIterator& boundary_face,
            const VectorSpace<float>& space, int var,
            std::vector<doffem::DofID>& dof_ids,
            std::vector<std::vector<float> >& dof_points );

} // namespace hiflow
