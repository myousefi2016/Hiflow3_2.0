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

#ifndef HIFLOW_DIRICHLET_BOUNDARY_CONDITIONS_H
#    define HIFLOW_DIRICHLET_BOUNDARY_CONDITIONS_H

#    include <mpi.h>

#    include <string>
#    include <vector>
#    include <cassert>

#    include "dof/degree_of_freedom.h"
#    include "mesh/iterator.h"
#    include "mesh/mesh.h"
#    include "space/vector_space.h"
#    include "mesh/types.h"
#    include "mesh/writer.h"

/// @author Eva Ketelaer, Staffan Ronnas

namespace hiflow
{

    /// \brief Finds the dofs and the coords on a given (boundary) facet.

    /// \details Serves as a helper function for the function
    /// compute_dirichlet_dofs_and_values.

    /// \param [in] boundary_face     iterator on the boundary facets.
    /// \param [in] space             space object containing the mesh and the FE
    ///                               approximation.
    /// \param [in] var               variable for which the dirichlet dofs should
    ///                               be set. The function can be called several
    ///                               times for vector-valued problems.
    /// \param [out] dof_ids          the id:s of the dofs on the facet.
    /// \param [out] dof_points       the coordinates of the dof id:s.
    template<class DataType>
    void find_dofs_and_coords_on_face ( const mesh::EntityIterator& boundary_face,
                                        const VectorSpace<DataType>& space,
                                        int var,
                                        std::vector<doffem::DofID>& dof_ids,
                                        std::vector<std::vector<DataType> >& dof_points );

    /// \brief Locate and evaluate Dirichlet boundary conditions.

    /// \details The function loops over all boundary facets and calls a
    /// user-provided function object (a functor) to obtain the values for
    /// the Dirichlet dofs on that boundary facet. This functor has to
    /// provide at least one function called evaluate, with the facet
    /// entity and a vector with the coordinates on that facet, that
    /// returns a std::vector<double> with the Dirichlet values for the
    /// dofs on this boundary facet. If the dofs of the boundary should
    /// not be constrained, an empty vector should be returned.
    /// \param [in] dirichlet_eval    user-provided function object which computes
    ///                               the dirichlet values.
    /// \param [in] space             space object containing the mesh and the FE
    ///                               approximation.
    /// \param [in] var               variable for which the dirichlet dofs should
    ///                               be set. The function can be called several
    ///                               times for vector-valued problems.
    /// \param [out] dirichlet_dofs   vector to which the indices of the dirichlet
    ///                               dofs are appended.
    /// \param [out] dirichlet_values vector which the values of the dirichlet dofs
    ///                               are appended.

    template<class DirichletEvaluator, class DataType>
    void compute_dirichlet_dofs_and_values ( DirichletEvaluator& dirichlet_eval,
                                             const VectorSpace<DataType>& space, int var,
                                             std::vector<doffem::DofID>& dirichlet_dofs,
                                             std::vector<DataType>& dirichlet_values )
    {
        using namespace mesh;

        typedef std::vector<DataType> Coord;

        const Mesh& mesh = space.mesh ( );
        const TDim tdim = mesh.tdim ( );

        // extract boundary of mesh
        MeshPtr boundary_mesh = mesh.extract_boundary_mesh ( );

        int rank = -1, size = -1;
        MPI_Comm_rank ( space.get_mpi_comm ( ), &rank );
        MPI_Comm_size ( space.get_mpi_comm ( ), &size );

        const bool is_sequential = ( size == 1 );
        if ( !is_sequential )
        {
            assert ( mesh.has_attribute ( "_sub_domain_", tdim ) );
        }

        // Loop over all faces which belong to the boundary
        for ( EntityIterator it_boundary = boundary_mesh->begin ( tdim - 1 );
              it_boundary != boundary_mesh->end ( tdim - 1 );
              ++it_boundary )
        {
            // get id of boundary face
            const Id boundary_id = it_boundary->id ( );

            // check if the boundary face exists and get the location
            // where the entity number should be stored
            int face_number;
            const bool check = mesh.find_entity ( tdim - 1, boundary_id, &face_number );
            assert ( check );

            // Get the face to be able to access to the data associated with the face
            Entity face = mesh.get_entity ( tdim - 1, face_number );

            // reset the cell because it was changed in the assert above
            IncidentEntityIterator cell = face.begin_incident ( tdim );

            std::vector<doffem::DofID> dofs_on_face;
            std::vector<Coord> coords_on_face;

            // Get the dofs indices and coords on the face
            find_dofs_and_coords_on_face ( it_boundary, space, var, dofs_on_face,
                                           coords_on_face );

            assert ( dofs_on_face.size ( ) == coords_on_face.size ( ) );

            // Evaluate user's function
            const std::vector<DataType> values_on_face =
                    dirichlet_eval.evaluate ( *it_boundary, coords_on_face );

            if ( !values_on_face.empty ( ) )
            {
                // If non-empty vector was returned, insert into output vectors
                // vectors to filter out only dofs that belong to our subdomain
                std::vector<doffem::DofID> dofs_on_face_checked;
                std::vector<DataType> values_on_face_checked;

                dofs_on_face_checked.reserve ( dofs_on_face.size ( ) );
                values_on_face_checked.reserve ( dofs_on_face.size ( ) );

                int k = 0;
                for ( std::vector<doffem::DofID>::iterator dof_it = dofs_on_face.begin ( );
                      dof_it != dofs_on_face.end ( ); ++dof_it )
                {
                    if ( space.dof ( ).owner_of_dof ( *dof_it ) == rank )
                    {
                        dofs_on_face_checked.push_back ( *dof_it );
                        values_on_face_checked.push_back ( values_on_face.at ( k ) );
                    }
                    ++k;
                }

                dirichlet_dofs.insert ( dirichlet_dofs.end ( ), dofs_on_face_checked.begin ( ),
                                        dofs_on_face_checked.end ( ) );
                dirichlet_values.insert ( dirichlet_values.end ( ),
                                          values_on_face_checked.begin ( ),
                                          values_on_face_checked.end ( ) );
            }
        }
    }

} // namespace hiflow

#endif
