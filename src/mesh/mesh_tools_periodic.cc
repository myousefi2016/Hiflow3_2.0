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

/// @author Teresa Beck

#include "mesh_tools_periodic.h"

#include <mpi.h>

// public headers
#include "common/log.h"
#include "attributes.h"
#include "entity.h"
#include "iterator.h"
#include "mesh_db_view.h"
#include "partitioning.h"
#include "reader.h"

// private headers
#include "communication.h"
#include "mpi_communication.h"

const int DEBUG_LEVEL = 1;

namespace hiflow
{
    namespace mesh
    {
        /// \details The ghost cells are computed in two steps. First, a
        /// search to identify the vertices shared with other processes is
        /// performed, and the results are stored in the shared_verts
        /// object. Second, the cells containing shared vertices are
        /// exchanged between all neighboring processes (i.e. processes
        /// that share at least one vertex).
        ///
        /// The returned mesh object contains an integer attribute named
        /// "__sub_domain__" which contains the owner process number for
        /// all cells. Another integer attribute "__remote_index__"
        /// contains the index of the cell for all ghost cells, and -1 for
        /// all local cells.
        ///
        /// This function is limited to work only with meshes
        /// whose concrete type is PeriodicMesh.
        ///
        /// \param local_mesh            The mesh on the current process.
        /// \param comm                  MPI communicator over which communication is performed.
        /// \param[in,out] shared_verts  Table with information about shared vertices.
        /// \return  Shared pointer to mesh object containing ghost cells.

        MeshPtr compute_ghost_cells_periodic ( const Mesh & local_mesh, MPI_Comm& comm, SharedVertexTable& shared_verts )
        {
            int rank = -1;
            MPI_Comm_rank ( comm, &rank );

            const PeriodicMesh& mesh_db_view = dynamic_cast < const PeriodicMesh& > ( local_mesh );

            // exchange shared vertices with other processes
            update_shared_vertex_table ( mesh_db_view, comm, shared_verts );

            // objects needed for ghost cell communication
            const TDim tdim = local_mesh.tdim ( );
            std::vector<SubDomainId> sub_domains ( local_mesh.num_entities ( tdim ), rank );

            std::vector<EntityNumber> remote_indices ( local_mesh.num_entities ( tdim ), -1 );

            PeriodicMeshBuilder builder ( mesh_db_view );

            // communicate ghost cells
            communicate_ghost_cells ( mesh_db_view, comm, shared_verts, builder, sub_domains, remote_indices );
            communicate_ghost_cell_facets ( mesh_db_view, comm, shared_verts, builder );
            MeshPtr mesh_with_ghost_cells ( builder.build ( ) );

            // add attributes to new mesh
            assert ( mesh_with_ghost_cells != 0 );
            AttributePtr sub_domain_attr = AttributePtr ( new IntAttribute ( sub_domains ) );
            mesh_with_ghost_cells->add_attribute ( "_sub_domain_", tdim, sub_domain_attr );
            AttributePtr remote_index_attr = AttributePtr ( new IntAttribute ( remote_indices ) );
            mesh_with_ghost_cells->add_attribute ( "_remote_index_", tdim, remote_index_attr );

            return mesh_with_ghost_cells;
        }

    }
}
