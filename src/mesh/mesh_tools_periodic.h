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

#ifndef _MESH_TOOLS_PERIODIC_H_
#    define _MESH_TOOLS_PERIODIC_H_

#    include "mesh.h"
#    include "periodicity_tools.h"
#    include <string>

#    include <mpi.h>

/// @brief High-level tools for dealing with periodic meshes
/// @author Teresa Beck

class MasterSlave;

namespace hiflow
{
    namespace mesh
    {
        class GraphPartitioner;
        struct SharedVertexTable;

        /// \brief Create new mesh containing the ghost cells from neighboring sub-domains.
        MeshPtr compute_ghost_cells_periodic ( const Mesh & local_mesh, MPI_Comm& comm,
                                               SharedVertexTable& shared_verts );

    }
}

#endif /* _MESH_TOOLS_PERIODIC_H_ */
