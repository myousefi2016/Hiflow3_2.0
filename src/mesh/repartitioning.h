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

#ifndef _MESH_REPARTITIONING_H_
#    define _MESH_REPARTITIONING_H_

#    include "mesh/mesh.h"

#    include <string>
#    include <sstream>

#    include <mpi.h>
#    include <vector>

#    include "mesh/boundary_domain_descriptor.h"
#    include "mesh/partitioning.h"
#    include "mesh/mesh_tools.h"
#    include "mesh/periodicity_tools.h"
#    include "mesh/communication.h"

/// @author Simon Gawlok
namespace hiflow
{
    namespace mesh
    {
        class GraphPartitioner;
        struct SharedVertexTable;

#    ifdef WITH_PARMETIS
        /// \brief Scattering communicator

        class MpiAllToAll : public Communicator
        {
          public:
            MpiAllToAll ( const MPI_Comm& mpi_communicator );

            void communicate ( EntityPackage* sent_entities,
                               EntityPackage* received_entities ) const
            {
                std::vector<EntityPackage> recv_ent;
                recv_ent.push_back ( *received_entities );
            }

            void communicate ( std::vector<EntityPackage>& sent_entities,
                               std::vector<EntityPackage>& received_entities ) const;

          private:
            int rank ( ) const;

            MPI_Comm mpi_comm_;
        };

        /// \brief Compute dual graph of an already distributed mesh in parallel
        void compute_dual_graph_distributed ( const mesh::MeshPtr mesh,
                                              MPI_Comm& comm,
                                              mesh::Graph* graph,
                                              const std::vector<int>& vtxdist );

        /// \brief Repartition an already distributed mesh in parallel
        mesh::MeshPtr repartition_mesh ( const mesh::MeshPtr mesh,
                                         MPI_Comm& comm,
                                         const mesh::ParMetisGraphPartitioner* gpartitioner );

#    endif
    }
}

#endif /*MESH_REPARTITIONING_H_*/
