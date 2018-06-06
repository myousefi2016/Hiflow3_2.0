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

#ifndef HIFLOW_DISTRIBUTED_BUILDER_H
#    define HIFLOW_DISTRIBUTED_BUILDER_H

#    include <vector>

#    include "mesh/mesh_builder.h"

#    include <mpi.h>

namespace hiflow
{
    namespace mesh
    {

        class DistributedBuilder : public MeshBuilder
        {
          public:
            typedef MeshBuilder::VertexHandle VertexHandle;
            typedef MeshBuilder::EntityHandle EntityHandle;

            DistributedBuilder ( const MPI_Comm& communicator,
                                 MeshBuilder& local_builder,
                                 int sender, const std::vector<int>* distribution = 0 );
            ~DistributedBuilder ( );

            // add vertices
            virtual VertexHandle add_vertex ( const std::vector<Coordinate>& coordinates );
            virtual std::vector<VertexHandle> add_vertices ( const std::vector<Coordinate>& coordinates );

            // add entities
            virtual EntityHandle add_entity ( TDim tdim, const std::vector<VertexHandle>& vertices );
            virtual std::vector<EntityHandle> add_entities ( TDim tdim,
                                                             const std::vector<VertexHandle>& vertices,
                                                             const std::vector<int>& sizes );

            virtual void set_material_number ( TDim tdim, EntityHandle entity, MaterialNumber material );

            virtual void clear ( );

            // create mesh
            virtual MeshPtr build ( );

          private:
            // Member functions called on receiving processes
            bool is_sender ( ) const;
            int rank ( ) const;
            int num_processes ( ) const;
            VertexHandle add_local_vertex ( int vertex_index );

            // Member functions called on sending process only
            void broadcast_vertices ( const std::vector<Coordinate>& coords, EntityCount num_new_vertices );

            void compute_entity_permutation ( int num_entities,
                                              const std::vector<int>& entity_offsets_for_proc,
                                              std::vector<int>& entity_permutation ) const;

            void compute_num_entities_on_proc ( int num_entities,
                                                std::vector<int>& num_entities_on_proc ) const;

            void compute_num_vertices_on_proc ( const std::vector<int>& permuted_entity_sizes,
                                                const std::vector<int>& entity_offsets_on_proc,
                                                std::vector<int>& num_vertices_on_proc ) const;

            void convert_sizes_to_offsets ( const std::vector<int>& sizes,
                                            std::vector<int>& offsets ) const;

            void permute_vector ( const std::vector<int>& source,
                                  const std::vector<int>& permutation,
                                  std::vector<int>& target ) const;

            void compute_permuted_vertices ( const std::vector<int>& entity_sizes,
                                             const std::vector<int>& num_entities_on_proc,
                                             const std::vector<int>& entity_permutation,
                                             const std::vector<VertexHandle>& vertices,
                                             std::vector<VertexHandle>& permuted_vertices ) const;

            // Members relevant to receiving processes (includes sending process)
            MPI_Comm communicator_;
            MeshBuilder& local_builder_;
            const int sender_;
            MPI_Datatype point_type_;
            std::vector<Coordinate> coordinates_;
            std::vector<VertexHandle> local_vertices_;

            // Members relevant to sending process only

            // destination process for each added entity in the order that they are added
            const std::vector<int>* partition_;

            // number of entities that have already been sent
            VertexHandle num_sent_entities_;
        };
    }
} // namespace hiflow
#endif /* _DISTRIBUTED_BUILDER_H_ */
