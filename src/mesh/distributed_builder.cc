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

#include "distributed_builder.h"

#include <algorithm>
#include <numeric>
#include <mpi.h>

#include "common/log.h"
#include "types.h"

const int DEBUG_LEVEL = 3;

namespace hiflow
{
    namespace mesh
    {
        typedef DistributedBuilder::VertexHandle VertexHandle;
        typedef DistributedBuilder::EntityHandle EntityHandle;

        DistributedBuilder::DistributedBuilder ( const MPI_Comm& communicator,
                                                 MeshBuilder& local_builder,
                                                 int sender,
                                                 const std::vector<int>* partition )
        : MeshBuilder ( local_builder.tdim ( ), local_builder.gdim ( ) ),
        communicator_ ( communicator ),
        local_builder_ ( local_builder ),
        sender_ ( sender ),
        partition_ ( partition ),
        num_sent_entities_ ( 0 )
        {
            // setup type for points
            MPI_Type_contiguous ( this->gdim ( ), MPI_DOUBLE, &point_type_ );
            MPI_Type_commit ( &point_type_ );
        }

        DistributedBuilder::~DistributedBuilder ( )
        {
            MPI_Type_free ( &point_type_ );
        }

        // NB coordinates is empty except for on sender

        VertexHandle DistributedBuilder::add_vertex ( const std::vector<Coordinate>& coordinates )
        {
            assert ( coordinates.empty ( ) || is_sender ( ) );
            if ( is_sender ( ) )
            {
                assert ( !coordinates.empty ( ) );
                assert ( static_cast < int > ( coordinates.size ( ) ) == gdim ( ) );
            }

            // broadcast vertex to all processes
            broadcast_vertices ( coordinates, 1 );

            // index of vertex in common coordinate structure
            return coordinates_.size ( ) % gdim ( );
        }

        std::vector<VertexHandle> DistributedBuilder::add_vertices ( const std::vector<Coordinate>& coordinates )
        {
            assert ( coordinates.empty ( ) || is_sender ( ) );
            if ( is_sender ( ) )
            {
                assert ( !coordinates.empty ( ) );
            }

            const EntityCount num_vertices_before = coordinates_.size ( ) / gdim ( );

            // broadcast number of vertices to receiving processes
            EntityCount num_new_vertices = coordinates.size ( ) / gdim ( );
            MPI_Bcast ( &num_new_vertices, 1, MPI_INT, sender_, communicator_ );

            // broadcast vertices to all processes
            broadcast_vertices ( coordinates, num_new_vertices );

            std::vector<VertexHandle> vertices ( num_new_vertices );
            for ( int i = 0; i < num_new_vertices; ++i )
            {
                vertices[i] = num_vertices_before + i;
            }

            return vertices;
        }

        EntityHandle DistributedBuilder::add_entity ( TDim tdim,
                                                      const std::vector<VertexHandle>& vertices )
        {
            assert ( vertices.empty ( ) || is_sender ( ) );
            if ( is_sender ( ) )
            {
                assert ( !vertices.empty ( ) );
            }

            // Broadcast dimension of entities
            MPI_Bcast ( &tdim, 1, MPI_INT, sender_, communicator_ );

            // Do point-to-point communication between sender and receiver
            int entity_size = vertices.size ( );

            // Find receiving process, and broadcast to all
            int receiver;
            if ( is_sender ( ) )
            {
                assert ( partition_ != 0 );
                assert ( num_sent_entities_ >= 0 );
                assert ( num_sent_entities_ < static_cast < int > ( partition_->size ( ) ) );

                receiver = ( *partition_ )[num_sent_entities_];
            }
            MPI_Bcast ( &receiver, 1, MPI_INT, sender_, communicator_ );

            const bool is_receiver = ( rank ( ) == receiver );

            // if sender_ == receiver, no communication is necessary
            std::vector<VertexHandle> received_vertices;
            if ( sender_ != receiver )
            {
                if ( is_sender ( ) )
                {
                    MPI_Send ( &tdim, 1, MPI_INT, receiver, 1, communicator_ );
                    MPI_Send ( &entity_size, 1, MPI_INT, receiver, 2, communicator_ );
                    MPI_Send ( const_cast < VertexHandle* > ( vec2ptr ( vertices ) ),
                               entity_size, MPI_INT, receiver, 3, communicator_ );
                }
                else if ( is_receiver )
                {
                    MPI_Recv ( &tdim, 1, MPI_INT, sender_, 1, communicator_, MPI_STATUS_IGNORE );
                    MPI_Recv ( &entity_size, 1, MPI_INT, sender_, 2, communicator_, MPI_STATUS_IGNORE );
                    received_vertices.resize ( entity_size );
                    MPI_Recv ( vec2ptr ( received_vertices ), entity_size, MPI_INT, sender_, 3, communicator_, MPI_STATUS_IGNORE );
                }
            }
            else
            {
                received_vertices = vertices;
            }

            if ( is_sender ( ) )
            {
                ++num_sent_entities_;
            }

            EntityHandle hdl = -1;
            if ( is_receiver )
            {
                std::vector<VertexHandle> entity ( entity_size, -1 );

                for ( int v = 0; v < entity_size; ++v )
                {
                    const int vertex_index = received_vertices[v];
                    entity[v] = add_local_vertex ( vertex_index );
                }

                hdl = local_builder_.add_entity ( tdim, entity );
            }

            return hdl;
        }

        std::vector<EntityHandle> DistributedBuilder::add_entities ( TDim tdim,
                                                                     const std::vector<VertexHandle>& vertices,
                                                                     const std::vector<int>& sizes )
        {
            assert ( vertices.empty ( ) || is_sender ( ) );
            assert ( sizes.empty ( ) || is_sender ( ) );
            if ( is_sender ( ) )
            {
                assert ( !vertices.empty ( ) );
                assert ( !sizes.empty ( ) );
            }

            // --- Setup

            // data to be sent (sender only)
            std::vector<int> num_entities_on_proc; // number of entities on each process (length = num processes)
            std::vector<int> permuted_entity_sizes; // sizes of entities in permuted order (length = num of entities)
            std::vector<int> num_vertices_on_proc; // number of vertices on each process (length = num processes)
            std::vector<VertexHandle> permuted_vertices; // vertices of entities in permuted order (length = num vertices)

            // extra data needed for sender
            std::vector<int> entity_offsets_on_proc; // offset into permuted_entity_sizes (length = num processes)
            std::vector<int> vertex_offsets_on_proc; // offset into permuted_vertices (length = num processes)

            if ( is_sender ( ) )
            {
                assert ( partition_ != 0 );
                assert ( num_sent_entities_ >= 0 );
                assert ( num_sent_entities_ < static_cast < int > ( partition_->size ( ) ) );
                assert ( num_sent_entities_ + sizes.size ( ) <= partition_->size ( ) );

                const int num_entities = sizes.size ( );

                compute_num_entities_on_proc ( num_entities, num_entities_on_proc );

                // Compute offsets into entity permutation for each proc
                convert_sizes_to_offsets ( num_entities_on_proc, entity_offsets_on_proc );

                // Compute entity permutation vector
                // This maps permuted entity indices to original entity indices
                std::vector<int> entity_permutation;
                compute_entity_permutation ( num_entities, entity_offsets_on_proc, entity_permutation );

                // Get entity sizes on each process by permuting the sizes vector
                permute_vector ( sizes, entity_permutation, permuted_entity_sizes );

                // Compute the number of vertices for each proc
                compute_num_vertices_on_proc ( permuted_entity_sizes, entity_offsets_on_proc, num_vertices_on_proc );

                // Compute offsets into permuted_vertices for each proc
                convert_sizes_to_offsets ( num_vertices_on_proc, vertex_offsets_on_proc );

                // Compute vertex offsets for permuted entities
                std::vector<int> permuted_entity_offsets;
                convert_sizes_to_offsets ( permuted_entity_sizes, permuted_entity_offsets );

                // Compute permuted vertex vector
                compute_permuted_vertices ( sizes, num_entities_on_proc, entity_permutation, vertices, permuted_vertices );
                // Update number of entities sent. Functions using this
                // member variable should not be called on sending process after this point.
                num_sent_entities_ += num_entities;
            }

            // --- Communication

            // Broadcast dimension of entities
            MPI_Bcast ( &tdim, 1, MPI_INT, sender_, communicator_ );

            // Scatter number of entities to receiving processes
            int num_received_entities;
            MPI_Scatter ( vec2ptr ( num_entities_on_proc ), 1, MPI_INT,
                          &num_received_entities, 1, MPI_INT, sender_, communicator_ );

            // Scatter entity sizes to receiving processes
            std::vector<int> received_entity_sizes ( num_received_entities );
            MPI_Scatterv ( vec2ptr ( permuted_entity_sizes ), vec2ptr ( num_entities_on_proc ),
                           vec2ptr ( entity_offsets_on_proc ), MPI_INT,
                           vec2ptr ( received_entity_sizes ), num_received_entities, MPI_INT, sender_, communicator_ );

            // Scatter vertices to receiving processes
            const int num_received_vertices = std::accumulate ( received_entity_sizes.begin ( ),
                                                                received_entity_sizes.end ( ), 0 );
            std::vector<int> received_vertices ( num_received_vertices );
            MPI_Scatterv ( vec2ptr ( permuted_vertices ), vec2ptr ( num_vertices_on_proc ),
                           vec2ptr ( vertex_offsets_on_proc ), MPI_INT,
                           vec2ptr ( received_vertices ), num_received_vertices, MPI_INT, sender_, communicator_ );

            // --- Compute offsets into received_vertices
            std::vector<int> received_entity_offsets;
            convert_sizes_to_offsets ( received_entity_sizes, received_entity_offsets );

            // --- Add vertices and entities locally
            std::vector<EntityHandle> entities ( num_received_entities, -1 );
            for ( int e = 0; e < num_received_entities; ++e )
            {
                const int entity_offset = received_entity_offsets[e];
                const int entity_size = received_entity_sizes[e];
                std::vector<VertexHandle> entity ( entity_size, -1 );

                for ( int v = 0; v < entity_size; ++v )
                {
                    const int vertex_index = received_vertices[entity_offset + v];
                    entity[v] = add_local_vertex ( vertex_index );
                }

                entities[e] = local_builder_.add_entity ( tdim, entity );
            }

            // NB: Only the locally created entities are returned. It
            // would be possible to send back all remote entity id:s to
            // the sender process if this is useful.
            return entities;
        }

        void DistributedBuilder::set_material_number ( TDim tdim, EntityHandle entity, MaterialNumber material )
        {
            // TODO: implement me!
        }

        void DistributedBuilder::clear ( )
        {
            coordinates_.clear ( );
            local_vertices_.clear ( );
            num_sent_entities_ = 0;
            local_builder_.clear ( );
        }

        MeshPtr DistributedBuilder::build ( )
        {
            return local_builder_.build ( );
        }

        bool DistributedBuilder::is_sender ( ) const
        {
            return rank ( ) == sender_;
        }

        int DistributedBuilder::rank ( ) const
        {
            int rank = -1;
            MPI_Comm_rank ( communicator_, &rank );
            return rank;
        }

        int DistributedBuilder::num_processes ( ) const
        {
            int size = -1;
            MPI_Comm_size ( communicator_, &size );
            return size;
        }

        // NB coordinates is empty except for on sender, but all processes
        // must now the number of vertices

        void DistributedBuilder::broadcast_vertices ( const std::vector<Coordinate>& coords,
                                                      EntityCount num_new_vertices )
        {
            assert ( num_new_vertices != 0 );
            if ( is_sender ( ) )
            {
                assert ( !coords.empty ( ) );
                assert ( static_cast < int > ( coords.size ( ) ) == num_new_vertices * gdim ( ) );
            }

            // allocate memory
            const int pos = coordinates_.size ( );
            if ( is_sender ( ) )
            {
                assert ( num_new_vertices * gdim ( ) == static_cast < int > ( coords.size ( ) ) );
                coordinates_.insert ( coordinates_.end ( ), coords.begin ( ), coords.end ( ) );
            }
            else
            {
                coordinates_.insert ( coordinates_.end ( ), num_new_vertices * gdim ( ), 0.0 );
            }

            // Broadcast from process sender_, putting data into local coordinates_ vector
            Coordinate* buf = &coordinates_.front ( ) + pos;
            MPI_Bcast ( buf, num_new_vertices, point_type_, sender_, communicator_ );

            local_vertices_.insert ( local_vertices_.end ( ), num_new_vertices, -1 );
        }

        void DistributedBuilder::compute_entity_permutation ( int num_entities,
                                                              const std::vector<int>& entity_offsets_for_proc,
                                                              std::vector<int>& entity_permutation ) const
        {
            assert ( partition_ != 0 );
            assert ( static_cast < int > ( entity_offsets_for_proc.size ( ) ) == num_processes ( ) );

            entity_permutation.clear ( );
            entity_permutation.resize ( num_entities );

            // store relative offsets within block belonging to each proc
            std::vector<int> relative_offsets ( num_processes ( ), 0 );

            // compute permuted position for each entity
            for ( int i = 0; i < num_entities; ++i )
            {
                const int ind = num_sent_entities_ + i;
                const int p = ( *partition_ )[ind];
                const int new_pos = entity_offsets_for_proc[p] + relative_offsets[p];
                entity_permutation[new_pos] = i;
                ++relative_offsets[p];
            }
        }

        // Computes the number of entities that will

        void DistributedBuilder::compute_num_entities_on_proc ( int num_entities,
                                                                std::vector<int>& num_entities_on_proc ) const
        {
            assert ( partition_ != 0 );
            num_entities_on_proc.clear ( );
            num_entities_on_proc.resize ( num_processes ( ) );

            for ( int i = 0; i < num_entities; ++i )
            {
                const int ind = num_sent_entities_ + i;
                const int proc = ( *partition_ )[ind];
                ++num_entities_on_proc[proc];
            }
        }

        void DistributedBuilder::compute_num_vertices_on_proc ( const std::vector<int>& permuted_entity_sizes,
                                                                const std::vector<int>& entity_offsets_on_proc,
                                                                std::vector<int>& num_vertices_on_proc ) const
        {
            const int num_processes = entity_offsets_on_proc.size ( );
            num_vertices_on_proc.clear ( );
            num_vertices_on_proc.resize ( num_processes );

            const int num_entities = permuted_entity_sizes.size ( );

            for ( int p = 0; p < num_processes; ++p )
            {
                const int begin_entity = entity_offsets_on_proc[p];
                const int end_entity = ( p + 1 == num_processes ) ?
                        num_entities : entity_offsets_on_proc[p + 1];

                // sum entity
                num_vertices_on_proc[p] =
                        std::accumulate ( permuted_entity_sizes.begin ( ) + begin_entity,
                                          permuted_entity_sizes.begin ( ) + end_entity, 0 );
            }
        }

        // Convert a vector of sizes of elements of a list to a vector of
        // offsets into the list.

        void DistributedBuilder::convert_sizes_to_offsets ( const std::vector<int>& sizes,
                                                            std::vector<int>& offsets ) const
        {
            offsets.clear ( );
            offsets.resize ( sizes.size ( ) );
            offsets[0] = 0;
            std::partial_sum ( sizes.begin ( ), sizes.end ( ) - 1, offsets.begin ( ) + 1 );
        }

        void DistributedBuilder::permute_vector ( const std::vector<int>& source,
                                                  const std::vector<int>& permutation,
                                                  std::vector<int>& target ) const
        {
            assert ( source.size ( ) == permutation.size ( ) );

            target.clear ( );
            target.resize ( source.size ( ) );

            for ( size_t i = 0; i != target.size ( ); ++i )
            {
                target[i] = source[permutation[i]];
            }
        }

        void DistributedBuilder::compute_permuted_vertices ( const std::vector<int>& entity_sizes,
                                                             const std::vector<int>& num_entities_on_proc,
                                                             const std::vector<int>& entity_permutation,
                                                             const std::vector<VertexHandle>& vertices,
                                                             std::vector<VertexHandle>& permuted_vertices ) const
        {
            permuted_vertices.clear ( );
            permuted_vertices.resize ( vertices.size ( ) );

            // Compute vertex offsets for each original entity
            std::vector<int> entity_offsets;
            convert_sizes_to_offsets ( entity_sizes, entity_offsets );

            // NB: sequential loop over indices in target (permuted) vector, copying
            // in the vertex from non-sequential blocks in the source (original) vector
            int target_offset = 0;
            int permuted_entity = 0;
            for ( int p = 0; p < num_processes ( ); ++p )
            {
                for ( int e = 0; e < num_entities_on_proc[p]; ++e )
                {
                    const int src_entity = entity_permutation[permuted_entity];
                    const int entity_size = entity_sizes[src_entity];
                    const int src_offset = entity_offsets[src_entity];
                    std::copy ( vertices.begin ( ) + src_offset,
                                vertices.begin ( ) + src_offset + entity_size,
                                permuted_vertices.begin ( ) + target_offset );
                    target_offset += entity_size;
                    ++permuted_entity;
                }
            }
        }

        // Adds the vertex with index vertex_index in the coordinates
        // vector, if it has not already been visited. Updates the
        // local_vertices_ structure with the returned VertexHandle, and
        // returns this handle to the caller.

        VertexHandle DistributedBuilder::add_local_vertex ( int vertex_index )
        {
            assert ( vertex_index >= 0 );
            assert ( vertex_index < static_cast < int > ( local_vertices_.size ( ) ) );
            typedef std::vector<Coordinate>::const_iterator CoordPtr;

            if ( local_vertices_[vertex_index] == -1 )
            {
                CoordPtr begin = coordinates_.begin ( ) + gdim ( ) * vertex_index;
                CoordPtr end = coordinates_.begin ( ) + gdim ( ) * ( vertex_index + 1 );
                VertexHandle v = local_builder_.add_vertex (
                                                             std::vector<Coordinate>( begin, end ) );
                local_vertices_[vertex_index] = v;
            }
            return local_vertices_[vertex_index];
        }

#if 0

        void permute_entity_vertices ( const std::vector<VertexHandle>& vertices,
                                       const std::vector<int>& sizes )
        {
            const EntityCount num_entities = sizes.size ( );
            std::vector<int> num_vertices_on_proc ( this->num_processes ( ), 0 );

            for ( int i = 0; i < num_entities; ++i )
            {
                const int ind = num_sent_entities_ + i;
                const int proc = partition_[ind];
                num_vertices_on_proc[proc] += sizes[i];
            }

            // compute offsets into the permuted vertex vector
            std::vector<int> process_offsets ( num_processes ( ) );
            process_offsets[0] = 0;
            std::partial_sum ( num_vertices_on_proc.begin ( ), num_vertices_on_proc.end ( ) - 1,
                               process_offsets.begin ( ) + 1 );

            // offset relative to process_offsets[proc] for position to store
            // next entity in permuted_vertices
            std::vector<int> relative_offset ( num_processes ( ), 0 );

            // copy data over into permuted vertex vector
            std::vector<VertexHandle> permuted_vertices ( vertices.size ( ) );
            int src_offset = 0;
            for ( int i = 0; i < num_entities; ++i )
            {
                const int ind = num_sent_entities_ + i;
                const int proc = partition_[ind];
                const int entity_size = sizes[i];
                const int target_offset = process_offsets[proc] + relative_offset[proc];
                std::copy ( vertices.begin ( ) + src_offset,
                            vertices.begin ( ) + src_offset + entity_size,
                            permuted_vertices.begin ( ) + target_offset );
                src_offset += entity_size;
                relative_offset[proc] += entity_size;
            }
        }
#endif

    }
} // namespace hiflow
