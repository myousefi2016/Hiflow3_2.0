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

#include "mpi_communication.h"

#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>

#include "common/log.h"

const int DEBUG_LEVEL = 1;

namespace hiflow
{
    namespace mesh
    {

        //////////////// MpiBroadcast ////////////////

        MpiBroadcast::MpiBroadcast ( const MPI_Comm& mpi_communicator, int root )
        : root_ ( root ), rank_ ( -1 ), mpi_comm_ ( mpi_communicator )
        {
            MPI_Comm_rank ( mpi_comm_, &rank_ );
        }

        void MpiBroadcast::communicate ( /*const*/ EntityPackage* sent_entities,
                                         EntityPackage* received_entities ) const
        {
            assert ( received_entities != 0 );
            if ( is_root ( ) )
            {
                assert ( sent_entities != 0 );
            }

            // Pack integer data into one structure for communication.
            std::vector<int> integer_data ( 5, -1 );
            if ( is_root ( ) )
            {
                integer_data[0] = sent_entities->tdim;
                integer_data[1] = sent_entities->gdim;
                integer_data[2] = sent_entities->coords.size ( );
                integer_data[3] = sent_entities->offsets.size ( );
                integer_data[4] = sent_entities->connections.size ( );
            }

            MPI_Bcast ( vec2ptr ( integer_data ), 5, MPI_INT, root_, mpi_comm_ );

            // Prepare received_entities
            received_entities->tdim = integer_data[0];
            received_entities->gdim = integer_data[1];
            const int num_coords = received_entities->gdim * integer_data[2];
            received_entities->coords.resize ( num_coords );
            const int num_entities = integer_data[3] - 1;
            received_entities->offsets.resize ( num_entities + 1, -1 );
            const int num_connections = integer_data[4];
            received_entities->connections.resize ( num_connections, -1 );
            received_entities->material_numbers.resize ( num_entities, -10 );

            // Swap data into received_entities.

            // sent_entities is only temporarily modified, and will be
            // reset to its original state before end of function (assuming no exceptions occur).
            if ( is_root ( ) )
            {
                received_entities->coords = sent_entities->coords;
                received_entities->offsets = sent_entities->offsets;
                received_entities->connections = sent_entities->connections;
                received_entities->material_numbers = sent_entities->material_numbers;
            }

            // Broadcast vectors.
            MPI_Bcast ( vec2ptr ( received_entities->coords ), num_coords, MPI_DOUBLE, root_, mpi_comm_ );
            MPI_Bcast ( vec2ptr ( received_entities->offsets ), num_entities + 1, MPI_INT, root_, mpi_comm_ );
            MPI_Bcast ( vec2ptr ( received_entities->connections ), num_connections, MPI_INT, root_, mpi_comm_ );
            MPI_Bcast ( vec2ptr ( received_entities->material_numbers ), num_entities, MPI_INT, root_, mpi_comm_ );

            MPI_Barrier ( mpi_comm_ );
        }

        bool MpiBroadcast::is_root ( ) const
        {
            return rank ( ) == root_;
        }

        int MpiBroadcast::rank ( ) const
        {
            return rank_;
        }

        //////////////// MpiScatter ////////////////

        MpiScatter::MpiScatter ( const MPI_Comm& mpi_communicator, int sender,
                                 const std::vector<EntityCount>& num_entities_on_proc )
        : sender_ ( sender ), mpi_comm_ ( mpi_communicator ), num_entities_on_proc_ ( num_entities_on_proc )
        {
            if ( is_sender ( ) )
            {
                int size = -1;
                MPI_Comm_size ( mpi_comm_, &size );
                if ( static_cast < int > ( num_entities_on_proc.size ( ) ) != size )
                {
                    LOG_DEBUG ( 2, "num_entities on proc = " << num_entities_on_proc.size ( ) << ", num procs = " << size );
                }
                assert ( static_cast < int > ( num_entities_on_proc.size ( ) ) == size );
            }
            else
            {
                // should be empty on non-sending procs
                assert ( num_entities_on_proc.empty ( ) );
            }
        }

        void MpiScatter::communicate ( /*const*/EntityPackage* sent_entities,
                                       EntityPackage* received_entities ) const
        {
            LOG_DEBUG ( 2, "Communication starting on process " << rank ( ) );
            assert ( received_entities != 0 );
            if ( is_sender ( ) )
            {
                assert ( sent_entities != 0 );
            }

            // broadcast dimensions
            if ( is_sender ( ) )
            {
                received_entities->tdim = sent_entities->tdim;
            }
            LOG_DEBUG ( 2, "Broadcasting tdim" );
            MPI_Bcast ( &( received_entities->tdim ), 1, MPI_INT, sender_, mpi_comm_ );

            if ( is_sender ( ) )
            {
                received_entities->gdim = sent_entities->gdim;
            }
            LOG_DEBUG ( 2, "Broadcasting gdim" );
            MPI_Bcast ( &( received_entities->gdim ), 1, MPI_INT, sender_, mpi_comm_ );
            LOG_DEBUG ( 2, "GDim = " << received_entities->gdim << ", TDim = " << received_entities->tdim << " on proc " << rank ( ) );

            // broadcast coordinates vector size
            int num_coordinates;
            if ( is_sender ( ) )
            {
                num_coordinates = sent_entities->coords.size ( );
            }
            LOG_DEBUG ( 2, "Broadcasting num_coordinates" );
            MPI_Bcast ( &num_coordinates, 1, MPI_INT, sender_, mpi_comm_ );
            received_entities->coords.resize ( num_coordinates );
            LOG_DEBUG ( 3, "On proc " << rank ( ) << ", " << num_coordinates << " coordinates will be received" );

            // broadcast coordinates vector
            LOG_DEBUG ( 2, "Broadcasting coordinates" );
            if ( is_sender ( ) )
            {
                // TODO(Staffan): can we avoid copy here, while still keeping the code clean?
                received_entities->coords = sent_entities->coords;
            }
            MPI_Bcast ( vec2ptr ( received_entities->coords ), num_coordinates, MPI_DOUBLE, sender_, mpi_comm_ );

            // scatter number of received entities for each proc
            int num_received_entities;
            EntityCount* num_ent_on_proc = const_cast < EntityCount* > ( vec2ptr ( num_entities_on_proc_ ) ); // work-around C weakness
            LOG_DEBUG ( 2, "Scattering number of received entities" );
            MPI_Scatter ( num_ent_on_proc, 1, MPI_INT,
                          &num_received_entities, 1, MPI_INT, sender_, mpi_comm_ );
            received_entities->offsets.resize ( num_received_entities );

            // scatter the offset vector
            std::vector<int> entity_offset_on_proc;
            if ( is_sender ( ) )
            {
                convert_sizes_to_offsets ( num_entities_on_proc_, entity_offset_on_proc );
            }
            LOG_DEBUG ( 2, "Scattering offset vector" );
            MPI_Scatterv ( vec2ptr ( sent_entities->offsets ),
                           num_ent_on_proc,
                           vec2ptr ( entity_offset_on_proc ), MPI_INT,
                           vec2ptr ( received_entities->offsets ), num_received_entities, MPI_INT,
                           sender_, mpi_comm_ );

            // scatter the material_numbers vector
            received_entities->material_numbers.resize ( num_received_entities );
            MPI_Scatterv ( vec2ptr ( sent_entities->material_numbers ),
                           num_ent_on_proc,
                           vec2ptr ( entity_offset_on_proc ), MPI_INT,
                           vec2ptr ( received_entities->material_numbers ),
                           num_received_entities, MPI_INT, sender_, mpi_comm_ );

            // Compute number of connections, and scatter to each
            // process. Simultaneously, compute offset to be subtracted
            // from the received offsets on each proc
            std::vector<int> process_offsets;
            std::vector<EntityCount> num_connections_on_proc;

            if ( is_sender ( ) )
            {
                int num_processes = -1;
                MPI_Comm_size ( mpi_comm_, &num_processes );
                process_offsets.resize ( num_processes );
                num_connections_on_proc.resize ( num_processes );
                int curr_entity = 0;
                int total_connections = 0;

                for ( int p = 0; p < num_processes; ++p )
                {
                    process_offsets[p] = total_connections;
                    const EntityCount num_entities = num_entities_on_proc_[p];
                    const int next_entity = curr_entity + num_entities;
                    assert ( next_entity <= static_cast < int > ( sent_entities->offsets.size ( ) ) );
                    if ( next_entity == static_cast < int > ( sent_entities->offsets.size ( ) ) )
                    {
                        num_connections_on_proc[p] =
                                sent_entities->connections.size ( ) - sent_entities->offsets[curr_entity];
                    }
                    else
                    {
                        num_connections_on_proc[p] =
                                sent_entities->offsets[next_entity] - sent_entities->offsets[curr_entity];
                    }

                    LOG_DEBUG ( 3, "On process p = " << p << " num_entities = " << num_entities <<
                                ", next_entity = " << next_entity << " offsets.size() = " << sent_entities->offsets.size ( ) );

                    total_connections += num_connections_on_proc[p];
                    curr_entity = next_entity;
                }
            }

            int proc_offset;
            LOG_DEBUG ( 2, "Scattering process offsets" );
            MPI_Scatter ( vec2ptr ( process_offsets ), 1, MPI_INT,
                          &proc_offset, 1, MPI_INT, sender_, mpi_comm_ );

            int num_received_connections;
            LOG_DEBUG ( 2, "Scattering number of received connections" );
            MPI_Scatter ( vec2ptr ( num_connections_on_proc ), 1, MPI_INT,
                          &num_received_connections, 1, MPI_INT, sender_, mpi_comm_ );
            received_entities->connections.resize ( num_received_connections );
            LOG_DEBUG ( 3, "On process " << rank ( ) << ", " << num_received_connections << " connections will be received." );

            // subtract proc_offset from all received offsets
            LOG_DEBUG ( 3, "proc_offset = " << proc_offset << " on process " << rank ( ) );
            std::transform ( received_entities->offsets.begin ( ), received_entities->offsets.end ( ), received_entities->offsets.begin ( ),
                             std::bind2nd ( std::minus<int>( ), proc_offset ) );

            // set ending offset
            received_entities->offsets.push_back ( num_received_connections );

            LOG_DEBUG ( 3, "Received offset vector on process " << rank ( )
                        << " \n" << string_from_range ( received_entities->offsets.begin ( ), received_entities->offsets.end ( ) ) );

            // scatter connections vector
            std::vector<int> connection_offset_on_proc;
            if ( is_sender ( ) )
            {
                convert_sizes_to_offsets ( num_connections_on_proc, connection_offset_on_proc );
            }
            LOG_DEBUG ( 2, "Scattering connections" );
            MPI_Scatterv ( vec2ptr ( sent_entities->connections ), vec2ptr ( num_connections_on_proc ),
                           vec2ptr ( connection_offset_on_proc ), MPI_INT,
                           vec2ptr ( received_entities->connections ), num_received_connections, MPI_INT,
                           sender_, mpi_comm_ );

            LOG_DEBUG ( 2, "Communication done on process " << rank ( ) );
        }

        void MpiScatter::communicate ( std::vector<EntityPackage>& sent_entities,
                                       EntityPackage& received_entities ) const
        {
            LOG_DEBUG ( 2, "Communication starting on process " << rank ( ) );
            int num_procs = -1;
            MPI_Comm_size ( mpi_comm_, &num_procs );

            // broadcast dimensions
            if ( is_sender ( ) )
            {
                received_entities.tdim = sent_entities[0].tdim;
            }
            LOG_DEBUG ( 2, "Broadcasting tdim" );
            MPI_Bcast ( &( received_entities.tdim ), 1, MPI_INT, sender_, mpi_comm_ );

            if ( is_sender ( ) )
            {
                received_entities.gdim = sent_entities[0].gdim;
            }
            LOG_DEBUG ( 2, "Broadcasting gdim" );
            MPI_Bcast ( &( received_entities.gdim ), 1, MPI_INT, sender_, mpi_comm_ );
            LOG_DEBUG ( 2, "GDim = " << received_entities.gdim << ", TDim = " << received_entities.tdim << " on proc " << rank ( ) );

            // broadcast coordinates vector size
            std::vector<int> num_coordinates_sent ( num_procs, 0 );

            if ( is_sender ( ) )
            {
                for ( int i = 0; i < num_procs; ++i )
                {
                    num_coordinates_sent[i] = sent_entities[i].coords.size ( );
                }
            }

            int num_coordinates_recv = -1;
            LOG_DEBUG ( 2, "Scatter num_coordinates_sent" );
            MPI_Scatter ( vec2ptr ( num_coordinates_sent ),
                          1, MPI_INT,
                          &num_coordinates_recv,
                          1, MPI_INT, sender_, mpi_comm_ );

            received_entities.coords.resize ( num_coordinates_recv );
            LOG_DEBUG ( 3, "On proc " << rank ( ) << ", "
                        << num_coordinates_recv
                        << " coordinates will be received" );

            // broadcast coordinates vector
            LOG_DEBUG ( 2, "Scatter coordinates" );

            // Create data structures for sending
            std::vector<Coordinate> coordinates_sent;
            std::vector<int> sdispl ( num_procs + 1, 0 );

            if ( is_sender ( ) )
            {
                for ( int i = 0; i < num_procs; ++i )
                {
                    coordinates_sent.insert ( coordinates_sent.end ( ),
                                              sent_entities[i].coords.begin ( ),
                                              sent_entities[i].coords.end ( ) );
                    sdispl[i + 1] = sdispl[i] + sent_entities[i].coords.size ( );
                }
            }

            MPI_Scatterv ( vec2ptr ( coordinates_sent ),
                           vec2ptr ( num_coordinates_sent ),
                           vec2ptr ( sdispl ),
                           MPI_DOUBLE,
                           vec2ptr ( received_entities.coords ),
                           num_coordinates_recv,
                           MPI_DOUBLE,
                           sender_,
                           mpi_comm_ );

            num_coordinates_sent.clear ( );
            coordinates_sent.clear ( );
            sdispl.clear ( );

            // scatter number of received entities for each proc
            LOG_DEBUG ( 2, "Scatter offset sizes" );
            std::vector<int> offset_sizes_sent ( num_procs, 0 );
            int offset_sizes_recv = -1;

            if ( is_sender ( ) )
            {
                for ( int i = 0; i < num_procs; ++i )
                {
                    offset_sizes_sent[i] = sent_entities[i].offsets.size ( );
                }
            }

            MPI_Scatter ( vec2ptr ( offset_sizes_sent ), 1, MPI_INT,
                          &offset_sizes_recv, 1, MPI_INT,
                          sender_,
                          mpi_comm_ );

            received_entities.offsets.resize ( offset_sizes_recv );
            LOG_DEBUG ( 3, "On proc " << rank ( ) << ", "
                        << offset_sizes_recv
                        << " offsets will be received" );

            LOG_DEBUG ( 2, "Scatter offset vector" );

            sdispl.resize ( num_procs + 1, 0 );
            std::vector<int> sent_offsets;
            if ( is_sender ( ) )
            {
                for ( int i = 0; i < num_procs; ++i )
                {
                    sent_offsets.insert ( sent_offsets.end ( ),
                                          sent_entities[i].offsets.begin ( ),
                                          sent_entities[i].offsets.end ( ) );
                    sdispl[i + 1] = sdispl[i] + sent_entities[i].offsets.size ( );

                }
            }

            MPI_Scatterv ( vec2ptr ( sent_offsets ),
                           vec2ptr ( offset_sizes_sent ),
                           vec2ptr ( sdispl ), MPI_INT,
                           vec2ptr ( received_entities.offsets ),
                           offset_sizes_recv,
                           MPI_INT,
                           sender_,
                           mpi_comm_ );
            offset_sizes_sent.clear ( );
            sent_offsets.clear ( );

            // scatter number of received entities for each proc
            LOG_DEBUG ( 2, "Scatter material_number sizes" );
            std::vector<int> mat_sizes_sent ( num_procs, 0 );
            int mat_sizes_recv = -1;

            if ( is_sender ( ) )
            {
                for ( int i = 0; i < num_procs; ++i )
                {
                    mat_sizes_sent[i] = sent_entities[i].material_numbers.size ( );
                }
            }

            MPI_Scatter ( vec2ptr ( mat_sizes_sent ), 1, MPI_INT,
                          &mat_sizes_recv, 1, MPI_INT,
                          sender_,
                          mpi_comm_ );

            received_entities.material_numbers.resize ( mat_sizes_recv );
            LOG_DEBUG ( 3, "On proc " << rank ( ) << ", "
                        << mat_sizes_recv
                        << " material numbers will be received" );

            LOG_DEBUG ( 2, "Scatter material_number vector" );

            sdispl.clear ( );
            sdispl.resize ( num_procs + 1, 0 );
            std::vector<int> sent_mat;
            if ( is_sender ( ) )
            {
                for ( int i = 0; i < num_procs; ++i )
                {
                    sent_mat.insert ( sent_mat.end ( ),
                                      sent_entities[i].material_numbers.begin ( ),
                                      sent_entities[i].material_numbers.end ( ) );
                    sdispl[i + 1] = sdispl[i] + sent_entities[i].material_numbers.size ( );

                }
            }

            MPI_Scatterv ( vec2ptr ( sent_mat ),
                           vec2ptr ( mat_sizes_sent ),
                           vec2ptr ( sdispl ), MPI_INT,
                           vec2ptr ( received_entities.material_numbers ),
                           mat_sizes_recv,
                           MPI_INT,
                           sender_,
                           mpi_comm_ );
            mat_sizes_sent.clear ( );
            sent_mat.clear ( );

            // scatter number of received entities for each proc
            LOG_DEBUG ( 2, "Scatter connections sizes" );
            std::vector<int> conn_sizes_sent ( num_procs, 0 );
            int conn_sizes_recv = -1;

            if ( is_sender ( ) )
            {
                for ( int i = 0; i < num_procs; ++i )
                {
                    conn_sizes_sent[i] = sent_entities[i].connections.size ( );
                }
            }

            MPI_Scatter ( vec2ptr ( conn_sizes_sent ), 1, MPI_INT,
                          &conn_sizes_recv, 1, MPI_INT,
                          sender_,
                          mpi_comm_ );

            received_entities.connections.resize ( conn_sizes_recv );
            LOG_DEBUG ( 3, "On proc " << rank ( ) << ", "
                        << conn_sizes_recv
                        << " connections will be received" );

            LOG_DEBUG ( 2, "Scatter connections vector" );

            sdispl.clear ( );
            sdispl.resize ( num_procs + 1, 0 );
            std::vector<int> sent_conn;
            if ( is_sender ( ) )
            {
                for ( int i = 0; i < num_procs; ++i )
                {
                    sent_conn.insert ( sent_conn.end ( ),
                                       sent_entities[i].connections.begin ( ),
                                       sent_entities[i].connections.end ( ) );
                    sdispl[i + 1] = sdispl[i] + sent_entities[i].connections.size ( );

                }
            }

            MPI_Scatterv ( vec2ptr ( sent_conn ),
                           vec2ptr ( conn_sizes_sent ),
                           vec2ptr ( sdispl ), MPI_INT,
                           vec2ptr ( received_entities.connections ),
                           conn_sizes_recv,
                           MPI_INT,
                           sender_,
                           mpi_comm_ );
            conn_sizes_sent.clear ( );
            sent_conn.clear ( );

            LOG_DEBUG ( 2, "Communication done on process " << rank ( ) );
        }

        bool MpiScatter::is_sender ( ) const
        {
            return rank ( ) == sender_;
        }

        int MpiScatter::rank ( ) const
        {
            int rank = -1;
            MPI_Comm_rank ( mpi_comm_, &rank );
            return rank;
        }

        //////////////// MpiNonBlockingPointToPoint ////////////////

        MpiNonBlockingPointToPoint::MpiNonBlockingPointToPoint ( const MPI_Comm& communicator,
                                                                 int sender,
                                                                 int receiver )
        : sender_ ( sender ), receiver_ ( receiver ), comm_ ( communicator ), is_communicating_ ( false )
        {
        }

        void MpiNonBlockingPointToPoint::communicate ( /*const*/ EntityPackage* sent_entities,
                                                       EntityPackage* received_entities ) const
        {
            // tags
            const int LOCAL_BUF_TAG = 100;
            const int COORDS_TAG = 101;
            const int OFFSETS_TAG = 102;
            const int CONNECTIONS_TAG = 103;
            const int MATERIALS_TAG = 104;

            assert ( is_sender ( ) || is_receiver ( ) );
            assert ( !is_communicating_ );

            is_communicating_ = true;

            if ( is_sender ( ) )
            {
                assert ( sent_entities != 0 );

                // send entities

                // copy dimensions and sizes into local array
                // to send all at once
                local_buf_[0] = sent_entities->tdim;
                local_buf_[1] = sent_entities->gdim;
                local_buf_[2] = sent_entities->coords.size ( );
                local_buf_[3] = sent_entities->offsets.size ( );
                local_buf_[4] = sent_entities->connections.size ( );

                LOG_DEBUG ( 2, "sending local buf from " << sender_ << " to " << receiver_ );
                MPI_Isend ( &local_buf_[0], 5, MPI_INT, receiver_, LOCAL_BUF_TAG, comm_, &local_buf_request_ );

                // send coords
                LOG_DEBUG ( 2, "sending coords from " << sender_ << " to " << receiver_ );
                MPI_Isend ( vec2ptr ( sent_entities->coords ),
                            sent_entities->coords.size ( ),
                            MPI_DOUBLE, receiver_, COORDS_TAG, comm_, &coords_request_ );

                // send offsets
                LOG_DEBUG ( 2, "sending offsets from " << sender_ << " to " << receiver_ );
                MPI_Isend ( vec2ptr ( sent_entities->offsets ),
                            sent_entities->offsets.size ( ),
                            MPI_INT, receiver_, OFFSETS_TAG, comm_, &offsets_request_ );

                // send connections
                LOG_DEBUG ( 2, "sending connections from " << sender_ << " to " << receiver_ );
                MPI_Isend ( vec2ptr ( sent_entities->connections ),
                            sent_entities->connections.size ( ),
                            MPI_INT,
                            receiver_,
                            CONNECTIONS_TAG, comm_, &connections_request_ );

                // send materials
                LOG_DEBUG ( 2, "sending material numbers from " << sender_ << " to " << receiver_ );
                MPI_Isend ( vec2ptr ( sent_entities->material_numbers ),
                            sent_entities->material_numbers.size ( ),
                            MPI_INT,
                            receiver_,
                            MATERIALS_TAG, comm_, &materials_request_ );

            }
            else
            {
                assert ( received_entities != 0 );
                // receive entities

                // Receive first dimensions and sizes. This receive is
                // blocking, since it must be finished before other receives can be posted.
                LOG_DEBUG ( 2, "Receiving local buf from " << sender_ << " on " << receiver_ << " (blocking)" );
                MPI_Recv ( &local_buf_, 5, MPI_INT, sender_, LOCAL_BUF_TAG, comm_, &status_ );

                // set dimensions
                received_entities->tdim = local_buf_[0];
                received_entities->gdim = local_buf_[1];

                // resize vectors
                received_entities->coords.resize ( local_buf_[2] );
                received_entities->offsets.resize ( local_buf_[3] );
                received_entities->connections.resize ( local_buf_[4] );
                const int num_entities = local_buf_[3] - 1;
                received_entities->material_numbers.resize ( num_entities );

                // Receive four vectors
                LOG_DEBUG ( 2, "Receiving coords from " << sender_ << " on " << receiver_ );
                MPI_Irecv ( vec2ptr ( received_entities->coords ),
                            received_entities->coords.size ( ),
                            MPI_DOUBLE, sender_, COORDS_TAG, comm_, &coords_request_ );

                LOG_DEBUG ( 2, "Receiving offsets from " << sender_ << " on " << receiver_ );
                MPI_Irecv ( vec2ptr ( received_entities->offsets ),
                            received_entities->offsets.size ( ),
                            MPI_INT, sender_, OFFSETS_TAG, comm_, &offsets_request_ );

                LOG_DEBUG ( 2, "Receiving connections from " << sender_ << " on " << receiver_ );
                MPI_Irecv ( vec2ptr ( received_entities->connections ),
                            received_entities->connections.size ( ),
                            MPI_INT, sender_, CONNECTIONS_TAG, comm_, &connections_request_ );

                LOG_DEBUG ( 2, "Receiving materials from " << sender_ << " on " << receiver_ );
                MPI_Irecv ( vec2ptr ( received_entities->material_numbers ),
                            num_entities, MPI_INT, sender_, MATERIALS_TAG, comm_, &materials_request_ );
            }
        }

        void MpiNonBlockingPointToPoint::wait ( ) const
        {
            assert ( is_communicating_ );

            assert ( is_sender ( ) || is_receiver ( ) );

            if ( is_sender ( ) )
            {
                LOG_DEBUG ( 2, "Waiting on send on " << sender_ );
                MPI_Wait ( &local_buf_request_, &status_ );
                MPI_Wait ( &coords_request_, &status_ );
                MPI_Wait ( &offsets_request_, &status_ );
                MPI_Wait ( &connections_request_, &status_ );
                MPI_Wait ( &materials_request_, &status_ );
            }
            else
            {
                LOG_DEBUG ( 2, "Waiting on receive on " << receiver_ );
                MPI_Wait ( &coords_request_, &status_ );
                MPI_Wait ( &offsets_request_, &status_ );
                MPI_Wait ( &connections_request_, &status_ );
                MPI_Wait ( &materials_request_, &status_ );
            }
            LOG_DEBUG ( 2, "Wait done!" );

            is_communicating_ = false;
        }

        int MpiNonBlockingPointToPoint::rank ( ) const
        {
            int rank = -1;
            MPI_Comm_rank ( comm_, &rank );
            return rank;
        }

        bool MpiNonBlockingPointToPoint::is_sender ( ) const
        {
            return rank ( ) == sender_;
        }

        bool MpiNonBlockingPointToPoint::is_receiver ( ) const
        {
            return rank ( ) == receiver_;
        }

    }
} // namespace hiflow
