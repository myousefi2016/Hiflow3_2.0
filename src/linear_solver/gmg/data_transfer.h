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

/// \author Aksel Alpay, Martin Wlotzka

#ifndef DATATRANSFER_H
#    define DATATRANSFER_H

#    include "tools/mpi_tools.h"
#    include "gmg_level.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            /// Holds information about a fine and a coarse grid and their processes,
            /// such as number of processes on both grids, and process lists.

            template<class LAD>
            class DataTransferInformation
            {
              public:
                /// Initializes the object. Collective on the communicator of
                /// the finer grid.
                /// @param coarse_level A pointer to a coarser grid of the MultiLevelHierarchy
                /// @param fine_level A pointer to a finer grid of the MultiLevelHierarchy
                /// While the two grids do not neccessarily have to stem from the
                /// same MultiLevelHierarchy object, the processes of the coarse grid
                /// communicator have to be subset of the processes in the communicator
                /// of the fine grid.

                DataTransferInformation ( const BasicLevel<LAD>* coarse_level,
                                          const BasicLevel<LAD>* fine_level );

                /// @return How many processes exist in the communicator of the coarse grid.

                int num_coarse_procs ( void ) const
                {
                    return num_coarse_procs_;
                }

                /// @return How many processes exist in the communicator of the fine grid.

                int num_fine_procs ( void ) const
                {
                    return num_fine_procs_;
                }

                /// @return The rank of the calling process in the communicator of the
                /// fine grid

                int fine_rank ( void ) const
                {
                    return fine_rank_;
                }

                /// @return The rank of the calling process in the communicator of the
                /// coarse grid

                int coarse_rank ( void ) const
                {
                    return coarse_rank_;
                }

                /// Returns a vector of the ranks in the fine communicator of the coarse processes
                /// (only working if the calling process is at least in the fine grid)
                /// @return the fine ranks of all coarse processes

                const std::vector<int>& coarse_procs_list ( void ) const
                {
                    return coarse_procs_list_;
                }

                /// Returns a vector of the ranks in the fine communicator of the fine processes
                /// (only working if the calling process is at least in the fine grid)
                /// @return the fine ranks of all fine processes

                const std::vector<int>& fine_procs_list ( void ) const
                {
                    return fine_procs_list_;
                }

                /// Translates the fine rank of a process to its coarse rank (if existing)
                /// @param fine_rank The fine rank to be translated
                /// @return The coarse rank corresponding to the fine rank, or -1 if
                /// the supplied rank does not belong to a coarse process

                int fine_to_coarse_rank ( const int fine_rank ) const
                {
                    assert ( fine_rank < coarse_procs_list_.size ( ) );
                    return coarse_procs_ranks_[fine_rank];
                }

                /// @return A pointer to the fine grid that has been used to
                /// construct the DataTransferInformation object

                const BasicLevel<LAD>* fine_level ( void ) const
                {
                    return fine_level_;
                }

                /// @return A pointer to the coarse grid that has been used to
                /// construct the DataTransferInformation object

                const BasicLevel<LAD>* coarse_level ( void ) const
                {
                    return coarse_level_;
                }

                /// @return Whether the calling process is in the communicator of the fine
                /// grid

                bool process_is_in_fine_level ( void ) const
                {
                    return fine_level_->is_scheduled_to_this_process ( );
                }

                /// @return Whether the calling process is in the communicator of the coarse
                /// grid

                bool process_is_in_coarse_level ( void ) const
                {
                    return coarse_level_->is_scheduled_to_this_process ( );
                }

              private:
                /// Fetches the ranks of all coarse processes, build the coarse process list
                /// stored in coarse_procs_list_ and counts the number of coarse processes
                /// which will be stored in num_coarse_procs_
                /// Collective on the communicator of the fine grid.

                void fetch_coarse_process_information ( void );

                int num_coarse_procs_;
                int num_fine_procs_;

                int fine_rank_;
                int coarse_rank_;

                std::vector<int> coarse_procs_ranks_;
                std::vector<int> coarse_procs_list_;

                std::vector<int> fine_procs_list_;

                const BasicLevel<LAD>* fine_level_;
                const BasicLevel<LAD>* coarse_level_;
            };

            /// This class implements data transfers via an array/map which specifies
            /// to which processes there is data to send.

            class AsymmetricDataTransfer
            {
              public:
                /// Initializes the object.
                /// @param comm The communicator to use for the transfer
                /// @param source_process_list A list of process ranks in the supplied
                /// communicator that will send data
                /// @param dest_process_list A list of process ranks in the supplied
                /// communicator that will receive data

                AsymmetricDataTransfer ( MPI_Comm comm,
                                         const std::vector<int>& source_process_list,
                                         const std::vector<int>& dest_process_list );

                virtual ~AsymmetricDataTransfer ( )
                {
                }

                /// transfers vectors of data from source to target processes.
                /// Collective on the union of the receiving and sending process groups.
                /// @param in The array of the input data. in[i] will be sent to the
                /// process with rank i (in the fine communicator). The sizes of the
                /// vectors stored in this argument may differ from each other.
                /// @param out The output array where the received values will be stored
                /// The i-th entry will be the data received from process i

                template<typename T>
                void transfer_vectors ( const std::vector<std::vector<T> >& in,
                                        std::vector<std::vector<T> >& out )
                {
                    if ( process_is_participating_ )
                    {
                        std::vector<int> input_sizes ( in.size ( ), 0 );
                        for ( unsigned i = 0; i < in.size ( ); ++i )
                            input_sizes[i] = static_cast < int > ( in[i].size ( ) );

                        transfer_vectors<T>( input_sizes, in, out );
                    }
                }

                /// Transfers vectors of data from source to target processes.
                /// Collective on the union of the receiving and sending process groups.
                /// @param in The array of the input data. in[i] will be sent to the
                /// process with rank i (in the fine communicator). The sizes of the
                /// vectors stored in this argument may differ from each other.
                /// @param out The output array where the received values will be stored
                /// The i-th entry will be the data received from process i
                /// @param input_sizes A vector whose i-th entry is the size of in[i].
                /// If this argument is not supplied, the transfer_vectors function
                /// calculates the sizes itself. If the the size data is however
                /// already available, it can be a (small) performance advantage
                /// to pass the existing data as argument to this function.

                template<typename T>
                void transfer_vectors ( const std::vector<int>& input_sizes,
                                        const std::vector<std::vector<T> >& in,
                                        std::vector<std::vector<T> >& out )
                {
                    if ( process_is_participating_ )
                    {
                        transfer_data<int>( input_sizes, transfer_sizes_ );
                        sizes_are_known_ = true;

                        perform_vector_transfer<T>( in, out );
                    }
                }

                /// Transfers vectors of data from source to target processes.
                /// If any of the transfer_vectors methods has been called before,
                /// this function will reuse the transfer sizes of that previous
                /// transfer. This avoids the need to send the sizes of the transfer
                /// to each of the receiving processes, which may be a performance
                /// gain.
                /// Collective on the union of the receiving and sending process groups.
                /// @param in The array of the input data. in[i] will be sent to the
                /// process with rank i (in the fine communicator). The sizes of the
                /// vectors stored in this argument may differ from each other.
                /// @param out The output array where the received values will be stored
                /// The i-th entry will be the data received from process i

                template<typename T>
                void transfer_vectors_by_reusing_sizes ( const std::vector<std::vector<T> >& in,
                                                         std::vector<std::vector<T> >& out )
                {
                    if ( sizes_are_known_ )
                        perform_vector_transfer<T>( in, out );
                    else
                        //unlike perform_vector_transfer, transfer_vectors will distribute the
                        //expected transfer sizes
                        transfer_vectors<T>( in, out );
                }

                /// Transfers vectors of data from source to target processes.
                /// If any of the transfer_vectors methods has been called before,
                /// this function will reuse the transfer sizes of that previous
                /// transfer. This avoids the need to send the sizes of the transfer
                /// to each of the receiving processes, which may be a performance
                /// gain.
                /// Collective on the union of the receiving and sending process groups.
                /// @param in The array of the input data. in[i] will be sent to the
                /// process with rank i (in the fine communicator). The sizes of the
                /// vectors stored in this argument may differ from each other.
                /// @param out The output array where the received values will be stored
                /// The i-th entry will be the data received from process i.
                /// @param input_sizes A vector whose i-th entry is the size of in[i].
                /// If this argument is not supplied, the transfer_vectors function
                /// calculates the sizes itself, if the sizes have not been sent before.
                /// If the the size data is however already available, it can be a (small)
                /// performance advantage to pass the existing data as argument to this function.

                template<typename T>
                void transfer_vectors_by_reusing_sizes ( const std::vector<int>& input_sizes,
                                                         const std::vector<std::vector<T> >& in,
                                                         std::vector<std::vector<T> >& out )
                {
                    if ( sizes_are_known_ )
                        perform_vector_transfer<T>( in, out );
                    else
                        transfer_vectors<T>( input_sizes, in, out );
                }

                /// Transfers plain data from source to target processes
                /// Collective on the union of the receiving and sending process groups.
                /// @param in The array of the input data. in[i] will be sent to the
                /// process with rank i (in the fine communicator)
                /// @param out The output array where the received values will be stored.
                /// The i-th entry will be the data received from process i.

                template<typename T>
                void transfer_data ( const std::vector<T>& in,
                                     std::vector<T>& out ) const
                {
                    if ( process_is_participating_ )
                    {
                        std::vector<MPI_Request> send_requests (
                                                                 receiving_processes_.size ( ),
                                                                 MPI_Request ( ) );

                        if ( process_is_sender_ )
                        {
                            assert ( in.size ( ) == receiving_processes_.size ( ) );
                            //send
                            for ( unsigned i = 0; i < receiving_processes_.size ( ); ++i )
                            {
                                int dest_process = receiving_processes_[i];

                                //           MPI_Isend(reinterpret_cast<void*> (const_cast<T*> (&(in[dest_process]))),
                                //                     1, mpi_data_type<T>::get_type(),
                                //                     dest_process, 0, comm_,
                                //                     &(send_requests[i]));
                                MPI_Isend ( reinterpret_cast < void* > ( const_cast < T* > ( &( in[i] ) ) ),
                                            1, mpi_data_type<T>::get_type ( ),
                                            dest_process, 0, comm_,
                                            &( send_requests[i] ) );
                            }
                        }

                        //receive
                        if ( process_is_receiver_ )
                        {
                            //received data will be stored here
                            out = std::vector<T>( sending_processes_.size ( ) );

                            std::vector<MPI_Request> recv_requests (
                                                                     sending_processes_.size ( ),
                                                                     MPI_Request ( ) );

                            for ( unsigned i = 0; i < sending_processes_.size ( ); ++i )
                            {
                                int source_process = sending_processes_[i];

                                MPI_Irecv ( &( out[i] ), 1, mpi_data_type<T>::get_type ( ),
                                            source_process, 0, comm_,
                                            &( recv_requests[i] ) );
                            }

                            std::vector<MPI_Status> recv_statuses (
                                                                    recv_requests.size ( ),
                                                                    MPI_Status ( ) );

                            MPI_Waitall ( recv_requests.size ( ),
                                          util::raw_array ( recv_requests ),
                                          util::raw_array ( recv_statuses ) );
                        }

                        if ( process_is_sender_ )
                        {
                            std::vector<MPI_Status> send_statuses (
                                                                    send_requests.size ( ),
                                                                    MPI_Status ( ) );

                            MPI_Waitall ( send_requests.size ( ),
                                          util::raw_array ( send_requests ),
                                          util::raw_array ( send_statuses ) );
                        }
                    }
                    if ( !process_is_receiver_ )
                    {
                        //For non-receiving processes, return empty array
                        out = std::vector<T>( );
                    }
                }

              private:
                /// The internal function used for the vector transfer.
                /// Requires that the transfer sizes have already been sent to the
                /// receiving processes and are stored in transfer_sizes_.
                /// Collective on the union of the receiving and sending process groups.
                /// @param in The array of the input data. in[i] will be sent to the
                /// process with rank i (in the fine communicator). The sizes of the
                /// vectors stored in this argument may differ from each other.
                /// @param out The output array where the received values will be stored
                /// The i-th entry will be the data received from process i

                template<typename T>
                void perform_vector_transfer ( const std::vector<std::vector<T> >& in,
                                               std::vector<std::vector<T> >& out ) const
                {
                    if ( process_is_participating_ )
                    {
                        //       assert(in.size() >= receiving_processes_.size());

                        std::vector<MPI_Request> send_requests;
                        send_requests.reserve ( receiving_processes_.size ( ) );

                        unsigned snd_ctr = 0;
                        if ( process_is_sender_ )
                        {
                            assert ( in.size ( ) == receiving_processes_.size ( ) );
                            //send
                            for ( unsigned i = 0; i < receiving_processes_.size ( ); ++i )
                            {
                                int dest_process = receiving_processes_[i];
                                if ( in[i].size ( ) != 0 )
                                    //           if (in[i].size() != 0)
                                {
                                    ++snd_ctr;

                                    send_requests.push_back ( MPI_Request ( ) );
                                    /*
                                                MPI_Isend(reinterpret_cast<void*> (const_cast<T*> (util::raw_array(in[dest_process]))),
                                                          in[dest_process].size(),
                                                          mpi_data_type<T>::get_type(),
                                                          dest_process, 0, comm_,
                                                          &(send_requests.back()));*/
                                    MPI_Isend ( reinterpret_cast < void* > ( const_cast < T* > ( util::raw_array ( in[i] ) ) ),
                                                in[i].size ( ),
                                                mpi_data_type<T>::get_type ( ),
                                                dest_process, 0, comm_,
                                                &( send_requests.back ( ) ) );
                                }
                            }
                        }
                        assert ( snd_ctr == send_requests.size ( ) );

                        if ( process_is_receiver_ )
                        {
                            //receive
                            //         out = std::vector<std::vector<T> >(
                            //                 comm_size_,
                            //                 std::vector<T>());
                            out = std::vector<std::vector<T> >(
                                    sending_processes_.size ( ),
                                    std::vector<T>( ) );

                            std::vector<MPI_Request> recv_requests;

                            recv_requests.reserve ( sending_processes_.size ( ) );

                            //         for (unsigned source_process = 0;
                            //                 source_process < sending_processes_.size();
                            //                 ++source_process)
                            //         {
                            //           out[source_process] =
                            //                   std::vector<T>(transfer_sizes_[source_process], 0);
                            //
                            //           if (transfer_sizes_[source_process] != 0)
                            //           {
                            //             recv_requests.push_back(MPI_Request());
                            //
                            //             MPI_Irecv(util::raw_array(out[source_process]),
                            //                       transfer_sizes_[source_process], mpi_data_type<T>::get_type(),
                            //                       source_process, 0, comm_,
                            //                       &(recv_requests.back()));
                            //           }
                            //         }
                            unsigned rcv_ctr = 0;
                            for ( unsigned i = 0; i < sending_processes_.size ( ); ++i )
                            {
                                out[i] = std::vector<T>( transfer_sizes_[i] );

                                if ( transfer_sizes_[i] > 0 )
                                {
                                    ++rcv_ctr;
                                    int source_process = sending_processes_[i];
                                    recv_requests.push_back ( MPI_Request ( ) );

                                    MPI_Irecv ( util::raw_array ( out[i] ),
                                                transfer_sizes_[i], mpi_data_type<T>::get_type ( ),
                                                source_process, 0, comm_,
                                                &( recv_requests.back ( ) ) );
                                }
                            }
                            assert ( rcv_ctr == recv_requests.size ( ) );

                            std::vector<MPI_Status> recv_statuses ( rcv_ctr, MPI_Status ( ) );

                            MPI_Waitall ( rcv_ctr,
                                          util::raw_array ( recv_requests ),
                                          util::raw_array ( recv_statuses ) );
                        }

                        if ( process_is_sender_ )
                        {
                            std::vector<MPI_Status> send_statuses ( snd_ctr,
                                                                    MPI_Status ( ) );

                            MPI_Waitall ( snd_ctr,
                                          util::raw_array ( send_requests ),
                                          util::raw_array ( send_statuses ) );
                        }
                    }
                    if ( !process_is_receiver_ )
                    {
                        out = std::vector<std::vector<T> >( );
                    }
                }

                MPI_Comm comm_;
                int rank_;
                int comm_size_;

                bool process_is_sender_;
                bool process_is_receiver_;
                bool process_is_participating_;

                std::vector<int> sending_processes_;
                std::vector<int> receiving_processes_;

                // for send_vectors_by_reusing_sizes
                bool sizes_are_known_;
                std::vector<int> transfer_sizes_;
            };

            /// Implements the data transfer from a finer to a coarser grid.

            template<class LAD>
            class DataTransferFineToCoarse : public AsymmetricDataTransfer
            {
              public:
                /// Initializes the object.
                /// @param info An DataTransferInformation object that has been constructed
                /// with the coarse and and fine grids that shall be used by the transfer.

                explicit DataTransferFineToCoarse ( const DataTransferInformation<LAD>& info )
                : AsymmetricDataTransfer ( info.fine_level ( )->comm ( ),
                                           info.fine_procs_list ( ),
                                           info.coarse_procs_list ( ) )
                {
                }
            };

            /// Implements the data transfer from a coarser to a finer grid.

            template<class LAD>
            class DataTransferCoarseToFine : public AsymmetricDataTransfer
            {
              public:
                /// Initializes the object.
                /// @param info An DataTransferInformation object that has been constructed
                /// with the coarse and and fine grids that shall be used by the transfer.

                explicit DataTransferCoarseToFine ( const DataTransferInformation<LAD>& info )
                : AsymmetricDataTransfer ( info.fine_level ( )->comm ( ),
                                           info.coarse_procs_list ( ),
                                           info.fine_procs_list ( ) )
                {
                }
            };

        } // namespace gmg
    } // namespace la
} // namespace hiflow

#endif
