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

#include "vector_transfer.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            /// Enables the transfer of a vector from a fine grid to a coarser grid
            /// and vice versa

            template<typename LAD>
            VectorTransfer<LAD>::VectorTransfer (
                                                  const boost::shared_ptr<const DofIdentification<LAD> >& dof_ident,
                                                  const boost::shared_ptr<const DataTransferInformation<LAD> >& info,
                                                  const boost::shared_ptr<VectorType>& fine_vector,
                                                  const boost::shared_ptr<VectorType>& coarse_vector,
                                                  const boost::shared_ptr<DataTransferCoarseToFine<LAD> >& to_fine,
                                                  const boost::shared_ptr<DataTransferFineToCoarse<LAD> >& to_coarse )
            : dof_ident_ ( dof_ident ), info_ ( info ), fine_vector_ ( fine_vector ), coarse_vector_ ( coarse_vector )
            {
                assert ( info != NULL );
                assert ( dof_ident != NULL );
                assert ( to_fine != NULL );
                assert ( to_coarse != NULL );

                transfer_to_fine_ = to_fine;

                transfer_to_coarse_ = to_coarse;

                is_coarse_process_ = info_->coarse_level ( )->is_scheduled_to_this_process ( );
                is_fine_process_ = info_->fine_level ( )->is_scheduled_to_this_process ( );

                if ( is_fine_process_ )
                {
                    assert ( fine_vector != NULL );
                    assert ( dof_ident_->get_fine_dofs_to_transfer ( ).size ( ) ==
                             info_->num_coarse_procs ( ) );
                }

                if ( is_coarse_process_ )
                {
                    assert ( coarse_vector != NULL );
                    assert ( dof_ident_->get_coarse_dofs_to_transfer ( ).size ( ) ==
                             info_->num_fine_procs ( ) );
                }
            }

            /// The internal method that is used for the actual transfer of the vectors.
            /// Collective on the union of the sender and receiver processes.
            /// @param is_sender Whether this process is among the processes that send
            /// the vector
            /// @param is_receiver Whether this process is among the processes that
            /// receive the vector
            /// @param sender_process_list A list of the ranks in the fine communicator
            /// of all processes that send the vector
            /// @param receiver_process_list A list of the ranks in the fine communicator
            /// of all processes that receive the vector
            /// @param sending_vector The vector that is to be sent. This only needs
            /// to be a valid pointer, if the calling process is among the sending
            /// processes.
            /// @param receiving_vector Where the received vector will be stored.
            /// This only needs to be a valid pointer, if the calling process is among
            /// the receiving processes.
            /// @param dofs_at_sender The vector obtained from the DofIdentification
            /// that specifies which DoFs have to be transferred to which receiving process.
            /// @param dofs_at_receiver The vector obtained from the DofIdentification
            /// that specifies which DoFs are received from which sending process.
            /// @param data_transfer The data transfer object used for the transfer

            template<typename LAD>
            void VectorTransfer<LAD>::transfer ( bool is_sender,
                                                 bool is_receiver,
                                                 const std::vector<int>& sender_process_list,
                                                 const std::vector<int>& receiver_process_list,
                                                 const boost::shared_ptr<VectorType>& sending_vector,
                                                 const boost::shared_ptr<VectorType>& receiving_vector,
                                                 const std::vector<std::vector<int> >& dofs_at_sender,
                                                 const std::vector<std::vector<int> >& dofs_at_receiver,
                                                 AsymmetricDataTransfer& data_transfer )
            {
                if ( is_sender || is_receiver )
                {
                    std::vector<std::vector<DataType> > received_data;

                    //     std::vector<std::vector<DataType> > data_to_send(info_->num_fine_procs(), std::vector<DataType>());
                    std::vector<std::vector<DataType> > data_to_send;

                    if ( is_sender )
                    {
                        data_to_send = std::vector<std::vector<DataType> >( receiver_process_list.size ( ), std::vector<DataType>( ) );
                        // For each process that receives data (i.e., a process
                        // to which data needs to be sent), obtain the data that
                        // has to be transferred from the vector
                        for ( unsigned i = 0; i < receiver_process_list.size ( ); ++i )
                        {
                            int receiving_process = receiver_process_list[i];

                            int num_dofs = dofs_at_sender[i].size ( );

                            data_to_send[i] = std::vector<DataType>( num_dofs, DataType ( ) );

                            // Make sure we have data to collect (it is possible, that
                            // we do not share any dofs with the receiving_process).
                            // A debug assertion in CoupledVector::GetValues fails
                            // if num_dofs == 0
                            if ( num_dofs != 0 )
                            {
                                // We know which DoFs to send to the process of
                                // rank receiving_process from the DofIdentification
                                // result stored in dofs_at_sender
                                sending_vector->GetValues (
                                                            util::raw_array ( dofs_at_sender[i] ),
                                                            num_dofs,
                                                            util::raw_array ( data_to_send[i] ) );
                            }
                        }

                    }

                    data_transfer.transfer_vectors_by_reusing_sizes<DataType>(
                            data_to_send,
                            received_data );

                    if ( is_receiver )
                    {
                        assert ( received_data.size ( ) == sender_process_list.size ( ) );
                        // Empty the receiving vector
                        receiving_vector->Zeros ( );

                        // For each process that sends, put the received
                        // data into the receiving_vector
                        for ( unsigned i = 0; i < sender_process_list.size ( ); ++i )
                        {
                            int sender_process = sender_process_list[i];

                            int num_dofs = received_data[i].size ( );
                            assert ( num_dofs == dofs_at_receiver[i].size ( ) );
                            if ( num_dofs != 0 )
                            {
                                // We know where to put the received values
                                // from the results of the DofIdentification
                                // stored in dofs_at_receiver.
                                receiving_vector->SetValues (
                                                              util::raw_array ( dofs_at_receiver[i] ),
                                                              num_dofs,
                                                              util::raw_array ( received_data[i] ) );
                            }
                        }
                    }
                }
            }

            template class VectorTransfer<LADescriptorCoupledD>;
            template class VectorTransfer<LADescriptorCoupledS>;

        } // namespace gmg
    } // namespace la
} // namespace hiflow
