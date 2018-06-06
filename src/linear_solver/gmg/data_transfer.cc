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

#include "data_transfer.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            /// Initializes the object. Collective on the communicator of
            /// the finer grid.
            /// @param coarse_level A pointer to a coarser grid of the MultiLevelHierarchy
            /// @param fine_level A pointer to a finer grid of the MultiLevelHierarchy
            /// While the two grids do not neccessarily have to stem from the
            /// same MultiLevelHierarchy object, the processes of the coarse grid
            /// communicator have to be subset of the processes in the communicator
            /// of the fine grid.

            template<class LAD>
            DataTransferInformation<LAD>::DataTransferInformation ( const BasicLevel<LAD>* coarse_level,
                                                                    const BasicLevel<LAD>* fine_level )
            : num_fine_procs_ ( 0 ), num_coarse_procs_ ( 0 ), fine_rank_ ( -1 ), coarse_rank_ ( -1 ),
            coarse_procs_ranks_ ( ), coarse_procs_list_ ( ), fine_procs_list_ ( ),
            coarse_level_ ( coarse_level ), fine_level_ ( fine_level )
            {
                assert ( fine_level_ != NULL );
                assert ( coarse_level_ != NULL );

                if ( fine_level_->is_scheduled_to_this_process ( ) )
                {
                    MPI_Comm_size ( fine_level_->comm ( ), &num_fine_procs_ );
                    MPI_Comm_rank ( fine_level_->comm ( ), &fine_rank_ );
                }

                if ( coarse_level_->is_scheduled_to_this_process ( ) )
                {
                    assert ( fine_level_->is_scheduled_to_this_process ( ) );

                    MPI_Comm_size ( coarse_level_->comm ( ), &num_coarse_procs_ );
                    MPI_Comm_rank ( coarse_level_->comm ( ), &coarse_rank_ );

                    assert ( num_coarse_procs_ <= num_fine_procs_ );
                }

                fetch_coarse_process_information ( );

                //fill fine_process_list
                if ( fine_level_->is_scheduled_to_this_process ( ) )
                {
                    fine_procs_list_ = std::vector<int>( num_fine_procs_, 0 );
                    for ( int i = 0; i < num_fine_procs_; ++i )
                    {
                        fine_procs_list_[i] = i;
                    }
                }
                assert ( num_fine_procs_ == fine_procs_list_.size ( ) );
            }

            /// Fetches the ranks of all coarse processes, build the coarse process list
            /// stored in coarse_procs_list_ and counts the number of coarse processes
            /// which will be stored in num_coarse_procs_
            /// Collective on the communicator of the fine grid.

            template<class LAD>
            void DataTransferInformation<LAD>::fetch_coarse_process_information ( )
            {
                if ( fine_level_->is_scheduled_to_this_process ( ) )
                {
                    coarse_procs_ranks_ = std::vector<int>( num_fine_procs_, -1 );

                    int transmitted_rank = -1;
                    if ( coarse_level_->is_scheduled_to_this_process ( ) )
                        //transmitting the coarse rank will enable us to translate
                        //between coarse and fine processes
                        transmitted_rank = coarse_rank_;

                    MPI_Allgather ( &transmitted_rank, 1, MPI_INT, util::raw_array ( coarse_procs_ranks_ ),
                                    1, MPI_INT, fine_level_->comm ( ) );

                    for ( int i = 0; i < num_fine_procs_; ++i )
                    {
                        // coarse processes have transmitted their rank in the
                        // fine communicator, while processes that are not
                        // in the coarse communicator have transferred -1
                        if ( coarse_procs_ranks_[i] != -1 )
                        {
                            //i is the fine rank of the process
                            coarse_procs_list_.push_back ( i );
                        }
                    }

                    if ( !coarse_level_->is_scheduled_to_this_process ( ) )
                    {
                        num_coarse_procs_ = coarse_procs_list_.size ( );
                    }
                }
                assert ( num_coarse_procs_ == coarse_procs_list_.size ( ) );
            }

            template class DataTransferInformation<LADescriptorCoupledD>;
            template class DataTransferInformation<LADescriptorCoupledS>;

            AsymmetricDataTransfer::AsymmetricDataTransfer ( MPI_Comm comm,
                                                             const std::vector<int>& source_process_list,
                                                             const std::vector<int>& dest_process_list )
            : comm_ ( comm ), rank_ ( -1 ), comm_size_ ( 0 ),
            process_is_sender_ ( false ), process_is_receiver_ ( false ), process_is_participating_ ( false ),
            sending_processes_ ( source_process_list ),
            receiving_processes_ ( dest_process_list ),
            sizes_are_known_ ( false ),
            transfer_sizes_ ( )
            {
                if ( comm != MPI_COMM_NULL )
                {
                    MPI_Comm_rank ( comm_, &rank_ );
                    MPI_Comm_size ( comm_, &comm_size_ );

                    std::vector<int>::iterator
                    it = std::find ( sending_processes_.begin ( ), sending_processes_.end ( ), rank_ );

                    process_is_sender_ = ( it != sending_processes_.end ( ) );

                    it = std::find ( receiving_processes_.begin ( ), receiving_processes_.end ( ), rank_ );

                    process_is_receiver_ = ( it != receiving_processes_.end ( ) );

                    assert ( process_is_sender_ || process_is_receiver_ );
                }

                process_is_participating_ = process_is_receiver_ || process_is_sender_;
            }

            template class DataTransferFineToCoarse<LADescriptorCoupledD>;
            template class DataTransferFineToCoarse<LADescriptorCoupledS>;

            template class DataTransferCoarseToFine<LADescriptorCoupledD>;
            template class DataTransferCoarseToFine<LADescriptorCoupledS>;

        } // namespace gmg
    } // namespace la
} // namespace hiflow
