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

/// @author Martin Wlotzka

#include <cassert>
#include <cmath>
#include <cstdlib>
#include "pp_vector.h"
#include "linear_algebra/la_descriptor.h"
#include "common/log.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        PpVector<LAD>::PpVector ( )
        {
            this->initialized_ = false;
            this->comm_ = MPI_COMM_NULL;
            this->dof_partition_ = NULL;
            this->vector_ = NULL;
            this->my_rank_ = -1;
            this->nb_procs_ = -1;
            this->dof_ids_.clear ( );
            this->offsets_.clear ( );
            this->values_.clear ( );
            this->nb_send_.clear ( );
            this->nb_recv_.clear ( );
            this->send_dofs_.clear ( );
            this->send_offsets_.clear ( );
            this->send_values_.clear ( );
        }

        template<class LAD>
        PpVector<LAD>::~PpVector ( )
        {
            this->Clear ( );
            int is_finalized;
            MPI_Finalized ( &is_finalized );
            if ( !is_finalized )
            {
                if ( this->comm_ != MPI_COMM_NULL )
                    MPI_Comm_free ( &this->comm_ );
            }
        }

        template<class LAD>
        void PpVector<LAD>::Init ( const MPI_Comm& comm, const doffem::DofPartition<typename LAD::DataType>& dp, const VectorType& vec )
        {
            assert ( comm != MPI_COMM_NULL );

            int info = MPI_Comm_size ( comm, &( this->nb_procs_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->nb_procs_ >= 1 );

            info = MPI_Comm_rank ( comm, &( this->my_rank_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );

            info = MPI_Comm_split ( comm, 0, this->my_rank_, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            this->dof_partition_ = &dp;
            assert ( this->dof_partition_ != NULL );

            this->vector_ = &vec;
            assert ( this->vector_ != NULL );
        }

        template<class LAD>
        void PpVector<LAD>::InitStructure ( )
        {
            assert ( this->comm_ != MPI_COMM_NULL );
            assert ( this->dof_partition_ != NULL );

            // clean up
            this->Clear ( );

            // get DoF-IDs and offsets for communication
            this->dof_partition_->GetPostProcessingStructure ( this->dof_ids_, this->offsets_ );
            assert ( this->offsets_.size ( ) == this->nb_procs_ + 1 );

            this->sorted_dof_ids_ = SortedArray<DofID>( this->dof_ids_.begin ( ), this->dof_ids_.end ( ) );

            this->values_.resize ( this->dof_ids_.size ( ), static_cast < DataType > ( 0. ) );

            int tag = 1;
            std::vector<MPI_Request> mpi_req_send ( this->nb_procs_ );
            std::vector<MPI_Request> mpi_req_recv ( this->nb_procs_ );

            // how much is to be sent to other processes
            this->nb_send_.resize ( nb_procs_, 0 );

            // how much do I request from others
            this->nb_recv_.resize ( nb_procs_, 0 );

            // exchange sizes
            int info = 0;
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( id != this->my_rank_ )
                {
                    this->nb_recv_[id] = this->offsets_[id + 1] - this->offsets_[id];
                    info = MPI_Isend ( &( this->nb_recv_[id] ), 1, MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( this->nb_send_[id] ), 1, MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // make sure all data has been sent and received
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Wait ( &( mpi_req_recv[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Wait ( &( mpi_req_send[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // allocate memory for DoF-IDs
            int nb_send_total = 0;

            for ( int i = 0; i < this->nb_procs_; i++ )
                nb_send_total += this->nb_send_[i];

            this->send_dofs_.resize ( nb_send_total, 0 );
            this->send_offsets_.resize ( this->nb_procs_ + 1, 0 );
            this->send_values_.resize ( nb_send_total, static_cast < DataType > ( 0. ) );

            for ( int i = 1; i < this->nb_procs_ + 1; i++ )
                this->send_offsets_[i] = this->send_offsets_[i - 1] + this->nb_send_[i - 1];

            // exchange DofIDs
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Isend ( &( this->dof_ids_[0] ) + this->offsets_[id],
                                       this->nb_recv_[id], MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( this->send_dofs_[0] ) + this->send_offsets_[id],
                                       this->nb_send_[id], MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // make sure all data has been sent and received
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Wait ( &( mpi_req_recv[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Wait ( &( mpi_req_send[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
            }

            this->initialized_ = true;
        }

        template<>
        void PpVector<LADescriptorCoupledD>::UpdateValues ( )
        {
            assert ( initialized_ );

            int tag = 1;
            int info;
            std::vector<MPI_Request> mpi_req_send ( nb_procs_ );
            std::vector<MPI_Request> mpi_req_recv ( nb_procs_ );

            // extract own values
            int n = this->offsets_[my_rank_ + 1] - this->offsets_[my_rank_];
            if ( n > 0 )
            {
                this->vector_->GetValues ( &( this->dof_ids_[this->offsets_[this->my_rank_]] ),
                                           n, &( this->values_[this->offsets_[this->my_rank_]] ) );
            }

            // extract values to send
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( ( id != this->my_rank_ ) && ( this->nb_send_[id] > 0 ) )
                {
                    this->vector_->GetValues ( &( this->send_dofs_[this->send_offsets_[id]] ),
                                               this->nb_send_[id],
                                               &( this->send_values_[this->send_offsets_[id]] ) );
                }
            }

            // exchange values
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Isend ( &( this->send_values_[0] ) + this->send_offsets_[id],
                                       this->nb_send_[id], MPI_DOUBLE,
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( this->values_[0] ) + this->offsets_[id],
                                       this->nb_recv_[id], MPI_DOUBLE,
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // make sure all data has been sent and received
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Wait ( &( mpi_req_recv[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Wait ( &( mpi_req_send[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
            }
        }

        template<>
        void PpVector<LADescriptorCoupledS>::UpdateValues ( )
        {
            assert ( initialized_ );

            int tag = 1;
            int info;
            std::vector<MPI_Request> mpi_req_send ( nb_procs_ );
            std::vector<MPI_Request> mpi_req_recv ( nb_procs_ );

            // extract own values
            int n = this->offsets_[my_rank_ + 1] - this->offsets_[my_rank_];
            if ( n > 0 )
            {
                this->vector_->GetValues ( &( this->dof_ids_[this->offsets_[this->my_rank_]] ),
                                           n, &( this->values_[this->offsets_[this->my_rank_]] ) );
            }

            // extract values to send
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( ( id != this->my_rank_ ) && ( this->nb_send_[id] > 0 ) )
                {
                    this->vector_->GetValues ( &( this->send_dofs_[this->send_offsets_[id]] ),
                                               this->nb_send_[id],
                                               &( this->send_values_[this->send_offsets_[id]] ) );
                }
            }

            // exchange values
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Isend ( &( this->send_values_[0] ) + this->send_offsets_[id],
                                       this->nb_send_[id], MPI_FLOAT,
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( this->values_[0] ) + this->offsets_[id],
                                       this->nb_recv_[id], MPI_FLOAT,
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // make sure all data has been sent and received
            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Wait ( &( mpi_req_recv[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Wait ( &( mpi_req_send[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
            }
        }

        template<class LAD>
        void PpVector<LAD>::GetDofsAndValues ( std::vector<DofID>& id, std::vector<DataType>& val ) const
        {
            id = this->dof_ids_;
            val = this->values_;
        }

        template<class LAD>
        void PpVector<LAD>::GetValues ( std::vector<DataType>& val ) const
        {
            val = this->values_;
        }

        template<class LAD>
        void PpVector<LAD>::GetValues ( const int* global_dof_ids, const int num_values, typename LAD::DataType* out_values ) const
        {
            int local_dof_id = -1;

            for ( int k = 0; k < num_values; ++k )
            {
                const bool found = this->sorted_dof_ids_.find ( *global_dof_ids, &local_dof_id );

                if ( !found )
                {
                    std::cerr << "Dof " << *global_dof_ids << " not found!\n";
                }

                assert ( found );

                *out_values = this->values_.at ( local_dof_id );
                global_dof_ids++;
                out_values++;
            }
        }

        template<class LAD>
        typename LAD::DataType PpVector<LAD>::at ( const int i ) const
        {
            assert ( this->initialized_ );
            assert ( i >= 0 );
            assert ( i < this->values_.size ( ) );
            return this->values_[i];
        }

        template<class LAD>
        void PpVector<LAD>::Clear ( )
        {
            this->initialized_ = false;
            this->my_rank_ = -1;
            this->nb_procs_ = -1;
            this->dof_ids_.clear ( );
            this->offsets_.clear ( );
            this->values_.clear ( );
            this->nb_send_.clear ( );
            this->nb_recv_.clear ( );
            this->send_dofs_.clear ( );
            this->send_offsets_.clear ( );
            this->send_values_.clear ( );
        }

        template<class LAD>
        void PpVector<LAD>::Print ( std::ostream& os ) const
        {
            os << "{ PpVector, rank = " << this->my_rank_ << "\n"
                    << "\tdof_ids = " << string_from_range ( this->dof_ids_.begin ( ), this->dof_ids_.end ( ) ) << "\n"
                    << "\tvalues = " << string_from_range ( this->values_.begin ( ), this->values_.end ( ) ) << "\n}\n";
        }

        /// template instantiation
        template class PpVector<LADescriptorCoupledD>;
        template class PpVector<LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow
