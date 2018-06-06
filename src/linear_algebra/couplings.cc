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

/// @author Chandramowli Subramanian, Martin Wlotzka

#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <set>

#include "couplings.h"
#include "dof/dof_partition.h"
#include "common/macros.h"
#include "mpi.h"

namespace hiflow
{
    namespace la
    {

        template<class DataType>
        Couplings<DataType>::Couplings ( )
        {
            this->Clear ( );
            this->comm_ = MPI_COMM_NULL;
            this->dof_partition_ = NULL;
            this->my_rank_ = -1;
        }

        template<class DataType>
        Couplings<DataType>::~Couplings ( )
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

        template<class DataType>
        void Couplings<DataType>::Init ( const MPI_Comm& comm, const typename doffem::DofPartition<DataType>& dp )
        {
            assert ( comm != MPI_COMM_NULL );
            int info = MPI_Comm_rank ( comm, &( this->my_rank_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->my_rank_ >= 0 );
            info = MPI_Comm_split ( comm, 0, this->my_rank_, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            this->dof_partition_ = &dp;
            assert ( this->dof_partition_ != NULL );
            assert ( this->my_rank_ == this->dof_partition_->get_my_subdomain ( ) );

            this->Clear ( );
        }

        template<class DataType>
        void Couplings<DataType>::InitializeCouplings ( const std::vector<int>& rows_offdiag,
                                                        const std::vector<int>& cols_offdiag )
        {
            assert ( this->dof_partition_ != NULL );
            assert ( this->my_rank_ == this->dof_partition_->get_my_subdomain ( ) );
            assert ( this->comm_ != MPI_COMM_NULL );
            assert ( rows_offdiag.size ( ) == cols_offdiag.size ( ) );

            this->Clear ( );

            // retrieve number of processes
            int nb_procs;
            int info = MPI_Comm_size ( this->comm_, &nb_procs );
            assert ( info == MPI_SUCCESS );
            assert ( nb_procs > this->my_rank_ );

            // get global offsets
            this->global_offsets_.resize ( nb_procs + 1 );
            this->global_offsets_[0] = 0;
            for ( int id = 1; id <= nb_procs; ++id )
            {
                this->global_offsets_[id] = this->global_offsets_[id - 1] +
                        this->dof_partition_->ndofs_on_sd ( id - 1 );
            }

            // a set of rows for every process id : vec[id] = { row_id_1, row_id_2, ...}
            std::vector< std::set<int> > cols_per_id ( nb_procs );

            for ( int i = 0; i < cols_offdiag.size ( ); ++i )
            {
                cols_per_id[this->dof_partition_->owner_of_dof ( cols_offdiag[i] )].insert ( cols_offdiag[i] );
            }
            // now the offdiag column indices are collected in cols_per_id

            // compute ghost offsets and indices (sorted per process by global column id)
            this->ghost_offsets_.resize ( nb_procs + 1 );
            int nb_ghost = 0;
            for ( int id = 0; id < nb_procs; ++id )
            {
                this->ghost_offsets_[id] = nb_ghost;
                nb_ghost += cols_per_id[id].size ( );
            }
            this->ghost_offsets_[nb_procs] = nb_ghost;

            std::vector<int> ghost_indices ( nb_ghost );
            std::set<int>::iterator set_it;
            int ind;
            for ( int id = 0; id < nb_procs; ++id )
            {
                ind = 0;
                set_it = cols_per_id[id].begin ( );
                while ( set_it != cols_per_id[id].end ( ) )
                {
                    ghost_indices[this->ghost_offsets_[id] + ind] = *set_it;
                    ind++;
                    set_it++;
                }
            }
            // now ghost_indices_ contains the global ghost indices

            // send border data as other processes need these informations for ghost layout
            int tag = 1;
            std::vector<MPI_Request> mpi_req_send ( nb_procs );
            std::vector<MPI_Request> mpi_req_recv ( nb_procs );

            // how much is received per process
            std::vector<int> nb_ghost_values ( nb_procs, 0 );
            std::vector<int> nb_border_values ( nb_procs, 0 );

            // send sizes
            for ( int id = 0; id < nb_procs; ++id )
            {
                if ( id != this->my_rank_ )
                {
                    nb_ghost_values[id] = this->ghost_offsets_[id + 1] - this->ghost_offsets_[id];
                    info = MPI_Isend ( &( nb_ghost_values[id] ), 1, MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( nb_border_values[id] ), 1, MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // make sure all data has been sended and received
            for ( int id = 0; id < nb_procs; ++id )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Wait ( &( mpi_req_send[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Wait ( &( mpi_req_recv[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // compute border offsets
            this->border_offsets_.resize ( nb_procs + 1 );
            this->border_offsets_[0] = 0;
            for ( int id = 1; id <= nb_procs; ++id )
            {
                this->border_offsets_[id] = this->border_offsets_[id - 1] + nb_border_values[id - 1];
            }
            this->border_indices_.resize ( this->border_offsets_[nb_procs] );

            // send ghost indices
            for ( int id = 0; id < nb_procs; ++id )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Isend ( &( ghost_indices[this->ghost_offsets_[id]] ),
                                       nb_ghost_values[id], MPI_INT, id, tag,
                                       this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( this->border_indices_[this->border_offsets_[id]] ),
                                       nb_border_values[id], MPI_INT, id, tag,
                                       this->comm_, &( mpi_req_recv[id] ) );
                }
            }

            // make sure all data has been sended and received
            for ( int id = 0; id < nb_procs; ++id )
            {
                if ( id != my_rank_ )
                {
                    info = MPI_Wait ( &( mpi_req_send[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Wait ( &( mpi_req_recv[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
            }
            // now border_indices_ contains the global border indices

            for ( int i = 0; i < this->border_indices_.size ( ); ++i )
            {
                ind = this->border_indices_[i];
                this->dof_partition_->global2local ( ind, &( this->border_indices_[i] ) );
            }

            // create global2offdiag and offdiag2global mappings
            this->CompressOffdiagonal ( cols_offdiag );

            // hopefully all things worked out fine
            this->initialized_ = true;
        }

        template class Couplings<double>;
        template class Couplings<float>;

    } // namespace la
} // namespace hiflow
