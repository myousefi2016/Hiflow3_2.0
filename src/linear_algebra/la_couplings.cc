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

#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <set>
#include <algorithm>
#include "cxx98_prettyprint.h"
#include "la_couplings.h"

namespace hiflow
{
    namespace la
    {

        LaCouplings::LaCouplings ( )
        {
            this->Clear ( );
            this->comm_ = MPI_COMM_NULL;
            this->my_rank_ = -1;
            this->initialized_ = false;
        }

        LaCouplings::~LaCouplings ( )
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

        void LaCouplings::Init ( const MPI_Comm& comm )
        {
            assert ( comm != MPI_COMM_NULL );
            int info = MPI_Comm_rank ( comm, &( this->my_rank_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->my_rank_ >= 0 );
            info = MPI_Comm_split ( comm, 0, this->my_rank_, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            this->Clear ( );
        }

        void LaCouplings::InitializeCouplings ( const std::vector<int>& global_offsets,
                                                const std::vector<int>& offdiag_cols,
                                                const std::vector<int>& offdiag_offsets )
        {
            // check input vectors
            int nb_procs = -1;
            int info = MPI_Comm_size ( this->comm_, &nb_procs );
            assert ( info == MPI_SUCCESS );
            assert ( nb_procs > this->my_rank_ );
            assert ( this->comm_ != MPI_COMM_NULL );
            assert ( static_cast < int > ( global_offsets.size ( ) ) == nb_procs + 1 );
            assert ( static_cast < int > ( offdiag_offsets.size ( ) ) == nb_procs + 1 );
            assert ( static_cast < int > ( offdiag_cols.size ( ) ) == offdiag_offsets.back ( ) );

            // clear content and begin
            this->Clear ( );

            // store global offsets
            this->global_offsets_.resize ( nb_procs + 1 );
            this->global_offsets_[0] = 0;
            for ( int i = 1; i <= nb_procs; ++i )
            {
                this->global_offsets_[i] = global_offsets[i];
            }

            // ghost variable  = variable which I need but another process owns it
            //                 = border variable of that other process
            // border variable = variable which I own but (several) other processes need it
            //                 = ghost variable of these other processes

            std::vector< std::set<int> > ghost_per_id ( nb_procs );

            // collect global ghost indices (i.e. the column indices in the offdiagonal part);
            // offdiag_cols may contain duplicates, these are skipped by std::map::insert
            for ( int id = 0; id < nb_procs; ++id )
            {
                for ( int j = offdiag_offsets[id], e_j = offdiag_offsets[id + 1]; j < e_j; j++ )
                {
                    ghost_per_id[id].insert ( offdiag_cols[j] );
                }
            }

            // compute ghost offsets
            this->ghost_offsets_.resize ( nb_procs + 1 );
            int nb_ghost = 0;
            for ( int id = 0; id < nb_procs; ++id )
            {
                this->ghost_offsets_[id] = nb_ghost;
                nb_ghost += ghost_per_id[id].size ( );
            }
            this->ghost_offsets_[nb_procs] = nb_ghost;
            // now ghost_offsets_.back() == nb. of ghost dofs

            // extract global ghost indices into a std::vector
            std::vector<int> ghost_global_ind ( this->ghost_offsets_.back ( ) );
            std::set<int>::iterator set_it;
            int ind;
            for ( int id = 0; id < nb_procs; ++id )
            {
                ind = 0;
                set_it = ghost_per_id[id].begin ( );
                while ( set_it != ghost_per_id[id].end ( ) )
                {
                    ghost_global_ind[this->ghost_offsets_[id] + ind] = *set_it;
                    ind++;
                    set_it++;
                }
            }

            // prepare for communication
            int tag = 1;
            std::vector<MPI_Request> mpi_req_send ( nb_procs );
            std::vector<MPI_Request> mpi_req_recv ( nb_procs );

            // how many border variables do other processes need
            std::vector<int> nb_border_values ( nb_procs, 0 );
            nb_border_values[this->my_rank_] = 0;

            // how many ghost variables are needed from other processes
            std::vector<int> nb_ghost_values ( nb_procs, 0 );
            nb_ghost_values[this->my_rank_] = this->ghost_offsets_[this->my_rank_ + 1] -
                    this->ghost_offsets_[this->my_rank_];
            assert ( nb_ghost_values[this->my_rank_] == 0 );

            // send number of ghost variables per process as this is the number of
            // border variables in view of the other processes
            // and recieve the number of ghost variables from the other processes as this
            // is the number of border variables in my view
            for ( int id = 0; id < nb_procs; ++id )
            {
                if ( id != this->my_rank_ )
                {
                    nb_ghost_values[id] = this->ghost_offsets_[id + 1] -
                            this->ghost_offsets_[id];
                    info = MPI_Isend ( &( nb_ghost_values[id] ), 1, MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( nb_border_values[id] ), 1, MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // make sure all data has been sent and received
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

            // std::vector to recieve global border indices
            //std::vector<int> border_global_ind(this->border_indices_.size());

            // send my global ghost indices as these are the global border indices in view
            // of the other processes
            // and recieve global ghost indices from the other processes as these are
            // the global border indices in my view
            for ( int id = 0; id < nb_procs; ++id )
            {
                if ( id != this->my_rank_ )
                {
                    info = MPI_Isend ( &( ghost_global_ind[this->ghost_offsets_[id]] ),
                                       nb_ghost_values[id], MPI_INT, id, tag,
                                       this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( this->border_indices_[this->border_offsets_[id]] ),
                                       nb_border_values[id], MPI_INT, id, tag,
                                       this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // make sure all data has been sent and received
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

            // shift global border indices to local numbering
            for ( size_t i = 0, e_i = this->border_indices_.size ( ); i != e_i; ++i )
            {
                assert ( this->border_indices_[i] >= this->global_offsets_[this->my_rank_] );
                assert ( this->border_indices_[i] < this->global_offsets_[this->my_rank_ + 1] );
                this->border_indices_[i] -= this->global_offsets_[this->my_rank_];
            }

            // create global2offdiag mapping
            this->CompressOffdiagonal ( offdiag_cols );

            // hopefully all things worked fine
            this->initialized_ = true;
        }

        void LaCouplings::Clear ( )
        {
            this->global_offsets_.clear ( );
            this->border_indices_.clear ( );
            this->border_offsets_.clear ( );
            this->ghost_offsets_.clear ( );
            this->global2offdiag_.clear ( );
            this->offdiag2global_.clear ( );
            this->initialized_ = false;
        }

        void LaCouplings::CompressOffdiagonal ( const std::vector<int>& offdiag_cols )
        {
            // sort by global column id that order matches with order of border indices
            std::vector<int> offdiag_cols_sorted ( offdiag_cols.size ( ) );
            std::copy ( offdiag_cols.begin ( ), offdiag_cols.end ( ),
                        offdiag_cols_sorted.begin ( ) );
            std::sort ( offdiag_cols_sorted.begin ( ), offdiag_cols_sorted.end ( ) );

            int ind = 0;
            for ( size_t k = 0, e_k = offdiag_cols_sorted.size ( ); k != e_k; ++k )
            {
                // insert column index to map if not done so far
                if ( this->global2offdiag_.find ( offdiag_cols_sorted[k] ) == this->global2offdiag_.end ( ) )
                {
                    this->global2offdiag_.insert ( std::pair<int, int>( offdiag_cols_sorted[k], ind ) );
                    this->offdiag2global_.insert ( std::pair<int, int>( ind, offdiag_cols_sorted[k] ) );
                    ind++;
                }
            }
            assert ( ind == this->size_ghost ( ) );
        }

        std::ostream& operator<< ( std::ostream& os, LaCouplings const& laCouplings )
        {
            return os << "global2offdiag: " << laCouplings.global2offdiag_ << "\n"
                    << "offdiag2global: " << laCouplings.offdiag2global_ << "\n"
                    << "global_offsets: " << laCouplings.global_offsets_ << "\n"
                    << "border_indices: " << laCouplings.border_indices_ << "\n"
                    << "border_offsets: " << laCouplings.border_offsets_ << "\n"
                    << "ghost_offsets: " << laCouplings.ghost_offsets_ << "\n"
                    << "initialized: " << laCouplings.initialized_ << "\n"
                    << "comm: " << laCouplings.comm_ << "\n"
                    << "my_rank: " << laCouplings.my_rank_ << "\n"
                    << std::endl;
        }

    } // namespace la
} // namespace hiflow
