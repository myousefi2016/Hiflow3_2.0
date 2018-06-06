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

/// @author Chandramowli Subramanian, Nico Trost, Dimitar Lukarski, Martin Wlotzka, Simon Gawlok

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "config.h"

#include "linear_algebra/coupled_vector.h"
#include "lmp/init_vec_mat.h"
#include "common/pointers.h"

#ifdef WITH_HDF5
#    include "hdf5.h"
#    include "common/hdf5_tools.h"
#endif

#include "common/log.h"
#include "tools/mpi_tools.h"

const int DEBUG_LEVEL = 0;
namespace hiflow
{
    namespace la
    {

        template<class DataType>
        CoupledVector<DataType>::CoupledVector ( )
        {
            this->la_couplings_ = NULL;
            this->comm_ = MPI_COMM_NULL;
            this->nb_procs_ = -1;
            this->my_rank_ = -1;
            this->interior_ = NULL;
            this->ghost_ = NULL;
            this->ownership_begin_ = -1;
            this->ownership_end_ = -1;
            this->pp_data_ = NULL;
            this->checked_for_dof_partition_ = false;
        }

        template<class DataType>
        CoupledVector<DataType>::~CoupledVector ( )
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
        void CoupledVector<DataType>::Init ( const MPI_Comm& comm, const LaCouplings& cp,
                                             PLATFORM plat, IMPLEMENTATION impl )
        {
            if ( this->comm_ != MPI_COMM_NULL )
                MPI_Comm_free ( &this->comm_ );

            assert ( comm != MPI_COMM_NULL );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( comm, &( this->nb_procs_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( comm, &( this->my_rank_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );

            //info = MPI_Comm_split ( comm, 0, this->my_rank_, &( this->comm_ ) );
            info = MPI_Comm_dup ( comm, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            // couplings
            this->la_couplings_ = &cp;
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );

            // clear old data
            this->Clear ( );

            // init interior
            this->interior_ = init_vector<DataType>( 0, "interior", plat, impl );
            assert ( this->interior_ != NULL );
            // init ghost
            this->ghost_ = init_vector<DataType>( 0, "ghost", plat, impl );
            assert ( this->ghost_ != NULL );
        }

        template<class DataType>
        void CoupledVector<DataType>::Init ( const MPI_Comm& comm, const LaCouplings& cp,
                                             PLATFORM plat, IMPLEMENTATION impl,
                                             const SYSTEM& my_system )
        {
            if ( this->comm_ != MPI_COMM_NULL )
                MPI_Comm_free ( &this->comm_ );

            assert ( comm != MPI_COMM_NULL );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( comm, &( this->nb_procs_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( comm, &( this->my_rank_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );

            //info = MPI_Comm_split ( comm, 0, this->my_rank_, &( this->comm_ ) );
            info = MPI_Comm_dup ( comm, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            // couplings
            this->la_couplings_ = &cp;
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );

            // clear old data
            this->Clear ( );

            if ( plat == OPENCL )
            {
                // init interior
                this->interior_ = init_vector<DataType>( 0, "interior", plat, impl, my_system );
                assert ( this->interior_ != NULL );
                // init ghost
                this->ghost_ = init_vector<DataType>( 0, "ghost", plat, impl, my_system );
                assert ( this->ghost_ != NULL );
            }
            else
            {
                // init interior
                this->interior_ = init_vector<DataType>( 0, "interior", plat, impl );
                assert ( this->interior_ != NULL );
                // init ghost
                this->ghost_ = init_vector<DataType>( 0, "ghost", plat, impl );
                assert ( this->ghost_ != NULL );
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::Init ( const MPI_Comm& comm,
                                             const LaCouplings& cp )
        {
            if ( this->comm_ != MPI_COMM_NULL )
                MPI_Comm_free ( &this->comm_ );

            assert ( comm != MPI_COMM_NULL );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( comm, &( this->nb_procs_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( comm, &( this->my_rank_ ) );
            assert ( info == MPI_SUCCESS );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );

            //info = MPI_Comm_split ( comm, 0, this->my_rank_, &( this->comm_ ) );
            info = MPI_Comm_dup ( comm, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            // couplings
            this->la_couplings_ = &cp;
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );

            this->Init_la_system ( CPU, NAIVE );
            this->InitStructure ( );
        }

        template<class DataType>
        void CoupledVector<DataType>::Init_la_system ( PLATFORM plat,
                                                       IMPLEMENTATION impl )
        {
            // first clear old data
            this->Clear ( );
            // init interior
            this->interior_ = init_vector<DataType>( 0, "interior", plat, impl );
            assert ( this->interior_ != NULL );
            // init ghost
            this->ghost_ = init_vector<DataType>( 0, "ghost", plat, impl );
            assert ( this->ghost_ != NULL );
        }

        template<class DataType>
        void CoupledVector<DataType>::InitStructure ( )
        {
            assert ( this->comm_ != MPI_COMM_NULL );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );
            assert ( this->interior_ != NULL );
            assert ( this->ghost_ != NULL );

            // compute ownership range
            this->ComputeOwnershipRange ( );

            // init structure of interior
            this->interior_->Clear ( );
            this->interior_->Init ( this->la_couplings_->nb_dofs ( this->my_rank_ ),
                                    "interior" );

            // set border indices
            this->interior_->set_indexset ( this->la_couplings_->border_indices ( ),
                                            this->la_couplings_->size_border_indices ( ) );

            // init structure of ghost
            if ( this->la_couplings_->size_ghost ( ) > 0 )
            {
                this->ghost_->Clear ( );
                this->ghost_->Init ( this->la_couplings_->size_ghost ( ),
                                     "ghost" );
            }

            // prepare for communication
            this->nb_sends_ = 0;
            this->nb_recvs_ = 0;
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->la_couplings_->border_offsets ( id + 1 ) -
                     this->la_couplings_->border_offsets ( id ) > 0 )
                {
                    this->nb_sends_++;
                }
                if ( this->la_couplings_->ghost_offsets ( id + 1 ) -
                     this->la_couplings_->ghost_offsets ( id ) > 0 )
                {
                    this->nb_recvs_++;
                }
            }
            this->mpi_req_.resize ( this->nb_sends_ + this->nb_recvs_ );
            this->mpi_stat_.resize ( this->nb_sends_ + this->nb_recvs_ );
            this->border_val_.resize ( this->la_couplings_->size_border_indices ( ) );
            this->ghost_val_.resize ( this->la_couplings_->size_ghost ( ) );
        }

        template<class DataType>
        void CoupledVector<DataType>::Zeros ( )
        {
            assert ( this->interior_ != NULL );
            assert ( this->ghost_ != NULL );

            this->interior_->Zeros ( );
            this->ghost_->Zeros ( );

            if ( this->pp_data_ != NULL )
            {
                for ( typename std::vector<DataType>::iterator iter = this->pp_data_->values.begin ( ),
                      e_iter = this->pp_data_->values.end ( );
                      iter != e_iter; ++iter )
                {
                    *iter = 0.;
                }
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::Add ( int global_dof_id, DataType val )
        {
            assert ( this->interior_ != NULL );
            assert ( this->ownership_begin_ <= global_dof_id );
            assert ( global_dof_id < this->ownership_end_ );

            interior_->add_value ( global_dof_id - this->ownership_begin_, val );
        }

        template<class DataType>
        void CoupledVector<DataType>::Add ( const int* indices, const int size_indices, const DataType* values )
        {
            assert ( this->interior_ != NULL );

            std::vector<int> shifted_indices ( size_indices );
            for ( int i = 0; i < size_indices; ++i )
            {
                assert ( this->ownership_begin_ <= indices[i] );
                assert ( indices[i] < this->ownership_end_ );
                shifted_indices[i] = indices[i] - this->ownership_begin_;
            }
            interior_->add_values ( vec2ptr ( shifted_indices ), size_indices, values );
        }

        template<class DataType>
        DataType CoupledVector<DataType>::Dot ( const Vector<DataType>& vec ) const
        {

            const CoupledVector<DataType> *cv_vec;

            cv_vec = dynamic_cast < const CoupledVector<DataType>* > ( &vec );

            if ( cv_vec != 0 )
            {
                return this->Dot ( *cv_vec );
            }
            else
            {
                LOG_ERROR ( "Called CoupledVector::Dot with incompatible argument type." );
                exit ( -1 );
                return 0.0;
            }
        }

        template<class DataType>
        DataType CoupledVector<DataType>::Dot ( const CoupledVector<DataType>& vec ) const
        {
            assert ( this->comm_ != MPI_COMM_NULL );
            assert ( this->interior_ != NULL );

            // local dot product
            DataType dot_local = this->interior_->Dot ( vec.interior ( ) );

            // now sum up
            DataType dot_global;
            MPI_Allreduce ( &dot_local, &dot_global, 1, mpi_data_type<DataType>::get_type ( ), MPI_SUM, this->comm_ );

            return dot_global;
        }

        template<class DataType>
        void CoupledVector<DataType>::Axpy ( const Vector<DataType>& vec, const DataType alpha )
        {
            const CoupledVector<DataType> *cv_vec;

            cv_vec = dynamic_cast < const CoupledVector<DataType>* > ( &vec );

            if ( cv_vec != 0 )
            {
                this->Axpy ( *cv_vec, alpha );
            }
            else
            {
                LOG_ERROR ( "Called CoupledVector::Axpy with incompatible argument type." );
                exit ( -1 );
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::Axpy ( const CoupledVector<DataType>& vec,
                                             const DataType alpha )
        {
            assert ( this->interior_ != NULL );

            this->interior_->Axpy ( vec.interior ( ), alpha );
        }

        template<class DataType>
        void CoupledVector<DataType>::ScaleAdd ( const Vector<DataType>& vec, const DataType alpha )
        {
            const CoupledVector<DataType> *cv_vec;

            cv_vec = dynamic_cast < const CoupledVector<DataType>* > ( &vec );

            if ( cv_vec != 0 )
            {
                this->ScaleAdd ( *cv_vec, alpha );
            }
            else
            {
                LOG_ERROR ( "Called CoupledVector::ScaleAdd with incompatible argument type." );
                exit ( -1 );
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::ScaleAdd ( const CoupledVector<DataType>& vec,
                                                 const DataType alpha )
        {
            assert ( this->interior_ != NULL );

            this->interior_->ScaleAdd ( alpha, vec.interior ( ) );
        }

        template<class DataType>
        void CoupledVector<DataType>::Scale ( const DataType alpha )
        {
            assert ( this->interior_ != NULL );

            this->interior_->Scale ( alpha );
        }

        template<class DataType>
        DataType CoupledVector<DataType>::Norm1 ( ) const
        {
            DataType norm1_local = this->interior_->Norm1 ( );
            assert ( norm1_local >= 0.0 );
            DataType norm1_global;
            int info = MPI_Allreduce ( &norm1_local, &norm1_global, 1, mpi_data_type<DataType>::get_type ( ), MPI_SUM, this->comm_ );
            assert ( info == MPI_SUCCESS );
            assert ( norm1_global >= 0.0 );
            return norm1_global;
        }

        template<class DataType>
        DataType CoupledVector<DataType>::NormMax ( ) const
        {
            DataType norm_max_local = this->interior_->NormMax ( );
            assert ( norm_max_local >= 0.0 );
            DataType norm_max_global;
            int info = MPI_Allreduce ( &norm_max_local, &norm_max_global, 1, mpi_data_type<DataType>::get_type ( ), MPI_MAX, this->comm_ );
            assert ( info == MPI_SUCCESS );
            assert ( norm_max_global >= 0.0 );
            return norm_max_global;
        }

        template<class DataType>
        DataType CoupledVector<DataType>::Norm2 ( ) const
        {
            return sqrt ( this->Dot ( *this ) );
        }

        template<class DataType>
        void CoupledVector<DataType>::GetValues ( const int* indices,
                                                  const int size_indices,
                                                  DataType* values ) const
        {
            assert ( this->interior_ != NULL );
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );

            if ( size_indices > 0 )
            {
                assert ( indices != 0 );
                assert ( values != 0 );

                // extract interior and ghost indices (and if possible pp inidices)
                std::vector<int> indices_interior ( size_indices );
                int ind_interior = 0;

                std::vector<int> indices_ghost ( size_indices );
                int ind_ghost = 0;

                std::vector<int> indices_pp ( size_indices );
                int ind_pp = 0;

                // transform indices according to local numbering of interior and ghost
                for ( int i = 0; i < size_indices; i++ )
                {
                    // interior
                    if ( this->ownership_begin_ <= indices[i] && indices[i] < this->ownership_end_ )
                    {
                        indices_interior[ind_interior] = indices[i] - this->ownership_begin_;
                        ind_interior++;
                    }

                        // ghost
                    else if ( this->la_couplings_->global2offdiag ( ).find ( indices[i] ) !=
                              this->la_couplings_->global2offdiag ( ).end ( ) )
                    {
                        indices_ghost[ind_ghost] = this->la_couplings_->Global2Offdiag ( indices[i] );
                        ind_ghost++;
                    }

                        // pp values
                    else if ( ( this->pp_data_ != NULL ) && ( this->pp_data_->sorted_dof_ids.find ( indices[i], &( indices_pp[ind_pp] ) ) ) )
                    {
                        ind_pp++;
                    }
                }

                // extract values
                std::vector<DataType> values_interior ( ind_interior );
                if ( ind_interior > 0 )
                {
                    this->interior_->GetValues ( &indices_interior.front ( ), ind_interior,
                                                 &values_interior.front ( ) );
                }

                std::vector<DataType> values_ghost ( ind_ghost );

                if ( ind_ghost > 0 )
                {
                    this->ghost_->GetValues ( &indices_ghost.front ( ), ind_ghost, &values_ghost.front ( ) );
                }

                std::vector<DataType> values_pp ( ind_pp );
                if ( this->pp_data_ != NULL )
                {
                    for ( int i = 0; i < ind_pp; ++i )
                    {
                        values_pp[i] = this->pp_data_->values[indices_pp[i]];
                    }
                }
                // put values together
                for ( int i = 0; i < size_indices; i++ )
                {
                    int cont;
                    // interior
                    if ( this->ownership_begin_ <= indices[size_indices - i - 1] &&
                         indices[size_indices - i - 1] < this->ownership_end_ )
                    {
                        values[size_indices - i - 1] = values_interior[ind_interior - 1];
                        --ind_interior;
                    }

                        // ghost
                    else if ( this->la_couplings_->global2offdiag ( ).find ( indices[size_indices - i - 1] ) !=
                              this->la_couplings_->global2offdiag ( ).end ( ) )
                    {
                        values[size_indices - i - 1] = values_ghost[ind_ghost - 1];
                        --ind_ghost;
                    }

                        // pp
                    else if ( ( this->pp_data_ != NULL ) &&
                              ( this->pp_data_->sorted_dof_ids.find ( indices[size_indices - i - 1], &cont ) ) )
                    {
                        values[size_indices - i - 1] = values_pp[ind_pp - 1];
                        --ind_pp;
                    }
                        // zeros
                    else
                    {
                        values[size_indices - i - 1] = static_cast < DataType > ( 0. );
                    }
                }

                assert ( ind_interior == 0 );
                assert ( ind_ghost == 0 );
                assert ( ind_pp == 0 );
            }
        }

        template<class DataType>
        DataType CoupledVector<DataType>::GetValue ( const int index ) const
        {
            assert ( index >= this->ownership_begin_ );
            assert ( index < this->ownership_end_ );
            DataType result;
            this->GetValues ( &index, 1, &result );
            return result;
        }

        template<class DataType>
        void CoupledVector<DataType>::GetValues ( DataType* values ) const
        {
            assert ( this->interior_ != NULL );

            this->interior_->GetBlockValues ( 0, this->size_local ( ), values );
        }

        template<class DataType>
        void CoupledVector<DataType>::SetValue ( const int index, const DataType value )
        {
            assert ( index >= this->ownership_begin_ );
            assert ( index < this->ownership_end_ );

            const int local_index = index - this->ownership_begin_;
            this->interior_->SetValues ( &local_index, 1, &value );
        }

        template<class DataType>
        void CoupledVector<DataType>::SetValues ( const int* indices,
                                                  const int size_indices,
                                                  const DataType* values )
        {
            assert ( this->interior_ != NULL );

#ifndef NDEBUG
            int own_beg, own_end;
            this->GetOwnershipRange ( &own_beg, &own_end );
            for ( int i = 0; i < size_indices; i++ )
            {
                assert ( own_beg <= indices[i] && indices[i] < own_end );
                //LOG_DEBUG(2, "index " << indices[i] << " begin " << this->ownership_begin_ <<" end " << this->ownership_end_);
                //assert( this->ownership_begin_ <= indices[i] && indices[i] < this->ownership_end_  );
            }
#endif  // NDEBUG

            // shift indices to local numbering
            std::vector<int> indices_shifted ( size_indices );

            for ( int i = 0; i < size_indices; ++i )
            {
                indices_shifted[i] = indices[i] - this->ownership_begin_;
            }

            this->interior_->SetValues ( &indices_shifted.front ( ), size_indices, values );
        }

        template<class DataType>
        void CoupledVector<DataType>::SetValues ( const DataType* values )
        {
            assert ( this->interior_ != NULL );

            this->interior_->SetBlockValues ( 0, this->size_local ( ), values );
        }

        template<class DataType>
        void CoupledVector<DataType>::SetGhostValues ( const DataType* values )
        {
            assert ( this->ghost_ != NULL );

            this->ghost_->SetBlockValues ( 0, this->size_local_ghost ( ), values );
        }

        template<class DataType>
        void CoupledVector<DataType>::Gather ( int recv_id, DataType* values ) const
        {
            assert ( this->interior_ != NULL );
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );

            if ( this->interior ( ).get_platform ( ) == CPU )
            {
                // dynamic cast to CPU vector in order to access buffer
                const CPU_lVector<DataType>* casted_vec =
                        dynamic_cast < const CPU_lVector<DataType>* > ( this->interior_ );
                assert ( casted_vec != NULL );

                int comm_size;
                int info = MPI_Comm_size ( this->comm_, &comm_size );
                assert ( info == MPI_SUCCESS );
                int tag = 1;

                // receive
                if ( this->my_rank_ == recv_id )
                {

                    int recv_index_begin;
                    int recv_index_end = 0;

                    for ( int id = 0; id < comm_size; ++id )
                    {
                        recv_index_begin = recv_index_end;
                        recv_index_end += this->la_couplings_->nb_dofs ( id );

                        if ( id != this->my_rank_ )
                        {
                            MPI_Status status;
                            info = MPI_Recv ( &values[recv_index_begin],
                                              recv_index_end - recv_index_begin,
                                              mpi_data_type<DataType>::get_type ( ), id, tag, this->comm_, &status );
                            assert ( info == MPI_SUCCESS );
                        }

                        else
                        {
                            for ( int i = recv_index_begin; i < recv_index_end; i++ )
                            {
                                values[i] = casted_vec->buffer[i - recv_index_begin];
                            }
                        }
                    }
                }

                    // send
                else
                {
                    info = MPI_Send ( casted_vec->buffer, this->size_local ( ),
                                      mpi_data_type<DataType>::get_type ( ), recv_id, tag, this->comm_ );
                    assert ( info == MPI_SUCCESS );
                }

            }

            else
            {
                LOG_ERROR ( "CoupledVector::Gather: not supported on non-CPU platforms." );
                exit ( -1 );
            }
        }

        template<class DataType>
        Vector<DataType>* CoupledVector<DataType>::Clone ( ) const
        {
            CoupledVector<DataType> * clone = new CoupledVector<DataType>( );

            clone->CloneFrom ( *this );

            return clone;
        }

        template<class DataType>
        void CoupledVector<DataType>::CloneFrom ( const CoupledVector<DataType>& vec )
        {
            if ( this != &vec )
            {
                // clone vector and structure
                this->CloneFromWithoutContent ( vec );

                // copy entries
                this->interior_->CopyFrom ( vec.interior ( ) );
                this->ghost_->CopyFrom ( vec.ghost ( ) );

                if ( vec.HasPpData ( ) )
                {
                    this->pp_data_ = vec.pp_data ( ).Clone ( );
                }
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::CopyFrom ( const CoupledVector<DataType>& vec )
        {
            if ( this != &vec )
            {
                assert ( this->nb_procs_ == vec.nb_procs ( ) );
                assert ( this->my_rank_ == vec.my_rank ( ) );
                assert ( this->ownership_begin_ == vec.ownership_begin ( ) );
                assert ( this->ownership_end_ == vec.ownership_end ( ) );

                this->interior_->CopyFrom ( vec.interior ( ) );
                this->ghost_->CopyFrom ( vec.ghost ( ) );

                if ( vec.HasPpData ( ) )
                {
                    if ( this->pp_data_ != NULL )
                    {
                        *( this->pp_data_ ) = vec.pp_data ( );
                    }
                    else
                    {
                        this->pp_data_ = vec.pp_data ( ).Clone ( );
                    }
                }
                else
                {
                    if ( this->pp_data_ != NULL )
                    {
                        delete this->pp_data_;
                    }
                    this->pp_data_ = NULL;
                }
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::CastInteriorFrom ( const CoupledVector<double>& other )
        {
            assert ( this->nb_procs_ == other.nb_procs ( ) );
            assert ( this->my_rank_ == other.my_rank ( ) );
            assert ( this->ownership_begin_ == other.ownership_begin ( ) );
            assert ( this->ownership_end_ == other.ownership_end ( ) );

            this->interior_->CastFrom ( other.interior ( ) );

            return;
        }

        template<class DataType>
        void CoupledVector<DataType>::CastInteriorFrom ( const CoupledVector<float>& other )
        {
            assert ( this->nb_procs_ == other.nb_procs ( ) );
            assert ( this->my_rank_ == other.my_rank ( ) );
            assert ( this->ownership_begin_ == other.ownership_begin ( ) );
            assert ( this->ownership_end_ == other.ownership_end ( ) );

            this->interior_->CastFrom ( other.interior ( ) );

            return;
        }

        template<class DataType>
        void CoupledVector<DataType>::CastInteriorTo ( CoupledVector<double>& other ) const
        {
            assert ( this->nb_procs_ == other.nb_procs ( ) );
            assert ( this->my_rank_ == other.my_rank ( ) );
            assert ( this->ownership_begin_ == other.ownership_begin ( ) );
            assert ( this->ownership_end_ == other.ownership_end ( ) );

            this->interior_->CastTo ( other.interior ( ) );

            return;
        }

        template<class DataType>
        void CoupledVector<DataType>::CastInteriorTo ( CoupledVector<float>& other ) const
        {
            assert ( this->nb_procs_ == other.nb_procs ( ) );
            assert ( this->my_rank_ == other.my_rank ( ) );
            assert ( this->ownership_begin_ == other.ownership_begin ( ) );
            assert ( this->ownership_end_ == other.ownership_end ( ) );

            this->interior_->CastTo ( other.interior ( ) );

            return;
        }

        template<class DataType>
        void CoupledVector<DataType>::CopyTo ( CoupledVector<DataType>& vec ) const
        {
            if ( this != &vec )
            {
                assert ( this->nb_procs_ == vec.nb_procs ( ) );
                assert ( this->my_rank_ == vec.my_rank ( ) );
                assert ( this->ownership_begin_ == vec.ownership_begin ( ) );
                assert ( this->ownership_end_ == vec.ownership_end ( ) );

                this->interior_->CopyTo ( vec.interior ( ) );
                this->ghost_->CopyTo ( vec.ghost ( ) );

                if ( this->HasPpData ( ) )
                {
                    if ( vec.pp_data_ != NULL )
                    {
                        *( vec.pp_data_ ) = this->pp_data ( );
                    }
                    else
                    {
                        vec.pp_data_ = this->pp_data ( ).Clone ( );
                    }
                }
                else
                {
                    if ( vec.pp_data_ != NULL )
                    {
                        delete vec.pp_data_;
                    }
                    vec.pp_data_ = NULL;
                }
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::CopyInteriorFrom ( const CoupledVector<DataType>& vec )
        {
            if ( this != &vec )
            {
                assert ( this->nb_procs_ == vec.nb_procs ( ) );
                assert ( this->my_rank_ == vec.my_rank ( ) );
                assert ( this->ownership_begin_ == vec.ownership_begin ( ) );
                assert ( this->ownership_end_ == vec.ownership_end ( ) );

                this->interior_->CopyFrom ( vec.interior ( ) );
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::CloneFromWithoutContent
        ( const CoupledVector<DataType>& vec )
        {
            if ( this != &vec )
            {
                this->Clear ( );
                int info = 0;
                if ( this->comm_ != MPI_COMM_NULL )
                {
                    info = MPI_Comm_free ( &this->comm_ );
                    assert ( info == MPI_SUCCESS );
                }
                this->nb_procs_ = vec.nb_procs ( );
                this->my_rank_ = vec.my_rank ( );

                //info = MPI_Comm_split ( vec.comm ( ), 0, this->my_rank_, &( this->comm_ ) );
                info = MPI_Comm_dup ( vec.comm ( ), &( this->comm_ ) );

                assert ( info == MPI_SUCCESS );
                this->la_couplings_ = &( vec.la_couplings ( ) );
                assert ( this->la_couplings_ != NULL );
                assert ( this->la_couplings_->initialized ( ) );
                this->ownership_begin_ = vec.ownership_begin ( );
                this->ownership_end_ = vec.ownership_end ( );
                this->nb_sends_ = vec.nb_sends ( );
                this->nb_recvs_ = vec.nb_recvs ( );
                this->mpi_req_ = vec.mpi_req ( );
                this->mpi_stat_ = vec.mpi_stat ( );
                this->border_val_ = vec.border_val ( );
                this->ghost_val_ = vec.ghost_val ( );

                // clone the vectors
                this->interior_ = vec.interior ( ).CloneWithoutContent ( );
                this->ghost_ = vec.ghost ( ).CloneWithoutContent ( );

                // copy the structure
                this->interior_->CopyStructureFrom ( vec.interior ( ) );
                this->ghost_->CopyStructureFrom ( vec.ghost ( ) );

                this->interior_->set_indexset ( this->la_couplings_->border_indices ( ),
                                                this->la_couplings_->size_border_indices ( ) );
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::CopyStructureFrom
        ( const CoupledVector<DataType>& vec )
        {
            if ( this != &vec )
            {

                // no Clear() !
                this->nb_procs_ = vec.nb_procs ( );
                this->my_rank_ = vec.my_rank ( );
                int info = 0;
                if ( this->comm_ != MPI_COMM_NULL )
                {
                    info = MPI_Comm_free ( &this->comm_ );
                    assert ( info == MPI_SUCCESS );
                }
                //info = MPI_Comm_split ( vec.comm ( ), 0, this->my_rank_, &( this->comm_ ) );
                info = MPI_Comm_dup ( vec.comm ( ), &( this->comm_ ) );

                assert ( info == MPI_SUCCESS );
                this->la_couplings_ = &( vec.la_couplings ( ) );
                assert ( this->la_couplings_ != NULL );
                assert ( this->la_couplings_->initialized ( ) );
                this->ownership_begin_ = vec.ownership_begin ( );
                this->ownership_end_ = vec.ownership_end ( );
                this->nb_sends_ = vec.nb_sends ( );
                this->nb_recvs_ = vec.nb_recvs ( );
                this->mpi_req_ = vec.mpi_req ( );
                this->mpi_stat_ = vec.mpi_stat ( );
                this->border_val_ = vec.border_val ( );
                this->ghost_val_ = vec.ghost_val ( );
                this->interior_->CopyStructureFrom ( vec.interior ( ) );
                this->ghost_->CopyStructureFrom ( vec.ghost ( ) );

                this->interior_->set_indexset ( this->la_couplings_->border_indices ( ),
                                                this->la_couplings_->size_border_indices ( ) );
                if ( this->pp_data_ != NULL )
                {
                    delete this->pp_data_;
                }
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::Clear ( )
        {
            // clear interior
            if ( this->interior_ != NULL )
            {
                delete this->interior_;
            }
            this->interior_ = NULL;

            // clear ghost
            if ( this->ghost_ != NULL )
            {
                delete this->ghost_;
            }
            this->ghost_ = NULL;

            // clear post processing data
            if ( this->pp_data_ != NULL )
            {
                delete this->pp_data_;
            }
            this->pp_data_ = NULL;

            this->mpi_req_.clear ( );
            this->mpi_stat_.clear ( );
            this->border_val_.clear ( );
            this->ghost_val_.clear ( );
            this->checked_for_dof_partition_ = false;
        }

        template<class DataType>
        void CoupledVector<DataType>::SendBorder ( )
        {
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );
            assert ( this->interior_ != NULL );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );

            this->interior_->GetIndexedValues ( &( this->border_val_[0] ) );

            int tag = 1;
            int ctr = 0;

            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( this->la_couplings_->border_offsets ( id + 1 ) -
                     this->la_couplings_->border_offsets ( id ) > 0 )
                {
                    int info = MPI_Isend ( &( this->border_val_[0] ) + this->la_couplings_->border_offsets ( id ),
                                           this->la_couplings_->border_offsets ( id + 1 ) -
                                           this->la_couplings_->border_offsets ( id ),
                                           mpi_data_type<DataType>::get_type ( ), id, tag, this->comm_,
                                           &( this->mpi_req_[this->nb_recvs_ + ctr] ) );
                    assert ( info == MPI_SUCCESS );
                    ctr++;
                }
            }
            assert ( ctr == this->nb_sends_ );
        }

        template<class DataType>
        void CoupledVector<DataType>::ReceiveGhost ( )
        {
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );

            int tag = 1;
            int ctr = 0;

            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->la_couplings_->ghost_offsets ( id + 1 ) -
                     this->la_couplings_->ghost_offsets ( id ) > 0 )
                {
                    int info = MPI_Irecv ( &( this->ghost_val_[0] ) + this->la_couplings_->ghost_offsets ( id ),
                                           this->la_couplings_->ghost_offsets ( id + 1 ) -
                                           this->la_couplings_->ghost_offsets ( id ),
                                           mpi_data_type<DataType>::get_type ( ), id, tag, this->comm_,
                                           &( this->mpi_req_[ctr] ) );
                    assert ( info == MPI_SUCCESS );
                    ctr++;
                }
            }
            assert ( ctr == this->nb_recvs_ );
        }

        template<class DataType>
        void CoupledVector<DataType>::WaitForSend ( )
        {
            int info = MPI_Waitall ( this->nb_sends_,
                                     &( this->mpi_req_[this->nb_recvs_] ),
                                     &( this->mpi_stat_[this->nb_recvs_] ) );
            assert ( info == MPI_SUCCESS );
        }

        template<class DataType>
        void CoupledVector<DataType>::WaitForRecv ( )
        {
            int info = MPI_Waitall ( this->nb_recvs_,
                                     &( this->mpi_req_[0] ),
                                     &( this->mpi_stat_[0] ) );
            assert ( info == MPI_SUCCESS );

            this->SetGhostValues ( &( this->ghost_val_[0] ) );
        }

        template<class DataType>
        void CoupledVector<DataType>::begin_update ( )
        {
            assert ( this->ghost_ != NULL );

            this->ReceiveGhost ( );
            this->SendBorder ( );
        }

        template<class DataType>
        void CoupledVector<DataType>::end_update ( )
        {
            assert ( this->ghost_ != NULL );

            this->WaitForRecv ( );
            this->WaitForSend ( );

            // update post processing if possible.
            if ( this->pp_data_ != NULL )
            {
                this->UpdatePpValues ( );
            }
            else if ( !checked_for_dof_partition_ )
            {
                const Couplings<DataType>* couplings = dynamic_cast < const Couplings<DataType>* > ( la_couplings_ );
                if ( couplings != NULL )
                {
                    this->InitializePostProcessing ( couplings );
                    this->UpdatePpValues ( );
                }
                else
                {
                    exit ( -1 );
                }
                checked_for_dof_partition_ = true;
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::UpdateCouplings ( )
        {

            this->UpdateGhost ( );

            // update post processing if possible.
            if ( this->pp_data_ != NULL )
            {
                this->UpdatePpValues ( );
            }
            else if ( !checked_for_dof_partition_ )
            {
                const Couplings<DataType>* couplings = dynamic_cast < const Couplings<DataType>* > ( la_couplings_ );
                if ( couplings != NULL )
                {
                    this->InitializePostProcessing ( couplings );
                    this->UpdatePpValues ( );
                }
                else
                {
                    exit ( -1 );
                }
                checked_for_dof_partition_ = true;
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::UpdateGhost ( )
        {
            assert ( this->ghost_ != NULL );

            this->ReceiveGhost ( );
            this->SendBorder ( );

            this->WaitForRecv ( );
            this->WaitForSend ( );
        }

        template<class DataType>
        int CoupledVector<DataType>::size_global ( ) const
        {
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );
            return this->la_couplings_->nb_total_dofs ( );
        }

        template<class DataType>
        void CoupledVector<DataType>::Print ( std::ostream &out ) const
        {
            this->interior_->print ( out );
            this->ghost_->print ( out );
        }

        template<class DataType>
        void CoupledVector<DataType>::ComputeOwnershipRange ( )
        {
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );
            assert ( this->my_rank_ >= 0 );

            this->ownership_begin_ = this->la_couplings_->dof_offset ( this->my_rank_ );
            this->ownership_end_ = this->ownership_begin_ +
                    this->la_couplings_->nb_dofs ( this->my_rank_ );
        }

        template<class DataType>
        void CoupledVector<DataType>::WriteHDF5 ( const std::string& filename,
                                                  const std::string& groupname,
                                                  const std::string& datasetname )
        {
#ifdef WITH_HDF5
            DataType* data;
            // Define Data in memory
            data = ( DataType * ) malloc ( sizeof (DataType ) * ( this->size_local ( ) ) );
            this->interior_->GetBlockValues ( 0, size_local ( ), data );

            H5FilePtr file_ptr ( new H5File ( filename, "w", this->comm_ ) );
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname, "w" ) );
            H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, this->size_global ( ),
                                                       datasetname, "w", data ) );
            dataset_ptr->write ( this->size_local ( ), this->ownership_begin_, data );

            free ( data );
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif

        }

        template<class DataType>
        void CoupledVector<DataType>::ReadHDF5 ( const std::string& filename,
                                                 const std::string& groupname,
                                                 const std::string& datasetname )
        {
#ifdef WITH_HDF5
            DataType* buffer;
            buffer = ( DataType * ) malloc ( sizeof (DataType ) * ( this->size_local ( ) ) );

            H5FilePtr file_ptr ( new H5File ( filename, "r", this->comm_ ) );
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname, "r" ) );
            H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, this->size_global ( ),
                                                       datasetname, "r", buffer ) );
            dataset_ptr->read ( this->size_local ( ), this->ownership_begin_, buffer );

            std::vector<int> ind ( this->size_local ( ) );
            for ( int i = 0; i < size_local ( ); ++i ) ind[i] = i;
            this->interior_->SetValues ( vec2ptr ( ind ), this->size_local ( ), buffer );

            // Update
            this->UpdateGhost ( );

            free ( buffer );
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        template<class DataType>
        void CoupledVector<DataType>::InitializePostProcessing ( const Couplings<DataType>* couplings )
        {
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );

            if ( this->pp_data_ != NULL )
            {
                delete this->pp_data_;
            }
            this->pp_data_ = new PpData<DataType>;
            // dof_ids and offsets are only temporarily needed since
            // we have most of the DoFs in interior.
            std::vector<int> dof_ids;
            std::vector<int> offsets;

            const typename doffem::DofPartition<DataType>* dp = couplings->dof_partition ( );
            assert ( dp != NULL );
            dp->GetPostProcessingStructure ( dof_ids, offsets );

            LOG_DEBUG ( 2, "[" << dp->get_my_subdomain ( ) << "] PostProcessingStructure: dof_ids "
                        << string_from_range ( dof_ids.begin ( ), dof_ids.end ( ) ) );

            LOG_DEBUG ( 2, "[" << dp->get_my_subdomain ( ) << "] PostProcessingStructure: offsets "
                        << string_from_range ( offsets.begin ( ), offsets.end ( ) ) );

            int tag = 1;
            std::vector<MPI_Request> mpi_req_send ( this->nb_procs_ );
            std::vector<MPI_Request> mpi_req_recv ( this->nb_procs_ );

            // How much is to be sent to other processes
            this->pp_data_->nb_sends.clear ( );
            this->pp_data_->nb_sends.resize ( this->nb_procs_, 0 );

            // How much do I request from others
            this->pp_data_->nb_recvs.clear ( );
            this->pp_data_->nb_recvs.resize ( this->nb_procs_, 0 );

            this->pp_data_->dof_ids.clear ( );
            this->pp_data_->dof_ids_offsets.resize ( this->nb_procs_ + 1, 0 );

            this->pp_data_->dof_ids.clear ( );

            // Extract the Dofs you additionally need besides the ones from
            // interior and ghost
            int counter = 0;
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( id != this->my_rank_ )
                {
                    for ( int k = offsets[id], e_k = offsets[id + 1]; k < e_k; ++k )
                    {
                        if ( this->la_couplings_->global2offdiag ( ).find ( dof_ids[k] ) ==
                             this->la_couplings_->global2offdiag ( ).end ( ) )
                        {
                            this->pp_data_->dof_ids.push_back ( dof_ids[k] );
                            ++counter;
                        }
                    }
                }
                this->pp_data_->dof_ids_offsets[id + 1] = counter;
            }
            this->pp_data_->values.resize ( this->pp_data_->dof_ids.size ( ), 0. );

            int info;
            // Exchange sizes
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( id != this->my_rank_ )
                {
                    this->pp_data_->nb_recvs[id] = this->pp_data_->dof_ids_offsets[id + 1] - this->pp_data_->dof_ids_offsets[id];
                    info = MPI_Isend ( &( this->pp_data_->nb_recvs[id] ), 1, MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );

                    info = MPI_Irecv ( &( this->pp_data_->nb_sends[id] ), 1, MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }
            // Make sure all data has been sent and received
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

            int nb_send_total = 0;

            for ( int i = 0; i < this->nb_procs_; ++i )
                nb_send_total += this->pp_data_->nb_sends[i];

            int nb_recv_total = 0;

            for ( int i = 0; i < this->nb_procs_; ++i )
                nb_recv_total += this->pp_data_->nb_recvs[i];

            this->pp_data_->sorted_dof_ids.clear ( );
            this->pp_data_->sorted_dof_ids = SortedArray<int>( this->pp_data_->dof_ids.begin ( ), this->pp_data_->dof_ids.end ( ) );

            // Allocate memory for the DoFs to be send.
            // Remark: We don't save the send values.
            this->pp_data_->send_dofs.clear ( );
            this->pp_data_->send_dofs.resize ( nb_send_total, 0 );
            this->pp_data_->send_offsets.clear ( );
            this->pp_data_->send_offsets.resize ( this->nb_procs_ + 1, 0 );

            for ( int i = 1; i < this->nb_procs_ + 1; ++i )
                this->pp_data_->send_offsets[i] = this->pp_data_->send_offsets[i - 1]
                    + this->pp_data_->nb_sends[i - 1];

            // Exchange DofIDs.
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->pp_data_->nb_recvs[id] > 0 )
                {
                    info = MPI_Isend ( &( this->pp_data_->sorted_dof_ids[this->pp_data_->dof_ids_offsets[id]] ),
                                       this->pp_data_->nb_recvs[id], MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );
                }

                if ( this->pp_data_->nb_sends[id] > 0 )
                {
                    info = MPI_Irecv ( &( this->pp_data_->send_dofs[this->pp_data_->send_offsets[id]] ),
                                       this->pp_data_->nb_sends[id], MPI_INT,
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // Make sure all data has been sent and received
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->pp_data_->nb_recvs[id] > 0 )
                {
                    info = MPI_Wait ( &( mpi_req_recv[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
                if ( this->pp_data_->nb_sends[id] > 0 )
                {
                    info = MPI_Wait ( &( mpi_req_send[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
            }
        }

        template<class DataType>
        bool CoupledVector<DataType>::HasPpData ( ) const
        {
            return (this->pp_data_ != NULL );
        }

        template<class DataType>
        void CoupledVector<DataType>::UpdatePpValues ( )
        {
            assert ( interior_ != NULL );
            assert ( pp_data_ != NULL );
            int tag = 1;
            int info;
            std::vector<MPI_Request> mpi_req_send ( nb_procs_ );
            std::vector<MPI_Request> mpi_req_recv ( nb_procs_ );

            //prepare the values for sending
            std::vector<DataType> send_pp_values ( this->pp_data_->send_dofs.size ( ) );
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->pp_data_->nb_sends[id] > 0 )
                {
                    this->GetValues ( &( this->pp_data_->send_dofs[this->pp_data_->send_offsets[id]] ),
                                      this->pp_data_->nb_sends[id],
                                      &( send_pp_values[this->pp_data_->send_offsets[id]] ) );
                }
            }

            // exchange values
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->pp_data_->nb_sends[id] > 0 )
                {
                    info = MPI_Isend ( &( send_pp_values[this->pp_data_->send_offsets[id]] ),
                                       this->pp_data_->nb_sends[id], mpi_data_type<DataType>::get_type ( ),
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
                if ( this->pp_data_->nb_recvs[id] > 0 )
                {
                    info = MPI_Irecv ( &( this->pp_data_->values[this->pp_data_->dof_ids_offsets[id]] ),
                                       this->pp_data_->nb_recvs[id], mpi_data_type<DataType>::get_type ( ),
                                       id, tag, this->comm_, &( mpi_req_send[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
            }

            // make sure all data has been sent and received
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->pp_data_->nb_sends[id] > 0 )
                {
                    info = MPI_Wait ( &( mpi_req_recv[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
                if ( this->pp_data_->nb_recvs[id] > 0 )
                {
                    info = MPI_Wait ( &( mpi_req_send[id] ), MPI_STATUS_IGNORE );
                    assert ( info == MPI_SUCCESS );
                }
            }
        }

        template<class DataType>
        void CoupledVector<DataType>::GetAllDofsAndValues ( std::vector<int>& id, std::vector<DataType>& val ) const
        {
            assert ( this->la_couplings_ != NULL );
            assert ( this->la_couplings_->initialized ( ) );
            assert ( this->interior_ != NULL );
            assert ( this->ghost_ != NULL );

            // Temporary containers for interior and ghost values
            DataType* values = new DataType[this->size_local ( )];
            DataType* ghost_values = new DataType[this->size_local_ghost ( )];

            // Combine interior, ghost and pp_data values
            this->interior_->GetBlockValues ( 0, this->size_local ( ), values );
            this->ghost_->GetBlockValues ( 0, this->size_local_ghost ( ), ghost_values );

            int total_size = this->size_local ( ) + this->size_local_ghost ( );

            if ( this->pp_data_ != NULL )
            {
                total_size += this->pp_data_->values.size ( );
            }

            id.resize ( total_size );
            val.resize ( total_size );
            // First: the DoFs and values from the interior
            for ( int i = 0; i < this->size_local ( ); ++i )
            {
                id[i] = this->ownership_begin_ + i;
                val[i] = values[i];
            }
            delete[] values;
            // Second: the DoFs and values from ghost
            int tmp_offset = this->size_local ( );
            for ( int i = 0; i < this->size_local_ghost ( ); ++i )
            {
                id[i + tmp_offset] = this->la_couplings_->Offdiag2Global ( i );
                val[i + tmp_offset] = ghost_values[i];
            }
            delete[] ghost_values;
            // Last: the DoFs and values additionaly needed for post processing
            if ( this->pp_data_ != NULL )
            {
                tmp_offset = this->size_local ( ) + this->size_local_ghost ( );
                for ( int i = 0, e_i = this->pp_data_->dof_ids_offsets[this->nb_procs_]; i < e_i; ++i )
                {
                    id[i + tmp_offset] = this->pp_data_->sorted_dof_ids[i];
                    val[i + tmp_offset] = this->pp_data_->values[i];
                }
            }
        }

        /// template instantiation
        template class CoupledVector<double>;
        template class CoupledVector<float>;

    } // namespace la
} // namespace hiflow
