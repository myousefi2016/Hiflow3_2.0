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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-11-17

#include "common/log.h"
#include "common/pointers.h"
#include "linear_algebra/petsc_vector.h"
#include "lmp/init_vec_mat.h"
#include "petsc.h"
#include "petsc_environment.h"
#include "petsc_vector_interface.h"
#include <vector>
#include <cmath>
#include <cstdlib>

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace la
    {

        template <class DataType>
        PETScVector<DataType>::PETScVector ( )
        : comm_ ( MPI_COMM_NULL ),
        ilower_ ( -1 ),
        iupper_ ( -1 ),
        cp_ ( NULL ),
        ghost_ ( init_vector<DataType>( 0, "ghost", CPU, NAIVE ) ),
        nb_procs_ ( -1 ),
        my_rank_ ( -1 ),
        nb_sends_ ( -1 ),
        nb_recvs_ ( -1 ),
        mpi_req_ ( ),
        mpi_stat_ ( ),
        border_val_ ( ),
        ghost_val_ ( ),
        pp_data_ ( NULL ),
        border_indices_ ( ),
        checked_for_dof_partition_ ( false ),
        initialized_ ( false ),
        ptr_vec_wrapper_ ( new petsc::Vec_wrapper )
        {
        }

        template <class DataType>
        PETScVector<DataType>::~PETScVector ( )
        {
            this->Clear ( );

            // clear ghost
            if ( this->ghost_ != NULL )
            {
                this->ghost_->Clear ( );
                delete this->ghost_;
            }
            this->ghost_ = NULL;
        }

        template <class DataType>
        Vector<DataType>* PETScVector<DataType>::Clone ( ) const
        {
            LOG_ERROR ( "Called PETScVector::Clone not yet implemented!!!" );
            exit ( -1 );
            return NULL;
        }

        template <class DataType>
        void PETScVector<DataType>::Init ( MPI_Comm& comm )
        {
            if ( this->comm_ != MPI_COMM_NULL ) MPI_Comm_free ( &this->comm_ );
            assert ( comm != MPI_COMM_NULL );

            // determine nb. of processes
            int info = MPI_Comm_size ( comm, &nb_procs_ );
            assert ( info == MPI_SUCCESS );
            assert ( nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( comm, &my_rank_ );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank_ >= 0 );
            assert ( my_rank_ < nb_procs_ );

            info = MPI_Comm_split ( comm, 0, my_rank_, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );
        }

        template <class DataType>
        void PETScVector<DataType>::Init ( MPI_Comm& comm, const LaCouplings& cp )
        {
            // clear possibly existing DataType
            if ( initialized_ ) Clear ( );

            Init ( comm );
            cp_ = &cp;

            // Compute indices range of this process
            ilower_ = cp.dof_offset ( my_rank_ );
            iupper_ = ilower_ + cp.nb_dofs ( my_rank_ ) - 1;

            // Prepare PETSc MPI interface
            PETScEnvironment::initialize ( );

            // Create PETSC Vector
            VecCreateMPI ( comm_, cp_->nb_dofs ( my_rank_ ), cp_->nb_total_dofs ( ), &ptr_vec_wrapper_->vec_ );

            // Set initialized_ flag
            initialized_ = true;
        }

        template <class DataType>
        void PETScVector<DataType>::Clear ( )
        {
            if ( initialized_ )
            {
                VecDestroy ( &ptr_vec_wrapper_->vec_ );

                // clear ghost
                if ( this->ghost_ != NULL )
                {
                    this->ghost_->Clear ( );
                    delete this->ghost_;
                }
                this->ghost_ = NULL;

                // clear post processing data
                if ( this->pp_data_ != NULL )
                {
                    delete this->pp_data_;
                }
                this->pp_data_ = NULL;

                this->cp_ = NULL;

                mpi_req_.clear ( );
                mpi_stat_.clear ( );
                border_val_.clear ( );
                ghost_val_.clear ( );
                border_indices_.clear ( );
                this->checked_for_dof_partition_ = false;
                ghost_ = init_vector<DataType>( 0, "ghost", CPU, NAIVE );
                assert ( this->ghost_ != NULL );
                initialized_ = false;
            }
        }

        template <class DataType>
        void PETScVector<DataType>::CloneFromWithoutContent (
                                                              const PETScVector<DataType>& vec )
        {
        }

        template <class DataType>
        void PETScVector<DataType>::CopyFrom ( const PETScVector<DataType>& vec )
        {
        }

        template <class DataType>
        int PETScVector<DataType>::size_local ( ) const
        {
            int size;
            CHKERRQ ( VecGetLocalSize ( ptr_vec_wrapper_->vec_, &size ) );
            return size;
        }

        template <class DataType>
        int PETScVector<DataType>::size_global ( ) const
        {
            int size;
            VecGetSize ( ptr_vec_wrapper_->vec_, &size );
            return size;
        }

        template <class DataType>
        void PETScVector<DataType>::Zeros ( )
        {
            VecSet ( ptr_vec_wrapper_->vec_, 0 );
        }

        template <class DataType>
        DataType PETScVector<DataType>::GetValue ( int index ) const
        {
            assert ( initialized_ );
            DataType value;
            VecGetValues ( ptr_vec_wrapper_->vec_, 1, &index, &value );
            return value;
        }

        template <class DataType>
        void PETScVector<DataType>::GetValues ( const int* indices,
                                                int length,
                                                DataType* values ) const
        {
            assert ( initialized_ );
            if ( length <= 0 ) return;
            assert ( indices != NULL );
            assert ( values != NULL );

            // extract interior indices
            std::vector<int> indices_interior;
            indices_interior.reserve ( length );

            // extract ghost indices
            std::vector<int> indices_ghost;
            indices_interior.reserve ( length );

            // extract pp indices
            std::vector<int> indices_pp;
            indices_pp.reserve ( length );

            // transform indices according to local numbering of interior and ghost
            for ( int i = 0; i < length; i++ )
            {
                // interior
                if ( this->ilower_ <= indices[i] && indices[i] <= this->iupper_ )
                {
                    indices_interior.push_back ( indices[i] );
                }

                    // ghost
                else if ( this->cp_->global2offdiag ( ).find ( indices[i] ) !=
                          this->cp_->global2offdiag ( ).end ( ) )
                {
                    indices_ghost.push_back ( this->cp_->Global2Offdiag ( indices[i] ) );
                }

                    // pp values
                else if ( this->pp_data_ )
                {
                    int pos;
                    this->pp_data_->sorted_dof_ids.find ( indices[i], &pos );
                    indices_pp.push_back ( pos );
                }
            }

            // extract values
            std::vector<DataType> values_interior ( indices_interior.size ( ) );
            if ( !indices_interior.empty ( ) )
            {
                VecSetValues ( ptr_vec_wrapper_->vec_, indices_interior.size ( ), vec2ptr ( indices_interior ), vec2ptr ( values_interior ), INSERT_VALUES );
                VecAssemblyBegin ( ptr_vec_wrapper_->vec_ );
                VecAssemblyEnd ( ptr_vec_wrapper_->vec_ );
            }

            std::vector<DataType> values_ghost ( indices_ghost.size ( ) );
            if ( !indices_ghost.empty ( ) )
            {
                this->ghost_->GetValues ( vec2ptr ( indices_ghost ), indices_ghost.size ( ), vec2ptr ( values_ghost ) );
            }

            std::vector<DataType> values_pp ( indices_pp.size ( ) );
            if ( this->pp_data_ )
            {
                for ( int i = 0; i < indices_pp.size ( ); ++i )
                {
                    values_pp[i] = this->pp_data_->values[indices_pp[i]];
                }
            }

            typename std::vector<DataType>::const_iterator iter_values_interior_cur = values_interior.begin ( );
            typename std::vector<DataType>::const_iterator iter_values_ghost_cur = values_ghost.begin ( );
            typename std::vector<DataType>::const_iterator iter_values_pp_cur = values_pp.begin ( );

            // put values together
            for ( int i = 0; i < length; i++ )
            {
                int cont;
                // interior
                if ( this->ilower_ <= indices[i] and indices[i] <= this->iupper_ )
                {
                    values[i] = *iter_values_interior_cur++;
                }
                    // ghost
                else if ( this->cp_->global2offdiag ( ).find ( indices[length - i - 1] ) !=
                          this->cp_->global2offdiag ( ).end ( ) )
                {
                    values[i] = *iter_values_ghost_cur++;
                }
                    // pp
                else if ( ( this->pp_data_ != NULL ) &&
                          ( this->pp_data_->sorted_dof_ids.find ( indices[length - i - 1], &cont ) ) )
                {
                    values[i] = *iter_values_pp_cur++;
                }
                    // zeros
                else
                {
                    values[i] = 0;
                }
            }
        }

        template <class DataType>
        DataType PETScVector<DataType>::Norm2 ( ) const
        {
            DataType value;
            VecNorm ( ptr_vec_wrapper_->vec_, NORM_2, &value );
            return value;
        }

        template <class DataType>
        DataType PETScVector<DataType>::Norm1 ( ) const
        {
            DataType value;
            VecNorm ( ptr_vec_wrapper_->vec_, NORM_1, &value );
            return value;
        }

        template <class DataType>
        DataType PETScVector<DataType>::NormMax ( ) const
        {
            DataType value;
            VecNorm ( ptr_vec_wrapper_->vec_, NORM_INFINITY, &value );
            return value;
        }

        template <class DataType>
        DataType PETScVector<DataType>::Dot ( const Vector<DataType>& vec ) const
        {
            const PETScVector<DataType>* hv = dynamic_cast < const PETScVector<DataType>* > ( &vec );
            if ( !hv )
            {
                LOG_ERROR ( "Called PETScVector::Dot with incompatible vector type." );
                exit ( -1 );
            }
            return this->Dot ( *hv );
        }

        template <class DataType>
        DataType PETScVector<DataType>::Dot ( const PETScVector<DataType>& vec ) const
        {
            DataType value;
            VecDot ( ptr_vec_wrapper_->vec_, vec.ptr_vec_wrapper_->vec_, &value );
            return value;
        }

        template <class DataType>
        void PETScVector<DataType>::Add ( int index, DataType value )
        {
            VecSetValues ( ptr_vec_wrapper_->vec_, 1, &index, &value, ADD_VALUES );
            VecAssemblyBegin ( ptr_vec_wrapper_->vec_ );
            VecAssemblyEnd ( ptr_vec_wrapper_->vec_ );
        }

        template <class DataType>
        void PETScVector<DataType>::Add ( const int* indices,
                                          int length,
                                          const DataType* values )
        {
            VecSetValues ( ptr_vec_wrapper_->vec_, length, indices, values, ADD_VALUES );
            VecAssemblyBegin ( ptr_vec_wrapper_->vec_ );
            VecAssemblyEnd ( ptr_vec_wrapper_->vec_ );
        }

        template <class DataType>
        void PETScVector<DataType>::SetValue ( int index, DataType value )
        {
            VecSetValues ( ptr_vec_wrapper_->vec_, 1, &index, &value, INSERT_VALUES );
            VecAssemblyBegin ( ptr_vec_wrapper_->vec_ );
            VecAssemblyEnd ( ptr_vec_wrapper_->vec_ );
        }

        template <class DataType>
        void PETScVector<DataType>::SetValues ( const int* indices,
                                                const int length,
                                                const DataType* values )
        {
            VecSetValues ( ptr_vec_wrapper_->vec_, length, indices, values, INSERT_VALUES );
            VecAssemblyBegin ( ptr_vec_wrapper_->vec_ );
            VecAssemblyEnd ( ptr_vec_wrapper_->vec_ );
        }

        template <class DataType>
        void PETScVector<DataType>::Axpy ( const Vector<DataType>& vecx,
                                           DataType alpha )
        {
            const PETScVector<DataType>* hv = dynamic_cast < const PETScVector<DataType>* > ( &vecx );
            if ( !hv )
            {
                LOG_ERROR ( "Called PETScVector::Axpy with incompatible vector type." );
                exit ( -1 );
            }
            this->Axpy ( *hv, alpha );
        }

        template <class DataType>
        void PETScVector<DataType>::Axpy ( const PETScVector<DataType>& vecx,
                                           DataType alpha )
        {
            VecAXPY ( ptr_vec_wrapper_->vec_, alpha, vecx.ptr_vec_wrapper_->vec_ );
        }

        template <class DataType>
        void PETScVector<DataType>::ScaleAdd ( const Vector<DataType>& vecx,
                                               DataType alpha )
        {
            const PETScVector<DataType>* hv = dynamic_cast < const PETScVector<DataType>* > ( &vecx );
            if ( !hv )
            {
                LOG_ERROR ( "Called PETScVector::Axpy with incompatible vector type." );
                exit ( -1 );
            }
            this->ScaleAdd ( *hv, alpha );
        }

        template <class DataType>
        void PETScVector<DataType>::ScaleAdd ( const PETScVector<DataType>& vecx,
                                               DataType alpha )
        {
            VecScale ( ptr_vec_wrapper_->vec_, alpha );
            VecAXPY ( ptr_vec_wrapper_->vec_, 1, vecx.ptr_vec_wrapper_->vec_ );
        }

        template <class DataType>
        void PETScVector<DataType>::Scale ( DataType alpha )
        {
            VecScale ( ptr_vec_wrapper_->vec_, alpha );
        }

        template <class DataType>
        void PETScVector<DataType>::SendBorder ( )
        {
            assert ( this->cp_->initialized ( ) );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );

            this->GetValues ( vec2ptr ( border_indices_ ), border_indices_.size ( ),
                              vec2ptr ( border_val_ ) );

            int tag = 1;
            int ctr = 0;

            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( this->cp_->border_offsets ( id + 1 ) - this->cp_->border_offsets ( id ) > 0 )
                {
                    int info = MPI_Isend (
                                           &( this->border_val_[0] ) + this->cp_->border_offsets ( id ),
                                           this->cp_->border_offsets ( id + 1 ) - this->cp_->border_offsets ( id ),
                                           mpi_data_type<DataType>::get_type ( ), id, tag, this->comm_,
                                           &( this->mpi_req_[this->nb_recvs_ + ctr] ) );
                    assert ( info == MPI_SUCCESS );
                    ctr++;
                }
            }
            assert ( ctr == this->nb_sends_ );
        }

        template <class DataType>
        void PETScVector<DataType>::ReceiveGhost ( )
        {
            assert ( this->cp_->initialized ( ) );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );

            int tag = 1;
            int ctr = 0;

            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->cp_->ghost_offsets ( id + 1 ) - this->cp_->ghost_offsets ( id ) > 0 )
                {
                    int info = MPI_Irecv (
                                           &( this->ghost_val_[0] ) + this->cp_->ghost_offsets ( id ),
                                           this->cp_->ghost_offsets ( id + 1 ) - this->cp_->ghost_offsets ( id ),
                                           mpi_data_type<DataType>::get_type ( ), id, tag, this->comm_,
                                           &( this->mpi_req_[ctr] ) );
                    assert ( info == MPI_SUCCESS );
                    ctr++;
                }
            }
            assert ( ctr == this->nb_recvs_ );
        }

        template <class DataType>
        void PETScVector<DataType>::WaitForSend ( )
        {
            int info = MPI_Waitall ( this->nb_sends_, &( this->mpi_req_[this->nb_recvs_] ),
                                     &( this->mpi_stat_[this->nb_recvs_] ) );
            assert ( info == MPI_SUCCESS );
        }

        template <class DataType>
        void PETScVector<DataType>::WaitForRecv ( )
        {
            int info =
                    MPI_Waitall ( this->nb_recvs_, &( this->mpi_req_[0] ), &( this->mpi_stat_[0] ) );
            assert ( info == MPI_SUCCESS );

            this->SetGhostValues ( &( this->ghost_val_[0] ) );
        }

        template <class DataType>
        void PETScVector<DataType>::Update ( )
        {
            this->UpdateGhost ( );

            // update post processing if possible.
            if ( this->pp_data_ != NULL )
            {
                this->UpdatePpValues ( );
            }
            else if ( !checked_for_dof_partition_ )
            {
                const Couplings<DataType>* couplings =
                        dynamic_cast < const Couplings<DataType>* > ( cp_ );
                if ( couplings != NULL )
                {
                    this->InitializePostProcessing ( couplings );
                    this->UpdatePpValues ( );
                }
                checked_for_dof_partition_ = true;
            }
        }

        template <class DataType>
        void PETScVector<DataType>::UpdateGhost ( )
        {
            assert ( this->ghost_ != NULL );

            this->ReceiveGhost ( );
            this->SendBorder ( );

            this->WaitForRecv ( );
            this->WaitForSend ( );
        }

        template <class DataType>
        void PETScVector<DataType>::InitializePostProcessing (
                                                               const Couplings<DataType>* couplings )
        {
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );

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
                        if ( this->cp_->global2offdiag ( ).find ( dof_ids[k] ) ==
                             this->cp_->global2offdiag ( ).end ( ) )
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

        template <class DataType>
        void PETScVector<DataType>::UpdatePpValues ( )
        {
            assert ( pp_data_ != NULL );
            int tag = 1;
            int info;
            std::vector<MPI_Request> mpi_req_send ( nb_procs_ );
            std::vector<MPI_Request> mpi_req_recv ( nb_procs_ );

            // prepare the values for sending
            std::vector<DataType> send_pp_values ( this->pp_data_->send_dofs.size ( ) );
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->pp_data_->nb_sends[id] > 0 )
                {
                    this->GetValues (
                                      &( this->pp_data_->send_dofs[this->pp_data_->send_offsets[id]] ),
                                      this->pp_data_->nb_sends[id],
                                      &( send_pp_values[this->pp_data_->send_offsets[id]] ) );
                }
            }

            // exchange values
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->pp_data_->nb_sends[id] > 0 )
                {
                    info = MPI_Isend ( &( send_pp_values[0] ) + this->pp_data_->send_offsets[id],
                                       this->pp_data_->nb_sends[id],
                                       mpi_data_type<DataType>::get_type ( ), id, tag,
                                       this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
                if ( this->pp_data_->nb_recvs[id] > 0 )
                {
                    info = MPI_Irecv (
                                       &( this->pp_data_->values[0] ) + this->pp_data_->dof_ids_offsets[id],
                                       this->pp_data_->nb_recvs[id], mpi_data_type<DataType>::get_type ( ), id,
                                       tag, this->comm_, &( mpi_req_send[id] ) );
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

        template <class DataType>
        void PETScVector<DataType>::SetGhostValues ( const DataType* values )
        {
            assert ( this->ghost_ != NULL );

            this->ghost_->SetBlockValues ( 0, this->size_local_ghost ( ), values );
        }

        template <class DataType>
        void PETScVector<DataType>::GetAllDofsAndValues ( std::vector<int>& id,
                                                          std::vector<DataType>& val )
        {
            assert ( this->ghost_ != NULL );

            // Temporary containers for interior and ghost values
            DataType* values = new DataType[this->size_local ( )];
            DataType* ghost_values = new DataType[this->size_local_ghost ( )];

            int total_size = this->size_local ( ) + this->size_local_ghost ( );
            if ( this->pp_data_ != NULL )
            {
                total_size += this->pp_data_->values.size ( );
            }
            id.resize ( total_size );

            // First: the DoFs from the interior
            for ( int i = 0; i < this->size_local ( ); ++i )
            {
                id[i] = ilower_ + i;
            }

            // Combine interior, ghost and pp_data values
            this->GetValues ( vec2ptr ( id ), this->size_local ( ), values );
            this->ghost_->GetBlockValues ( 0, this->size_local_ghost ( ), ghost_values );

            val.resize ( total_size );
            // First: values from the interior
            for ( int i = 0; i < this->size_local ( ); ++i )
            {
                val[i] = values[i];
            }
            delete[] values;
            // Second: the DoFs and values from ghost
            int tmp_offset = this->size_local ( );
            for ( int i = 0; i < this->size_local_ghost ( ); ++i )
            {
                id[i + tmp_offset] = this->cp_->Offdiag2Global ( i );
                val[i + tmp_offset] = ghost_values[i];
            }
            delete[] ghost_values;
            // Last: the DoFs and values additionaly needed for post processing
            if ( this->pp_data_ != NULL )
            {
                tmp_offset = this->size_local ( ) + this->size_local_ghost ( );
                for ( int i = 0, e_i = this->pp_data_->dof_ids_offsets[this->nb_procs_];
                      i < e_i; ++i )
                {
                    id[i + tmp_offset] = this->pp_data_->dof_ids[i];
                    val[i + tmp_offset] = this->pp_data_->values[i];
                }
            }
        }

        // template instantiation
        template class PETScVector<double>;

    } // namespace la
} // namespace hiflow
