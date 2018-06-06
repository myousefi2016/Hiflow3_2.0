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

/// @author Simon Gawlok

#include "linear_algebra/hypre_vector.h"
#include "lmp/init_vec_mat.h"
#include "common/pointers.h"
#include "common/log.h"
#include <vector>
#include <cmath>
#include <cstdlib>

#ifdef WITH_HYPRE
#    include "_hypre_parcsr_mv.h"
#endif

#ifdef WITH_HDF5
#    include "hdf5.h"
#    include "common/hdf5_tools.h"
#endif

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace la
    {

        template<class DataType>
        HypreVector<DataType>::HypreVector ( )
        {
            this->initialized_ = false;
            this->cp_ = NULL;
            this->comm_ = MPI_COMM_NULL;
            this->nb_procs_ = -1;
            this->my_rank_ = -1;
            this->ghost_ = NULL;
            this->ilower_ = -1;
            this->iupper_ = -1;
            this->pp_data_ = NULL;
            this->checked_for_dof_partition_ = false;
            this->ghost_ = init_vector<DataType>( 0, "ghost", CPU, NAIVE );
            assert ( this->ghost_ != NULL );
            this->global_indices_ = NULL;
        }

        template<class DataType>
        HypreVector<DataType>::~HypreVector ( )
        {
            this->Clear ( );

            // clear ghost
            if ( this->ghost_ != NULL )
            {
                this->ghost_->Clear ( );
                delete this->ghost_;
            }
            this->ghost_ = NULL;

            int is_finalized;
            MPI_Finalized ( &is_finalized );
            if ( !is_finalized )
            {
                if ( this->comm_ != MPI_COMM_NULL )
                {
                    MPI_Comm_free ( &this->comm_ );
                    assert ( this->comm_ == MPI_COMM_NULL );
                }

            }
            if ( this->global_indices_ != NULL )
                delete [] this->global_indices_;
            this->global_indices_ = NULL;
        }

        template<class DataType>
        Vector<DataType>* HypreVector<DataType>::Clone ( ) const
        {
            LOG_ERROR ( "Called HypreVector::Clone not yet implemented!!!" );
            exit ( -1 );
            return NULL;
        }

        template<class DataType>
        void HypreVector<DataType>::Init ( const MPI_Comm &comm )
        {
            // clear possibly existing DataType
            if ( initialized_ )
            {
                this->Clear ( );
            }

            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
                assert ( this->comm_ == MPI_COMM_NULL );
            }

            assert ( comm != MPI_COMM_NULL );

            MPI_Comm_dup ( comm, &this->comm_ );
            assert ( this->comm_ != MPI_COMM_NULL );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( this->comm_, &nb_procs_ );
            assert ( info == MPI_SUCCESS );
            assert ( nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( this->comm_, &my_rank_ );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank_ >= 0 );
            assert ( my_rank_ < nb_procs_ );
        }

        template<class DataType>
        void HypreVector<DataType>::Init ( const MPI_Comm &comm, const LaCouplings &cp )
        {
            // clear possibly existing DataType
            if ( initialized_ )
            {
                this->Clear ( );
            }

            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
                assert ( this->comm_ == MPI_COMM_NULL );
            }

            assert ( comm != MPI_COMM_NULL );

            MPI_Comm_dup ( comm, &this->comm_ );
            assert ( this->comm_ != MPI_COMM_NULL );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( this->comm_, &nb_procs_ );
            assert ( info == MPI_SUCCESS );
            assert ( nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( this->comm_, &my_rank_ );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank_ >= 0 );
            assert ( my_rank_ < nb_procs_ );
#ifdef WITH_HYPRE

            this->cp_ = &cp;
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );

            // Get rank of current process
            MPI_Comm_rank ( comm_, &my_rank_ );

            // Compute indices range of this process
            ilower_ = this->cp_->dof_offset ( my_rank_ );
            iupper_ = ilower_ + this->cp_->nb_dofs ( my_rank_ ) - 1;

            // Create HYPRE Vector
            HYPRE_IJVectorCreate ( comm_, ilower_, iupper_, &x_ );

            // Use parallel csr format
            HYPRE_IJVectorSetObjectType ( x_, HYPRE_PARCSR );

            HYPRE_IJVectorSetPrintLevel ( x_, 100 );

            // Tell HYPRE that no vector entries need to be communicated to other processors
            HYPRE_IJVectorSetMaxOffProcElmts ( x_, 0 );

            // Initialize
            HYPRE_IJVectorInitialize ( x_ );

            // Initialize exact structure of vector. To achieve this, we set every element to zero.
            const int local_size = iupper_ - ilower_ + 1;

            this->global_indices_ = new int[local_size ];

            std::vector<DataType> val ( local_size, static_cast < DataType > ( 0 ) );

            const int N = local_size;

            assert ( N > 0 );

            /* Version that uses loop unrolling by an unroll-factor of 5*/
            // compute overhead to unroll factor
            const int M = N % 5;

            // if N is a multiple of 5
            if ( M == 0 )
            {
#    pragma clang loop vectorize(enable)
                for ( int i = 0; i < N; i += 5 )
                {
                    this->global_indices_[i] = ilower_ + i;
                    this->global_indices_[i + 1] = ilower_ + i + 1;
                    this->global_indices_[i + 2] = ilower_ + i + 2;
                    this->global_indices_[i + 3] = ilower_ + i + 3;
                    this->global_indices_[i + 4] = ilower_ + i + 4;
                }
            }
            else
            {
                // result for overhead to unroll factor
#    pragma clang loop vectorize(enable)
                for ( int i = 0; i < M; ++i )
                {
                    this->global_indices_[i] = ilower_ + i;
                }

                // result for rest of vectors if length is greater than the unroll factor
                if ( N > 5 )
                {
#    pragma clang loop vectorize(enable)
                    for ( int i = M; i < N; i += 5 )
                    {
                        this->global_indices_[i] = ilower_ + i;
                        this->global_indices_[i + 1] = ilower_ + i + 1;
                        this->global_indices_[i + 2] = ilower_ + i + 2;
                        this->global_indices_[i + 3] = ilower_ + i + 3;
                        this->global_indices_[i + 4] = ilower_ + i + 4;
                    }
                }
            }

            int nvals = local_size;

            HYPRE_IJVectorSetValues ( x_, nvals, this->global_indices_, vec2ptr ( val ) );

            // Finalize initialization of vector
            HYPRE_IJVectorAssemble ( x_ );
            HYPRE_IJVectorGetObject ( x_, ( void ** ) &parcsr_x_ );

            // set border indices
            const size_t size_border = this->cp_->size_border_indices ( );
            border_indices_.resize ( size_border );
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < size_border; ++i )
            {
                this->border_indices_[i] = this->cp_->border_indices ( )[i] + ilower_;
            }

            // Initialize ghost part
            // init structure of ghost
            if ( this->cp_->size_ghost ( ) > 0 )
            {
                this->ghost_->Clear ( );
                this->ghost_->Init ( this->cp_->size_ghost ( ),
                                     "ghost" );
            }

            // prepare for communication
            this->nb_sends_ = 0;
            this->nb_recvs_ = 0;
            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->cp_->border_offsets ( id + 1 ) -
                     this->cp_->border_offsets ( id ) > 0 )
                {
                    this->nb_sends_++;
                }
                if ( this->cp_->ghost_offsets ( id + 1 ) -
                     this->cp_->ghost_offsets ( id ) > 0 )
                {
                    this->nb_recvs_++;
                }
            }
            this->mpi_req_.resize ( this->nb_sends_ + this->nb_recvs_ );
            this->mpi_stat_.resize ( this->nb_sends_ + this->nb_recvs_ );
            this->border_val_.resize ( this->cp_->size_border_indices ( ) );
            this->ghost_val_.resize ( this->cp_->size_ghost ( ) );

            // Set initialized_ flag
            initialized_ = true;

            val.clear ( );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreVector<DataType>::Clear ( )
        {
#ifdef WITH_HYPRE
            if ( initialized_ )
            {
                HYPRE_IJVectorDestroy ( x_ );
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

                mpi_req_.clear ( );
                mpi_stat_.clear ( );
                border_val_.clear ( );
                ghost_val_.clear ( );
                border_indices_.clear ( );
                this->checked_for_dof_partition_ = false;
                ghost_ = init_vector<DataType>( 0, "ghost", CPU, NAIVE );
                assert ( this->ghost_ != NULL );
            }
            initialized_ = false;
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreVector<DataType>::CloneFromWithoutContent ( const HypreVector<DataType> &vec )
        {
#ifdef WITH_HYPRE
            if ( this != &vec )
            {
                this->Clear ( );
                int info = 0;
                if ( this->comm_ != MPI_COMM_NULL )
                {
                    MPI_Comm_free ( &this->comm_ );
                    assert ( this->comm_ == MPI_COMM_NULL );
                }

                assert ( vec.comm ( ) != MPI_COMM_NULL );

                MPI_Comm_dup ( vec.comm ( ), &this->comm_ );
                assert ( this->comm_ != MPI_COMM_NULL );
                // MPI communicator

                // determine nb. of processes
                info = MPI_Comm_size ( this->comm_, &nb_procs_ );
                assert ( info == MPI_SUCCESS );
                assert ( nb_procs_ > 0 );

                // retrieve my rank
                info = MPI_Comm_rank ( this->comm_, &my_rank_ );
                assert ( info == MPI_SUCCESS );
                assert ( my_rank_ >= 0 );
                assert ( my_rank_ < nb_procs_ );

                this->cp_ = &( vec.la_couplings ( ) );
                assert ( this->cp_ != NULL );
                assert ( this->cp_->initialized ( ) );

                // Compute indices range of this process
                ilower_ = vec.la_couplings ( ).dof_offset ( my_rank_ );
                iupper_ = ilower_ + vec.la_couplings ( ).nb_dofs ( my_rank_ ) - 1;

                // Create HYPRE Vector
                HYPRE_IJVectorCreate ( comm_, ilower_, iupper_, &x_ );

                // Use parallel csr format
                HYPRE_IJVectorSetObjectType ( x_, HYPRE_PARCSR );

                // Tell HYPRE that no vector entries need to be communicated to other processors
                HYPRE_IJVectorSetMaxOffProcElmts ( x_, 0 );

                // Initialize
                HYPRE_IJVectorInitialize ( x_ );

                // Initialize exact structure of vector. To achieve this, we set every element to zero.
                const int local_size = iupper_ - ilower_ + 1;

                this->global_indices_ = new int[local_size ];

                std::vector<DataType> val ( local_size, static_cast < DataType > ( 0 ) );

                const int N = local_size;

                assert ( N > 0 );

                /* Version that uses loop unrolling by an unroll-factor of 5*/
                // compute overhead to unroll factor
                const int M = N % 5;

                // if N is a multiple of 5
                if ( M == 0 )
                {
#    pragma clang loop vectorize(enable)
                    for ( int i = 0; i < N; i += 5 )
                    {
                        this->global_indices_[i] = ilower_ + i;
                        this->global_indices_[i + 1] = ilower_ + i + 1;
                        this->global_indices_[i + 2] = ilower_ + i + 2;
                        this->global_indices_[i + 3] = ilower_ + i + 3;
                        this->global_indices_[i + 4] = ilower_ + i + 4;
                    }
                }
                else
                {
                    // result for overhead to unroll factor
#    pragma clang loop vectorize(enable)
                    for ( int i = 0; i < M; ++i )
                    {
                        this->global_indices_[i] = ilower_ + i;
                    }

                    // result for rest of vectors if length is greater than the unroll factor
                    if ( N > 5 )
                    {
#    pragma clang loop vectorize(enable)
                        for ( int i = M; i < N; i += 5 )
                        {
                            this->global_indices_[i] = ilower_ + i;
                            this->global_indices_[i + 1] = ilower_ + i + 1;
                            this->global_indices_[i + 2] = ilower_ + i + 2;
                            this->global_indices_[i + 3] = ilower_ + i + 3;
                            this->global_indices_[i + 4] = ilower_ + i + 4;
                        }
                    }
                }

                int nvals = local_size;

                HYPRE_IJVectorSetValues ( x_, nvals, this->global_indices_, vec2ptr ( val ) );

                // Finalize initialization of vector
                HYPRE_IJVectorAssemble ( x_ );
                HYPRE_IJVectorGetObject ( x_, ( void ** ) &parcsr_x_ );

                // prepare for communication
                this->nb_sends_ = 0;
                this->nb_recvs_ = 0;
                for ( int id = 0; id < this->nb_procs_; ++id )
                {
                    if ( this->cp_->border_offsets ( id + 1 ) -
                         this->cp_->border_offsets ( id ) > 0 )
                    {
                        this->nb_sends_++;
                    }
                    if ( this->cp_->ghost_offsets ( id + 1 ) -
                         this->cp_->ghost_offsets ( id ) > 0 )
                    {
                        this->nb_recvs_++;
                    }
                }
                this->mpi_req_.resize ( this->nb_sends_ + this->nb_recvs_ );
                this->mpi_stat_.resize ( this->nb_sends_ + this->nb_recvs_ );
                this->border_val_.resize ( this->cp_->size_border_indices ( ) );
                this->ghost_val_.resize ( this->cp_->size_ghost ( ) );

                // set border indices
                border_indices_.resize ( this->cp_->size_border_indices ( ) );
#    pragma clang loop vectorize(enable)
                for ( size_t i = 0; i < border_indices_.size ( ); ++i )
                {
                    border_indices_[i] = this->cp_->border_indices ( )[i] + ilower_;
                }

                this->ghost_ = vec.ghost ( ).CloneWithoutContent ( );
                this->ghost_->CopyStructureFrom ( vec.ghost ( ) );

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

                initialized_ = true;

                val.clear ( );
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreVector<DataType>::print_statistics ( ) const
        {
#ifdef WITH_HYPRE
            int my_rank, num_procs;
            MPI_Comm_rank ( comm_, &my_rank );
            MPI_Comm_size ( comm_, &num_procs );

            int jlower, jupper;
            HYPRE_IJVectorGetLocalRange ( x_, &jlower, &jupper );

            // print statistics
            for ( int i = 0; i < num_procs; ++i )
            {
                MPI_Barrier ( comm_ );
                if ( i == my_rank )
                {
                    std::cout << "HypreVector on process " << my_rank << ":" << std::endl;
                    // print size information
                    std::cout << "\t jlower: " << jlower << std::endl;
                    std::cout << "\t jupper: " << jupper << std::endl;
                }
                MPI_Barrier ( comm_ );
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        int HypreVector<DataType>::size_local ( ) const
        {
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );
            return this->cp_->nb_dofs ( my_rank_ );
        }

        template<class DataType>
        int HypreVector<DataType>::size_global ( ) const
        {
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );
            return this->cp_->nb_total_dofs ( );
        }

        template<class DataType>
        void HypreVector<DataType>::Zeros ( )
        {
#ifdef WITH_HYPRE
            HYPRE_ParVectorSetConstantValues ( parcsr_x_, static_cast < DataType > ( 0 ) );
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
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        DataType HypreVector<DataType>::GetValue ( const int index ) const
        {
            DataType val;
            this->GetValues ( &index, 1, &val );
            return val;
        }

        template<class DataType>
        void HypreVector<DataType>::GetValues ( const int* indices, const int size_indices, DataType* values ) const
        {
#ifdef WITH_HYPRE
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );

            if ( size_indices > 0 )
            {
                assert ( indices != 0 );
                assert ( values != 0 );

                // extract interior and ghost indices (and if possible pp indices)
                std::vector<int> indices_interior ( size_indices );
                int ind_interior = 0;

                std::vector<int> indices_ghost ( size_indices );
                int ind_ghost = 0;

                std::vector<int> indices_pp ( size_indices );
                int ind_pp = 0;

                const size_t local_length = this->iupper_ - this->ilower_;

                // transform indices according to local numbering of interior and ghost
                for ( int i = 0; i < size_indices; i++ )
                {
                    // interior
                    if ( ( unsigned ) ( indices[i] - this->ilower_ ) <= local_length )
                    {
                        indices_interior[ind_interior] = indices[i];
                        ++ind_interior;
                    }

                        // ghost
                    else if ( this->cp_->global2offdiag ( ).find ( indices[i] ) !=
                              this->cp_->global2offdiag ( ).end ( ) )
                    {
                        indices_ghost[ind_ghost] = this->cp_->Global2Offdiag ( indices[i] );
                        ++ind_ghost;
                    }

                        // pp values
                    else if ( ( this->pp_data_ != NULL ) && ( this->pp_data_->sorted_dof_ids.find ( indices[i], &( indices_pp[ind_pp] ) ) ) )
                    {
                        ++ind_pp;
                    }
                }

                // extract values
                std::vector<DataType> values_interior ( ind_interior );
                if ( ind_interior > 0 )
                {
                    HYPRE_IJVectorGetValues ( x_, ind_interior, vec2ptr ( indices_interior ), vec2ptr ( values_interior ) );
                }

                std::vector<DataType> values_ghost ( ind_ghost );

                if ( ind_ghost > 0 )
                {
                    this->ghost_->GetValues ( &indices_ghost.front ( ), ind_ghost, &values_ghost.front ( ) );
                }

                std::vector<DataType> values_pp ( ind_pp );
                if ( this->pp_data_ != NULL )
                {
#    pragma clang loop vectorize(enable)
                    for ( int i = 0; i < ind_pp; ++i )
                    {
                        values_pp[i] = this->pp_data_->values[indices_pp[i]];
                    }
                }
                // put values together
                for ( int i = 0; i < size_indices; i++ )
                {
                    int cont;
                    const size_t current_index = size_indices - i - 1;
                    const size_t index_current = indices[current_index];
                    // interior
                    if ( ( unsigned ) ( index_current - this->ilower_ ) <= local_length )
                    {
                        values[current_index] = values_interior[--ind_interior];
                    }

                        // ghost
                    else if ( this->cp_->global2offdiag ( ).find ( index_current ) !=
                              this->cp_->global2offdiag ( ).end ( ) )
                    {
                        values[current_index] = values_ghost[--ind_ghost];
                    }

                        // pp
                    else if ( ( this->pp_data_ != NULL ) &&
                              ( this->pp_data_->sorted_dof_ids.find ( index_current, &cont ) ) )
                    {
                        values[current_index] = values_pp[--ind_pp];
                    }
                        // zeros
                    else
                    {
                        values[current_index] = static_cast < DataType > ( 0. );
                    }
                }

                assert ( ind_interior == 0 );
                assert ( ind_ghost == 0 );
                assert ( ind_pp == 0 );
                indices_interior.clear ( );
                indices_ghost.clear ( );
                indices_pp.clear ( );
                values_interior.clear ( );
                values_ghost.clear ( );
                values_pp.clear ( );
            }
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        DataType HypreVector<DataType>::Norm2 ( ) const
        {
            return std::sqrt ( this->Dot ( *this ) );
        }

        template<class DataType>
        DataType HypreVector<DataType>::Norm1 ( ) const
        {
            LOG_ERROR ( "Called HypreVector::Norm1 not yet implemented!!!" );
            exit ( -1 );
            return -1.;
        }

        template<class DataType>
        DataType HypreVector<DataType>::NormMax ( ) const
        {
            LOG_ERROR ( "Called HypreVector::NormMax not yet implemented!!!" );
            exit ( -1 );
            return -1.;
        }

        template<class DataType>
        DataType HypreVector<DataType>::Dot ( const Vector<DataType>& vec ) const
        {
            const HypreVector<DataType> *hv = dynamic_cast < const HypreVector<DataType>* > ( &vec );

            if ( hv != 0 )
            {
                return this->Dot ( *hv );
            }
            else
            {
                LOG_ERROR ( "Called HypreVector::Dot with incompatible vector type." );
                exit ( -1 );
                return -1.;
            }
            return -1.;
        }

        template<class DataType>
        DataType HypreVector<DataType>::Dot ( const HypreVector<DataType>& vec ) const
        {
#ifdef WITH_HYPRE
            DataType result = 0.;
            HYPRE_ParVectorInnerProd ( parcsr_x_, *( vec.GetParVector ( ) ), &result );
            return result;
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
            return -1.;
#endif
        }

        template<class DataType>
        void HypreVector<DataType>::Add ( int index, DataType scalar )
        {
            this->Add ( &index, 1, &scalar );
        }

        template<class DataType>
        void HypreVector<DataType>::Add ( const int* indices, int length, const DataType* values )
        {
#ifdef WITH_HYPRE
            HYPRE_IJVectorAddToValues ( x_, length, indices, values );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreVector<DataType>::SetValue ( const int index, const DataType value )
        {
            this->SetValues ( &index, 1, &value );
        }

        template<class DataType>
        void HypreVector<DataType>::SetValues ( const int* indices, const int size_indices, const DataType* values )
        {
#ifdef WITH_HYPRE
            HYPRE_IJVectorSetValues ( x_, size_indices, indices, values );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreVector<DataType>::Axpy ( const Vector<DataType>& vecx, const DataType alpha )
        {
            const HypreVector<DataType> *hv = dynamic_cast < const HypreVector<DataType>* > ( &vecx );

            if ( hv != 0 )
            {
                this->Axpy ( *hv, alpha );
            }
            else
            {
                LOG_ERROR ( "Called HypreVector::Axpy with incompatible vector type." );
                exit ( -1 );
            }
        }

        template<class DataType>
        void HypreVector<DataType>::Axpy ( const HypreVector<DataType>& vecx, const DataType alpha )
        {
#ifdef WITH_HYPRE
            HYPRE_ParVectorAxpy ( alpha, *( vecx.GetParVector ( ) ), parcsr_x_ );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreVector<DataType>::ScaleAdd ( const Vector<DataType>& vecx, const DataType alpha )
        {
            this->Scale ( alpha );
            this->Axpy ( vecx, static_cast < DataType > ( 1. ) );
        }

        template<class DataType>
        void HypreVector<DataType>::Scale ( const DataType alpha )
        {
#ifdef WITH_HYPRE
            HYPRE_ParVectorScale ( alpha, parcsr_x_ );
#else
            LOG_ERROR ( "HiFlow was not compiled with HYPRE support!" );
            exit ( -1 );
#endif
        }

        template<class DataType>
        void HypreVector<DataType>::SendBorder ( )
        {
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );
            assert ( this->comm_ != MPI_COMM_NULL );

            this->GetValues ( vec2ptr ( border_indices_ ), border_indices_.size ( ), vec2ptr ( border_val_ ) );

            int tag = 1;
            int ctr = 0;

            for ( int id = 0; id < this->nb_procs_; id++ )
            {
                if ( this->cp_->border_offsets ( id + 1 ) -
                     this->cp_->border_offsets ( id ) > 0 )
                {
                    int info = MPI_Isend ( &( this->border_val_[0] ) + this->cp_->border_offsets ( id ),
                                           this->cp_->border_offsets ( id + 1 ) -
                                           this->cp_->border_offsets ( id ),
                                           mpi_data_type<DataType>::get_type ( ), id, tag, this->comm_,
                                           &( this->mpi_req_[this->nb_recvs_ + ctr] ) );
                    assert ( info == MPI_SUCCESS );
                    ctr++;
                }
            }
            assert ( ctr == this->nb_sends_ );
        }

        template<class DataType>
        void HypreVector<DataType>::ReceiveGhost ( )
        {
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );
            assert ( this->my_rank_ >= 0 );
            assert ( this->my_rank_ < this->nb_procs_ );
            assert ( this->comm_ != MPI_COMM_NULL );

            int tag = 1;
            int ctr = 0;

            for ( int id = 0; id < this->nb_procs_; ++id )
            {
                if ( this->cp_->ghost_offsets ( id + 1 ) -
                     this->cp_->ghost_offsets ( id ) > 0 )
                {
                    int info = MPI_Irecv ( &( this->ghost_val_[0] ) + this->cp_->ghost_offsets ( id ),
                                           this->cp_->ghost_offsets ( id + 1 ) -
                                           this->cp_->ghost_offsets ( id ),
                                           mpi_data_type<DataType>::get_type ( ), id, tag, this->comm_,
                                           &( this->mpi_req_[ctr] ) );
                    assert ( info == MPI_SUCCESS );
                    ctr++;
                }
            }
            assert ( ctr == this->nb_recvs_ );
        }

        template<class DataType>
        void HypreVector<DataType>::WaitForSend ( )
        {
            int info = MPI_Waitall ( this->nb_sends_,
                                     &( this->mpi_req_[this->nb_recvs_] ),
                                     &( this->mpi_stat_[this->nb_recvs_] ) );
            assert ( info == MPI_SUCCESS );
        }

        template<class DataType>
        void HypreVector<DataType>::WaitForRecv ( )
        {
            int info = MPI_Waitall ( this->nb_recvs_,
                                     &( this->mpi_req_[0] ),
                                     &( this->mpi_stat_[0] ) );
            assert ( info == MPI_SUCCESS );

            this->SetGhostValues ( &( this->ghost_val_[0] ) );
        }

        template<class DataType>
        void HypreVector<DataType>::Update ( )
        {

            this->UpdateGhost ( );

            // update post processing if possible.
            if ( this->pp_data_ != NULL )
            {
                this->UpdatePpValues ( );
            }
            else if ( !checked_for_dof_partition_ )
            {
                const Couplings<DataType>* couplings = dynamic_cast < const Couplings<DataType>* > ( cp_ );
                if ( couplings != NULL )
                {
                    this->InitializePostProcessing ( couplings );
                    this->UpdatePpValues ( );
                }
                checked_for_dof_partition_ = true;
            }
        }

        template<class DataType>
        void HypreVector<DataType>::UpdateGhost ( )
        {
            assert ( this->ghost_ != NULL );

            this->ReceiveGhost ( );
            this->SendBorder ( );

            this->WaitForRecv ( );
            this->WaitForSend ( );
        }

        template<class DataType>
        void HypreVector<DataType>::InitializePostProcessing ( const Couplings<DataType>* couplings )
        {
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );
            assert ( this->comm_ != MPI_COMM_NULL );

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

        template<class DataType>
        void HypreVector<DataType>::UpdatePpValues ( )
        {
            assert ( pp_data_ != NULL );
            assert ( this->comm_ != MPI_COMM_NULL );

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
                    info = MPI_Isend ( &( send_pp_values[0] ) + this->pp_data_->send_offsets[id],
                                       this->pp_data_->nb_sends[id], mpi_data_type<DataType>::get_type ( ),
                                       id, tag, this->comm_, &( mpi_req_recv[id] ) );
                    assert ( info == MPI_SUCCESS );
                }
                if ( this->pp_data_->nb_recvs[id] > 0 )
                {
                    info = MPI_Irecv ( &( this->pp_data_->values[0] ) + this->pp_data_->dof_ids_offsets[id],
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

            mpi_req_send.clear ( );
            mpi_req_recv.clear ( );
            send_pp_values.clear ( );
        }

        template<class DataType>
        void HypreVector<DataType>::SetGhostValues ( const DataType* values )
        {
            assert ( this->ghost_ != NULL );

            this->ghost_->SetBlockValues ( 0, this->size_local_ghost ( ), values );
        }

        template<class DataType>
        void HypreVector<DataType>::GetAllDofsAndValues ( std::vector<int>& id, std::vector<DataType>& val ) const
        {
            assert ( this->cp_ != NULL );
            assert ( this->cp_->initialized ( ) );
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
#pragma clang loop vectorize(enable)
            for ( int i = 0; i < this->size_local ( ); ++i )
            {
                id[i] = ilower_ + i;
            }

            // Combine interior, ghost and pp_data values
            this->GetValues ( vec2ptr ( id ), this->size_local ( ), values );
            this->ghost_->GetBlockValues ( 0, this->size_local_ghost ( ), ghost_values );

            val.resize ( total_size );
            // First: values from the interior
#pragma clang loop vectorize(enable)
            for ( int i = 0; i < this->size_local ( ); ++i )
            {
                val[i] = values[i];
            }
            delete[] values;
            // Second: the DoFs and values from ghost
            int tmp_offset = this->size_local ( );
#pragma clang loop vectorize(enable)
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
#pragma clang loop vectorize(enable)
                for ( int i = 0, e_i = this->pp_data_->dof_ids_offsets[this->nb_procs_]; i < e_i; ++i )
                {
                    id[i + tmp_offset] = this->pp_data_->dof_ids[i];
                    val[i + tmp_offset] = this->pp_data_->values[i];
                }
            }
        }

#ifdef WITH_HYPRE

        template<class DataType>
        HYPRE_ParVector* HypreVector<DataType>::GetParVector ( )
        {
            return &parcsr_x_;
        }

        template<class DataType>
        const HYPRE_ParVector* HypreVector<DataType>::GetParVector ( ) const
        {
            return &parcsr_x_;
        }
#endif

        template<class DataType>
        void HypreVector<DataType>::WriteHDF5 ( const std::string& filename,
                                                const std::string& groupname,
                                                const std::string& datasetname )
        {
            assert ( this->comm_ != MPI_COMM_NULL );
#ifdef WITH_HDF5
            // Define Data in memory
            const size_t local_size = this->size_local ( );
            std::vector<DataType> data ( local_size, 0. );

            this->GetValues ( this->global_indices_, local_size, vec2ptr ( data ) );

            H5FilePtr file_ptr ( new H5File ( filename, "w", this->comm_ ) );
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname, "w" ) );
            H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, this->size_global ( ),
                                                       datasetname, "w", vec2ptr ( data ) ) );
            dataset_ptr->write ( this->size_local ( ), this->ilower_, vec2ptr ( data ) );

            data.clear ( );
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif

        }

        template<class DataType>
        void HypreVector<DataType>::ReadHDF5 ( const std::string& filename,
                                               const std::string& groupname,
                                               const std::string& datasetname )
        {
            assert ( this->comm_ != MPI_COMM_NULL );
#ifdef WITH_HDF5
            DataType* buffer;
            buffer = new DataType[this->size_local ( )];

            H5FilePtr file_ptr ( new H5File ( filename, "r", this->comm_ ) );
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname, "r" ) );
            H5DatasetPtr dataset_ptr ( new H5Dataset ( group_ptr, this->size_global ( ),
                                                       datasetname, "r", buffer ) );
            dataset_ptr->read ( this->size_local ( ), this->ilower_, buffer );

            const size_t local_size = this->size_local ( );

            this->SetValues ( this->global_indices_, local_size, buffer );

            // Update
            this->UpdateGhost ( );

            delete buffer;
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        // template instantiation
        template class HypreVector<double>;

    } // namespace la
} // namespace hiflow
