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

/// \author Michael Schick

#include "polynomial_chaos/pc_galerkin_vector.h"
#include "linear_algebra/lmp/lvector.h"
#include <cmath>

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class DataType>
        PCGalerkinVector<DataType>::PCGalerkinVector ( )
        {
            nmodes_ = -1;
            noff_modes_ = -1;
            nmodes_total_ = -1;
            comm_ = MPI_COMM_NULL;
        }

        template<class DataType>
        PCGalerkinVector<DataType>::PCGalerkinVector ( PCGalerkinVector<DataType> const& vec )
        {
            this->Clear ( );

            nmodes_ = vec.NModes ( );
            noff_modes_ = vec.NOffModes ( );
            nmodes_total_ = nmodes_ + noff_modes_;

            modes_.clear ( );

            modes_.resize ( nmodes_total_ );

            for ( int n = 0; n < nmodes_total_; ++n )
            {
                modes_[n] = new la::CoupledVector<DataType>;
                modes_[n]->CloneFrom ( *( vec.Mode ( n ) ) );
            }

            SetMPIComm ( vec.comm ( ) );
            Initialize ( );
        }

        template<class DataType>
        PCGalerkinVector<DataType>::~PCGalerkinVector<DataType> ( )
        {
            this->Clear ( );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::ModeAxpy ( int mode, la::CoupledVector<DataType> const& vec, const DataType alpha )
        {
            assert ( mode < nmodes_total_ && nmodes_total_ > -1 );
            modes_[mode]->Axpy ( vec, alpha );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::Axpy ( PCGalerkinVector<DataType> const& vec, DataType const alpha )
        {
            int min = nmodes_;
            if ( vec.NModes ( ) < min ) min = vec.NModes ( );

            for ( int n = 0; n < min; ++n )
            {
                modes_[n]->Axpy ( *( vec.Mode ( n ) ), alpha );
            }
        }

        template<class DataType>
        DataType PCGalerkinVector<DataType>::Dot ( const PCGalerkinVector<DataType>& vec ) const
        {
            int min = nmodes_;
            if ( vec.NModes ( ) < min ) min = vec.NModes ( );

            DataType dot = 0.0;
            for ( int n = 0; n < min; ++n )
            {
                dot += modes_[n]->Dot ( *( vec.Mode ( n ) ) );
            }

            return dot;
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::Zeros ( )
        {
            for ( int i = 0; i < nmodes_; ++i )
            {
                modes_[i]->Zeros ( );
            }
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::ScaleAdd ( const PCGalerkinVector<DataType>& vec, const DataType alpha )
        {
            int min = nmodes_;
            if ( vec.NModes ( ) < min ) min = vec.NModes ( );

            for ( int n = 0; n < min; ++n )
            {
                modes_[n]->ScaleAdd ( *( vec.Mode ( n ) ), alpha );
            }
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::Scale ( const DataType alpha )
        {
            for ( int n = 0; n < nmodes_; ++n )
            {
                modes_[n]->Scale ( alpha );
            }
        }

        template<class DataType>
        DataType PCGalerkinVector<DataType>::Norm2 ( ) const
        {
            return sqrt ( this->Dot ( *this ) );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::SynchronizeModes ( GalerkinTensor* pctensor )
        {
            int my_rank = pctensor->MyRank ( );
            int numproc = pctensor->NumProc ( );
            int nsend = 0;
            for ( int rank = 0; rank < numproc; ++rank )
            {
                if ( rank != my_rank ) nsend += pctensor->SendHostileModes ( rank ).size ( );
            }

            std::vector<MPI_Request> request ( nsend );
            std::vector<MPI_Status> status ( nsend );
            int cntr = 0;

            // Who wants something from me -> MPI Send
            for ( int rank = 0; rank < numproc; ++rank )
            {
                if ( rank != my_rank )
                {
                    for ( int n = 0; n < pctensor->SendHostileModes ( rank ).size ( ); ++n )
                    {
                        //Extract Data
                        int glo_idx = pctensor->SendHostileModes ( rank )[n];
                        int loc_idx = pctensor->G2L ( glo_idx );
                        int tag = 4;
                        int info = MPI_Isend ( this->Mode ( loc_idx )->interior ( ).GetBuffer ( ), this->DSizeLocal ( ), MPI_DOUBLE,
                                               rank, tag, pctensor->Comm ( ), &request[cntr] );
                        assert ( info == MPI_SUCCESS );
                        cntr++;
                    }

                    for ( int n = 0; n < pctensor->RecvHostileModes ( rank ).size ( ); ++n )
                    {
                        int glo_idx = pctensor->RecvHostileModes ( rank )[n];
                        int loc_idx = pctensor->G2L ( glo_idx );
                        int tag = 4;

                        MPI_Status status;
                        MPI_Recv ( this->Mode ( loc_idx )->interior ( ).GetBuffer ( ), this->DSizeLocal ( ), MPI_DOUBLE, rank, tag, pctensor->Comm ( ), &status );
                    }
                }
            }
            int info_send = MPI_Waitall ( nsend, &request.front ( ), &status.front ( ) );
            assert ( info_send == MPI_SUCCESS );

            MPI_Barrier ( pctensor->Comm ( ) );

            SynchronizeModesGhost ( pctensor );
            UpdateCouplings ( );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::SynchronizeModesGhost ( GalerkinTensor* pctensor )
        {
            int my_rank = pctensor->MyRank ( );
            int numproc = pctensor->NumProc ( );
            int nsend = 0;
            for ( int rank = 0; rank < numproc; ++rank )
            {
                if ( rank != my_rank ) nsend += pctensor->SendHostileModes ( rank ).size ( );
            }

            std::vector<MPI_Request> request ( nsend );
            std::vector<MPI_Status> status ( nsend );
            int cntr = 0;

            // Who wants something from me -> MPI Send
            for ( int rank = 0; rank < numproc; ++rank )
            {
                if ( rank != my_rank )
                {
                    for ( int n = 0; n < pctensor->SendHostileModes ( rank ).size ( ); ++n )
                    {
                        //Extract Data
                        int glo_idx = pctensor->SendHostileModes ( rank )[n];
                        int loc_idx = pctensor->G2L ( glo_idx );
                        int tag = 5;
                        int info = MPI_Isend ( this->Mode ( loc_idx )->ghost ( ).GetBuffer ( ), this->Mode ( loc_idx )->size_local_ghost ( ), MPI_DOUBLE,
                                               rank, tag, pctensor->Comm ( ), &request[cntr] );
                        assert ( info == MPI_SUCCESS );
                        cntr++;
                    }

                    for ( int n = 0; n < pctensor->RecvHostileModes ( rank ).size ( ); ++n )
                    {
                        int glo_idx = pctensor->RecvHostileModes ( rank )[n];
                        int loc_idx = pctensor->G2L ( glo_idx );
                        int tag = 5;

                        MPI_Status status;
                        MPI_Recv ( this->Mode ( loc_idx )->ghost ( ).GetBuffer ( ), this->Mode ( loc_idx )->size_local_ghost ( ), MPI_DOUBLE, rank, tag, pctensor->Comm ( ), &status );
                    }
                }
            }
            int info_send = MPI_Waitall ( nsend, &request.front ( ), &status.front ( ) );
            assert ( info_send == MPI_SUCCESS );
            MPI_Barrier ( pctensor->Comm ( ) );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::CreateCoarseVector ( GalerkinTensor* pctensor, PCGalerkinVector<DataType>& coarse_vec ) const
        {
            //1. Allocate data in coarse vec
            int n_coarse_modes = pctensor->SizeLocal ( );
            int n_off_coarse_modes = pctensor->SizeOffModes ( );
            std::vector<la::CoupledVector<DataType>* > coarse_modes ( n_coarse_modes + n_off_coarse_modes );
            for ( int i = 0; i < coarse_modes.size ( ); ++i )
            {
                coarse_modes[i] = new la::CoupledVector<DataType>;
                coarse_modes[i]->CloneFrom ( *this->Mode ( 0 ) ); //Use dummy mode instead -> otherwise bug for low number of rvs?
            }
            coarse_vec.SetModes ( coarse_modes, n_coarse_modes, n_off_coarse_modes );
            coarse_vec.SetMPIComm ( this->comm ( ) );
            for ( int i = 0; i < coarse_modes.size ( ); ++i )
            {
                coarse_modes[i]->Clear ( );
                delete coarse_modes[i];
            }

            //2. Copy data
            int coarse_level = pctensor->GetLevel ( );
            int my_rank = pctensor->MyRank ( );

            for ( int lid_coarse = 0; lid_coarse < n_coarse_modes; ++lid_coarse )
            {
                int gid_coarse = pctensor->L2G ( lid_coarse, coarse_level );
                if ( pctensor->Ownership ( gid_coarse, coarse_level + 1 ) == my_rank )
                {
                    coarse_vec.Mode ( lid_coarse )->CopyFrom ( *this->Mode ( pctensor->G2L ( gid_coarse, coarse_level + 1 ) ) );
                }
            }
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::AssignToRefinedVector ( GalerkinTensor* pctensor, const PCGalerkinVector<DataType>& coarse_vec )
        {
            // 1. Allocate data in coarse vec
            //int n_coarse_modes = coarse_vec.NModes();

            int fine_level = pctensor->GetLevel ( );
            int my_rank = pctensor->MyRank ( );

            for ( int lid_fine = 0; lid_fine < nmodes_; ++lid_fine )
            {
                int gid_fine = pctensor->L2G ( lid_fine, fine_level );
                if ( gid_fine < pctensor->Size ( fine_level - 1 ) )
                {
                    if ( pctensor->Ownership ( gid_fine, fine_level - 1 ) == my_rank )
                    {
                        this->Mode ( lid_fine )->CopyFrom ( *coarse_vec.Mode ( pctensor->G2L ( gid_fine, fine_level - 1 ) ) );
                    }
                }
            }
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::Initialize ( )
        {
            if ( modes_.size ( ) > 0 )
                size_local_deterministic_ = modes_[0]->size_local ( );
            else
                size_local_deterministic_ = 0;
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::Clear ( )
        {
            if ( nmodes_total_>-1 )
            {
                for ( int i = 0; i < nmodes_total_; ++i )
                {
                    modes_[i]->Clear ( );
                    delete modes_[i];
                }
            }
            nmodes_total_ = -1;
        }

        template<class DataType>
        la::CoupledVector<DataType>* PCGalerkinVector<DataType>::Mode ( int n ) const
        {
            assert ( n < nmodes_total_ );
            return modes_[n];
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::CopyMode ( int mode, la::CoupledVector<DataType> const& vec )
        {
            assert ( mode < nmodes_total_ && nmodes_total_ > -1 );
            modes_[mode]->CopyFrom ( vec );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::SetMPIComm ( const MPI_Comm& comm )
        {
            assert ( comm != MPI_COMM_NULL );
            if ( this->comm_ != MPI_COMM_NULL )
                MPI_Comm_free ( &this->comm_ );
            int my_rank;
            int info = MPI_Comm_rank ( comm, &my_rank );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank >= 0 );
            int numproc;
            info = MPI_Comm_size ( comm, &numproc );
            assert ( info == MPI_SUCCESS );
            assert ( numproc > 0 );
            assert ( my_rank < numproc );

            info = MPI_Comm_split ( comm, 0, my_rank, &comm_ );
            assert ( info == MPI_SUCCESS );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::SetModes ( std::vector<la::CoupledVector<DataType>* > const& modes, int nmodes, int noff_modes )
        {
            if ( nmodes_total_ > -1 ) this->Clear ( );

            nmodes_ = nmodes;
            noff_modes_ = noff_modes;

            nmodes_total_ = modes.size ( );
            assert ( nmodes_total_ == nmodes_ + noff_modes_ );

            modes_.resize ( nmodes_total_ );
            for ( int i = 0; i < nmodes_total_; ++i )
            {
                modes_[i] = new la::CoupledVector<DataType>;
                modes_[i]->CloneFrom ( *modes[i] );
            }
            Initialize ( );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::CloneFrom ( const PCGalerkinVector<DataType>& vec )
        {
            this->Clear ( );

            nmodes_ = vec.NModes ( );
            noff_modes_ = vec.NOffModes ( );

            nmodes_total_ = nmodes_ + noff_modes_;

            modes_.clear ( );
            modes_.resize ( nmodes_total_ );
            for ( int i = 0; i < nmodes_total_; ++i )
            {
                modes_[i] = new la::CoupledVector<DataType>;
                modes_[i]->CloneFrom ( *( vec.Mode ( i ) ) );
            }

            SetMPIComm ( vec.comm ( ) );
            Initialize ( );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::CopyFrom ( const PCGalerkinVector<DataType>& vec )
        {
            assert ( nmodes_ == vec.NModes ( ) && noff_modes_ == vec.NOffModes ( ) );
            for ( int i = 0; i < nmodes_total_; ++i )
            {
                modes_[i]->CopyFrom ( *( vec.Mode ( i ) ) );
            }
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::CloneFromWithoutContent ( const PCGalerkinVector<DataType>& vec )
        {
            this->Clear ( );
            nmodes_ = vec.NModes ( );
            noff_modes_ = vec.NOffModes ( );
            nmodes_total_ = nmodes_ + noff_modes_;

            modes_.resize ( nmodes_total_ );
            for ( int i = 0; i < nmodes_total_; ++i )
            {
                modes_[i] = new la::CoupledVector<DataType>;
                modes_[i]->CloneFromWithoutContent ( *( vec.Mode ( i ) ) );
            }

            SetMPIComm ( vec.comm ( ) );
            Initialize ( );
        }

        template class PCGalerkinVector<double>;
    }
}
