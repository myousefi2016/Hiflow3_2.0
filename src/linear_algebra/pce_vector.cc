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

#include "pce_vector.h"

namespace hiflow
{
    namespace la
    {

        // constructor

        template<class DataType>
        PCEVector<DataType>::PCEVector ( )
        {
            this->nmode_ = -1;
            this->totlevel_ = -1;
            this->initialized_ = false;
        }

        // deconstructor

        template<class DataType>
        PCEVector<DataType>::~PCEVector ( )
        {

            this->nmode_ = -1;
            this->totlevel_ = -1;
            this->initialized_ = false;
            this->Clear ( );
        }

        // initialize

        template<class DataType>
        void PCEVector<DataType>::Init ( const PCTensor& pctensor, const HypreVector<DataType>& mean_vector )
        {

            // assign pctensor_
            this->pctensor_ = pctensor;

            // calculate nb of modes
            this->nmode_ = pctensor.Size ( );

            // calculate total level
            this->totlevel_ = pctensor_.GetLevel ( );

            // allocate nb of modes
            this->mode_vector_.resize ( this->nmode_ );

            // initialize with mean_vector
            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->mode_vector_[mode].CloneFromWithoutContent ( mean_vector );
            }

            this->initialized_ = true;
        }

        template<class DataType>
        void PCEVector<DataType>::Init ( const PCTensor& pctensor, const MPI_Comm& comm, const LaCouplings& cp )
        {

            // assign pctensor_
            this->pctensor_ = pctensor;

            // calculate nb of modes
            this->nmode_ = pctensor_.Size ( );

            // calculate total level
            this->totlevel_ = pctensor.GetLevel ( );

            // allocate nb of modes
            this->mode_vector_.resize ( this->nmode_ );

            // tmp vec
            HypreVector<DataType> tmp;
            tmp.Init ( comm, cp );

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->mode_vector_[mode].CloneFromWithoutContent ( tmp );
            }

            this->initialized_ = true;
        }

        // access the member of mode_vector_

        template<class DataType>
        HypreVector<DataType>& PCEVector<DataType>::Mode ( const int& mode )
        {
            assert ( this->initialized_ == true );
            return this->mode_vector_[mode];
        }

        template<class DataType>
        const HypreVector<DataType>& PCEVector<DataType>::GetMode ( const int& mode ) const
        {
            assert ( this->initialized_ == true );
            return this->mode_vector_[mode];
        }

        // return number of modes

        template<class DataType>
        int PCEVector<DataType>::nb_mode ( ) const
        {
            assert ( this->initialized_ == true );
            return this->nmode_;
        }

        // return total level

        template<class DataType>
        int PCEVector<DataType>::total_level ( ) const
        {
            assert ( this->initialized_ == true );
            return this->totlevel_;
        }

        // return pctensor

        template<class DataType>
        PCTensor PCEVector<DataType>::GetPCTensor ( ) const
        {
            assert ( this->initialized_ == true );
            return this->pctensor_;
        }

        // Update

        template<class DataType>
        void PCEVector<DataType>::Update ( )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->mode_vector_[mode].Update ( );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::Update ( const int& mode )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( mode >= 0 );

            this->mode_vector_[mode].Update ( );
        }

        // Clear

        template<class DataType>
        void PCEVector<DataType>::Clear ( )
        {
            if ( this->initialized_ == true )
            {
                for ( int mode = 0; mode != this->nmode_; ++mode )
                {
                    this->mode_vector_[mode].Clear ( );
                }
            }
        }

        // Zeros

        template<class DataType>
        void PCEVector<DataType>::Zeros ( )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->mode_vector_[mode].Zeros ( );
            }
        }

        // Axpy

        template<class DataType>
        void PCEVector<DataType>::Axpy ( const PCEVector<DataType>& vec, const DataType& alpha )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                Axpy ( mode, vec.GetMode ( mode ), alpha );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::Axpy ( const PCEVector<DataType>& vec, const DataType& alpha, const int& l )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( l >= 0 );
            assert ( this->totlevel_ >= 0 );

            for ( int mode = 0; mode != this->pctensor_.Size ( l ); ++mode )
            {
                Axpy ( mode, vec.GetMode ( mode ), alpha );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::Axpy ( const int& mode, const PCEVector<DataType>& vec, const DataType& alpha )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( this->totlevel_ >= 0 );
            assert ( mode >= 0 );
            this->mode_vector_[mode].Axpy ( vec.GetMode ( mode ), alpha );
        }

        template<class DataType>
        void PCEVector<DataType>::Axpy ( const int& mode, const HypreVector<DataType>& vec, const DataType& alpha )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( this->totlevel_ >= 0 );
            assert ( mode >= 0 );
            this->mode_vector_[mode].Axpy ( vec, alpha );
        }

        // Dot

        template<class DataType>
        DataType PCEVector<DataType>::Dot ( const PCEVector<DataType>& vec ) const
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );

            DataType dot = 0.0;
            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                dot += this->Dot ( mode, vec.GetMode ( mode ) );
            }

            return dot;
        }

        template<class DataType>
        DataType PCEVector<DataType>::Dot ( const int& mode, const HypreVector<DataType>& vec ) const
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( mode >= 0 );

            return this->mode_vector_[mode].Dot ( vec );
        }

        // ScaleAdd

        template<class DataType>
        void PCEVector<DataType>::ScaleAdd ( const PCEVector<DataType>& vec, const DataType& alpha )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->ScaleAdd ( mode, vec.GetMode ( mode ), alpha );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::ScaleAdd ( const int& mode, const HypreVector<DataType>& vec, const DataType& alpha )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( mode >= 0 );

            this->mode_vector_[mode].ScaleAdd ( vec, alpha );
        }

        // Scale

        template<class DataType>
        void PCEVector<DataType>::Scale ( const DataType& alpha )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->mode_vector_[mode].Scale ( alpha );
            }
        }

        // size

        template<class DataType>
        int PCEVector<DataType>::size_local ( ) const
        {
            assert ( this->initialized_ == true );
            return (this->nmode_ * this->size_local ( 0 ) );
        }

        template<class DataType>
        int PCEVector<DataType>::size_local ( const int& mode ) const
        {
            assert ( this->initialized_ == true );
            return this->mode_vector_[mode].size_local ( );
        }

        template<class DataType>
        int PCEVector<DataType>::size_global ( ) const
        {
            assert ( this->initialized_ == true );
            return (this->nmode_ * this->size_global ( 0 ) );
        }

        template<class DataType>
        int PCEVector<DataType>::size_global ( const int& mode ) const
        {
            assert ( this->initialized_ == true );
            return this->mode_vector_[mode].size_global ( );
        }

        template<class DataType>
        int PCEVector<DataType>::size_local_ghost ( ) const
        {
            assert ( this->initialized_ == true );
            return (this->nmode_ * this->size_local_ghost ( 0 ) );
        }

        template<class DataType>
        int PCEVector<DataType>::size_local_ghost ( const int& mode ) const
        {
            assert ( this->initialized_ == true );
            return this->mode_vector_[mode].size_local_ghost ( );
        }

        // Norm2

        template<class DataType>
        DataType PCEVector<DataType>::Norm2 ( ) const
        {
            assert ( this->initialized_ == true );
            return sqrt ( this->Dot ( *this ) );
        }

        // CloneFrom

        template<class DataType>
        void PCEVector<DataType>::CloneFrom ( const PCEVector<DataType>& vec )
        {
            if ( this->initialized_ == false )
            {
                this->Init ( vec.GetPCTensor ( ), vec.GetMode ( 0 ) );
            }

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->CloneFrom ( mode, vec.GetMode ( mode ) );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::CloneFrom ( const PCEVector<DataType>& vec, const int& l )
        {
            if ( this->initialized_ == false )
            {
                PCTensor pctensor = vec.GetPCTensor ( );
                pctensor.SetLevel ( l );

                // restrict to own size is smaller than input
                assert ( pctensor.Size ( ) < vec.nb_mode ( ) );
                this->Init ( pctensor, vec.GetMode ( 0 ) );
            }

            for ( int mode = 0; mode != this->pctensor_.Size ( l ); ++mode )
            {
                this->CloneFrom ( mode, vec.GetMode ( mode ) );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::CloneFrom ( const int& mode, const HypreVector<DataType>& vec )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( mode >= 0 );
            this->mode_vector_[mode].CloneFrom ( vec );
        }

        // CopyFrom

        template<class DataType>
        void PCEVector<DataType>::CopyFrom ( const PCEVector<DataType>& vec )
        {
            if ( this->initialized_ == false )
            {
                this->Init ( vec.GetPCTensor ( ), vec.GetMode ( 0 ) );
            }

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->CopyFrom ( mode, vec.GetMode ( mode ) );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::CopyFrom ( const PCEVector<DataType>& vec, const int& l )
        {
            if ( this->initialized_ == false )
            {
                PCTensor pctensor = vec.GetPCTensor ( );
                pctensor.SetLevel ( l );

                // restrict to own size is smaller than input
                assert ( pctensor.Size ( ) < vec.nb_mode ( ) );
                this->Init ( vec.GetPCTensor ( ), vec.GetMode ( 0 ) );
            }

            for ( int mode = 0; mode != this->pctensor_.Size ( l ); ++mode )
            {
                this->CopyFrom ( mode, vec.GetMode ( mode ) );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::CopyFrom ( const int& mode, const HypreVector<DataType>& vec )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( mode >= 0 );
            this->mode_vector_[mode].CopyFrom ( vec );

        }

        // CopyFromWithoutGhost

        template<class DataType>
        void PCEVector<DataType>::CopyFromWithoutGhost ( const PCEVector<DataType>& vec )
        {
            if ( this->initialized_ == false )
            {
                this->Init ( vec.GetPCTensor ( ), vec.GetMode ( 0 ) );
            }

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->CopyFromWithoutGhost ( mode, vec.GetMode ( mode ) );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::CopyFromWithoutGhost ( const PCEVector<DataType>& vec, const int& l )
        {
            if ( this->initialized_ == false )
            {
                PCTensor pctensor = vec.GetPCTensor ( );
                pctensor.SetLevel ( l );

                // restrict to own size is smaller than input
                assert ( pctensor.Size ( ) < vec.nb_mode ( ) );
                this->Init ( vec.GetPCTensor ( ), vec.GetMode ( 0 ) );
            }

            for ( int mode = 0; mode != this->pctensor_.Size ( l ); ++mode )
            {
                this->CopyFromWithoutGhost ( mode, vec.GetMode ( mode ) );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::CopyFromWithoutGhost ( const int& mode, const HypreVector<DataType>& vec )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( mode >= 0 );
            this->mode_vector_[mode].CopyFromWithoutGhost ( vec );
        }

        // CloneFromWithoutContent

        template<class DataType>
        void PCEVector<DataType>::CloneFromWithoutContent ( const PCEVector<DataType>& vec )
        {
            if ( this->initialized_ == false )
            {
                this->Init ( vec.GetPCTensor ( ), vec.GetMode ( 0 ) );
            }

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                this->CloneFromWithoutContent ( mode, vec.GetMode ( mode ) );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::CloneFromWithoutContent ( const PCEVector<DataType>& vec, const int& l )
        {
            if ( this->initialized_ == false )
            {
                PCTensor pctensor = vec.GetPCTensor ( );
                pctensor.SetLevel ( l );

                // restrict to own size is smaller than input
                assert ( pctensor.Size ( ) < vec.nb_mode ( ) );

                this->Init ( vec.GetPCTensor ( ), vec.GetMode ( 0 ) );
            }

            for ( int mode = 0; mode != this->pctensor_.Size ( l ); ++mode )
            {
                this->CloneFromWithoutContent ( mode, vec.GetMode ( mode ) );
            }
        }

        template<class DataType>
        void PCEVector<DataType>::CloneFromWithoutContent ( const int& mode, const HypreVector<DataType>& vec )
        {
            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );
            assert ( mode >= 0 );

            this->mode_vector_[mode].CloneFromWithoutContent ( vec );
        }

        // Write and Read vector content to/from HDF5 file

        template<class DataType>
        void PCEVector<DataType>::WriteHDF5 ( const std::string& filename,
                                              const std::string& groupname,
                                              const std::string& datasetname )
        {

            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );

#ifdef WITH_HDF5

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                std::ostringstream groupname_mode;
                groupname_mode << groupname << "_mode" << mode;

                std::ostringstream datasetname_mode;
                datasetname_mode << datasetname << ".Mode(" << mode << ")";

                this->mode_vector_[mode].WriteHDF5 ( filename, groupname_mode.str ( ), datasetname_mode.str ( ) );
            }

#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        template<class DataType>
        void PCEVector<DataType>::ReadHDF5 ( const std::string& filename,
                                             const std::string& groupname,
                                             const std::string& datasetname )
        {

            assert ( this->initialized_ == true );
            assert ( this->nmode_ > 0 );

#ifdef WITH_HDF5

            for ( int mode = 0; mode != this->nmode_; ++mode )
            {
                std::ostringstream groupname_mode;
                groupname_mode << groupname << "_mode" << mode;

                std::ostringstream datasetname_mode;
                datasetname_mode << datasetname << ".Mode(" << mode << ")";

                this->mode_vector_[mode].ReadHDF5 ( filename, groupname_mode.str ( ), datasetname_mode.str ( ) );
            }

#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        template class PCEVector<double>;

    }
}
