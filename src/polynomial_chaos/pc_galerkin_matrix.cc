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

#include "polynomial_chaos/pc_galerkin_matrix.h"
#include "common/timer.h"
#include <cassert>
#ifdef WITH_OPENMP
#    include <omp.h>
#endif

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class DataType>
        PCGalerkinMatrix<DataType>::PCGalerkinMatrix ( la::CoupledVector<DataType> const& dummy )
        {
            numthreads_ = 1;
            modes_allocated_ = false;
#ifdef WITH_OPENMP
            omp_set_num_threads ( numthreads_ );
#endif
            tmp_vec_.resize ( 1 );
            tmp_vec_[0] = new la::CoupledVector<DataType>;
            tmp_vec_[0]->CloneFrom ( dummy );
        }

        template<class DataType>
        PCGalerkinMatrix<DataType>::~PCGalerkinMatrix ( )
        {
            for ( int thread = 0; thread < numthreads_; ++thread )
            {
                delete tmp_vec_[thread];
            }
            tmp_vec_.clear ( );
            if ( modes_allocated_ )
            {
                for ( int i = 0; i < modes_.size ( ); ++i )
                    delete modes_[i];
            }
            modes_.clear ( );
        }

        template<class DataType>
        void PCGalerkinMatrix<DataType>::SetNumThreads ( int numthreads )
        {
            la::CoupledVector<DataType> dummy;
            dummy.CloneFrom ( *tmp_vec_[0] );

            for ( int thread = 0; thread < numthreads_; ++thread )
            {
                delete tmp_vec_[thread];
            }
            tmp_vec_.clear ( );

            numthreads_ = numthreads;
#ifdef WITH_OPENMP
            omp_set_num_threads ( numthreads_ );
#endif

            tmp_vec_.resize ( numthreads_ );
            for ( int thread = 0; thread < numthreads_; ++thread )
            {
                tmp_vec_[thread] = new la::CoupledVector<DataType>;
                tmp_vec_[thread]->CloneFrom ( dummy );
                tmp_vec_[thread]->Zeros ( );
            }
        }

        template<class DataType>
        void PCGalerkinMatrix<DataType>::SetMatrix ( std::vector<la::CoupledMatrix<DataType>*> modes )
        {
            nmodes_ = modes.size ( );
            modes_.resize ( nmodes_ );
            for ( int i = 0; i < nmodes_; ++i )
                modes_[i] = modes[i];
        }

        template<class DataType>
        void PCGalerkinMatrix<DataType>::SetMatrixAllocate ( std::vector<la::CoupledMatrix<DataType>*> modes )
        {
            nmodes_ = modes.size ( );
            modes_.resize ( nmodes_ );
            for ( int i = 0; i < nmodes_; ++i )
            {
                modes_[i] = new la::CoupledMatrix<DataType>;
                modes_[i]->CloneFrom ( *modes[i] );
            }
            modes_allocated_ = true;
        }

        template<class DataType>
        void PCGalerkinMatrix<DataType>::SetTensor ( GalerkinTensor* pctensor )
        {
            pctensor_ = pctensor;
        }

        template<class DataType>
        void PCGalerkinMatrix<DataType>::VectorMult ( PCGalerkinVector<DataType>& in,
                                                      PCGalerkinVector<DataType>* out ) const
        {
            out->Zeros ( );

            std::vector<int> breakpoints = pctensor_->k_breakpoints ( );

#ifdef WITH_OPENMP
            if ( numthreads_ == 1 )
            {
                in.SynchronizeModes ( pctensor_ );
                omp_set_num_threads ( 1 );
                for ( int loc_idx = 0; loc_idx < pctensor_->SizeLocal ( ); ++loc_idx )
                {
                    int bp = pctensor_->L2G ( loc_idx ) + 1;
                    for ( int pos = breakpoints[bp - 1]; pos < breakpoints[bp]; ++pos )
                    {
                        modes_[pctensor_->Indices ( pos )[0]]->VectorMult ( *in.Mode ( pctensor_->Indices ( pos )[1] ), tmp_vec_[0] );
                        out->ModeAxpy ( pctensor_->Indices ( pos )[2], *tmp_vec_[0], pctensor_->Val ( pos ) );
                    }
                }
            }
            else
            {
                omp_set_num_threads ( numthreads_ );

                in.ReceiveGhost ( );
                in.SendBorder ( );

#    pragma omp parallel for
                for ( int loc_idx = 0; loc_idx < pctensor_->SizeLocal ( ); ++loc_idx )
                {
                    int bp = pctensor_->L2G ( loc_idx ) + 1;
                    for ( int pos = breakpoints[bp - 1]; pos < breakpoints[bp]; ++pos )
                    {
                        modes_[pctensor_->Indices ( pos )[0]]->diagonal ( ).VectorMult ( in.Mode ( pctensor_->Indices ( pos )[1] )->interior ( ), &( tmp_vec_[omp_get_thread_num ( )]->interior ( ) ) );
                        out->Mode ( pctensor_->Indices ( pos )[2] )->interior ( ).Axpy ( tmp_vec_[omp_get_thread_num ( )]->interior ( ), pctensor_->Val ( pos ) );
                    }
                }

                in.WaitForRecv ( );

#    pragma omp parallel for
                for ( int loc_idx = 0; loc_idx < pctensor_->SizeLocal ( ); ++loc_idx )
                {
                    int bp = pctensor_->L2G ( loc_idx ) + 1;
                    for ( int pos = breakpoints[bp - 1]; pos < breakpoints[bp]; ++pos )
                    {
                        tmp_vec_[omp_get_thread_num ( )]->Zeros ( );
                        modes_[pctensor_->Indices ( pos )[0]]->offdiagonal ( ).VectorMult ( in.Mode ( pctensor_->Indices ( pos )[1] )->ghost ( ), &( tmp_vec_[omp_get_thread_num ( )]->interior ( ) ) );
                        out->Mode ( pctensor_->Indices ( pos )[2] )->interior ( ).Axpy ( tmp_vec_[omp_get_thread_num ( )]->interior ( ), pctensor_->Val ( pos ) );
                    }
                }

                in.WaitForSend ( );
            }
#else
            in.SynchronizeModes ( pctensor_ );
            for ( int loc_idx = 0; loc_idx < pctensor_->SizeLocal ( ); ++loc_idx )
            {
                int bp = pctensor_->L2G ( loc_idx ) + 1;
                for ( int pos = breakpoints[bp - 1]; pos < breakpoints[bp]; ++pos )
                {
                    modes_[pctensor_->Indices ( pos )[0]]->VectorMult ( *in.Mode ( pctensor_->Indices ( pos )[1] ), tmp_vec_[0] );
                    out->ModeAxpy ( pctensor_->Indices ( pos )[2], *tmp_vec_[0], pctensor_->Val ( pos ) );
                }
            }
#endif
        }

        template<class DataType>
        void PCGalerkinMatrix<DataType>::Zeros ( )
        {
            for ( int i = 0; i < nmodes_; ++i )
                modes_[i]->Zeros ( );
        }

        template<class DataType>
        void PCGalerkinMatrix<DataType>::CopyFrom ( const PCGalerkinMatrix<DataType>& mat )
        {
            this->SetMatrix ( *mat.GetModes ( ) );
            this->SetTensor ( mat.GetTensor ( ) );
        }

        template<class DataType>
        std::vector<la::CoupledMatrix<DataType>*> const* PCGalerkinMatrix<DataType>::GetModes ( ) const
        {
            return &modes_;
        }

        template<class DataType>
        typename PCGalerkinMatrix<DataType>::GalerkinTensor* PCGalerkinMatrix<DataType>::GetTensor ( ) const
        {
            return pctensor_;
        }

        template class PCGalerkinMatrix<double>;

    }
}
