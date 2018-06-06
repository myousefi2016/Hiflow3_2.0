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

#include "pce_matrix.h"

namespace hiflow
{
    namespace la
    {

        // constructor

        template<class DataType>
        PCEMatrix<DataType>::PCEMatrix ( )
        {

            this->basis_matrix_.clear ( );
            this->nbasis_ = -1;

            // continue filling
        }

        // destructor

        template<class DataType>
        PCEMatrix<DataType>::~PCEMatrix ( )
        {

            this->Clear ( );
            basis_matrix_.clear ( );

        }

        // Inititialize

        template<class DataType>
        void PCEMatrix<DataType>::Init ( PCTensor& pctensor, const MPI_Comm& comm, const LaCouplings& cp )
        {
            // assign pctensor
            this->pctensor_ = pctensor;

            // calculate size of basis
            std::vector<int> bk = pctensor_.k_breakpoints ( );
            this->nbasis_ = bk[1] - bk[0];
            assert ( this->nbasis_ > 0 );

            // initialize the size of basis_matrix_
            this->basis_matrix_.resize ( this->nbasis_ );

            // initialize each HypreMatrix in basis_matrix_
            // currently, only use the mean_matrix to initialize the set of matrices,
            // it can be extended with a set of matrices by further implementation
            for ( int i = 0; i != this->nbasis_; ++i )
            {
                this->basis_matrix_[i].Init ( comm, cp );
            }

        }

        // accessing the member of basis_matrix_

        template<class DataType>
        HypreMatrix<DataType>& PCEMatrix<DataType>::BasisMode ( const int& i )
        {
            return this->basis_matrix_[i];
        }

        template<class DataType>
        const HypreMatrix<DataType>& PCEMatrix<DataType>::GetBasisMode ( const int& i ) const
        {
            return this->basis_matrix_[i];
        }

        // number of basis

        template<class DataType>
        int PCEMatrix<DataType>::nb_basis ( ) const
        {
            assert ( this->nbasis_ >= 0 );
            return this->nbasis_;
        }

        // Zeros

        template<class DataType>
        void PCEMatrix<DataType>::Zeros ( )
        {
            assert ( this->nbasis_ > 0 );
            for ( int i = 0; i != this->nbasis_; ++i )
            {
                this->basis_matrix_[i].Zeros ( );
            }
        }

        template<class DataType>
        void PCEMatrix<DataType>::Zeros ( const int& i )
        {
            assert ( this->nbasis_ > 0 );
            this->basis_matrix_[i].Zeros ( );
        }

        // Clear

        template<class DataType>
        void PCEMatrix<DataType>::Clear ( )
        {
            assert ( this->nbasis_ > 0 );
            for ( int i = 0; i != this->nbasis_; ++i )
            {
                this->basis_matrix_[i].Clear ( );
            }
        }

        // VectorMult

        template<class DataType>
        void PCEMatrix<DataType>::VectorMult ( PCEVector<DataType>& in, PCEVector<DataType> *out ) const
        {
            assert ( this->nbasis_ > 0 );
            assert ( in.nb_mode ( ) == out->nb_mode ( ) );
            assert ( this->pctensor_.Size ( ) == in.nb_mode ( ) );

            out->Zeros ( );

            std::vector<int> bk = this->pctensor_.k_breakpoints ( );

            HypreVector<DataType> vec_tmp;
            vec_tmp.CloneFromWithoutContent ( in.Mode ( 0 ) );
            vec_tmp.Zeros ( );

            for ( int mode = 0; mode != this->pctensor_.Size ( ); ++mode )
            {

                for ( int pos = bk[mode]; pos != bk[mode + 1]; ++pos )
                {

                    vec_tmp.Zeros ( );

                    std::vector<int> mode_idx = this->pctensor_.IndicesGlobal ( pos );

                    this->VectorMult ( mode_idx[0], in.Mode ( mode_idx[1] ), &vec_tmp );

                    out->Mode ( mode ).Axpy ( vec_tmp, this->pctensor_.Val ( pos ) );

                }

            }

        }

        template<class DataType>
        void PCEMatrix<DataType>::VectorMult ( PCEVector<DataType>& in, PCEVector<DataType> *out, const int& l ) const
        {
            assert ( this->nbasis_ > 0 );
            assert ( in.nb_mode ( ) == out->nb_mode ( ) );

            out->Zeros ( );

            const int mode_total = this->pctensor_.Size ( l );

            assert ( in.nb_mode ( ) == mode_total );
            assert ( out->nb_mode ( ) == mode_total );

            std::vector<int> bk = this->pctensor_.k_breakpoints ( );

            HypreVector<DataType> vec_tmp;
            vec_tmp.CloneFromWithoutContent ( in.Mode ( 0 ) );
            vec_tmp.Zeros ( );

            for ( int mode = 0; mode != this->pctensor_.Size ( ); ++mode )
            {

                for ( int pos = bk[mode]; pos != bk[mode + 1]; ++pos )
                {

                    std::vector<int> mode_idx = this->pctensor_.IndicesGlobal ( pos );

                    if ( mode_idx[0] < mode_total && mode_idx[1] < mode_total && mode_idx[2] < mode_total )
                    {

                        vec_tmp.Zeros ( );

                        this->VectorMult ( mode_idx[0], in.Mode ( mode_idx[1] ), &vec_tmp );

                        out->Mode ( mode ).Axpy ( vec_tmp, this->pctensor_.Val ( pos ) );

                    }

                }

            }

        }

        template<class DataType>
        void PCEMatrix<DataType>::VectorMult ( const int& i, HypreVector<DataType>& in, HypreVector<DataType>* out ) const
        {
            assert ( this->nbasis_ > 0 );
            assert ( this->nbasis_ > i );

            this->basis_matrix_[i].VectorMult ( in, out );
        }

        template class PCEMatrix<double>;

    }
}
