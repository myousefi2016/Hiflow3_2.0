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

#include "GPUcublas2_CSR_lMatrix.h"
#include "../lmp_log.h"

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        GPUcublas2_CSR_lMatrix<ValueType>::GPUcublas2_CSR_lMatrix ( )
        {
#ifdef WITH_CUDA
            this->implementation_name_ = "cublas2";
            this->implementation_id_ = CUBLAS2;
            cusparseCreate ( &this->handle_ );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        GPUcublas2_CSR_lMatrix<ValueType>::GPUcublas2_CSR_lMatrix (
                                                                    int init_nnz, int init_num_row, int init_num_col,
                                                                    std::string const &init_name )
        {
#ifdef WITH_CUDA
            this->Init ( init_nnz, init_num_row, init_num_col, init_name );
            this->implementation_name_ = "cublas2";
            this->implementation_id_ = CUBLAS2;
            cusparseCreate ( &this->handle_ );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        GPUcublas2_CSR_lMatrix<ValueType>::~GPUcublas2_CSR_lMatrix ( )
        {
#ifdef WITH_CUDA
            cusparseDestroy ( this->handle_ );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        void GPUcublas2_CSR_lMatrix<ValueType>::VectorMult (
                                                             const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
        {
#ifdef WITH_CUDA
            const GPU_lVector<ValueType> *casted_invec =
                    dynamic_cast < const GPU_lVector<ValueType> * > ( &invec );
            GPU_lVector<ValueType> *casted_outvec =
                    dynamic_cast < GPU_lVector<ValueType> * > ( outvec );

            if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
            {
                LOG_ERROR (
                            "ERROR GPUcublas2_CSR_lMatrix<ValueType>::VectorMult unsupported in or "
                            "out vector" );
                this->print ( );
                invec.print ( );
                outvec->print ( );
                exit ( -1 );
            }

            this->VectorMult ( *casted_invec, casted_outvec );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_CSR_lMatrix<double>::VectorMult (
                                                          const GPU_lVector<double> &invec, GPU_lVector<double> *outvec ) const
        {
#ifdef WITH_CUDA
            // Init matrix properties
            cusparseMatDescr_t descr = 0;
            cusparseCreateMatDescr ( &descr );
            cusparseSetMatType ( descr, CUSPARSE_MATRIX_TYPE_GENERAL );
            cusparseSetMatIndexBase ( descr, CUSPARSE_INDEX_BASE_ZERO );

            const double alpha = 1.0;
            const double beta = 0.0;

            /* cublasStatus_t status */ cusparseDcsrmv (
                                                         this->handle_, // cusparseHandle_ handle
                                                         CUSPARSE_OPERATION_NON_TRANSPOSE, // cusparseOperation_t transa
                                                         this->num_col_, // int m
                                                         this->num_row_, // int n
                                                         this->nnz_, // int nnz
                                                         &alpha, // const double *alpha
                                                         descr, // const cusparseMatDescr_t descrA
                                                         this->matrix.val, // const double *csrValA
                                                         this->matrix.row, // const int *csrRowPtrA
                                                         this->matrix.col, // const int *csrColIndA
                                                         invec.buffer, // const double *x
                                                         &beta, // const double *beta
                                                         outvec->buffer // double *y
                                                         );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <>
        void GPUcublas2_CSR_lMatrix<float>::VectorMult (
                                                         const GPU_lVector<float> &invec, GPU_lVector<float> *outvec ) const
        {
#ifdef WITH_CUDA
            // Init matrix properties
            cusparseMatDescr_t descr = 0;
            cusparseCreateMatDescr ( &descr );
            cusparseSetMatType ( descr, CUSPARSE_MATRIX_TYPE_GENERAL );
            cusparseSetMatIndexBase ( descr, CUSPARSE_INDEX_BASE_ZERO );

            const float alpha = 1.0;
            const float beta = 0.0;

            /* cublasStatus_t status */ cusparseScsrmv (
                                                         this->handle_, // cusparseHandle_ handle
                                                         CUSPARSE_OPERATION_NON_TRANSPOSE, // cusparseOperation_t transa
                                                         this->num_col_, // int m
                                                         this->num_row_, // int n
                                                         this->nnz_, // int nnz
                                                         &alpha, // const double *alpha
                                                         descr, // const cusparseMatDescr_t descrA
                                                         this->matrix.val, // const double *csrValA
                                                         this->matrix.row, // const int *csrRowPtrA
                                                         this->matrix.col, // const int *csrColIndA
                                                         invec.buffer, // const double *x
                                                         &beta, // const double *beta
                                                         outvec->buffer // double *y
                                                         );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        void GPUcublas2_CSR_lMatrix<ValueType>::VectorMultAdd (
                                                                const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
        {
#ifdef WITH_CUDA
            const GPU_lVector<ValueType> *casted_invec =
                    dynamic_cast < const GPU_lVector<ValueType> * > ( &invec );
            GPU_lVector<ValueType> *casted_outvec =
                    dynamic_cast < GPU_lVector<ValueType> * > ( outvec );

            if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
            {
                LOG_ERROR (
                            "ERROR GPUcublas2_CSR_lMatrix<ValueType>::VectorMultAdd unsupported in "
                            "or out vector" );
                this->print ( );
                invec.print ( );
                outvec->print ( );
                exit ( -1 );
            }

            this->VectorMultAdd ( *casted_invec, casted_outvec );
#else
            LOG_ERROR ( "No CUDA support" );
            exit ( -1 );
#endif
        }

        template <typename ValueType>
        void GPUcublas2_CSR_lMatrix<ValueType>::VectorMultAdd (
                                                                const GPU_lVector<ValueType> &invec, GPU_lVector<ValueType> *outvec ) const
        {
            LOG_ERROR ( "ERROR GPUcublas2_CSR_lMatrix<ValueType>::VectorMultAdd not implemented yet." );
        }

        template <typename ValueType>
        void GPUcublas2_CSR_lMatrix<ValueType>::VectorMultNoDiag (
                                                                   const lVector<ValueType> &in, lVector<ValueType> *out ) const
        {
            LOG_ERROR ( "ERROR GPUcublas2_CSR_lMatrix<ValueType>::VectorMultNoDiag not implemented yet." );
        }

        template class GPUcublas2_CSR_lMatrix<float>;
        template class GPUcublas2_CSR_lMatrix<double>;

    } // namespace la
} // namespace hiflow
