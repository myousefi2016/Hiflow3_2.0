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

#include "umfpack_solver.h"
#include <cassert>

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        UmfpackSolver<LAD>::UmfpackSolver ( )
        : LinearSolver<LAD>( )
        {
            Control = new DataType[UMFPACK_CONTROL];
            Info = new DataType[UMFPACK_INFO];
            umfpack_di_defaults ( Control );
            Control[UMFPACK_STRATEGY] = 3;
            Ap_b_ = false;
            Ai_b_ = false;
            Ax_b_ = false;

            Ap_ = NULL;
            Ax_ = NULL;
            Ai_ = NULL;

            symbolic_ = NULL;
            numeric_ = NULL;

            is_symb_factorized_ = false;
            is_num_factorized_ = false;
        }

        template<class LAD>
        UmfpackSolver<LAD>::~UmfpackSolver ( )
        {
            Clear ( );
            delete[] Control;
            delete[] Info;
        }

        template<class LAD>
        void UmfpackSolver<LAD>::InitParameter ( )
        {

        }

        template<class LAD>
        void UmfpackSolver<LAD>::SetupOperator ( OperatorType& op )
        {
            bool reuse = this->reuse_;
            Clear ( );
            this->SetReuse ( reuse );

            this->op_ = &op;
            this->umf_n_ = this->op_->nrows_local ( );
            this->umf_nnz_ = this->op_->nnz_local ( );

            this->SetModifiedOperator ( true );
        }

        template<class LAD>
        void UmfpackSolver<LAD>::Build ( )
        {
            assert ( this->op_ != NULL );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Build Solver", 1 );
            }

            Ap_ = new int[umf_n_ + 1];
            Ap_b_ = true;
            Ai_ = new int[umf_nnz_];
            Ai_b_ = true;
            Ax_ = new DataType [umf_nnz_];
            Ax_b_ = true;

            for ( int k = 0; k < umf_n_ + 1; ++k )
                Ap_[k] = 0;
            for ( int k = 0; k < umf_nnz_; ++k )
            {
                Ai_[k] = 0;
                Ax_[k] = 0.0;
            }

            this->op_->ExtractCSR ( Ap_, Ai_, Ax_ );

            if ( !is_symb_factorized_ )
                this->FactorizeSymbolic ( );
            if ( !is_num_factorized_ )
                this->FactorizeNumeric ( );

            this->SetState ( true );
            this->SetModifiedOperator ( true );
        }

        template<class LAD>
        void UmfpackSolver<LAD>::FactorizeSymbolic ( )
        {
            assert ( umf_n_ > 0 );
            assert ( Ap_ != NULL );
            assert ( Ai_ != NULL );
            assert ( Ax_ != NULL );

            int status;

            status = umfpack_di_symbolic ( umf_n_, umf_n_, Ap_, Ai_, Ax_, &symbolic_, Control, Info );

            umfpack_di_report_status ( Control, status );

            assert ( status == UMFPACK_OK );
            assert ( symbolic_ != NULL );

            is_symb_factorized_ = true;
        }

        template<class LAD>
        void UmfpackSolver<LAD>::FactorizeNumeric ( )
        {
            assert ( umf_n_ > 0 );
            assert ( Ap_ != NULL );
            assert ( Ai_ != NULL );
            assert ( Ax_ != NULL );

            int status;

            status = umfpack_di_numeric ( Ap_, Ai_, Ax_, symbolic_, &numeric_, Control, Info );

            umfpack_di_report_status ( Control, status );

            assert ( status == UMFPACK_OK );
            assert ( numeric_ != NULL );

            is_num_factorized_ = true;
        }

        template<class LAD>
        LinearSolverState UmfpackSolver<LAD>::SolveFactorized ( const VectorType& b, VectorType* x )
        {
            assert ( umf_n_ > 0 );
            assert ( b.size_local ( ) == umf_n_ );
            assert ( x->size_local ( ) == umf_n_ );

            assert ( is_symb_factorized_ );
            assert ( is_num_factorized_ );

            int status;

            DataType* b_ptr = new DataType [umf_n_];
            b.GetValues ( b_ptr );
            DataType* x_ptr = new DataType [umf_n_];
            for ( int k = 0; k < umf_n_; ++k ) x_ptr[k] = 0.0;

            status = umfpack_di_solve ( UMFPACK_Aat, Ap_, Ai_, Ax_, x_ptr, b_ptr, numeric_, Control, Info );

            umfpack_di_report_status ( Control, status );
            umfpack_di_report_control ( Control );
            umfpack_di_report_info ( Control, Info );

            x->SetValues ( x_ptr );

            delete[] x_ptr;
            delete[] b_ptr;

            if ( status == UMFPACK_OK )
                return kSolverSuccess;

            return kSolverError;
        }

        template<class LAD>
        LinearSolverState UmfpackSolver<LAD>::SolveFactorized ( const std::vector<DataType>& b, std::vector<DataType>* x )
        {
            assert ( umf_n_ > 0 );
            assert ( b.size ( ) == umf_n_ );
            assert ( x->size ( ) == umf_n_ );

            assert ( is_symb_factorized_ );
            assert ( is_num_factorized_ );

            int status;

            DataType* b_ptr = new DataType [umf_n_];
            for ( int k = 0; k < umf_n_; ++k ) b_ptr[k] = b[k];
            DataType* x_ptr = new DataType [umf_n_];
            for ( int k = 0; k < umf_n_; ++k ) x_ptr[k] = 0.0;

            status = umfpack_di_solve ( UMFPACK_Aat, Ap_, Ai_, Ax_, x_ptr, b_ptr, numeric_, Control, Info );

            umfpack_di_report_status ( Control, status );
            umfpack_di_report_control ( Control );
            umfpack_di_report_info ( Control, Info );

            for ( int k = 0; k < umf_n_; ++k ) ( *x )[k] = x_ptr[k];

            delete[] x_ptr;
            delete[] b_ptr;

            if ( status == UMFPACK_OK )
                return kSolverSuccess;

            return kSolverError;
        }

        template<class LAD>
        LinearSolverState UmfpackSolver<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }
            return this->SolveFactorized ( b, x );
        }

        template<class LAD>
        LinearSolverState UmfpackSolver<LAD>::Solve ( const std::vector<DataType>& b, std::vector<DataType>* x )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }
            return this->SolveFactorized ( b, x );
        }

        template<class LAD>
        void UmfpackSolver<LAD>::Clear ( )
        {
            if ( Ap_b_ ) delete[] Ap_;
            if ( Ai_b_ ) delete[] Ai_;
            if ( Ax_b_ ) delete[] Ax_;

            Ap_b_ = false;
            Ai_b_ = false;
            Ax_b_ = false;

            if ( is_symb_factorized_ )
                umfpack_di_free_symbolic ( &symbolic_ );
            if ( is_num_factorized_ )
                umfpack_di_free_numeric ( &numeric_ );

            is_symb_factorized_ = false;
            is_num_factorized_ = false;

            symbolic_ = NULL;
            numeric_ = NULL;

            LinearSolver<LAD>::Clear ( );
        }

        template class UmfpackSolver<LADescriptorCoupledD>;
    }
}
