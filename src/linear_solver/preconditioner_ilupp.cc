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

/// \author Chandramowli Subramanian, Hendryk Bockelmann

#include "preconditioner_ilupp.h"

#include <cassert>

#include "common/log.h"
#include "linear_algebra/la_descriptor.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        PreconditionerIlupp<LAD>::PreconditionerIlupp ( )
        : PreconditionerBlockJacobi<LAD>( )
        {
        }

        template<class LAD>
        PreconditionerIlupp<LAD>::~PreconditionerIlupp ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void PreconditionerIlupp<LAD>::InitParameter ( int prepro_type, int precond_no,
                                                       int max_levels, double mem_factor,
                                                       double threshold, double min_pivot )
        {
            assert ( prepro_type >= 0 );
            assert ( precond_no >= 0 );
            assert ( max_levels >= 1 );
            assert ( mem_factor > 0. );
            assert ( threshold > 0. );
            assert ( min_pivot >= 0. );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Preprocessing Type", prepro_type );
                LOG_INFO ( "Precondition Number", precond_no );
                LOG_INFO ( "Maximum level", max_levels );
                LOG_INFO ( "Memory factor", mem_factor );
                LOG_INFO ( "Threshold", threshold );
                LOG_INFO ( "Minimal pivot", min_pivot );
            }

#ifdef WITH_ILUPP
            // set preprocessing_type
            switch ( prepro_type )
            {
                case 0: this->ilupp_preproc_.set_normalize ( );
                    break;
                case 1: this->ilupp_preproc_.set_PQ ( );
                    break;
                case 2: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING ( );
                    break;
                case 3: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_UNIT_DIAG ( );
                    break;
                case 4: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_DD_MOV_COR_IM ( );
                    break;
                case 5: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_UNIT_DIAG_DD_MOV_COR_IM ( );
                    break;
                case 6: this->ilupp_preproc_.set_SPARSE_FIRST_MAX_WEIGHTED_MATCHING_ORDERING ( );
                    break;
                case 7: this->ilupp_preproc_.set_SPARSE_FIRST_MAX_WEIGHTED_MATCHING_ORDERING_UNIT_DIAG ( );
                    break;
                case 8: this->ilupp_preproc_.set_SPARSE_FIRST_MAX_WEIGHTED_MATCHING_ORDERING_DD_MOV_COR_IM ( );
                    break;
                case 9: this->ilupp_preproc_.set_SPARSE_FIRST_MAX_WEIGHTED_MATCHING_ORDERING_UNIT_DIAG_DD_MOV_COR_IM ( );
                    break;
                case 10: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING ( );
                    break;
                case 11: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_DD_MOV_COR_IM ( );
                    break;
                case 12: this->ilupp_preproc_.set_SPARSE_FIRST ( );
                    break;
                case 13: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_PQ ( );
                    break;
                case 14: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_PQ ( );
                    break;
                case 15: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_MOVE_CORNER ( );
                    break;
                case 16: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_MOVE_CORNER ( );
                    break;
                case 17: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_MOVE_CORNER_IM ( );
                    break;
                case 18: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_MOVE_CORNER_IM ( );
                    break;
                case 19: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_WGT_MOV_COR ( );
                    break;
                case 20: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_WGT_MOV_COR ( );
                    break;
                case 21: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_WGT_MOV_COR_IM ( );
                    break;
                case 22: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_WGT_MOV_COR_IM ( );
                    break;
                case 23: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_WGT2_MOV_COR ( );
                    break;
                case 24: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_WGT2_MOV_COR ( );
                    break;
                case 25: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_WGT2_MOV_COR_IM ( );
                    break;
                case 26: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_WGT2_MOV_COR_IM ( );
                    break;
                case 27: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_SYM_PQ ( );
                    break;
                case 28: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_SYM_PQ ( );
                    break;
                case 29: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_SYMB_MOVE_CORNER ( );
                    break;
                case 30: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_SYMB_MOVE_CORNER ( );
                    break;
                case 31: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_SYMB_MOVE_CORNER_IM ( );
                    break;
                case 32: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_SYMB_MOVE_CORNER_IM ( );
                    break;
                case 33: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_SP_MOVE_CORNER ( );
                    break;
                case 34: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_SP_MOVE_CORNER ( );
                    break;
                case 35: this->ilupp_preproc_.set_MAX_WEIGHTED_MATCHING_ORDERING_SP_MOVE_CORNER_IM ( );
                    break;
                case 36: this->ilupp_preproc_.set_NORM_MAX_WEIGHTED_MATCHING_ORDERING_SP_MOVE_CORNER_IM ( );
                    break;

                default: this->ilupp_preproc_.set_normalize ( );
                    break;
            }

            // set precond_parameters for ILU++
            this->ilupp_param_.init ( this->ilupp_preproc_, precond_no, "some comment" ); // ACHTUNG: setzt default-parameters

            // set_threshold(tau) sets dropping threshold to 10^{-tau}
            this->ilupp_param_.set_threshold ( threshold );

            // set minimal pivot
            this->ilupp_param_.set_MEM_FACTOR ( mem_factor );
            this->ilupp_param_.set_MIN_PIVOT ( min_pivot );
            this->ilupp_param_.set_MAX_LEVELS ( max_levels );
#else
            printf ( "PreconditionerIlupp::InitParameter: No ILU++ support.\n" );
#endif
        }

        template<class LAD>
        void PreconditionerIlupp<LAD>::Build ( )
        {
            assert ( this->op_ != NULL );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Build Solver", 1 );
            }
            this->CreateLocalMatrix ( );
            this->Factorize ( );

            this->SetModifiedOperator ( false );
            this->SetState ( true );
        }

        template<class LAD>
        void PreconditionerIlupp<LAD>::Factorize ( )
        {
#ifdef WITH_ILUPP
            assert ( !this->ia_.empty ( ) );
            assert ( !this->ja_.empty ( ) );
            assert ( !this->val_.empty ( ) );

            this->ilupp_precond_.setup ( this->val_, this->ja_, this->ia_, iluplusplus::ROW, this->ilupp_param_ );
#else
            printf ( "PreconditionerIlupp::Factorize: No ILU++ support.\n" );
#endif
        }

        template<class LAD>
        LinearSolverState PreconditionerIlupp<LAD>::ApplyPreconditioner ( const VectorType& b, VectorType* x )
        {
#ifdef WITH_ILUPP
            assert ( this->ilupp_precond_.exists ( ) );

            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            DataType* x_local = new DataType[b.size_local ( )];
            b.GetValues ( x_local );

            this->ilupp_precond_.apply_preconditioner ( x_local, b.size_local ( ) );

            x->SetValues ( x_local );

            delete[] x_local;

            return kSolverSuccess;
#else
            printf ( "PreconditionerIlupp::ApplyPreconditioner: No ILU++ support.\n" );
            return kSolverError;
#endif
        }

        template<class LAD>
        void PreconditionerIlupp<LAD>::Clear ( )
        {
            this->ia_.clear ( ),
                    this->ja_.clear ( );
            this->val_.clear ( );
            Preconditioner<LAD>::Clear ( );
        }

        template<class LAD>
        void PreconditionerIlupp<LAD>::CreateLocalMatrix ( )
        {
            assert ( this->op_ != NULL );

            this->ia_.resize ( this->op_->nrows_local ( ) + 1 );
            this->ja_.resize ( this->op_->nnz_local_diag ( ) );
            this->val_.resize ( this->op_->nnz_local_diag ( ) );

            this->op_->ExtractDiagonalCSR ( &( this->ia_.front ( ) ), &( this->ja_.front ( ) ), &( this->val_.front ( ) ) );
        }

        /// template instantiation
        template class PreconditionerIlupp<LADescriptorCoupledD>;
        // template class PreconditionerIlupp<LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow
