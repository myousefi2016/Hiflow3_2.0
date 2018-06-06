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
/// @date 2015-12-21

#include "petsc.h"
#include "petsc_general_ksp.h"
#include "../linear_algebra/petsc_matrix_interface.h"
#include "../linear_algebra/petsc_vector_interface.h"

namespace hiflow
{
    namespace la
    {
        namespace petsc
        {

            struct KSP_wrapper
            {
                ::KSP ksp_;
            };

            /// Convert enumerator into macro name for KSP type

            const char* get_type ( KSPType ksp_type )
            {
                switch ( ksp_type )
                {
                    case CG:
                    {
                        return KSPCG;
                    }
                    case GMRES:
                    {
                        return KSPGMRES;
                    }
                    default:
                    {
                        LOG_ERROR ( "Unkown PETSc KSPType." );
                        exit ( 1 );
                    }
                }
            }

        } // namespace petsc

        template <class LAD>
        PETScGeneralKSP<LAD>::PETScGeneralKSP ( const MPI_Comm& comm, const OperatorType& op, petsc::KSPType ksp_type )
        : initialized_ ( true ),
        comm_ ( comm ),
        op_ ( op ),
        PETScLinearSolver<LAD>( new petsc::KSP_wrapper )
        {
            KSPCreate ( comm, &this->ptr_ksp_wrapper_->ksp_ );
            op_.Assembly ( );
            KSPCreate ( comm, &this->ptr_ksp_wrapper_->ksp_ );
            KSPSetType ( this->ptr_ksp_wrapper_->ksp_, petsc::get_type ( ksp_type ) );
            KSPSetInitialGuessNonzero ( this->ptr_ksp_wrapper_->ksp_, PETSC_TRUE );
            KSPSetOperators ( this->ptr_ksp_wrapper_->ksp_, op_.ptr_mat_wrapper_->mat_, op_.ptr_mat_wrapper_->mat_ );
            KSPSetUp ( this->ptr_ksp_wrapper_->ksp_ );
        }

        template <class LAD>
        PETScGeneralKSP<LAD>::~PETScGeneralKSP ( )
        {
            Clear ( );
        }

        template <class LAD>
        void PETScGeneralKSP<LAD>::InitControl ( int maxits, double abstol, double reltol, double divtol )
        {
            KSPSetTolerances ( this->ptr_ksp_wrapper_->ksp_, reltol, abstol, divtol, maxits );
        }

        template <class LAD>
        void PETScGeneralKSP<LAD>::InitControl ( int maxits, DataType abstol, DataType reltol )
        {
            DataType cur_reltol, cur_abstol, cur_dtol;
            int cur_maxits;
            KSPGetTolerances ( this->ptr_ksp_wrapper_->ksp_, &cur_reltol, &cur_abstol, &cur_dtol, &cur_maxits );
            KSPSetTolerances ( this->ptr_ksp_wrapper_->ksp_, reltol, abstol, cur_dtol, maxits );
        }

        template <class LAD>
        void PETScGeneralKSP<LAD>::SetRelativeTolerance ( DataType reltol )
        {
            DataType cur_reltol, cur_abstol, cur_dtol;
            int cur_maxits;
            KSPGetTolerances ( this->ptr_ksp_wrapper_->ksp_, &cur_reltol, &cur_abstol, &cur_dtol, &cur_maxits );
            KSPSetTolerances ( this->ptr_ksp_wrapper_->ksp_, reltol, cur_abstol, cur_dtol, cur_maxits );
        }

        template <class LAD>
        void PETScGeneralKSP<LAD>::SetupOperator ( OperatorType& op )
        {
            KSPSetOperators ( this->ptr_ksp_wrapper_->ksp_, op_.ptr_mat_wrapper_->mat_, op_.ptr_mat_wrapper_->mat_ );
        }

        template <class LAD>
        LinearSolverState PETScGeneralKSP<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
            PetscErrorCode ierr = KSPSolve ( this->ptr_ksp_wrapper_->ksp_, b.ptr_vec_wrapper_->vec_, x->ptr_vec_wrapper_->vec_ );
            return ierr ? kSolverError : kSolverSuccess;
        }

        template <class LAD>
        void PETScGeneralKSP<LAD>::Clear ( )
        {
            if ( initialized_ )
            {
                KSPDestroy ( &this->ptr_ksp_wrapper_->ksp_ );
                initialized_ = false;
            }
        }

        template class PETScGeneralKSP<LADescriptorPETSc>;

    } // namespace la
} // namespace hiflow
