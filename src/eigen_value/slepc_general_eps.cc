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

/// @author Philipp Gerstner

#include "slepc_general_eps.h"
#include "slepc.h"
#include <cstdlib>
#include "common/log.h"
#include "../linear_algebra/petsc_wrapper.h"
#include "slepc_wrapper.h"

namespace hiflow
{
    namespace la
    {

        namespace slepc
        {

            /// Convert enumerator into macro name for EPS type

            const EPSType get_type ( EpsType eps_type )
            {
                switch ( eps_type )
                {
                    case POWER:
                    {
                        return EPSPOWER;
                    }
                    case SUBSPACE:
                    {
                        return EPSSUBSPACE;
                    }
                    case ARNOLDI:
                    {
                        return EPSARNOLDI;
                    }
                    case LANCZOS:
                    {
                        return EPSLANCZOS;
                    }
                    case KRYLOVSCHUR:
                    {
                        return EPSKRYLOVSCHUR;
                    }
                    case GD:
                    {
                        return EPSGD;
                    }
                    case JD:
                    {
                        return EPSJD;
                    }
                    case RQCG:
                    {
                        return EPSRQCG;
                    }
                    case ARPACK:
                    {
                        return EPSARPACK;
                    }

                    default:
                    {
                        LOG_ERROR ( "Unkown SLEPc EPSType." );
                        exit ( 1 );
                    }
                }
            }

            /// Convert enumerator into macro name for EPS problem type

            const EPSProblemType get_problem_type ( EpsProblemType eps_problem_type )
            {
                switch ( eps_problem_type )
                {
                    case HEP:
                    {
                        return EPS_HEP;
                    }
                    case NHEP:
                    {
                        return EPS_NHEP;
                    }
                    case GHEP:
                    {
                        return EPS_GHEP;
                    }
                    case GHIEP:
                    {
                        return EPS_GHIEP;
                    }
                    case GNHEP:
                    {
                        return EPS_GNHEP;
                    }
                    case PGNHEP:
                    {
                        return EPS_PGNHEP;
                    }

                    default:
                    {
                        LOG_ERROR ( "Unkown SLEPc Problem Type." );
                        exit ( 1 );
                    }
                }
            }

            /// Convert enumerator into macro name for EPS eigenvalue type

            const EPSWhich get_eigenvalue_type ( EpsWhich eigenvalue_type )
            {
                switch ( eigenvalue_type )
                {
                    case L_MAG:
                    {
                        return EPS_LARGEST_MAGNITUDE;
                    }
                    case S_MAG:
                    {
                        return EPS_SMALLEST_MAGNITUDE;
                    }
                    case L_REAL:
                    {
                        return EPS_LARGEST_REAL;
                    }
                    case S_REAL:
                    {
                        return EPS_SMALLEST_REAL;
                    }
                    case L_IMAG:
                    {
                        return EPS_LARGEST_IMAGINARY;
                    }
                    case S_IMAG:
                    {
                        return EPS_SMALLEST_IMAGINARY;
                    }

                    default:
                    {
                        LOG_ERROR ( "Unkown SLEPc Eigenvalue Type." );
                        exit ( 1 );
                    }
                }
            }

            /// Convert enumerator into macro name for EPS convergence criterion

            const EPSConv get_conv_type ( EpsConv conv_type )
            {
                switch ( conv_type )
                {
                    case ABS:
                    {
                        return EPS_CONV_ABS;
                    }
                    case REL:
                    {
                        return EPS_CONV_REL;
                    }
                    case NORM:
                    {
                        return EPS_CONV_NORM;
                    }

                    default:
                    {
                        LOG_ERROR ( "Unkown SLEPc convergence Type." );
                        exit ( 1 );
                    }
                }
            }

            /// Convert enumerator into macro name for EPS error estimation

            const EPSErrorType get_error_type ( EpsConv conv_type )
            {
                switch ( conv_type )
                {
                    case ABS:
                    {
                        return EPS_ERROR_ABSOLUTE;
                    }
                    case REL:
                    {
                        return EPS_ERROR_RELATIVE;
                    }
                    case NORM:
                    {
                        return EPS_ERROR_BACKWARD;
                    }

                    default:
                    {
                        LOG_ERROR ( "Unkown SLEPc error Type." );
                        exit ( 1 );
                    }
                }
            }

            /// Convert enumerator into macro name for EPS target 

            const EPSWhich get_target_type ( EpsMeasure meas_type )
            {
                switch ( meas_type )
                {
                    case MAG:
                    {
                        return EPS_TARGET_MAGNITUDE;
                    }
                    case REAL:
                    {
                        return EPS_TARGET_REAL;
                    }
                    case IMAG:
                    {
                        return EPS_TARGET_IMAGINARY;
                    }

                    default:
                    {
                        LOG_ERROR ( "Unkown SLEPc measure Type." );
                        exit ( 1 );
                    }
                }
            }

            /// Convert enumerator into macro name for ST Type

            const STType get_st_type ( StType type )
            {
                switch ( type )
                {
                    case SHIFT:
                    {
                        return STSHIFT;
                    }
                    case SINVERT:
                    {
                        return STSINVERT;
                    }
                    case CAYLEY:
                    {
                        return STCAYLEY;
                    }
                    case PRECOND:
                    {
                        return STPRECOND;
                    }

                    default:
                    {
                        LOG_ERROR ( "Unkown SLEPc ST Type." );
                        exit ( 1 );
                    }
                }
            }

        } // namespace petsc

        template <class LAD>
        SLEPcGeneralEPS<LAD>::SLEPcGeneralEPS ( slepc::EpsType eps_type )
        : SLEPcEigenSolver<LAD>( new slepc::EPS_wrapper ),
        eps_type_ ( slepc::KRYLOVSCHUR ),
        problem_type_ ( slepc::NHEP ),
        eval_type_ ( slepc::L_REAL ),
        conv_type_ ( slepc::ABS ),
        error_type_ ( slepc::ABS ),
        ptr_st_wrapper_ ( new slepc::ST_wrapper ),
        ptr_st_ksp_wrapper_ ( new petsc::KSP_wrapper ),
        ptr_st_pc_wrapper_ ( new petsc::PC_wrapper )
        {
            this->eps_type_ = eps_type;
        }

        template <class LAD>
        SLEPcGeneralEPS<LAD>::~SLEPcGeneralEPS ( )
        {
            SLEPcGeneralEPS<LAD>::Clear ( );
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::Init ( const MPI_Comm& comm )
        {
            if ( this->initialized_ )
                SLEPcGeneralEPS<LAD>::Clear ( );

            if ( this->comm_ != MPI_COMM_NULL )
            {
                MPI_Comm_free ( &this->comm_ );
            }

            assert ( comm != MPI_COMM_NULL );
            // MPI communicator

            // determine nb. of processes
            int info = MPI_Comm_size ( comm, &this->nb_procs_ );
            assert ( info == MPI_SUCCESS );
            assert ( this->nb_procs_ > 0 );

            // retrieve my rank
            info = MPI_Comm_rank ( comm, &this->rank_ );
            assert ( info == MPI_SUCCESS );
            assert ( this->rank_ >= 0 );
            assert ( this->rank_ < this->nb_procs_ );

            info = MPI_Comm_split ( comm, 0, this->rank_, &( this->comm_ ) );
            assert ( info == MPI_SUCCESS );

            EPSCreate ( this->comm_, &this->ptr_eps_wrapper_->eps_ );
            EPSSetType ( this->ptr_eps_wrapper_->eps_, slepc::get_type ( eps_type_ ) );
            EPSGetST ( this->ptr_eps_wrapper_->eps_, &this->ptr_st_wrapper_->st_ );
            STGetKSP ( this->ptr_st_wrapper_->st_, &this->ptr_st_ksp_wrapper_->ksp_ );
            KSPGetPC ( this->ptr_st_ksp_wrapper_->ksp_, &this->ptr_st_pc_wrapper_->pc_ );

            this->initialized_ = true;
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::InitControl ( int nev, int ncv, int mpd, int maxits, double tol, slepc::EpsConv conv_type, slepc::EpsConv error_type, bool true_res_norm )
        {

            EPSSetTolerances ( this->ptr_eps_wrapper_->eps_, tol, maxits );
            EPSSetConvergenceTest ( this->ptr_eps_wrapper_->eps_, slepc::get_conv_type ( conv_type ) );

            if ( this->eval_type_ != slepc::ALL )
            {
                if ( ncv == 0 )
                {
                    if ( mpd == 0 )
                        EPSSetDimensions ( this->ptr_eps_wrapper_->eps_, nev, PETSC_DEFAULT, PETSC_DEFAULT );
                    else
                        EPSSetDimensions ( this->ptr_eps_wrapper_->eps_, nev, PETSC_DEFAULT, mpd );
                }
                else
                {
                    if ( mpd == 0 )
                        EPSSetDimensions ( this->ptr_eps_wrapper_->eps_, nev, ncv, PETSC_DEFAULT );
                    else
                        EPSSetDimensions ( this->ptr_eps_wrapper_->eps_, nev, ncv, mpd );
                }
            }

            this->nev_ = nev;
            this->ncv_ = ncv;
            this->mpd_ = mpd;
            this->error_type_ = error_type;
            this->conv_type_ = conv_type;
            this->control_set_ = true;

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "Number of Eigenvalues to compute: ", nev );
                LOG_INFO ( "Maximum size of subspace:         ", ncv );
                LOG_INFO ( "Maximum size of projected matrix: ", mpd );
                LOG_INFO ( "Maximum number of iterations:     ", maxits );
                LOG_INFO ( "Tolerance:                        ", tol );
                LOG_INFO ( "Convergence test type: 			", conv_type );
                LOG_INFO ( "Error estimation type:			", error_type );
                LOG_INFO ( "Compute true residual:			", true_res_norm );
            }
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::PassOperatorsToSolver ( )
        {
            assert ( this->mat_A_set_ );

            if ( this->mat_A_set_ && this->mat_B_set_ )
                EPSSetOperators ( this->ptr_eps_wrapper_->eps_, this->ptr_A_wrapper_->mat_, this->ptr_B_wrapper_->mat_ );

            if ( this->mat_A_set_ && !this->mat_B_set_ )
                EPSSetOperators ( this->ptr_eps_wrapper_->eps_, this->ptr_A_wrapper_->mat_, NULL );

            this->operator_set_ = true;
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::SetProblemType ( slepc::EpsProblemType type )
        {
            EPSSetProblemType ( this->ptr_eps_wrapper_->eps_, slepc::get_problem_type ( type ) );
            this->problem_set_ = true;
            if ( this->print_level_ > 0 )
                LOG_INFO ( "Eigenvalue problem type: 			", type );
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::ComputeExtremalEigenValue ( slepc::EpsWhich type )
        {
            EPSSetWhichEigenpairs ( this->ptr_eps_wrapper_->eps_, slepc::get_eigenvalue_type ( type ) );
            this->eval_type_ = type;
            this->eigenvalue_set_ = true;
            if ( this->print_level_ > 0 )
                LOG_INFO ( "Compute extremal Eigenvalue: 		", type );
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::ComputeTargetEigenValue ( DataType target, slepc::EpsMeasure type )
        {
            EPSSetWhichEigenpairs ( this->ptr_eps_wrapper_->eps_, slepc::get_target_type ( type ) );
            EPSSetTarget ( this->ptr_eps_wrapper_->eps_, target );
            this->eval_type_ = slepc::TARGET;
            this->eigenvalue_set_ = true;
            if ( this->print_level_ > 0 )
                LOG_INFO ( "Compute target Eigenvalue: 		", target );
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::ComputeIntervalEigenValue ( DataType a, DataType b )
        {
            EPSSetWhichEigenpairs ( this->ptr_eps_wrapper_->eps_, EPS_ALL );
            EPSSetInterval ( this->ptr_eps_wrapper_->eps_, a, b );
            this->eval_type_ = slepc::ALL;
            this->eigenvalue_set_ = true;

            if ( this->ncv_ == 0 )
            {
                if ( this->mpd_ == 0 )
                    EPSKrylovSchurSetDimensions ( this->ptr_eps_wrapper_->eps_, this->nev_, PETSC_DEFAULT, PETSC_DEFAULT );
                else
                    EPSKrylovSchurSetDimensions ( this->ptr_eps_wrapper_->eps_, this->nev_, PETSC_DEFAULT, this->mpd_ );
            }
            else
            {
                if ( this->mpd_ == 0 )
                    EPSKrylovSchurSetDimensions ( this->ptr_eps_wrapper_->eps_, this->nev_, this->ncv_, PETSC_DEFAULT );
                else
                    EPSKrylovSchurSetDimensions ( this->ptr_eps_wrapper_->eps_, this->nev_, this->ncv_, this->mpd_ );
            }

            if ( this->print_level_ > 0 )
            {
                LOG_INFO ( "Compute interval Eigenvalue: 		", a );
                LOG_INFO ( "Compute interval Eigenvalue: 		", b );
            }
        }

        // TODO check whether transposed operator is needed

        template <class LAD>
        int SLEPcGeneralEPS<LAD>::CheckSolver ( )
        {
            if ( this->eps_type_ == slepc::POWER )
            {
                if ( this->eval_type_ != slepc::L_MAG )
                {
                    LOG_ERROR ( "Eigensolver: Power method only supports EigenvalueType L_MAG" );
                    return -1;
                }
            }
            if ( this->eps_type_ == slepc::SUBSPACE )
            {
                if ( this->eval_type_ != slepc::L_MAG )
                {
                    LOG_ERROR ( "Eigensolver: SubSpace method only supports EigenvalueType L_MAG" );
                    return -1;
                }
            }
            if ( this->eps_type_ == slepc::LANCZOS )
            {
                if ( this->problem_type_ != slepc::HEP && this->problem_type_ != slepc::GHEP )
                {
                    LOG_ERROR ( "Eigensolver: Lanczos method only supports problem type HEP and GHEP" );
                    return -1;
                }
                /*
                if (this->mat_type == petsc::MPI_MATFREE) 
                {
                        if ( !this->mat_At_set_ ) 
                        {
                                LOG_ERROR("Eigensolver: Lanczos method requires transposed operator A");
                                return -1;
                        }
                        if ( this->problem_type != HEP && this->problem_type != NHEP ) 
                        {
                                if ( !this->mat_Bt_set_ ) 
                                {	
                                        LOG_ERROR("Eigensolver: Lanczos method requires transposed operator B");
                                        return -1;
                                }
                        }
                }
                 */
            }
            if ( this->eps_type_ == slepc::RQCG )
            {
                if ( this->eval_type_ != slepc::S_REAL )
                {
                    LOG_ERROR ( "Eigensolver: RQCG method only supports EigenvalueType S_REAL" );
                    return -1;
                }
                if ( this->problem_type_ != slepc::HEP && this->problem_type_ != slepc::GHEP )
                {
                    LOG_ERROR ( "Eigensolver: RQCG method only supports problem type HEP and GHEP" );
                    return -1;
                }
            }
            return 0;
        }

        template <class LAD>
        EigenSolverState SLEPcGeneralEPS<LAD>::Solve ( )
        {

            if ( !this->control_set_ )
            {
                LOG_ERROR ( "Eigensolver: InitControl not called" );
                return EigenSolverError;
            }
            if ( !this->operator_set_ )
            {
                LOG_ERROR ( "Eigensolver: SetupOperators not called" );
                return EigenSolverError;
            }
            if ( !this->eigenvalue_set_ )
            {
                LOG_ERROR ( "Eigensolver: Desired Eigenvalues not specified " );
                return EigenSolverError;
            }
            if ( !this->problem_set_ )
            {
                LOG_ERROR ( "Eigensolver SetProblem not called" );
                return EigenSolverError;
            }

            // Check whether problem type is supported by chosen solver 
            if ( this->CheckSolver ( ) != 0 )
            {
                return EigenSolverWrong;
            }

            // Solve eigenvalue problem 
            PetscErrorCode ierr;
            ierr = EPSSolve ( this->ptr_eps_wrapper_->eps_ );

            // Get number of converged eigenvalues and iterations 
            EPSGetConverged ( this->ptr_eps_wrapper_->eps_, &this->n_conv_ );
            EPSGetIterationNumber ( this->ptr_eps_wrapper_->eps_, &this->iter_ );

            // Get residuals 
            this->res_.resize ( this->n_conv_, 0. );
            for ( int j = 0; j<this->n_conv_; j++ )
            {
                EPSComputeError ( this->ptr_eps_wrapper_->eps_, j, slepc::get_error_type ( this->error_type_ ), &this->res_[j] );
            }

            // Get eigenvalues
            this->eig_val_.resize ( this->n_conv_ );
#ifdef WITH_COMPLEX_PETSC
            this->ptr_v_wrapper_.resize ( this->n_conv_ );
#else 
            this->ptr_vr_wrapper_.resize ( this->n_conv_ );
            this->ptr_vi_wrapper_.resize ( this->n_conv_ );
#endif

            for ( int j = 0; j<this->n_conv_; j++ )
            {
#ifdef WITH_COMPLEX_PETSC
                this->ptr_v_wrapper_[j] = new (petsc::Vec_wrapper );
                MatCreateVecs ( this->ptr_A_wrapper_->mat_, NULL, &this->ptr_v_wrapper_[j]->vec_ );
                EPSGetEigenpair ( this->ptr_eps_wrapper_->eps_, j, &this->eig_val_[j].val_, NULL, this->ptr_v_wrapper_[j]->vec_, NULL );

                this->eig_val_[j].real_ = PetscRealPart ( this->eig_val_[j].val_ );
                this->eig_val_[j].imag_ = PetscImaginaryPart ( this->eig_val_[j].val_ );
#else  
                this->ptr_vr_wrapper_[j] = new (petsc::Vec_wrapper );
                this->ptr_vi_wrapper_[j] = new (petsc::Vec_wrapper );
                MatCreateVecs ( this->ptr_A_wrapper_->mat_, NULL, &this->ptr_vr_wrapper_[j]->vec_ );
                MatCreateVecs ( this->ptr_A_wrapper_->mat_, NULL, &this->ptr_vi_wrapper_[j]->vec_ );

                EPSGetEigenpair ( this->ptr_eps_wrapper_->eps_, j, &this->eig_val_[j].real_, &this->eig_val_[j].imag_, this->ptr_vr_wrapper_[j]->vec_, this->ptr_vi_wrapper_[j]->vec_ );
#endif  
                if ( this->print_level_ > 0 )
                {
                    LOG_INFO ( "Computed Eigenvalue real part:	", this->eig_val_[j].real_ );
                    LOG_INFO ( "Computed Eigenvalue imag part:	", this->eig_val_[j].imag_ );
                    LOG_INFO ( "Computed Eigenvalue residual:		", this->res_[j] );
                }
            }

            if ( this->n_conv_ > 0 )
                return EigenSolverSuccess;
            else
                return EigenSolverExceeded;
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::SetSTShift ( double shift )
        {
            STSetShift ( this->ptr_st_wrapper_->st_, shift );
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::SetSTType ( slepc::StType type )
        {
            STSetType ( this->ptr_st_wrapper_->st_, slepc::get_st_type ( type ) );
        }

        template <class LAD>
        petsc::KSP_wrapper* SLEPcGeneralEPS<LAD>::GetSTKSP ( )
        {
            return this->ptr_st_ksp_wrapper_;
        }

        template <class LAD>
        petsc::PC_wrapper* SLEPcGeneralEPS<LAD>::GetSTPC ( )
        {
            return this->ptr_st_pc_wrapper_;
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::SetSTKSP ( slepc::KSPType ksp_type, int max_iter, int max_size, double abs_tol, double rel_tol )
        {
            double div_tol = 1e6;
            KSPSetTolerances ( this->ptr_st_ksp_wrapper_->ksp_, rel_tol, abs_tol, div_tol, max_iter );

            switch ( ksp_type )
            {
                case (slepc::GMRES ):
                    KSPSetType ( this->ptr_st_ksp_wrapper_->ksp_, KSPGMRES );
                    KSPGMRESSetRestart ( this->ptr_st_ksp_wrapper_->ksp_, max_size );
                    break;
                case (slepc::FGMRES ):
                    KSPSetType ( this->ptr_st_ksp_wrapper_->ksp_, KSPFGMRES );
                    KSPGMRESSetRestart ( this->ptr_st_ksp_wrapper_->ksp_, max_size );
                    break;
                case (slepc::CG ):
                    KSPSetType ( this->ptr_st_ksp_wrapper_->ksp_, KSPCG );
                    break;
                case (slepc::PREONLY ):
                    KSPSetType ( this->ptr_st_ksp_wrapper_->ksp_, KSPPREONLY );
                    break;
            }
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::SetSTPC ( slepc::PcType pc_type, slepc::SolverPackage package )
        {
            switch ( pc_type )
            {
                case (slepc::NONE ):
                    PCSetType ( this->ptr_st_pc_wrapper_->pc_, PCNONE );
                    break;
                case (slepc::LU ):
                    PCSetType ( this->ptr_st_pc_wrapper_->pc_, PCLU );
                    switch ( package )
                    {
                        case (slepc::MUMPS ):
                            PCFactorSetMatSolverPackage ( this->ptr_st_pc_wrapper_->pc_, MATSOLVERMUMPS );
                            break;
                        case (slepc::UMFPACK ):
                            PCFactorSetMatSolverPackage ( this->ptr_st_pc_wrapper_->pc_, MATSOLVERUMFPACK );
                            break;
                        case (slepc::NOPACKAGE ):
                            break;
                    }
                    break;
                case (slepc::ILU ):
                    PCSetType ( this->ptr_st_pc_wrapper_->pc_, PCILU );
                    break;
                case (slepc::JACOBI ):
                    PCSetType ( this->ptr_st_pc_wrapper_->pc_, PCJACOBI );
                    break;
                case (slepc::BJACOBI ):
                    PCSetType ( this->ptr_st_pc_wrapper_->pc_, PCBJACOBI );
                    break;
            }
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::SetupSTSolver ( )
        {
            PCSetUp ( this->ptr_st_pc_wrapper_->pc_ );
            KSPSetUp ( this->ptr_st_ksp_wrapper_->ksp_ );
        }

        template <class LAD>
        void SLEPcGeneralEPS<LAD>::Clear ( )
        {
            if ( this->initialized_ )
            {
                EPSDestroy ( &this->ptr_eps_wrapper_->eps_ );
                this->initialized_ = false;
            }
            this->Clear ( );
        }

#ifdef WITH_PETSC   
        template class SLEPcGeneralEPS<LADescriptorPETSc>;
#endif
#ifdef WITH_HYPRE
        template class SLEPcGeneralEPS<LADescriptorHypreD>;
#endif
        template class SLEPcGeneralEPS<LADescriptorCoupledD>;

    } // namespace la
} // namespace hiflow
