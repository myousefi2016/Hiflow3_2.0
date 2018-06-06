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

#include "polynomial_chaos/pc_multilevel.h"

#include <cassert>
#include <cmath>

#include "common/log.h"
#include "common/macros.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class LAD>
        PCMultilevelSolver<LAD>::PCMultilevelSolver ( const PropertyTree &config )
        : la::LinearSolver<LAD>( ),
        config_ ( &config )
        {
            config["Multilevel"]["Nu1"].read<int>( nu1_ );
            config["Multilevel"]["Nu2"].read<int>( nu2_ );
            config["Multilevel"]["Mu"].read<int>( mu_ );
            config["Multilevel"]["Smoothing"].read<std::string>( smoother_type_ );
            config["Multilevel"]["MLType"].read<std::string>( mltype_ );
            MPI_Comm_rank ( MPI_COMM_WORLD, &global_mpi_rank_ );
            this->SetModifiedSmoother ( false );
            this->smoother_ = NULL;
            this->res_assembler_ = NULL;
            this->mean_solver_ = NULL;
            this->pctensor_ = NULL;
        }

        template<class LAD>
        PCMultilevelSolver<LAD>::~PCMultilevelSolver ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::Clear ( )
        {
            la::LinearSolver<LAD>::Clear ( );

            this->SetModifiedSmoother ( false );
            this->smoother_ = NULL;
            this->res_assembler_ = NULL;
            this->mean_solver_ = NULL;
            this->mean_solver_cg_.Clear ( );
            this->mean_solver_gmres_.Clear ( );
            this->pctensor_ = NULL;
            this->global_mpi_rank_ = -1;
            this->nu1_ = 0;
            this->nu2_ = 0;
            this->mu_ = 0;
            this->smoother_type_ = "";
            this->mltype_ = "Free";
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::SetupOperator ( OperatorType& op )
        {
            this->op_ = &op;
            SetPCTensor ( this->op_->GetTensor ( ) );
            smoother_ = this->op_->GetModes ( )->at ( 0 );
            if ( mltype_ == "Free" ) std::cout << "Warning, using Multilevel with matrix in matrix free mode\n";

            this->SetModifiedOperator ( true );
            this->SetModifiedSmoother ( true );
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::SetupSmoother ( la::CoupledMatrix<DataType>& smoother )
        {
            smoother_ = &smoother;
            this->SetModifiedSmoother ( true );
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::SetupResidualAssembler ( ResidualAssembler<LAD>& assembler )
        {
            res_assembler_ = &assembler;
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::Build ( )
        {
            assert ( this->smoother_ != NULL );

            if ( !this->GetModifiedSmoother ( ) )
            {
                this->SetModifiedOperator ( false );
                this->SetModifiedSmoother ( false );
                this->SetState ( true );
                return;
            }

            if ( smoother_type_ == "Umfpack" )
            {
#ifdef WITH_UMFPACK
                mean_solver_umf_.SetupOperator ( *smoother_ );
                mean_solver_umf_.Build ( );
                mean_solver_ = &mean_solver_umf_;
#else
                interminable_assert ( 0 );
#endif
            }
            else if ( smoother_type_ == "PGMRES" )
            {
#ifdef WITH_ILUPP
                ilupp_.InitParameter ( 0, 11, 20, 0.8, 1.0, 0.05 );
                ilupp_.SetupOperator ( *smoother_ );

                mean_solver_gmres_.InitControl ( 1000, 1.e-14, 1.e-1, 1.e6 );
                mean_solver_gmres_.InitParameter ( 50, "RightPreconditioning" );
                mean_solver_gmres_.SetupOperator ( *smoother_ );
                mean_solver_gmres_.SetupPreconditioner ( ilupp_ );
                mean_solver_gmres_.Build ( );
                mean_solver_ = &mean_solver_gmres_;
#else
                interminable_assert ( 0 );
#endif
            }
            else if ( smoother_type_ == "Inexact_GMRES" )
            {
                mean_solver_gmres_.InitControl ( 1000, 1.e-14, 1.e-1, 1.e6 );
                mean_solver_gmres_.InitParameter ( 50, "NoPreconditioning" );
                mean_solver_gmres_.SetupOperator ( *smoother_ );
                mean_solver_gmres_.Build ( );
                mean_solver_ = &mean_solver_gmres_;
            }
            else if ( smoother_type_ == "GMRES" )
            {
                mean_solver_gmres_.InitControl ( 1000, 1.e-14, 1.e-12, 1.e6 );
                mean_solver_gmres_.InitParameter ( 50, "NoPreconditioning" );
                mean_solver_gmres_.SetupOperator ( *smoother_ );
                mean_solver_gmres_.Build ( );
                mean_solver_ = &mean_solver_gmres_;
            }
            else if ( smoother_type_ == "InexactPCG" )
            {
                det_precond_.Init_SSOR ( 1.0 );
                det_precond_.SetupOperator ( *smoother_ );

                mean_solver_cg_.InitControl ( 1000, 1.e-14, 1.e-1, 1.e6 );
                mean_solver_cg_.InitParameter ( "Preconditioning" );
                mean_solver_cg_.SetupOperator ( *smoother_ );
                mean_solver_cg_.SetupPreconditioner ( det_precond_ );
                mean_solver_cg_.Build ( );
                mean_solver_ = &mean_solver_cg_;
            }
            else if ( smoother_type_ == "InexactCG" )
            {
                mean_solver_cg_.InitControl ( 1000, 1.e-14, 1.e-1, 1.e6 );
                mean_solver_cg_.InitParameter ( "NoPreconditioning" );
                mean_solver_cg_.SetupOperator ( *smoother_ );
                mean_solver_cg_.Build ( );
                mean_solver_ = &mean_solver_cg_;
            }
            else
            {
                std::cout << "Unkown Smoother type for PC Multilevel";
                exit ( -1 );
            }
            this->SetModifiedOperator ( false );
            this->SetModifiedSmoother ( false );
            this->SetState ( true );
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::UpdateControl ( int maxit, double atol, double rtol, double dtol )
        {
            mean_solver_->InitControl ( maxit, atol, rtol, dtol );
        }

        template<class LAD>
        la::LinearSolverState PCMultilevelSolver<LAD>::ApplyPreconditioner ( const VectorType& b,
                                                                             VectorType* x )
        {
            if ( global_mpi_rank_ == 0 )
            {
                if ( mu_ == 1 )
                    std::cout << "Solving one V-cycle of Multilevel........\n";
                if ( mu_ == 2 )
                    std::cout << "Solving one W-cycle of Multilevel........\n";
            }
            ML ( b, x, pctensor_->GetLevel ( ) );

            if ( global_mpi_rank_ == 0 )
            {
                std::cout << "Multilevel preconditioner completed\n";
            }
            return la::kSolverSuccess;
        }

        template<class LAD>
        la::LinearSolverState PCMultilevelSolver<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            if ( global_mpi_rank_ == 0 )
            {
                if ( mu_ == 1 )
                    std::cout << "Solving V-cycle of Multilevel........\n";
                if ( mu_ == 2 )
                    std::cout << "Solving W-cycle of Multilevel........\n";
            }

            int iter = 0;
            IterateControl::State conv = IterateControl::kIterate;

            VectorType res;
            res.CloneFrom ( *x );
            res.Zeros ( );

            if ( mltype_ == "Free" )
            {
                res_assembler_->assemble_residual ( b, x, res );
            }
            else
            {
                this->op_->VectorMult ( *x, &res );
                res.ScaleAdd ( b, -1.0 );
            }

            DataType ressquared = res.Dot ( res );
            this->res_ = sqrt ( ressquared );
            conv = this->control ( ).Check ( iter, this->res ( ) );

            if ( global_mpi_rank_ == 0 )
            {
                std::cout << "Multilevel starts with residual: " << this->res_ << "\n";
            }

            while ( conv == IterateControl::kIterate )
            {
                iter++;

                ML ( b, x, pctensor_->GetLevel ( ) );

                res.Zeros ( );

                if ( mltype_ == "Free" )
                {
                    res_assembler_->assemble_residual ( b, x, res );
                }
                else
                {
                    this->op_->VectorMult ( *x, &res );
                    res.ScaleAdd ( b, -1.0 );
                }

                ressquared = res.Dot ( res );
                this->res_ = sqrt ( ressquared );

                conv = this->control ( ).Check ( iter, this->res ( ) );
                if ( conv != IterateControl::kIterate )
                {
                    break;
                }
                if ( global_mpi_rank_ == 0 )
                {
                    std::cout << "Multilevel current residual: " << this->res_ << " in iteration: " << iter << "\n";
                }
            }
            if ( global_mpi_rank_ == 0 )
            {
                std::cout << "Multilevel final residual: " << this->res_ << " finished in " << iter << " iterations\n";
            }

            if ( conv == IterateControl::kFailureDivergenceTol ||
                 conv == IterateControl::kFailureMaxitsExceeded )
                return la::kSolverExceeded;
            else
                return la::kSolverSuccess;
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::Solve_Mean_Problem ( VectorType const&b, VectorType* x )
        {
            Solve_Mode_Problem ( 0, b, x );
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::Solve_Mode_Problem ( int mode, VectorType const&b, VectorType* x )
        {
            mean_solver_->Solve ( *b.Mode ( mode ), x->Mode ( mode ) );
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::ML ( VectorType const& b, VectorType* x, int l )
        {
            if ( l == 0 )
            {
                if ( pctensor_->MyRank ( ) == pctensor_->Ownership ( 0, 0 ) )
                    Solve_Mean_Problem ( b, x );
            }
            else
            {
                //1. Apply pre-smoothing nu1 times on fine grid
                Smoothing ( b, x, nu1_ );

                //2. Compute residual on fine grid
                VectorType res_fine;
                res_fine.CloneFrom ( *x );
                res_fine.Zeros ( );

                if ( mltype_ == "Free" )
                {
                    res_assembler_->assemble_residual ( b, x, res_fine );
                }
                else
                {
                    this->op_->VectorMult ( *x, &res_fine );
                    res_fine.ScaleAdd ( b, -1.0 );
                }

                //3. Restrict vector to coarse grid
                res_fine.SynchronizeModes ( pctensor_ );
                x->SynchronizeModes ( pctensor_ );

                pctensor_->SetLevel ( l - 1 );
                VectorType x_coarse;
                VectorType c_coarse;
                VectorType res_coarse;

                VectorType u_fine;
                VectorType u_coarse;

                if ( mltype_ == "Free" )
                {
                    if ( res_assembler_->nonlinear_problem ( ) )
                    {
                        res_assembler_->get_linearization_point ( &u_fine );
                        u_fine.CreateCoarseVector ( pctensor_, u_coarse );
                        res_assembler_->set_linearization_point ( u_coarse );
                    }
                }

                res_fine.CreateCoarseVector ( pctensor_, res_coarse );
                x->CreateCoarseVector ( pctensor_, x_coarse );

                c_coarse.CloneFrom ( res_coarse );
                c_coarse.Zeros ( );

                //4. Correction on coarse grid
                for ( int mu = 0; mu < mu_; ++mu )
                    ML ( res_coarse, &c_coarse, l - 1 );

                //5. Apply correction
                c_coarse.SynchronizeModes ( pctensor_ );
                x_coarse.Axpy ( c_coarse, 1.0 );

                //6. Prolongation
                pctensor_->SetLevel ( l );
                x->AssignToRefinedVector ( pctensor_, x_coarse );

                if ( mltype_ == "Free" )
                {
                    if ( res_assembler_->nonlinear_problem ( ) )
                    {
                        u_fine.AssignToRefinedVector ( pctensor_, u_coarse );
                        res_assembler_->set_linearization_point ( u_fine );
                    }
                }

                //7. Apply post-smoothing nu2 times on fine grid
                Smoothing ( b, x, nu2_ );
            }
        }

        template<class LAD>
        void PCMultilevelSolver<LAD>::Smoothing ( VectorType const& b, VectorType* x, int nu )
        {
            VectorType res_fine, cor;
            res_fine.CloneFrom ( *x );
            res_fine.Zeros ( );
            cor.CloneFrom ( res_fine );

            if ( mltype_ == "Free" )
            {
                res_assembler_->assemble_residual ( b, x, res_fine );
            }
            else
            {
                this->op_->VectorMult ( *x, &res_fine );
                res_fine.ScaleAdd ( b, -1.0 );
            }

            for ( int i = 0; i < nu; ++i )
            {
                cor.Zeros ( );
                for ( int mode = 0; mode < pctensor_->SizeLocal ( ); ++mode )
                {
                    Solve_Mode_Problem ( mode, res_fine, &cor );
                    x->Mode ( mode )->Axpy ( *cor.Mode ( mode ), 1.0 );
                }
                if ( i != nu - 1 )
                {
                    res_fine.Zeros ( );
                    if ( mltype_ == "Free" )
                    {
                        res_assembler_->assemble_residual ( b, x, res_fine );
                    }
                    else
                    {
                        this->op_->VectorMult ( *x, &res_fine );
                        res_fine.ScaleAdd ( b, -1.0 );
                    }
                }
            }
        }

        /// template instantiation
        template class PCMultilevelSolver<la::LADescriptorPolynomialChaosD>;

    } // namespace polynomialchaos
} // namespace hiflow
