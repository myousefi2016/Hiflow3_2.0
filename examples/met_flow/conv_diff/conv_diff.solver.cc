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

#include "conv_diff.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>

/// ***************************************************************************
/// SOLVER
/// ***************************************************************************

void MetFlowConvDiffApp::prepare_linear_solver ( bool ic )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare primal linear solver" << std::endl;
    const int max_iter = this->base_params_["PrimalLinearSolver"]["MaximumIterations"].get<int>( );
    const double abs_tol = this->base_params_["PrimalLinearSolver"]["AbsoluteTolerance"].get<double>( );
    const double rel_tol = this->base_params_["PrimalLinearSolver"]["RelativeTolerance"].get<double>( );
    const double div_tol = this->base_params_["PrimalLinearSolver"]["DivergenceLimit"].get<double>( );
    const int basis_size = this->base_params_["PrimalLinearSolver"]["BasisSize"].get<int>( );

#ifdef USE_HYPRE
    T_precond_.Clear ( );
    T_precond_.SetPreconditioningParameters ( );
    T_precond_.SetNumFunctions ( 1 );
    T_precond_.SetCycleType ( params_["AMG"]["CycleType"].get<int>( ) );
    T_precond_.InitControl ( params_["AMG"]["MaxIterations"].get<int>( ), 0, 0 );
    T_precond_.SetRelaxType ( params_["AMG"]["RelaxType"].get<int>( ) );
    T_precond_.SetRelaxWt ( params_["AMG"]["RelaxWeight"].get<DATATYPE>( ) );
    T_precond_.SetInterpType ( params_["AMG"]["InterpolationType"].get<int>( ) );
    T_precond_.SetStrongThreshold ( params_["AMG"]["StrongThreshold"].get<DATATYPE>( ) );
    T_precond_.SetAggNumLevels ( params_["AMG"]["AggNumLevels"].get<int>( ) );
    T_precond_.SetCoarsenType ( params_["AMG"]["CoarsenType"].get<int>( ) );
    T_precond_.SetCycleNumSweeps ( params_["AMG"]["NumDownSweeps"].get<int>( ), 1 );
    T_precond_.SetCycleNumSweeps ( params_["AMG"]["NumUpSweeps"].get<int>( ), 2 );
    T_precond_.SetSmoothType ( params_["AMG"]["SmoothType"].get<int>( ) );
    T_precond_.SetSmoothNumLevels ( params_["AMG"]["SmoothNumLevels"].get<int>( ) );
    T_precond_.SetMaxCoarseSize ( params_["AMG"]["MaxCoarseSize"].get<int>( ) );
    T_precond_.SetMaxLevels ( params_["AMG"]["MaxLevels"].get<int>( ) );
    //T_precond_.SetNodal(1);
    //T_precond_.SetNodalDiag(1);
    T_precond_.Init ( );
    T_linear_solver_->InitControl ( max_iter, 1e-12, 1e-6, 1e6 );

    T_linear_solver_->SetupPreconditioner ( T_precond_ );
    T_linear_solver_->InitParameter ( basis_size, "RightPreconditioning" );
    T_linear_solver_->SetupOperator ( T_matrix_dual_ );
#else
    T_linear_solver_->InitControl ( max_iter, abs_tol, rel_tol, div_tol );

#    ifdef WITH_ILUPP
    T_linear_solver_->InitParameter ( basis_size, "RightPreconditioning" );
    // prepare ILUPP
    const int prepro_type = this->base_params_["ILUPP"]["PreprocessingType"].get<int>( );
    const int precond_no = this->base_params_["ILUPP"]["PreconditionerNumber"].get<int>( );
    const int max_levels = this->base_params_["ILUPP"]["MaxMultilevels"].get<int>( );
    const double mem_factor = this->base_params_["ILUPP"]["MemFactor"].get<double>( );
    const double threshold = this->base_params_["ILUPP"]["PivotThreshold"].get<double>( );
    const double min_pivot = this->base_params_["ILUPP"]["MinPivot"].get<double>( );

    T_precond_.InitParameter ( prepro_type, precond_no, max_levels, mem_factor, threshold, min_pivot );
    T_linear_solver_->SetupPreconditioner ( this->T_precond_ );
#    else
    T_linear_solver_->InitParameter ( basis_size, "NoPreconditioning" );
#    endif
    // set the matrix to be used as the operator
    T_linear_solver_->SetupOperator ( *matrix_ );

#endif
}

void MetFlowConvDiffApp::prepare_linear_solver_dual ( )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Prepare dual linear solver" << std::endl;

    const int max_iter = this->base_params_["DualLinearSolver"]["MaximumIterations"].get<int>( );
    const double abs_tol = this->base_params_["DualLinearSolver"]["AbsoluteTolerance"].get<double>( );
    const double rel_tol = this->base_params_["DualLinearSolver"]["RelativeTolerance"].get<double>( );
    const double div_tol = this->base_params_["DualLinearSolver"]["DivergenceLimit"].get<double>( );
    const int basis_size = this->base_params_["DualLinearSolver"]["BasisSize"].get<int>( );

#ifdef USE_HYPRE
    T_precond_dual_.Clear ( );
    T_precond_dual_.SetPreconditioningParameters ( );
    T_precond_dual_.SetNumFunctions ( 1 );
    T_precond_dual_.SetCycleType ( params_["AMG"]["CycleType"].get<int>( ) );
    T_precond_dual_.InitControl ( params_["AMG"]["MaxIterations"].get<int>( ), 0, 0 );
    T_precond_dual_.SetRelaxType ( params_["AMG"]["RelaxType"].get<int>( ) );
    T_precond_dual_.SetRelaxWt ( params_["AMG"]["RelaxWeight"].get<DATATYPE>( ) );
    T_precond_dual_.SetInterpType ( params_["AMG"]["InterpolationType"].get<int>( ) );
    T_precond_dual_.SetStrongThreshold ( params_["AMG"]["StrongThreshold"].get<DATATYPE>( ) );
    T_precond_dual_.SetAggNumLevels ( params_["AMG"]["AggNumLevels"].get<int>( ) );
    T_precond_dual_.SetCoarsenType ( params_["AMG"]["CoarsenType"].get<int>( ) );
    T_precond_dual_.SetCycleNumSweeps ( params_["AMG"]["NumDownSweeps"].get<int>( ), 1 );
    T_precond_dual_.SetCycleNumSweeps ( params_["AMG"]["NumUpSweeps"].get<int>( ), 2 );
    T_precond_dual_.SetSmoothType ( params_["AMG"]["SmoothType"].get<int>( ) );
    T_precond_dual_.SetSmoothNumLevels ( params_["AMG"]["SmoothNumLevels"].get<int>( ) );
    T_precond_dual_.SetMaxCoarseSize ( params_["AMG"]["MaxCoarseSize"].get<int>( ) );
    T_precond_dual_.SetMaxLevels ( params_["AMG"]["MaxLevels"].get<int>( ) );
    //T_precond_dual_.SetNodal(1);
    //T_precond_dual_.SetNodalDiag(1);
    T_precond_dual_.Init ( );
    T_linear_solver_dual_->InitControl ( max_iter, 1e-12, 1e-6, 1e6 );

    T_linear_solver_dual_->SetupPreconditioner ( T_precond_ );
    T_linear_solver_dual_->InitParameter ( basis_size, "RightPreconditioning" );
    T_linear_solver_dual_->SetupOperator ( T_matrix_dual_ );
#else

    T_linear_solver_dual_->InitControl ( max_iter, abs_tol, rel_tol, 1e6 );

#    ifdef WITH_ILUPP
    this->T_linear_solver_dual_->InitParameter ( basis_size, "RightPreconditioning" );

    // prepare ILUPP
    const int prepro_type = this->base_params_["ILUPP"]["PreprocessingType"].get<int>( );
    const int precond_no = this->base_params_["ILUPP"]["PreconditionerNumber"].get<int>( );
    const int max_levels = this->base_params_["ILUPP"]["MaxMultilevels"].get<int>( );
    const double mem_factor = this->base_params_["ILUPP"]["MemFactor"].get<double>( );
    const double threshold = this->base_params_["ILUPP"]["PivotThreshold"].get<double>( );
    const double min_pivot = this->base_params_["ILUPP"]["MinPivot"].get<double>( );

    this->T_precond_dual_.InitParameter ( prepro_type, precond_no, max_levels, mem_factor, threshold, min_pivot );
    this->T_linear_solver_dual_->SetupPreconditioner ( this->T_precond_dual_ );
#    else
    this->T_linear_solver_dual_->InitParameter ( basis_size, "NoPreconditioning" );
#    endif
    // set the matrix to be used as the operator
    this->T_linear_solver_dual_->SetupOperator ( *this->matrix_dual_ );
#endif
}

void MetFlowConvDiffApp::update_preconditioner_dual ( LAD::MatrixType* DF )
{
#ifndef USE_HYPRE
#    ifdef WITH_ILUPP
    this->T_precond_dual_.SetupOperator ( *DF );
#    endif
#endif
}

void MetFlowConvDiffApp::update_preconditioner ( const LAD::VectorType& u, LAD::MatrixType* DF )
{
#ifndef USE_HYPRE
#    ifdef WITH_ILUPP
    this->T_precond_.SetupOperator ( *DF );
#    endif
#endif
}
