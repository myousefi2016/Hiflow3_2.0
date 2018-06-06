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

#ifndef MY_NEWTON_H
#    define MY_NEWTON_H

///
/// \file my_newton.h
/// \brief Adjusted Newton solver for use in FE space on an h-adapted mesh. Hanging nodes on the adapted mesh lead to constraints in global solution vector.
///
/// \author Teresa Beck, Philipp Gerstner
///
#    include <sstream>
#    include <string>
#    include <fstream>
#    include "hiflow.h"
#    include "nonlinear/newton.h"
#    include "common/timer.h"
#    include "../tmp_config/met_flow_vars.h"

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;
using namespace hiflow::la;

template<class LAD>
class MyNewton : public Newton<LAD>
{
  public:
    typedef typename LAD::MatrixType MatrixType;
    typedef typename LAD::VectorType VectorType;
    typedef typename LAD::DataType DataType;

    MyNewton ( VectorType* residual, MatrixType* matrix ) : res_ ( residual ), jac_ ( matrix )
    {
        Newton<LAD>( residual, matrix );
        Newton<LAD>::InitParameter ( residual, matrix );
    };

    virtual NonlinearSolverState Solve ( const VectorType& y, VectorType* x )
    {
        std::cout << "DONT USE THIS FUNCTION! \nFUNCTION IS DEPRECATED IN ADAPTIVE MODE! \n\n";
        interminable_assert ( 0 );
    };

    inline NonlinearSolverState Solve ( const VectorType& rhs, VectorType* x, const VectorSpace<double>& space );

    void SetOutputFilename ( std::string filename )
    {
        this->filename_linear_solver_ = filename;
    }

  private:
    VectorType* res_;
    MatrixType* jac_;
    std::string filename_linear_solver_;
};

template<class LAD>
NonlinearSolverState MyNewton<LAD>::Solve ( const VectorType& rhs, VectorType* x, const VectorSpace<double>& space )
{
    assert ( this->op_ != NULL );
    assert ( this->linsolve_ != NULL );
    assert ( x != NULL );

    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    //Init
    LinearSolverState LinSolState = kSolverSuccess;
    IterateControl::State conv = IterateControl::kIterate;
    this->res_->Clear ( );
    this->res_->CloneFromWithoutContent ( rhs );

    VectorType* cor = new VectorType ( );
    cor->Clear ( );
    cor->CloneFromWithoutContent ( rhs );

    VectorType* sol = new VectorType ( );
    sol->Clear ( );
    if ( this->InitialSolution_ == this->NewtonInitialSolutionOwn )
    {
        sol->CloneFrom ( *x );
    }
    else if ( this->InitialSolution_ == this->NewtonInitialSolution0 )
    {
        sol->CloneFromWithoutContent ( rhs );
        sol->Zeros ( );
    }
    else
    {
        return kNonlinearSolverInitError;
    }

    //Step 0
    double resnorm;
    double resnorm_old;

    this->iter_ = 0;
    this->op_->Reinit ( );
    this->ComputeResidual ( *sol, rhs, this->res_ );
    resnorm = this->GetResidual ( );
    resnorm_old = resnorm;

    if ( rank == 0 ) std::cout << "  [" << this->iter_ << "] Newton residual norm: " << resnorm << std::endl;

    conv = this->control ( ).Check ( this->iter ( ), this->GetResidual ( ) );
    if ( this->ForcingStrategy_ == this->NewtonForcingStrategyOwn )
    {
        assert ( this->ForcingStratObject_ != NULL );
        this->ForcingStratObject_->SetResidualNorm ( this->res_->Norm2 ( ) );
    }

    while ( conv == IterateControl::kIterate )
    {
        //NextStep
        this->iter_++;
        this->ComputeJacobian ( *sol, this->jac_ );

        Timer timer;
        timer.start ( );
        LinSolState = this->SolveJacobian ( *this->jac_, *this->res_, cor );
        timer.stop ( );

        if ( LinSolState == kSolverError ) break;

        this->UpdateSolution ( *cor, rhs, sol );
        interpolate_constrained_vector ( space, *sol );
        sol->Update ( );

        resnorm_old = resnorm;
        resnorm = this->res_->Norm2 ( );
        double duration = timer.get_duration ( );
        double lin_iter = this->linsolve_->iter ( );
        double lin_res = this->linsolve_->res ( );
        double newton_res = this->GetResidual ( );
        double newton_iter = this->iter_;
        double cor_norm = cor->Norm2 ( );

        if ( rank == 0 )
        {
            std::cout << "  [" << newton_iter << "] Newton abs. ResNorm: " << newton_res << ", iter: " << lin_iter
                    << ", CPU: " << duration << ", abs. ResNorm: " << lin_res << ", rel. ResNorm: " << lin_res / resnorm_old
                    << ", cor norm " << cor_norm << std::endl;

            std::string path = this->filename_linear_solver_;
            std::ofstream out;
            out.open ( path.c_str ( ), std::ios::out | std::ios::app );
            out.precision ( 6 );
            out << std::scientific;
            out << newton_iter << " " << newton_res << " " << lin_iter << " " << duration << " " << lin_res << " " << lin_res / resnorm_old << "\n";
            out.close ( );
        }

        conv = this->control ( ).Check ( this->iter ( ), this->GetResidual ( ) );
        if ( this->ForcingStrategy_ == this->NewtonForcingStrategyOwn )
        {
            assert ( this->ForcingStratObject_ != NULL );
            this->ForcingStratObject_->ComputeForcingTerm ( this->res_->Norm2 ( ), this->linsolve_->res ( ) );
        }
    }
    delete cor;

    if ( LinSolState == kSolverError )
    {
        sol->Clear ( );
        delete sol;
        return kNonlinearSolverError;
    }
    else if ( conv == IterateControl::kFailureDivergenceTol && conv == IterateControl::kFailureMaxitsExceeded )
    {
        sol->Clear ( );
        delete sol;
        return kNonlinearSolverExceeded;
    }
    else
    {
        x->CloneFrom ( *sol );
        sol->Clear ( );
        delete sol;
        return kNonlinearSolverSuccess;
    }
}

#endif //MY_NEWTON_H
