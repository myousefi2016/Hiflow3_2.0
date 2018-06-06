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

#include "common/macros.h"
#include "common/log.h"
#include "newton.h"

#include <cassert>
#include <vector>

#include "linear_algebra/la_descriptor.h"

#include "damping_strategy.h"
#include "forcing_strategy.h"
#include "nonlinear_problem.h"

/// @author Tobias Hahn, Michael Schick

using namespace hiflow::la;

namespace hiflow
{

    /// Sets up paramaters for initial solution, forcing and damping

    template<class LAD>
    NonlinearSolverState Newton<LAD>::InitParameter ( NonlinearSolverParameter param )
    {
        if ( param == NewtonDampingStrategyNone || param == NewtonDampingStrategyOwn )
        {
            this->DampingStrategy_ = param;
        }
        else if ( param == NewtonForcingStrategyConstant || param == NewtonForcingStrategyOwn )
        {
            this->ForcingStrategy_ = param;
        }
        else if ( param == NewtonInitialSolution0 || param == NewtonInitialSolutionOwn )
        {
            this->InitialSolution_ = param;
        }
        else
        {
            return kNonlinearSolverInitError;
        }

        return kNonlinearSolverSuccess;
    }

    template<class LAD>
    void Newton<LAD>::InitParameter ( VectorType* residual, MatrixType* matrix )
    {
        res_ = residual;
        jac_ = matrix;
        assert ( res_ != NULL );
        assert ( jac_ != NULL );
    }

    /// Sets up the damping strategy
    /// @param dampstrat DampingStrategy

    template<class LAD>
    void Newton<LAD>::SetDampingStrategy ( DampingStrategy<LAD>& dampstrat )
    {
        this->DampStratObject_ = &dampstrat;
        this->DampingStrategy_ = NewtonDampingStrategyOwn;
    }

    /// Sets up the forcing strategy
    /// @param forcingstrat ForcingStrategy

    template<class LAD>
    void Newton<LAD>::SetForcingStrategy ( ForcingStrategy<LAD>& forcingstrat )
    {
        this->ForcingStratObject_ = &forcingstrat;
        this->ForcingStrategy_ = NewtonForcingStrategyOwn;
    }

    /// Get the current forcing term
    /// returns zero if no forcing is activated

    template<class LAD>
    void Newton<LAD>::GetForcingTerm ( DataType& forcing ) const
    {
        if ( this->ForcingStrategy_ == NewtonForcingStrategyOwn )
        {
            forcing = this->ForcingStratObject_->GetCurrentForcingTerm ( );
        }
        else
        {
            forcing = 0.;
        }
    }

    /// Set a forcing term, necessary in combination with damping
    /// only if forcing is activated
    /// @param forcing DataType

    template<class LAD>
    void Newton<LAD>::SetForcingTerm ( DataType forcing )
    {
        if ( this->ForcingStrategy_ == NewtonForcingStrategyOwn )
        {
            this->ForcingStratObject_->SetForcingTerm ( forcing );
        }
    }

    /// Provides information is forcing is used
    /// necessary for damping

    template<class LAD>
    bool Newton<LAD>::Forcing ( ) const
    {
        if ( this->ForcingStrategy_ == NewtonForcingStrategyOwn )
        {
            return true;
        }
        else
            return false;
    }

    /// Updates solution vector using correction and possible damping/
    /// forcing strategies that use in turn a right-hand side
    /// @param cor correction vector
    /// @param rhs right-hand side vector
    /// @param sol solution vector

    template<class LAD>
    NonlinearSolverState Newton<LAD>::UpdateSolution ( const VectorType& cor, const VectorType& rhs, VectorType* sol )
    {

        if ( this->DampingStrategy_ == NewtonDampingStrategyNone )
        {
            assert ( sol->size_local ( ) == cor.size_local ( ) );
            assert ( sol->size_global ( ) == cor.size_global ( ) );
            sol->Axpy ( cor, static_cast < DataType > ( -1. ) );
            sol->Update ( );
            this->op_->ApplyFilter ( *sol );
            sol->Update ( );
            if ( non_const_mode_ )
                this->ComputeResidualNonConst ( *sol, rhs, this->res_ );
            else
                this->ComputeResidual ( *sol, rhs, this->res_ );
            LOG_INFO ( "Newton residual norm", this->GetResidual ( ) );
            return kNonlinearSolverSuccess;
        }
        else if ( this->DampingStrategy_ == NewtonDampingStrategyOwn )
        {
            assert ( DampStratObject_ != NULL );
            DampingState state = DampStratObject_->Update ( cor, rhs, this->res_, sol, this );
            this->residual_ = DampStratObject_->GetResidual ( );
            this->resids_.push_back ( this->residual_ );
            LOG_INFO ( "Newton residual norm", this->GetResidual ( ) );
            if ( state != 0 ) return kNonlinearSolverError;
        }
        else
            return kNonlinearSolverError;

        return kNonlinearSolverSuccess;
    }

    /// Updates solution vector using correction and possible damping/
    /// forcing strategies
    /// @param cor correction vector
    /// @param sol solution vector

    template<class LAD>
    NonlinearSolverState Newton<LAD>::UpdateSolution ( const VectorType& cor, VectorType* sol )
    {
        VectorType* rhs = new VectorType ( );
        rhs->Clear ( );
        NonlinearSolverState state = this->UpdateSolution ( cor, *rhs, sol );
        delete rhs;
        return state;
    }

    /// Returns jacobian matrix J of nonlinear problem F at x
    /// @param x point of evaluation
    /// @param jacobian at x

    template<class LAD>
    void Newton<LAD>::ComputeJacobian ( const VectorType& x, MatrixType* jacobian )
    {
        this->op_->EvalGrad ( x, jacobian );
    }

    /// Returns jacobian matrix J of nonlinear problem F at x
    /// @param x point of evaluation
    /// @param jacobian at x

    template<class LAD>
    void Newton<LAD>::ComputeJacobianNonConst ( VectorType& x, MatrixType* jacobian )
    {
        this->op_->EvalGradNonConst ( x, jacobian );
    }

    /// Solves linear problem J*c=r
    /// If Forcing is activated, then the system
    /// is solved in an inexact way with
    /// relative tolerance determined by forcing terms
    /// @param jacobian jacobian matrix J
    /// @param residual residual vector r
    /// @param correction correction vector c

    template<class LAD>
    LinearSolverState Newton<LAD>::SolveJacobian ( const MatrixType& jacobian, const VectorType& residual, VectorType* correction )
    {
        assert ( correction != NULL );

        // Reset correction if needed
        if ( ( residual.size_local ( ) != correction->size_local ( ) ) || ( residual.size_global ( ) != correction->size_global ( ) ) )
        {
            correction->CloneFromWithoutContent ( residual );
        }

        // start vector
        correction->Zeros ( );

        correction->Update ( );

        if ( this->ForcingStrategy_ == NewtonForcingStrategyOwn )
        {
            this->linsolve_->SetRelativeTolerance ( this->ForcingStratObject_->GetCurrentForcingTerm ( ) );
            LOG_INFO ( "Forcing term", this->ForcingStratObject_->GetCurrentForcingTerm ( ) );
        }

        // solve
        LinearSolverState state = this->linsolve_->Solve ( residual, correction );

        correction->Update ( );

        return state; //solve jacobian
    }

    /// Computes residual vector F(sol)-rhs for non-linear problem F with
    /// right-hand side rhs
    /// @param sol solution vector
    /// @param rhs right-hand side vector
    /// @param res residual vector

    template<class LAD>
    void Newton<LAD>::ComputeResidual ( const VectorType& sol, const VectorType& rhs,
                                        VectorType* res )
    {
        assert ( res != NULL );
        // Reset residual if needed
        if ( ( res->size_local ( ) != sol.size_local ( ) ) ||
             ( res->size_global ( ) != sol.size_global ( ) ) )
        {
            res->CloneFromWithoutContent ( sol );
            res->Zeros ( );
        }

        // Compute residual
        this->op_->EvalFunc ( sol, res );
        if ( ( rhs.size_local ( ) == res->size_local ( ) ) &&
             ( rhs.size_global ( ) == res->size_global ( ) ) )
        {
            res->Axpy ( rhs, static_cast < DataType > ( -1. ) );
        }

        //res->Update( );

        // Compute new residual norm
        this->residual_ = res->Norm2 ( );

        // Store it in vector
        this->resids_.push_back ( this->residual_ );
    }

    /// Computes residual vector F(sol)-rhs for non-linear problem F with
    /// right-hand side rhs
    /// @param sol solution vector
    /// @param rhs right-hand side vector
    /// @param res residual vector

    template<class LAD>
    void Newton<LAD>::ComputeResidualNonConst ( VectorType& sol, const VectorType& rhs,
                                                VectorType* res )
    {
        assert ( res != NULL );
        // Reset residual if needed
        if ( ( res->size_local ( ) != sol.size_local ( ) ) ||
             ( res->size_global ( ) != sol.size_global ( ) ) )
        {
            res->CloneFromWithoutContent ( sol );
            res->Zeros ( );
        }

        // Compute residual
        this->op_->EvalFuncNonConst ( sol, res );
        if ( ( rhs.size_local ( ) == res->size_local ( ) ) &&
             ( rhs.size_global ( ) == res->size_global ( ) ) )
        {
            res->Axpy ( rhs, static_cast < DataType > ( -1. ) );
        }

        //res->Update( );

        // Compute new residual norm
        this->residual_ = res->Norm2 ( );

        // Store it in vector
        this->resids_.push_back ( this->residual_ );
    }

    /// Computes residual vector F(sol) for non-linear problem F
    /// @param sol solution vector
    /// @param res residual vector

    template<class LAD>
    void Newton<LAD>::ComputeResidual ( const VectorType& sol, VectorType* res )
    {
        VectorType* rhs = new VectorType ( );
        rhs->Clear ( );
        this->ComputeResidual ( sol, *rhs, res );
        delete rhs;
    }

    /// Solves F(x)=y
    /// @param y right hand side vectorNewtonDampingStrategyArmijo
    /// @param x solution vector
    /// @return status if solver succeeded

    template<class LAD>
    NonlinearSolverState Newton<LAD>::Solve ( const VectorType& rhs, VectorType* x )
    {
        assert ( this->res_ != NULL );
        assert ( this->jac_ != NULL );
        assert ( this->op_ != NULL );
        assert ( this->linsolve_ != NULL );
        assert ( x != NULL );

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
        if ( InitialSolution_ == NewtonInitialSolutionOwn )
            sol->CloneFrom ( *x );
        else if ( InitialSolution_ == NewtonInitialSolution0 )
        {
            sol->CloneFromWithoutContent ( rhs );
            sol->Zeros ( );
        }
        else
            return kNonlinearSolverInitError;

        //Step 0
        this->iter_ = 0;
        this->op_->Reinit ( );

        sol->Update ( );

        if ( non_const_mode_ )
            this->ComputeResidualNonConst ( *sol, rhs, this->res_ );
        else
            this->ComputeResidual ( *sol, rhs, this->res_ );

        conv = this->control ( ).Check ( this->iter ( ), this->GetResidual ( ) );
        if ( this->ForcingStrategy_ == NewtonForcingStrategyOwn )
        {
            assert ( this->ForcingStratObject_ != NULL );
            this->ForcingStratObject_->SetResidualNorm ( this->res_->Norm2 ( ) );
        }

        LOG_INFO ( "Newton starts with residual norm", this->GetResidual ( ) );
        if ( this->ForcingStrategy_ == NewtonForcingStrategyOwn )
        {
            LOG_INFO ( "Newton starts with forcing term", this->ForcingStratObject_->GetCurrentForcingTerm ( ) );
        }

        while ( conv == IterateControl::kIterate )
        {
            //NextStep
            this->iter_++;
            // sol->Update();
            if ( non_const_mode_ )
                this->ComputeJacobianNonConst ( *sol, this->jac_ );
            else
                this->ComputeJacobian ( *sol, this->jac_ );
            LinSolState = this->SolveJacobian ( *this->jac_, *this->res_, cor );
            if ( LinSolState == kSolverError ) break;
            this->UpdateSolution ( *cor, rhs, sol );
            //sol->Update( );

            conv = this->control ( ).Check ( this->iter ( ), this->GetResidual ( ) );
            if ( this->ForcingStrategy_ == NewtonForcingStrategyOwn )
            {
                assert ( this->ForcingStratObject_ != NULL );
                this->ForcingStratObject_->ComputeForcingTerm ( this->res_->Norm2 ( ), this->linsolve_->res ( ) );
                LOG_INFO ( "New forcing term", this->ForcingStratObject_->GetCurrentForcingTerm ( ) );
            }
        }
        delete cor;

        LOG_INFO ( "Newton finished with residual norm", this->GetResidual ( ) );

        if ( LinSolState == kSolverError )
        {
            delete sol;
            return kNonlinearSolverError;
        }
        else if ( conv == IterateControl::kFailureDivergenceTol ||
                  conv == IterateControl::kFailureMaxitsExceeded )
        {
            delete sol;
            return kNonlinearSolverExceeded;
        }
        else
        {
            x->CopyFrom ( *sol );
            x->Update ( );
            delete sol;
            return kNonlinearSolverSuccess;
        }
    }

    /// Solves F(x)=0
    /// @param x solution vector
    /// @return status if solver succeeded

    template<class LAD>
    NonlinearSolverState Newton<LAD>::Solve ( VectorType* x )
    {
        VectorType* rhs = new VectorType ( );
        rhs->Clear ( );
        rhs->CloneFromWithoutContent ( *x );
        rhs->Zeros ( );
        NonlinearSolverState state = this->Solve ( *rhs, x );
        delete rhs;
        return state;
    }

    /// Plain constructor that does not initialize the solver.

    template<class LAD>
    Newton<LAD>::Newton ( )
    {
        this->InitialSolution_ = NewtonInitialSolution0;
        this->DampingStrategy_ = NewtonDampingStrategyNone;
        this->ForcingStrategy_ = NewtonForcingStrategyConstant;
        this->non_const_mode_ = false;
    }

    /// Standard constructor requiring pointers to user reserved space
    /// for residual vector and jacobian. Sets no forcing and no damping and
    /// uses initial solution zero.

    template<class LAD>
    Newton<LAD>::Newton ( VectorType* residual, MatrixType* matrix )
    : res_ ( residual ), jac_ ( matrix )
    {
        assert ( res_ != NULL );
        assert ( jac_ != NULL );

        this->InitialSolution_ = NewtonInitialSolution0;
        this->DampingStrategy_ = NewtonDampingStrategyNone;
        this->ForcingStrategy_ = NewtonForcingStrategyConstant;
        this->non_const_mode_ = false;
    }

    /// Standard destructor

    template<class LAD>
    Newton<LAD>::~Newton ( )
    {
        this->res_ = NULL;
        this->jac_ = NULL;
        this->linsolve_ = NULL;
        this->op_ = NULL;
        this->DampStratObject_ = NULL;
        this->ForcingStratObject_ = NULL;
    }

    /// template instantiation
    template class Newton<LADescriptorCoupledD>;
    template class Newton<LADescriptorCoupledS>;
    template class Newton<LADescriptorHypreD>;
    template class Newton<LADescriptorPolynomialChaosD>;
    template class Newton<LADescriptorPolynomialChaosExpansionD>;

} // namespace hiflow
