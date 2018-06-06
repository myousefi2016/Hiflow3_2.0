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

#ifndef HIFLOW_NONLINEAR_NEWTON_H_
#    define HIFLOW_NONLINEAR_NEWTON_H_

#    include "common/property_tree.h"
#    include "nonlinear/nonlinear_solver.h"
#    include "nonlinear/nonlinear_solver_creator.h"

namespace hiflow
{

    template<class LAD> class DampingStrategy;
    template<class LAD> class ForcingStrategy;

    /// @brief Newton nonlinear solver
    /// @author Tobias Hahn, Michael Schick
    ///
    /// Newton solver for generic nonlinear problems.
    /// Requires allocated space for jacobian and residual.
    /// May use forcing, damping and custom initial solution, but does not
    /// by default.

    template<class LAD>
    class Newton : public NonlinearSolver<LAD>
    {
      public:
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;
        typedef typename LAD::DataType DataType;

        enum NonlinearSolverParameter
        {
            NewtonInitialSolutionOwn,
            NewtonInitialSolution0,
            NewtonDampingStrategyNone,
            NewtonDampingStrategyOwn,
            NewtonForcingStrategyConstant,
            NewtonForcingStrategyOwn
        };

        Newton ( );
        Newton ( VectorType* residual, MatrixType* matrix );
        ~Newton ( );

        void InitParameter ( VectorType* residual, MatrixType* matrix );
        NonlinearSolverState InitParameter ( NonlinearSolverParameter param );
        void SetDampingStrategy ( DampingStrategy<LAD>& dampstrat );
        void SetForcingStrategy ( ForcingStrategy<LAD>& forcingstrat );
        bool Forcing ( ) const;

        void GetForcingTerm ( DataType& forcing ) const;
        void SetForcingTerm ( DataType forcing );

        void ActivateNonConstMode ( )
        {
            non_const_mode_ = true;
        }

        bool GetNonConstMode ( )
        {
            return non_const_mode_;
        }

        NonlinearSolverState Solve ( const VectorType& y, VectorType* x );
        NonlinearSolverState Solve ( VectorType* x );

        void ComputeJacobian ( const VectorType& x, MatrixType* jacobian );
        void ComputeJacobianNonConst ( VectorType& x, MatrixType* jacobian );
        la::LinearSolverState SolveJacobian ( const MatrixType& jacobian, const VectorType& residual, VectorType* correction );

        NonlinearSolverState UpdateSolution ( const VectorType& cor, VectorType* sol );
        NonlinearSolverState UpdateSolution ( const VectorType& cor, const VectorType& rhs, VectorType* sol );

        void ComputeResidual ( const VectorType& sol, VectorType* res );
        void ComputeResidual ( const VectorType& sol, const VectorType& rhs, VectorType* res );
        void ComputeResidualNonConst ( VectorType& sol, const VectorType& rhs, VectorType* res );

      protected:
        VectorType* res_;
        MatrixType* jac_;
        DampingStrategy<LAD>* DampStratObject_;
        ForcingStrategy<LAD>* ForcingStratObject_;
        std::vector<DataType> resids_;

        NonlinearSolverParameter InitialSolution_;
        NonlinearSolverParameter DampingStrategy_;
        NonlinearSolverParameter ForcingStrategy_;

        bool non_const_mode_;
    };

    /// @brief Newton creator class
    /// @author Tobias Hahn

    template<class LAD>
    class Newtoncreator : public NonlinearSolverCreator<LAD>
    {
      public:
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;

        NonlinearSolver<LAD>* params ( VectorType* residual, MatrixType* matrix, const PropertyTree& c )
        {
            Newton<LAD>* newNewton = new Newton<LAD>( residual, matrix );
            if ( c.contains ( "MaxIterations" ) && c.contains ( "AbsTolerance" ) &&
                 c.contains ( "RelTolerance" ) && c.contains ( "DivTolerance" ) )
                newNewton->InitControl ( c["MaxIterations"].template get<int>( ),
                                         c["AbsTolerance"].template get<double>( ),
                                         c["RelTolerance"].template get<double>( ),
                                         c["DivTolerance"].template get<double>( ) );
            return newNewton;
        }

        NonlinearSolver<LAD>* params ( const PropertyTree& c )
        {
            return new Newton<LAD>( );
        }
    };

} // namespace hiflow

#endif  // HIFLOW_NONLINEAR_NEWTON_H_
