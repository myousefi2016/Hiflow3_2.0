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

#ifndef HIFLOW_NONLINEAR_NONLINEAR_PROBLEM_H_
#    define HIFLOW_NONLINEAR_NONLINEAR_PROBLEM_H_

namespace hiflow
{

    /// @brief Newton nonlinear solver
    /// @author Tobias Hahn
    ///
    /// Base class for nonlinear problems as required by the Nonlinear
    /// Solver class. Basically provides a function evaluation routine and
    /// if applicable also for the gradient

    template<class LAD>
    class NonlinearProblem
    {
      public:
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;

        NonlinearProblem ( )
        {
        }

        virtual ~NonlinearProblem ( )
        {
        }

        /// Optional user supplied init function

        virtual void Reinit ( )
        {
        }

        /// Evaluates function at given point

        virtual void EvalFunc ( const VectorType& in, VectorType* out )
        {
        }

        /// Evaluates function at given point

        virtual void EvalFuncNonConst ( VectorType& in, VectorType* out )
        {
        }

        /// Compute gradient at given point

        virtual void EvalGrad ( const VectorType& in, MatrixType* out )
        {
        }

        /// Compute gradient at given point

        virtual void EvalGradNonConst ( VectorType& in, MatrixType* out )
        {
        }

        /// Compute pressure filter or something other custom after each Newton Step

        virtual void ApplyFilter ( VectorType& u )
        {
        }

    };

} // namespace hiflow

#endif
