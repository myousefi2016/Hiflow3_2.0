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

/// \author Tobias Hahn

#ifndef HIFLOW_NONLINEAR_NONLINEAR_PROBLEM_FROM_APP_H_
#    define HIFLOW_NONLINEAR_NONLINEAR_PROBLEM_FROM_APP_H_

#    include "nonlinear_problem.h"
#    include "application.h"

namespace hiflow
{

    /// Simplification, waiting for final integration routine

    template<class LAD>
    class NonlinearProblemFromApp : public NonlinearProblem<LAD>
    {
      public:
        typedef typename LAD::MatrixType DataType;
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;

        NonlinearProblemFromApp ( );

        NonlinearProblemFromApp ( Application<DataType> &newapp )
        {
            app_ = &newapp;
        }
        ~NonlinearProblemFromApp ( app_ = NULL; );

        /// user supplied init function

        void Reinit ( Application<DataType> &newapp )
        {
            app_ = &newapp;
        }

        /// evaluate function

        void EvalFunc ( const VectorType& in, VectorType* out )
        {
            //Integrate<Laplace<LAD::DataType>, LAD, UnitIntegrator2d>(solvec, quadrature, laplace, &Laplace<LAD::DataType>::rhs2d_1, &b);
        }

        /// compute gradient

        void EvalGrad ( const VectorType& in, MatrixType* out )
        {
            //assemble jacobian , allocate if needed

        }

      private:
        Application<DataType> *app_;
    };

} // namespace hiflow

#endif
