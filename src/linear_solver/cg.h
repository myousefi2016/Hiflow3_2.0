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

/// @author Chandramowli Subramanian

#ifndef HIFLOW_LINEARSOLVER_CG_H_
#    define HIFLOW_LINEARSOLVER_CG_H_

#    include <string>
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/linear_solver_creator.h"

namespace hiflow
{
    namespace la
    {

        /// @brief CG solver
        ///
        /// CG solver for symmetric linear systems Ax=b.

        template<class LAD>
        class CG : public LinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            CG ( );
            virtual ~CG ( );

            void InitParameter ( std::string method );

            LinearSolverState Solve ( const VectorType& b, VectorType* x );

          private:

            LinearSolverState SolveNoPrecond ( const VectorType& b, VectorType* x );

            LinearSolverState SolvePrecond ( const VectorType& b, VectorType* x );

        };

        /// @brief CG creator class
        /// @author Tobias Hahn

        template<class LAD>
        class CGcreator : public LinearSolverCreator<LAD>
        {
          public:

            LinearSolver<LAD>* params ( const PropertyTree& c )
            {
                CG<LAD>* newCG = new CG<LAD>( );
                if ( c.contains ( "Method" ) )
                    newCG->InitParameter ( c["Method"].template get<std::string>( ).c_str ( ) );
                if ( c.contains ( "MaxIterations" ) && c.contains ( "AbsTolerance" ) &&
                     c.contains ( "RelTolerance" ) && c.contains ( "DivTolerance" ) )
                    newCG->InitControl ( c["MaxIterations"].template get<int>( ),
                                         c["AbsTolerance"].template get<double>( ),
                                         c["RelTolerance"].template get<double>( ),
                                         c["DivTolerance"].template get<double>( ) );
                return newCG;
            }
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_CG_H_
