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

/// @author Simon Gawlok

#ifndef HIFLOW_LINEARSOLVER_BICGSTAB_H_
#    define HIFLOW_LINEARSOLVER_BICGSTAB_H_

#    include <string>
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/linear_solver_creator.h"

namespace hiflow
{
    namespace la
    {

        /// @brief BiCGSTAB solver
        ///
        /// BiCGSTAB solver for regular linear systems Ax=b.

        template<class LAD>
        class BiCGSTAB : public LinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            BiCGSTAB ( );
            virtual ~BiCGSTAB ( );

            void InitParameter ( std::string method );
            LinearSolverState Solve ( const VectorType& b, VectorType* x );

            /// Sets the relative tolerance.
            /// Needed by Inexact Newton Methods
            /// @param reltol relative tolerance of residual to converge

            void SetRelativeTolerance ( double reltol )
            {
                int maxits = this->control_.maxits ( );
                double atol = this->control_.absolute_tol ( );
                double dtol = this->control_.divergence_tol ( );
                this->control_.Init ( maxits, atol, reltol, dtol );
            }

          private:
            LinearSolverState SolveNoPrecond ( const VectorType& b, VectorType* x );
            LinearSolverState SolvePrecondRight ( const VectorType& b, VectorType* x );
            LinearSolverState SolvePrecondLeft ( const VectorType& b, VectorType* x );

        };

        /// @brief BiCGSTAB creator class
        /// @author Simon Gawlok

        template<class LAD>
        class BiCGSTABcreator : public LinearSolverCreator<LAD>
        {
          public:

            LinearSolver<LAD>* params ( const PropertyTree& c )
            {
                BiCGSTAB<LAD>* newBiCGSTAB = new BiCGSTAB<LAD>( );
                if ( c.contains ( "Method" ) )
                    newBiCGSTAB->InitParameter ( c["Method"].template get<std::string>( ).c_str ( ) );
                if ( c.contains ( "MaxIterations" ) && c.contains ( "AbsTolerance" ) &&
                     c.contains ( "RelTolerance" ) && c.contains ( "DivTolerance" ) )
                    newBiCGSTAB->InitControl ( c["MaxIterations"].template get<int>( ),
                                               c["AbsTolerance"].template get<double>( ),
                                               c["RelTolerance"].template get<double>( ),
                                               c["DivTolerance"].template get<double>( ) );
                return newBiCGSTAB;
            }
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_BICGSTAB_H_
