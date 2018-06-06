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

/// @author Hendryk Bockelmann, Chandramowli Subramanian

#ifndef HIFLOW_LINEARSOLVER_GMRES_H_
#    define HIFLOW_LINEARSOLVER_GMRES_H_

#    include <string>
#    include <vector>
#    include <cmath>
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/linear_solver_creator.h"

namespace hiflow
{
    namespace la
    {

        template<class DataType> class SeqDenseMatrix;

        /// @brief GMRES solver
        ///
        /// GMRES solver for linear systems Ax=b with left, right or no preconditioning.
        /// (not flexible!)

        template<class LAD>
        class GMRES : public LinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            GMRES ( );
            virtual ~GMRES ( );

            virtual void InitParameter ( int size_basis, std::string method );

            virtual LinearSolverState Solve ( const VectorType& b, VectorType* x );

            virtual int size_basis ( ) const
            {
                return this->size_basis_;
            }

            /// Sets the relative tolerance.
            /// Needed by Inexact Newton Methods
            /// @param reltol relative tolerance of residual to converge

            virtual void SetRelativeTolerance ( double reltol )
            {
                int maxits = this->control_.maxits ( );
                double atol = this->control_.absolute_tol ( );
                double dtol = this->control_.divergence_tol ( );
                this->control_.Init ( maxits, atol, reltol, dtol );
            }

          protected:
            /// Applies Givens rotation.
            /// @param cs cos(phi)
            /// @param sn sin(phi)
            /// @param dx first coordinate
            /// @param dy second coordinate

            inline virtual void ApplyPlaneRotation ( const DataType& cs, const DataType& sn,
                                                     DataType* dx, DataType* dy ) const
            {
                const DataType temp = cs * ( *dx ) + sn * ( *dy );
                *dy = -sn * ( *dx ) + cs * ( *dy );
                *dx = temp;
            }

            /// Generates Givens rotation.
            /// @param dx first coordinate
            /// @param dy second coordinate
            /// @param cs cos(phi)
            /// @param sn sin(phi)

            inline virtual void GeneratePlaneRotation ( const DataType& dx, const DataType& dy,
                                                        DataType* cs, DataType* sn ) const
            {
                const DataType beta = std::sqrt ( dx * dx + dy * dy );
                *cs = dx / beta;
                *sn = dy / beta;
            }
            virtual LinearSolverState SolveNoPrecond ( const VectorType& b, VectorType* x );
            virtual LinearSolverState SolveLeft ( const VectorType& b, VectorType* x );
            virtual LinearSolverState SolveRight ( const VectorType& b, VectorType* x );

            virtual void UpdateSolution ( const VectorType * const * V,
                                          const SeqDenseMatrix<DataType>& H,
                                          const std::vector<DataType>& g,
                                          int k,
                                          VectorType* x ) const;

            /// max size of the Krylov subspace basis
            int size_basis_;

        };

        /// @brief GMRES creator class
        /// @author Tobias Hahn

        template<class LAD>
        class GMREScreator : public LinearSolverCreator<LAD>
        {
          public:

            LinearSolver<LAD>* params ( const PropertyTree& c )
            {
                GMRES<LAD>* newGMRES = new GMRES<LAD>( );
                if ( c.contains ( "Method" ) && c.contains ( "SizeBasis" ) )
                {
                    newGMRES->InitParameter ( c["SizeBasis"].template get<int>( ), c["Method"].template get<std::string>( ).c_str ( ) );
                }
                if ( c.contains ( "MaxIterations" ) && c.contains ( "AbsTolerance" ) &&
                     c.contains ( "RelTolerance" ) && c.contains ( "DivTolerance" ) )
                {
                    newGMRES->InitControl ( c["MaxIterations"].template get<int>( ),
                                            c["AbsTolerance"].template get<double>( ),
                                            c["RelTolerance"].template get<double>( ),
                                            c["DivTolerance"].template get<double>( ) );
                }
                return newGMRES;
            }
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_GMRES_H_
