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

/// @author Philipp Gerstner

#ifndef HIFLOW_EIGEN_VALUE_EIGENSOLVER_H_
#    define HIFLOW_EIGEN_VALUE_EIGENSOLVER_H_

#    include <iostream>
#    include "linear_algebra/la_descriptor.h"

namespace hiflow
{
    namespace la
    {

        /// Enumerator @em EigenSolverState as return value for the eigen value solvers

        enum EigenSolverState
        {
            EigenSolverSuccess = 0,
            EigenSolverExceeded = 1,
            EigenSolverError = 2,
            EigenSolverWrong = 3

        };

        /// @brief Base class for all eigenvalue solvers in HiFlow.

        template<class LAD>
        class EigenSolver
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// Constructor

            EigenSolver ( )
            : n_conv_ ( 0 ),
            print_level_ ( 0 ),
            nev_ ( 0 ),
            ncv_ ( 0 ),
            maxits_ ( 0 ),
            tol_ ( 0. ),
            initialized_ ( false ),
            control_set_ ( false ),
            eigenvalue_set_ ( false ),
            problem_set_ ( false ),
            operator_set_ ( false ),
            complex_A_ ( false ),
            complex_B_ ( false ),
            opA_ ( NULL ),
            opB_ ( NULL ),
            opAt_ ( NULL ),
            opBt_ ( NULL ),
            opAi_ ( NULL ),
            opBi_ ( NULL ),
            opAit_ ( NULL ),
            opBit_ ( NULL )
            {
                res_.resize ( 1, 0. );
            }

            /// Destructor

            virtual ~EigenSolver ( )
            {
            };

            /// Initialize convergence control
            /// @param[in] nev	      number of eigenvalues to be computed
            /// @param[in] ncv        largest dimension of working subspace, default: 2*nev
            /// @param[in] maxits     maximum number of iterations
            /// @param[in] tol        residual tolerance

            void InitControl ( int nev, int ncv, int maxits, double tol )
            {
                maxits_ = maxits;
                tol_ = tol;
                nev_ = nev;
                ncv_ = ncv;
            }

            /// Setup operator for Eigenvalue / SingaluarValue problem
            /// @param[in] opA real part of left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] opAi (optional) imaginary part of left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] cp Coupling object of considered matrix
            /// @param[in] d_rows SparsityStrute.diagonal_rows
            /// @param[in] d_cols SparsityStrute.diagonal_cols
            /// @param[in] o_rows SparsityStrute.off_diagonal_rows
            /// @param[in] o_cols SparsityStrute.off_diagonal_cols

            virtual void SetupOperatorA ( const OperatorType* opA, const OperatorType* opAi,
                                          const Couplings<DataType>& cp, const std::vector<int>& d_rows, const std::vector<int>& d_cols, const std::vector<int>& o_rows, const std::vector<int>& o_cols )
            {
            }

            /// Setup operator for Eigenvalue / SingaluarValue problem with matrix free implementation.
            /// Note: not every eigensolver needs the transposed operator. When using such a solver, use opAt = opA.
            /// @param[in] opA real part of left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] opAt real part of transposed left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] opAi (optional) imag part of left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] opAit (optional) imag part of transposed left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD

            virtual void SetupOperatorA ( const OperatorType* opA, const OperatorType* opAt, const OperatorType* opAi, const OperatorType* opAit, const Couplings<DataType>& cp )
            {
            }

            /// Setup operator for Eigenvalue / SingaluarValue problem
            /// @param[in] opB real part of right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] opBi (optional) imaginary part of right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] cp Coupling object of considered matrix
            /// @param[in] d_rows SparsityStrute.diagonal_rows
            /// @param[in] d_cols SparsityStrute.diagonal_cols
            /// @param[in] o_rows SparsityStrute.off_diagonal_rows
            /// @param[in] o_cols SparsityStrute.off_diagonal_cols

            virtual void SetupOperatorB ( const OperatorType* opB, const OperatorType* opBi,
                                          const Couplings<DataType>& cp, const std::vector<int>& d_rows, const std::vector<int>& d_cols, const std::vector<int>& o_rows, const std::vector<int>& o_cols )
            {
            }

            /// Setup operator for Eigenvalue / SingaluarValue problem with matrix free implementation.
            /// Note: not every eigensolver needs the transposed operator. When using such a solver, use opAt = opA.
            /// @param[in] opB real part of right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] opBt real part of transposed right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] opBi (optional) imag part of right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD
            /// @param[in] opBit (optional) imag part of transposed right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD

            virtual void SetupOperatorB ( const OperatorType* opB, const OperatorType* opBt, const OperatorType* opBi, const OperatorType* opBit, const Couplings<DataType>& cp )
            {
            }

            /// Call this functions after setting all operators defining the eigenvalue problem
            virtual void PassOperatorsToSolver ( ) = 0;

            /// Apply the solver
            virtual EigenSolverState Solve ( ) = 0;

            /// Return residual of specific eigenvalue
            /// @param[in] j index of converged eigenvalue

            DataType Res ( int j ) const
            {
                return this->res_[j];
            }

            /// Return number of iterations

            int Iter ( ) const
            {
                return this->iter_;
            }

            /// Return number of converged eigenvalues

            int ConvEigValues ( ) const
            {
                return this->n_conv_;
            }

            /// Erase possibly allocated data.

            void Clear ( )
            {
                opA_ = NULL;
                opB_ = NULL;
                opAt_ = NULL;
                opBt_ = NULL;
                opAi_ = NULL;
                opBi_ = NULL;
                opAit_ = NULL;
                opBit_ = NULL;
                n_conv_ = 0;
                print_level_ = 0;
                operator_set_ = false;
                control_set_ = false;
                problem_set_ = false;
                eigenvalue_set_ = false;
                initialized_ = false;
                complex_A_ = false;
                complex_B_ = false;
                nev_ = 0;
                ncv_ = 0;
                maxits_ = 0;
                tol_ = 0.;
                iter_ = 0.;
                res_.clear ( );
            }

            /// Print possibly info

            virtual void Print ( std::ostream &out = std::cout ) const
            {
            };

            /// Set Print level
            /// @param[in] level print level

            void SetPrintLevel ( int level )
            {
                print_level_ = level;
            }

          protected:
            /// Left system matrix (real part)
            OperatorType* opA_;

            /// Right system matrix in case of generalized eigenvalue problem (real part)
            OperatorType* opB_;

            /// Transposed left system matrix (real part)
            OperatorType* opAt_;

            /// Transposed right system matrix in case of generalized eigenvalue problem (real part)
            OperatorType* opBt_;

            /// Left system matrix, imaginary part
            OperatorType* opAi_;

            /// Right system matrix in case of generalized eigenvalue problem, imaginary part
            OperatorType* opBi_;

            /// Transposed left system matrix, imaginary part
            OperatorType* opAit_;

            /// Transposed right system matrix in case of generalized eigenvalue problem, imaginary part
            OperatorType* opBit_;

            /// residual norm for each eigenpair
            std::vector<DataType> res_;

            /// Maximum number of iterations
            int maxits_;

            /// Number of eigenvalues to compute
            int nev_;

            /// Maximum dimension of Krylov subspace
            int ncv_;

            /// Maximum dimension of projection matrix
            int mpd_;

            /// Tolerance
            DataType tol_;

            /// Iteration counter
            int iter_;

            /// Number of computed eigenpairs
            int n_conv_;

            /// Flag if solver is initialized
            bool initialized_;

            /// Flag if problem type is set
            bool problem_set_;

            /// Flag if eigenvalue type is set
            bool eigenvalue_set_;

            /// Flag if control is set
            bool control_set_;

            /// Flag if operators are set
            bool operator_set_;

            /// Print level for log file
            int print_level_;

            /// Flag if operator A is complex valued
            bool complex_A_;

            /// Flag if operator B is complex valued
            bool complex_B_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PRECONDITIONER_H_
