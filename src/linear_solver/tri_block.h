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

/// \author Philipp Gerstner

#ifndef HIFLOW_LINEARSOLVER_TRI_BLOCK_H_
#    define HIFLOW_LINEARSOLVER_TRI_BLOCK_H_

#    include <mpi.h>
#    include <vector>
#    include <map>

#    include "config.h"
#    include "mesh/types.h"
#    include "mesh/iterator.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_algebra/seq_dense_matrix.h"
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/schur_complement.h"
#    include "space/vector_space.h"
#    include "common/log.h"
#    include "common/sorted_array.h"
#    include "nonlinear/nonlinear_problem.h"
#    include "common/timer.h"

namespace hiflow
{
    namespace la
    {

        /// This class provides functionality to solve a block triangular or 
        /// block diagonal system of equations of the form: \n
        /// \f$\left(\begin{array}{cc} A & B \\ 0 & D \end{array}\right)
        /// \left(\begin{array}{cc} x \\ y \end{array}\right)
        /// = \left(\begin{array}{cc} f \\ g \end{array}\right) \f$,
        /// or \n
        /// \f$\left(\begin{array}{cc} A & 0 \\ C & D \end{array}\right)
        /// \left(\begin{array}{cc} x \\ y \end{array}\right)
        /// = \left(\begin{array}{cc} f \\ g \end{array}\right) \f$,
        /// or \n
        /// \f$\left(\begin{array}{cc} A & 0 \\ 0 & D \end{array}\right)
        /// \left(\begin{array}{cc} x \\ y \end{array}\right)
        /// = \left(\begin{array}{cc} f \\ g \end{array}\right) \f$. \n
        /// In the upper triangular case, the corresponding system is solved by \n
        /// 1. \f$Dy = g\f$ \n
        /// 2. \f$r  = f-By\f$ \n
        /// 3. \f$Ax = r\f$. \n \n
        ///
        ///
        /// Functionalities to set up the four submatrices and the realisation of the 
        /// defined solution algorithm are provided in this class.
        /// The user has to specify matrix operators from which each submatrix should be extracted and 
        /// linear solvers for block A and block D. It is possible to provide different operators for
        /// each submatrix. However, all matrix operators must have the same coupling structure involving all 
        /// occuring variables, as it is the case for the Schur complement class.\n
        /// In order to setup the solvers, the user can choose, to either 
        /// - provide a solver with already defined operator
        /// - provide a solver with not defined operator. \n
        /// In the later case, the defined solvers are initialized for the corresponding submatrices and the user has 
        /// to call SetupOperator(matrix_op, block_id) for the corresponding subblock (0=A, 1=B, 2=C, 3=D).
        /// If no operator for block B or C is set, they are assumed to be zero. If the user sets an operator for
        /// block A or D, the solvers are initialized for the submatrices extracted from these operators.

        ///
        /// Furthermore, the solver can be set to SIMPLE mode, where \f$A^{-1}\f$ is
        /// replaced by \f$diag(A)^{-1}\f$ and \f$D^{-1}\f$ is
        /// replaced by \f$diag(D)^{-1}\f$ \n\n
        ///
        /// Finally, it is possible to obtain a solver \f$D^{-1} = eE^{-1}-sQ^{-1} Fp H^{-1}\f$. To this end, the user needs
        /// to define scaling factors \f$e,s \f$, solvers \f$E^{-1}, Q^{-1}, H^{-1}\f$ and the matrix \f$Fp\f$.\n \n 
        ///
        /// Currently, the implementation has two limitations: \n
        /// 1. The setup of the block matrices only works correctly in the case of 
        /// standard FEM, i.e., hp-FEM is currently not supported. \n
        /// 2. The internal solvers are hard-coded and set to those provided by 
        /// the HYPRE software packages. As a consequence, the HYPRE implementation 
        /// of the linear algebra is the only LAD that is supported currently.
        ///

        /// \author Philipp Gerstner
        ///
        /// \brief Block Triangular solver interface

        template<class LAD>
        class TriBlock : public SchurComplement<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType ValueType;

            /// standard constructor
            TriBlock ( );
            /// destructor
            virtual ~TriBlock ( );

            /// Init the Schur complement solver. Generates internally the needed 
            /// sparsity structures for the submatrices and the mappings of the system
            /// DoFs to the new block-wise enumeration and vice versa.
            /// \param space Finite element space
            /// \param block_one_variables Variable numbers that belong to the first block
            /// \param block_two_variables Variable numbers that belong to the second block
            virtual void Init ( const hiflow::VectorSpace<ValueType>& space,
                                const std::vector<int>& block_one_variables,
                                const std::vector<int>& block_two_variables );

            /// Set the global matrix operator from which submatrices should be extracted.
            /// Allows the use of different operators for different submatrices.
            /// Non-set B,C - blocks are assumed to be zero. For non-set A,D-blocks, the 
            /// operator of the assigned linear solver is used. 
            /// @param op Matrix operator
            /// @param block ID of submatrix (A=0,B=1,C=2,D=3)
            virtual void SetupOperator ( OperatorType& op, const int block );

            /// Following functions of base class SchurComplement should not be used by TriBlock class

            virtual void SetupOperator ( OperatorType& op )
            {
                std::cout << "TriBlock: Don't use SetupOperator function " << std::endl;
                exit ( -1 );
            }

            virtual void SetupPreconditioner ( Preconditioner<LAD>& precond )
            {
                std::cout << "TriBlock: Don't use SetupPreconditioner function " << std::endl;
                exit ( -1 );
            }

            virtual void SetupSolverPrecondBlockTwo ( LinearSolver<LAD>& solver_precond_block_two )
            {
                std::cout << "TriBlock: Don't use SetupSolverPrecondBlockTwo function " << std::endl;
                exit ( -1 );
            }

            virtual void SetupBlockTwoMatrix ( const OperatorType& op )
            {
                std::cout << "TriBlock: Don't use SetupBlockTwoMatrix function " << std::endl;
                exit ( -1 );
            }

            virtual void SetBlockTwoMatrixScaling ( ValueType scaling )
            {
                std::cout << "TriBlock: Don't use SetBlockTwoMatrixScaling function " << std::endl;
                exit ( -1 );
            }

            /// Setup solver for operator E
            /// @param solver Solver object to solve operator E
            virtual void SetupSolverE ( LinearSolver<LAD>& solver );

            /// Setup operator E
            /// @param op Matrix operator for pressure stabilization 
            virtual void SetupOperatorE ( OperatorType& op );

            /// Set scaling factor for approximate Schur complement
            /// @param scaling scaling factor

            virtual void SetPressConvDiffScaling ( const ValueType scaling_E, const ValueType scaling_S )
            {
                this->scaling_E_ = scaling_E;
                this->scaling_S_ = scaling_S;
            }

            /// Set flag indicating whether to use approximate Schur solver as block D
            /// @param flag 

            virtual void UsePressConvDiff ( bool flag )
            {
                this->use_press_conv_diff_ = flag;
                if ( flag && block_op_[3] > 0 )
                {
                    std::cout << "TriBlock: Don't SetApproximateSchur if Operator D is already defined. " << std::endl;
                    exit ( -1 );
                }
            }

            /// Build the preconditioner, i.e. pass the operators to the subsolvers and build the subsolvers
            virtual void Build ( );

            /// Applies the Schur complement solver.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if preconditioning succeeded
            virtual LinearSolverState Solve ( const VectorType& b, VectorType* x );

            /// Clears allocated data.
            virtual void Clear ( );

          protected:

            /// Solve linear system with submatrix of index block
            /// \param[in] b right hand side vector
            /// \param[out] x solution vector
            virtual void SolveBlock ( const VectorType &b, VectorType &x, const int block );

            /// Apply the filter specified in nonlinear operator to the global vector y. 
            virtual void FilterVector ( VectorType* y );

            /// Apply the filter specified in nonlinear operator to the subvector y. 
            virtual void FilterSubVector ( VectorType* y, VectorType* rhs_temp, std::vector<int>* indexset_two );

            /// Apply the operator Sp^{-1} = scaling_E * E^{-1} - scaling_S * Q^{-1} * Fp * H^{-1} which is an approximation to the Navier-Stokes Schur complement,
            /// if Q = pressure mass matrix, F_p = pressure mass-diffusion-convection matrix, H = pressure diffusion matrix 
            /// \param[in] b right hand side vector
            /// \param[out] x solution vector
            void virtual ApplyPressureConvDiff ( const VectorType &b, VectorType* x );

            /// Operator \f$E\f$ in approximate Schur inverse \f$D^{-1} = eE^{-1}-sQ^{-1} Fp H^{-1}\f$
            OperatorType E_;

            /// Flag if operator E is initialized
            bool E_is_initialized_;

            /// Flag if operator E is set up
            bool E_is_set_;

            /// Flag if operator D has been modified
            bool E_modified_;
            /// Flag if operator D has been passed to respective linear solver
            bool E_passed2solver_;

            /// Flags to indicate whether operators for all blocks are set 
            std::vector<int> block_op_;

            /// Diagonal of matrix A in case of SIMPLE preconditioning
            std::vector<ValueType> diag_D_;

            /// Solver E
            LinearSolver<LAD>* solver_E_;

            /// Scaling
            ValueType scaling_E_;
            ValueType scaling_S_;

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_TRI_BLOCK_H_

