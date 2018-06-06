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

/// \author Simon Gawlok

#ifndef HIFLOW_LINEARSOLVER_SCHUR_COMPLEMENT_H_
#    define HIFLOW_LINEARSOLVER_SCHUR_COMPLEMENT_H_

#    include <mpi.h>
#    include <vector>
#    include <map>

#    include "config.h"
#    include "mesh/types.h"
#    include "mesh/iterator.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_algebra/seq_dense_matrix.h"
#    include "linear_solver/linear_solver.h"
#    include "space/vector_space.h"
#    include "common/log.h"
#    include "common/sorted_array.h"
#    include "nonlinear/nonlinear_problem.h"
#    include "common/timer.h"

namespace hiflow
{
    namespace la
    {

        /// This class provides functionality to precondition a two-block linear
        /// system of equations by the Schur complement method: \n
        /// Given \f$\left(\begin{array}{cc} A & B \\ C & D \end{array}\right)
        /// \left(\begin{array}{cc} x \\ y \end{array}\right)
        /// = \left(\begin{array}{cc} f \\ g \end{array}\right) \f$,
        /// the system is solved by \n
        /// 1. \f$Sy = g-CA^{-1}f\f$ \n
        /// 2. \f$Ax = f - By\f$ \n
        /// with the Schur complement matrix \f$S = D-CA^{-1}B\f$.
        ///
        /// Functionalities to set up the four submatrices and to solve the Schur
        /// complement equation (1.) with the FGMRES method are provided in this class.
        /// For the Schur complement equation, the user can choose, to either
        /// - set up a general preconditioner or
        /// - to define a solver/preconditioner using a matrix, such as the mass or
        ///   the diffusion matrix of the second block of variables.
        ///
        /// Furthermore, the solver can be set to SIMPLE mode, where \f$A^{-1}\f$ is
        /// replaced by \f$diag(A)^{-1}\f$.
        ///
        /// Finally, it is possible to obtain an approximate Schur inverse \f$S^{-1} = s_D*D^{-1} - s_S*Q^{-1} F H^{-1}\f$ for the (Navier)-Stokes equations, 
        /// which is based on a Pressure Convection-Diffusion Preconditioner. \n
        /// To this end, the user needs  to define scaling factors \f$s_D, s_S \f$ and provide solvers \f$D^{-1}, Q^{-1}, H^{-1}\f$ and the matrices \f$F,Q,H\f$.\n 
        /// In order to obtain a good approximation, \f$Q\f$ should be set to the pressure mass matrix, \f$H\f$ should be set to the pressure Laplace matrix with 
        /// appropriate boundary conditions (Purely Neumann Zero BC + regularization Matrix for enclosed flow; Dirichlet Zero BC at outflow boundaries) and 
        /// \f$F\f$ should be set to a pressure mass-Laplace-convection matrix discretizing \f$\partial_t p + \delta_1 \Delta p + \delta_2 u \cdot \nabla p\f$
        /// with convection field \f$u\f$ being the last Newton iterate for the velocity field. 
        ///
        /// Currently, the implementation has two limitations: \n
        /// 1. The setup of the block matrices only works correctly in the case of 
        /// standard FEM, i.e., hp-FEM is currently not supported. \n
        /// 2. The internal solvers are hard-coded and set to those provided by 
        /// the HYPRE software packages. As a consequence, the HYPRE implementation 
        /// of the linear algebra is the only LAD that is supported currently.
        ///

        /// \author Simon Gawlok
        ///
        /// \brief Schur complement solver interface

        template<class LAD>
        class SchurComplement : public LinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType ValueType;

            /// standard constructor
            SchurComplement ( );
            /// destructor
            virtual ~SchurComplement ( );

            /// Init the Schur complement solver. Generates internally the needed 
            /// sparsity structures for the submatrices and the mappings of the system
            /// DoFs to the new block-wise enumeration and vice versa.
            /// \param space Finite element space
            /// \param block_one_variables Variable numbers that belong to the first block
            /// \param block_two_variables Variable numbers that belong to the second block
            virtual void Init ( const hiflow::VectorSpace<ValueType>& space,
                                const std::vector<int>& block_one_variables,
                                const std::vector<int>& block_two_variables );

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

            /// Inits the Schur complement operator i.e. sets up the submatrices.
            /// @param op linear operator to be preconditioned
            virtual void SetupOperator ( OperatorType& op );

            /// Setup operator F
            /// @param op Matrix operator for pressure mass-diffusion-convection
            virtual void SetupOperatorF ( const OperatorType& op );

            /// Setup operator Q
            /// @param op Matrix operator for pressure mass 
            virtual void SetupOperatorQ ( const OperatorType& op );

            /// Setup operator H
            /// @param op Matrix operator for pressure diffusion 
            virtual void SetupOperatorH ( const OperatorType& op );

            /// Setup solver for submatrix A_
            /// @param solver_A Solver object to solve with submatrix A_
            virtual void SetupSolverBlockA ( LinearSolver<LAD>& solver_A );

            /// Setup solver for operator Q
            /// @param solver Solver object to solve operator Q
            virtual void SetupSolverQ ( LinearSolver<LAD>& solver );

            /// Setup solver for operator H
            /// @param solver Solver object to solve operator H
            virtual void SetupSolverH ( LinearSolver<LAD>& solver );

            /// Setup solver for operator D
            /// @param solver Solver object to solve operator D
            virtual void SetupSolverD ( LinearSolver<LAD>& solver );

            /// Setup a preconditioner object for the Schur complement equation.
            /// A potential matrix operator for this preconditioner has to be
            /// provided by the user.
            /// If a general preconditioner is set via this function the block
            /// two matrix based solver/preconditioner is not used.
            /// @param precond Preconditioner object
            virtual void SetupPreconditioner ( Preconditioner<LAD>& precond );

            /// Setup a solver/preconditioner object for the Schur complement
            /// equation, which uses \f$\epsilon M + D\f$ as operator.
            /// The matrix M of the second block has to be provided by
            /// SetupBlockTwoMatrix. This object is not used, if a general
            /// preconditioner is provided by SetupPreconditioner.
            /// @param solver_precond_block_two Solver object to solve matrix solver_precond_block_two_
            virtual void SetupSolverPrecondBlockTwo ( LinearSolver<LAD>& solver_precond_block_two );

            /// Setup matrix for second block of variables; used as preconditioner
            /// in FGMRES iterations for the Schur Complement.
            /// \param[in] op Matrix of global equation system
            virtual void SetupBlockTwoMatrix ( const OperatorType& op );

            /// Setup a SIMPLE preconditioner, i.e., instead of using \f$A^{-1}\f$,
            /// \f$diag(A)^{-1}\f$ is used

            virtual void SetSimplePreconditioning ( bool is_simple_precond )
            {
                this->is_simple_precond_ = is_simple_precond;
            }

            /// Set preconditioning matrix scaling factor for second block
            /// \param[in] scaling Scaling factor for matrix of second block

            virtual void SetBlockTwoMatrixScaling ( ValueType scaling )
            {
                this->scaling_block_two_matrix_ = scaling;
            }

            /// Set scaling factors for approximate Schur complement
            /// @param[in] scaling_D scaling factor for inverse D
            /// @param[in] scaling_S scaling factor for pressure convection-diffusion operator

            virtual void SetPressConvDiffScaling ( const ValueType scaling_D, const ValueType scaling_S )
            {
                this->scaling_D_ = scaling_D;
                this->scaling_S_ = scaling_S;
            }

            /// Set flag indicating whether to use approximate Schur operator based on pressure convection diffusion
            /// @param flag 

            virtual void UsePressConvDiff ( bool flag )
            {
                this->use_press_conv_diff_ = flag;
                if ( flag )
                    this->standard_precond_ = false;
                if ( flag && this->precond_ != NULL )
                {
                    std::cout << "SchurComp: Don't UsePressureConvDiff if preconditioner for Schur complement is already defined. " << std::endl;
                    exit ( -1 );
                }
            }

            /// Set a level ID in case of nested schur complements to make LOG file better readable

            virtual void SetLevel ( int level )
            {
                this->nested_level_ = level;
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

            /// Destroy the solver objects (in case of external libraries). Reason: MPI communicator problems when reusing HYpre linear solvers
            virtual void DestroySolver ( );

          protected:

            /// Evaluate multiplication of Schur complement with a vector
            /// \param[in] x Input vector for matrix-vector-multiplication
            /// \param[out] r Result vector of matrix-vector multiplication
            virtual void SchurVectorMult ( VectorType &x, VectorType &r );

            /// Apply preconditioner for schur complement equation
            /// \param[in] b right hand side vector
            /// \param[out] x solution vector
            virtual void SchurPrecond ( const VectorType &b, VectorType &x, VectorType &rhs_temp, VectorType &sol_temp );

            /// Solve linear system with submatrix A
            /// \param[in] b right hand side vector
            /// \param[out] x solution vector
            virtual void SolveBlockA ( const VectorType &b, VectorType &x );

            /// Applies Givens rotation.
            /// @param cs cos(phi)
            /// @param sn sin(phi)
            /// @param dx first coordinate
            /// @param dy second coordinate

            inline virtual void ApplyPlaneRotation ( const ValueType& cs, const ValueType& sn,
                                                     ValueType* dx, ValueType* dy ) const
            {
                const ValueType temp = cs * ( *dx ) + sn * ( *dy );
                *dy = -sn * ( *dx ) + cs * ( *dy );
                *dx = temp;
            }

            /// Generates Givens rotation.
            /// @param dx first coordinate
            /// @param dy second coordinate
            /// @param cs cos(phi)
            /// @param sn sin(phi)

            inline virtual void GeneratePlaneRotation ( const ValueType& dx, const ValueType& dy,
                                                        ValueType* cs, ValueType* sn ) const
            {
                const ValueType beta = std::sqrt ( dx * dx + dy * dy );
                *cs = dx / beta;
                *sn = dy / beta;
            }

            /// Updates solution: \f$x = x + Vy\f$ with \f$y\f$ solution of least squares problem.
            /// @param V Krylov subspace basis
            /// @param H Hessenberg matrix
            /// @param g rhs of least squares problem
            /// @param k iteration step
            /// @param x solution vector
            void UpdateSolution ( const SeqDenseMatrix<ValueType>& H,
                                  const std::vector<ValueType>& g,
                                  int k,
                                  VectorType* x ) const;

            /// Apply the operator Sp^{-1} = scaling_D * D^{-1} - scaling_S * Q^{-1} * Fp * H^{-1} which is an approximation to the Navier-Stokes Schur complement,
            /// if Q = pressure mass matrix, F_p = pressure mass-diffusion-convection matrix, H = pressure diffusion matrix and D corresponds to block D
            /// \param[in] b right hand side vector
            /// \param[out] x solution vector
            void virtual ApplyPressureConvDiff ( const VectorType &b, VectorType* x );

            /// Apply filtering to a vector, e.g., pressure filtering, in the block 
            /// of the Schur complement equation
            /// \param[inout] vec Vector to be filtered
            /// \param[in] temp Vector for temporary data. Must be initialized to
            /// system vector size
            /// \param[in] indexset Indices of entries in vec in block numbering, 
            /// which belong to part that shall be filtered

            inline virtual void FilterVector ( VectorType& vec, VectorType& temp )
            {
                if ( this->filter_solution_ )
                {
                    Timer timer;
                    timer.reset ( );
                    timer.start ( );

                    vec.GetValues ( vec2ptr ( this->indexset_two_ ), this->indexset_two_.size ( ), vec2ptr ( this->val_temp_two_ ) );
                    temp.SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( this->val_temp_two_ ) );

                    temp.Update ( );

                    this->non_lin_op_->ApplyFilter ( temp );

                    temp.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( this->val_temp_two_ ) );

                    vec.SetValues ( vec2ptr ( this->indexset_two_ ), this->indexset_two_.size ( ), vec2ptr ( this->val_temp_two_ ) );

                    timer.stop ( );
                    if ( this->print_level_ > 2 ) LOG_INFO ( "[" << this->nested_level_ << "] Solution filtering: CPU time                         ", timer.get_duration ( ) );
                }
            }

            /// MPI communicator
            MPI_Comm comm_;
            /// Rank of current process
            int my_rank_;
            /// Global number of processes
            int nb_procs_;

            /// Variable for monitoring the current iteration number in the FGMRES 
            /// solver routine
            int num_iter_;

            /// Number of variables in the first block
            int num_var_one_;
            /// Number of variables in the second block
            int num_var_two_;

            /// Couplings of submatrix \f$A\f$ in system numbering
            std::map<int, SortedArray<int> > couplings_A_;
            /// Couplings of submatrix \f$B\f$ in system numbering
            std::map<int, SortedArray<int> > couplings_B_;
            /// Couplings of submatrix \f$C\f$ in system numbering
            std::map<int, SortedArray<int> > couplings_C_;
            /// Couplings of submatrix \f$D\f$ in system numbering
            std::map<int, SortedArray<int> > couplings_D_;

            /// Sparsity structure of submatrix \f$A\f$ in block numbering
            std::vector<std::vector<int> > sparsity_A_;
            /// Sparsity structure of submatrix \f$B\f$ in block numbering
            std::vector<std::vector<int> > sparsity_B_;
            /// Sparsity structure of submatrix \f$C\f$ in block numbering
            std::vector<std::vector<int> > sparsity_C_;
            /// Sparsity structure of submatrix \f$D\f$ in block numbering
            std::vector<std::vector<int> > sparsity_D_;

            /// Offsets of the processes in the block numbering of the first block
            std::vector<int> offsets_block_one_;
            /// Offsets of the processes in the block numbering of the second block
            std::vector<int> offsets_block_two_;
            /// Mapping block->system numbering of first block
            std::vector<int>mapb2s_one_;
            /// Mapping block->system numbering of second block
            std::vector<int>mapb2s_two_;
            /// Mapping system->block numbering of first block
            std::map<int, int> maps2b_one_;
            /// Mapping system->block numbering of second block
            std::map<int, int> maps2b_two_;

            /// LaCouplings for first block
            LaCouplings la_c_one_;
            /// LaCouplings for second block
            LaCouplings la_c_two_;

            /// Global indices of DoFs in first block
            std::vector<int> indexset_one_;
            /// Global indices of DoFs in second block
            std::vector<int> indexset_two_;

            /// Vector for temporary values in first block
            std::vector<ValueType> val_temp_one_;
            /// Vector for temporary values in second block
            std::vector<ValueType> val_temp_two_;

            /// Submatrix \f$A\f$
            OperatorType A_;
            /// Submatrix \f$B\f$
            OperatorType B_;
            /// Submatrix \f$C\f$
            OperatorType C_;
            /// Submatrix \f$D\f$
            OperatorType D_;

            /// Flag if operator A has been modified
            bool A_modified_;
            /// Flag if operator A has been passed to respective linear solver
            bool A_passed2solver_;
            /// Flag if operator Q has been modified
            bool Q_modified_;
            /// Flag if operator Q has been passed to respective linear solver
            bool Q_passed2solver_;
            /// Flag if operator H has been modified
            bool H_modified_;
            /// Flag if operator H has been passed to respective linear solver
            bool H_passed2solver_;
            /// Flag if operator D has been modified
            bool D_modified_;
            /// Flag if operator D has been passed to respective linear solver
            bool D_passed2solver_;

            /// Operator \f$Q\f$ in Pressure Convection Diffusion \f$S^{-1} = scaling_D*D^{-1} - scaling_S*Q^{-1} F H^{-1}\f$
            OperatorType Q_;
            bool Q_is_initialized_;
            bool Q_is_set_;

            /// Operator \f$F\f$ in Pressure Convection Diffusion \f$S^{-1} = scaling_D*D^{-1} - scaling_S*Q^{-1} F H^{-1}\f$
            OperatorType F_;
            bool F_is_initialized_;
            bool F_is_set_;

            /// Operator \f$H\f$ in Pressure Convection Diffusion \f$S^{-1} = scaling_D*D^{-1} - scaling_S*Q^{-1} F H^{-1}\f$
            OperatorType H_;
            bool H_is_initialized_;
            bool H_is_set_;

            /// Preconditioning matrix \f$M\f$ of second block
            OperatorType block_two_matrix_;
            bool block_two_matrix_is_initialized_;
            bool block_two_matrix_is_set_;

            /// Preconditioning matrix \f$\epsilon M + D\f$ of second block
            OperatorType precond_two_;
            bool precond_two_is_initialized_;
            bool precond_two_is_set_;
            /// Flag if operator precond_two_ has been modified
            bool precond_two_modified_;
            /// Flag if operator precond_two_ has been passed to respective linear solver
            bool precond_two_passed2solver_;

            /// Block of right hand side belonging to the first block of variables
            VectorType f_;
            /// Block of right hand side belonging to the second block of variables
            VectorType g_;

            /// Block of solution vector belonging to the first block of variables
            VectorType x_;
            /// Block of solution vector belonging to the second block of variables
            VectorType y_;

            //*****************************************************************
            // Auxiliary vectors
            //*****************************************************************
            // In block one space
            VectorType h_one_1_, rhs_x_, h1_, h2_;

            // In block two space
            VectorType h_two_1_, w_, h3_;

            // Krylov subspace bases
            std::vector<VectorType*> V_;
            std::vector<VectorType*> Z_;

            /// General preconditioner for the Schur complement equation
            Preconditioner<LAD> *precond_;
            /// If no general preconditioner is set, standard block two matrix-related
            /// preconditioning is done
            bool standard_precond_;
            /// Diagonal of matrix A in case of SIMPLE preconditioning
            std::vector<ValueType> diag_A_;
            /// Use as preconditioner with SIMPLE preconditioning
            bool is_simple_precond_;

            /// If set, ApplyPressureConvDiff is called instead of solver_D_
            bool use_press_conv_diff_;

            /// If set, Schur complement equation is not solved by preconditioend FGMRES. Instead, only the preconditioner is applied 
            bool use_schur_precond_only_;

            /// Pointer to nonlinear problem in order to apply solution filtering
            NonlinearProblem<LAD> *non_lin_op_;
            /// Flag if filtering of solution shall be done
            bool filter_solution_;

            /// Scaling factor for preconditioning matrix of second block
            ValueType scaling_block_two_matrix_;

            /// Flag for level of printing information
            int nested_level_;

            /// Solver for matrix A_
            LinearSolver<LAD>* solver_A_;

            /// Matrix-related solver for the Schur complement equation
            LinearSolver<LAD>* solver_precond_block_two_;

            /// Solver for matrix D_
            LinearSolver<LAD>* solver_D_;

            /// Solver Q
            LinearSolver<LAD>* solver_Q_;

            /// Solver H
            LinearSolver<LAD>* solver_H_;

            /// Scaling
            ValueType scaling_D_;
            ValueType scaling_S_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_SCHUR_COMPLEMENT_H_

