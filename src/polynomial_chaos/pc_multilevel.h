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

#ifndef HIFLOW_POLYNOMIALCHAOS_PCMULTILEVEL_H_
#    define HIFLOW_POLYNOMIALCHAOS_PCMULTILEVEL_H_

/// \file pc_multilevel.h
/// \brief Linear solver and preconditioner for multilevel Polynomial Chaos.
/// It uses the natural hierarchy provided by the total polynomial degree
/// of Polynomial Chaos to define a multilevel approach similar to multigrid
/// in spatial variables.
/// \author Michael Schick

#    include <vector>

#    include "config.h"
#    include "common/property_tree.h"
#    include "linear_solver/linear_solver.h"
#    include "linear_solver/preconditioner.h"
#    include "linear_solver/cg.h"
#    include "linear_solver/gmres.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_algebra/coupled_matrix.h"
#    include "polynomial_chaos/pc_tensor.h"
#    include "linear_solver/preconditioner_bjacobi_standard.h"
#    ifdef WITH_ILUPP
#        include "linear_solver/preconditioner_ilupp.h"
#    endif
#    ifdef WITH_UMFPACK
#        include "linear_solver/umfpack_solver.h"
#    endif
#    include "polynomial_chaos/pc_residual_assembler.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class LAD>
        class PCMultilevelSolver : public la::LinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;
            typedef PCTensor GalerkinTensor;

            /// Default constructor
            PCMultilevelSolver ( const PropertyTree &config );
            /// Destructor
            virtual ~PCMultilevelSolver ( );
            /// Initializes the solvers paramaters.
            virtual void Build ( );
            ///Clear function
            virtual void Clear ( );
            /// Set up operator (matrix) which needs to be solved
            virtual void SetupOperator ( OperatorType& op );
            /// Set up smoother explicitly, provided by a matrix
            virtual void SetupSmoother ( la::CoupledMatrix<DataType>& smoother );
            /// Set up residual assembler for matrix free approach
            virtual void SetupResidualAssembler ( ResidualAssembler<LAD>& assembler );
            /// Set up tensor of Polynomial Chaos

            virtual void SetPCTensor ( GalerkinTensor* pctensor )
            {
                pctensor_ = pctensor;
            }
            /// Update control of iterative solution parameters
            virtual void UpdateControl ( int maxit, double atol, double rtol, double dtol );
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

            /// Applies the multilevel solver as a preconditioner.
            virtual la::LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );
            /// Apply multilevel solver on solution vector x and rhs b to solve Ax = b
            /// where the matrix A is provided either by SetupOperator or by SetupResidualAssembler
            virtual la::LinearSolverState Solve ( const VectorType& b, VectorType* x );

            /// Set Flag if smoothing operator has changed

            virtual void SetModifiedSmoother ( bool flag )
            {
                this->modified_smoother_ = flag;
                if ( flag )
                    this->SetState ( false );
            }

            /// Get flag if smoothing operator has changed

            virtual bool GetModifiedSmoother ( )
            {
                return this->modified_smoother_;
            }

          private:
            /// Apply Multilevel scheme
            virtual void ML ( VectorType const& b, VectorType* x, int l );
            /// Solve deterministic mean problem associated to zero index mode
            virtual void Solve_Mean_Problem ( VectorType const&b, VectorType* x );
            /// Solve deterministic mean problem for a specific mode
            virtual void Solve_Mode_Problem ( int mode, VectorType const&b, VectorType* x );
            /// Apply smoothing method
            virtual void Smoothing ( VectorType const& b, VectorType* x, int nu );
            /// Smoothing parameters
            int nu1_, nu2_, mu_;
            /// Smoother type
            std::string smoother_type_;
            /// Matrix free or not
            std::string mltype_;
            /// Smoothing operator
            la::CoupledMatrix<DataType>* smoother_;
            /// Residual assembler
            ResidualAssembler<LAD>* res_assembler_;
            /// Smoother
            la::LinearSolver<la::LADescriptorCoupledD>* mean_solver_;
            /// Smoother CG
            la::CG<la::LADescriptorCoupledD> mean_solver_cg_;
            /// Smoother GMRES
            la::GMRES<la::LADescriptorCoupledD> mean_solver_gmres_;
#    ifdef WITH_UMFPACK
            /// Smoother Umfpack
            la::UmfpackSolver<la::LADescriptorCoupledD> mean_solver_umf_;
#    endif
#    ifdef WITH_ILUPP
            /// Smoother preconditioner incomplete LU
            la::PreconditionerIlupp<la::LADescriptorCoupledD> ilupp_;
#    endif

            /// Tensor of Polynomial Chaos
            GalerkinTensor* pctensor_;
            /// Deterministic block Jacobi preconditioner
            la::PreconditionerBlockJacobiStand<la::LADescriptorCoupledD> det_precond_;

            /// Configutation parameters
            const PropertyTree* config_;

            /// Globak MPI rank of MPI_COMM_WORLD
            int global_mpi_rank_;

            /// Flag if smoother operator has changed
            bool modified_smoother_;
        };

    } // namespace polynomialchaos
} // namespace hiflow

#endif  // HIFLOW_POLYNOMIALCHAOS_PCMULTILEVEL_H_
