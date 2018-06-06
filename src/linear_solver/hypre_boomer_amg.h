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

#ifndef HIFLOW_LINEARSOLVER_HYPRE_BOOMER_AMG_H_
#    define HIFLOW_LINEARSOLVER_HYPRE_BOOMER_AMG_H_

#    include <cstdlib>
#    include <map>

#    include "common/log.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_solver/hypre_linear_solver.h"

namespace hiflow
{
    namespace la
    {
        /// @author Simon Gawlok

        /// @brief Wrapper class for BoomerAMG implementation of Hypre
        /// A linear solver is in particular a preconditioner.

        template<class LAD>
        class HypreBoomerAMG : public HypreLinearSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            HypreBoomerAMG ( );

            ~HypreBoomerAMG ( );

            /// Sets suitable tolerances for usage as preconditioner
            void SetPreconditioningParameters ( );

            /// Solves a linear system.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if solver succeeded
            LinearSolverState Solve ( const VectorType& b, VectorType* x );

            /// Clear allocated data
            void Clear ( );

            /// Initialize solver/preconditioner after parameters have been 
            /// set.
            void Init ( );

#    ifdef WITH_HYPRE
            /// Get pointer to solve function of preconditioner

            HYPRE_PtrToSolverFcn get_solve_function ( )
            {
                return ( HYPRE_PtrToSolverFcn ) HYPRE_BoomerAMGSolve;
            }

            /// Get pointer to setup function of preconditioner

            HYPRE_PtrToSolverFcn get_setup_function ( )
            {
                return ( HYPRE_PtrToSolverFcn ) HYPRE_BoomerAMGSetup;
            }

            /// Get hypre preconditioner object

            HYPRE_Solver& get_solver ( )
            {
                if ( !( this->initialized_ ) || this->modified_param_ )
                {
                    this->Init ( );
                }
                return this->solver_;
            }

            /// call init and setup

            void Build ( );

            /// Destroy solver object 

            void DestroySolver ( )
            {
                HYPRE_BoomerAMGDestroy ( this->solver_ );
                this->SetInitialized ( false );
            }
#    endif

            /// Sets maximum size of coarsest grid. The default is 9.

            void SetMaxCoarseSize ( int max_coarse_size )
            {
                this->max_coarse_size_ = max_coarse_size;
                this->SetModifiedParam ( true );
            }

            /// Sets minimum size of coarsest grid. The default is 1.

            void SetMinCoarseSize ( int min_coarse_size )
            {
                this->min_coarse_size_ = min_coarse_size;
                this->SetModifiedParam ( true );
            }

            /// Sets maximum number of multigrid levels

            void SetMaxLevels ( int max_levels )
            {
                this->max_levels_ = max_levels;
                this->SetModifiedParam ( true );
            }

            /// Defines the number of levels of aggressive coarsening.

            void SetAggNumLevels ( int agg_num_levels )
            {
                this->agg_num_levels_ = agg_num_levels;
                this->SetModifiedParam ( true );
            }

            /// Defines the interpolation type on levels of aggressive coarsening

            void SetAggInterpType ( int agg_interp_type )
            {
                this->agg_interp_type_ = agg_interp_type;
                this->SetModifiedParam ( true );
            }

            /// Defines which parallel coarsening algorithm is used

            void SetCoarsenType ( int coarsen_type )
            {
                this->coarsen_type_ = coarsen_type;
                this->SetModifiedParam ( true );
            }

            /// Sets the size of the system of PDEs, if using the systems version. The default is 1, i.e. a scalar system.

            void SetNumFunctions ( int num_functions )
            {
                this->num_functions_ = num_functions;
                this->SetModifiedParam ( true );
            }

            /// Sets the mapping that assigns the function to each variable, if using the systems version. If no
            /// assignment is made and the number of functions is k > 1, the mapping generated is (0,1,...,k-1,0,1,...,k-1,...).

            void SetDofFunc ( const std::vector<int>& dof_func )
            {
                this->dof_func_.clear ( );
                this->dof_func_ = dof_func;
                this->SetModifiedParam ( true );
            }

            /// Sets AMG strength threshold. The default is 0.25. For 2d Laplace operators, 0.25 is a good value,
            /// for 3d Laplace operators, 0.5 or 0.6 is a better value. For elasticity problems, a large strength threshold,
            /// such as 0.9, is often better.

            void SetStrongThreshold ( DataType strong_threshold )
            {
                this->strong_threshold_ = strong_threshold;
                this->SetModifiedParam ( true );
            }

            /// Set multigrid cycle type
            /// \param[in] cycle_type 1: V-cycle, 2: W-cycle

            void SetCycleType ( int cycle_type )
            {
                this->cycle_type_ = cycle_type;
                this->SetModifiedParam ( true );
            }

            /// Sets whether to use the nodal systems coarsening. Should be used for linear systems generated
            /// from systems of PDEs.
            /// \param[in] nodal type of coarsening algorithm

            void SetNodal ( int nodal )
            {
                this->nodal_ = nodal;
                this->SetModifiedParam ( true );
            }

            /// Sets whether to give special treatment to diagonal elements in the nodal systems version. The
            /// default is 0. If set to 1, the diagonal entry is set to the negative sum of all off diagonal entries. If set to 2,
            /// the signs of all diagonal entries are inverted.

            void SetNodalDiag ( int nodal_diag )
            {
                this->nodal_diag_ = nodal_diag;
                this->SetModifiedParam ( true );
            }

            /// Defines which parallel interpolation operator is used.

            void SetInterpType ( int interp_type )
            {
                this->interp_type_ = interp_type;
                this->SetModifiedParam ( true );
            }

            /// Defines a truncation factor for the interpolation. The default is 0.

            void SetTruncFactor ( DataType trunc_factor )
            {
                this->trunc_factor_ = trunc_factor;
                this->SetModifiedParam ( true );
            }

            /// Defines the maximal number of elements per row for the interpolation. The default is 0.

            void SetPMaxElmts ( int P_max_elem )
            {
                this->P_max_elem_ = P_max_elem;
                this->SetModifiedParam ( true );
            }

            /// Defines the smoother to be used. It uses the given smoother on the fine grid, the up and the
            /// down cycle and sets the solver on the coarsest level to Gaussian elimination.

            void SetRelaxType ( int relax_type )
            {
                this->relax_type_ = relax_type;
                this->SetModifiedParam ( true );
            }

            /// Defines the relaxation weight for smoothed Jacobi and hybrid SOR on all levels.

            void SetRelaxWt ( DataType relax_weight )
            {
                this->relax_weight_ = relax_weight;
                this->SetModifiedParam ( true );
            }

            /// Sets the number of sweeps at a specified cycle.

            void SetCycleNumSweeps ( int num_sweeps, int k )
            {
                this->cycle_num_sweeps_[k] = num_sweeps;
                this->SetModifiedParam ( true );
            }

            /// Defines the smoother at a given cycle.

            void SetCycleRelaxType ( int relax_type, int k )
            {
                this->cycle_relax_type_[k] = relax_type;
                this->SetModifiedParam ( true );
            }

            /// Defines the outer relaxation weight for hybrid SOR and SSOR on all levels.

            void SetOuterWt ( DataType omega )
            {
                this->omega_ = omega;
                this->SetModifiedParam ( true );
            }

            /// Enables the use of more complex smoothers.

            void SetSmoothType ( int smooth_type )
            {
                this->smooth_type_ = smooth_type;
                this->SetModifiedParam ( true );
            }

            /// Set Variant of Schwarz smoother

            void SetVariant ( int variant )
            {
                this->variant_ = variant;
                this->SetModifiedParam ( true );
            }

            /// Set overlap Schwarz smoother

            void SetOverlap ( int overlap )
            {
                this->overlap_ = overlap;
                this->SetModifiedParam ( true );
            }

            /// (Optional) Defines the type of domain used for the Schwarz method

            void SetDomainType ( int domain_type )
            {
                this->domain_type_ = domain_type;
                this->SetModifiedParam ( true );
            }

            /// Sets the number of levels for more complex smoothers.
            /// The smoothers, as defined by SetSmoothType, will be used on level 0 (the finest level) through level
            /// smooth num levels-1. The default is 0, i.e. no complex smoothers are used.

            void SetSmoothNumLevels ( int smooth_num_levels )
            {
                this->smooth_num_levels_ = smooth_num_levels;
                this->SetModifiedParam ( true );
            }

            /// Sets the number of sweeps for more complex smoothers. The default is 1.

            void SetSmoothNumSweeps ( int smooth_num_sweeps )
            {
                this->smooth_num_sweeps_ = smooth_num_sweeps;
                this->SetModifiedParam ( true );
            }

            /// (Optional) Indicates that the aggregates may not be SPD for the 
            /// Schwarz method. The following options exist for use nonsymm: 
            /// 0 - assume SPD (default); 1 - assume non-symmetric

            void SetSchwarzUseNonSymm ( int use_nonsymm )
            {
                this->use_nonsymm_ = use_nonsymm;
                this->SetModifiedParam ( true );
            }

            /// Defines symmetry for ParaSAILS

            void SetSym ( int sym )
            {
                this->sym_ = sym;
                this->SetModifiedParam ( true );
            }

            /// Defines number of levels for ParaSAILS

            void SetLevel ( int level )
            {
                this->level_ = level;
                this->SetModifiedParam ( true );
            }

            /// Defines threshold for ParaSAILS.

            void SetThreshold ( DataType threshold )
            {
                this->threshold_ = threshold;
                this->SetModifiedParam ( true );
            }

            /// Defines filter for ParaSAILS

            void SetFilter ( DataType filter )
            {
                this->filter_ = filter;
                this->SetModifiedParam ( true );
            }

            /// Defines drop tolerance for PILUT.

            void SetDropTol ( DataType drop_tol )
            {
                this->drop_tol_ = drop_tol;
                this->SetModifiedParam ( true );
            }

            /// Defines maximal number of nonzeros for PILUT.

            void SetMaxNzPerRow ( int max_nz_per_row )
            {
                this->max_nz_per_row_ = max_nz_per_row;
                this->SetModifiedParam ( true );
            }

            /// Defines name of an input file for Euclid parameters. For further explanation see description of Euclid.

            void SetEuclidFile ( char* euclidfile )
            {
                this->euclidfile_ = euclidfile;
                this->SetModifiedParam ( true );
            }

            /// Defines number of levels for ILU(k) in Euclid. For further explanation see description of Euclid.

            void SetEuLevel ( int eu_level )
            {
                this->eu_level_ = eu_level;
                this->SetModifiedParam ( true );
            }

            /// Defines filter for ILU(k) for Euclid. For further explanation see description of Euclid.

            void SetEuSparseA ( DataType eu_sparse_A )
            {
                this->eu_sparse_A_ = eu_sparse_A;
                this->SetModifiedParam ( true );
            }

            /// Defines use of block jacobi ILUT for Euclid. For further explanation see description of Euclid.

            void SetEuBJ ( int eu_bj )
            {
                this->eu_bj_ = eu_bj;
                this->SetModifiedParam ( true );
            }

            /// Requests automatic printing of setup and solve information.

            void SetPrintLevel ( int print_level )
            {
                this->print_level_ = print_level;
                this->SetModifiedParam ( true );
            }

            /// If set to 1, the local interpolation transposes will be saved to use more efficient matvecs instead
            /// of matvecTs

            void SetKeepTranspose ( int keepTranspose )
            {
                this->keepTranspose_ = keepTranspose;
                this->SetModifiedParam ( true );
            }

          private:
            /// Maximum number of dofs on coarsest level
            int max_coarse_size_;
            /// Minium number of dofs on coarsest level
            int min_coarse_size_;
            /// Maximum number of multigrid levels
            int max_levels_;
            /// Number of levels of aggressive coarsening
            int agg_num_levels_;
            /// Interpolation on levels of aggressive coarsening
            int agg_interp_type_;
            /// Coarsening operator
            int coarsen_type_;
            /// Number of functions in PDE system
            int num_functions_;
            /// Mapping of dofs to functions
            std::vector<int> dof_func_;
            /// Mapping of dofs to functions (temporary auxiliary vector)
            int* dof_func_temp_;
            /// Strong threshold
            DataType strong_threshold_;
            /// Type of multigrid cycle
            int cycle_type_;
            /// Sets type of nodal treatment
            int nodal_;
            /// Special treatment for diagonal elements in systems version
            int nodal_diag_;
            /// Parallel interpolation operator
            int interp_type_;
            /// Truncation factor for interpolation
            DataType trunc_factor_;
            /// Defines type of relaxation operation
            int relax_type_;
            /// Maximum number of elements per row for interpolation
            int P_max_elem_;
            /// Relaxation weight
            DataType relax_weight_;
            /// Outer relaxation weight
            DataType omega_;
            /// Information for sweeping
            std::map<int, int> cycle_num_sweeps_;
            /// Information for relax type in specific cycles
            std::map<int, int> cycle_relax_type_;
            /// Enables use of more complex smoothers
            int smooth_type_;
            /// Set Variant of Schwarz smoother
            int variant_;
            /// Overlap in Schwarz smoother
            int overlap_;
            /// Type of domain in Schwarz method
            int domain_type_;
            /// Number of levels for more complex smoothers
            int smooth_num_levels_;
            /// Number of sweeps for more complex smoothers. The default is 1.
            int smooth_num_sweeps_;
            /// Indicates that the aggregates may not be SPD for the Schwarz method
            int use_nonsymm_;
            /// Set symmetry for ParaSails
            int sym_;
            /// Number of levels for ParaSails
            int level_;
            /// Threshold for ParaSails
            DataType threshold_;
            /// Filter for ParaSails
            DataType filter_;
            /// Drop tolerance for PILUT
            DataType drop_tol_;
            /// Maximum number of non-zero elements per row for PILUT
            int max_nz_per_row_;
            /// Parameter file for Euclid
            char* euclidfile_;
            /// Number of levels for Euclid
            int eu_level_;
            /// Defines filter for ILU(k) for Euclid
            DataType eu_sparse_A_;
            /// Define use of block jacobi ILUT for Euclid
            int eu_bj_;
            /// Define to keep local interpolation transposes
            int keepTranspose_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_HYPRE_BOOMER_AMG_H_
