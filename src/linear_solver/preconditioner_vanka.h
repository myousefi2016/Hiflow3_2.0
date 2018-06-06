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

#ifndef HIFLOW_LINEARSOLVER_PRECONDITIONER_VANKA_H_
#    define HIFLOW_LINEARSOLVER_PRECONDITIONER_VANKA_H_

#    include <mpi.h>
#    include <vector>

#    include "config.h"
#    include "mesh/types.h"
#    include "mesh/iterator.h"
#    include "linear_algebra/la_descriptor.h"
#    include "linear_solver/preconditioner_bjacobi.h"
#    include "linear_solver/preconditioner_ilupp.h"
#    include "linear_algebra/seq_dense_matrix.h"
#    include "space/vector_space.h"

namespace hiflow
{
    namespace la
    {
        /// \author Simon Gawlok
        /// \brief Vanka preconditioner interface

        template<class LAD>
        class PreconditionerVanka : public PreconditionerBlockJacobi<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// standard constructor
            PreconditionerVanka ( );
            /// destructor
            virtual ~PreconditionerVanka ( );

            /// Inits parameters for Vanka preconditioner.
            /// \param space Finite element space
            /// \param damping_param Damping parameter
            /// \param num_iter Preconditioner is applied num_iter times; 
            /// this results in a Jacobi iteration over subdomains
            /// \param use_block_of_cells Uses for each cell the block 
            /// of all adjacent cells
            void InitParameter ( const hiflow::VectorSpace<DataType>& space,
                                 const DataType damping_param = 1.,
                                 const int num_iter = 1,
                                 const bool use_block_of_cells = false,
                                 const bool use_preconditioner = false
                                 );

            /// Set configuration parameters for additional ILU++ preconditioner 
            /// param prepro_type type of preprocessing
            /// \param precond_no number of preconditioner
            /// \param max_levels maximum number of multilevels
            /// \param mem_factor see ILU++ manual
            /// \param threshold see ILU++ manual
            /// \param min_pivot see ILU++ manual
            void InitIluppPrecond ( int prepro_type,
                                    int precond_no,
                                    int max_levels,
                                    DataType mem_factor,
                                    DataType threshold,
                                    DataType min_pivot );

            /// Setup the local operator for the local preconditioner
            void SetupOperator ( OperatorType& op );

            /// Build the preconditioner, i.e. perform the factorization and additional time consuming stuff 
            void Build ( );

            /// Applies the Vanka preconditioner.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if preconditioning succeeded
            LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );

            /// Clears allocated data.
            void Clear ( );

          protected:
            /// Creates local matrices
            void CreateLocalMatrices ( );
            /// Computes the  LU factorizations of local matrices
            void Factorize ( );

            const hiflow::VectorSpace<DataType>* space_; // Finite element space

            DataType damping_param_;

            bool use_block_of_cells_; // Uses block of all adjacent cells to "diagonal" cell

            bool use_preconditioner_; // Uses ILU++ preconditioner to get already improved solution
            hiflow::la::PreconditionerIlupp<LAD> precond_; // ILU++ preconditioner

            std::vector<hiflow::la::SeqDenseMatrix<DataType>* > local_mat_diag_; // local matrices (diagonal part)
            hiflow::la::SeqDenseMatrix<DataType> local_mat_diag_block_mode_; // local matrix in block mode

            std::vector<std::vector<int> > sorted_dofs_diag_; // local dofs, sorted for each cell

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PRECONDITIONER_VANKA_H_
