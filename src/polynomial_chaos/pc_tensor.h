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

#ifndef HIFLOW_POLYNOMIALCHAOS_PCTENSOR_H_
#    define HIFLOW_POLYNOMIALCHAOS_PCTENSOR_H_

/// \file pc_tensor.h
/// \brief Tensor classes for integrals of Polynomial Chaos.
///
/// \author Michael Schick

#    include "polynomial_chaos/pc_basis.h"
#    include <vector>
#    include <cassert>
#    include <iostream>
#    include <mpi.h>
#    include <map>

namespace hiflow
{
    namespace polynomialchaos
    {

        /// 3rd order Tensor for Polynomial Chaos c=c_{ijk].
        /// Supports sparsity in storage and computation, as well as
        /// distributed memory computation based on MPI.
        /// Also supports multilevel structure of Polynomial Chaos
        /// with respect to the polynomial degree.

        class PCTensor
        {
          public:
            /// Constructor
            PCTensor ( );
            /// Destructor
            ~PCTensor ( );

            /// Set MPI communicator for distributed memory computation
            void SetMPIComm ( const MPI_Comm& comm );
            /// Set threshold for assembly of tensor values (tensor entries
            /// less than threshold are not assembled, default is 1.e-15)

            void SetTensorEps ( double tensor_eps )
            {
                tensor_eps_ = tensor_eps;
            }
            /// Precompute and store values of the tensor
            void ComputeTensor ( int input_degree, PCBasis const& pcbasis );

            /// Returns number of modes on current level

            int Size ( ) const
            {
                return pc_size_[level_];
            }
            /// Returns number of modes on specific level

            int Size ( int level ) const
            {
                assert ( level <= level_ );
                return pc_size_[level];
            }
            /// Returns number of modes owned by a process on the current level

            int SizeLocal ( ) const
            {
                return nmodes_local_[level_];
            }
            /// Returns number of ghost layer modes on the current level

            int SizeOffModes ( ) const
            {
                return noff_modes_[level_];
            }

            /// Return value of tensor at specific position
            /// To be used either with IndicesNL(int pos) or Indices(int pos)

            inline double Val ( int pos ) const
            {
                return val_[level_][pos];
            }
            /// Return value of tensor at specific position (special case linear random input)
            /// To be used together with the function IndicesL(int pos)

            inline double ValL ( int pos ) const
            {
                return val_[level_][pos_linear_[level_][pos]];
            }

            /// Return indices for summation of modes on current level
            /// Local indices in case of parallel computing
            /// It is assumed that the input parameter of the model PDE
            /// exhibits a linear dependence on the uncertainty, which
            /// is approximated by some Polynomial Chaos expansion.

            inline std::vector<int> Indices ( int pos )
            {
                std::vector<int> idx ( 3 );
                idx[0] = glo2loc_[level_][pos][0];
                idx[1] = g2l_[level_].at ( glo2loc_[level_][pos][1] );
                idx[2] = g2l_[level_].at ( glo2loc_[level_][pos][2] );
                return idx;
            }
            /// Return indices for summation of modes on current level
            /// Local indices in case of parallel computing
            /// It is assumed that the input parameter of the model PDE
            /// exhibits a nonlinear dependence on the uncertainty expressed
            /// by some Polynomial Chaos expansion.

            inline std::vector<int> IndicesNL ( int pos )
            {
                std::vector<int> idx ( 3 );
                idx[0] = g2l_[level_].at ( glo2loc_[level_][pos][0] );
                idx[1] = g2l_[level_].at ( glo2loc_[level_][pos][1] );
                idx[2] = g2l_[level_].at ( glo2loc_[level_][pos][2] );
                return idx;
            }
            /// Return indices for summation of modes on current level
            /// Local indices in case of parallel computing (mixed nonlinear case)
            /// It is assumed that two uncertain inputs arise in the
            /// model PDE: one exhibits a linear dependence, the other one a
            /// noninear dependence on the uncertainty. Both expressed by Polynomial
            /// Chaos expansions. The first component of the return vector corresponds
            /// to the global index of the Chaos Polynomial of the linear input, whereas
            /// the other indices refer to the nonlinear case (ansatz and test functions).

            inline std::vector<int> IndicesL ( int pos )
            {
                std::vector<int> idx ( 3 );
                idx[0] = glo2loc_[level_][pos_linear_[level_][pos]][0];
                idx[1] = g2l_[level_].at ( glo2loc_[level_][pos_linear_[level_][pos]][1] );
                idx[2] = g2l_[level_].at ( glo2loc_[level_][pos_linear_[level_][pos]][2] );
                return idx;
            }

            /// Return indices for summation of modes on current level
            /// Global indices in case of parallel computing
            /// In contrast: the function IndicesNL(int pos) or Indices(int pos)
            /// return local indices in parallel computations. For sequential
            /// version, i.e., using only one MPI processor (for Polynomial Chaos)
            /// local and global indices coincide.

            inline std::vector<int> IndicesGlobal ( int pos ) const
            {
                return glo2loc_[level_][pos];
            }

            /// Direct access to L^2 - squared norm of a Polynomial Chaos basis function

            double SquareNorm ( int mode ) const
            {
                return square_norm_[mode];
            }

            /// Returns end index of PCTensor on current level

            inline int End ( ) const
            {
                return nnz_[level_];
            }
            /// Returns k-breakpoints for seperation of PC test functions

            inline std::vector<int> k_breakpoints ( ) const
            {
                return k_breakpoints_[level_];
            }
            /// Returns k-breakpoints for seperation of PC test functions (special case for linear random inputs)

            inline std::vector<int> k_breakpointsL ( ) const
            {
                return k_linear_breakpoints_[level_];
            }

            /// Set level = total polynomial degree

            inline void SetLevel ( int level )
            {
                level_ = level;
            }
            /// Get current level of tensor

            inline int GetLevel ( ) const
            {
                return level_;
            }

            /// Get MPI rank of processor

            inline int MyRank ( ) const
            {
                return my_rank_;
            }
            /// Get number of MPI processors

            inline int NumProc ( ) const
            {
                return numproc_;
            }
            /// Mapping global to local index on current level

            inline int G2L ( int glo_idx ) const
            {
                return g2l_[level_].at ( glo_idx );
            }
            /// Mapping global to local index on specified level

            inline int G2L ( int glo_idx, int level ) const
            {
                return g2l_[level].at ( glo_idx );
            }
            /// Mapping local to global index on current level

            inline int L2G ( int loc_idx ) const
            {
                return l2g_[level_][loc_idx];
            }
            /// Mapping local to global index on specified level

            inline int L2G ( int loc_idx, int level ) const
            {
                return l2g_[level][loc_idx];
            }
            /// Returns ownership (MPI rank) of a global index mode on specified level

            inline int Ownership ( int glo_idx, int level ) const
            {
                return ownership_[level].at ( glo_idx );
            }

            /// MPI communication for sending ghost layer modes on current level

            const std::vector<int>& SendHostileModes ( int rank ) const
            {
                return hostile_modes_send_[level_][rank];
            }
            /// MPI communication for receiving ghost layer modes on current level

            const std::vector<int>& RecvHostileModes ( int rank ) const
            {
                return hostile_modes_needed_[level_][rank];
            }
            /// Access to MPI communicator

            const MPI_Comm& Comm ( ) const
            {
                return comm_;
            }

          protected:
            /// Compute ownership, i.e., determining which MPI processor owns specific
            /// indices corresponding to the Polynomial Chaos (global index)
            void ComputeOwnership ( PCBasis const& pcbasis );
            /// Load balancing of ghost layers (minimize redundant copies, which are
            /// needed in case of parallel computations)
            void ComputeDistributedModesTable ( PCBasis const& pcbasis );

            /// Number of nonzero entries of tensor
            std::vector<int> nnz_;
            /// Size of PC expansion on each level
            std::vector<int> pc_size_;
            /// Polynomial degree of random input
            int q_;
            /// Threshold value for determining if c_{ijk}/c_{0kk} = 0 or not. Default value is 1.e-15
            double tensor_eps_;
            /// Stored nonzero values of tensor
            std::vector<std::vector<double> > val_;
            /// L^2 - square norms of PC basis functions
            std::vector<double> square_norm_;
            /// Mapping global 2 local for indexing scheme
            std::vector<std::vector<std::vector<int> > > glo2loc_;
            /// k-breakpoints, sorting the iteration over the Tensor entries by the testfunctions
            /// global index (the index k of c_{ijk})
            std::vector<std::vector<int> > k_breakpoints_, k_linear_breakpoints_;
            /// position vector for internal access to breakpoints
            std::vector<std::vector<double> > pos_linear_;

            /// Current level, MPI rank and number of processors
            int level_, my_rank_, numproc_;
            /// Number of modes owned by this process, offsets for this process and ghost layer size
            std::vector<int> nmodes_local_, my_mode_offset_, noff_modes_;
            /// Container for ownerships
            std::vector<std::map<int, int> > ownership_;

            /// Global to local mapping
            std::vector<std::map<int, int> > g2l_;
            /// Local to global numbering
            std::vector<std::vector<int> > l2g_;
            /// MPI communicator
            MPI_Comm comm_;
            /// MPI distribution table (required modes)
            std::vector<std::vector<std::vector<int> > > hostile_modes_needed_;
            /// MPI distribution table (sended modes)
            std::vector<std::vector<std::vector<int> > > hostile_modes_send_;
            /// Offset for mode indexing
            std::vector<int> mode_offset_;
        };

    }
}

#endif
