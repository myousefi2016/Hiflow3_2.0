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

/// @author Dimitar Lukarski

#ifndef __LPRECONDITIONER_MC_H
#    define __LPRECONDITIONER_MC_H

#    include <iostream>
#    include <stdlib.h>

#    include "lvector.h"
#    include "lmatrix.h"
#    include "lpreconditioner.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Multi-Coloring Preconditioners for LU types
        /// @author Dimitar Lukarski
        ///
        /// Provides full toolchain for analysing, building and solving LU types of preconditioner

        template <typename ValueType>
        class lPreconditioner_MultiColoring : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_MultiColoring ( );
            virtual ~lPreconditioner_MultiColoring ( );

            /// Default Init/Clear function
            virtual void Init ( void );
            virtual void Clear ( void );

            /// Setup the local vector for permutation if local permutation is used
            /// @param x - Set the unknown vector x for permutation
            /// @param rhs - Set the rhs vector for permutation
            virtual void SetupPermVectors ( hiflow::la::lVector<ValueType> *x,
                                            hiflow::la::lVector<ValueType> *rhs );

            /// Setup the local vector for analyzing
            /// @param Set the unknown vector x for analysis
            virtual void SetupVector ( const hiflow::la::lVector<ValueType> *x );

            /// Setup the local matrix for permutation if local permutation is used
            /// @param Set the matrix for local permutation
            virtual void SetupPermOperator ( hiflow::la::lMatrix<ValueType> &op );

            /// Re-setup the local matrix after external permutation is used
            /// @param Re-setup the matrix
            virtual void reInitOperator ( const hiflow::la::lMatrix<ValueType> &op );

            /// backward permutation for lMatrix
            virtual void PermuteBack ( hiflow::la::lMatrix<ValueType> *op ) const;
            /// backward permutation for lVector
            virtual void PermuteBack ( hiflow::la::lVector<ValueType> *vec ) const;

            /// Set up relaxation parameter
            virtual void SetRelax ( const ValueType omega_relax );

            /// Build function: Analyse() + (permute) + factorize + dropoff + LS + Prepare_LU
            virtual void Build ( void );

            /// Contains build_aux_matrix + multicoloring + allocate_LU
            virtual void Analyse ( void );
            /// Build auxilary matrix for MC
            virtual void build_aux_matrix ( void ) = 0;
            /// Apply multicoloring to the auxilary matrix pattern
            virtual void multicoloring ( void );
            /// Allocate LU matrix
            virtual void allocate_LU ( void ) = 0;

            /// Do a local permutation
            virtual void permute ( void );
            /// Permute Operator + LU matrix + x and rhs
            virtual void permute_all ( void );

            /// Contains factorize + dropoffs + levelscheduling
            virtual void Build_LU ( void );
            /// Factorization
            virtual void factorize ( void ) = 0;
            /// Delete non-diagonal entries in the diagonal blocks
            virtual void dropoffs ( void );
            /// Perform LS analysis
            virtual void levelscheduling ( void );

            /// Contains extract_mat + allocate_vec + convert_Dmat2vec + scale_D + check_for_dropoffs + cleanup
            virtual void Prepare_LU ( void );
            /// Extract D[], L[], R[] sub matrices from LU (wrt to the specific platform of Operator)
            virtual void extract_mat ( void );
            /// Allocate sub vectors
            virtual void allocate_vec ( void );
            /// Convert D[] to vectors
            virtual void convert_Dmat2vec ( void );
            /// Scale diago entries (if need it)
            virtual void scale_D ( void ) = 0;
            /// check for drop-off entries
            virtual void check_for_dropoffs ( void );
            /// Clean up some temporally (internal) data
            virtual void cleanup ( void );

            /// Solve output=(L+D)*D*(U+D)*input
            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

            /// Split input vector into vector chuncks (this->output_chunks_)
            virtual void splitvectors ( const hiflow::la::lVector<ValueType> &input );
            /// Solve the backwrd step
            virtual void backwardstep ( void ) = 0;
            /// Diagonal scaling step
            virtual void diagonalstep ( void ) = 0;
            /// Solve the forward step
            virtual void forwardstep ( void ) = 0;
            /// concatenate output_chunks_ vectors to output vector
            virtual void concatenatevectors ( hiflow::la::lVector<ValueType> *output );

            /// If a global permutation is necessary to be done (instead of local), this function
            /// return the permutation mapping which is allocated from the multicoloring() function.
            /// @return Permutation vector
            virtual int* get_permut ( void );

          protected:

            // color info
            bool flag_mc_; // should mc permutation be used
            bool flag_dropoff_mc_; // dropoff after mc permutation
            bool flag_ls_; // should ls permutation be used

            int ncolors_; // number of colors (levels)
            int *color_sizes_; // color (levels) array sizes

            int *permut_index_; // permutation vector

            // Relaxation Parameter (for SOR, SSOR)
            ValueType omega_relax_;

            // D,L,R matrices (array of matrices)
            lMatrix<ValueType> **L_; // L matrices
            lMatrix<ValueType> **R_; // R matrices
            lMatrix<ValueType> **D_; // D matrices

            // Diagonal Values (array of vectors)
            lVector<ValueType> **Dv_; // Diagonal values
            lVector<ValueType> **iDv_; // Inverse-diagonal values

            // Auxiliary Matrix
            lMatrix<ValueType> *aux_mat_analysis_cpu_;

            // LU Matrix
            lMatrix<ValueType> *LU_cpu_;
            lMatrix<ValueType> *LU_; // LU matrix
            // or pointer to the original matrix for Gauss-Seidel type of preconditioners

            // pointers to the x and rhs (for permutation)
            lVector<ValueType> *x_;
            lVector<ValueType> *rhs_;

            // pointers to x only for analyzing
            const lVector<ValueType> *x_const_;

            // input and output chunks lvectors
            lVector<ValueType> **output_chunks_;

            // Operator for permutation
            hiflow::la::lMatrix<ValueType> *Operator_perm_;

        };

    } // namespace la
} // namespace hiflow

#endif
