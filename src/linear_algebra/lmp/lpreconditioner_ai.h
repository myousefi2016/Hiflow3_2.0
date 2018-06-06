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

#ifndef __LPRECONDITIONER_AI_H
#    define __LPRECONDITIONER_AI_H

#    include <iostream>
#    include <stdlib.h>

#    include "lvector.h"
#    include "lmatrix.h"
#    include "lmatrix_csr_cpu.h"
#    include "lpreconditioner.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Approximate Inverse Preconditioners
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_ApproximateInverse : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_ApproximateInverse ( );
            virtual ~lPreconditioner_ApproximateInverse ( );

            virtual void Init ( void );
            virtual void Clear ( void );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

          protected:

            lMatrix<ValueType> *AI_;

        };

        /// @brief Approximate Inverse Preconditioners - FSAI
        /// @author Dimitar Lukarski, Niels Wegh
        ///
        /// Factorized Sparse Approximate Inverse Preconditioner - FSAI
        /// min Fob norm of $|I - GA|$, where $G=G_L^T G_L$

        /*
        @ARTICLE{FSAI,
        author = {L.Y. Kolotilina and A.Y. Yeremin},
        title = {Factorized sparse approximate inverse preconditionings, {I}:
                  theory},
        journal = {SIAM J. Matrix Anal. Appl.},
        volume = {14},
        issue = {1},
        month = {January},
        year = {1993},
        issn = {0895-4798},
        pages = {45--58},
        numpages = {14},
        url = {http://dx.doi.org/10.1137/0614004},
        doi = {http://dx.doi.org/10.1137/0614004},
        acmid = {155939},
        publisher = {Society for Industrial and Applied Mathematics},
        address = {Philadelphia, PA, USA},
        keywords = {convergent splittings, preconditioning, sparse approximate
                    inverse},
        }
         */

        template <typename ValueType>
        class lPreconditioner_ApproximateInverse_FSAI : public hiflow::la::lPreconditioner_ApproximateInverse<ValueType>
        {
          public:

            lPreconditioner_ApproximateInverse_FSAI ( );
            virtual ~lPreconditioner_ApproximateInverse_FSAI ( );

            virtual void Init ( void );
            virtual void Init ( const int solver_max_iter,
                                const ValueType solver_rel_eps,
                                const ValueType solver_abs_eps,
                                const int matrix_power );

            /// set up the power matrix pattern (i.e. patter = |matrix|^matrix_power )
            virtual void set_matrix_power ( const int matrix_power );

            /// set the matrix patten via lMatrix
            virtual void set_ext_matrix_pattern ( const hiflow::la::lMatrix<ValueType> &mat );

            // set up the lvector types
            virtual void SetupVector ( const hiflow::la::lVector<ValueType> *x );

            // set up build the approximate invese
            virtual void Build ( void );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

          protected:

            int solver_max_iter_; // max number of iteration for the interal solver (CG)
            ValueType solver_rel_eps_; // relative tol for the interal solver (CG)
            ValueType solver_abs_eps_; // abs for the interal solver (CG)

            int matrix_power_; // matrix power

            lVector<ValueType> *mult_tmp_; // internal lvector for applying the splitting preconditioner

            lMatrix<ValueType> *AI_L_; // L
            lMatrix<ValueType> *AI_Lt_; // L^t

            bool ext_matrix_pat_;
            CPU_CSR_lMatrix<ValueType> *matrix_pat_; // matrix pattern

        };

    } // namespace la
} // namespace hiflow

#endif
