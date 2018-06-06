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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-09-29

#ifndef SRC_LINEAR_SOLVER_AMG_TRIPLEPRODUCTIMPL_H_
#    define SRC_LINEAR_SOLVER_AMG_TRIPLEPRODUCTIMPL_H_

#    include "common/smart_pointers.h"
#    include "linear_algebra/coupled_matrix.h"

namespace hiflow
{
    namespace AMG
    {

        /// For local matrices

        template <class MatrixType>
        struct TripleProductImpl
        {
            typedef hiflow::shared_ptr<MatrixType> PtrMatrixType;

            void operator() ( PtrMatrixType &ptr_Ac, MatrixType const& R, MatrixType const& Af, MatrixType const& P ) const
            {
                // Ac = R * Af * P;
                PtrMatrixType tmp ( static_cast < MatrixType* > ( R.MatrixMult ( *( static_cast < const hiflow::la::lMatrix<typename MatrixType::value_type>* > ( &Af ) ) ) ) );
                ptr_Ac = PtrMatrixType ( static_cast < MatrixType* > ( tmp->MatrixMult ( *( static_cast < const hiflow::la::lMatrix<typename MatrixType::value_type>* > ( &P ) ) ) ) );
            }
        };

        /// For coupled matrices

        template <class T>
        struct TripleProductImpl<hiflow::la::CoupledMatrix<T> >
        {
            typedef hiflow::la::CoupledMatrix<T> MatrixType;
            typedef hiflow::shared_ptr<MatrixType> PtrMatrixType;

            void operator() ( PtrMatrixType &ptr_coarse_matrix, MatrixType const& restriction_matrix, MatrixType const& Af, MatrixType const& interpolation_matrix ) const
            {
                // Ac = R * Af * P;
            }
        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_TRIPLEPRODUCTIMPL_H_ */
