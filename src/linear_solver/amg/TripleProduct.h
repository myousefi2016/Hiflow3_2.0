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

#ifndef SRC_LINEAR_SOLVER_AMG_TRIPLEPRODUCT_H_
#    define SRC_LINEAR_SOLVER_AMG_TRIPLEPRODUCT_H_

#    include "typedefs.h"
#    include "TripleProductImpl.h"

namespace hiflow
{
    namespace AMG
    {

        /// Implementation of the CoarseningMatrixConstructionMethod as
        /// matrix product Ac = R * Af * P.
        /// R: restriction matrix
        /// P: interpolation matrix
        /// Af: fine matrix
        /// Ac: coarse matrix

        template <class LevelGenerator>
        class TripleProduct
        {
          public:

            typedef typename LevelGenerator::MatrixType MatrixType;
            typedef typename LevelGenerator::PtrMatrixType PtrMatrixType;
            typedef typename LevelGenerator::ResultType ResultType;
            typedef NoSettings Settings;
            typedef NoOutput Output;

            TripleProduct ( Settings const& settings = Settings ( ) )
            : settings_ ( settings )
            {
            }

            void operator() ( ResultType &result, MatrixType const& Af ) const
            {
                TripleProductImpl<MatrixType>( )( result.ptr_coarse_matrix, *result.ptr_restriction_matrix, Af, *result.ptr_interpolation_matrix );
            }

          private:

            Settings settings_;

        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_TRIPLEPRODUCT_H_ */
