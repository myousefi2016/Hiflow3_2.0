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
/// @date 2015-10-07

#ifndef SRC_LINEAR_SOLVER_AMG_PRETENDCOARSENING_H_
#    define SRC_LINEAR_SOLVER_AMG_PRETENDCOARSENING_H_

#    include "StandardCoarseningInterpolationConnection.h"
#    include "typedefs.h"

namespace hiflow
{
    namespace AMG
    {

        /**
         * The coarse variables are manually pretended.
         * Main purpose is for testing.
         */
        template <class LevelGenerator>
        class PretendCoarsening
        {
          public:

            typedef typename LevelGenerator::MatrixType MatrixType;
            typedef typename LevelGenerator::PtrMatrixType PtrMatrixType;
            typedef typename LevelGenerator::ResultType ResultType;
            typedef StandardCoarseningInterpolationConnection ConnectionType;
            typedef StandardCoarseningInterpolationConnection Settings;
            typedef NoOutput Output;

            PretendCoarsening ( Settings const& settings )
            : settings_ ( settings )
            {
            }

            ConnectionType operator() ( ResultType &result, MatrixType const& Af ) const
            {
                return ConnectionType (
                                        settings_.coarse_variables,
                                        settings_.fine_variables,
                                        settings_.strong_connections
                                        );
            }

          private:

            Settings settings_;

        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_PRETENDCOARSENING_H_ */
