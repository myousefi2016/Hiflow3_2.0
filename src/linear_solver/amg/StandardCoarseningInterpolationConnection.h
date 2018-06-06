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
/// @date 2015-10-06

#ifndef SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENINGINTERPOLATIONCONNECTION_H_
#    define SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENINGINTERPOLATIONCONNECTION_H_

#    include "typedefs.h"

namespace hiflow
{
    namespace AMG
    {

        /**
         * \brief Standard type of the C/F-Splitting.
         *
         * Connection of CoarseningMethod::Output and InterpolationMethod::Input.
         * Only the coarse variables (C) will be stored. The fine variables (F) will be
         * the complement of all variables (Omega) with the fine variables:
         *
         *   F = Omega \ C.
         */
        struct StandardCoarseningInterpolationConnection
        {

            StandardCoarseningInterpolationConnection (
                                                        IndexSet const& coarse_variables = IndexSet ( ),
                                                        IndexSet const& fine_variables = IndexSet ( ),
                                                        AdjacencyList const& strong_connections = AdjacencyList ( )
                                                        ) :
            coarse_variables ( coarse_variables ),
            fine_variables ( fine_variables ),
            strong_connections ( strong_connections )
            {
            }

            /// List of coarse variables
            IndexSet coarse_variables;

            /// List of fine variables in the order of coarse variables
            IndexSet fine_variables;

            /// Strong connection matrix as adjacency list
            AdjacencyList strong_connections;

        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENINGINTERPOLATIONCONNECTION_H_ */
