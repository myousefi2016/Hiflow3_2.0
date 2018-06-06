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
/// @date 2015-10-20

#ifndef SRC_LINEAR_SOLVER_AMG_TYPEDEFS_H_
#    define SRC_LINEAR_SOLVER_AMG_TYPEDEFS_H_

namespace hiflow
{
    namespace AMG
    {

        typedef std::vector<size_t> IndexVector;
        typedef std::set<size_t> IndexSet;
        typedef std::vector<IndexSet> AdjacencyList;

        /// Empty class if no settings are needed

        struct NoSettings
        {
        };

        /// Empty class if no output shall be defined

        struct NoOutput
        {
        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_TYPEDEFS_H_ */
