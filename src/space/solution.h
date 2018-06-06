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

#ifndef HIFLOW_SOLUTION_H_
#    define HIFLOW_SOLUTION_H_

#    include "space/vector_space.h"

namespace hiflow
{

    /// @brief Solution structure which represents a vector and its space.
    /// @author Martin Baumann, Chandramowli Subramanian
    ///
    /// A solution consists of a vector and the space in which it is contained.

    template<class DataType>
    class VectorSpace;

    template<class VectorType, class DataType>
    struct Solution
    {
      public:
        const VectorSpace<DataType>* space_;
        const VectorType* vector_;
    };

} // namespace hiflow

#endif  // HIFLOW_SOLUTION_H_
