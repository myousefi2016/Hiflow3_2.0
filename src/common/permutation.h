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

#ifndef HIFLOW_COMMON_PERMUTATION_H
#    define HIFLOW_COMMON_PERMUTATION_H

#    include <algorithm>
#    include <cassert>
#    include <vector>

/// \file This file contains a set of template functions to help
/// working with permutations of std::vector.

/// \author Staffan Ronnas

namespace hiflow
{
    /// Comparison operator based on values in a vector. Instead of
    /// comparing indices i and j, it compares values[i] and
    /// values[j]. This can be used e.g. to compute the permutation which
    /// will sort the vector values. To do this, simply call sort on the vector
    /// [0, N), where N is the length of values.

    struct VectorMapCmp
    {

        VectorMapCmp ( const std::vector<int>& values ) : values_ ( values )
        {
        }

        bool operator() ( size_t i, size_t j ) const
        {
            return values_[i] < values_[j];
        }
        const std::vector<int>& values_;
    };

    inline void compute_sorting_permutation ( const std::vector<int>& values,
                                              std::vector<int>& permutation )
    {
        const size_t N = values.size ( );
        permutation.resize ( N );
        for ( size_t i = 0; i < N; ++i )
        {
            permutation[i] = i;
        }

        VectorMapCmp cmp ( values );
        std::stable_sort ( permutation.begin ( ), permutation.end ( ), cmp );
    }

    /// \brief Permute a vector.
    // NB: in_values cannot be equal to out_values

    template<class T>
    void permute_vector ( const std::vector<int>& permutation, const std::vector<T>& in_values,
                          std::vector<T>& out_values )
    {
        const size_t N = in_values.size ( );
        out_values.resize ( N );

        assert ( *std::min_element ( permutation.begin ( ), permutation.end ( ) ) == 0 );
        assert ( *std::max_element ( permutation.begin ( ), permutation.end ( ) ) == N - 1 );

        for ( size_t i = 0; i < N; ++i )
        {
            out_values[i] = in_values[permutation[i]];
        }
    }
}

#endif
