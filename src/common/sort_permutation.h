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

#ifndef HIFLOW_SORT_PERMUTATION_H
#    define HIFLOW_SORT_PERMUTATION_H

#    include <vector>
#    include <algorithm>
#    include <iostream>

/// \author Simon Gawlok

namespace hiflow
{

    /// @brief Structure for sorting values and get the corresponding permutation.
    /// To achieve this, the comparison for two indices is performed on the 
    /// corresponding values

    template<class T>
    struct CmpPairs
    {
        /// Constructor

        CmpPairs ( const std::vector<T> &v ) : v_ ( v )
        {
        }

        /// Vector of values
        std::vector<T> v_;

        /// Comparison operator

        bool operator() ( int a, int b )
        {
            return v_[a] < v_[b];
        }
    };

    /// @brief Creator function for CmpPairs

    template<class T>
    CmpPairs<T> CreateCmpPairs ( const std::vector<T> & v )
    {
        return CmpPairs<T>( v );
    }

    /// @brief Function to get the corresponding permutation for sorting of
    /// values vector. Permutation vector is allocated inside the function
    /// @param[in] values Vector of values to get the sorting permutation for
    /// @param[out] v Permutation vector

    template<class T>
    void sortingPermutation ( const std::vector<T>& values, std::vector<int>& v )
    {
        size_t size = values.size ( );
        v.clear ( );
        v.reserve ( size );
        for ( size_t i = 0; i != size; ++i )
            v.push_back ( i );

        std::stable_sort ( v.begin ( ), v.end ( ), CreateCmpPairs ( values ) );
        ;

    }

}

#endif
