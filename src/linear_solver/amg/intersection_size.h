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
/// @date 2015-10-22

#ifndef SRC_LINEAR_SOLVER_AMG_INTERSECTION_SIZE_H_
#    define SRC_LINEAR_SOLVER_AMG_INTERSECTION_SIZE_H_

#    include <algorithm>
#    include <map>

namespace hiflow
{
    namespace AMG
    {

        struct Counter
        {
            typedef size_t const_reference;

            Counter ( ) : count ( 0 )
            {
            }

            struct value_type
            {

                template <typename T> value_type ( const T& )
                {
                }
            };

            void push_back ( const value_type& )
            {
                ++count;
            }
            size_t count;
        };

        /// Return the intersection size of two containers using std::intersection.
        /// Both containers s1, s2 shall be sorted.

        template <typename T1, typename T2>
        size_t intersection_size ( T1 const& s1, T2 const& s2 )
        {
            Counter c;
            std::set_intersection ( s1.begin ( ), s1.end ( ), s2.begin ( ), s2.end ( ), std::back_inserter ( c ) );
            return c.count;
        }

        /// Functor returning true if value of v1 is equal to first value of v2 or vice versa.
        /// Should be a lambda function using C++11.

        struct MapFirstLessThan
        {

            template <class T1, class T2, class T3>
            bool operator() ( T1 const& v1, std::pair<T2, T3> const& v2 ) const
            {
                return v1 < v2.first;
            }

            template <class T1, class T2, class T3>
            bool operator() ( std::pair<T1, T2> const& v1, T3 const& v2 ) const
            {
                return v1.first < v2;
            }
        };

        /// Return the intersection size of two containers using std::intersection.
        /// s1 containers shall be sorted.
        /// s2 is from type std::map.

        template <typename T1, typename T2, typename T3>
        size_t intersection_size ( T1 const& s1, std::map<T2, T3> const& s2 )
        {
            Counter c;
            std::set_intersection ( s1.begin ( ), s1.end ( ), s2.begin ( ), s2.end ( ), std::back_inserter ( c ), MapFirstLessThan ( ) );
            return c.count;
        }

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_INTERSECTION_SIZE_H_ */
