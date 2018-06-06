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

#ifndef SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENING_H_
#    define SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENING_H_

#    include "intersection_size.h"
#    include "StandardCoarseningInterpolationConnection.h"
#    include "StandardCoarseningSettings.h"
#    include "typedefs.h"
#    include "cxx98_prettyprint.h"

//#define PRINT_DEBUG

namespace hiflow
{
    namespace AMG
    {

        /**
         * \brief Standard coarsening
         *
         * Reference: K. Stueben, Algebraic Multigrid (AMG): An Introduction with Applications
         *            GMD-Report 70, November 1999
         */
        template <class LevelGenerator>
        class StandardCoarsening
        {
          public:

            typedef typename LevelGenerator::MatrixType MatrixType;
            typedef typename LevelGenerator::PtrMatrixType PtrMatrixType;
            typedef typename LevelGenerator::ResultType ResultType;
            typedef typename MatrixType::value_type value_type;
            typedef StandardCoarseningInterpolationConnection ConnectionType;
            typedef StandardCoarseningSettings Settings;
            typedef NoOutput Output;

            StandardCoarsening ( Settings const& settings = Settings ( ) )
            : settings_ ( settings )
            {
            }

            ConnectionType operator() ( ResultType &result, MatrixType const& Af ) const
            {
                ConnectionType connection;
                size_t num_variables = Af.get_num_row ( );

                create_strong_connections ( connection.strong_connections, Af );

#    ifdef PRINT_DEBUG
                std::cout << "S = " << connection.strong_connections << std::endl;
#    endif

                AdjacencyList ST = transpose ( connection.strong_connections );

#    ifdef PRINT_DEBUG
                std::cout << "ST = " << ST << std::endl;
#    endif

                // Map undefined variables with lambda value
                std::map<size_t, size_t> U;
                for ( size_t i = 0; i != num_variables; ++i )
                {
                    U.insert ( std::make_pair ( i, ST[i].size ( ) ) );
                }

                // Iterate until all undefined variables are attached
                while ( !U.empty ( ) )
                {
#    ifdef PRINT_DEBUG
                    std::cout << "U = " << U << std::endl;
#    endif

                    size_t next_C = std::max_element ( U.begin ( ), U.end ( ), SecondLessThan ( ) )->first;

#    ifdef PRINT_DEBUG
                    std::cout << "next C = " << next_C << std::endl;
#    endif

                    connection.coarse_variables.insert ( next_C );
                    U.erase ( next_C );

                    // Store indices which must be updated in lambda
                    IndexSet update_list;

                    for ( IndexSet::const_iterator iter_cur ( ST[next_C].begin ( ) ),
                          iter_end ( ST[next_C].end ( ) ); iter_cur != iter_end; ++iter_cur )
                    {
                        if ( U.find ( *iter_cur ) == U.end ( ) ) continue;
#    ifdef PRINT_DEBUG
                        std::cout << "next F = " << *iter_cur << std::endl;
#    endif
                        connection.fine_variables.insert ( *iter_cur );
                        U.erase ( *iter_cur );
                        update_list.insert ( ST[*iter_cur].begin ( ), ST[*iter_cur].end ( ) );
                    }

                    // Assigned indices shall not be updated
                    for ( IndexSet::iterator iter_cur ( update_list.begin ( ) ); iter_cur != update_list.end ( ); )
                    {
                        if ( U.find ( *iter_cur ) == U.end ( ) ) update_list.erase ( iter_cur++ );
                        else ++iter_cur;
                    }

#    ifdef PRINT_DEBUG
                    std::cout << "update_list = " << update_list << std::endl;
#    endif

                    // Update lambda
                    for ( IndexSet::const_iterator iter_cur ( update_list.begin ( ) ), iter_end ( update_list.end ( ) );
                          iter_cur != iter_end; ++iter_cur )
                    {
                        U[*iter_cur] = intersection_size ( ST[*iter_cur], U )
                                + 2 * intersection_size ( ST[*iter_cur], connection.fine_variables );
                    }
                }

#    ifdef PRINT_DEBUG
                IndexVector is;
                std::set_intersection ( connection.coarse_variables.begin ( ), connection.coarse_variables.end ( ),
                                        connection.fine_variables.begin ( ), connection.fine_variables.end ( ),
                                        std::back_inserter ( is ) );
                std::cout << "is: " << is << std::endl;
#    endif

                assert ( num_variables == connection.fine_variables.size ( ) + connection.coarse_variables.size ( ) );

                return connection;
            }

          private:

            /// Determine strong connection matrix as adjacency list

            void create_strong_connections ( AdjacencyList& S, MatrixType const& Af ) const
            {
                size_t num_variables = Af.get_num_row ( );
                S.resize ( num_variables );
                for ( size_t i = 0; i != num_variables; ++i )
                {
                    // Get absolute value of maximal negative element of row i != j
                    value_type abs_max_neg_value ( 0 );
                    for ( size_t p = Af.matrix.row[i]; p != Af.matrix.row[i + 1]; ++p )
                    {
                        if ( Af.matrix.col[p] != i and Af.matrix.val[p] < 0 and std::abs ( Af.matrix.val[p] ) > abs_max_neg_value )
                            abs_max_neg_value = std::abs ( Af.matrix.val[p] );
                    }
                    value_type alpha = abs_max_neg_value * settings_.strength_threshold;

                    // Store strongly connected variables as adjacent list
                    for ( size_t p = Af.matrix.row[i]; p != Af.matrix.row[i + 1]; ++p )
                    {
                        if ( Af.matrix.col[p] != i and -Af.matrix.val[p] >= alpha )
                        {
                            S[i].insert ( Af.matrix.col[p] );
                        }
                    }
                }
            }

            /// Transposition of an adjacency list

            AdjacencyList transpose ( AdjacencyList const& S ) const
            {
                AdjacencyList ST ( S.size ( ) );
                for ( AdjacencyList::const_iterator iter_row_cur ( S.begin ( ) ),
                      iter_row_end ( S.end ( ) ); iter_row_cur != iter_row_end; ++iter_row_cur )
                {
                    for ( AdjacencyList::value_type::const_iterator iter_col_cur ( iter_row_cur->begin ( ) ),
                          iter_col_end ( iter_row_cur->end ( ) ); iter_col_cur != iter_col_end; ++iter_col_cur )
                    {
                        ST[*iter_col_cur].insert ( std::distance ( S.begin ( ), iter_row_cur ) );
                    }
                }
                return ST;
            }

            /// Functor returning true if absolute value of v1 is lesser than absolute value of v2.
            /// Should be a lambda function using C++11.

            struct AbsLessThan
            {

                bool operator() ( double v1, double v2 ) const
                {
                    return std::abs ( v1 ) < std::abs ( v2 );
                }
            };

            /// Functor returning true if second value of v1 is lesser than the second value of v2.
            /// Should be a lambda function using C++11.

            struct SecondLessThan
            {

                template <class T1, class T2>
                bool operator() ( const std::pair<T1, T2>& v1, const std::pair<T1, T2>& v2 ) const
                {
                    return v1.second < v2.second;
                }
            };

            Settings settings_;

        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_STANDARDCOARSENING_H_ */
