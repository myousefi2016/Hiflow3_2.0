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

#include "connectivity.h"

#include <cassert>

#include "common/log.h"

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace mesh
    {

        typedef Connectivity::ConnectionIterator ConnectionIterator;
        typedef Connectivity::ConstConnectionIterator ConstConnectionIterator;

        bool is_subset ( const ConstConnectionIterator& begin_sub,
                         const ConstConnectionIterator& end_sub,
                         const ConstConnectionIterator& begin_sup,
                         const ConstConnectionIterator& end_sup );

        EntityCount Connectivity::num_entities ( ) const
        {
            return connections_.size ( );
        }

        ConstConnectionIterator Connectivity::begin ( int id ) const
        {
            assert ( id >= 0 );
            assert ( id < static_cast < int > ( connections_.size ( ) ) );
            return connections_[id].begin ( );
        }

        ConnectionIterator Connectivity::begin ( int id )
        {
            assert ( id >= 0 );
            assert ( id < static_cast < int > ( connections_.size ( ) ) );
            return connections_[id].begin ( );
        }

        ConstConnectionIterator Connectivity::end ( int id ) const
        {
            assert ( id >= 0 );
            assert ( id < static_cast < int > ( connections_.size ( ) ) );
            return connections_[id].end ( );
        }

        ConnectionIterator Connectivity::end ( int id )
        {
            assert ( id >= 0 );
            assert ( id < static_cast < int > ( connections_.size ( ) ) );
            return connections_[id].end ( );
        }

        void Connectivity::clear ( )
        {
            connections_.resize ( 0 );
        }

        EntityCount Connectivity::num_connections ( int id ) const
        {
            assert ( id >= 0 );
            assert ( id < static_cast < int > ( connections_.size ( ) ) );
            return connections_[id].size ( );
        }

        void Connectivity::add_connections ( const std::vector<int>& connections )
        {
            connections_.push_back ( connections );
        }

        void Connectivity::transpose ( const EntityCount num_d1_entities,
                                       Connectivity& d1_d2_connectivity ) const
        {
            LOG_DEBUG ( 1, "Transposing full connectivity" );
            std::vector< std::vector<int> >& d1_d2_conn = d1_d2_connectivity.connections_;
            d1_d2_conn.clear ( );
            d1_d2_conn.resize ( num_d1_entities );

            for ( size_t d2_id = 0; d2_id != connections_.size ( ); ++d2_id )
            {
                for ( ConstConnectionIterator d1_it = begin ( d2_id ); d1_it != end ( d2_id ); ++d1_it )
                {
                    d1_d2_conn[*d1_it].push_back ( d2_id );
                }
            }
        }

        void Connectivity::transpose_subset ( const std::vector<int>& d2_indices,
                                              const SortedArray<int>& d1_indices,
                                              Connectivity& d1_d2_connectivity ) const
        {
            LOG_DEBUG ( 1, "Transposing subset of connectivity" );
            typedef std::vector<int>::const_iterator IndexIterator;

            std::vector< std::vector<int> >& d1_d2_conn = d1_d2_connectivity.connections_;

            // create all needed vectors
            d1_d2_conn.clear ( );
            d1_d2_conn.resize ( d1_indices.size ( ) );

            for ( IndexIterator d2_it = d2_indices.begin ( ); d2_it != d2_indices.end ( ); ++d2_it )
            {
                const int d2_id = *d2_it;

                for ( ConstConnectionIterator d1_it = begin ( d2_id ); d1_it != end ( d2_id ); ++d1_it )
                {
                    const int d1_id = *d1_it;
                    int d1_index;
                    if ( d1_indices.find ( d1_id, &d1_index ) )
                    {
                        d1_d2_conn[d1_index].push_back ( d2_id );
                    }
                }
            }
        }

#if 0

        void ConnectivityAlgorithms::intersect_equal_dimensions ( const Connectivity& d_zero_connectivity,
                                                                  const Connectivity& zero_d_connectivity,
                                                                  Connectivity& d_d_connectivity )
        {
            LOG_DEBUG ( 1, "Intersecting (equal dimensions)" );

            const EntityCount num_d_entities = d_zero_connectivity.num_entities ( );

            for ( int d_id = 0; d_id < num_d_entities; ++d_id )
            {
                std::vector<int> connected_d_entities;

                for ( ConstConnectionIterator v_it = d_zero_connectivity.begin ( d_id );
                      v_it != d_zero_connectivity.end ( d_id ); ++v_it )
                {
                    const int v_id = *v_it;
                    for ( ConstConnectionIterator d_it = zero_d_connectivity.begin ( v_id );
                          d_it != zero_d_connectivity.end ( v_id ); ++d_it )
                    {
                        const int connected_d_id = *d_it;

                        const bool is_already_added = find ( connected_d_entities.begin ( ),
                                                             connected_d_entities.end ( ),
                                                             connected_d_id ) != connected_d_entities.end ( );

                        if ( connected_d_id != d_id &&
                             !is_already_added )
                        {
                            connected_d_entities.push_back ( connected_d_id );
                        }
                    }
                }

                d_d_connectivity.add_connections ( connected_d_entities );
            }
        }
#endif

        void ConnectivityAlgorithms::intersect_subset_equal_dimensions (
                                                                         const SortedArray<int>& d_indices,
                                                                         const Connectivity& d_zero_connectivity,
                                                                         const Connectivity& zero_d_connectivity,
                                                                         Connectivity& d_d_connectivity )
        {
            LOG_DEBUG ( 1, "Intersecting subset (equal dimensions)" );
            typedef SortedArray<int>::const_iterator IndexIterator;

            for ( IndexIterator d_it = d_indices.begin ( ); d_it != d_indices.end ( ); ++d_it )
            {
                // NB d_zero_connectivity is assumed to be complete, i.e. contain all entities.
                const int d_id = *d_it;
                std::vector<int> connected_d_entities;

                for ( ConstConnectionIterator v_it = d_zero_connectivity.begin ( d_id );
                      v_it != d_zero_connectivity.end ( d_id ); ++v_it )
                {
                    const int v_id = *v_it;
                    // NB zero_d_connectivity is assumed to be complete - i.e. contain all vertices.
                    // The v_id is then directly the index to use in this connectivity.
                    for ( ConstConnectionIterator cd_it = zero_d_connectivity.begin ( v_id );
                          cd_it != zero_d_connectivity.end ( v_id ); ++cd_it )
                    {
                        // add the connected entity only if it is in d_indices
                        const int connected_d_id = *cd_it;

                        const bool is_already_added = find ( connected_d_entities.begin ( ),
                                                             connected_d_entities.end ( ),
                                                             connected_d_id ) != connected_d_entities.end ( );

                        if ( connected_d_id != d_id &&
                             !is_already_added &&
                             d_indices.find ( connected_d_id, 0 ) )
                        {
                            connected_d_entities.push_back ( connected_d_id );
                        }
                    }
                }
                d_d_connectivity.add_connections ( connected_d_entities );
            }
        }

#if 0

        void ConnectivityAlgorithms::intersect_unequal_dimensions (
                                                                    const Connectivity& d1_zero_connectivity,
                                                                    const Connectivity& d2_zero_connectivity,
                                                                    const Connectivity& zero_d2_connectivity,
                                                                    Connectivity& d1_d2_connectivity )
        {
            LOG_DEBUG ( 1, "Intersecting (unequal dimensions)" );
            const EntityCount num_d1_entities = d1_zero_connectivity.num_entities ( );

            for ( int d1_id = 0; d1_id < num_d1_entities; ++d1_id )
            {
                std::vector<int> connected_d2_entities;
                for ( ConstConnectionIterator v_it = d1_zero_connectivity.begin ( d1_id );
                      v_it != d1_zero_connectivity.end ( d1_id ); ++v_it )
                {
                    const int v_id = *v_it;
                    for ( ConstConnectionIterator d2_it = zero_d2_connectivity.begin ( v_id );
                          d2_it != zero_d2_connectivity.end ( v_id ); ++d2_it )
                    {
                        const int d2_id = *d2_it;

                        const bool is_connected = is_subset ( d2_zero_connectivity.begin ( d2_id ),
                                                              d2_zero_connectivity.end ( d2_id ),
                                                              d1_zero_connectivity.begin ( d1_id ),
                                                              d1_zero_connectivity.end ( d1_id ) );

                        const bool is_already_added = find ( connected_d2_entities.begin ( ),
                                                             connected_d2_entities.end ( ),
                                                             d2_id ) != connected_d2_entities.end ( );

                        if ( is_connected && !is_already_added )
                        {
                            connected_d2_entities.push_back ( d2_id );
                        }
                    }
                }
                d1_d2_connectivity.add_connections ( connected_d2_entities );
            }
        }
#endif

        void ConnectivityAlgorithms::intersect_subset_unequal_dimensions (
                                                                           const SortedArray<int>& d1_indices,
                                                                           const SortedArray<int>& d2_indices,
                                                                           const Connectivity& d1_zero_connectivity,
                                                                           const Connectivity& d2_zero_connectivity,
                                                                           const Connectivity& zero_d2_connectivity,
                                                                           Connectivity& d1_d2_connectivity )
        {
            LOG_DEBUG ( 1, "Intersecting subset (unequal dimensions)" );
            typedef SortedArray<int>::const_iterator IndexIterator;
            std::vector< std::vector<int> >& d1_d2_conn = d1_d2_connectivity.connections_;

            d1_d2_conn.clear ( );

            for ( IndexIterator d1_it = d1_indices.begin ( ); d1_it != d1_indices.end ( ); ++d1_it )
            {
                const int d1_id = *d1_it;

                std::vector<int> connected_d2_entities;

                for ( ConstConnectionIterator v_it = d1_zero_connectivity.begin ( d1_id );
                      v_it != d1_zero_connectivity.end ( d1_id ); ++v_it )
                {
                    const int v_id = *v_it;

                    for ( ConstConnectionIterator d2_it = zero_d2_connectivity.begin ( v_id );
                          d2_it != zero_d2_connectivity.end ( v_id ); ++d2_it )
                    {
                        const int d2_id = *d2_it;

                        const bool is_connected = is_subset ( d2_zero_connectivity.begin ( d2_id ),
                                                              d2_zero_connectivity.end ( d2_id ),
                                                              d1_zero_connectivity.begin ( d1_id ),
                                                              d1_zero_connectivity.end ( d1_id ) );

                        const bool is_already_added = find ( connected_d2_entities.begin ( ),
                                                             connected_d2_entities.end ( ),
                                                             d2_id ) != connected_d2_entities.end ( );

                        if ( is_connected &&
                             !is_already_added &&
                             d2_indices.find ( d2_id, 0 ) )
                        {
                            connected_d2_entities.push_back ( d2_id );
                        }
                    }
                }
                d1_d2_connectivity.add_connections ( connected_d2_entities );
            }
        }

        bool is_subset ( const ConstConnectionIterator& begin_sub,
                         const ConstConnectionIterator& end_sub,
                         const ConstConnectionIterator& begin_sup,
                         const ConstConnectionIterator& end_sup )
        {
            for ( ConstConnectionIterator it = begin_sub; it != end_sub; ++it )
            {
                if ( find ( begin_sup, end_sup, *it ) == end_sup )
                {
                    return false;
                }
            }
            return true;
        }
    }
} // namespace hiflow
