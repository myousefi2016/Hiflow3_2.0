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

/// \author Aksel Alpay, Martin Wlotzka

#ifndef GMG_HIERARCHY_H
#    define GMG_HIERARCHY_H

#    include "level_connection.h"
#    include "basic_hierarchy.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            /// A MultiLevelHierarchy that, additionally to the BasicHierarchy,
            /// creates InterGridConnection objects between the levels.

            template<class ConnectedLevelType>
            class ConnectedHierarchy : public BasicHierarchy<ConnectedLevelType>
            {
              public:

                typedef ConnectedLevelType LevelType;
                typedef typename ConnectedLevelType::Connection ConnectionType;

                IMPORT_FROM_BASECLASS ( BasicHierarchy<ConnectedLevelType>, IteratorFromFinest );
                IMPORT_FROM_BASECLASS ( BasicHierarchy<ConnectedLevelType>, IteratorFromCoarsest );
                IMPORT_FROM_BASECLASS ( BasicHierarchy<ConnectedLevelType>, ConstIteratorFromFinest );
                IMPORT_FROM_BASECLASS ( BasicHierarchy<ConnectedLevelType>, ConstIteratorFromCoarsest );

                /// Initializes the hierarchy. Additionally, creates InterGridConnection
                /// objects between the levels of the hierarchy. Collective on the root
                /// communicator that has been supplied to the communicator hierarchy
                /// generator.
                /// @param master_rank The rank of the master process relative to the
                /// root communicator.
                /// @param levels How many levels the hierarchy is supposed to have
                /// @param master_mesh The master mesh that is to be refined and distributed.
                /// This parameter only needs to be a valid pointer on the master process.
                /// @param fe_degrees Specifies the fe-degrees. The i-th entry of the
                /// vector is the polynomial degree of the i-th variable of the problem
                /// @param settings The settings that are to be used by the linear algebra
                /// routines
                /// @param comm_generator A pointer to a communicator hierarchy generator.

                ConnectedHierarchy ( unsigned master_rank,
                                     unsigned levels,
                                     const MeshPtr& master_mesh,
                                     const std::vector<int>& fe_degrees,
                                     const Settings& settings,
                                     communicator_hierarchy_generators::BasicGenerator* comm_generator )
                : BasicHierarchy<ConnectedLevelType>( master_rank, levels, master_mesh,
                fe_degrees, settings, comm_generator )
                {
                    //we want the last and first connection to remain NULL,
                    //since these levels have no coarser/finer grids.
                    //Thus, our connection vector holds two additional NULL
                    //connections.
                    std::vector<boost::shared_ptr<ConnectionType> > connections ( levels + 1, boost::shared_ptr<ConnectionType>( ) );

                    if ( levels > 1 )
                    {
                        //The first non-NULL connection (from the rear) is at size - 2
                        unsigned position = connections.size ( ) - 2;

                        IteratorFromFinest coarser_grid = this->begin_from_finest_level ( ) + 1;
                        IteratorFromFinest finer_grid = this->begin_from_finest_level ( );
                        for (; coarser_grid != this->end_at_coarsest_level ( );
                              ++coarser_grid, ++finer_grid )
                        {
                            //Initialize connections
                            connections[position] =
                                    boost::shared_ptr<ConnectionType>(
                                    new ConnectionType ( *coarser_grid, *finer_grid ) );

                            --position;
                        }
                    }

                    typename std::vector<boost::shared_ptr<ConnectionType> >::const_iterator
                    connection_iterator = connections.begin ( );

                    for ( IteratorFromCoarsest level_iterator = this->begin_from_coarsest_level ( );
                          level_iterator != this->end_at_finest_level ( );
                          ++level_iterator, ++connection_iterator )
                    {
                        level_iterator->set_connections ( *connection_iterator,
                                                          *( connection_iterator + 1 ) );
                    }
                }
            };

            template<class LAD>
            struct HierarchyTypes
            {
                typedef ConnectedHierarchy<GMGLevel<LAD, GMGConnection<LAD> > > GMGHierarchy;
            };

        } // namespace gmg
    } // namespace la
} // namespace hiflow

#endif
