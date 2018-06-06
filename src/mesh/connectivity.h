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

/// \author Staffan Ronnas

#ifndef HIFLOW_MESH_CONNECTIVITY_H
#    define HIFLOW_MESH_CONNECTIVITY_H

#    include <vector>

#    include "common/sorted_array.h"
#    include "mesh/types.h"
//#include <mpi.h>
//#include "common/hdf5_tools.h"

namespace hiflow
{
    namespace mesh
    {

        /// Represents the connectivity between the entities of two
        /// dimensions. This implementation uses a directed adjacency
        /// list.

        class Connectivity
        {
          public:
            typedef std::vector<int>::const_iterator ConstConnectionIterator;
            typedef std::vector<int>::iterator ConnectionIterator;

            /// Returns number of entities in this container
            EntityCount num_entities ( ) const;

            /// Returns iterator to first entity connected to entity index
            ConstConnectionIterator begin ( int index ) const;
            ConnectionIterator begin ( int index );

            /// Returns iterator to one past the last entity connected to entity index
            ConstConnectionIterator end ( int index ) const;
            ConnectionIterator end ( int index );

            /// clear Connectivity
            void clear ( );
            /// Returns number of entities connected to entity index
            EntityCount num_connections ( int index ) const;

            /// Adds connections for a new entity
            void add_connections ( const std::vector<int>& connections );

            /// Computes the transposed connectivity with the given Index vector
            void transpose ( const EntityCount num_d1_entities, Connectivity& d1_d2_connectivity ) const;
            void transpose_subset ( const std::vector<int>& d2_indices,
                                    const SortedArray<int>& d1_indices,
                                    Connectivity& d1_d2_connectivity ) const;

            friend class ConnectivityAlgorithms;

            // TODO(Staffan): make this private again if VertexConnectivity is removed
          protected:
            std::vector< std::vector<int> > connections_;
        };

        class ConnectivityAlgorithms
        {
          public:

#    if 0
            // This function is not currently in use. It has never been tested.
            static void intersect_equal_dimensions ( const Connectivity& d_zero_connectivity,
                                                     const Connectivity& zero_d_connectivity,
                                                     Connectivity& d_d_connectivity );
#    endif

            static void intersect_subset_equal_dimensions ( const SortedArray<int>& d_indices,
                                                            const Connectivity& d_zero_connectivity,
                                                            const Connectivity& zero_d_connectivity,
                                                            Connectivity& d_d_connectivity );
#    if 0
            // This function is not currently in use. It has never been tested.
            static void intersect_unequal_dimensions ( const Connectivity& d1_zero_connectivity,
                                                       const Connectivity& d2_zero_connectivity,
                                                       const Connectivity& zero_d2_connectivity,
                                                       Connectivity& d1_d2_connectivity );
#    endif
            static void intersect_subset_unequal_dimensions ( const SortedArray<int>& d1_indices,
                                                              const SortedArray<int>& d2_indices,
                                                              const Connectivity& d1_zero_connectivity,
                                                              const Connectivity& d2_zero_connectivity,
                                                              const Connectivity& zero_d2_connectivity,
                                                              Connectivity& d1_d2_connectivity );
        };
    }
} // namespace hiflow
#endif /* _CONNECTIVITY_H_ */
