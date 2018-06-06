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

/// \author Thomas Gengenbach and Staffan Ronnas

#ifndef HIFLOW_MESH_ITERATOR_H
#    define HIFLOW_MESH_ITERATOR_H

#    include <boost/iterator/iterator_facade.hpp>
#    include <vector>

#    include "mesh/entity.h"
#    include "mesh/types.h"

namespace hiflow
{
    namespace mesh
    {

        class Mesh;

        /// \brief Provides different iterators to get entities.
        /// \details EntityIterator: Forward Traversal Iterator
        /// Implemented through boost::iterator_facade

        class EntityIterator
        : public boost::iterator_facade <EntityIterator,
        const Entity,
        boost::forward_traversal_tag>
        {
          public:

            EntityIterator ( )
            : entity_ ( ), index_ ( -1 )
            {
            }

            EntityIterator ( const Entity& entity, EntityNumber index )
            : entity_ ( entity ), index_ ( index )
            {
            }

          private:
            // provide access for boost::iterator_facade
            friend class boost::iterator_core_access;

            // interface for boost::iterator_facade

            bool equal ( const EntityIterator& other ) const
            {
                return index_ == other.index_;
            }

            void increment ( )
            {
                ++index_;
            }

            const Entity& dereference ( ) const
            {
                entity_.set_index ( index_ );
                return entity_;
            }

            mutable Entity entity_;
            EntityNumber index_;
        };

        /// \brief Provides iterator to get entities that are incident to
        // another entity..
        /// \details IncidentEntityIterator: Forward Traversal Iterator
        /// Implemented through boost::iterator_facade

        class IncidentEntityIterator
        : public boost::iterator_facade <IncidentEntityIterator,
        const Entity,
        boost::forward_traversal_tag>
        {
          public:

            IncidentEntityIterator ( )
            {
            }

            IncidentEntityIterator ( const Entity& entity, std::vector<EntityNumber>::const_iterator index_it )
            : entity_ ( entity ), index_it_ ( index_it )
            {
            }

            const Entity& get_entity ( ) const
            {
                return entity_;
            }

            std::vector<EntityNumber>::const_iterator get_incidence_iterator ( ) const
            {
                return index_it_;
            }
          private:
            // provide access for boost::iterator_facade
            friend class boost::iterator_core_access;

            // interface for boost::iterator_facade

            bool equal ( const IncidentEntityIterator& other ) const
            {
                return index_it_ == other.index_it_;
            }

            void increment ( )
            {
                ++index_it_;
            }

            const Entity& dereference ( ) const
            {
                entity_.set_index ( *index_it_ );
                return entity_;
            }

            mutable Entity entity_;
            std::vector<EntityNumber>::const_iterator index_it_;
        };
    }
} // namespace hiflow
#endif /* _MESH_ITERATOR_H_ */
