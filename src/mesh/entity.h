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

#ifndef HIFLOW_MESH_ENTITY_H
#    define HIFLOW_MESH_ENTITY_H

#    include <iosfwd>
#    include <vector>

#    include "mesh/attributes.h"
#    include "mesh/mesh.h"
#    include "mesh/types.h"

namespace hiflow
{
    namespace mesh
    {

        class CellType;
        class Mesh;
        class IncidentEntityIterator;

        typedef std::vector<Id>::const_iterator ChildrenIdIterator;

        /// \brief View of one entity of a mesh.

        class Entity
        {
          public:

            /// \brief Default constructor. Creates invalid Entity object.
            Entity ( );

            /// \brief Constructor.
            Entity ( ConstMeshPtr mesh, TDim tdim, EntityNumber index );

            // use default copy ctor, assignment operator and dtor

            /// \brief Access topological dimension of the entity.
            TDim tdim ( ) const;

            /// \brief Access geometrical dimension of the entity.
            GDim gdim ( ) const;

            /// \brief Access the id of the entity.
            Id id ( ) const;

            /// \brief Access the index of the entity.
            EntityNumber index ( ) const;

            /// \brief Access the mesh to which the entity belongs.
            ConstMeshPtr mesh ( ) const;

            /// \brief Access the cell type of the entity.
            const CellType& cell_type ( ) const;

            /// \brief Create an iterator to the first incident entity of dimension incident_dim.
            IncidentEntityIterator begin_incident ( TDim incident_dim ) const;

            /// \brief Create an iterator to one-past-the-last incident entity of dimension incident_dim.
            IncidentEntityIterator end_incident ( TDim incident_dim ) const;

            /// \brief Get the number of incident entities of dimension incident_dim.
            EntityCount num_incident_entities ( TDim incident_dim ) const;

            /// \brief Create an iterator to the id of the first vertex.
            VertexIdIterator begin_vertex_ids ( ) const;

            /// \brief Create an iterator to the id of one-past-the-last-vertex.
            VertexIdIterator end_vertex_ids ( ) const;

            /// \brief Get the number of vertices contained in the entity.
            EntityCount num_vertices ( ) const;

            /// \brief Get the id of a vertex in the entity.
            Id vertex_id ( int vertex_number ) const;

            /// \brief Get the parent cell of a cell entity.
            Entity parent ( ) const;

            /// \brief Check whether a cell entity has a parent.
            bool has_parent ( ) const;

            /// \brief Create an iterator to the id of the first child cell of a cell entity.
            ChildrenIdIterator begin_children_ids ( ) const;

            /// \brief Create an iterator to the id of one-past-the-last child cell of a cell entity.
            ChildrenIdIterator end_children_ids ( ) const;

            /// \brief Get the number of children cells of a cell entity.
            EntityCount num_children ( ) const;

            /// \brief Get the id of a child cell of a cell entity
            Id child_id ( int child_number ) const;

            /// \brief Get the coordinates in interleaved format (x0 y0 z0
            /// x1 y1 z1 ...) of the vertices contained in the entity.
            //template<typename T> const std::vector<T> get_coordinates() const;
            template<typename T> void get_coordinates ( std::vector<T>& coords ) const;

            /// \brief Get the coordinates of one vertex of the entity.
            //template<typename T> const std::vector<T> get_coordinates(int vertex_number) const;
            template<typename T> void get_coordinates ( std::vector<T>& coords, int vertex_number ) const;

            /// \brief Get the material number of an entity.
            MaterialNumber get_material_number ( ) const;

            /// \brief Get an attribute of the entity
            template <typename T> void get ( const std::string& name, T* value ) const;

            // \brief Set an attribute of the entity
            template<typename T> void set ( const std::string& name, const T& value ) const;

            // Modify the index of the Entity. This function should primarily
            // be used by the iterator.
            //
            // TODO(Staffan): consider making this private
            // and providing special friend access for iterators

            /// \brief Set the index of an entity. (Internal use only)
            void set_index ( EntityNumber index );

          private:
            /// \brief Fetches coordinates from the associated Mesh object.
            void compute_coordinates ( ) const;

            /// \brief Fetches vertex id:s from the associated Mesh object.
            void compute_vertices ( ) const;

            /// \brief Fetches children cell id:s from the associated Mesh object.
            void compute_children ( ) const;

            /// \brief Check if the Entity object is valid.
            bool is_valid ( ) const;

            /// \brief Mesh object to which Entity belongs.
            ConstMeshPtr mesh_;

            /// \brief Topological dimension of the entity.
            TDim tdim_;

            /// \brief Index of the Entity in mesh_ .
            EntityNumber index_;

            /// \brief Cached coordinate array
            mutable std::vector<Coordinate> coordinates_;

            /// \brief Cached vertex id array
            mutable std::vector<Id> vertices_;

            /// \brief Cached children cell array
            mutable std::vector<Id> children_;
        };

        template<typename T>
        void Entity::get ( const std::string& name, T* value ) const
        {
            mesh_->get_attribute_value ( name, tdim_, index_, value );
        }

        template<typename T>
        void Entity::set ( const std::string& name, const T& value ) const
        {
            mesh_->set_attribute_value ( name, tdim_, index_, value );
        }

        template<typename T>
        void Entity::get_coordinates ( std::vector<T>& coords ) const
        {
            assert ( is_valid ( ) );
            compute_coordinates ( );

            coords.resize ( coordinates_.size ( ) );
            std::copy ( coordinates_.begin ( ), coordinates_.end ( ), coords.begin ( ) );
        }

        template<typename T>
        void Entity::get_coordinates ( std::vector<T>& coords, int vertex_number ) const
        {
            assert ( is_valid ( ) );
            this->compute_coordinates ( );
            const GDim d = gdim ( );

            std::vector<Coordinate>::iterator begin_sub = coordinates_.begin ( ) + d * vertex_number;
            std::vector<Coordinate>::iterator end_sub = coordinates_.begin ( ) + d * ( vertex_number + 1 );

            coords = std::vector<Coordinate>( begin_sub, end_sub );
        }

        /// \brief Output operator for Entity.
        std::ostream& operator<< ( std::ostream& os, const Entity& entity );
    }
} // namespace hiflow

#endif /* _ENTITY_H_ */
