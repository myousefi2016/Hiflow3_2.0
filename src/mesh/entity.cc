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

#include "entity.h"

#include <cassert>
#include <iostream>
#include <algorithm>

#include "cell_type.h"
#include "mesh.h"
#include "iterator.h"

namespace hiflow
{
    namespace mesh
    {

        Entity::Entity ( )
        : mesh_ ( 0 ), tdim_ ( -1 ), index_ ( -1 )
        {
        }

        Entity::Entity ( ConstMeshPtr mesh, TDim tdim, EntityNumber index )
        : mesh_ ( mesh ), tdim_ ( tdim ), index_ ( index )
        {
            if ( tdim == 0 && is_valid ( ) )
            {
                vertices_.push_back ( id ( ) );
            }
        }

        TDim Entity::tdim ( ) const
        {
            assert ( is_valid ( ) );
            return tdim_;
        }

        GDim Entity::gdim ( ) const
        {
            assert ( is_valid ( ) );
            return mesh_->gdim ( );
        }

        Id Entity::id ( ) const
        {
            return mesh_->get_id ( tdim_, index_ );
        }

        EntityNumber Entity::index ( ) const
        {
            assert ( is_valid ( ) );
            return index_;
        }

        ConstMeshPtr Entity::mesh ( ) const
        {
            assert ( is_valid ( ) );
            return mesh_;
        }

        const CellType& Entity::cell_type ( ) const
        {
            assert ( is_valid ( ) );
            return CellType::get_instance ( tdim_, num_vertices ( ) );
        }

        IncidentEntityIterator Entity::begin_incident ( TDim incident_dim ) const
        {
            // Incident_dim has to be smaller or equal than the current
            // entities dimension. If it has to be bigger, use the
            // begin_incident function that is provided in the mesh.
            assert ( is_valid ( ) );
            return mesh_->begin_incident ( *this, incident_dim );
        }

        IncidentEntityIterator Entity::end_incident ( TDim incident_dim ) const
        {
            // Incident_dim has to be smaller or equal than the current
            // entities dimension. If it has to be bigger, use the
            // end_incident function that is provided in the mesh.
            assert ( is_valid ( ) );
            return mesh_->end_incident ( *this, incident_dim );
        }

        EntityCount Entity::num_incident_entities ( TDim incident_dim ) const
        {
            assert ( is_valid ( ) );
            return mesh_->num_incident_entities ( *this, incident_dim );
        }

        VertexIdIterator Entity::begin_vertex_ids ( ) const
        {
            assert ( is_valid ( ) );
            compute_vertices ( );
            return vertices_.begin ( );
        }

        VertexIdIterator Entity::end_vertex_ids ( ) const
        {
            assert ( is_valid ( ) );
            compute_vertices ( );
            return vertices_.end ( );
        }

        EntityCount Entity::num_vertices ( ) const
        {
            assert ( is_valid ( ) );

            if ( tdim_ == 0 )
            {
                return 1;
            }
            else
            {
                compute_vertices ( );
                return vertices_.size ( );
            }
        }

        Id Entity::vertex_id ( int vertex_number ) const
        {
            assert ( is_valid ( ) );

            if ( tdim_ == 0 )
            {
                return id ( );
            }
            else
            {
                compute_vertices ( );

                assert ( vertex_number >= 0 );
                assert ( vertex_number < num_vertices ( ) );
                return vertices_[vertex_number];
            }
        }

        Entity Entity::parent ( ) const
        {
            assert ( tdim_ == mesh_->tdim ( ) );
            assert ( mesh_->cell_has_parent ( index_ ) );
            return mesh_->get_parent_cell ( index_ );
        }

        bool Entity::has_parent ( ) const
        {
            assert ( tdim_ == mesh_->tdim ( ) );
            return mesh_->cell_has_parent ( index_ );
        }

        ChildrenIdIterator Entity::begin_children_ids ( ) const
        {
            assert ( is_valid ( ) );
            assert ( tdim_ == mesh_->tdim ( ) );
            compute_children ( );
            return children_.begin ( );
        }

        ChildrenIdIterator Entity::end_children_ids ( ) const
        {
            assert ( is_valid ( ) );
            assert ( tdim_ == mesh_->tdim ( ) );
            compute_children ( );
            return children_.end ( );
        }

        EntityCount Entity::num_children ( ) const
        {
            assert ( is_valid ( ) );
            assert ( tdim_ == mesh_->tdim ( ) );
            compute_children ( );
            return children_.size ( );
        }

        Id Entity::child_id ( int child_number ) const
        {
            assert ( is_valid ( ) );
            assert ( tdim_ == mesh_->tdim ( ) );
            compute_children ( );

            assert ( child_number >= 0 );
            assert ( child_number < static_cast < int > ( children_.size ( ) ) );

            return children_[child_number];
        }

        MaterialNumber Entity::get_material_number ( ) const
        {
            assert ( is_valid ( ) );
            return mesh_->get_material_number ( tdim_, index_ );
        }

        void Entity::compute_coordinates ( ) const
        {
            assert ( is_valid ( ) );
            if ( coordinates_.empty ( ) )
            {
                coordinates_ = mesh_->get_coordinates ( tdim_, index_ );
            }
        }

        void Entity::compute_vertices ( ) const
        {
            assert ( is_valid ( ) );
            if ( vertices_.empty ( ) )
            {
                if ( tdim_ != 0 )
                {
                    vertices_ = mesh_->get_vertex_ids ( tdim_, index_ );
                }
                else
                {
                    vertices_.push_back ( id ( ) );
                }
            }
        }

        void Entity::compute_children ( ) const
        {
            assert ( is_valid ( ) );
            assert ( tdim_ == mesh_->tdim ( ) );
            if ( children_.empty ( ) )
            {
                children_ = mesh_->get_children_cell_ids ( index_ );
            }
        }

        void Entity::set_index ( EntityNumber index )
        {
            assert ( mesh_ != 0 );

            // clear cached coordinates and vertices
            coordinates_.clear ( );
            vertices_.clear ( );
            children_.clear ( );

            // set the index
            index_ = index;
        }

        bool Entity::is_valid ( ) const
        {
            if ( mesh_ == 0 )
            {
                return false;
            }

            return tdim_ >= 0 && tdim_ <= mesh_->tdim ( ) &&
                    index_ >= 0 && index_ < mesh_->num_entities ( tdim_ );
        }

        std::ostream& operator<< ( std::ostream& os, const Entity& entity )
        {
            os << "Entity " << entity.index ( ) << " [id = " << entity.id ( ) << ", tdim = " << entity.tdim ( ) << ", gdim = " << entity.gdim ( ) << "]\n";
            os << "\tVertices = { ";
            std::copy ( entity.begin_vertex_ids ( ), entity.end_vertex_ids ( ), std::ostream_iterator<Id>( os, " " ) );
            os << "}\n";

            os << "\tVertex coordinates = ";
            for ( int v = 0; v < entity.num_vertices ( ); ++v )
            {
                os << "( ";
                std::vector<Coordinate> coords;
                entity.get_coordinates ( coords, v );
                std::copy ( coords.begin ( ), coords.end ( ), std::ostream_iterator<Coordinate>( os, " " ) );
                os << ")";
            }
            os << "\n";
            return os;
        }
    }
} // namespace hiflow
