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

#ifndef HIFLOW_MESH_MESH_BUILDER_H
#    define HIFLOW_MESH_MESH_BUILDER_H

#    include <vector>

#    include "mesh/types.h"
#    include "mesh/periodicity_tools.h"

namespace hiflow
{
    namespace mesh
    {

        class Mesh;

        /// \brief Abstract interface for mesh construction.

        class MeshBuilder
        {
          public:
            /// \brief Type for identifiers of vertices that have been added to the MeshBuilder
            typedef int VertexHandle;

            /// \brief Type for identifiers of entities that have been
            /// added to the MeshBuilder. This later becomes the Id values of the entities.
            typedef int EntityHandle;

            /// \brief Constructor for MeshBuilder object.
            /// \param tdim   The topological dimension of the Mesh to be constructed
            /// \param gdim   The geometrical dimension of the Mesh to be constructed
            /// \param period Describing periodicity. Default: No periodicity
            inline MeshBuilder ( TDim tdim, GDim gdim, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Destructor
            inline virtual ~MeshBuilder ( );

            /// \brief Access the topological dimension of the Mesh under construction.
            inline TDim tdim ( ) const;

            /// \brief Access the geometrical dimension of the Mesh under construction.
            inline GDim gdim ( ) const;

            /// \brief Add a vertex to the Mesh under construction.
            /// \param coordinates   the coordinates of the vertex
            /// \return the identifier of the new vertex
            virtual VertexHandle add_vertex ( const std::vector<Coordinate>& coordinates ) = 0;

            /// \brief Add a list of vertices to the Mesh under construction
            /// \param coordinates   the coordinates of the vertices in interleaved format (x0 y0 z0 x1 y1 z1 ...)
            /// \return a list with the identifiers of the new vertices
            virtual std::vector<VertexHandle> add_vertices ( const std::vector<Coordinate>& coordinates ) = 0;

            /// \brief Add an entity to the Mesh under construction
            /// \param tdim       the topological dimension of the added entity
            /// \param vertices   the list of vertices in the entity, identified through the VertexHandle:s returned from add_vertex or add_vertices.
            /// \return  the identifier of the new entity
            virtual EntityHandle add_entity ( TDim tdim, const std::vector<VertexHandle>& vertices ) = 0;

            /// \brief Add a list of entities to the Mesh under construction
            /// \param tdim       the topological dimension of the added entities
            /// \param vertices   the concatenated list of vertices of the entities, identified through their VertexHandle:s.
            /// \param sizes      a list of sizes (number of vertices) of the entities. This length of this list equals the number of entities to be added.
            /// \return a list with the identifiers of the new entities
            virtual std::vector<EntityHandle> add_entities ( TDim tdim,
                                                             const std::vector<VertexHandle>& vertices,
                                                             const std::vector<int>& sizes ) = 0;

            /// \brief Set the material number for an entity.
            /// \param tdim       the topological dimension of the entity
            /// \param entity     identifier for the entity
            /// \param material   new material number for the entity
            virtual void set_material_number ( TDim tdim, EntityHandle entity, MaterialNumber material ) = 0;

            /// \brief Get periodicity
            /// \return Vector containing the information about the periodicity
            /// of the mesh
            inline std::vector<MasterSlave> get_period ( ) const;

            /// \brief Clean the MeshBuilder
            virtual void clear ( ) = 0;

            /// \brief Build the Mesh.
            ///
            /// The indices of the cells in the new mesh are guaranteed to
            /// be the same as that in which they were added to the
            /// mesh. No guarantees are made about entities of lower
            /// dimension, including vertices.
            ///
            /// \return a new Mesh object containing the added
            /// entities. Ownership is transferred to the caller.
            virtual MeshPtr build ( ) = 0;

          protected:
            TDim tdim_;
            GDim gdim_;
            std::vector<MasterSlave> period_;

        };

        //////////////// MeshBuilder implementation ////////////////

        MeshBuilder::MeshBuilder ( TDim tdim, GDim gdim, std::vector<MasterSlave> period )
        : tdim_ ( tdim ), gdim_ ( gdim ), period_ ( period )
        {
        }

        MeshBuilder::~MeshBuilder ( )
        {
        }

        TDim MeshBuilder::tdim ( ) const
        {
            return tdim_;
        }

        GDim MeshBuilder::gdim ( ) const
        {
            return gdim_;
        }

        std::vector<MasterSlave> MeshBuilder::get_period ( ) const
        {
            return period_;
        }
    }
} // namespace hiflow

#endif
