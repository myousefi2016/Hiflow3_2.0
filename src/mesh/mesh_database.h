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

#ifndef MESH_DATABASE_H
#    define MESH_DATABASE_H

#    include <vector>
#    include <map>
#    include <iostream>
#    include <string>

#    include "common/sorted_array.h"
#    include "communication.h"
#    include "connectivity.h"
#    include "types.h"
#    include <mpi.h>
#    include "common/hdf5_tools.h"

namespace hiflow
{
    namespace mesh
    {

        class VertexConnectivity;
        class VertexSearchTable;
        class MeshDatabase;

        typedef SharedPtr<MeshDatabase>::Type MeshDatabasePtr;
        typedef std::vector<int>::const_iterator ConstConnectionIterator;
        typedef std::vector<int>::iterator ConnectionIterator;

        ///
        /// \brief Container class for mesh entities.
        ///
        /// The MeshDatabase class stores topological and geometrical data
        /// for mesh entities. It manages two types of objects: vertices,
        /// which are mesh entities of topological dimension 0 together
        /// with the position of the vertex in geometrical space; and
        /// entities, which encompass edges, faces (in 3d) and cells. Each
        /// entity is defined by the ordered sequence of vertices that it
        /// contains.
        ///
        /// Each vertex and entity is assigned a id (of type Id). The id
        /// is always non-negative, but otherwise no guarantee is given as
        /// to the order in which id:s are assigned. The id is unique for
        /// entities of a given dimension.
        ///
        /// The mesh database has a topological and a geometrical
        /// dimension. The topological dimension of the mesh database is
        /// the maximum topological dimension of the entities that it can
        /// contain. The geometrical dimension is the number of
        /// coordinates of the vertices. The geometrical dimension is
        /// always greater or equal to the topological dimension, but it
        /// is possible to have the mesh database contain the entities of
        /// e.g. a two-dimensional mesh in three-dimensional space.
        ///
        ///
        /// \author Staffan Ronnas

        class MeshDatabase : public NonCopyable
        {
          public:
            /// \brief Constructor
            MeshDatabase ( TDim tdim, GDim gdim );

            ~MeshDatabase ( )
            {
            }

            /// \brief Returns the topological dimension of the mesh database
            virtual TDim tdim ( ) const;

            /// \brief Returns the geometrical dimension of the mesh database
            virtual GDim gdim ( ) const;

            /// \brief Adds a vertex, creating it if it does not already exist
            virtual Id add_vertex ( const Coordinate* point );

            /// \brief Replace a vertex
            virtual void replace_vertex ( const Coordinate* destination, Id vertex_id );

            /// \brief Move a vertex
            virtual void move_vertex ( const Coordinate* displacement, Id vertex_id );

            /// \brief Check if vertex with given coordinates exists
            virtual bool has_vertex ( const Coordinate* point, Id& id ) const;

            /// \brief Access the coordinates of a vertex
            virtual std::vector<Coordinate> get_coordinates ( Id vertex_id ) const;

            /// \brief Access the coordinates of a list of vertices
            virtual std::vector<Coordinate> get_coordinates ( const std::vector<Id>& vertex_ids ) const;

            /// \brief Get point that exists in database that is closest
            /// to x within an epsilon range eps.
            virtual bool find_closest_vertex ( const Coordinate* x, double eps, Id* closest_vertex ) const;

            /// \brief Access the number of vertices
            virtual EntityCount num_vertices ( ) const;

            /// \brief Adds an entity, creating it if it does not already exist
            virtual Id add_entity ( TDim entity_dim, const std::vector<Id>& entity_vertices );

            /// \brief Check if entity with given vertices exists
            virtual bool has_entity ( TDim entity_dim, const std::vector<Id>& entity_vertices, Id& id ) const;

            /// \brief Check if entity with given Id exists
            virtual bool has_entity ( TDim entity_dim, Id id ) const;

            /// \brief Access the number of entities in the database
            virtual EntityCount num_entities ( TDim entity_dim ) const;

            /// \brief Obtain iterator to first vertex connected to an entity
            virtual VertexIdIterator begin_vertices ( TDim entity_dim, Id id ) const;

            /// \brief Obtain iterator to one-past-the-last vertex connected to an entity
            virtual VertexIdIterator end_vertices ( TDim entity_dim, Id id ) const;

            /// \brief Obtain the number of vertices connected to an entity
            virtual EntityCount entity_size ( TDim entity_dim, Id id ) const;

            /// \brief Creates sub-entities of cells.
            /// 0 < d < D
            virtual void build ( TDim d, const SortedArray<Id>& cells,
                                 SortedArray<Id>& d_entities,
                                 Connectivity& cell_d_connections );

            /// \brief Computes incidence connections between existing entities
            /// 0 < d2 <= d1 <= D
            virtual void compute_incident_entities ( TDim d, const SortedArray<Id>& d_entities,
                                                     Connectivity& d_d_connectivity ) const;

            virtual void compute_incident_entities ( TDim d1, TDim d2,
                                                     const SortedArray<Id>& d1_entities,
                                                     const SortedArray<Id>& d2_entities,
                                                     Connectivity& d1_d2_connectivity ) const;

            /// \brief Computes incidence 0 -> d, 0 < d <= D
            virtual void vertex_entity_incidence ( TDim d, const SortedArray<Id>& vertices,
                                                   const SortedArray<Id>& d_entities,
                                                   Connectivity& zero_d_connectivity ) const;

            /// \brief Sets the material number of an entity
            /// \pre d == tdim() or d == tdim() - 1
            virtual void set_material_number ( TDim d, Id id, MaterialNumber material );

            /// \brief Returns the material number of an entity
            virtual MaterialNumber get_material_number ( TDim d, Id id ) const;

            //virtual void reset_coordinates(Id vertex_id, std::vector<Coordinate> coords);

            /// \brief write content of database into file
            virtual void save ( std::string filename, const MPI_Comm& comm ) const;

            /// \brief read and setup database from file
            virtual void load ( std::string filename, const MPI_Comm& comm );

            /// \brief creates an EntityPackage containing all entities of given topological dimension and Ids
            virtual void create_entity_package ( int tdim, const std::vector<Id>& ids, EntityPackage* entities ) const;

            /// \brief make a copy of complete database and return pointer to copy
            virtual void copy_from ( const MeshDatabasePtr db );

            /// \brief make a copy of complete database and return pointer to copy.
            /// Compared to copy_from, this function copies the mesh objects stored in mesh_history, rather than copying the corresponding MeshPtrs only.
            virtual void deep_copy_from ( const MeshDatabasePtr db );

          protected:
            /// \brief Creates a new vertex
            Id create_vertex ( const Coordinate* point );

            /// \brief Creates a new entity
            Id create_entity ( TDim entity_dim, const std::vector<Id>& entity_vertices );

            /// Helper function to get connectivtity for a given dimension
            const Connectivity& entity_vertex_connectivity ( TDim entity_dim ) const;
            const Connectivity& vertex_entity_connectivity ( TDim entity_dim ) const;

            /// Maximal topological dimension of entities in database
            /// TDim is essentially part of the type of MeshDatabase
            TDim tdim_;

            /// Geometrical dimension of vertices in database
            /// GDim is essentially part of the type of MeshDatabase
            GDim gdim_;

            /// Vector containing vertex coordinates stored interleaved
            /// (size is gdim_ * number of vertices)
            std::vector<Coordinate> coordinates_;

            /// Array containing d->0 connectivities for 0 < d <= tdim_
            std::vector<Connectivity> entity_vertex_connectivities_;

            /// Array containing 0->d connectivities for 0 < d <= tdim_
            std::vector<VertexConnectivity> vertex_entity_connectivities_;

            /// Vertex search table - used for looking up vertex by coordinates
            ScopedPtr<VertexSearchTable>::Type vertex_search_table_;

            /// Material numbers
            std::vector<Id> material_numbers_[2];
        };

        //////////////// VertexConnectivity ////////////////
        /// Contains connectivity for vertices (0->d). This class is
        /// different from Connectivity in that the connections for a
        /// given vertex can be edited, and that they are stored in sorted
        /// order. This class is only meant to be used by the MeshDatabase class

        class VertexConnectivity : public Connectivity
        {
            // TODO(Staffan): maybe this class can be refactored to use Connectivity instead
            // by moving functions outside
          public:
            typedef std::vector<Id>::const_iterator ConstConnectionIterator;
            void add_vertex ( );
            void add_vertex_connection ( Id vertex_id, Id connected_entity_id );
            bool find_entity ( const std::vector<Id>& entity_vertices,
                               const TDim tdim, const MeshDatabase* db, Id& id ) const;
        };

        //////////////// VertexSearchTable ////////////////
        /// Contains search table used to find vertices from their
        /// coordinates. This class is only meant to be used by the
        /// MeshDatabase class.

        class VertexSearchTable
        {
          public:
            VertexSearchTable ( GDim gdim, const std::vector<Coordinate>& coordinates );
            void add_vertex ( const Coordinate* point, Id id );
            bool find_vertex ( const Coordinate* point, Id& id ) const;
            Coordinate distance ( const Coordinate* point ) const;
            Coordinate distance ( const Coordinate* point1, const Coordinate* point2 ) const;
            Coordinate distance ( Id id ) const;
            std::vector<Id> vertices ( ) const;
            void set_vertices ( std::vector<Id>& vertices );
            const Coordinate* get_coordinates ( Id id ) const;
            void update_vertex ( Id id );
          private:
            bool compare_distances ( Coordinate d1, Coordinate d2 ) const;
            bool compare_points ( const Coordinate* p1, const Coordinate* p2 ) const;

            std::vector<Id> vertices_; // vector of vertex Id:s sorted by distance from origin
            const GDim gdim_;
            const std::vector<Coordinate>& coordinates_; // reference to coordinates_ vector in database
        };

    }
} // namespace hiflow
#endif /* _MESH_DB_H_ */
