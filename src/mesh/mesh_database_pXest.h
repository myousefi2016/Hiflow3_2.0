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

/// \author Jonathan Schwegler, Philipp Gerstner

#ifndef MESH_DATABASE_PXEST_H
#    define MESH_DATABASE_PXEST_H

#    include <vector>
#    include <map>
#    include <iostream>
#    include <string>
#    include "common/sorted_array.h"
#    include "mesh_database.h"
#    include "config.h"
#    include "pXest_utils.h"
#    include "pXest_structs.h"

#    ifdef WITH_P4EST
#        include "p4est.h"
#        include "p8est.h"
#        include "p4est_connectivity.h"
#        include "p8est_connectivity.h"
#        include "p4est_ghost.h"
#        include "p8est_ghost.h"
#    endif

namespace hiflow
{
    namespace mesh
    {

        class MeshPXestDatabase;
        typedef SharedPtr< MeshPXestDatabase >::Type MeshPXestDatabasePtr;

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

        class MeshPXestDatabase : public MeshDatabase
        {
          public:
            /// \brief Constructor
            MeshPXestDatabase ( TDim tdim, GDim gdim );
            ~MeshPXestDatabase ( );

            virtual void build ( TDim dim, const SortedArray<Id>& cells, SortedArray<Id>& d_entities, Connectivity& cell_d_connections );

            // STRUCT_TODO move additional functions to protected
#    ifdef WITH_P4EST
            /// \brief set p4est connectivity for coarse mesh

            virtual void set_p4est_conn ( p4est_connectivity_t* conn );

            virtual void set_p8est_conn ( p8est_connectivity_t* conn );

            /// \brief set p4est main object: local part of distributed forest

            virtual inline void set_p4est_forest ( p4est_t* forest )
            {
                p4est_forest_ = forest;
            }

            virtual inline void set_p8est_forest ( p8est_t* forest )
            {
                p8est_forest_ = forest;
            }

            /// \brief set p4est ghost cells of forest

            virtual inline void set_p4est_ghost ( p4est_ghost_t* ghost )
            {
                p4est_ghost_ = ghost;
            }

            virtual inline void set_p8est_ghost ( p8est_ghost_t* ghost )
            {
                p8est_ghost_ = ghost;
            }

            virtual inline p4est_connectivity_t* get_p4est_conn ( ) const
            {
                return this->p4est_conn_;
            }

            virtual inline p8est_connectivity_t* get_p8est_conn ( ) const
            {
                return this->p8est_conn_;
            }

            virtual inline p4est_t* get_p4est_forest ( ) const
            {
                return this->p4est_forest_;
            }

            virtual inline p8est_t* get_p8est_forest ( ) const
            {
                return this->p8est_forest_;
            }

            virtual inline p4est_ghost_t* get_p4est_ghost ( ) const
            {
                return this->p4est_ghost_;
            }

            virtual inline p8est_ghost_t* get_p8est_ghost ( ) const
            {
                return this->p8est_ghost_;
            }

            virtual inline void set_connectivity_type ( p4est_connect_type_t ctype )
            {
                this->ctype4_ = ctype;
            }

            virtual inline void set_connectivity_type ( p8est_connect_type_t ctype )
            {
                this->ctype8_ = ctype;
            }
#    endif

            virtual inline void set_quadrants_sorted ( bool flag )
            {
                this->quadrant_arrays_sorted_ = flag;
            }

            virtual inline bool quadrants_sorted ( ) const
            {
                return this->quadrant_arrays_sorted_;
            }

            /// \brief Add new entity_to_quad_map
            /// @param[in] tdim topological dimension of entities
            /// @param[in] entity_id Id of entity
            /// @param[in] coord coordinates of p4est quadrant corresponding to entity
            /// @param[in] which_mapping flag indicating whether map should be added to coarse maps (0) or distributed maps(1)
            virtual void add_entity_to_quad_coord ( int tdim, Id entity_id, QuadCoord& coord, int which_mapping );

            virtual void permute_entity_to_quad_map ( const std::vector<size_t>& perm, int tdim, int which_mapping );

            virtual void clear_entity_to_quad_map ( int tdim, int which_mapping );

            /// \brief Set complete entity_to_quad_mapping for specified topological dimension
            /// @param[in] tdim topological dimension
            /// @param[in] which_mapping flag indicating whether map should be see as coarse maps (0) or distributed maps(1)
            /// @param[in] maps vector of entity_to_quad maps to be set
            virtual void set_entity_to_quad_map ( int tdim, int which_mapping, EntityToQuadMap& maps );

            virtual QuadCoord get_entity_to_quad_coord ( int tdim, Id entity_id, int which_mapping );

            virtual EntityToQuadMap* get_entity_to_quad_map ( int tdim, int which_mapping );

            /// \brief set data structures used for building p4est connectivity for coarse mesh
            virtual void set_conn_data ( int num_vertices, int num_cells, std::vector<double>& vertices, std::vector<int>& tree_to_vertices );

            virtual inline std::vector<double>* get_conn_vertices ( )
            {
                return &this->conn_vertices_;
            }

            virtual inline std::vector<int>* get_tree_to_vertices ( )
            {
                return &this->tree_to_vertices_;
            }

            virtual inline int get_num_conn_vertices ( ) const
            {
                return this->num_conn_vertices_;
            }

            virtual inline int get_num_conn_cells ( ) const
            {
                return this->num_conn_cells_;
            }

            virtual inline int get_layer_width ( )
            {
                return this->layer_width_;
            }

            virtual inline void set_layer_width ( int width )
            {
                this->layer_width_ = width;
            }

            /// \brief get meshpointer for given history index
            /// @param[in] history_index
            /// @return Pointer to mesh
            // RESTRUCTURING -> Database
            // STRUCT_TODO make const
            virtual MeshPtr get_mesh ( int history_index );

            virtual std::map< int, MeshPtr> get_mesh_history ( ) const
            {
                return mesh_history_;
            }

            /// \brief set mesh pointer with given history index
            // RESTRUCTURING -> Database
            // STRUCT_TODO const pointer?
            virtual void set_mesh ( MeshPtr mesh_ptr, int history_index );

            /// \brief add mesh history indices for given cells
            // RESTRUCTURING -> Database
            virtual void add_cell_to_mesh_maps ( const std::vector<Id>& cells, std::vector<int>& history_indices );

            /// \brief return all mesh history indices for a given cell
            // RESTRUCTURING -> Database
            // STRUCT_TODO make const
            virtual int get_last_mesh_history_index ( Id cell_id );

            /// \brief write content of database into file
            virtual void save ( std::string filename, const MPI_Comm& comm ) const;

            /// \brief read and setup database from file
            virtual void load ( std::string filename, const MPI_Comm& comm );

            /// \brief make a copy of complete database and return pointer to copy
            virtual void copy_from ( const MeshPXestDatabasePtr db );

            /// \brief make a copy of complete database and return pointer to copy.
            /// Compared to copy_from, this function copies the mesh objects stored in mesh_history, rather than copying the corresponding MeshPtrs only.
            virtual void deep_copy_from ( const MeshPXestDatabasePtr db );

          protected:

#    ifdef WITH_P4EST
            // Connectivity object for coarse mesh
            p4est_connectivity_t* p4est_conn_;
            p8est_connectivity_t* p8est_conn_;

            // Forest object for grid hierarchy
            p4est_t* p4est_forest_;
            p8est_t* p8est_forest_;

            // Ghost layer for forest
            p4est_ghost_t* p4est_ghost_;
            p8est_ghost_t* p8est_ghost_;

            p4est_connect_type_t ctype4_;
            p8est_connect_type_t ctype8_;
#    endif

            // data for building p4est connectivity
            std::vector<double> conn_vertices_;
            std::vector<int> tree_to_vertices_;
            int num_conn_vertices_;
            int num_conn_cells_;

            // Mapping from hiflow entity to quadrant in forest for distributed mesh
            EntityToQuadMapping entity_to_quad_;

            // Mapping from hiflow entity to quadrant in forest for nondistributed coarse mesh
            EntityToQuadMapping coarse_entity_to_quad_;

            // RESTRUCTURING -> Database
            std::map< Id, int > cell_to_mesh_mapping_;
            std::map< int, MeshPtr> mesh_history_;

            // Width of ghost layer in number of cells
            int layer_width_;

            // flag indicating whether quadrant arrays in forest are sorted
            bool quadrant_arrays_sorted_;
        };

    }
} // namespace hiflow
#endif /* _MESH_DB_H_ */
