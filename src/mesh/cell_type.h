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

#ifndef HIFLOW_MESH_CELL_TYPE_H
#    define HIFLOW_MESH_CELL_TYPE_H

#    include <iosfwd>
#    include <utility>
#    include <vector>

#    include "mesh/types.h"

namespace hiflow
{
    namespace mesh
    {

        class RefinementTree;

        /// \brief Representation of a type of cell, including its
        /// lower-dimensional constituent entities and its possible
        /// refinements.

        class CellType
        {
          public:
            /// \brief Enumeration of available cell types.

            enum Tag
            {
                POINT = 0,
                LINE,
                TRIANGLE,
                QUADRILATERAL,
                TETRAHEDRON,
                HEXAHEDRON,
                PYRAMID,
                NUM_CELL_TYPES
            };

            /// \brief Access to CellType instance corresponding to given tag.
            static const CellType& get_instance ( Tag tag );

            /// \brief Access to CellType instance for entities of given
            /// dimension and number of vertices.
            static const CellType& get_instance ( TDim dim, int num_vertices );

            /// \brief Destructor.
            virtual ~CellType ( );

            /// \brief Access the tag corresponding to the cell type.
            Tag tag ( ) const;

            /// \brief Access the topological dimension of the cell type.
            TDim tdim ( ) const;

            /// \brief Access the number of regular sub-entities of the cell type.
            EntityCount num_regular_entities ( TDim dim ) const;

            /// \brief Access the total number of sub-entities of the cell type.
            EntityCount num_entities ( TDim dim ) const;

            /// \brief Access a refined vertex.
            const std::vector<int>& refined_vertex ( int i ) const;

            /// \brief Access the total number of vertices.
            int num_vertices ( ) const;

            /// \brief Access local vertex numbers of an entity.
            const std::vector<int>& local_vertices_of_entity ( TDim d, int i ) const;

            /// \brief Map local vertices of an entity to vertex id:s.
            std::vector<int> vertices_of_entity ( TDim d, int i, const std::vector<int>& cell_vertices ) const;

            /// \brief Access the entity numbers of a (refined) cell.
            const std::vector<int>& sub_entities_of_cell ( TDim d, int cell ) const;

            /// \brief Access the same-dimension parents of an entity.
            std::vector<int> parents_of_entity ( TDim d, int i ) const;

            /// \brief Access the regular parent of an entity.
            int regular_parent ( TDim d, int i ) const;

            /// \brief Access the refinement tree for a given refinement type.
            const RefinementTree* refinement_tree ( int refinement_type ) const;

            /// \brief Check geometry of an instance of this CellType. If
            /// not overridden in the sub-class, this function simply
            /// returns true, without performing any check.
            template<typename T>
            bool check_cell_geometry ( const std::vector<T>& coords, GDim gdim ) const;

          private:
            ////////////////////////////////////////////////////////////
            /// Callback interface to be implemented in sub-classes ////
            ////////////////////////////////////////////////////////////

            /// \brief Callback in which the regular entities of the
            /// CellType should be created through calls to add_regular_entity.
            virtual void create_regular_entities ( ) = 0;

            /// \brief Callback in which the refined vertices should be
            /// created through calls to add_refined_vertex.
            virtual void create_refined_vertices ( ) = 0;

            /// \brief Callback in which the refined cells should be
            /// created through calls to add_refined_cell.
            virtual void create_refined_cells ( ) = 0;

            /// \brief Callback in which refinements should be created
            /// through calls to add_refinement.
            virtual void create_refinements ( ) = 0;

          protected:
            ////////////////////////////////////////////////////////////
            /// Service interface for sub-classes to use.           ////
            ////////////////////////////////////////////////////////////

            /// \brief Base class constructor. The sub-class should call
            /// this with the correct values for the parameters.
            CellType ( Tag tag, TDim dim, int number_vertices );

            /// \brief Add a regular sub-entity.
            int add_regular_entity ( TDim dim, const std::vector<int>& vertices );

            /// \brief Add a refined vertex.
            int add_refined_vertex ( const std::vector<int>& super_vertices );

            /// \brief Add a refined cell.
            int add_refined_cell ( const std::vector<int>& vertices );

            /// \brief Add a refinement.
            int add_refinement ( const std::vector<int>& sub_cell_numbers );

          private:
            typedef std::vector< std::vector<int> > ConnectionList;
            friend std::ostream& operator<< ( std::ostream&, const CellType& );

            /// \brief Initializes all CellTypes.
            static void initialize_cell_types ( );

            /// \brief Returns true if CellTypes were already initialized.
            static bool is_initialized ( );

            /// \brief Add an entity.
            int add_entity ( TDim dim, const std::vector<int>& vertices );

            /// \brief Checks whether a vertex is regular.
            bool is_regular_vertex ( int v ) const;

            /// \brief Checks whether a sub-entity is regular
            bool is_regular_entity ( TDim dim, int e ) const;

            /// \brief Checks whether all vertices exist.
            bool vertices_exist ( const std::vector<int>& vertices ) const;

            /// \brief Searches for a sub-entity based on its vertices.
            bool find_entity ( TDim dim, const std::vector<int>& vertices,
                               int* entity_number ) const;

            /// \brief Checks if one entity is contained in another
            /// entity of the same dimension.
            bool is_entity_contained ( TDim dim, int entity, int contained_entity ) const;

            /// \brief Computes the hierarchies of all sub-entities, once
            /// they have been initialized in the subclasses.
            void compute_sub_entity_hierarchy ( );

            /// \brief Creates all sub-entities of a refined cell.
            ///
            /// \pre The regular sub-entities of the CellType for the
            /// sub-cell have already been created.
            void create_sub_entities_of_refined_cell ( int cell );

            /// \brief Container instances of all CellTypes.
            ///
            /// \details A SharedPtr is used since a ScopedPtr has no copy constructor
            static std::vector< SharedPtr<CellType>::Type > cell_types_;

            /// \brief Access sub_entities_ for a given dimension.
            ConnectionList& entities ( TDim dim );
            const ConnectionList& entities ( TDim dim ) const;

            /// \brief Access to sub_cell_sub_entities_ for a given dimension.
            ConnectionList& cell_sub_entities ( TDim dim );
            const ConnectionList& cell_sub_entities ( TDim dim ) const;

            /// \brief Access sub_entity_parents_ for a given dimension.
            ConnectionList& entity_parents ( TDim dim );
            const ConnectionList& entity_parents ( TDim dim ) const;

            /// \brief Topological dimension of the CellType.
            const TDim tdim_;

            /// \brief Tag of the CellType.
            const Tag tag_;

            /// \brief Description of regular and refined vertices.
            ConnectionList vertices_;

            /// \brief Vertices of all entities (regular and refined) 0 < d <= D.
            std::vector<ConnectionList> entities_;

            /// \brief Cell -> entity connectivity 0 < d < D
            std::vector<ConnectionList> cell_sub_entities_;

            /// \brief Same-dimension parents of all entities 0 < d <= D.
            std::vector<ConnectionList> entity_parents_;

            /// \brief List of RefinementTree:s for the CellType.
            std::vector< SharedPtr<RefinementTree>::Type > refinement_trees_;
        };

        std::ostream& operator<< ( std::ostream& os, const CellType& cell_type );

        ////////////////////////////////////////////////////
        //////////////// Concrete CellTypes ////////////////
        ////////////////////////////////////////////////////

        //////////////// Point ////////////////

        class Point : public CellType
        {
            friend class CellType;
          private:
            Point ( );

            virtual void create_regular_entities ( );
            virtual void create_refined_vertices ( );
            virtual void create_refined_cells ( );
            virtual void create_refinements ( );
        };

        //////////////// Line ////////////////

        class Line : public CellType
        {
            friend class CellType;
          private:
            Line ( );

            virtual void create_regular_entities ( );
            virtual void create_refined_vertices ( );
            virtual void create_refined_cells ( );
            virtual void create_refinements ( );
        };

        //////////////// Triangle ////////////////

        class Triangle : public CellType
        {
          public:

            enum
            {
                TRI_4TRI_REFINEMENT_TYPE = 0,
                TRI_2TRI_V0_REFINEMENT_TYPE,
                TRI_2TRI_V1_REFINEMENT_TYPE,
                TRI_2TRI_V2_REFINEMENT_TYPE
            };

            friend class CellType;

            virtual bool check_cell_geometry ( const std::vector<Coordinate>& coords, GDim gdim ) const;
          private:
            Triangle ( );

            virtual void create_regular_entities ( );
            virtual void create_refined_vertices ( );
            virtual void create_refined_cells ( );
            virtual void create_refinements ( );
        };

        //////////////// Quadrilateral ////////////////

        class Quadrilateral : public CellType
        {

            enum
            {
                QUAD_4QUAD_REFINEMENT_TYPE = 0,
                QUAD_2QUAD_HORIZONTAL_REFINEMENT_TYPE,
                QUAD_2QUAD_VERTICAL_REFINEMENT_TYPE,
                QUAD_2TRI_V0_V2_DIAGONAL_REFINEMENT_TYPE,
                QUAD_2TRI_V1_V3_DIAGONAL_REFINEMENT_TYPE,
                QUAD_3TRI_V0_V1_CORNERS_REFINEMENT_TYPE,
                QUAD_3TRI_V1_V2_CORNERS_REFINEMENT_TYPE,
                QUAD_3TRI_V2_V3_CORNERS_REFINEMENT_TYPE,
                QUAD_3TRI_V3_V0_CORNERS_REFINEMENT_TYPE
            };

            friend class CellType;
          private:
            Quadrilateral ( );
            virtual void create_regular_entities ( );
            virtual void create_refined_vertices ( );
            virtual void create_refined_cells ( );
            virtual void create_refinements ( );
        };

        //////////////// Tetrahedron ////////////////

        class Tetrahedron : public CellType
        {
          public:

            enum
            {
                TET_8TET_REFINEMENT_TYPE = 0,
                TET_4HEX_REFINEMENT_TYPE
            };
            friend class CellType;

            virtual bool check_cell_geometry ( const std::vector<Coordinate>& coords, GDim gdim ) const;
          private:
            Tetrahedron ( );
            virtual void create_regular_entities ( );
            virtual void create_refined_vertices ( );
            virtual void create_refined_cells ( );
            virtual void create_refinements ( );
        };

        //////////////// Hexahedron ////////////////

        class Hexahedron : public CellType
        {
          public:

            enum
            {
                HEX_8HEX_REFINEMENT_TYPE
            };
            friend class CellType;
          private:
            Hexahedron ( );

            virtual void create_regular_entities ( );
            virtual void create_refined_vertices ( );
            virtual void create_refined_cells ( );
            virtual void create_refinements ( );
        };

        //////////////// Pyramid ////////////////

        class Pyramid : public CellType
        {
          public:

            enum
            {
                PYR_6PYR4TET_REFINEMENT_TYPE
            };
            friend class CellType;
          private:
            Pyramid ( );

            virtual void create_regular_entities ( );
            virtual void create_refined_vertices ( );
            virtual void create_refined_cells ( );
            virtual void create_refinements ( );
        };

        //////////////// Helper functions ////////////////

        /// \brief Compares two vectors as if they were sets, i.e. without
        /// respect to order.
        bool compare_as_sets ( const std::vector<int>& v1, const std::vector<int>& v2 );

    } // namespace mesh
} // namespace hiflow

#endif
