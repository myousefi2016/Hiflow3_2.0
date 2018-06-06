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

#include "cell_type.h"

#include "common/log.h"
#include "refinement.h"

#include <algorithm>
#include <iostream>
#include <stack>

namespace hiflow
{
    namespace mesh
    {

        //////////////// PUBLIC INTERFACE ////////////////

        /// \param tag  tag for cell type.
        /// \returns reference to the unique CellType object corresponding to tag.

        const CellType& CellType::get_instance ( CellType::Tag tag )
        {
            if ( !is_initialized ( ) )
            {
                initialize_cell_types ( );
            }
            return *( cell_types_[tag].get ( ) );
        }

        /// \param dim              topological dimension of the cell type.
        /// \param num_vertices     number of (regular) vertices in the cell type.
        /// \returns reference to the unique CellType object corresponding to (dim, num_vertices)

        const CellType& CellType::get_instance ( TDim dim, int num_vertices )
        {
            Tag tag = CellType::NUM_CELL_TYPES;
            switch ( dim )
            {
                case 0:
                    tag = CellType::POINT;
                    break;
                case 1:
                    tag = CellType::LINE;
                    break;
                case 2:
                    switch ( num_vertices )
                    {
                        case 3:
                            tag = CellType::TRIANGLE;
                            break;
                        case 4:
                            tag = CellType::QUADRILATERAL;
                            break;
                        default:
                            assert ( false && "No 2d cell type found" );
                    }
                    break;
                case 3:
                    switch ( num_vertices )
                    {
                        case 4:
                            tag = CellType::TETRAHEDRON;
                            break;
                        case 5:
                            tag = CellType::PYRAMID;
                            break;
                        case 8:
                            tag = CellType::HEXAHEDRON;
                            break;
                        default:
                            assert ( false && "No 3d cell type found" );
                    }
                    break;
                default:
                    assert ( false && "No cell type found" );
            }

            assert ( tag != CellType::NUM_CELL_TYPES );
            return get_instance ( tag );
        }

        CellType::~CellType ( )
        {
            // Empty destructor - only here since CellType has virtual functions.
        }

        /// \returns the tag corresponding to the cell type.

        CellType::Tag CellType::tag ( ) const
        {
            return tag_;
        }

        /// \returns the topological dimension of the cell type.

        TDim CellType::tdim ( ) const
        {
            return tdim_;
        }

        /// \param dim   toplogical dimension of the the sub-entities.
        /// \pre  0 <= dim <= CellType::tdim() .
        /// \returns  the number of regular sub-entities of the cell type.

        EntityCount CellType::num_regular_entities ( TDim dim ) const
        {
            // Regular entities are those which are direct sub-entities of
            // the regular cell, which is always the first cell (number 0).
            // Here we check how many vertices / sub-entities this
            // cell has.

            // There is only one regular cell. This check is done to avoid
            // problems for 0d and 1d CellTypes.
            if ( dim == tdim_ )
            {
                return 1;
            }
            else
            {
                if ( dim == 0 )
                {
                    // number of vertices is length of entities(tdim)-vector
                    return entities ( tdim_ )[0].size ( );
                }
                else
                {
                    // number of sub-entities is length of cell_sub_entities-vector for root-cell
                    return cell_sub_entities ( dim )[0].size ( );
                }
            }
        }

        /// \param dim   toplogical dimension of the the sub-entities.
        /// \pre  0 <= dim <= CellType::tdim() .
        /// \returns  the total number of sub-entities of dimension dim for all refined cells.

        EntityCount CellType::num_entities ( TDim dim ) const
        {
            if ( dim == 0 )
            {
                return vertices_.size ( );
            }
            else
            {
                return entities ( dim ).size ( );
            }
        }

        /// \param  i   the local number of the vertex.
        /// \returns    vector of parent vertex number of vertex i.

        const std::vector<int>& CellType::refined_vertex ( int i ) const
        {
            range_check ( vertices_, i );
            return vertices_[i];
        }

        /// \returns  total number of vertices, including both regular and refined vertices.

        int CellType::num_vertices ( ) const
        {
            return vertices_.size ( );
        }

        /// \param d  dimension of the entity.
        /// \param i  local entity number of the entity.
        /// \pre 0 <= d <= CellType::tdim() .
        /// \pre 0 <= i < CellType::num_entities(d) .
        /// \returns  a vector containing the local vertex numbers of the entity.

        const std::vector<int>& CellType::local_vertices_of_entity ( TDim d, int i ) const
        {
            assert ( d >= 0 );
            assert ( d <= tdim_ );
            assert ( i >= 0 );
            assert ( i < num_entities ( d ) );
            if ( d == 0 )
            {
                return vertices_[i];
            }
            else
            {
                return entities ( d )[i];
            }
        }

        /// \details this function maps provides a way to "construct" the
        /// id:s of a local sub-entity of a cell, given the id:s of the
        /// cell's vertices. For instance, if the second (i = 1)
        /// one-dimensional sub-entity of a two-dimensional CellType is
        /// (1, 2), the global vertices (101, 102) of that sub-entity in a
        /// cell (100, 101, 102, 103) can be the found by a call to
        /// vertices_of_entity(1, 1, vector<int>({100, 101, 102, 103})) .
        ///
        /// \details The function can be used either for sub-entities of
        /// the regular cell, in which only vertex id:s for the regular
        /// vertices must be provided; or for any refined cell, in which
        /// the vertex id:s of all vertices must be provided.
        ///
        /// \param d              dimension of the entity.
        /// \param i              local entity number of the entity.
        /// \param cell_vertices  vertex id:s for the vertices of the entity.
        /// \pre 0 <= d <= CellType::tdim() .
        /// \pre 0 <= i < CellType::num_entities(d) .
        /// \pre 0 <= cell_vertices[j] < CellType::num_regular_entities(0) if i < CellType::num_regular_entites(i) .
        /// \pre 0 <= cell_vertices[j] < CellType::num_vertices if i >= CellType::num_regular_entites(i) .

        std::vector<int> CellType::vertices_of_entity ( TDim d, int i,
                                                        const std::vector<int>& cell_vertices ) const
        {
            // NB: some code will call this with the cell_vertices only
            // containing entries for the regular vertices, and other code
            // with call it with all vertices. Both should work, as long
            // as those vertex numbers which are contained in the entity,
            // are smaller than the size of cell_vertices. This means in
            // practice, that if only regular vertices are provided, then
            // only regular entities (that are part of the root cell) can
            // be queried.

            assert ( d <= tdim_ );

            const int num_vertices = ( d == 0 ) ? 1 : entities ( d )[i].size ( );

            std::vector<int> entity_vertices ( num_vertices, -1 );

            if ( d == 0 )
            {
                // if dimension is 0, we simply return the i:th vertex in the vector
                entity_vertices[0] = cell_vertices[i];
            }
            else
            {
                // Map cell_vertices -> entity_vertices via the entities(d)[i] array.
                const int* p = vec2ptr ( entities ( d )[i] );
                for ( int k = 0; k < num_vertices; ++k, ++p )
                {
                    entity_vertices[k] = cell_vertices[*p];
                }
            }
            return entity_vertices;
        }

        /// \param d     topological dimension of the sub-entity.
        /// \param cell  local number of the refined cell (0 for the regular cell).
        /// \pre 0 < d < CellType::tdim()
        /// \returns     reference to vector containing the numbers of the sub-entities of dimension d of the cell, in the current CellType.

        const std::vector<int>& CellType::sub_entities_of_cell ( TDim d, int cell ) const
        {
            return cell_sub_entities ( d )[cell];
        }

        /// \details A sub-entity can have one or several parents. One
        /// sub-entity (dim, p) is the parent of another sub-entity (dim,
        /// c) if for each vertex v of (dim, c), (dim, p) contains v or,
        /// if v is a refined vertex, all of v:s parents, recursively.
        ///
        /// \param d   topological dimension of the sub-entity.
        /// \param i   local number of the sub-entity in the CellType.
        /// \pre   0 < d <= CellType::tdim().
        /// \pre   0 <= i < CellType::num_entities(d) .
        ///
        /// \returns vector with the local numbers of all parents of the
        /// entity in this CellType. If the entity has no parent, an empty
        /// vector is returned.

        std::vector<int> CellType::parents_of_entity ( TDim d, int i ) const
        {
            return entity_parents ( d )[i];
        }

        /// \details A regular entity is an entity that is part of the
        /// root cell, and hence only contains regular vertices. Although
        /// a sub-entity of a refined cell can have several parent
        /// entities, only one of them will be a regular entity. Regular
        /// entities have no parents themselves.
        ///
        /// \param d   topological dimension of the sub-entity.
        /// \param i   local number of the sub-entity in the CellType.
        /// \pre   0 < d <= CellType::tdim() .
        /// \pre   0 <= i < CellType::num_entities(d) .
        ///
        /// \returns the unique regular parent entity of the entity, or -1 if no parent exists.

        int CellType::regular_parent ( TDim d, int i ) const
        {
            const std::vector<int> parents = parents_of_entity ( d, i );

            // return -1 if no parents exist
            if ( parents.empty ( ) )
            {
                return -1;
            }

            // find parent entity in the range [0, num_regular_entities(d))
            const int max_regular_entity = num_regular_entities ( d );

            // check that there is one and only one regular parent
            assert ( std::count_if ( parents.begin ( ), parents.end ( ),
                                     std::bind2nd ( std::less<int>( ), max_regular_entity ) ) == 1 );

            // find the regular parent
            std::vector<int>::const_iterator parent_it =
                    std::find_if ( parents.begin ( ), parents.end ( ),
                                   std::bind2nd ( std::less<int>( ), max_regular_entity ) );
            return *parent_it;

        }

        /// \param refinement_type the type of refinement.
        ///
        /// \returns pointer to a refinement tree object for the
        /// refinement. This object should not be deleted by the caller.

        const RefinementTree* CellType::refinement_tree ( int refinement_type ) const
        {
            assert ( refinement_type >= 0 );
            assert ( refinement_type < static_cast < int > ( refinement_trees_.size ( ) ) );
            return refinement_trees_[refinement_type].get ( );
        }

        /// \details This function can be overridden in the various
        /// CellType sub-classes to provide the possibility to verify if
        /// the geometry of a cell of the CellType is correct. If it is
        /// not overridden, true is returned by default.
        ///
        /// \param coords   the coordinates of the cell in interleaved form.
        /// \param gdim     the geometrical dimension of the cell (number of coordinates per point in the coords vector).
        /// \pre   coords.size() / gdim = CellType::num_regular_entities(0).
        /// \returns  true if the geometry of the cell is correct.

        template<typename T>
        bool CellType::check_cell_geometry ( const std::vector<T>& coords, GDim gdim ) const
        {
            return true;
        }

        //////////////// END PUBLIC INTERFACE ////////////////

        //////////////// PROTECTED INTERFACE ////////////////

        /// \details The constructor initializes the CellType with
        /// number_vertices regular vertices, which do not have any parents.
        ///
        /// \param tag                 tag for CellType of subclass
        /// \param tdim                topological dimension of subclass
        /// \param number_vertices     number of vertices of subclass

        CellType::CellType ( Tag tag, TDim dim, int number_vertices )
        : tdim_ ( dim ),
        tag_ ( tag )
        {
            // create all regular vertices
            vertices_.resize ( number_vertices );
            for ( int v = 0; v < number_vertices; ++v )
            {
                vertices_[v] = std::vector<int>( 1, v );
            }

            // initialize remaining structures
            if ( dim > 0 )
            {
                entities_.resize ( dim );
                cell_sub_entities_.resize ( dim - 1 );
                entity_parents_.resize ( dim );

                // create vector of vertices of the regular cell
                std::vector<int> cell_vertices ( number_vertices, -1 );
                for ( int v = 0; v < number_vertices; ++v )
                {
                    cell_vertices[v] = v;
                }

                // add to entities_ structure, as cell 0
                entities ( tdim_ ).push_back ( cell_vertices );
            }
        }

        /// \details Adds a regular entity with the given topological
        /// dimension and vertices. It is verified that the vertices
        /// corresponding to the indices given in the vertices vector
        /// are all regular.
        ///
        /// \param tdim      topological dimension of the entity
        /// \param vertices  local vertex numbers defining the entity
        /// \returns         the local entity number

        int CellType::add_regular_entity ( TDim dim, const std::vector<int>& vertices )
        {
            assert ( dim > 0 );
            assert ( dim < tdim_ );
            assert ( vertices_exist ( vertices ) );

            // check that all vertices are regular
            for ( std::vector<int>::const_iterator it = vertices.begin ( );
                  it != vertices.end ( ); ++it )
            {
                assert ( is_regular_vertex ( *it ) );
            }

            // check that the same entity is not added twice
            assert ( !find_entity ( dim, vertices, 0 ) );

            const int entity_number = add_entity ( dim, vertices );
            cell_sub_entities ( dim )[0].push_back ( entity_number );
            return entity_number;
        }

        /// \param super_vertices   the vertices from which this vertex is computed.
        /// \returns                the number of the refined vertex

        int CellType::add_refined_vertex ( const std::vector<int>& super_vertices )
        {
            assert ( vertices_exist ( super_vertices ) );
            const int number = vertices_.size ( );
            vertices_.push_back ( super_vertices );
            return number;
        }

        /// \param vertices         the vertices defining the refined cell.
        /// \returns                the number of the new cell

        int CellType::add_refined_cell ( const std::vector<int>& vertices )
        {
            return add_entity ( tdim_, vertices );
        }

        /// \param sub_cell_numbers   sub-cells to be created in the refinement.
        /// \returns                  the number of the new refinement

        int CellType::add_refinement ( const std::vector<int>& sub_cell_numbers )
        {
            // Check that sub-cells exist, and that the regular cell
            // (number 0) is not among them.
            assert ( *std::min_element ( sub_cell_numbers.begin ( ), sub_cell_numbers.end ( ) ) > 0 );
            assert ( *std::max_element ( sub_cell_numbers.begin ( ), sub_cell_numbers.end ( ) ) <
                     num_entities ( tdim ( ) ) );

            // Create new tree object.
            SharedPtr<RefinementTree>::Type tree =
                    SharedPtr<RefinementTree>::Type ( new RefinementTree ( tdim ( ), *this ) );
            assert ( tree != 0 );

            // Create new sub-tree of root of the tree with the given
            // sub-cell numbers.
            tree->add_subtree ( 0, sub_cell_numbers );

            // Add tree to list of refinements, and return its number.
            const int refinement_number = refinement_trees_.size ( );
            refinement_trees_.push_back ( tree );
            return refinement_number;
        }

        //////////////// END PROTECTED INTERFACE ////////////////

        //////////////// PRIVATE FUNCTIONS ////////////////

        void CellType::initialize_cell_types ( )
        {
            if ( !is_initialized ( ) )
            {
                cell_types_.resize ( CellType::NUM_CELL_TYPES );

                // Create one instance for each CellType::Tag.
                for ( int t = 0; t < CellType::NUM_CELL_TYPES; ++t )
                {
                    CellType* ptr = 0;

                    switch ( t )
                    {
                        case POINT:
                            ptr = new Point ( );
                            break;
                        case LINE:
                            ptr = new Line ( );
                            break;
                        case TRIANGLE:
                            ptr = new Triangle ( );
                            break;
                        case QUADRILATERAL:
                            ptr = new Quadrilateral ( );
                            break;
                        case TETRAHEDRON:
                            ptr = new Tetrahedron ( );
                            break;
                        case HEXAHEDRON:
                            ptr = new Hexahedron ( );
                            break;
                        case PYRAMID:
                            ptr = new Pyramid ( );
                            break;
                        default:
                            assert ( false );
                    };

                    assert ( ptr != 0 );

                    assert ( ptr->tag ( ) == t );

                    cell_types_[t] = SharedPtr<CellType>::Type ( ptr );

                    // Add a cell->sub_entity connectivity entry for the regular cell
                    // for each entity dimension.
                    for ( int d = 1; d < ptr->tdim ( ); ++d )
                    {
                        ptr->cell_sub_entities ( d ).resize ( 1 );
                    }
                }

                // Initialize regular entities for all cell types.
                for ( int t = 0; t < CellType::NUM_CELL_TYPES; ++t )
                {
                    CellType* ptr = cell_types_[t].get ( );
                    ptr->create_regular_entities ( );
                }

                // Initialize refined vertices and refined cells. This
                // must be done after all the regular entities of all
                // CellTypes have been created, since the sub-entities of
                // the refined cells can be of any CellType.
                for ( int t = 0; t < CellType::NUM_CELL_TYPES; ++t )
                {
                    CellType* ptr = cell_types_[t].get ( );

                    if ( ptr->tdim ( ) == 0 )
                    {
                        continue;
                    }

                    ptr->create_refined_vertices ( );
                    ptr->create_refined_cells ( );

                    ptr->compute_sub_entity_hierarchy ( );

                    ptr->create_refinements ( );
                }
            }
        }

        bool CellType::is_initialized ( )
        {
            return cell_types_.size ( ) == CellType::NUM_CELL_TYPES;
        }

        int CellType::add_entity ( TDim dim, const std::vector<int>& vertices )
        {
            assert ( dim > 0 );
            assert ( dim <= tdim_ );
            assert ( vertices_exist ( vertices ) );

            const int number = entities ( dim ).size ( );
            entities ( dim ).push_back ( vertices );
            return number;
        }

        bool CellType::is_regular_vertex ( int v ) const
        {
            assert ( v >= 0 );
            assert ( v < static_cast < int > ( vertices_.size ( ) ) );

            // A vertex is regular if it is connected only to itself.
            return vertices_[v].size ( ) == 1 && vertices_[v][0] == v;
        }

        bool CellType::is_regular_entity ( TDim dim, int e ) const
        {
            assert ( dim > 0 );
            assert ( dim <= tdim_ );

            const std::vector<int>& vertices = entities ( dim )[e];

            // An entity is regular if all its vertices are regular.
            for ( std::vector<int>::const_iterator it = vertices.begin ( );
                  it != vertices.end ( ); ++it )
            {
                if ( !is_regular_vertex ( *it ) )
                {
                    return false;
                }
            }
            return true;
        }

        bool CellType::vertices_exist ( const std::vector<int>& vertices ) const
        {
            const int num_vertices = vertices_.size ( );

            // Check that all vertices exist.
            for ( std::vector<int>::const_iterator it = vertices.begin ( );
                  it != vertices.end ( ); ++it )
            {
                if ( ( *it < 0 ) || ( *it >= num_vertices ) )
                {
                    return false;
                }
            }
            return true;
        }

        /// \param dim       topological dimension of the sought entity.
        /// \param vertices   list of vertex numbers defining the sought entity.
        /// \param entity_number  output parameter for the number of the entity, if it is found.
        /// \returns   true if the sought entity was found.

        bool CellType::find_entity ( TDim dim, const std::vector<int>& vertices,
                                     int* entity_number ) const
        {
            assert ( dim > 0 );
            assert ( dim <= tdim_ );

            const ConnectionList& existing_entities = entities ( dim );

            int k = 0;
            for ( ConnectionList::const_iterator it = existing_entities.begin ( );
                  it != existing_entities.end ( ); ++it )
            {
                if ( compare_as_sets ( vertices, *it ) )
                {
                    if ( entity_number != 0 )
                    {
                        *entity_number = k;
                    }
                    return true;
                }
                ++k;
            }

            return false;
        }

        /// \param dim               topological dimension of the entities.
        /// \param entity             entity number for the containing entity.
        /// \param contained_entity   entity number for the contained entity.
        /// \returns    true if all vertices of contained_entity are either equal to,
        /// or children of, the vertices in entity.

        bool CellType::is_entity_contained ( TDim dim, int entity, int contained_entity ) const
        {

            // vertex numbers of containing entity
            const std::vector<int>& entity_vertices = entities ( dim )[entity];

            // vertex numbers of the contained entity
            const std::vector<int>& contained_entity_vertices = entities ( dim )[contained_entity];

            // check first if entities are equal - this does not count as containment
            if ( compare_as_sets ( entity_vertices, contained_entity_vertices ) )
            {
                return false;
            }

            // stack for vertex numbers which remain to be checked
            std::stack<int> vertex_numbers;

            // put vertex numbers of the contained entity on the stack
            for ( std::vector<int>::const_iterator it = contained_entity_vertices.begin ( );
                  it != contained_entity_vertices.end ( ); ++it )
            {
                vertex_numbers.push ( *it );
            }

            bool is_contained = true;

            while ( !vertex_numbers.empty ( ) )
            {
                const int v = vertex_numbers.top ( );
                vertex_numbers.pop ( );

                if ( find ( entity_vertices.begin ( ), entity_vertices.end ( ), v )
                     == entity_vertices.end ( ) )
                {
                    // vertex number not found in containing entity

                    if ( !is_regular_vertex ( v ) )
                    {
                        // add parent vertices of v to stack, in order to
                        // check them in the next iterations
                        const std::vector<int>& parent_vertices = vertices_[v];
                        for ( std::vector<int>::const_iterator it = parent_vertices.begin ( );
                              it != parent_vertices.end ( ); ++it )
                        {
                            vertex_numbers.push ( *it );
                        }
                    }
                    else
                    {
                        // vertex is regular, so contained_entity is not a
                        // subset of entity
                        is_contained = false;
                        break;
                    }
                }
            }

            return is_contained;
        }

        void CellType::compute_sub_entity_hierarchy ( )
        {
            // Check that first cell is regular.
            assert ( is_regular_entity ( tdim_, 0 ) );

            const int num_cells = entities ( tdim_ ).size ( );
            for ( int cell = 1; cell < num_cells; ++cell )
            {
                // Create all sub-entities of all dimensions (can be of different types).
                create_sub_entities_of_refined_cell ( cell );
            }

            // For all entities of all dimensions, create parent relation
            for ( int d = 1; d <= tdim_; ++d )
            {
                const int num_entities = this->num_entities ( d );
                entity_parents ( d ).resize ( num_entities );
                for ( int e = 0; e < num_entities; ++e )
                {
                    for ( int e2 = 0; e2 < num_entities; ++e2 )
                    {
                        if ( is_entity_contained ( d, e, e2 ) )
                        {
                            entity_parents ( d )[e2].push_back ( e );
                        }
                    }
                }
            }
        }

        void CellType::create_sub_entities_of_refined_cell ( int cell )
        {
            range_check ( entities ( tdim_ ), cell );

            // TODO Basically the same algorithm as build() in MeshDatabase -> can we merge?

            // This function should not be called for root cell.
            assert ( cell > 0 );

            std::vector<int>& cell_vertices = entities ( tdim_ )[cell];

            // Get cell type of sub-cell
            const CellType& cell_type = CellType::get_instance ( tdim_, cell_vertices.size ( ) );

            for ( int d = 1; d < tdim_; ++d )
            {

                // Allocate entry in cell_sub_entities for this cell if it does not exist
                if ( cell >= static_cast < int > ( cell_sub_entities ( d ).size ( ) ) )
                {
                    cell_sub_entities ( d ).resize ( cell + 1 );
                }

                // Get number of sub-entities to create.
                const int num_entities = cell_type.num_regular_entities ( d );

                // Access the cell_entity entry for this cell and dimension.
                std::vector<int>& cell_entity_connectivity = cell_sub_entities ( d )[cell];
                cell_entity_connectivity.clear ( );

                // Add all sub-entities of this cell to the cell type, as
                // well as to the cell->sub-entity connectivity.
                for ( int s = 0; s < num_entities; ++s )
                {
                    const std::vector<int> entity_vertices =
                            cell_type.vertices_of_entity ( d, s, cell_vertices );

                    int entity_number;

                    if ( !find_entity ( d, entity_vertices, &entity_number ) )
                    {
                        entity_number = add_entity ( d, entity_vertices );
                    }
                    cell_sub_entities ( d )[cell].push_back ( entity_number );
                }
            }
        }

        CellType::ConnectionList& CellType::entities ( TDim dim )
        {
            assert ( dim > 0 );
            assert ( dim - 1 < static_cast < int > ( entities_.size ( ) ) );
            return entities_[dim - 1];
        }

        const CellType::ConnectionList& CellType::entities ( TDim dim ) const
        {
            assert ( dim > 0 );
            assert ( dim - 1 < static_cast < int > ( entities_.size ( ) ) );
            return entities_[dim - 1];
        }

        CellType::ConnectionList& CellType::cell_sub_entities ( TDim dim )
        {
            assert ( dim > 0 );
            assert ( dim - 1 < static_cast < int > ( cell_sub_entities_.size ( ) ) );
            return cell_sub_entities_[dim - 1];
        }

        const CellType::ConnectionList& CellType::cell_sub_entities ( TDim dim ) const
        {
            assert ( dim > 0 );
            assert ( dim - 1 < static_cast < int > ( cell_sub_entities_.size ( ) ) );
            return cell_sub_entities_[dim - 1];
        }

        CellType::ConnectionList& CellType::entity_parents ( TDim dim )
        {
            assert ( dim > 0 );
            assert ( dim - 1 < static_cast < int > ( entity_parents_.size ( ) ) );
            return entity_parents_[dim - 1];
        }

        const CellType::ConnectionList& CellType::entity_parents ( TDim dim ) const
        {
            assert ( dim > 0 );
            assert ( dim - 1 < static_cast < int > ( entity_parents_.size ( ) ) );
            return entity_parents_[dim - 1];
        }

        //////////////// END PRIVATE FUNCTIONS ////////////////

        std::ostream& operator<< ( std::ostream& os, const CellType& cell_type )
        {
            os << "CellType: dim = " << cell_type.tdim ( )
                    << ", tag = " << cell_type.tag ( ) << "\n";

            for ( int d = 0; d < cell_type.tdim ( ); ++d )
            {
                os << "Num entities of dimension " << d
                        << " = " << cell_type.num_entities ( d )
                        << " of which " << cell_type.num_regular_entities ( d )
                        << " are regular\n";
            }

            for ( int v = 0; v < cell_type.num_entities ( 0 ); ++v )
            {
                os << "Vertex " << v;
                if ( cell_type.is_regular_vertex ( v ) )
                {
                    os << " is regular.\n";
                }
                else
                {
                    os << " has parents "
                            << string_from_range ( cell_type.vertices_[v].begin ( ),
                                                   cell_type.vertices_[v].end ( ) ) << "\n";
                }
            }

            for ( int d = 1; d <= cell_type.tdim ( ); ++d )
            {
                for ( int i = 0; i < cell_type.num_entities ( d ); ++i )
                {
                    const std::vector<int>& vertices = cell_type.entities ( d )[i];
                    os << "Entity (dim = " << d << ", number = " << i << ")\n";
                    os << "\tvertices = "
                            << string_from_range ( vertices.begin ( ), vertices.end ( ) ) << "\n";

                    const std::vector<int>& parents = cell_type.entity_parents ( d )[i];
                    os << "\tparents = "
                            << string_from_range ( parents.begin ( ), parents.end ( ) ) << "\n";
                    os << "\tregular parent = "
                            << cell_type.regular_parent ( d, i ) << "\n";
                }
            }

            for ( int i = 0; i < cell_type.num_entities ( cell_type.tdim ( ) ); ++i )
            {
                for ( int d = 1; d < cell_type.tdim ( ); ++d )
                {
                    const std::vector<int>& sub_entities = cell_type.cell_sub_entities ( d )[i];
                    os << "Cell " << i << " has sub-entities "
                            << string_from_range ( sub_entities.begin ( ), sub_entities.end ( ) )
                            << " of dimension " << d << "\n";
                }
            }
            return os;
        }

        // Definition of static container of CellType singleton objects.
        std::vector< SharedPtr<CellType>::Type > CellType::cell_types_;

        //////////////////////////////////////////////////
        //////////////// Helper functions ////////////////
        //////////////////////////////////////////////////

        ///
        /// \details The algorithm performs a simple linear search on v2
        /// for each element in v1. Its complexity is hence quadratic,
        /// which means that it will only be efficient for small vectors
        /// v1 and v2.
        ///
        /// \param v1  first vector
        /// \param v2  second vector
        /// \returns true if v1.size() == v2.size and all entries in v1 exist in v2.

        bool compare_as_sets ( const std::vector<int>& v1, const std::vector<int>& v2 )
        {
            // check if sets are the same size
            if ( v1.size ( ) != v2.size ( ) )
            {
                return false;
            }

            // check if all elements of v1 exist in v2, using linear search
            for ( std::vector<int>::const_iterator it = v1.begin ( ); it != v1.end ( ); ++it )
            {
                if ( find ( v2.begin ( ), v2.end ( ), *it ) == v2.end ( ) )
                {
                    return false;
                }
            }
            return true;
        }

        template bool CellType::check_cell_geometry ( const std::vector<double>& coords, GDim gdim ) const;
        template bool CellType::check_cell_geometry ( const std::vector<float>& coords, GDim gdim ) const;
    } // namespace mesh
} // namespace hiflow
