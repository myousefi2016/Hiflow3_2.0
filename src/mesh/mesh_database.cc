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

#include "mesh_database.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

#include "common/log.h"

#include "cell_type.h"

const int DEBUG_LEVEL = 0;

std::string HDF5_GROUP_NAME = "MeshDatabase";
//dataset base names
std::string HDF5_COORDS_NAME = "coordinates";
std::string HDF5_MATERIAL_NUM = "material_numbers";
std::string HDF5_ENTITY_VERTEX = "entity_vertex_connectivity";
std::string HDF5_ENTITY_VERTEX_OFFSETS = "entity_vertex_offsets";
std::string HDF5_VERTEX_ENTITY = "vertex_entity_connectivity";
std::string HDF5_VERTEX_ENTITY_OFFSETS = "vertex_entity_offsets";
std::string HDF5_VERTEX_SEARCH_TABLE = "vertex_search_table";
std::string HDF5_TDIM = "tdim";
std::string HDF5_GDIM = "gdim";

namespace hiflow
{
    namespace mesh
    {

        //////////////// BEGIN MeshDatabase implementation ////////////////

        /// \param tdim           the topological dimension of the database
        /// \param gdim           the geometrical dimension of the database

        MeshDatabase::MeshDatabase ( TDim tdim, GDim gdim )
        : tdim_ ( tdim ),
        gdim_ ( gdim ),
        entity_vertex_connectivities_ ( tdim ),
        vertex_entity_connectivities_ ( tdim ),
        vertex_search_table_ ( new VertexSearchTable ( gdim, coordinates_ ) )
        {

            assert ( tdim_ >= 0 );
            assert ( gdim_ >= 0 );
            //    assert(gdim_ >= tdim_);

        }

        TDim MeshDatabase::tdim ( ) const
        {
            return tdim_;
        }

        GDim MeshDatabase::gdim ( ) const
        {
            return gdim_;
        }

        /// \param point          the coordinates of the vertex (should point to an array of gdim coordinate values)
        /// \returns              the id of the vertex

        Id MeshDatabase::add_vertex ( const Coordinate* point )
        {
            Id id;
            if ( !has_vertex ( point, id ) )
            {
                id = create_vertex ( point );
            }
            return id;
        }

        /// \param destination    the vector of the destination (should point to an array of gdim coordinate values)
        /// \param id             the id of the vertex

        void MeshDatabase::replace_vertex ( const Coordinate* destination, Id vertex_id )
        {
            assert ( coordinates_.size ( ) >= gdim_ * vertex_id + gdim_ - 1 );
            for ( int c = 0; c < gdim_; ++c )
            {
                coordinates_[gdim_ * vertex_id + c] = destination[c];
            }
            vertex_search_table_->update_vertex ( vertex_id );
        }

        /// \param displacement   the vector of the displacement (should point to an array of gdim coordinate values)
        /// \param id             the id of the vertex

        void MeshDatabase::move_vertex ( const Coordinate* displacement, Id vertex_id )
        {
            assert ( coordinates_.size ( ) >= gdim_ * vertex_id + gdim_ - 1 );
            for ( int c = 0; c < gdim_; ++c )
            {
                coordinates_[gdim_ * vertex_id + c] += displacement[c];
            }
            vertex_search_table_->update_vertex ( vertex_id );
        }

        /// \param point          the coordinates of the vertex (should point to an array of gdim coordinate values)
        /// \param id             output parameter which contains the id of the vertex if it was found
        /// \returns              true if the database contains the vertex

        bool MeshDatabase::has_vertex ( const Coordinate* point, Id& id ) const
        {
            return vertex_search_table_->find_vertex ( point, id );
        }

        /// \param vertex_id      the id of the vertex
        /// \returns              a vector containing the coordinates of the vertex

        std::vector<Coordinate> MeshDatabase::get_coordinates ( Id vertex_id ) const
        {
            return std::vector<Coordinate>( &coordinates_[gdim_ * vertex_id],
                    &coordinates_[gdim_ * ( vertex_id + 1 )] );
        }

        /// \param vertex_ids     a vector with the id:s of the requested vertices
        /// \returns              a vector containing the coordinates of the vertices in interleaved form (in 3d: x1 y1 z1 x2 y2 z2, ... )

        std::vector<Coordinate> MeshDatabase::get_coordinates ( const std::vector<Id>& vertex_ids ) const
        {
            std::vector<Coordinate> coords ( gdim_ * vertex_ids.size ( ) );
            int k = 0;
            for ( std::vector<Id>::const_iterator it = vertex_ids.begin ( ); it != vertex_ids.end ( ); ++it )
            {
                for ( int c = 0; c < gdim_; ++c )
                {
                    coords[k++] = coordinates_[gdim_ * ( *it ) + c];
                }
            }
            return coords;
        }

        /// \param x [in] Coordinates of the point to query
        /// \param eps [in] An epsilon range to search in.
        /// \param vertex_id [out] Id of the closest point if it exists in
        /// the epsilon range
        /// \returns true if a closest point exists within the epsilon range

        bool MeshDatabase::find_closest_vertex ( const Coordinate* x, double eps, Id* closest_vertex ) const
        {
            // 1. get candidate points in vertex search table
            // 2. if points exist: compute distances |x - x_hat|
            //    else return false
            // 3. vertex_id = id of x_hat with minimum distance to x
            // 4. return true.
            double distance_x_origin = vertex_search_table_->distance ( x );
            LOG_DEBUG ( 2, "Distance of search point to origin: " << distance_x_origin );
            bool found = false;

            // LOG SEARCH: Scheme: |          |   | |||    |
            //                    min                x    max

            // Create a list of candidates (~500-1000 points), and chose
            // the one closest to the original point within a epsilon
            // range.
            int range = 1000;
            std::vector<Id> vertices_in_db = vertex_search_table_->vertices ( );
            EntityNumber number_of_vertices_in_db = vertices_in_db.size ( );

            assert ( number_of_vertices_in_db != 0 );
            LOG_DEBUG ( 3, "number_of_vertices_in_db == " << number_of_vertices_in_db );

            int min_index = 0;
            int max_index = number_of_vertices_in_db - 1;
            int middle_index = static_cast < int > ( ( min_index + max_index ) / 2 );
            while ( max_index - min_index > range )
            {
                LOG_DEBUG ( 4, "min_index == " << min_index );
                LOG_DEBUG ( 4, "max_index == " << max_index );
                middle_index = static_cast < int > ( ( min_index + max_index ) / 2 );
                // compute distance of middle_index vertex to origin

                double middle_index_distance = vertex_search_table_->distance ( vertices_in_db[middle_index] );
                if ( middle_index_distance < distance_x_origin )
                {
                    min_index = middle_index;
                }
                else
                {
                    max_index = middle_index;
                }
            }
            middle_index = static_cast < int > ( ( min_index + max_index ) / 2 );
            min_index = middle_index - range / 2 >= 0 ? middle_index - range / 2 : 0;
            max_index = middle_index + range / 2 <= number_of_vertices_in_db - 1 ? middle_index + range / 2 : number_of_vertices_in_db - 1;
            LOG_DEBUG ( 4, "(min middle max) == " << min_index << " " << middle_index << " " << max_index );

            // fill candidates list
            std::vector<Id> candidates;
            for ( int i = min_index; i <= max_index; ++i )
            {
                candidates.push_back ( vertices_in_db[i] );
            }
            assert ( !candidates.empty ( ) );

            // compute distance between x and each candidate!
            std::vector<Id>::iterator it = candidates.begin ( );
            *closest_vertex = *it;
            Coordinate minimum = vertex_search_table_->distance ( x, vertex_search_table_->get_coordinates ( *closest_vertex ) );

            for (; it != candidates.end ( ); ++it )
            {
                const Coordinate* cand_coords = vertex_search_table_->get_coordinates ( *it );
                double distance_of_points = vertex_search_table_->distance ( x, cand_coords );
                if ( distance_of_points < minimum )
                {
                    *closest_vertex = *it;
                    minimum = distance_of_points;
                }

                LOG_DEBUG ( 2, "Current minimum distance of point with id "
                            << *it
                            << " from my search point is: " << minimum );
            }
            if ( minimum < eps ) found = true;
            LOG_DEBUG ( 1, "closest_vertex == " << *closest_vertex );
            return found;
        }

        /// \returns              the number of vertices in the database

        EntityCount MeshDatabase::num_vertices ( ) const
        {
            return coordinates_.size ( ) / gdim_;
        }

        /// \param entity_dim    the topological dimension of the entity
        /// \param entity_vertices the vertices contained in the
        /// entity
        /// \param mat_num the material number of the entity, -1 by default.
        /// \returns             the id of the entity

        Id MeshDatabase::add_entity ( TDim entity_dim,
                                      const std::vector<Id>& entity_vertices )
        {

            assert ( entity_dim >= 1 );
            assert ( entity_dim <= tdim_ );

            // check if entity exists
            Id id;

            if ( !has_entity ( entity_dim, entity_vertices, id ) )
            {
                // create entity
                id = create_entity ( entity_dim, entity_vertices );
            }

            return id;
        }

        /// \param entity_dim    the topological dimension of the entity
        /// \param entity_vertices the vertices contained in the entity
        /// \param id            output parameter containing the id of the entity if it was found
        /// \returns             true if the database contains the entity with the given vertices

        bool MeshDatabase::has_entity ( TDim entity_dim,
                                        const std::vector<Id>& entity_vertices,
                                        Id& id ) const
        {
            assert ( entity_dim >= 1 );
            assert ( entity_dim <= tdim_ );

            return vertex_entity_connectivities_[entity_dim - 1].find_entity ( entity_vertices, entity_dim, this, id );
        }

        /// \param entity_dim    the topological dimension of the sought entity
        /// \param id            the id of the sought entity
        /// \returns             true if an entity (entity_dim, id) exists in the database

        bool MeshDatabase::has_entity ( TDim entity_dim, Id id ) const
        {
            assert ( entity_dim >= 0 );
            assert ( entity_dim <= tdim_ );
            return id < num_entities ( entity_dim ) && id >= 0;
        }

        /// \returns             the number of entities of dimension entity_dim

        EntityCount MeshDatabase::num_entities ( TDim entity_dim ) const
        {
            if ( entity_dim == 0 )
            {
                return coordinates_.size ( ) / gdim_;
            }
            else
            {
                assert ( entity_dim > 0 );
                assert ( entity_dim <= tdim_ );
                const Connectivity& connectivity = entity_vertex_connectivity ( entity_dim );
                return connectivity.num_entities ( );
            }
        }

        /// \param entity_dim    the topological dimension of the entity
        /// \param id            the id of the entity
        /// \returns             a VertexIdIterator referencing the first vertex of the entity

        VertexIdIterator MeshDatabase::begin_vertices ( TDim entity_dim, Id id ) const
        {
            const Connectivity& connectivity = entity_vertex_connectivity ( entity_dim );
            return connectivity.begin ( id );
        }

        /// \param entity_dim    the topological dimension of the entity
        /// \param id            the id of the entity
        /// \returns             a VertexIdIterator referencing one past the last vertex of the entity

        VertexIdIterator MeshDatabase::end_vertices ( TDim entity_dim, Id id ) const
        {
            const Connectivity& connectivity = entity_vertex_connectivity ( entity_dim );
            return connectivity.end ( id );
        }

        /// \param entity_dim    the topological dimension of the entity
        /// \param id            the id of the entity
        /// \returns             the number of vertices connected to the entity

        EntityCount MeshDatabase::entity_size ( TDim entity_dim, Id id ) const
        {
            if ( entity_dim > 0 )
            {
                const Connectivity& connectivity = entity_vertex_connectivity ( entity_dim );
                return connectivity.num_connections ( id );
            }
            else
            {
                assert ( id >= 0 );
                assert ( id < num_vertices ( ) );
                return 1;
            }
        }

        void MeshDatabase::build ( TDim dim, const SortedArray<Id>& cells,
                                   SortedArray<Id>& d_entities,
                                   Connectivity& cell_d_connections )
        {
            LOG_DEBUG ( 2, "Building entities of dimension " << dim );
            typedef SortedArray<Id>::const_iterator IdIterator;

            // An implementation of the build() algorithm
            // described in the paper by A. Logg.

            assert ( dim > 0 );
            assert ( dim < tdim ( ) );
            assert ( tdim ( ) > 1 ); // operation does not make sense for tdim() == 1

            for ( IdIterator cell_it = cells.begin ( ); cell_it != cells.end ( ); ++cell_it )
            {
                const Id cell_id = *cell_it;
                // TODO(Staffan) consider if this iteration should not be
                // replaced with EntityNumber/id access, since we are in
                // MeshDatabase and have direct access to entity data.
                const std::vector<Id> cell_vertices ( begin_vertices ( tdim_, cell_id ), end_vertices ( tdim_, cell_id ) );
                const CellType& cell_type = CellType::get_instance ( tdim_, cell_vertices.size ( ) );
                const EntityCount num_sub_entities = cell_type.num_regular_entities ( dim );

                std::vector<Id> connected_entities ( num_sub_entities );

                for ( int s = 0; s < num_sub_entities; ++s )
                {
                    // get vertices of sub-entity
                    std::vector<Id> sub_entity_vertices =
                            cell_type.vertices_of_entity ( dim, s, cell_vertices );

                    const Id sub_entity_id = add_entity ( dim, sub_entity_vertices );

                    // add it to the d_entities vector if it does not already exist
                    d_entities.find_insert ( sub_entity_id );

                    connected_entities[s] = sub_entity_id;
                }

                cell_d_connections.add_connections ( connected_entities );
            }
        }

        void MeshDatabase::compute_incident_entities ( TDim d, const SortedArray<Id>& d_entities,
                                                       Connectivity& d_d_connectivity ) const
        {
            LOG_DEBUG ( 2, "Computing connectivity " << d << " -> " << d );
            assert ( d > 0 );
            assert ( d <= tdim ( ) );
            const Connectivity& d_0_connectivity = entity_vertex_connectivity ( d );
            const Connectivity& zero_d_connectivity = vertex_entity_connectivity ( d );

            ConnectivityAlgorithms::intersect_subset_equal_dimensions ( d_entities,
                                                                        d_0_connectivity,
                                                                        zero_d_connectivity,
                                                                        d_d_connectivity );
        }

        void MeshDatabase::compute_incident_entities ( TDim d1, TDim d2,
                                                       const SortedArray<Id>& d1_entities,
                                                       const SortedArray<Id>& d2_entities,
                                                       Connectivity& d1_d2_connectivity ) const
        {
            LOG_DEBUG ( 2, "Computing connectivity " << d1 << " -> " << d2 );
            assert ( d1 > 0 );
            assert ( d1 <= tdim ( ) );
            assert ( d2 < d1 );

            const Connectivity& d1_0_connectivity = entity_vertex_connectivity ( d1 );
            const Connectivity& d2_0_connectivity = entity_vertex_connectivity ( d2 );
            const Connectivity& zero_d2_connectivity = vertex_entity_connectivity ( d2 );

            ConnectivityAlgorithms::intersect_subset_unequal_dimensions ( d1_entities,
                                                                          d2_entities,
                                                                          d1_0_connectivity,
                                                                          d2_0_connectivity,
                                                                          zero_d2_connectivity,
                                                                          d1_d2_connectivity );
        }

        void MeshDatabase::vertex_entity_incidence ( TDim d, const SortedArray<Id>& vertices,
                                                     const SortedArray<Id>& d_entities,
                                                     Connectivity& zero_d_connectivity ) const
        {
            LOG_DEBUG ( 2, "Computing connectivity 0 -> " << d );
            assert ( d > 0 );
            assert ( d <= tdim ( ) );
            const Connectivity& d_0_connectivity = entity_vertex_connectivity ( d );
            d_0_connectivity.transpose_subset ( d_entities, vertices, zero_d_connectivity );
        }

        /// \pre              the entity is a cell or facet
        /// \param d          the topological dimension of the entity
        /// \param id         id of the entity
        /// \returns          material number of the entity

        MaterialNumber MeshDatabase::get_material_number ( TDim d, Id id ) const
        {
            assert ( d == tdim_ || d == tdim_ - 1 );
            assert ( has_entity ( d, id ) );

            return material_numbers_[tdim_ - d][id];
        }

        /// \pre              the entity is a cell or facet
        /// \param d          the topological dimension of the entity
        /// \param id         id of the entity
        /// \param material   new material number for the entity

        void MeshDatabase::set_material_number ( TDim d, Id id, MaterialNumber material )
        {
            assert ( d == tdim_ || d == tdim_ - 1 );
            assert ( has_entity ( d, id ) );

            material_numbers_[tdim_ - d][id] = material;
        }

        /// \param point          the coordinates of the vertex (should point to an array of gdim coordinate values)
        /// \returns              the id of the newly created vertex

        Id MeshDatabase::create_vertex ( const Coordinate* point )
        {
            const Id id = num_vertices ( );
            coordinates_.insert ( coordinates_.end ( ), point, point + gdim_ );
            for ( int d = 0; d < tdim_; ++d )
            {
                vertex_entity_connectivities_[d].add_vertex ( );
            }
            vertex_search_table_->add_vertex ( point, id );
            return id;
        }

        /// \param entity_dim    the topological dimension of the entity
        /// \param entity_vertices the vertices contained in the
        /// entity
        /// \returns             the id of the newly created entity

        Id MeshDatabase::create_entity ( TDim entity_dim, const std::vector<Id>& entity_vertices )
        {
            assert ( entity_dim >= 1 );
            assert ( entity_dim <= tdim_ );

            // create entity
            Id id = entity_vertex_connectivities_[entity_dim - 1].num_entities ( );
            entity_vertex_connectivities_[entity_dim - 1].add_connections ( entity_vertices );
            for ( std::vector<Id>::const_iterator it = entity_vertices.begin ( );
                  it != entity_vertices.end ( ); ++it )
            {
                vertex_entity_connectivities_[entity_dim - 1].add_vertex_connection ( *it, id );
            }

            if ( entity_dim == tdim_ || entity_dim == tdim_ - 1 )
            {
                // Allocate space for material number for entity if it is a cell or facet
                // Default material number is -1.
                material_numbers_[tdim_ - entity_dim].push_back ( -1 );
            }

            LOG_DEBUG ( 2, "created entity (d=" << entity_dim << ", id=" << id << "): " );
            LOG_DEBUG ( 2, string_from_range ( entity_vertices.begin ( ), entity_vertices.end ( ) ) );

            return id;
        }

        const Connectivity& MeshDatabase::entity_vertex_connectivity ( TDim entity_dim ) const
        {
            assert ( entity_dim > 0 );
            assert ( entity_dim <= tdim_ );
            return entity_vertex_connectivities_[entity_dim - 1];
        }

        const Connectivity& MeshDatabase::vertex_entity_connectivity ( TDim entity_dim ) const
        {
            assert ( entity_dim > 0 );
            assert ( entity_dim <= tdim_ );
            return vertex_entity_connectivities_[entity_dim - 1];
        }

        void MeshDatabase::save ( std::string filename, const MPI_Comm& comm ) const
        {
#ifdef WITH_HDF5
            int rank, num_part;
            int master_rank = 0;
            MPI_Comm_rank ( comm, &rank );
            MPI_Comm_size ( comm, &num_part );

            H5FilePtr file_ptr ( new H5File ( filename, "w", comm ) );
            bool is_writer = ( master_rank == rank );
            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << HDF5_GROUP_NAME;

            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "w" ) );

            // WRITING GDIM AND TDIM
            H5DatasetPtr gdim_dataset_ptr ( new H5Dataset ( group_ptr, 1,
                                                            HDF5_GDIM, "w", &gdim_ ) );
            gdim_dataset_ptr->write ( is_writer, 0, &gdim_ );

            H5DatasetPtr tdim_dataset_ptr ( new H5Dataset ( group_ptr, 1,
                                                            HDF5_TDIM, "w", &tdim_ ) );
            tdim_dataset_ptr->write ( is_writer, 0, &tdim_ );

            // WRITING COORDINATES
            write_array_parallel<Coordinate>( group_ptr, HDF5_COORDS_NAME, coordinates_, comm );

            //WRITING VERTEX SEARCH TABLE
            //currently no way to avoid copying vertices Ids >.>
            std::vector<Id> vertex_table = vertex_search_table_->vertices ( );
            write_array_parallel<Id>( group_ptr, HDF5_VERTEX_SEARCH_TABLE, vertex_table, comm );

            //WRITING MATERIAL NUMBERS
            for ( int mnc = 0; mnc < 2; ++mnc )
            {
                std::stringstream mat_data_name;
                mat_data_name << HDF5_MATERIAL_NUM << mnc;
                write_array_parallel<Id>( group_ptr, mat_data_name.str ( ), material_numbers_[mnc], comm );
            }

            // WRITING CONNECTIVITY
            // FIRST ENTITY -> VERTEX
            assert ( tdim_ >= 2 );
            for ( int dim = 0; dim < tdim_; ++dim )
            {
                const Connectivity& conny = entity_vertex_connectivities_[dim];
                //convert connectivity to vector
                std::vector<int> conn_as_vec ( 0 );
                std::vector<int> offsets ( conny.num_entities ( ) );
                for ( int ent = 0; ent < conny.num_entities ( ); ++ent )
                {
                    int num_conns = 0;
                    for ( ConstConnectionIterator inc = conny.begin ( ent ); inc != conny.end ( ent ); ++inc )
                    {
                        conn_as_vec.push_back ( *inc );
                        ++num_conns;
                    }
                    offsets[ent] = num_conns;
                }

                std::stringstream conn_data_name;
                conn_data_name << HDF5_ENTITY_VERTEX << dim;
                write_array_parallel<int>( group_ptr, conn_data_name.str ( ), conn_as_vec, comm );
                std::stringstream offset_data_name;
                offset_data_name << HDF5_ENTITY_VERTEX_OFFSETS << dim;
                write_array_parallel<int>( group_ptr, offset_data_name.str ( ), offsets, comm );
            }

            // SECOND VERTEX -> ENTITY
            for ( int dim = 0; dim < tdim_; ++dim )
            {
                const VertexConnectivity& vertex_conny = vertex_entity_connectivities_[dim];
                //convert connectivity to vector
                std::vector<int> vertex_conn_as_vec ( 0 );
                std::vector<int> vertex_offsets ( vertex_conny.num_entities ( ) );
                for ( int ent = 0; ent < vertex_conny.num_entities ( ); ++ent )
                {
                    int num_conns = 0;
                    for ( ConstConnectionIterator inc = vertex_conny.begin ( ent ); inc != vertex_conny.end ( ent ); ++inc )
                    {
                        vertex_conn_as_vec.push_back ( *inc );
                        ++num_conns;
                    }
                    vertex_offsets[ent] = num_conns;
                }
                std::stringstream vertex_conn_data_name;
                vertex_conn_data_name << HDF5_VERTEX_ENTITY << dim;
                write_array_parallel<int>( group_ptr, vertex_conn_data_name.str ( ), vertex_conn_as_vec, comm );
                std::stringstream vertex_offset_data_name;
                vertex_offset_data_name << HDF5_VERTEX_ENTITY_OFFSETS << dim;
                write_array_parallel<int>( group_ptr, vertex_offset_data_name.str ( ), vertex_offsets, comm );
            }

#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        void MeshDatabase::load ( std::string filename, const MPI_Comm& comm )
        {
#ifdef WITH_HDF5
            //assert Database is empty?
            int rank, size;
            int master_rank = 0;
            MPI_Comm_rank ( comm, &rank );
            MPI_Comm_size ( comm, &size );
            H5FilePtr file_ptr ( new H5File ( filename, "r", comm ) );
            bool is_reader = ( master_rank == rank );
            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << HDF5_GROUP_NAME;

            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "r" ) );

            // READING GDIM AND TDIM
            H5DatasetPtr gdim_dataset_ptr ( new H5Dataset ( group_ptr, 1,
                                                            HDF5_GDIM, "r", &gdim_ ) );
            gdim_dataset_ptr->read ( is_reader, 0, &gdim_ );
            H5DatasetPtr tdim_dataset_ptr ( new H5Dataset ( group_ptr, 1,
                                                            HDF5_TDIM, "r", &tdim_ ) );
            tdim_dataset_ptr->read ( is_reader, 0, &tdim_ );
            MPI_Bcast ( &tdim_, 1, MPI_INT, master_rank, comm );
            MPI_Bcast ( &gdim_, 1, MPI_INT, master_rank, comm );

            // READING COORDINATES
            read_array_parallel<Coordinate>( group_ptr, HDF5_COORDS_NAME, coordinates_, comm );

            //READING VERTEX SEARCH TABLE
            //pointer so we can avoid copying data
            std::vector<Id>* vertex_table = new std::vector<Id>( 0 );
            read_array_parallel<Id>( group_ptr, HDF5_VERTEX_SEARCH_TABLE, *vertex_table, comm );
            this->vertex_search_table_->set_vertices ( *vertex_table );

            //READING MATERIAL NUMBERS
            for ( int mnc = 0; mnc < 2; ++mnc )
            {
                std::stringstream mat_data_name;
                mat_data_name << HDF5_MATERIAL_NUM << mnc;
                read_array_parallel<Id>( group_ptr, mat_data_name.str ( ), material_numbers_[mnc], comm );
            }

            // READING CONNECTIVITY
            // FIRST ENTITY -> VERTEX
            assert ( tdim_ >= 2 );
            for ( int dim = 0; dim < tdim_; ++dim )
            {
                std::vector<int> conn_as_vec ( 0 );
                std::vector<int> offsets ( 0 );
                std::stringstream conn_data_name;
                conn_data_name << HDF5_ENTITY_VERTEX << dim;
                read_array_parallel<int>( group_ptr, conn_data_name.str ( ), conn_as_vec, comm );
                std::stringstream offset_data_name;
                offset_data_name << HDF5_ENTITY_VERTEX_OFFSETS << dim;
                read_array_parallel<int>( group_ptr, offset_data_name.str ( ), offsets, comm );
                int curr_pos = 0;
                int next_pos = 0;
                for ( int i = 0; i < offsets.size ( ); ++i )
                {
                    next_pos += offsets[i];
                    entity_vertex_connectivities_[dim].add_connections ( std::vector<int>( conn_as_vec.begin ( ) + curr_pos, conn_as_vec.begin ( ) + next_pos ) );
                    curr_pos = next_pos;
                }
            }

            // SECOND VERTEX -> ENTITY
            for ( int dim = 0; dim < tdim_; ++dim )
            {
                std::vector<int> vertex_conn_as_vec ( 0 );
                std::vector<int> vertex_offsets ( 0 );
                std::stringstream vertex_conn_data_name;
                vertex_conn_data_name << HDF5_VERTEX_ENTITY << dim;
                read_array_parallel<int>( group_ptr, vertex_conn_data_name.str ( ), vertex_conn_as_vec, comm );
                std::stringstream vertex_offset_data_name;
                vertex_offset_data_name << HDF5_VERTEX_ENTITY_OFFSETS << dim;
                read_array_parallel<int>( group_ptr, vertex_offset_data_name.str ( ), vertex_offsets, comm );
                int curr_pos = 0;
                int next_pos = 0;
                for ( int i = 0; i < vertex_offsets.size ( ); ++i )
                {
                    next_pos += vertex_offsets[i];
                    vertex_entity_connectivities_[dim].add_connections ( std::vector<int>( vertex_conn_as_vec.begin ( ) + curr_pos, vertex_conn_as_vec.begin ( ) + curr_pos + vertex_offsets[i] ) );
                    curr_pos = next_pos;
                }
            }
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        void MeshDatabase::copy_from ( const MeshDatabasePtr db )
        {
            const int tdim = this->tdim ( );
            assert ( tdim >= 2 );

            // copy vertex coordinates
            assert ( db->coordinates_.size ( ) > 0 );
            this->coordinates_ = db->coordinates_;

            // copy material numbers
            for ( int l = 0; l < 2; ++l )
            {
                this->material_numbers_[l] = db->material_numbers_[l];
            }

            // copy entity_vertex connectivity
            this->entity_vertex_connectivities_.resize ( db->entity_vertex_connectivities_.size ( ) );

            for ( int dim = tdim - 2; dim < tdim; ++dim )
            {
                Connectivity& conny = db->entity_vertex_connectivities_[dim];
                //convert connectivity to vector
                std::vector<int> conn_as_vec ( 0 );
                std::vector<int> offsets ( conny.num_entities ( ) );
                for ( int ent = 0; ent < conny.num_entities ( ); ++ent )
                {
                    int num_conns = 0;
                    for ( ConnectionIterator inc = conny.begin ( ent ); inc != conny.end ( ent ); ++inc )
                    {
                        conn_as_vec.push_back ( *inc );
                        ++num_conns;
                    }
                    offsets[ent] = num_conns;
                }

                int curr_pos = 0;
                int next_pos = 0;
                for ( int i = 0; i < offsets.size ( ); ++i )
                {
                    next_pos += offsets[i];
                    this->entity_vertex_connectivities_[dim].add_connections ( std::vector<int>( conn_as_vec.begin ( ) + curr_pos, conn_as_vec.begin ( ) + next_pos ) );
                    curr_pos = next_pos;
                }
            }

            // copy vertex_entity connections
            this->vertex_entity_connectivities_.resize ( db->vertex_entity_connectivities_.size ( ) );
            for ( int dim = tdim_ - 2; dim < tdim_; ++dim )
            {
                VertexConnectivity& vertex_conny = db->vertex_entity_connectivities_[dim];
                //convert connectivity to vector
                std::vector<int> vertex_conn_as_vec ( 0 );
                std::vector<int> vertex_offsets ( vertex_conny.num_entities ( ) );
                for ( int ent = 0; ent < vertex_conny.num_entities ( ); ++ent )
                {
                    int num_conns = 0;
                    for ( ConnectionIterator inc = vertex_conny.begin ( ent ); inc != vertex_conny.end ( ent ); ++inc )
                    {
                        vertex_conn_as_vec.push_back ( *inc );
                        ++num_conns;
                    }
                    vertex_offsets[ent] = num_conns;
                }

                int curr_pos = 0;
                int next_pos = 0;
                for ( int i = 0; i < vertex_offsets.size ( ); ++i )
                {
                    next_pos += vertex_offsets[i];
                    this->vertex_entity_connectivities_[dim].add_connections ( std::vector<int>( vertex_conn_as_vec.begin ( ) + curr_pos, vertex_conn_as_vec.begin ( ) + curr_pos + vertex_offsets[i] ) );
                    curr_pos = next_pos;
                }
            }

            // copy vertex search table
            std::vector<Id> vertex_table = db->vertex_search_table_->vertices ( );
            this->vertex_search_table_->set_vertices ( vertex_table );
        }

        void MeshDatabase::deep_copy_from ( const MeshDatabasePtr db )
        {
            this->copy_from ( db );
        }

        void MeshDatabase::create_entity_package ( int entity_tdim, const std::vector<Id>& ids, EntityPackage* entities ) const
        {
            const GDim gdim = this->gdim ( );
            assert ( entity_tdim >= 0 );
            assert ( entities != 0 );

            entities->tdim = entity_tdim;
            entities->gdim = gdim;

            std::tr1::unordered_map<Id, int> vertex_map;
            for ( int j = 0; j < ids.size ( ); ++j )
            {
                entities->offsets.push_back ( entities->connections.size ( ) );

                std::vector<Id> vertex_id = std::vector<Id>( this->begin_vertices ( entity_tdim, ids[j] ),
                        this->end_vertices ( entity_tdim, ids[j] ) );
                MaterialNumber material_number = this->get_material_number ( entity_tdim, ids[j] );
                int num_vertices = vertex_id.size ( );

                for ( int v = 0; v < num_vertices; ++v )
                {
                    const Id v_id = vertex_id[v];
                    if ( vertex_map.find ( v_id ) == vertex_map.end ( ) )
                    {
                        std::vector<Coordinate> v_coords = this->get_coordinates ( v_id );
                        const int v_index = entities->coords.size ( ) / gdim;
                        vertex_map[v_id] = v_index;
                        entities->coords.insert ( entities->coords.end ( ), v_coords.begin ( ), v_coords.end ( ) );
                    }
                    entities->connections.push_back ( vertex_map[v_id] );
                }
                entities->material_numbers.push_back ( material_number );
            }
            // last offset
            entities->offsets.push_back ( entities->connections.size ( ) );
        }
        //////////////// END MeshDatabase implementation ////////////////

        //////////////// BEGIN VertexConnectivity implementation ////////////////

        void VertexConnectivity::add_vertex ( )
        {
            connections_.push_back ( std::vector<Id>( ) );
        }

        void VertexConnectivity::add_vertex_connection ( Id vertex_id, Id connected_entity_id )
        {
            assert ( vertex_id >= 0 );
            assert ( vertex_id < static_cast < int > ( connections_.size ( ) ) );

            // OLD IMPLEMENTATION
            // NB: no check that the vertex is not already connected to the entity is performed

            // add new connection to the end
            /*connections_[vertex_id].push_back(connected_entity_id);

            // bubble-sort in reverse direction to put it in sorted position
            std::vector<Id>::reverse_iterator it = connections_[vertex_id].rbegin();
            std::vector<Id>::reverse_iterator next = it;
            ++next;
            std::vector<Id>::reverse_iterator end = connections_[vertex_id].rend();

            while (next != end && *it < *next) {
                std::swap(*it, *next);
                ++it;
                ++next;
            }*/

            // NEW IMPLEMENTATION

            // check if entity is already connected
            std::vector<int>::iterator contained = std::find ( connections_[vertex_id].begin ( ), connections_[vertex_id].end ( ), connected_entity_id );
            if ( contained != connections_[vertex_id].end ( ) )
            {
                return;
            }
            else
            {
                connections_[vertex_id].push_back ( connected_entity_id );
                std::stable_sort ( connections_[vertex_id].begin ( ), connections_[vertex_id].end ( ) );
            }
        }

        bool VertexConnectivity::find_entity ( const std::vector<Id>& entity_vertices,
                                               const TDim tdim,
                                               const MeshDatabase* db,
                                               Id& id ) const
        {
            // typedef's for convenience
            typedef std::vector<Id>::iterator Iterator;
            typedef std::vector<Id>::const_iterator Citerator;

            // Check if the intersection of entities connected to each
            // vertex in entity_vertices is non-empty. If so, it should
            // have only one element, which is returned in id.
            assert ( entity_vertices.size ( ) >= 1 );

            const Id first_vertex = *entity_vertices.begin ( );
            assert ( first_vertex >= 0 );
            assert ( first_vertex < static_cast < int > ( connections_.size ( ) ) );

            // copy entities connected to first vertex
            // TODO(Staffan): return early if this is empty
            std::vector<Id> intersection = connections_[first_vertex];

            LOG_DEBUG ( 3, "first vertex = " << first_vertex << "\tstart intersection: " );
            LOG_DEBUG ( 3, string_from_range ( intersection.begin ( ), intersection.end ( ) ) );

            // intersect with entities connected to each subsequent vertex
            for ( Citerator v_it = entity_vertices.begin ( ) + 1; v_it != entity_vertices.end ( ); ++v_it )
            {
                const Id vertex = *v_it;
                Citerator connected_begin = connections_[vertex].begin ( );
                Citerator connected_end = connections_[vertex].end ( );

                LOG_DEBUG ( 3, "\tEntities connected to " << vertex << ": " );
                LOG_DEBUG ( 3, string_from_range ( connected_begin, connected_end ) );

                //copy all entities that are connected with this vertex into
                //a new vector. In the next iteration we'll use this vector.
                std::vector<Id> new_intersection ( 0 );
                new_intersection.reserve ( intersection.size ( ) );
                for ( Iterator target = intersection.begin ( ); target != intersection.end ( ); ++target )
                {
                    const bool found = std::binary_search ( connected_begin, connected_end, *target );
                    if ( found )
                    {
                        new_intersection.push_back ( *target );
                    }
                }
                intersection = new_intersection;
            }
            LOG_DEBUG ( 3, "Number of entity vertices: " << entity_vertices.size ( ) );
            if ( !intersection.empty ( ) )
            {
                LOG_DEBUG ( 3, "Number of vertices of entity to find: "
                            << db->entity_size ( tdim, intersection.front ( ) ) );
            }

            bool found_unique = false;
            if ( !intersection.empty ( ) )
            {
                id = intersection.front ( );
                if ( static_cast < int > ( entity_vertices.size ( ) ) == static_cast < int > ( db->entity_size ( tdim, id ) ) )
                {
                    found_unique = true;
                    if ( intersection.size ( ) != 1 )
                    {
                        LOG_ERROR ( "Intersection contains " << intersection.size ( ) << " elements:" );
                        LOG_ERROR ( string_from_range ( intersection.begin ( ), intersection.end ( ) ) );
                        assert ( intersection.size ( ) == 1 );
                    }
                }
                LOG_DEBUG ( 3, "Id: " << id );
            }
            return found_unique;
        }
        //////////////// END VertexConnectivity implementation ////////////////

        //////////////// BEGIN VertexSearchTable implementation ////////////////

        VertexSearchTable::VertexSearchTable ( GDim gdim, const std::vector<Coordinate>& coordinates )
        : vertices_ ( 0 ), gdim_ ( gdim ), coordinates_ ( coordinates )
        {
        }

        void VertexSearchTable::add_vertex ( const Coordinate* point, Id id )
        {
            if ( vertices_.empty ( ) )
            {
                vertices_.push_back ( id );
            }
            else
            {
                if ( vertices_.size ( ) == vertices_.capacity ( ) )
                {
                    vertices_.reserve ( static_cast < int > ( vertices_.capacity ( ) * 1.5 ) );
                }
                // special binary search with only one cmp per iteration
                // http://en.wikipedia.org/wiki/Binary_search_algorithm
                const Coordinate d_point = distance ( point );
                int lower = 0;
                int upper = vertices_.size ( );

                while ( lower < upper )
                {
                    const int pos = lower + static_cast < int > ( ( upper - lower ) / 2 );
                    const Coordinate d_cur = distance ( vertices_[pos] );
                    if ( d_cur < d_point )
                    {
                        lower = pos + 1;
                    }
                    else
                    {
                        upper = pos;
                    }
                }

                assert ( lower == upper );

                // Here lower == upper and d_cur >= d_point.
                vertices_.insert ( vertices_.begin ( ) + lower, id );
            }
        }

        /// \param id
        /// \return id

        bool VertexSearchTable::find_vertex ( const Coordinate* point, Id& id ) const
        {
            // handle special case where no vertices exist
            if ( vertices_.empty ( ) )
            {
                return false;
            }

            const Coordinate d_target = distance ( point );

            int lower = 0;
            int upper = vertices_.size ( );

            do
            {
                const int pos = lower + static_cast < int > ( ( upper - lower ) / 2 );
                const Coordinate d_cur = distance ( vertices_[pos] );
                // TODO(Thomas): See if more compares are needed -> this
                // was a really ugly bug fix!

                // The second part in the if-condition is necessary
                // because of the floating point comparison issue. The "<"
                // comparison is "too exact" and falls in the trap of
                // rounding mistakes, so it gives false
                // positives. Therefore we have to ensure that the two
                // points aren:t actually the same one by a less precise
                // and more controllable comparison.
                if ( d_cur < d_target && !compare_distances ( d_target, d_cur ) )
                {
                    lower = pos + 1;
                }
                else
                {
                    upper = pos;
                }
            }
            while ( lower < upper );

            assert ( lower == upper );
            assert ( lower >= 0 );
            assert ( lower <= static_cast < int > ( vertices_.size ( ) ) );

            // If lower == vertices_.size(), this means d_target >
            // distance() of all existing points. We could however have
            // the case that d_target is close to distance() of the last
            // point, so we let lower point to the last point, and
            // continue the comparison.
            //
            // TODO(Staffan): the effect of using the < operator with
            // floating point numbers must be investigated. As it is now,
            // we rely on the incorrect assumption that the
            // compare_distances function is associative. It should work
            // in most cases, but for points that are really close, we
            // might get in trouble.
            if ( lower == static_cast < int > ( vertices_.size ( ) ) )
            {
                --lower;
            }

            const Coordinate d_cur = distance ( vertices_[lower] );

            if ( compare_distances ( d_target, d_cur ) )
            {
                // Find neighbors which could also match point. These are
                // stored in ascending order of distance.
                // TODO(Staffan): could perhaps be simplified with
                // upper_bound and lower_bound functions
                std::vector<int> positions;
                int neighbor_pos = lower - 1;
                while ( neighbor_pos >= 0 &&
                        compare_distances ( d_target, distance ( vertices_[neighbor_pos] ) ) )
                {
                    positions.push_back ( neighbor_pos );
                    --neighbor_pos;
                }
                std::reverse ( positions.begin ( ), positions.end ( ) );

                positions.push_back ( lower );

                neighbor_pos = lower + 1;
                while ( neighbor_pos < static_cast < int > ( vertices_.size ( ) )
                        && compare_distances ( d_target, distance ( vertices_[neighbor_pos] ) ) )
                {
                    positions.push_back ( neighbor_pos );
                    ++neighbor_pos;
                }

                // compare all coordinates for candidate points
                for ( std::vector<int>::const_iterator it = positions.begin ( );
                      it != positions.end ( ); ++it )
                {
                    const Id neighbor_id = vertices_[*it];
                    const Coordinate* pt = get_coordinates ( neighbor_id );
                    if ( compare_points ( pt, point ) )
                    {
                        id = neighbor_id;
                        return true;
                    }
                }
            }
            return false;
        }

        Coordinate VertexSearchTable::distance ( const Coordinate* point ) const
        {
            Coordinate r = 0.0;

            for ( int i = 0; i < gdim_; ++i )
            {
                r += point[i] * point[i];
            }

            return r;
        }

        Coordinate VertexSearchTable::distance ( Id id ) const
        {
            const Coordinate* coords = get_coordinates ( id );
            return distance ( coords );
        }

        Coordinate VertexSearchTable::distance ( const Coordinate* point1, const Coordinate* point2 ) const
        {
            Coordinate dist = static_cast < Coordinate > ( 0 );
            for ( int i = 0; i < gdim_; ++i )
            {
                dist += ( point1[i] - point2[i] ) * ( point1[i] - point2[i] );
            }

            return dist;
        }

        std::vector<Id> VertexSearchTable::vertices ( ) const
        {
            return vertices_;
        }

        void VertexSearchTable::set_vertices ( std::vector<Id>& vertices )
        {
            vertices_ = vertices;
            return;
        }
        // Helper function for comparing floating point numbers

        bool compare_floating_point ( double x1, double x2, double abs_tol, double rel_tol )
        {
            // TODO(Staffan): consider replacing this with the two-complements
            // technique described on
            // http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
            if ( std::abs ( x1 - x2 ) < abs_tol )
            {
                return true;
            }
            else
            {
                return std::abs ( x1 - x2 ) < rel_tol * std::max ( std::abs ( x1 ), std::abs ( x2 ) );
            }
        }

        bool VertexSearchTable::compare_distances ( Coordinate d1, Coordinate d2 ) const
        {
            assert ( d1 >= 0.0 );
            assert ( d2 >= 0.0 );
            // TODO(Staffan): analyze impact of the tolerance values for
            // comparing distances.
            //
            // How can we guarantee
            // that returning false here implies returning false also in
            // compare_points() ?
            const Coordinate ABS_DIST_TOL = 1.0e-6;
            const Coordinate REL_DIST_TOL = 1.0e-6;
            return compare_floating_point ( d1, d2, ABS_DIST_TOL, REL_DIST_TOL );
        }

        bool VertexSearchTable::compare_points ( const Coordinate* p1, const Coordinate* p2 ) const
        {
            assert ( p1 != 0 );
            assert ( p2 != 0 );

            // TODO(Staffan): analyze the impact of the tolerance values
            // for comparing coordinates.
            const Coordinate ABS_COORD_TOL = 1.0e-6;
            const Coordinate REL_COORD_TOL = 1.0e-6;

            for ( int i = 0; i < gdim_; ++i )
            {
                if ( !compare_floating_point ( p1[i], p2[i], ABS_COORD_TOL, REL_COORD_TOL ) )
                {
                    return false;
                }
            }
            return true;
        }

        const Coordinate* VertexSearchTable::get_coordinates ( Id id ) const
        {
            assert ( id >= 0 );
            assert ( gdim_ * id < static_cast < int > ( coordinates_.size ( ) ) );
            return &coordinates_[gdim_ * id];
        }

        void VertexSearchTable::update_vertex ( Id id )
        {
            int pos = -1;
            int num_vert = vertices_.size ( );
            for ( int i = 0; i < num_vert; ++i )
            {
                if ( vertices_[i] == id )
                {
                    pos = i;
                    break;
                }
            }
            if ( pos != -1 )
            { // if vertex already exists
                const Coordinate dist = distance ( id );
                //check if its still in the right position.
                //if not: delete it and re-add it.
                if ( pos > 0 )
                {
                    const Coordinate dist_low = this->distance ( vertices_[pos - 1] );
                    if ( dist_low > dist )
                    {
                        vertices_.erase ( vertices_.begin ( ) + pos );
                        this->add_vertex ( &coordinates_[gdim_ * id], id );
                        return;
                    }
                }
                if ( pos < num_vert - 1 )
                {
                    const Coordinate dist_upp = this->distance ( vertices_[pos + 1] );
                    if ( dist_upp < dist )
                    {
                        vertices_.erase ( vertices_.begin ( ) + pos );
                        this->add_vertex ( &coordinates_[gdim_ * id], id );
                        return;
                    }
                }
                //do nothing if it is in the right position
                return;
            }
            //add vertex if it didn't exist before.
            this->add_vertex ( &coordinates_[gdim_ * id], id );
            return;
        }
        //////////////// END VertexSearchTable implementation ////////////////

    } // end namespace mesh
} // namespace hiflow
