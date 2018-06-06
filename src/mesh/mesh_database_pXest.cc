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

#include "mesh_database_pXest.h"
#include "mesh_pXest.h"
#include "common/log.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace mesh
    {

        //////////////// BEGIN MeshDatabase implementation ////////////////

        /// \param tdim the topological dimension of the database
        /// \param gdim the geometrical dimension of the database

        MeshPXestDatabase::MeshPXestDatabase ( TDim tdim, GDim gdim )
        : MeshDatabase ( tdim, gdim ), layer_width_(-1)
#ifdef WITH_P4EST
        ,
        p4est_conn_ ( NULL ),
        p8est_conn_ ( NULL ),
        p4est_forest_ ( NULL ),
        p8est_forest_ ( NULL ),
        p4est_ghost_ ( NULL ),
        p8est_ghost_ ( NULL ),
        quadrant_arrays_sorted_ ( false ),
        ctype4_ ( P4EST_CONNECT_FACE ),
        ctype8_ ( P8EST_CONNECT_FACE )
#endif
        {
            this->entity_to_quad_.resize ( tdim + 1 );
            this->coarse_entity_to_quad_.resize ( tdim + 1 );
        }

        void MeshPXestDatabase::build ( TDim dim, const SortedArray<Id>& cells,
                                        SortedArray<Id>& d_entities,
                                        Connectivity& cell_d_connections )
        {
            LOG_DEBUG ( 1, "Building entities of dimension " << dim );
            typedef SortedArray<Id>::const_iterator IdIterator;

            assert ( dim > 0 );
            assert ( dim < tdim ( ) );
            assert ( tdim ( ) > 1 ); // operation does not make sense for tdim() == 1

            // only works for quadrilaterals and hexahedrons
            int num_sub_entities = 0;
            if ( this->tdim ( ) == 2 )
            {
                num_sub_entities = 4;
            }
            else if ( this->tdim ( ) == 3 )
            {
                switch ( dim )
                {
                    case 1:
                        num_sub_entities = 12;
                        break;
                    case 2:
                        num_sub_entities = 6;
                        break;
                }
            }

            /* sort + unique vector approach */
            std::vector<Id> tmp_d_entities;
            tmp_d_entities.reserve ( d_entities.data ( ).size ( ) + cells.size ( ) * num_sub_entities );
            tmp_d_entities.insert ( tmp_d_entities.begin ( ), d_entities.data ( ).begin ( ), d_entities.data ( ).end ( ) );

            /*set approach*/
            /*
            std::set<Id> tmp_d_entities;
            for( unsigned i = 0; i < d_entities.data().size(); ++i )
            {
                tmp_d_entities.insert( d_entities.data()[i] );
            }
             */

            /* std approach */
            //d_entities.reserve (d_entities.data().size() + cells.size() * num_sub_entities);

            // An implementation of the build() algorithm
            // described in the paper by A. Logg.
            for ( IdIterator cell_it = cells.begin ( ); cell_it != cells.end ( ); ++cell_it )
            {
                const Id cell_id = *cell_it;
                // TODO(Staffan) consider if this iteration should not be
                // replaced with EntityNumber/id access, since we are in
                // MeshDatabase and have direct access to entity data.
                const std::vector<Id> cell_vertices ( begin_vertices ( tdim_, cell_id ), end_vertices ( tdim_, cell_id ) );
                const CellType& cell_type = CellType::get_instance ( tdim_, cell_vertices.size ( ) );
                //const EntityCount num_sub_entities = cell_type.num_regular_entities ( dim );

                std::vector<Id> connected_entities ( num_sub_entities );

                for ( int s = 0; s < num_sub_entities; ++s )
                {
                    // get vertices of sub-entity
                    std::vector<Id> sub_entity_vertices =
                            cell_type.vertices_of_entity ( dim, s, cell_vertices );

                    const Id sub_entity_id = add_entity ( dim, sub_entity_vertices );

                    // add it to the d_entities vector if it does not already exist

                    /* sort + unique vector approach */
                    tmp_d_entities.push_back ( sub_entity_id );

                    /*set approach*/
                    //tmp_d_entities.insert (sub_entity_id);

                    /* std approach */
                    //d_entities.find_insert ( sub_entity_id );

                    connected_entities[s] = sub_entity_id;
                }

                cell_d_connections.add_connections ( connected_entities );
            }
            /* sort + unique vector approach */
            std::sort ( tmp_d_entities.begin ( ), tmp_d_entities.end ( ) );
            tmp_d_entities.erase ( unique ( tmp_d_entities.begin ( ), tmp_d_entities.end ( ) ), tmp_d_entities.end ( ) );
            d_entities.data ( ) = tmp_d_entities;

            /*set approach*/
            //d_entities.data().assign(tmp_d_entities.begin(), tmp_d_entities.end());
        }

        void MeshPXestDatabase::set_conn_data ( int num_vertices, int num_cells, std::vector<double>& vertices, std::vector<int>& tree_to_vertices )
        {
            assert ( vertices.size ( ) > 0 );
            assert ( tree_to_vertices.size ( ) > 0 );

            this->conn_vertices_.resize ( 0 );
            this->conn_vertices_.insert ( this->conn_vertices_.end ( ), vertices.begin ( ), vertices.end ( ) );

            this->tree_to_vertices_.resize ( 0 );
            this->tree_to_vertices_.insert ( this->tree_to_vertices_.end ( ), tree_to_vertices.begin ( ), tree_to_vertices.end ( ) );

            this->num_conn_vertices_ = num_vertices;
            this->num_conn_cells_ = num_cells;
        }

#ifdef WITH_P4EST

        void MeshPXestDatabase::set_p4est_conn ( p4est_connectivity_t* conn )
        {
            this->p4est_conn_ = conn;

            // update conn data structs in database
            this->num_conn_vertices_ = conn->num_vertices;
            this->num_conn_cells_ = conn->num_trees;
            this->conn_vertices_ .resize ( 2 * this->num_conn_vertices_ );
            this->tree_to_vertices_ .resize ( 4 * this->num_conn_cells_ );

            for ( int l = 0; l<this->conn_vertices_.size ( ); ++l )
            {
                this->conn_vertices_[l] = conn->vertices[l];
            }
            for ( int l = 0; l<this->tree_to_vertices_.size ( ); ++l )
            {
                this->tree_to_vertices_[l] = conn->tree_to_vertex[l];
            }
        }

        void MeshPXestDatabase::set_p8est_conn ( p8est_connectivity_t* conn )
        {
            this->p8est_conn_ = conn;

            // update conn data structs in database
            this->num_conn_vertices_ = conn->num_vertices;
            this->num_conn_cells_ = conn->num_trees;
            this->conn_vertices_ .resize ( 3 * this->num_conn_vertices_ );
            this->tree_to_vertices_ .resize ( 8 * this->num_conn_cells_ );

            for ( int l = 0; l<this->conn_vertices_.size ( ); ++l )
            {
                this->conn_vertices_[l] = conn->vertices[l];
            }
            for ( int l = 0; l<this->tree_to_vertices_.size ( ); ++l )
            {
                this->tree_to_vertices_[l] = conn->tree_to_vertex[l];
            }
        }
#endif

        MeshPXestDatabase::~MeshPXestDatabase ( )
        {
#ifdef WITH_P4EST
            if ( this->p4est_conn_ != NULL )
            {
                p4est_connectivity_destroy ( this->p4est_conn_ );
                p4est_destroy ( this->p4est_forest_ );
            }
            if ( this->p8est_conn_ != NULL )
            {
                p8est_connectivity_destroy ( p8est_conn_ );
                p8est_destroy ( p8est_forest_ );
            }
#endif
        }

        void MeshPXestDatabase::add_cell_to_mesh_maps ( const std::vector<Id>& cells, std::vector<int>& history_indices )
        {
            assert ( cells.size ( ) == history_indices.size ( ) );
            LOG_DEBUG ( 2, "Add the following cell id to mesh index maps: " << string_from_range ( cells.begin ( ), cells.end ( ) ) << "  =>  "
                        << string_from_range ( history_indices.begin ( ), history_indices.end ( ) ) );

            for ( int l = 0; l < cells.size ( ); ++l )
            {
                Id cell_id = cells[l];
                int ind = history_indices[l];

                std::map<Id, int>::iterator it = this->cell_to_mesh_mapping_.find ( cell_id );
                if ( it == cell_to_mesh_mapping_.end ( ) )
                    cell_to_mesh_mapping_[cell_id] = ind;
                else
                {
                    it->second = ind;
                }
            }
        }

        int MeshPXestDatabase::get_last_mesh_history_index ( Id cell_id )
        {
            if ( this->cell_to_mesh_mapping_.count ( cell_id ) == 0 )
                return -1;
            return ( this->cell_to_mesh_mapping_[cell_id] );
        }

        void MeshPXestDatabase::set_mesh ( MeshPtr mesh_ptr, int history_index )
        {
            assert ( history_index >= 0 );
            LOG_DEBUG ( 2, "Existing MeshPtr at index " << history_index << " is " << this->mesh_history_[history_index] << ", to be inserted MeshPtr is " << mesh_ptr );

            std::pair < std::map<int, MeshPtr >::iterator, bool> ret = this->mesh_history_.insert ( std::make_pair ( history_index, mesh_ptr ) );

            // if element already exists -> update
            if ( !ret.second )
            {
                this->mesh_history_[history_index] = mesh_ptr;
            }
            LOG_DEBUG ( 2, "New MeshPtr at index " << history_index << " is " << this->mesh_history_[history_index] );
        }

        MeshPtr MeshPXestDatabase::get_mesh ( int history_index )
        {
            assert ( history_index >= 0 );
            assert ( history_index <= mesh_history_.size ( ) );
            assert ( this->mesh_history_[history_index] != NULL );
            return this->mesh_history_[history_index];
        }

        void MeshPXestDatabase::add_entity_to_quad_coord ( int tdim, Id entity_id, QuadCoord& coord, int which_mapping )
        {
            assert ( tdim >= 0 );
            assert ( tdim <= this->tdim ( ) );
            assert ( entity_id >= 0 );
            assert ( coord.is_valid ( ) );

            EntityToQuadMap* map;
            switch ( which_mapping )
            {
                case 0:
                {
                    map = &this->coarse_entity_to_quad_[tdim];
                    break;
                }
                case 1:
                {
                    map = &this->entity_to_quad_[tdim];
                    break;
                }
                default:
                {
                    map = NULL;
                    return;
                    break;
                }
            }
            std::pair < EntityToQuadMap::iterator, bool> ret = map->insert ( std::make_pair ( entity_id, coord ) );
            // if element alreayd exists -> update
            if ( !ret.second )
            {
                ( *map )[entity_id] = coord;
            }
        }

        void MeshPXestDatabase::clear_entity_to_quad_map ( int tdim, int which_mapping )
        {
            assert ( tdim >= 0 );
            assert ( tdim <= this->tdim ( ) );

            if ( which_mapping == 0 )
            {
                this->coarse_entity_to_quad_[tdim].clear ( );
            }
            else
            {
                this->entity_to_quad_[tdim].clear ( );
            }
        }

        void MeshPXestDatabase::permute_entity_to_quad_map ( const std::vector<size_t>& perm, int tdim, int which_mapping )
        {
            // so far, only this arguments have been tested
            assert ( which_mapping == 0 );
            assert ( tdim == this->tdim_ );

            EntityToQuadMap* map;
            switch ( which_mapping )
            {
                case 0:
                {
                    map = &this->coarse_entity_to_quad_[tdim];
                    assert ( map->size ( ) == this->num_conn_cells_ );
                    break;
                }
                case 1:
                {
                    map = &this->entity_to_quad_[tdim];
                    break;
                }
                default:
                {
                    map = NULL;
                    return;
                    break;
                }
            }

            for ( EntityToQuadMap::iterator it = map->begin ( ); it != map->end ( ); ++it )
            {
                int32_t old_tree = it->second.tree;
                it->second.tree = perm[old_tree];
                assert ( perm[old_tree] >= 0 );
                assert ( perm[old_tree] < this->num_conn_cells_ );
            }
        }

        QuadCoord MeshPXestDatabase::get_entity_to_quad_coord ( int tdim, Id entity_id, int which_mapping )
        {
            assert ( tdim >= 0 );
            assert ( tdim <= this->tdim ( ) );
            assert ( entity_id >= 0 );

            switch ( which_mapping )
            {
                case 0:
                {
                    assert ( this->coarse_entity_to_quad_[tdim].find ( entity_id ) != this->coarse_entity_to_quad_[tdim].end ( ) );
                    return this->coarse_entity_to_quad_[tdim][entity_id];
                    break;
                }
                case 1:
                {
                    assert ( this->entity_to_quad_[tdim].find ( entity_id ) != this->entity_to_quad_[tdim].end ( ) );
                    return this->entity_to_quad_[tdim][entity_id];
                    break;
                }
                default:
                {
                    assert ( this->entity_to_quad_[tdim].find ( entity_id ) != this->entity_to_quad_[tdim].end ( ) );
                    return this->entity_to_quad_[tdim][entity_id];
                    break;
                }
            }
        }

        EntityToQuadMap* MeshPXestDatabase::get_entity_to_quad_map ( int tdim, int which_mapping )
        {
            assert ( tdim >= 0 );
            assert ( tdim <= this->tdim ( ) );
            switch ( which_mapping )
            {
                case 0:
                {
                    return &this->coarse_entity_to_quad_[tdim];
                    break;
                }
                case 1:
                {
                    return &this->entity_to_quad_[tdim];
                    break;
                }
                default:
                {
                    return NULL;
                    break;
                }
            }
        }

        void MeshPXestDatabase::set_entity_to_quad_map ( int tdim, int which_mapping, EntityToQuadMap& map )
        {
            assert ( tdim >= 0 );
            assert ( tdim <= this->tdim ( ) );
            switch ( which_mapping )
            {
                case 0:
                {
                    this->coarse_entity_to_quad_[tdim] = map;
                    break;
                }
                case 1:
                {
                    this->entity_to_quad_[tdim] = map;
                    break;
                }
                default:
                {
                    break;
                }
            }
        }

        void MeshPXestDatabase::save ( std::string filename, const MPI_Comm& comm ) const
        {
#ifdef WITH_HDF5
            const int tdim = this->tdim ( );
            MeshDatabase::save ( filename, comm );

#    ifdef WITH_P4EST
            std::size_t pos_suffix = filename.find_last_of ( "." );

            std::string filename_without_suffix = filename.substr ( 0, pos_suffix );

            assert ( !filename_without_suffix.empty ( ) );

            std::string forest_filename = filename_without_suffix + "_forest.pXest";
            std::string connectivity_filename = filename_without_suffix + "_connectivity.pXest";

            if ( tdim == 2 )
            {
                p4est_save ( forest_filename.c_str ( ), p4est_forest_, 1 );
            }

            if ( tdim == 3 )
            {
                p8est_save ( forest_filename.c_str ( ), p8est_forest_, 1 );
            }
#    endif
            int rank, num_part;
            MPI_Comm_rank ( comm, &rank );
            MPI_Comm_size ( comm, &num_part );

            H5FilePtr file_ptr ( new H5File ( filename, "w", comm ) );
            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "MeshPXestDatabase";

            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "w" ) );

            // save EntityToQuadMappings entity_to_quad and coarse_entity_to_quad
            for ( int i = 0; i < tdim + 1; ++i )
            {
                const EntityToQuadMap& m = entity_to_quad_[i];
                std::vector<int> keys;
                std::vector<int32_t> values;
                keys.reserve ( m.size ( ) );
                values.reserve ( 8 * m.size ( ) );
                for ( std::map<int, QuadCoord>::const_iterator it = m.begin ( ); it != m.end ( ); ++it )
                {
                    keys.push_back ( it->first );
                    QuadCoord temp_coord = it->second;
                    values.push_back ( temp_coord.tree );
                    values.push_back ( ( int32_t ) temp_coord.level );
                    values.push_back ( temp_coord.x );
                    values.push_back ( temp_coord.y );
                    values.push_back ( temp_coord.z );
                    values.push_back ( ( int32_t ) temp_coord.localId );
                    values.push_back ( ( int32_t ) temp_coord.tdim );
                    values.push_back ( ( int32_t ) temp_coord.coarseId );
                }
                assert ( keys.size ( ) == m.size ( ) );
                assert ( values.size ( ) == 8 * m.size ( ) );
                std::stringstream keys_name;
                keys_name << "entity_to_quad" << i << "_keys";
                std::stringstream values_name;
                values_name << "entity_to_quad" << i << "_values";
                write_array_parallel<int>( group_ptr, keys_name.str ( ), keys, comm );
                write_array_parallel<int32_t>( group_ptr, values_name.str ( ), values, comm );

                const EntityToQuadMap& m2 = coarse_entity_to_quad_[i];
                keys.clear ( );
                values.clear ( );
                keys.reserve ( m2.size ( ) );
                values.reserve ( 8 * m2.size ( ) );
                for ( std::map<int, QuadCoord>::const_iterator it = m2.begin ( ); it != m2.end ( ); ++it )
                {
                    keys.push_back ( it->first );
                    const QuadCoord& temp_coord = it->second;
                    values.push_back ( temp_coord.tree );
                    values.push_back ( ( int32_t ) temp_coord.level );
                    values.push_back ( temp_coord.x );
                    values.push_back ( temp_coord.y );
                    values.push_back ( temp_coord.z );
                    values.push_back ( ( int32_t ) temp_coord.localId );
                    values.push_back ( ( int32_t ) temp_coord.tdim );
                    values.push_back ( ( int32_t ) temp_coord.coarseId );
                }
                assert ( keys.size ( ) == m2.size ( ) );
                assert ( values.size ( ) == 8 * m2.size ( ) );
                std::stringstream coarse_keys_name;
                coarse_keys_name << "coarse_entity_to_quad" << i << "_keys";
                std::stringstream coarse_values_name;
                coarse_values_name << "coarse_entity_to_quad" << i << "_values";
                write_array_parallel<int>( group_ptr, coarse_keys_name.str ( ), keys, comm );
                write_array_parallel<int32_t>( group_ptr, coarse_values_name.str ( ), values, comm );
            }

            // save p4est_connectivity related stuff
            write_array_parallel<double>( group_ptr, "conn_vertices", conn_vertices_, comm );
            write_array_parallel<int>( group_ptr, "tree_to_vertices", tree_to_vertices_, comm );
            write_array_parallel<int>( group_ptr, "num_conn_vertices", &num_conn_vertices_, 1, comm );
            write_array_parallel<int>( group_ptr, "num_conn_cells", &num_conn_cells_, 1, comm );

            //save cell_to_mesh_mapping_
            std::vector<int> mapping_values ( cell_to_mesh_mapping_.size ( ) );
            std::vector< Id > mapping_keys ( cell_to_mesh_mapping_.size ( ) );
            int counter = 0;
            for ( std::map< Id, int >::const_iterator it = cell_to_mesh_mapping_.begin ( ); it != cell_to_mesh_mapping_.end ( ); ++it )
            {
                mapping_keys[counter] = it->first;
                mapping_values[counter] = it->second;
                ++counter;
            }

            write_array_parallel<Id>( group_ptr, "cell_to_mesh_mapping_keys", mapping_keys, comm );
            write_array_parallel<int>( group_ptr, "cell_to_mesh_mapping_values", mapping_values, comm );

            //save individual values
            write_array_parallel<int>( group_ptr, "layer_width", &layer_width_, 1, comm );
            int quadrant_arrays_sorted_int = ( int ) quadrant_arrays_sorted_;
            write_array_parallel<int>( group_ptr, "quadrant_arrays_sorted", &quadrant_arrays_sorted_int, 1, comm );
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        void MeshPXestDatabase::load ( std::string filename, const MPI_Comm& comm )
        {
#ifdef WITH_HDF5
            const int tdim = this->tdim ( );
            MeshDatabase::load ( filename, comm );

#    ifdef WITH_P4EST
            std::size_t pos_suffix = filename.find_last_of ( "." );

            std::string filename_without_suffix = filename.substr ( 0, pos_suffix );

            assert ( !filename_without_suffix.empty ( ) );

            std::string forest_filename = filename_without_suffix + "_forest.pXest";
            if ( tdim == 2 )
            {
                p4est_forest_ = p4est_load ( forest_filename.c_str ( ), comm, QuadData_size ( tdim ), 1, NULL, &p4est_conn_ );
            }

            if ( tdim == 3 )
            {
                p8est_forest_ = p8est_load ( forest_filename.c_str ( ), comm, QuadData_size ( tdim ), 1, NULL, &p8est_conn_ );
            }
#    endif
            int rank, size;
            MPI_Comm_rank ( comm, &rank );
            MPI_Comm_size ( comm, &size );
            H5FilePtr file_ptr ( new H5File ( filename, "r", comm ) );
            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "MeshPXestDatabase";

            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "r" ) );

            // read EntityToQuadMappings entity_to_quad and coarse_entity_to_quad
            for ( int i = 0; i < tdim + 1; ++i )
            {
                std::vector<int> keys;
                std::vector<int32_t> values;
                std::stringstream keys_name;
                keys_name << "entity_to_quad" << i << "_keys";
                std::stringstream values_name;
                values_name << "entity_to_quad" << i << "_values";
                read_array_parallel<int>( group_ptr, keys_name.str ( ), keys, comm );
                read_array_parallel<int32_t>( group_ptr, values_name.str ( ), values, comm );
                assert ( 8 * keys.size ( ) == values.size ( ) );
                EntityToQuadMap& m1 = entity_to_quad_[i];
                m1.clear ( );
                int num_quads = keys.size ( );
                for ( int j = 0; j < num_quads; ++j )
                {
                    QuadCoord new_quad ( values[8 * j], values[8 * j + 1], values[8 * j + 2], values[8 * j + 3], values[8 * j + 4], values[8 * j + 5], values[8 * j + 6], values[8 * j + 7] );
                    m1.insert ( std::make_pair ( keys[j], new_quad ) );
                }

                keys.clear ( );
                values.clear ( );
                std::stringstream coarse_keys_name;
                coarse_keys_name << "coarse_entity_to_quad" << i << "_keys";
                std::stringstream coarse_values_name;
                coarse_values_name << "coarse_entity_to_quad" << i << "_values";
                read_array_parallel<int>( group_ptr, coarse_keys_name.str ( ), keys, comm );
                read_array_parallel<int32_t>( group_ptr, coarse_values_name.str ( ), values, comm );
                assert ( 8 * keys.size ( ) == values.size ( ) );
                EntityToQuadMap& m2 = coarse_entity_to_quad_[i];
                m2.clear ( );
                num_quads = keys.size ( );
                for ( int j = 0; j < num_quads; ++j )
                {
                    QuadCoord new_quad ( values[8 * j], values[8 * j + 1], values[8 * j + 2], values[8 * j + 3], values[8 * j + 4], values[8 * j + 5], values[8 * j + 6], values[8 * j + 7] );
                    m2.insert ( std::make_pair ( keys[j], new_quad ) );
                }
            }

            // read p4est_connectivity related stuff
            std::vector<int> temp_container;
            read_array_parallel<double>( group_ptr, "conn_vertices", conn_vertices_, comm );
            read_array_parallel<int>( group_ptr, "tree_to_vertices", tree_to_vertices_, comm );
            read_array_parallel<int>( group_ptr, "num_conn_vertices", temp_container, comm );
            num_conn_vertices_ = temp_container[0];
            read_array_parallel<int>( group_ptr, "num_conn_cells", temp_container, comm );
            num_conn_cells_ = temp_container[0];

            // read cell_to_mesh_mapping_
            std::vector<int> mapping_values ( 0 );
            std::vector<int> offsets ( 0 );
            std::vector< Id > mapping_keys ( 0 );

            read_array_parallel<Id>( group_ptr, "cell_to_mesh_mapping_keys", mapping_keys, comm );
            read_array_parallel<int>( group_ptr, "cell_to_mesh_mapping_values", mapping_values, comm );

            cell_to_mesh_mapping_.clear ( );
            for ( int j = 0; j < mapping_keys.size ( ); ++j )
            {
                cell_to_mesh_mapping_.insert ( std::make_pair ( mapping_keys[j], mapping_values[j] ) );
            }

            // read individual values
            int dummy;
            int* buffer;
            read_array_parallel<int>( group_ptr, "layer_width", buffer, dummy, comm );
            layer_width_ = *buffer;
            delete buffer;
            assert(dummy == 1);

            read_array_parallel<int>( group_ptr, "quadrant_arrays_sorted", buffer, dummy, comm );
            quadrant_arrays_sorted_ = (*buffer != 0);
            delete buffer;
            assert(dummy == 1);
#else
            throw "HiFlow was not compiled with HDF5 support!\n";
#endif
        }

        void MeshPXestDatabase::copy_from ( const MeshPXestDatabasePtr db )
        {
            const int tdim = this->tdim ( );

            // copy content of base class and cast pointer
            MeshDatabase::copy_from ( db );

#ifdef WITH_P4EST
            //boost::shared_ptr<MeshPXestDatabase> new_db_pXest = boost::static_pointer_cast<MeshPXestDatabase> ( new_db );

            // copy connectivity
            if ( tdim == 2 )
            {
                assert ( db->p4est_conn_ != NULL );

                int num_conn_vertices = db->p4est_conn_->num_vertices;
                int num_conn_cells = db->p4est_conn_->num_trees;

                std::vector<double> conn_vertices ( num_conn_vertices * 3 );
                for ( int i = 0; i < num_conn_vertices * 3; ++i )
                {
                    conn_vertices[i] = db->p4est_conn_->vertices[i];
                }

                std::vector<int> tree_to_vertices ( num_conn_cells * 4 );
                for ( int i = 0; i < num_conn_cells * 4; ++i )
                {
                    tree_to_vertices[i] = db->p4est_conn_->tree_to_vertex[i];
                }

                p4est_connectivity_t* conn;
                pXest_build_conn ( num_conn_vertices,
                                   num_conn_cells,
                                   conn_vertices,
                                   tree_to_vertices,
                                   conn );
                this->p4est_conn_ = conn;
            }

            if ( tdim == 3 )
            {
                assert ( db->p8est_conn_ != NULL );

                int num_conn_vertices = db->p8est_conn_->num_vertices;
                int num_conn_cells = db->p8est_conn_->num_trees;

                std::vector<double> conn_vertices ( num_conn_vertices * 3 );
                for ( int i = 0; i < num_conn_vertices * 3; ++i )
                {
                    conn_vertices[i] = db->p8est_conn_->vertices[i];
                }

                std::vector<int> tree_to_vertices ( num_conn_cells * 8 );
                for ( int i = 0; i < num_conn_cells * 8; ++i )
                {
                    tree_to_vertices[i] = db->p8est_conn_->tree_to_vertex[i];
                }

                p8est_connectivity_t* conn;
                pXest_build_conn ( num_conn_vertices,
                                   num_conn_cells,
                                   conn_vertices,
                                   tree_to_vertices,
                                   conn );
                this->p8est_conn_ = conn;
            }

            // copy forest and ghost
            this->layer_width_ = db->layer_width_;

            if ( tdim == 2 )
            {
                assert ( db->get_p4est_forest ( ) != NULL );
                p4est_t* forest = p4est_copy ( db->p4est_forest_, 1 );
                forest->connectivity = this->p4est_conn_;

                p4est_ghost_t* ghost;
                pXest_copy_ghost ( db->p4est_ghost_, ghost );
                //assert ( p4est_ghost_is_valid( forest, ghost ) );

                this->p4est_forest_ = forest;
                this->p4est_ghost_ = ghost;
            }
            if ( tdim == 3 )
            {
                assert ( db->get_p8est_forest ( ) != NULL );
                p8est_t* forest = p8est_copy ( db->p8est_forest_, 1 );
                forest->connectivity = this->p8est_conn_;

                p8est_ghost_t* ghost;
                pXest_copy_ghost ( db->p8est_ghost_, ghost );
                //assert ( p4est_ghost_is_valid( forest, ghost ) );

                this->p8est_forest_ = forest;
                this->p8est_ghost_ = ghost;
            }

            // copy mesh history stuff
            this->cell_to_mesh_mapping_ = db->cell_to_mesh_mapping_;
            this->mesh_history_ = db->mesh_history_;

            // copy entity to quad mappings
            this->entity_to_quad_ = db->entity_to_quad_;
            this->coarse_entity_to_quad_ = db->coarse_entity_to_quad_;
#endif
        }

        void MeshPXestDatabase::deep_copy_from ( const MeshPXestDatabasePtr db )
        {
            // standard copy
            this->copy_from ( db );

            MeshPXestDatabasePtr self_ptr ( this );

            // make copies of meshes
            for ( int jm = 0; jm < db->mesh_history_.size ( ); ++jm )
            {
                MeshPtr tmp_mesh ( new MeshPXest ( this->tdim ( ), this->gdim ( ) ) );
                tmp_mesh->copy_from ( db->mesh_history_[jm] );

                boost::intrusive_ptr<MeshPXest> mesh_pXest = boost::static_pointer_cast<MeshPXest> ( tmp_mesh );
                mesh_pXest->set_db ( self_ptr );

                this->mesh_history_[jm] = tmp_mesh;
            }
        }
    } // end namespace mesh
} // namespace hiflow
