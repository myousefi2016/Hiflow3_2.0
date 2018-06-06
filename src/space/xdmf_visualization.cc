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

#include "xdmf_visualization.h"

#include "common/log.h"
#include "fem/cell_transformation.h"

#include "mesh/entity.h"
#include "mesh/iterator.h"
#include "mesh/attributes.h"

const int DEBUG_LEVEL = 1;

const std::string DEFAULT_SOL_NAME = "_sol_";

namespace hiflow
{

    template<class DataType>
    XdmfVisualization<DataType>::XdmfVisualization ( MPI_Comm& comm, std::string filename, bool append ) :
    comm_ ( comm ), rank_ ( -1 ), num_part_ ( -1 ), xdmf_file_ ( "default" ), current_iteration_ ( -1 ), current_grid_ ( "NO_GRID" )
    {
        {
            //get location
            size_t last_slash_idx = filename.rfind ( '/' );
            if ( std::string::npos != last_slash_idx )
            {
                location_ = filename.substr ( 0, last_slash_idx + 1 );
            }
            LOG_DEBUG ( 3, "Location: " << location_ );

            //get filename
            // Remove directory if present.
            // Do this before extension removal incase directory has a period character.
            last_slash_idx = filename.find_last_of ( "/" );
            if ( std::string::npos != last_slash_idx )
            {
                filename.erase ( 0, last_slash_idx + 1 );
            }

            // Remove extension if present.
            const size_t period_idx = filename.rfind ( '.' );
            if ( std::string::npos != period_idx )
            {
                filename.erase ( period_idx );
            }
            filename_ = filename;
        }
        LOG_DEBUG ( 3, "Filename: " << filename_ );
        // open hdf5 file

        if ( append == true )
        {
            continue_from_xdmf ( );
        }
        else
        {
            current_iteration_ = -1;

            //header of the xdmf file
            TiXmlDeclaration decl ( "1.0", "", "" );

            xdmf_file_.InsertEndChild ( decl );

            TiXmlElement root ( "Xdmf" );
            root.SetAttribute ( "xmlns:xi", "http://www.w3.org/2001/XInclude" );
            root.SetAttribute ( "Version", "2.0" );

            TiXmlElement domain ( "Domain" );

            TiXmlElement grid_coll ( "Grid" );

            grid_coll.SetAttribute ( "GridType", "Collection" );
            grid_coll.SetAttribute ( "CollectionType", "temporal" );

            domain.InsertEndChild ( grid_coll );

            root.InsertEndChild ( domain );
            xdmf_file_.InsertEndChild ( root );
        }

        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_part_ );

    }

    template<class DataType>
    void XdmfVisualization<DataType>::open_data_file ( )
    {
        std::string complete_path = location_;
        complete_path += filename_;
        complete_path += ".h5";
        LOG_DEBUG ( 3, "Open file " << complete_path );
#ifdef WITH_HDF5
        file_ptr_.reset ( new H5File ( complete_path, "w", comm_ ) );
#endif
    }

    template<class DataType>
    void XdmfVisualization<DataType>::close_data_file ( )
    {
        LOG_DEBUG ( 3, "Close file." );
#ifdef WITH_HDF5
        file_ptr_.reset ( );
#endif
    }

    template<class DataType>
    void XdmfVisualization<DataType>::add_timestep ( double time )
    {
        current_iteration_ += 1;
        LOG_DEBUG ( 3, "Iteration: " << current_iteration_ << " Time: " << time );
        TiXmlNode* grid_coll = xdmf_file_.LastChild ( )->FirstChild ( )->FirstChild ( );

        // Remove all later and current timesteps
        if ( grid_coll->LastChild ( ) != NULL )
        {
            while ( atof ( grid_coll->LastChild ( )->FirstChildElement ( )->Attribute ( "Value" ) ) >= time )
            {
                grid_coll->RemoveChild ( grid_coll->LastChild ( ) );
            }
        }

        // this is the element that will be attached to grid_coll
        TiXmlElement new_grid ( "Grid" );

        //the rest is just filling up new_grid with the correct data!
        new_grid.SetAttribute ( "Type", "Uniform" );
        new_grid.SetAttribute ( "GridType", "Collection" );
        new_grid.SetAttribute ( "CollectionType", "spatial" );

        TiXmlElement time_data ( "Time" );

        time_data.SetAttribute ( "Type", "Single" );
        time_data.SetDoubleAttribute ( "Value", time );

        new_grid.InsertEndChild ( time_data );

        grid_coll->InsertEndChild ( new_grid );
    }

    template<class DataType>
    void XdmfVisualization<DataType>::write_coupled_vector ( la::CoupledVector<DataType>& cd_vr, std::string identifier )
    {
#ifdef WITH_HDF5
        if ( current_iteration_ == -1 )
        {
            LOG_ERROR ( "Use the function add_timestep first (even in stationary case)" );
            exit ( -1 );
            return;
        }
        if ( identifier == "" )
        {
            std::stringstream ss;
            ss << DEFAULT_SOL_NAME;
            identifier = ss.str ( );
        }

        std::string complete_path = location_;
        complete_path += filename_;
        complete_path += ".h5";
        std::stringstream groupname;
        groupname << "time_step" << current_iteration_;
        cd_vr.WriteHDF5 ( complete_path, groupname.str ( ), identifier );
#else
        ERROR;
#endif //WITH_HDF5
    }

    template<class DataType>
    void XdmfVisualization<DataType>::write_hypre_vector ( la::HypreVector<DataType>& cd_vr, std::string identifier )
    {
#ifdef WITH_HDF5
        if ( current_iteration_ == -1 )
        {
            LOG_ERROR ( "Use the function add_timestep first (even in stationary case)" );
            exit ( -1 );
            return;
        }
        if ( identifier == "" )
        {
            std::stringstream ss;
            ss << DEFAULT_SOL_NAME;
            identifier = ss.str ( );
        }

        std::string complete_path = location_;
        complete_path += filename_;
        complete_path += ".h5";
        std::stringstream groupname;
        groupname << "time_step" << current_iteration_;
        cd_vr.WriteHDF5 ( complete_path, groupname.str ( ), identifier );
#else
        ERROR;
#endif //WITH_HDF5
    }

    template<class DataType>
    void XdmfVisualization<DataType>::write_grid ( std::vector< VectorSpace<DataType>* > spaces, int num_intervals, std::string grid_name, std::vector<int> variables, std::vector<DataType> ( *geom_map )( std::vector<DataType> coord ) )
    {
#ifdef WITH_HDF5

        // TOL used to determine whether two points are coinciding
        DataType TOL = 1.0e-6;
        if ( typeid (DataType ) == typeid (double ) )
        {
            TOL = 1.0e-13;
        }

        //assuming all meshes are exactly the same (TODO: make consistency tests)
        const mesh::Mesh& local_mesh = spaces[0]->mesh ( );
        const mesh::TDim tdim = local_mesh.tdim ( );
        const mesh::GDim gdim = local_mesh.gdim ( );
        CellVisualizationGrids<DataType>* grids = new CellVisualizationGrids<DataType>( &local_mesh, num_intervals, 0., 1. );

        std::vector< DataType> coords ( 0 );
        coords.reserve ( gdim * grids->num_visu_points ( ) );
        std::vector<DataType> ref_pt ( gdim, 0. ), mapped_pt ( 3, 0. );

        int num_variables = variables.size ( );
        std::vector<std::vector<DofID > > dof_ids ( num_variables );
        for ( int v = 0; v < num_variables; v++ )
        {
            dof_ids[v].reserve ( grids->num_visu_points ( ) );
        }

        for ( mesh::EntityIterator it = local_mesh.begin ( tdim ), end_it = local_mesh.end ( tdim );
              it != end_it; ++it )
        {
            int remote_index;
            local_mesh.get_attribute_value ( "_remote_index_", local_mesh.tdim ( ),
                                             it->index ( ), &remote_index );
            // we skip the ghost cells.
            if ( remote_index != -1 ) continue;
            std::vector< std::vector< std::vector<DataType> > > cell_coords ( num_variables );
            std::vector< std::vector<DofID> > cell_dofs ( num_variables );

            for ( int v = 0; v < num_variables; v++ )
            {
                spaces[v]->dof ( ).get_coord_on_cell ( variables[v], it->index ( ), cell_coords[v] );
                spaces[v]->dof ( ).get_dofs_on_cell ( variables[v], it->index ( ), cell_dofs[v] );
            }
            // like cell_visualization:
            const doffem::CellTransformation<DataType>& cell_trans = spaces[0]->GetCellTransformation ( *it );

            // TODO: Presorting for different num_intervals!!
            // list of dofs we want to search for every variable
            std::vector< std::list<int> > search_lists ( num_variables );
            for ( int v = 0; v < num_variables; v++ )
            {
                for ( int i = 0; i < cell_dofs[v].size ( ); i++ )
                {
                    search_lists[v].push_back ( i );
                }
            }

            // TODO: this can be done with a single call now
            //      for (int p = 0, p_end = grid_.num_points(); p != p_end; ++p) {
            for ( int p = 0, p_end = grids->num_points ( it->cell_type ( ).tag ( ) ); p != p_end; ++p )
            {
                const int offset = gdim * p;
                for ( int c = 0; c < gdim; ++c )
                {
                    ref_pt[c] = grids->coords ( it->cell_type ( ).tag ( ) )[offset + c];
                    mapped_pt[c] = 0.; // reset
                }

                // TODO: Would be great if we could do this with one call, as
                // with Hp-transformations.
                if ( gdim >= 1 ) mapped_pt[0] = cell_trans.x ( ref_pt );
                if ( gdim >= 2 ) mapped_pt[1] = cell_trans.y ( ref_pt );
                if ( gdim >= 3 ) mapped_pt[2] = cell_trans.z ( ref_pt );

                //apply geometrical mapping
                std::vector<DataType> geo_mapped_pt = geom_map ( mapped_pt );
                coords.insert ( coords.end ( ), geo_mapped_pt.begin ( ), geo_mapped_pt.begin ( ) + gdim );

                //for each variable find a dof lying on the current phys coord
                //(TODO: Maybe we can use reference coords here)
                for ( int v = 0; v < num_variables; v++ )
                {
                    // flag to break the for-loop when the right DoF was found
                    bool found_dof = false;
                    // search only through the cell_coords that hasn't been associated
                    // with a grid point yet. This only works if there aren't two
                    // DoFs of one cell on one coordinate!
                    std::list<int>::iterator list_iter = search_lists[v].begin ( );
                    while ( list_iter != search_lists[v].end ( ) )
                    {
                        int i = *list_iter;
                        for ( int d = 0; d < gdim; ++d )
                        {
                            if ( std::abs ( mapped_pt[d] - cell_coords[v][i][d] ) > TOL )
                            {
                                list_iter++;
                                break;
                            }
                            if ( d == gdim - 1 )
                            {
                                // Delete the coord index here, so we don't search it twice
                                list_iter = search_lists[v].erase ( list_iter );
                                dof_ids[v].push_back ( cell_dofs[v][i] );
                                found_dof = true;
                            }
                        }
                        if ( found_dof ) break;
                    }
                    if ( !found_dof )
                    {
                        LOG_ERROR ( "num_intervals " << num_intervals << " not compatible with variable " << v << "!" );
                        exit ( -1 );
                        return;
                    }
                }
            }
        }

        H5GroupPtr group_ptr ( new H5Group ( file_ptr_, grid_name, "w" ) );

        //getting offset and total size data to correctly write down the coords and the "maps"
        int coords_num = grids->num_visu_points ( );
        assert ( coords.size ( ) == coords_num * gdim );
        int coords_offset = 0;
        if ( rank_ > 0 )
        {
            MPI_Status status;
            MPI_Recv ( &coords_offset, 1, MPI_INT, rank_ - 1, rank_ - 1, comm_, &status );
        }
        int coords_next_offset = coords_offset + coords_num;
        if ( rank_ < num_part_ - 1 )
        {
            MPI_Send ( &coords_next_offset, 1, MPI_INT, ( rank_ + 1 ), rank_, comm_ );
        }

        MPI_Bcast ( &coords_next_offset, 1, MPI_INT, num_part_ - 1, comm_ );
        // printing coords
        H5DatasetPtr dataset_ptr_crds ( new H5Dataset ( group_ptr, coords_next_offset * gdim,
                                                        "coordinates", "w", &coords[0] ) );
        dataset_ptr_crds->write ( coords_num * gdim, coords_offset* gdim, &coords[0] );
        coords.clear ( );
        //printing maps
        for ( int v = 0; v < variables.size ( ); v++ )
        {
            std::stringstream map_name;
            map_name << grid_name << v;
            H5DatasetPtr dataset_ptr_map ( new H5Dataset ( group_ptr, coords_next_offset,
                                                           map_name.str ( ), "w", &dof_ids[v][0] ) );
            dataset_ptr_map->write ( coords_num, coords_offset, &dof_ids[v][0] );
            dof_ids[v].clear ( );
        }
        std::vector<int> incidents ( 0 );

        int p_offset = 0;
        int pt_offset = coords_offset;

        // Connectivity
        static const int xdmf_cell_types[] = { 1, 2, 4, 5, 6, 9 }; //not sure about the fist two

        for ( mesh::EntityIterator it = local_mesh.begin ( tdim ); it != local_mesh.end ( tdim ); ++it )
        {
            int rem_ind = -100;
            it->get<int>( "_remote_index_", &rem_ind );
            if ( rem_ind != -1 ) continue;

            for ( int c = 0, c_end = grids->num_cells ( it->cell_type ( ).tag ( ) ); c != c_end; ++c )
            {
                const std::vector<int>& verts = grids->vertices_of_cell ( it->cell_type ( ).tag ( ), c );
                // getting cell type id for xdmf
                incidents.push_back ( xdmf_cell_types[static_cast < int > ( it->cell_type ( ).tag ( ) )] );
                for ( int v = 0, v_end = verts.size ( ); v != v_end; ++v )
                {
                    incidents.push_back ( verts[v] + p_offset + pt_offset );
                }
            }
            p_offset += grids->num_points ( it->cell_type ( ).tag ( ) );
        }
        // some communication to get the total length of the incid(ents) array
        int inc_size = incidents.size ( );
        int inc_offset = 0;
        if ( rank_ > 0 )
        {
            MPI_Status status;
            MPI_Recv ( &inc_offset, 1, MPI_INT, rank_ - 1, rank_ - 1, comm_, &status );
        }
        int inc_next_offset = inc_offset + inc_size;
        if ( rank_ < num_part_ - 1 )
        {
            MPI_Send ( &inc_next_offset, 1, MPI_INT, ( rank_ + 1 ), rank_, comm_ );
        }
        MPI_Bcast ( &inc_next_offset, 1, MPI_INT, num_part_ - 1, comm_ );

        // once communication is done we can write!
        H5DatasetPtr dataset_ptr_inc ( new H5Dataset ( group_ptr, inc_next_offset,
                                                       "incidents", "w", &incidents[0] ) );
        dataset_ptr_inc->write ( inc_size, inc_offset, &incidents[0] );

        incidents.clear ( );
        // more communication needed to get the number of entities in total!
        int num_ent = grids->num_visu_cells ( );
        int num_ent_offset = 0;
        if ( rank_ > 0 )
        {
            MPI_Status status;
            MPI_Recv ( &num_ent_offset, 1, MPI_INT, rank_ - 1, rank_ - 1, comm_, &status );
        }
        int num_ent_next_offset = num_ent_offset + num_ent;
        if ( rank_ < num_part_ - 1 )
        {
            MPI_Send ( &num_ent_next_offset, 1, MPI_INT, ( rank_ + 1 ), rank_, comm_ );
        }
        MPI_Bcast ( &num_ent_next_offset, 1, MPI_INT, num_part_ - 1, comm_ );

        GridData* local_grid = new GridData;
        local_grid->num_visu_cells_ = num_ent_next_offset;
        local_grid->num_incidents_ = inc_next_offset;
        local_grid->num_visu_points_ = coords_next_offset;
        local_grid->gdim_ = gdim;
        local_grid->name_ = grid_name;

        grid_info_[grid_name] = *local_grid;

        delete grids;
#else
        ERROR;
#endif //WITH_HDF5
    }

    template<class DataType>
    void XdmfVisualization<DataType>::add_xdmfattribute ( std::string att_type, std::string cv_id, std::string func_name, std::vector<int> variables, std::string grid_name, int num_dofs )
    {
        // writing the xmf file: we do this for all processes -> this might
        // be changed so that only the Master_rank does this
        if ( cv_id == "" )
        {
            std::stringstream ss;
            ss << DEFAULT_SOL_NAME;
            cv_id = ss.str ( );
        }
        //this is a pointer to the grid collection we want to write in.
        TiXmlNode* grid_coll = xdmf_file_.LastChild ( )->FirstChild ( )->FirstChild ( )->LastChild ( );

        TiXmlElement* grid = dynamic_cast < TiXmlElement* > ( grid_coll->FirstChild ( "Grid" ) );

        //we check if this grid exists already in the xdmf or has yet to be printed
        //maybe we can instead use a flag as member variable to check this
        bool grid_exists = false;
        while ( grid != NULL )
        {
            if ( grid->Attribute ( "Name" ) == grid_name )
            {
                grid_exists = true;
                break;
            }
            grid = dynamic_cast < TiXmlElement* > ( grid_coll->IterateChildren ( "Grid", grid ) );
        }
        if ( !grid_exists )
        {
            this->add_xdmfgrid ( grid_name );
            this->add_xdmfattribute ( att_type, cv_id, func_name, variables, grid_name, num_dofs );
            return;
        }
        //hdf5 file:
        std::string h5_filename = filename_;
        h5_filename += ".h5";
        //special treatment for 2d vectors (set 3rd component to zero)
        bool is_2d_vector = false;
        if ( att_type == "2dVector" )
        {
            is_2d_vector = true;
            att_type = "Vector";
        }

        GridData* local_grid = &( grid_info_[grid_name] ); //get_grid_info(grid_name);

        int map_total_size = local_grid->num_visu_points_;

        //collect dataitems for the XdmfFunction
        std::vector<XdmfBaseItem*> collected_vars;
        for ( int v = 0; v < variables.size ( ); v++ )
        {

            std::stringstream time_group;
            time_group << grid_name;
            std::stringstream map_location_name;
            map_location_name << h5_filename << ":/" << time_group.str ( ) << "/" << grid_name << variables[v];

            XdmfDataItem* map_data = new XdmfDataItem ( map_location_name.str ( ), map_total_size, "Int" );

            std::stringstream func_location_name;

            func_location_name << h5_filename << ":/time_step" << current_iteration_ << "/" << cv_id;

            // TODO get right number of dofs
            XdmfDataItem* sol_data = new XdmfDataItem ( func_location_name.str ( ), num_dofs, "Float", 8 );
            XdmfCoordinate* map_coord = new XdmfCoordinate ( map_data, sol_data );

            collected_vars.push_back ( map_coord );
        }
        if ( is_2d_vector )
        {

            std::stringstream time_group;
            time_group << grid_name;
            std::stringstream map_location_name;
            map_location_name << h5_filename << ":/" << time_group.str ( ) << "/" << "zerofunc";

            XdmfDataItem* map_data = new XdmfDataItem ( map_location_name.str ( ), map_total_size, "Float", 8 );
            collected_vars.push_back ( map_data );
        }
        XdmfFunction* vector_data = new XdmfFunction ( collected_vars, "JOIN" );
        XdmfAttribute att_map ( vector_data, func_name, att_type, "Node" );
        grid->LinkEndChild ( att_map.get_xdmf_element ( ) );
        att_map.clean_deletion ( );
    }

    template<class DataType>
    void XdmfVisualization<DataType>::add_xdmfgrid ( std::string grid_name )
    {
        std::string h5_filename = filename_;
        h5_filename += ".h5";
        GridData* grid_info = &( grid_info_[grid_name] );
        TiXmlNode* grid_coll = xdmf_file_.LastChild ( )->FirstChild ( )->FirstChild ( )->LastChild ( );

        // this is the element that will be attached to grid_coll
        TiXmlElement* new_grid = new TiXmlElement ( "Grid" );

        //the rest is just filling up new_grid with the correct data!
        new_grid->SetAttribute ( "Name", grid_name );
        new_grid->SetAttribute ( "Type", "Uniform" );

        std::stringstream inc_location_name;
        inc_location_name << h5_filename << ":/" << grid_name << "/incidents";

        XdmfDataItem incidents ( inc_location_name.str ( ), grid_info->num_incidents_, "Int" );
        XdmfTopology topology ( &incidents, grid_info->num_visu_cells_ );

        new_grid->LinkEndChild ( topology.get_xdmf_element ( ) );

        mesh::GDim gdim = grid_info->gdim_;

        std::stringstream coords_location_name;
        coords_location_name << h5_filename << ":/" << grid_name << "/coordinates";
        XdmfDataItem coords ( coords_location_name.str ( ), grid_info->num_visu_points_ * gdim, "Float", 8 );
        XdmfGeometry geometry ( &coords, gdim );

        new_grid->LinkEndChild ( geometry.get_xdmf_element ( ) );

        grid_coll->LinkEndChild ( new_grid );
    }

    template<class DataType>
    void XdmfVisualization<DataType>::add_zero_function ( std::string grid_name )
    {
#ifdef WITH_HDF5
        bool* has_func = &( grid_info_[grid_name].has_zero_func_ );
        //return early if we already have a null function for this grid
        if ( *has_func ) return;
        this->open_data_file ( );
        //get data about the size of this grid
        int num_pts = grid_info_[grid_name].num_visu_points_;
        int rest = num_pts % num_part_;
        int lcl_offset = 0;
        int lcl_num_pts = num_pts / num_part_;
        if ( rank_ == 0 )
        {
            lcl_num_pts += rest;
        }
        else
        {
            lcl_offset = rank_ * lcl_num_pts + rest;
        }
        /*int pt_offset = 0;
        if( rank_ > 0) {
          MPI_Status status;
          MPI_Recv(&pt_offset, 1, MPI_INT, rank_ - 1,  rank_ - 1, comm_, &status);
        }
        int pt_next_offset = pt_offset + pt_size;
        if( rank_ < num_part_ - 1) {
          MPI_Send(&pt_next_offset, 1, MPI_INT, (rank_ + 1), rank_, comm_);
        }

        MPI_Bcast(&pt_next_offset, 1, MPI_INT, num_part_ - 1, comm_);*/

        //write zero vector to the hdf5 file
        std::vector<DataType> zero_func ( lcl_num_pts, static_cast < DataType > ( 0. ) );

        H5GroupPtr group_ptr ( new H5Group ( file_ptr_, grid_name, "w" ) );
        std::stringstream func_name;
        func_name << "zerofunc";
        H5DatasetPtr dataset_ptr_map ( new H5Dataset ( group_ptr, num_pts,
                                                       func_name.str ( ), "w", &zero_func[0] ) );
        dataset_ptr_map->write ( lcl_num_pts, lcl_offset, &zero_func[0] );
        zero_func.clear ( );

        *has_func = true;
        this->close_data_file ( );
#else
        ERROR;
#endif
    }

    template<class DataType>
    void XdmfVisualization<DataType>::add_to_view ( VectorSpace<DataType>* space, int variable, std::string id, std::string func_name )
    {
        std::vector<int> variables ( 1, variable );
        this->add_to_view ( space, variables, id, func_name );
    }

    template<class DataType>
    void XdmfVisualization<DataType>::add_to_view ( VectorSpace<DataType>* space, std::vector<int> variables, std::string id, std::string func_name )
    {
        new_spaces_.push_back ( space );
        new_variables_.push_back ( variables );
        new_vec_ids_.push_back ( id );
        new_func_names_.push_back ( func_name );
    }

    template<class DataType>
    void XdmfVisualization<DataType>::write_view ( int num_intervals, std::string grid_name, std::vector<DataType> ( *geom_map )( std::vector<DataType> coord ) )
    {
        current_grid_ = grid_name;
        std::vector<int> variables ( 0 );
        std::vector< VectorSpace<DataType>* > spaces ( 0 );
        for ( int i = 0, n = new_variables_.size ( ); i < n; ++i )
        {
            for ( int j = 0, m = new_variables_[i].size ( ); j < m; ++j )
            {
                variables.push_back ( new_variables_[i][j] );
                spaces.push_back ( new_spaces_[i] );
            }
        }
        this->open_data_file ( );
        this->write_grid ( spaces, num_intervals, current_grid_, variables, geom_map );
        this->close_data_file ( );
        variables_ = new_variables_;
        func_names_ = new_func_names_;
        vec_ids_ = new_vec_ids_;

        new_spaces_.clear ( );
        new_variables_.clear ( );
        new_func_names_.clear ( );
        new_vec_ids_.clear ( );
    }

    template<class DataType>
    void XdmfVisualization<DataType>::write_solution_vectors ( std::vector< la::CoupledVector<DataType>* > solutions, std::vector< std::string> vec_ids )
    {

        assert ( solutions.size ( ) == vec_ids.size ( ) );
        std::map<std::string, la::CoupledVector<DataType>* > id_to_vec;
        for ( int i = 0, n = solutions.size ( ); i < n; ++i )
        {
            write_coupled_vector ( *( solutions[i] ), vec_ids[i] );
            id_to_vec[vec_ids[i]] = solutions[i];
        }
        int var_cntr = 0;
        std::vector<int> visu_vars ( 0 );
        for ( int i = 0, n = func_names_.size ( ); i < n; ++i )
        {
            std::string func_type;
            switch ( variables_[i].size ( ) )
            {
                case 1:
                    func_type = "Scalar";
                    break;
                case 2:
                    func_type = "2dVector";
                    this->add_zero_function ( current_grid_ );
                    break;
                case 3:
                    func_type = "Vector";
                    break;
                case 6:
                    func_type = "Tensor6";
                    break;
                case 9:
                    func_type = "Tensor";
                    break;
                default:
                    LOG_ERROR ( "Number of variables must be either 1 (scalar), 2 (2d vector), 3 (3d vector), 6 (symmetrical tensor) or 9  (full tensor)" );
                    exit ( -1 );
            }
            visu_vars.resize ( variables_[i].size ( ) );
            for ( int j = 0, m = variables_[i].size ( ); j < m; ++j )
            {
                visu_vars[j] = var_cntr;
                ++var_cntr;
            }

            this->add_xdmfattribute ( func_type, vec_ids_[i], func_names_[i], visu_vars, current_grid_, id_to_vec[vec_ids_[i]]->size_global ( ) );
        }
    }

    template<class DataType>
    void XdmfVisualization<DataType>::write_solution_vectors ( std::vector< la::HypreVector<DataType>* > solutions, std::vector< std::string> vec_ids )
    {

        assert ( solutions.size ( ) == vec_ids.size ( ) );
        std::map<std::string, la::HypreVector<DataType>* > id_to_vec;
        for ( int i = 0, n = solutions.size ( ); i < n; ++i )
        {
            write_hypre_vector ( *( solutions[i] ), vec_ids[i] );
            id_to_vec[vec_ids[i]] = solutions[i];
        }
        int var_cntr = 0;
        std::vector<int> visu_vars ( 0 );
        for ( int i = 0, n = func_names_.size ( ); i < n; ++i )
        {
            std::string func_type;
            switch ( variables_[i].size ( ) )
            {
                case 1:
                    func_type = "Scalar";
                    break;
                case 2:
                    func_type = "2dVector";
                    this->add_zero_function ( current_grid_ );
                    break;
                case 3:
                    func_type = "Vector";
                    break;
                case 6:
                    func_type = "Tensor6";
                    break;
                case 9:
                    func_type = "Tensor";
                    break;
                default:
                    LOG_ERROR ( "Number of variables must be either 1 (scalar), 2 (2d vector), 3 (3d vector), 6 (symmetrical tensor) or 9  (full tensor)" );
                    exit ( -1 );
            }
            visu_vars.resize ( variables_[i].size ( ) );
            for ( int j = 0, m = variables_[i].size ( ); j < m; ++j )
            {
                visu_vars[j] = var_cntr;
                ++var_cntr;
            }

            this->add_xdmfattribute ( func_type, vec_ids_[i], func_names_[i], visu_vars, current_grid_, id_to_vec[vec_ids_[i]]->size_global ( ) );
        }
    }

    template<class DataType>
    void XdmfVisualization<DataType>::write_solution_vector ( la::CoupledVector<DataType>* solution, std::string vec_id )
    {
        std::vector< la::CoupledVector<DataType>* > solutions ( 1, solution );
        std::vector< std::string > vec_ids ( 1, vec_id );
        this->write_solution_vectors ( solutions, vec_ids );
    }

    template<class DataType>
    void XdmfVisualization<DataType>::write_solution_vector ( la::HypreVector<DataType>* solution, std::string vec_id )
    {
        std::vector< la::HypreVector<DataType>* > solutions ( 1, solution );
        std::vector< std::string > vec_ids ( 1, vec_id );
        this->write_solution_vectors ( solutions, vec_ids );
    }

    template<class DataType>
    void XdmfVisualization<DataType>::continue_from_xdmf ( )
    {
        std::string complete_path = location_;
        complete_path += filename_;
        complete_path += ".xmf";
        bool status = xdmf_file_.LoadFile ( complete_path );
        if ( !status )
        {
            std::cerr << "Loading file failed!" << std::endl;
        }
        //load data from the xmf file
        // infos we need:
        //  --- current iteration
        //  --- all grid data
        //  --- informations about the functions that are visualized (attributes in xdmf)
        TiXmlHandle docHandle ( &xdmf_file_ );
        TiXmlHandle grid_time_coll = docHandle.FirstChild ( "Xdmf" ).FirstChild ( "Domain" ).FirstChild ( "Grid" );
        TiXmlElement* grid_spat_coll = grid_time_coll.FirstChild ( "Grid" ).ToElement ( );
        //get current iteration by counting the number of children
        while ( grid_spat_coll != NULL )
        {
            current_iteration_ += 1;
            grid_spat_coll = grid_spat_coll->NextSiblingElement ( );
        }
        //assert that it has been done one iteration before
        assert ( current_iteration_ > -1 );
        //now get the attributes from the xdmf file (of the last grid only!)
        //the only variable we don't get is "has_zero_func"
        TiXmlHandle last_grid = grid_time_coll.Child ( current_iteration_ ).FirstChild ( "Grid" );

        std::string grid_name = last_grid.ToElement ( )->Attribute ( "Name" );

        GridData* new_grid = new GridData;

        new_grid->name_ = grid_name;

        new_grid->num_visu_cells_ = std::atoi ( last_grid.FirstChild ( "Topology" ).ToElement ( )->Attribute ( "NumberOfElements" ) );
        new_grid->num_incidents_ = std::atoi ( last_grid.FirstChild ( "Topology" ).FirstChild ( "DataItem" ).ToElement ( )->Attribute ( "Dimensions" ) );
        std::string geo_dim = last_grid.FirstChild ( "Geometry" ).ToElement ( )->Attribute ( "Type" );
        if ( geo_dim == "XYZ" )
        {
            new_grid->gdim_ = 3;
        }
        else if ( geo_dim == "XY" )
        {
            new_grid->gdim_ = 2;
        }
        else
        {
            std::cerr << "Unkown Geometry Type" << std::endl;
        }

        new_grid->num_visu_points_ = std::atoi ( last_grid.FirstChild ( "Geometry" ).FirstChild ( "DataItem" ).ToElement ( )->Attribute ( "Dimensions" ) ) / new_grid->gdim_;
        new_grid->has_zero_func_ = false;

        //now get data about the functions to be visualized!

        TiXmlElement* attributes = last_grid.FirstChild ( "Attribute" ).ToElement ( );
        int node_counter = 0;
        while ( attributes != NULL )
        {
            TiXmlHandle att_handle ( attributes );
            std::string att_center = attributes->Attribute ( "Center" );
            if ( att_center == "Node" )
            {
                int num_vars = 0;
                std::string att_type = attributes->Attribute ( "AttributeType" );
                if ( att_type == "Scalar" )
                {
                    num_vars = 1;
                }
                else if ( att_type == "Vector" )
                {
                    if ( att_handle.FirstChild ( "DataItem" ).ToElement ( )->LastChild ( "DataItem" )->ToElement ( )->Attribute ( "ItemType" ) == NULL )
                    {
                        att_type = "2dVector";
                        new_grid->has_zero_func_ = true;
                        num_vars = 2;
                    }
                    else
                    {
                        num_vars = 3;
                    }
                }
                else if ( att_type == "Tensor6" )
                {
                    num_vars = 6;
                }
                else if ( att_type == "Tensor" )
                {
                    num_vars = 9;
                }
                std::vector<int> temp_vector ( 0 );
                for ( int i = 0; i < num_vars; ++i )
                {
                    temp_vector.push_back ( node_counter );
                    node_counter++;
                }
                variables_.push_back ( temp_vector );
                std::string att_name = attributes->Attribute ( "Name" );
                func_names_.push_back ( att_name );
                std::string location_text = att_handle.FirstChild ( "DataItem" ).FirstChild ( "DataItem" ).ChildElement ( "DataItem", 1 ).ToElement ( )->GetText ( );

                size_t last_slash_idx = location_text.find_last_of ( "/" );
                if ( std::string::npos != last_slash_idx )
                {
                    location_text.erase ( 0, last_slash_idx + 1 );
                }
                vec_ids_.push_back ( location_text );
            }
            attributes = attributes->NextSiblingElement ( "Attribute" );
        }
        grid_info_[grid_name] = *new_grid;
        current_grid_ = grid_name;

    }
    template class XdmfVisualization<double>;
    //template class XdmfVisualization<float>;
}
