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

#ifndef HIFLOW_SPACE_XDMF_VISUALIZATION
#    define HIFLOW_SPACE_XDMF_VISUALIZATION

/// \author Jonathan Schwegler
///
/// \brief Visualization of finite element functions.
///
/// This class provides functions to visualize FeFunctions using only
/// values saved in the corresponding CoupledVector. The generated outputs
/// are a HDF5 file storing the "heavy data" and a XDMF file that can be
/// opened with e.g. ParaView (www.paraview.org/).
/// It is not (yet) possible to use this class to append the visualization
/// data to an existing one instead it will overwrite the previous.

#    include "config.h"

#    include "mesh/mesh.h"

#    include <sstream>
#    include <string>
#    include <vector>
#    include <map>
#    include <list>
#    include <mpi.h>
#    include "common/log.h"
#    include <fstream>
#    include <iostream>
#    include <sys/stat.h>

#    include "space/vector_space.h"
#    include "space/cell_visualization.h" //for CellVisualizationGrids

#    ifdef WITH_HDF5
#        include "common/hdf5_tools.h"
#    else
#        define ERROR LOG_ERROR("HiFlow was not compiled with HDF5 support!");  exit(-1);
#    endif

#    include "linear_algebra/coupled_vector.h"

#    define TIXML_USE_STL
#    include <tinyxml/tinyxml.h>

namespace hiflow
{
    struct GridData;

    template<class DataType>
    std::vector< DataType> identity ( std::vector< DataType> coord )
    {
        return coord;
    }

    template<class DataType>
    class XdmfVisualization
    {
        typedef hiflow::doffem::DofID DofID;
      public:

        XdmfVisualization ( MPI_Comm& comm, std::string filename, bool append = false );

        ~XdmfVisualization ( )
        {
        };

        /// \brief Add a scalar valued finite element function to the view
        /// \param [in] space Pointer to the space the function is from
        /// \param [in] variable Variable to be visualized
        /// \param [in] vec_id Identifier for the solution vector
        /// \param [in] func_name Name that will be displayed in ParaView
        void add_to_view ( VectorSpace<DataType>* space, int variable, std::string vec_id, std::string func_name );

        /// \brief Add a vector/tensor valued finited element function to the view.
        /// If used two or three variables -> 2d/3d Vector
        /// If used six or nine variables  -> (symmetrical) tensor
        /// Other numbers of variables will result in an error.
        void add_to_view ( VectorSpace<DataType>* space, std::vector<int> variables, std::string vec_id, std::string func_name );

        /// \brief Write the view added with add_to_view
        /// \param [in] num_intervals Resolution of the view:
        /// e.g. num_intervals=1 for Q1/P1-view, =2 for Q2/P2, etc...
        /// The degree of all variables added before must be divisible by
        /// num_intervals. So num_intervals=1 will always work.
        /// \param [in] grid_name Name of the grid. Must be unique during
        /// the simulation process if more than one grid is used.
        /// \param [in] geom_map Function that maps the coordinates of the
        /// mesh for the visualization. Default is identity.
        void write_view ( int num_intervals, std::string grid_name, std::vector<DataType> ( *geom_map )( std::vector<DataType> coord ) = identity<DataType> );

        /// \brief Add a new time step in a instationary simulation. You also
        /// have to use this function once in a stationary case.
        /// \param [in] time Current time in the simulation process.
        void add_timestep ( double time = 0. );

        /// \brief Writes the solution vectors that are to be visualized.
        /// the vec_ids must correspond to the ids given in add_to_view.
        /// \param [in] solutions Containing the vectors to be visualized.
        void write_solution_vectors ( std::vector< la::CoupledVector<DataType>* > solutions, std::vector< std::string> vec_ids );

        /// \brief Writes the solution vectors that are to be visualized.
        /// the vec_ids must correspond to the ids given in add_to_view.
        /// \param [in] solutions Containing the vectors to be visualized.
        void write_solution_vectors ( std::vector< la::HypreVector<DataType>* > solutions, std::vector< std::string> vec_ids );

        /// \brief Alternative to write_solution_vectors. Instead of collecting
        /// the CoupledVectors into a vector, you can call this function for
        /// each CoupledVector separatly.
        void write_solution_vector ( la::CoupledVector<DataType>* solution, std::string vec_id );

        /// \brief Alternative to write_solution_vectors. Instead of collecting
        /// the CoupledVectors into a vector, you can call this function for
        /// each CoupledVector separatly.
        void write_solution_vector ( la::HypreVector<DataType>* solution, std::string vec_id );

        /// \brief This function completes the visualization process and
        /// writes a XDMF file that can be opened with e.g ParaView
        /// \param [in] print_rank The rank that writes the xdmf-file

        void print_xdmf ( int print_rank = 0 )
        {
            assert ( print_rank < num_part_ );
            if ( rank_ == print_rank )
            {
                std::string complete_path = location_;
                complete_path += filename_;
                complete_path += ".xmf";
                xdmf_file_.SaveFile ( complete_path );
            }

            this->close_data_file ( );
        }

        /// \brief Overwrite automatically detected iteration number by hand
        /// \param[in] iteration number of current iteration

        void set_current_iteration ( int iteration )
        {
            current_iteration_ = iteration;
        }

      private:
        //open the hdf5 file
        void open_data_file ( );
        //close the hdf5 file
        void close_data_file ( );
        //add a attribute to the internal xdmf document
        void add_xdmfattribute ( std::string att_type, std::string cv_id, std::string func_name, std::vector<int> variables, std::string grid_name, int num_dofs );
        //add a grid to the internal xdmf document
        void add_xdmfgrid ( std::string grid_name );
        //write down the grid to the hdf5 file
        void write_grid ( std::vector<VectorSpace<DataType>* > spaces, int num_intervals, std::string grid_name, std::vector<int> variables, std::vector<DataType> ( *geom_map )( std::vector<DataType> coord ) );
        //write down a CoupledVector to the hdf5 file
        void write_coupled_vector ( la::CoupledVector<DataType>& cd_vr, std::string identifier );
        //write down a HypreVector to the hdf5 file
        void write_hypre_vector ( la::HypreVector<DataType>& cd_vr, std::string identifier );
        //add a zero function to the grid
        void add_zero_function ( std::string grid_name );
        void continue_from_xdmf ( );

        const MPI_Comm& comm_;

        //destiny of the output
        std::string location_;
        //filename of the xdmf and hdf5 file (without extension)
        std::string filename_;
#    ifdef WITH_HDF5
        H5FilePtr file_ptr_; //pointer to the hdf5 file
#    endif
        int rank_;
        int num_part_;
        TiXmlDocument xdmf_file_; //manages the creation of the xdmf file
        std::map< std::string, GridData > grid_info_; //Not sure if it is necessary to save the GridData for every Grid...

        int current_iteration_;

        //NEW INTERFACE VARIABLES
        std::vector<VectorSpace<DataType>* > new_spaces_;
        std::vector< std::vector<int> > new_variables_;
        std::vector<std::string> new_func_names_;
        std::vector<std::string> new_vec_ids_;

        std::vector< std::vector<int> > variables_;
        std::vector<std::string> func_names_;
        std::vector<std::string> vec_ids_;
        std::string current_grid_;

    };

    // holds meta-data about the grids

    struct GridData
    {

        GridData ( ) : name_ ( "" ), has_zero_func_ ( false ), num_visu_cells_ ( -1 ), num_visu_points_ ( -1 ), num_incidents_ ( -1 ), gdim_ ( -1 )
        {
        };

        ~GridData ( )
        {
        };
        std::string name_;

        bool has_zero_func_;
        // following 4 "int" are there to create a xdmf file suiting this grid!
        int num_visu_cells_;
        int num_visu_points_;

        int num_incidents_;

        int gdim_;
    };

    //Structs to create XDMF-conform XmlElements

    struct XdmfBaseItem
    {
      public:

        XdmfBaseItem ( )
        {
        };
        virtual TiXmlElement* get_xdmf_element ( ) = 0;
        virtual void clean_deletion ( ) = 0;
        virtual std::vector<int> dimensions ( ) = 0;

        std::string dim_to_string ( )
        {
            std::vector<int> temp_dim = dimensions ( );
            return string_from_range ( temp_dim.begin ( ), temp_dim.end ( ) );
        };

        virtual ~XdmfBaseItem ( )
        {
        };
    };

    //Xdmf DataItem containing hdf5-Data

    struct XdmfDataItem : XdmfBaseItem
    {
      public:

        XdmfDataItem ( std::string location, int size, std::string number_type, int precision = -1 )
        : dimensions_ ( size ), number_type_ ( number_type ), location_ ( location )
        {
            assert ( ( number_type_ == "Int" ) || ( number_type_ == "Float" ) );
            if ( number_type_ == "Float" )
            {
                assert ( precision > 0 );
                precision_ = precision;
            }
            else
            {
                precision_ = -1;
            }
        }

        TiXmlElement* get_xdmf_element ( )
        {
            TiXmlElement* dataitem = new TiXmlElement ( "DataItem" );

            dataitem->SetAttribute ( "Dimensions", dim_to_string ( ) );
            dataitem->SetAttribute ( "NumberType", number_type_ );
            if ( precision_ > 0 )
                dataitem->SetAttribute ( "Precision", precision_ );
            dataitem->SetAttribute ( "Format", "HDF" );
            TiXmlText* location = new TiXmlText ( location_ );
            dataitem->LinkEndChild ( location );
            return dataitem;
        }

        void clean_deletion ( )
        {
        }

        std::vector<int> dimensions ( )
        {
            return std::vector<int>( 1, dimensions_ );
        }

      private:
        int dimensions_; //size of Data
        std::string number_type_; //Float or Int
        std::string location_; //location of data
        int precision_;
    };

    struct XdmfCoordinate : XdmfBaseItem
    {

        XdmfCoordinate ( XdmfBaseItem* map, XdmfBaseItem* data ) : map_ ( map ), data_ ( data )
        {
        };

        TiXmlElement* get_xdmf_element ( )
        {
            assert ( map_ != NULL );
            assert ( data_ != NULL );
            TiXmlElement* coordinate = new TiXmlElement ( "DataItem" );

            coordinate->SetAttribute ( "ItemType", "coordinates" );
            coordinate->SetAttribute ( "Dimensions", dim_to_string ( ) );

            coordinate->LinkEndChild ( map_->get_xdmf_element ( ) );
            coordinate->LinkEndChild ( data_->get_xdmf_element ( ) );

            return coordinate;
        }

        void clean_deletion ( )
        {
            map_->clean_deletion ( );
            delete map_;
            data_->clean_deletion ( );
            delete data_;
        }

        std::vector<int> dimensions ( )
        {
            return map_->dimensions ( );
        };

        XdmfBaseItem* map_;
        XdmfBaseItem* data_;
    };

    struct XdmfFunction : XdmfBaseItem
    {

        XdmfFunction ( std::vector<XdmfBaseItem*> data_items, std::string function ) : data_items_ ( data_items ), function_ ( function )
        {
            dimensions_.resize ( 2 );
            assert ( data_items_[0]->dimensions ( ).size ( ) == 1 ); // only support one dim DataItems for Functions
            dimensions_[0] = data_items_[0]->dimensions ( )[0];
            dimensions_[1] = data_items_.size ( );
        };

        TiXmlElement* get_xdmf_element ( )
        {
            assert ( data_items_.size ( ) != 0 );
            TiXmlElement* xdmf_function = new TiXmlElement ( "DataItem" );

            xdmf_function->SetAttribute ( "ItemType", "Function" );

            xdmf_function->SetAttribute ( "Dimensions", dim_to_string ( ) );

            std::stringstream function_text;
            function_text << function_ << "($0";
            for ( int v = 1; v < static_cast < int > ( data_items_.size ( ) ); v++ )
            {
                function_text << ", $" << v;
            }
            function_text << ")";

            xdmf_function->SetAttribute ( "Function", function_text.str ( ) );

            for ( int i = 0; i < static_cast < int > ( data_items_.size ( ) ); ++i )
            {
                xdmf_function->LinkEndChild ( data_items_[i]->get_xdmf_element ( ) );
            }

            return xdmf_function;
        }

        void clean_deletion ( )
        {
            for ( int i = 0; i < static_cast < int > ( data_items_.size ( ) ); i++ )
            {
                data_items_[i]->clean_deletion ( );
                delete data_items_[i];
            }
            data_items_.clear ( );
        }

        std::vector<int> dimensions ( )
        {
            return dimensions_;
        };

        std::vector<int> dimensions_;
        std::vector<XdmfBaseItem*> data_items_;
        std::string function_;
    };

    struct XdmfAttribute
    {

        XdmfAttribute ( XdmfBaseItem* data_item, std::string name, std::string att_type, std::string center )
        : data_item_ ( data_item ), name_ ( name ), att_type_ ( att_type ), center_ ( center )
        {
        };

        TiXmlElement* get_xdmf_element ( )
        {
            assert ( data_item_ != NULL );
            TiXmlElement* attribute = new TiXmlElement ( "Attribute" );

            attribute->SetAttribute ( "Name", name_ );
            attribute->SetAttribute ( "AttributeType", att_type_ );
            attribute->SetAttribute ( "Center", center_ );

            attribute->LinkEndChild ( data_item_->get_xdmf_element ( ) );
            return attribute;
        }

        void clean_deletion ( )
        {
            data_item_->clean_deletion ( );
            delete data_item_;
        }

        XdmfBaseItem* data_item_;
        std::string name_;
        std::string att_type_;
        std::string center_;
    };

    struct XdmfTopology
    {

        XdmfTopology ( XdmfBaseItem* data_item, int number_of_elements )
        : data_item_ ( data_item ), number_of_elements_ ( number_of_elements )
        {
        };

        TiXmlElement* get_xdmf_element ( )
        {
            assert ( data_item_ != NULL );
            TiXmlElement* topology = new TiXmlElement ( "Topology" );

            topology->SetAttribute ( "Type", "Mixed" );
            topology->SetAttribute ( "NumberOfElements", number_of_elements_ );

            topology->LinkEndChild ( data_item_->get_xdmf_element ( ) );

            return topology;
        }

        void clean_deletion ( )
        {
            data_item_->clean_deletion ( );
            delete data_item_;
        }

        XdmfBaseItem* data_item_;
        int number_of_elements_;
    };

    struct XdmfGeometry
    {

        XdmfGeometry ( XdmfBaseItem* data_item, int gdim )
        : data_item_ ( data_item ), gdim_ ( gdim )
        {
        };

        TiXmlElement* get_xdmf_element ( )
        {
            assert ( data_item_ != NULL );
            TiXmlElement* geometry = new TiXmlElement ( "Geometry" );

            std::string type;
            if ( gdim_ == 2 )
            {
                type = "XY";
            }
            else if ( gdim_ == 3 )
            {
                type = "XYZ";
            }
            else
            {
                std::cerr << "XdmfGeometry: Only gdim 2 or 3 supported." << std::endl;
            }
            geometry->SetAttribute ( "Type", type );

            geometry->LinkEndChild ( data_item_->get_xdmf_element ( ) );

            return geometry;
        }

        void clean_deletion ( )
        {
            data_item_->clean_deletion ( );
            delete data_item_;
        }

        XdmfBaseItem* data_item_;
        int gdim_;
    };
}
#endif
