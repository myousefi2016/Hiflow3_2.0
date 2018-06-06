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

#include "reader.h"

#include <iostream>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <tr1/unordered_map>
#include <vector>
#include <utility>

#include <tinyxml.h>

#include "common/log.h"

#include "attributes.h"
#include "types.h"
#include "iterator.h"

#include "mesh_db_view.h"
#include "mesh_database.h"
#include "mesh_pXest.h"
#include "mesh_database_pXest.h"

#ifdef WITH_P4EST
#    include "p4est.h"
#    include "p4est_connectivity.h"
#    include "p8est.h"
#    include "p8est_connectivity.h"
#endif

const int DEBUG_LEVEL = 1;

namespace hiflow
{
    namespace mesh
    {

        ///////////////////////////////////////////////////
        //////// Reader ///////////////////////////////////
        ///////////////////////////////////////////////////

        /// \details Uses the MeshBuilder to build up the mesh.
        /// \param mb MeshBuilder pointer to the MeshBuilder it works with.

        Reader::Reader ( MeshBuilder* mb )
        : mesh_builder_ ( mb )
        {
        }

        /// \details Reads entities from file using the MeshBuilder. A view
        /// with the added entities is returned in mesh, if it is
        /// non-zero. The caller is responsible for deleting this mesh
        /// object.
        /// \param [in]  filename The name of the mesh file to be read.
        /// \param [out] mesh     The view of the read-in mesh.

        void Reader::read ( const char* filename, MeshPtr& mesh )
        {
            try
            {
                // clean builder
                mesh_builder_->clear ( );
                this->read_file ( filename, mesh );

                LOG_INFO ( "Filename", filename );
                LOG_INFO ( "Points", mesh->num_entities ( 0 ) );
                LOG_INFO ( "Cells", mesh->num_entities ( mesh->tdim ( ) ) );
                LOG_INFO ( "Topological dimension", mesh->tdim ( ) );
                LOG_INFO ( "Geometrical dimension", mesh->gdim ( ) );
            }
            catch ( const ReadGridException& exc )
            {
                std::cerr << "Error reading grid from file " << filename << "\n";
                std::cerr << exc.what ( ) << "\n";
                exit ( 1 );
            }
        }

        MeshBuilder* Reader::mesh_builder ( )
        {
            return mesh_builder_;
        }

        ////////////////////////////////////////////////////
        //////// UCDReader /////////////////////////////////
        ////////////////////////////////////////////////////

        /// \brief Helper structure for UCDReader

        struct UcdEntityDescription
        {
            TDim dim;
            int num_verts;
            std::vector<Id> vertices;
        };

        /// \brief Creates UcdEntityDescription based on name of type.

        UcdEntityDescription create_entity_description ( std::string type )
        {
            UcdEntityDescription descr;
            if ( type == std::string ( "line" ) )
            {
                descr.dim = 1;
                descr.num_verts = 2;
            }
            else if ( type == std::string ( "tri" ) )
            {
                descr.dim = 2;
                descr.num_verts = 3;
            }
            else if ( type == std::string ( "quad" ) )
            {
                descr.dim = 2;
                descr.num_verts = 4;
            }
            else if ( type == std::string ( "tet" ) )
            {
                descr.dim = 3;
                descr.num_verts = 4;
            }
            else if ( type == std::string ( "hex" ) )
            {
                descr.dim = 3;
                descr.num_verts = 8;
            }
            else if ( type == std::string ( "pyr" ) )
            {
                descr.dim = 3;
                descr.num_verts = 5;
            }
            else if ( type == std::string ( "prism" ) )
            {
                descr.dim = 3;
                descr.num_verts = 6;
            }
            else
            {
                descr.dim = -1;
                descr.num_verts = 0;
            }
            return descr;
        }

        /// \brief Helper structure for UCDPXestReader -> need to restrict allowed  entities

        struct UcdPXestEntityDescription
        {
            TDim dim;
            int num_verts;
            std::vector<Id> vertices;
        };

        /// \brief Creates UcdEntityDescription based on name of type.

        UcdPXestEntityDescription create_entity_description_pXest ( std::string type )
        {
            UcdPXestEntityDescription descr;
            if ( type == std::string ( "line" ) )
            {
                descr.dim = 1;
                descr.num_verts = 2;
            }
            else if ( type == std::string ( "quad" ) )
            {
                descr.dim = 2;
                descr.num_verts = 4;
            }
            else if ( type == std::string ( "hex" ) )
            {
                descr.dim = 3;
                descr.num_verts = 8;
            }
            else
            {
                descr.dim = -1;
                descr.num_verts = 0;
            }
            return descr;
        }

        UcdReader::UcdReader ( MeshBuilder* mb )
        : Reader ( mb )
        {
        }

        void UcdReader::read_file ( const char* filename, MeshPtr& mesh )
        {
            std::ifstream s ( filename );
            if ( !s )
            {
                std::cerr << "\nCan't open the file: " << filename << std::endl;
                exit ( 1 );
            }

            // Current line of the grid
            std::string current_strline;

            // get topological and geometrical dimension of the mesh
            // database
            const TDim tdim = mesh_builder ( )->tdim ( );
            const GDim gdim = mesh_builder ( )->gdim ( );

            // Filter comments
            do
            {
                getline ( s, current_strline );
            }
            while ( current_strline.find ( "#" ) != std::string::npos );

            // Read first line
            std::istringstream first_line ( current_strline.c_str ( ) );

            EntityCount num_nodes, num_entities;

            // we need only the information on the number of points and cells
            // from the first line.
            first_line >> num_nodes >> num_entities;

            assert ( num_nodes > 0 );
            assert ( num_entities > 0 );

            // Read the vertices
            std::tr1::unordered_map<int, Id> vertex_map;
            std::vector<Coordinate> coords ( gdim );

            for ( EntityCount index = 0; index < num_nodes; ++index )
            {
                // read line
                getline ( s, current_strline );
                std::istringstream current_point ( current_strline.c_str ( ) );

                // read UCD vertex id and coordinates
                Id idUCD;
                current_point >> idUCD;
                for ( int c = 0; c < gdim; ++c )
                {
                    current_point >> coords[c];
                }

                // add to database and vertex_map
                const Id vertex_id = mesh_builder ( )->add_vertex ( coords );
                vertex_map.insert ( std::make_pair ( idUCD, vertex_id ) );
            }

            std::vector<Id> cells;

            // Read the entities (lines, trias, quads, tetras, hexas)
            for ( EntityCount index = 0; index < num_entities; ++index )
            {
                getline ( s, current_strline );
                std::istringstream current_entity ( current_strline.c_str ( ) );

                // Get Id and Material Number
                Id idUCD;
                MaterialNumber mat_num;
                current_entity >> idUCD >> mat_num;

                // Get type
                std::string type;
                current_entity >> type;

                // create entity description object
                UcdEntityDescription entity_description = create_entity_description ( type );
                if ( entity_description.dim == -1 )
                {
                    // unknown entity type
                    throw ReadGridException (
                                              std::string ( "Unknown cell type: " ) += type,
                                              std::string ( current_strline ) );
                }

                // check dimension
                if ( entity_description.dim > tdim )
                {
                    // incorrect cell type for dimension
                    std::ostringstream err_msg;
                    err_msg << "Incorrect cell type " << type << " for mesh of dimension " << tdim;
                    throw ReadGridException ( err_msg.str ( ), std::string ( current_strline ) );
                }

                // read vertices, map to correct id with vertex_map
                for ( int v = 0; v < entity_description.num_verts; ++v )
                {
                    Id v_id;
                    current_entity >> v_id;
                    entity_description.vertices.push_back ( vertex_map[v_id] );
                }

                if ( type == "pyr" )
                {
                    std::rotate ( entity_description.vertices.begin ( ), entity_description.vertices.begin ( ) + 1, entity_description.vertices.end ( ) );
                }

                LOG_DEBUG ( 3, "Read entity " << string_from_range ( entity_description.vertices.begin ( ),
                                                                     entity_description.vertices.end ( ) ) );

                // read cells and facets and their material numbers
                if ( entity_description.dim == tdim || entity_description.dim == tdim - 1 )
                {
                    const Id entity_id = mesh_builder ( )->add_entity ( entity_description.dim, entity_description.vertices );
                    mesh_builder ( )->set_material_number ( entity_description.dim, entity_id, mat_num );
                }
            }

            // Create mesh containing cells
            mesh = mesh_builder ( )->build ( );
        }

        ////////////////////////////////////////////////////
        //////// VtkReader /////////////////////////////////
        ////////////////////////////////////////////////////

        VtkReader::VtkReader ( MeshBuilder* mb )
        : Reader ( mb )
        {
        }
        // reads a vtkXMLUnstructuredGrid, uses tinyXML to create a DOM

        void VtkReader::read_file ( const char* filename, MeshPtr& mesh )
        {

            /*
              vtk cell types:
             */
            enum VTKTYPES
            {
                VERTEX = 1,
                POLY_VERTEX = 2,
                LINE = 3,
                POLY_LINE = 4,
                TRIANGLE = 5,
                TRIANGLE_STRIP = 6,
                POLYGON = 7,
                PIXEL = 8,
                QUAD = 9,
                TETRA = 10,
                VOXEL = 11,
                HEXAHEDRON = 12,
                WEDGE = 13,
                PYRAMID = 14
            };

            // get topological and geometrical dimension of the mesh
            // database
            const TDim tdim = mesh_builder ( )->tdim ( );
            const GDim gdim = mesh_builder ( )->gdim ( );

            EntityCount num_nodes, num_cells;

            // load XML file
            TiXmlDocument doc ( filename );
            doc.LoadFile ( );

            // file could not be loaded, or vtu is corrupt
            assert ( doc.FirstChild ( ) != 0 );

            TiXmlElement* grid_root;

            // Paraview generates vtk files without "<?xml version="1.0" ?>"
            if ( doc.FirstChild ( "VTKFile" ) != 0 )
                grid_root = doc.FirstChild ( "VTKFile" )->ToElement ( );

            else if ( doc.FirstChild ( )->NextSibling ( "VTKFile" ) != 0 )
                grid_root = doc.FirstChild ( )->NextSibling ( "VTKFile" )->ToElement ( );

            else
            {
                std::cerr << "Input file " << filename << " is corrupt!";
                exit ( -1 );
            }

            // check type of grid (only unstructured for the moment)
            assert ( std::string ( grid_root->ToElement ( )->Attribute ( "type" ) ) == "UnstructuredGrid" );
            // check first child to be of the correct type
            assert ( std::string ( grid_root->FirstChild ( )->Value ( ) ) == std::string ( grid_root->ToElement ( )->Attribute ( "type" ) ) );

            TiXmlHandle docHandle ( &doc );
            TiXmlElement* piece_handle = docHandle.FirstChild ( "VTKFile" ).FirstChild ( "UnstructuredGrid" ).FirstChild ( "Piece" ).Element ( );

            if ( piece_handle )
            {
                // check for the correct attribute names
                assert ( std::string ( piece_handle->FirstAttribute ( )->Name ( ) ) == "NumberOfPoints" );
                assert ( std::string ( piece_handle->LastAttribute ( )->Name ( ) ) == "NumberOfCells" );

                // get number of nodes and cells
                num_nodes = atoi ( piece_handle->FirstAttribute ( )->Value ( ) );
                num_cells = atoi ( piece_handle->LastAttribute ( )->Value ( ) );

                // check if data is useful
                assert ( num_nodes > 0 );
                assert ( num_cells > 0 );

                // get XML nodes of the specific data
                TiXmlNode* pointData = piece_handle->FirstChildElement ( "PointData" );
                TiXmlNode* cellData = piece_handle->FirstChildElement ( "CellData" );
                TiXmlNode* points = piece_handle->FirstChildElement ( "Points" );
                TiXmlNode* cells = piece_handle->FirstChildElement ( "Cells" );

                // possible use of openMP here!
                //      int proc_num = omp_get_num_procs();
                //      omp_set_num_threads(proc_num);

                // PointData ////////
                TiXmlElement* pointDataArray;
                std::vector<std::string> pointDataArrayNames;

                // TODO(Thomas): Use not only double data!
                std::vector<std::vector<double> > pointDataArrayVectors;
                if ( pointData != 0 )
                {
                    if ( pointData->FirstChild ( "DataArray" ) != 0 )
                    {
                        pointDataArray = pointData->FirstChild ( "DataArray" )->ToElement ( );
                        // loop through all pointDataArrays in pointData
                        for (; pointDataArray; pointDataArray = pointDataArray->NextSiblingElement ( ) )
                        {

                            const char * num_comps = pointDataArray->Attribute ( "NumberOfComponents" );
                            int number_of_comps = 1;
                            if ( num_comps != 0 )
                            {
                                number_of_comps = atoi ( num_comps );
                            }
                            // get data value
                            const char* data = pointDataArray->ToElement ( )->GetText ( );
                            if ( data == 0 ) continue;
                            std::istringstream current_data ( data );

                            // store name of data array
                            pointDataArrayNames.push_back ( std::string ( pointDataArray->Attribute ( "Name" ) ) );

                            // store data value in vector
                            std::vector<double> temp_data_vector ( number_of_comps * num_nodes );
                            for ( EntityCount index = 0;
                                  index < ( number_of_comps * num_nodes );
                                  ++index )
                            {
                                current_data >> temp_data_vector[index];
                            }

                            if ( number_of_comps > 1 )
                            {
                                pointDataArrayVectors.resize ( number_of_comps );
                                // TODO(Thomas): Resize Name array as
                                // well, and give names for the additional
                                // data arrays.
                                std::string tmp_name = pointDataArrayNames.back ( );
                                pointDataArrayNames.resize ( number_of_comps );
                                for ( int i = 0; i < number_of_comps; ++i )
                                {
                                    for ( std::vector<double>::const_iterator it = temp_data_vector.begin ( ) + i;
                                          it < temp_data_vector.end ( ); it = it + number_of_comps )
                                    {
                                        pointDataArrayVectors[i].push_back ( *it );
                                    }
                                    std::ostringstream index_str;
                                    index_str << i;
                                    pointDataArrayNames[i] = tmp_name + index_str.str ( );

                                }
                            }
                            else
                            {
                                pointDataArrayVectors.push_back ( temp_data_vector );
                            }

                        }
                        assert ( pointDataArrayVectors.size ( ) == pointDataArrayNames.size ( ) );
                    }
                }
                // CellData /////////

                // Vectors to store material numbers.
                std::vector<MaterialNumber> material_numbers ( num_cells );

                // Vectors to store the hierarchy of the grid if given.
                std::vector<Id> parent ( num_cells ), first_child ( num_cells ), num_child ( num_cells );
                TiXmlElement* cellDataArray;

                bool has_cell_data = true;
                if ( cellData != 0 )
                {
                    if ( cellData->FirstChild ( "DataArray" ) != 0 )
                    {
                        cellDataArray = cellData->FirstChild ( "DataArray" )->ToElement ( );
                    }
                    else
                    {
                        std::cerr << "Are you sure? You don't even provide a material id?!\n"
                                << "No good can come from this...\n";
                        std::cerr << "I'll set some material id for you. Will be all -1.\n";
                        for ( int i = 0; i < num_cells; ++i )
                            material_numbers.push_back ( -1 );
                        has_cell_data = false;
                    }
                    if ( has_cell_data )
                    {
                        // loop through all cellDataArrays in CellData
                        for (; cellDataArray; cellDataArray = cellDataArray->NextSiblingElement ( ) )
                        {
                            const char* data = cellDataArray->ToElement ( )->GetText ( );
                            std::istringstream current_data ( data );

                            if ( std::string ( cellDataArray->Attribute ( "Name" ) ) == "Material Id" )
                            {
                                // Vector to store material numbers on the mesh.
                                for ( EntityCount index = 0; index < num_cells; ++index )
                                {
                                    current_data >> material_numbers[index];
                                }
                            }
                            else if ( std::string ( cellDataArray->Attribute ( "Name" ) ) == "Parent" )
                            {
                                for ( EntityCount index = 0; index < num_cells; ++index )
                                {
                                    current_data >> parent[index];
                                }
                            }
                            else if ( std::string ( cellDataArray->Attribute ( "Name" ) ) == "NumberOfChildren" )
                            {
                                for ( EntityCount index = 0; index < num_cells; ++index )
                                {
                                    current_data >> num_child[index];
                                }
                            }
                            else if ( std::string ( cellDataArray->Attribute ( "Name" ) ) == "FirstChild" )
                            {
                                for ( EntityCount index = 0; index < num_cells; ++index )
                                {
                                    current_data >> first_child[index];
                                }
                            }
                            else
                            {
                                std::cerr << "Could not read Unknown CellData type \"" << cellDataArray->Attribute ( "Name" ) << "\"\t... skipping.\n";
                            }
                        }
                    }
                }

                // Points ///////////
                std::tr1::unordered_map<int, Id> vertex_map;
                std::vector<Coordinate> coords ( gdim );

                // get number of components per point
                const int num_comps_available = atoi ( points->FirstChild ( "DataArray" )->ToElement ( )->Attribute ( "NumberOfComponents" ) );
                const int num_comps = gdim;
                assert ( num_comps <= num_comps_available );

                // is there a better way to convert this in a vector<double>?
                const char* node_vec = points->FirstChild ( "DataArray" )->ToElement ( )->GetText ( );
                assert ( node_vec != 0 );
                std::istringstream current_point ( node_vec );

                for ( EntityCount index = 0; index < num_nodes; ++index )
                {
                    for ( int c = 0; c < num_comps; ++c )
                    {
                        current_point >> coords[c];
                    }
                    // read extra components
                    double dummy_comp;
                    for ( int c = 0; c < num_comps_available - num_comps; ++c )
                    {
                        current_point >> dummy_comp;
                    }

                    // add to database and vertex_map
                    const Id vertex_id = mesh_builder ( )->add_vertex ( coords );
                    vertex_map.insert ( std::make_pair ( index, vertex_id ) );
                }

                // Cells ////////////
                TiXmlNode* connectivities = cells->FirstChild ( "DataArray" );
                TiXmlNode* offsets = connectivities->NextSibling ( "DataArray" );
                TiXmlNode* types = offsets->NextSibling ( "DataArray" );

                assert ( std::string ( connectivities->ToElement ( )->Attribute ( "Name" ) ) == "connectivity" );
                assert ( std::string ( offsets->ToElement ( )->Attribute ( "Name" ) ) == "offsets" );
                assert ( std::string ( types->ToElement ( )->Attribute ( "Name" ) ) == "types" );

                // get offsets
                std::vector<int> offs ( num_cells + 1 );
                offs[0] = 0;
                const char* offs_char = offsets->ToElement ( )->GetText ( );
                std::istringstream current_offs ( offs_char );

                for ( EntityCount index = 1; index < num_cells + 1; ++index )
                {
                    current_offs >> offs[index];
                }

                // get types of elements
                std::vector<int> typs ( num_cells );
                const char* typs_char = types->ToElement ( )->GetText ( );
                std::istringstream current_types ( typs_char );
                for ( EntityCount index = 0; index < num_cells; ++index )
                {
                    current_types >> typs[index];
                }

                // get connectivities of elements with respect to type and
                // set the material id of the corresponding element.

                // create connectivity vector with a max of 8
                // (HEXAHEDRON)
                std::vector<Id> cons;
                cons.reserve ( 8 );

                const char* cons_char = connectivities->ToElement ( )->GetText ( );
                std::istringstream current_cons ( cons_char );
                Id con;
                Id entity_id;

                // Create Material Id structure
                //            ScopedArray<std::vector<MaterialNumber> >::Type material_id(new std::vector<MaterialNumber>[tdim]);

                ScopedArray< std::vector<Id> >::Type parent_id ( new std::vector<Id>[tdim] );
                ScopedArray< std::vector<Id> >::Type children_id ( new std::vector<Id>[tdim] );
                ScopedArray< std::vector<Id> >::Type first_child_id ( new std::vector<Id>[tdim] );

                for ( EntityCount index = 1; index < num_cells + 1; ++index )
                {
                    for ( int offset = 0; offset < offs[index] - offs[index - 1]; ++offset )
                    {
                        current_cons >> con;
                        cons.push_back ( vertex_map[con] );
                    }
                    TDim tdim_ent;
                    switch ( typs[index - 1] )
                    {
                        case LINE:
                            tdim_ent = 1;
                            assert ( tdim == 2 || tdim == 1 );
                            assert ( cons.size ( ) == 2 );
                            entity_id = mesh_builder ( )->add_entity ( tdim_ent, cons );
                            mesh_builder ( )->set_material_number ( tdim_ent, entity_id, material_numbers[index - 1] );
                            //                    material_id[tdim_ent - 1].push_back(material_numbers[entity_id]);
                            parent_id[tdim_ent - 1].push_back ( parent[entity_id] );
                            children_id[tdim_ent - 1].push_back ( num_child[entity_id] );
                            first_child_id[tdim_ent - 1].push_back ( first_child[entity_id] );
                            break;
                        case TRIANGLE:
                            tdim_ent = 2;
                            assert ( tdim > 1 );
                            assert ( cons.size ( ) == 3 );
                            entity_id = mesh_builder ( )->add_entity ( tdim_ent, cons );
                            mesh_builder ( )->set_material_number ( tdim_ent, entity_id, material_numbers[index - 1] );
                            //          if (material_numbers[entity_id] != 4)
                            //                    material_id[tdim_ent - 1].push_back(material_numbers[entity_id]);
                            parent_id[tdim_ent - 1].push_back ( parent[entity_id] );
                            children_id[tdim_ent - 1].push_back ( num_child[entity_id] );
                            first_child_id[tdim_ent - 1].push_back ( first_child[entity_id] );
                            break;
                        case QUAD:
                            tdim_ent = 2;
                            assert ( tdim > 1 );
                            assert ( cons.size ( ) == 4 );
                            entity_id = mesh_builder ( )->add_entity ( tdim_ent, cons );
                            mesh_builder ( )->set_material_number ( tdim_ent, entity_id, material_numbers[index - 1] );
                            //                    material_id[tdim_ent - 1].push_back(material_numbers[entity_id]);
                            parent_id[tdim_ent - 1].push_back ( parent[entity_id] );
                            children_id[tdim_ent - 1].push_back ( num_child[entity_id] );
                            first_child_id[tdim_ent - 1].push_back ( first_child[entity_id] );
                            break;
                        case TETRA:
                            tdim_ent = 3;
                            assert ( tdim == 3 );
                            assert ( cons.size ( ) == 4 );
                            entity_id = mesh_builder ( )->add_entity ( tdim_ent, cons );
                            mesh_builder ( )->set_material_number ( tdim_ent, entity_id, material_numbers[index - 1] );

                            //                    material_id[tdim_ent - 1].push_back(material_numbers[entity_id]);
                            parent_id[tdim_ent - 1].push_back ( parent[entity_id] );
                            children_id[tdim_ent - 1].push_back ( num_child[entity_id] );
                            first_child_id[tdim_ent - 1].push_back ( first_child[entity_id] );
                            break;
                        case VOXEL:
                        { // brackets needed because types are
                            // initialized in here

                            // Treat voxels like hexahedrons, but
                            // resort cons (connectivity) vector
                            // first!
                            assert ( cons.size ( ) == 8 );
                            std::vector<int> tmp_cons ( cons );
                            int permute[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
                            for ( int i = 0; i < static_cast < int > ( cons.size ( ) ); ++i )
                            {
                                cons[i] = tmp_cons[permute[i]];
                            }
                        }
                            // no break!
                        case HEXAHEDRON:
                            tdim_ent = 3;
                            assert ( cons.size ( ) == 8 );
                            assert ( tdim == 3 );
                            entity_id = mesh_builder ( )->add_entity ( tdim_ent, cons );
                            mesh_builder ( )->set_material_number ( tdim_ent, entity_id, material_numbers[index - 1] );
                            //                    material_id[tdim_ent - 1].push_back(material_numbers[entity_id]);
                            parent_id[tdim_ent - 1].push_back ( parent[entity_id] );
                            children_id[tdim_ent - 1].push_back ( num_child[entity_id] );
                            first_child_id[tdim_ent - 1].push_back ( first_child[entity_id] );
                            break;
                        case VERTEX:
                            std::cerr << "Vertices are currently not supported!\n";
                            NOT_YET_IMPLEMENTED;
                            break;
                        case POLY_VERTEX:
                            std::cerr << "Poly vertices are currently not supported!\n";
                            NOT_YET_IMPLEMENTED;
                            break;
                        case POLY_LINE:
                            std::cerr << "Poly lines are currently not supported!\n";
                            NOT_YET_IMPLEMENTED;
                            break;
                        case TRIANGLE_STRIP:
                            std::cerr << "Triangle strips are currently not supported!\n";
                            NOT_YET_IMPLEMENTED;
                            break;
                        case POLYGON:
                            std::cerr << "Polygons are currently not supported!\n";
                            NOT_YET_IMPLEMENTED;
                            break;
                        case PIXEL:
                            std::cerr << "Pixels are currently not supported!\n";
                            NOT_YET_IMPLEMENTED;
                            break;
                        case WEDGE:
                            std::cerr << "Wedges are currently not supported!\n";
                            NOT_YET_IMPLEMENTED;
                            break;
                        case PYRAMID:
                            std::cerr << "Pyramids are currently not supported!\n";
                            NOT_YET_IMPLEMENTED;
                            break;
                        default:
                            std::cerr << "This type of element is unknown!\n";
                            NOT_YET_IMPLEMENTED;
                    }
                    cons.clear ( );
                }
                // Create mesh containing cells
                mesh = mesh_builder ( )->build ( );

                // add point data to mesh
                assert ( pointDataArrayNames.size ( ) == pointDataArrayVectors.size ( ) );
                std::vector<std::vector<double> >::const_iterator point_data_it = pointDataArrayVectors.begin ( );
                for ( std::vector<std::string>::const_iterator names_it = pointDataArrayNames.begin ( ); names_it != pointDataArrayNames.end ( ); ++names_it )
                {
                    AttributePtr point_attribute ( new DoubleAttribute ( *point_data_it ) );
                    mesh->add_attribute ( *names_it, 0, point_attribute );
                    ++point_data_it;
                }

            }
            else
            {
                std::cerr << "Invalid vtk grid file!\n";
                exit ( -1 );
            }
        }

        ////////////////////////////////////////////////////
        //////// Parallel VtkReader ////////////////////////
        ////////////////////////////////////////////////////
        /// \param [in] mb       The MeshBuilder used to create the mesh.
        /// \param [in] mpi_comm An MPI Communicator.

        PVtkReader::PVtkReader ( MeshBuilder* mb, const MPI_Comm& mpi_comm )
        : Reader ( mb ), mpi_comm_ ( mpi_comm )
        {
        }

        /// \details Reads a parallel vtkXMLUnstructuredGrid, uses
        /// tinyXML to create a DOM and VtkReader to read the single
        /// files.

        void PVtkReader::read_file ( const char* filename, MeshPtr& mesh )
        {
            int rank = -1, size = -1;
            MPI_Comm_rank ( mpi_comm_, &rank );
            MPI_Comm_size ( mpi_comm_, &size );

            // load XML file
            TiXmlDocument doc ( filename );
            doc.LoadFile ( );

            // file could not be loaded, or vtu is corrupt
            assert ( doc.FirstChild ( ) != 0 );

            TiXmlElement* grid_root = doc.FirstChild ( )->NextSibling ( "VTKFile" )->ToElement ( );

            // check type of grid -> Parallel Unstructured Grid
            assert ( std::string ( grid_root->ToElement ( )->Attribute ( "type" ) ) == "PUnstructuredGrid" );
            // check first child to be of the correct type
            assert ( std::string ( grid_root->FirstChild ( )->Value ( ) ) == std::string ( grid_root->ToElement ( )->Attribute ( "type" ) ) );

            TiXmlHandle docHandle ( &doc );
            TiXmlElement* pUG_handle = docHandle.FirstChild ( "VTKFile" ).FirstChild ( "PUnstructuredGrid" ).Element ( );

            if ( pUG_handle )
            {

                // get XML nodes of the specific data
                TiXmlNode* piece = pUG_handle->FirstChildElement ( "Piece" );

                // get the correct filename including the path
                std::istringstream filename_root_dir ( filename );
                std::size_t dir = filename_root_dir.str ( ).find_last_of ( "/\\" );

                std::string path = filename_root_dir.str ( ).substr ( 0, dir );
                if ( path == filename )
                {
                    path = ".";
                }
                LOG_DEBUG ( 3, "Pathname to file: " << path );
                assert ( !path.empty ( ) );

                // Loop through names of actual data files an pass them to
                // the respective VtkReader.
                int i = 0;
                while ( piece )
                {
                    if ( rank == i )
                    {
                        const char* cstr_src_filename = piece->ToElement ( )->Attribute ( "Source" );
                        assert ( cstr_src_filename != 0 );

                        std::stringstream src_filename ( cstr_src_filename );
                        std::string str_src_filename = ( path + "/" + src_filename.str ( ) );
                        LOG_DEBUG ( 3, "Filename: " << str_src_filename );
                        assert ( !str_src_filename.empty ( ) );
                        VtkReader reader ( mesh_builder ( ) );
                        reader.read ( str_src_filename.c_str ( ), mesh );
                    }
                    else
                    {
                        piece = piece->NextSibling ( "Piece" );
                    }
                    ++i;
                }
            }
        }

        UcdPXestReader::UcdPXestReader ( MeshBuilder* mb )
        : Reader ( mb )
        {
        }

        void UcdPXestReader::read_file ( const char* filename, MeshPtr& mesh )
        {
#ifdef WITH_P4EST
            // typecast builder 
            MeshPXestBuilder* builder_pXest = static_cast < MeshPXestBuilder* > ( mesh_builder ( ) );

            // mesh builder must have an empty database
            const int tdim = builder_pXest->tdim ( );
            const int gdim = builder_pXest->gdim ( );

            assert ( builder_pXest->get_db ( )->num_entities ( 0 ) == 0 );
            assert ( builder_pXest->get_db ( )->num_entities ( tdim ) == 0 );

            std::ifstream s ( filename );
            if ( !s )
            {
                std::cerr << "\nCan't open the file: " << filename << std::endl;
                exit ( 1 );
            }

            // Current line of the grid
            std::string current_strline;

            // Filter comments
            do
            {
                getline ( s, current_strline );
            }
            while ( current_strline.find ( "#" ) != std::string::npos );

            // Read first line
            std::istringstream first_line ( current_strline.c_str ( ) );

            EntityCount num_nodes, num_entities;

            // we need only the information on the number of points and cells
            // from the first line.
            first_line >> num_nodes >> num_entities;

            assert ( num_nodes > 0 );
            assert ( num_entities > 0 );
            /// !
            // Read the vertices
            std::tr1::unordered_map<int, Id> vertex_map;
            std::vector<Coordinate> coords ( gdim );

            int v_shift = 0;
            for ( EntityCount index = 0; index < num_nodes; ++index )
            {
                // read line
                getline ( s, current_strline );
                std::istringstream current_point ( current_strline.c_str ( ) );

                // read UCD vertex id and coordinates
                Id idUCD;
                current_point >> idUCD;

                for ( int c = 0; c < gdim; ++c )
                {
                    current_point >> coords[c];
                }

                // add to database and vertex_map
                const Id vertex_id = builder_pXest->add_vertex ( coords );
                vertex_map.insert ( std::make_pair ( idUCD, vertex_id ) );
            }

            // store coordinates of vertices: Mind periodicity!
            int num_db_vertices = builder_pXest->get_db ( )->num_entities ( 0 );
            std::vector<double> pXest_vertices ( num_db_vertices * 3, 0 );

            for ( Id v_id = 0; v_id < num_db_vertices; ++v_id )
            {
                std::vector<Coordinate> coord = builder_pXest->get_db ( )->get_coordinates ( v_id );
                std::vector<Coordinate> v_coord = periodify ( coord, gdim, mesh_builder ( )->get_period ( ) );

                // 3rd coordinate stays 0 in case of 2d
                for ( int c = 0; c < gdim; ++c )
                {
                    pXest_vertices[3 * v_id + c] = v_coord[c];
                }
                if ( gdim == 2 )
                {
                    pXest_vertices[3 * v_id + 2] = 0.;
                }
            }

            LOG_DEBUG ( 1, "Number of vertices in inp file: " << num_nodes << ", number of vertices in database: " << num_db_vertices );

            std::vector<int> pXest_tree_to_vertices ( 0 );
            int num_cells = 0;
            int e_shift = 0;
            // Read the entities (lines, trias, quads, tetras, hexas)
            for ( EntityCount index = 0; index < num_entities; ++index )
            {
                getline ( s, current_strline );
                std::istringstream current_entity ( current_strline.c_str ( ) );

                // Get Id and Material Number
                Id idUCD;
                MaterialNumber mat_num;
                current_entity >> idUCD >> mat_num;

                // Get type
                std::string type;
                current_entity >> type;

                // create entity description object
                UcdPXestEntityDescription entity_description = create_entity_description_pXest ( type );
                if ( entity_description.dim == -1 )
                {
                    // unknown entity type
                    throw ReadGridException (
                                              std::string ( "Unknown cell type: " ) += type,
                                              std::string ( current_strline ) );
                }

                // check dimension
                if ( entity_description.dim > tdim )
                {
                    // incorrect cell type for dimension
                    std::ostringstream err_msg;
                    err_msg << "Incorrect cell type " << type << " for mesh of dimension " << tdim;
                    throw ReadGridException ( err_msg.str ( ), std::string ( current_strline ) );
                }

                std::vector<Id> v_ids ( entity_description.num_verts );
                // read vertices, map to correct id with vertex_map
                for ( int v = 0; v < entity_description.num_verts; ++v )
                {
                    Id v_id;
                    current_entity >> v_id;
                    entity_description.vertices.push_back ( vertex_map[v_id] );
                    v_ids[v] = vertex_map[v_id];
                }

                LOG_DEBUG ( 3, "Read entity " << string_from_range ( entity_description.vertices.begin ( ),
                                                                     entity_description.vertices.end ( ) ) );

                // read cells and facets and their material numbers
                Id entity_id = -1;
                if ( entity_description.dim == tdim || entity_description.dim == tdim - 1 )
                {
                    // create entity in mesh database
                    entity_id = builder_pXest->add_entity ( entity_description.dim, entity_description.vertices );
                    builder_pXest->set_material_number ( entity_description.dim, entity_id, mat_num );
                }

                // get cell
                if ( entity_description.dim == tdim )
                {
                    std::vector<Id> v ( entity_description.num_verts );
                    v[0] = v_ids[0];
                    v[1] = v_ids[1];
                    v[2] = v_ids[3];
                    v[3] = v_ids[2];
                    if ( type == "hex" )
                    {
                        v[4] = v_ids[4];
                        v[5] = v_ids[5];
                        v[6] = v_ids[7];
                        v[7] = v_ids[6];
                    }

                    for ( int n = 0; n < entity_description.num_verts; ++n )
                    {
                        pXest_tree_to_vertices.push_back ( v[n] );
                    }

                    // store entitiy(tdim)-to-quad map
                    // [cell_id, ref_level, x, y, z, (as topological coordinates inside of cell), localId]
                    assert ( num_cells == entity_id );

                    QuadCoord coord ( num_cells, 0, 0, 0, 0, 0, tdim, entity_id );
                    builder_pXest->add_entity_to_quad_coord ( tdim, entity_id, coord, 0 );

                    num_cells++;
                }

            }
            int num_verts = static_cast < int > ( std::pow ( static_cast < double > ( 2 ), tdim ) );
            assert ( num_cells * num_verts == pXest_tree_to_vertices.size ( ) );

            // set data needed for creating p4est connectivity
            builder_pXest->set_conn_data ( num_nodes, num_cells, pXest_vertices, pXest_tree_to_vertices );

            // Create mesh containing cells
            mesh = builder_pXest->build ( );

#endif
        }

    }
} // namespace hiflow
