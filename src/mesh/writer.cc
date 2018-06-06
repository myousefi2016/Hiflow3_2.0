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

#include "writer.h"

#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <tr1/unordered_map>

#include <mpi.h>
#include <tinyxml.h>

#include "common/log.h"

#include "attributes.h"
#include "cell_type.h"
#include "iterator.h"

const int DEBUG_LEVEL = 0;

namespace hiflow
{
    namespace mesh
    {

        class Entity;

        //////// Writer ///////////////////////

        /// \param [in] filename The name of the mesh file to be written.
        /// \param [in] mesh The view of the mesh to be written.

        void Writer::write ( const char* filename, const Mesh& mesh )
        {
            this->write_file ( filename, mesh );
        }

        /// \details Makes the writer aware of data that got attached to
        /// the mesh in order to be able to write it. If the specified
        /// data by means of its string and topological dimension doesn:t
        /// exist, it will print only -1:s for all entities.
        /// \param [in] array_name The name of the data array in the mesh.
        /// \param [in] tdim The topological dimension associated to the
        /// data in the mesh.

        void Writer::add_data_array ( std::string array_name, TDim tdim )
        {
            data_to_write_.push_back ( std::make_pair ( array_name, tdim ) );
        }

        /// \details Informs the writer that the attribute with the given
        /// name should be written to the output file.
        ///
        /// \param [in] attribute_name The name of the attribute to be
        /// written.
        /// \param [in] tdim The topological dimension of the
        /// entities for which the attribute is defined.

        void Writer::add_attribute ( const std::string& attribute_name, TDim tdim )
        {
            assert ( !attribute_name.empty ( ) );
            if ( find ( attributes_to_write_.begin ( ), attributes_to_write_.end ( ), std::make_pair ( attribute_name, tdim ) )
                 == attributes_to_write_.end ( ) )
            {
                attributes_to_write_.push_back ( std::make_pair ( attribute_name, tdim ) );
            }
        }

        void Writer::add_all_attributes ( const Mesh& mesh, bool internal_mesh_attributes )
        {
            for ( TDim i = 0; i <= mesh.tdim ( ); ++i )
            { // dimensions
                std::vector<std::string> names = mesh.get_attribute_names ( i );
                for ( std::vector<std::string>::const_iterator it = names.begin ( );
                      it != names.end ( ); ++it )
                { // attributes of dimension

                    // skip internal attributes (starting with '_') unless flag is set.
                    if ( !internal_mesh_attributes && ( *it )[0] == '_' )
                    {
                        continue;
                    }

                    // take out "Material Id" since this is added always,
                    // and paraview crashes if the same array exists
                    // twice.
                    if ( *it == "Material Id" )
                    {
                        continue;
                    }

                    // only add attributes once
                    if ( find ( attributes_to_write_.begin ( ), attributes_to_write_.end ( ), std::make_pair ( *it, i ) )
                         == attributes_to_write_.end ( ) )
                    {
                        attributes_to_write_.push_back ( std::make_pair ( *it, i ) );
                    }
                }
            }
        }

        //////// UcdWriter ////////////////////////

        void UcdWriter::write_file ( const char* filename, const Mesh& mesh ) const
        {
            std::ofstream file ( filename );
            if ( !file )
            {
                LOG_ERROR ( "Could not write mesh to file" << filename );
                return;
            }

            const TDim tdim = mesh.tdim ( );
            const EntityCount num_cells = mesh.num_entities ( tdim );
            const EntityCount num_verts = mesh.num_entities ( 0 );

            MeshPtr bdy_mesh = mesh.extract_boundary_mesh ( );
            const TDim bdim = bdy_mesh->tdim ( );
            const EntityCount num_bdy_cells = bdy_mesh->num_entities ( bdim );

            std::tr1::unordered_map<Id, EntityNumber> id_map;

            // Write header
            file << "# Ucd mesh file written using HiFlow3\n"
                    << num_verts << " " << num_cells + num_bdy_cells << " 0 0 0\n";

            // Write vertices
            for ( int v = 0; v < num_verts; ++v )
            {
                const Entity& vertex_entity = mesh.get_entity ( 0, v );
                std::vector<Coordinate> coords = mesh.get_coordinates ( 0, v );

                id_map.insert ( std::make_pair ( vertex_entity.id ( ), v ) );

                // make point three-dimensional if it is not already.
                if ( coords.size ( ) < 3 )
                {
                    coords.resize ( 3, 0. ); // resize does not modify existing elements
                }

                file << v << " " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
            }

            // Write cells
            std::vector<int> cell_verts;
            cell_verts.reserve ( 8 );
            for ( EntityIterator cell_it = mesh.begin ( tdim ), end_cell = mesh.end ( tdim );
                  cell_it != end_cell; ++cell_it )
            {

                cell_verts.clear ( );
                for ( IncidentEntityIterator vert_it = cell_it->begin_incident ( 0 ), end_vert = cell_it->end_incident ( 0 );
                      vert_it != end_vert; ++vert_it )
                {
                    cell_verts.push_back ( vert_it->index ( ) );
                }

                std::string cell_type = "unknown";
                if ( tdim == 2 )
                {
                    if ( cell_verts.size ( ) == 3 )
                    {
                        cell_type = "tri";
                    }
                    else if ( cell_verts.size ( ) == 4 )
                    {
                        cell_type = "quad";
                    }
                }
                else if ( tdim == 3 )
                {
                    if ( cell_verts.size ( ) == 4 )
                    {
                        cell_type = "tet";
                    }
                    else if ( cell_verts.size ( ) == 8 )
                    {
                        cell_type = "hex";
                    }
                    else if ( cell_verts.size ( ) == 5 )
                    {
                        cell_type = "pyr";
                    }
                }

                if ( cell_type == "unknown" )
                {
                    LOG_ERROR ( "Unknown cell type: tdim = " << tdim << ", num_verts = " << cell_verts.size ( ) );
                }

                const MaterialNumber mat_num = cell_it->get_material_number ( );

                file << cell_it->index ( ) << " " << mat_num << " " << cell_type << " ";
                if ( cell_type == "pyr" )
                {
                    file << " " << cell_verts[cell_verts.size ( ) - 1];
                    for ( int c = 0; c != static_cast < int > ( cell_verts.size ( ) - 1 ); ++c )
                    {
                        file << " " << cell_verts[c];
                    }

                }
                else
                {
                    for ( int c = 0; c != static_cast < int > ( cell_verts.size ( ) ); ++c )
                    {
                        file << " " << cell_verts[c];
                    }
                }

                file << "\n";
            }

            // Write boundary facets
            std::vector<int> facet_verts;
            facet_verts.reserve ( 8 );
            for ( EntityIterator facet_it = bdy_mesh->begin ( bdim ), end_facet = bdy_mesh->end ( bdim );
                  facet_it != end_facet; ++facet_it )
            {

                const MaterialNumber mat_num = facet_it->get_material_number ( );

                facet_verts.clear ( );
                facet_verts.insert ( facet_verts.end ( ),
                                     facet_it->begin_vertex_ids ( ), facet_it->end_vertex_ids ( ) );

                std::string facet_type = "unknown";
                if ( bdim == 1 )
                {
                    if ( facet_verts.size ( ) == 2 )
                    {
                        facet_type = "line";
                    }
                }
                else if ( bdim == 2 )
                {
                    if ( facet_verts.size ( ) == 3 )
                    {
                        facet_type = "tri";
                    }
                    else if ( facet_verts.size ( ) == 4 )
                    {
                        facet_type = "quad";
                    }
                }

                if ( facet_type == "unknown" )
                {
                    LOG_ERROR ( "Unknown facet type: bdim = " << bdim << ", num_verts = " << facet_verts.size ( ) );
                }

                file << facet_it->index ( ) << " " << mat_num << " " << facet_type << " ";
                for ( int c = 0; c != static_cast < int > ( facet_verts.size ( ) ); ++c )
                {
                    if ( id_map.count ( facet_verts[c] ) == 0 )
                    {
                        LOG_ERROR ( "No vertex with id " << facet_verts[c] << " found in id_map!" );
                        throw "Unknown vertex!";
                    }
                    file << " " << id_map[facet_verts[c]];
                }
                file << "\n";
            }
        }

        //////// VtkWriter ////////////////////////

        VtkWriter::VtkWriter ( ) : deformation_ ( 0 )
        {
        }

        void VtkWriter::set_deformation ( const std::vector<double>* deformation )
        {
            deformation_ = deformation;
        }

        void PVtkWriter::set_deformation ( const std::vector<double>* deformation )
        {
            deformation_ = deformation;
        }

        void VtkWriter::write_file ( const char* filename, const Mesh & mesh ) const
        {
            // mapping from the CellType::Type enum to the VTK cell types
            // see the VTK file formats description
            static const int vtk_cell_types[] = { 1, 3, 5, 9, 10, 12, 14 };

            const int tdim = mesh.tdim ( );

            TiXmlDocument doc;
            TiXmlDeclaration* decl = new TiXmlDeclaration ( "1.0", "", "" );
            doc.LinkEndChild ( decl );

            TiXmlElement* vtkFile = new TiXmlElement ( "VTKFile" );
            vtkFile->SetAttribute ( "type", "UnstructuredGrid" );
            vtkFile->SetAttribute ( "version", "0.1" );
            vtkFile->SetAttribute ( "byte_order", "LittleEndian" );
            vtkFile->SetAttribute ( "compressor", "vtkZLibDataCompressor" );
            doc.LinkEndChild ( vtkFile );

            TiXmlElement* umesh = new TiXmlElement ( "UnstructuredGrid" );
            vtkFile->LinkEndChild ( umesh );

            TiXmlElement* piece = new TiXmlElement ( "Piece" );
            piece->SetAttribute ( "NumberOfPoints", mesh.num_entities ( 0 ) );
            // TODO(Thomas): add boundary cells
            piece->SetAttribute ( "NumberOfCells", mesh.num_entities ( tdim ) );
            umesh->LinkEndChild ( piece );

            TiXmlElement* point_data = new TiXmlElement ( "PointData" );
            piece->LinkEndChild ( point_data );

            TiXmlElement* cell_data = new TiXmlElement ( "CellData" );
            // TODO(Thomas): set attributes scalars, vectors, etc.
            piece->LinkEndChild ( cell_data );

            // write point data
            TiXmlElement* data_array;

            std::stringstream os;

            // write attributes
            for ( DataVectorPair::const_iterator it = attributes_to_write_.begin ( );
                  it != attributes_to_write_.end ( ); ++it )
            {

                if ( !mesh.has_attribute ( it->first, it->second ) )
                {
                    LOG_ERROR ( "No attribute " << it->first << " for dimension " << it->second << " found in mesh" );
                    continue;
                }

                if ( ( it->second != 0 ) && ( it->second != tdim ) )
                {
                    LOG_ERROR ( "VTK does not support data for entities other than cells or vertices." );
                    continue;
                }

                data_array = new TiXmlElement ( "DataArray" );
                data_array->SetAttribute ( "Name", it->first );

                AttributePtr attr_ptr = mesh.get_attribute ( it->first, it->second );
                Attribute* attr = attr_ptr.get ( );
                IntAttribute* int_attr;
                DoubleAttribute* double_attr;
                InheritedAttribute* inh_attr;

                if ( ( int_attr = dynamic_cast < IntAttribute* > ( attr ) ) != 0 )
                {
                    // int attribute
                    data_array->SetAttribute ( "type", "Int32" );
                    for ( EntityIterator entity_it = mesh.begin ( it->second );
                          entity_it != mesh.end ( it->second ); ++entity_it )
                    {
                        os << int_attr->get_int_value ( entity_it->index ( ) ) << " ";
                    }
                    TiXmlText* data = new TiXmlText ( os.str ( ) );
                    data_array->LinkEndChild ( data );

                }
                else if ( ( double_attr = dynamic_cast < DoubleAttribute* > ( attr ) ) != 0 )
                {
                    // double attribute
                    data_array->SetAttribute ( "type", "Float64" );
                    for ( EntityIterator entity_it = mesh.begin ( it->second );
                          entity_it != mesh.end ( it->second ); ++entity_it )
                    {
                        os << double_attr->get_double_value ( entity_it->index ( ) ) << " ";

                    }
                    TiXmlText* data = new TiXmlText ( os.str ( ) );
                    data_array->LinkEndChild ( data );

                }
                else if ( ( inh_attr = dynamic_cast < InheritedAttribute* > ( attr ) ) != 0 )
                {
                    // TODO (Staffan): this code duplication was done under extreme stress. Clean up!
                    Attribute* data_attr = inh_attr->data_attribute ( ).get ( );
                    if ( ( int_attr = dynamic_cast < IntAttribute* > ( data_attr ) ) != 0 )
                    {
                        // int attribute
                        data_array->SetAttribute ( "type", "Int32" );
                        for ( EntityIterator entity_it = mesh.begin ( it->second );
                              entity_it != mesh.end ( it->second ); ++entity_it )
                        {
                            os << inh_attr->get_int_value ( entity_it->index ( ) ) << " ";
                        }
                        TiXmlText* data = new TiXmlText ( os.str ( ) );
                        data_array->LinkEndChild ( data );
                    }
                    else if ( ( double_attr = dynamic_cast < DoubleAttribute* > ( data_attr ) ) != 0 )
                    {
                        // double attribute
                        data_array->SetAttribute ( "type", "Float64" );
                        for ( EntityIterator entity_it = mesh.begin ( it->second );
                              entity_it != mesh.end ( it->second ); ++entity_it )
                        {
                            os << inh_attr->get_double_value ( entity_it->index ( ) ) << " ";

                        }
                        TiXmlText* data = new TiXmlText ( os.str ( ) );
                        data_array->LinkEndChild ( data );
                    }
                    else
                    {
                        // TODO Any attribute
                        LOG_ERROR ( "AnyAttributes are currently not supported" );
                        continue;
                    }
                }
                else
                {
                    // TODO Any attribute
                    LOG_ERROR ( "AnyAttributes are currently not supported" );
                    continue;
                }

                data_array->SetAttribute ( "format", "ascii" );

                // set data for the right tdim
                if ( it->second == 0 ) // -> is point_data
                    point_data->LinkEndChild ( data_array );
                if ( it->second == tdim ) // -> is cell_data
                    cell_data->LinkEndChild ( data_array );
                os.str ( "" );
                os.clear ( );

            }

            // write material ids
            data_array = new TiXmlElement ( "DataArray" );
            data_array->SetAttribute ( "type", "Int32" );
            data_array->SetAttribute ( "Name", "Material Id" );
            data_array->SetAttribute ( "format", "ascii" );

            MaterialNumber min_mat_num = 0; // mesh.begin(tdim)->get_material_number();
            MaterialNumber max_mat_num = 1000; // mesh.begin(tdim)->get_material_number();
            for ( EntityIterator it = mesh.begin ( tdim ); it != mesh.end ( tdim ); ++it )
            {
                // TODO(Thomas): Get real material number and determine
                // min and max range.

                const MaterialNumber mat_num = it->get_material_number ( );
                os << mat_num << " ";

                if ( mat_num < min_mat_num )
                {
                    min_mat_num = mat_num;
                }

                if ( mat_num > max_mat_num )
                {
                    max_mat_num = mat_num;
                }
            }
            data_array->SetAttribute ( "RangeMin", min_mat_num );
            data_array->SetAttribute ( "RangeMax", max_mat_num );

            cell_data->LinkEndChild ( data_array );
            TiXmlText* material_id = new TiXmlText ( os.str ( ) );
            data_array->LinkEndChild ( material_id );
            os.str ( "" );
            os.clear ( );

            // write points
            TiXmlElement* points = new TiXmlElement ( "Points" );
            piece->LinkEndChild ( points );
            data_array = new TiXmlElement ( "DataArray" );
            data_array->SetAttribute ( "type", "Float64" );
            data_array->SetAttribute ( "Name", "Array" );
            // always 3 comps, since vtk doesn:t handle 2D.
            data_array->SetAttribute ( "NumberOfComponents", "3" );
            data_array->SetAttribute ( "format", "ascii" );

            // output vertex coordinates separated by spaces
            std::map<Id, Id> vtk_id_map;
            Id index = 0;
            double min = 0.0, max = 0.0;

            for ( EntityIterator it = mesh.begin ( 0 ); it != mesh.end ( 0 ); ++it )
            {
                std::vector<Coordinate> vertex;
                it->get_coordinates ( vertex, 0 );
                // vtk doesn:t like 2D and 1D coordinates
                if ( mesh.gdim ( ) < 3 )
                {
                    vertex.push_back ( 0 );
                    if ( mesh.gdim ( ) < 2 )
                    {
                        vertex.push_back ( 0 );
                    }
                }

                if ( deformation_ )
                {
                    for ( int c = 0; c < 3; ++c )
                    {
                        vertex[c] += ( *deformation_ )[3 * it->index ( ) + c];
                    }
                }

                copy ( vertex.begin ( ), vertex.end ( ), std::ostream_iterator<double>( os, " " ) );

                // build up id map to have connectivities zero-based.
                vtk_id_map.insert ( std::pair<Id, Id>( it->id ( ), index ) );
                ++index;

                // keep largest and smallest value
                min = *min_element ( vertex.begin ( ), vertex.end ( ) ) < min ? *min_element ( vertex.begin ( ), vertex.end ( ) ) : min;
                max = *max_element ( vertex.begin ( ), vertex.end ( ) ) > max ? *max_element ( vertex.begin ( ), vertex.end ( ) ) : max;
            }

            TiXmlText* coords = new TiXmlText ( os.str ( ) );
            data_array->SetAttribute ( "RangeMin", min );
            data_array->SetAttribute ( "RangeMax", max );

            points->LinkEndChild ( data_array );
            data_array->LinkEndChild ( coords );
            os.str ( "" );
            os.clear ( );

            // write cells
            TiXmlElement* cells = new TiXmlElement ( "Cells" );
            piece->LinkEndChild ( cells );

            // connectivities
            Id id_min = 0, id_max = 0;
            data_array = new TiXmlElement ( "DataArray" );
            data_array->SetAttribute ( "type", "Int64" );
            data_array->SetAttribute ( "Name", "connectivity" );
            data_array->SetAttribute ( "format", "ascii" );
            for ( EntityIterator it = mesh.begin ( tdim ); it != mesh.end ( tdim ); ++it )
            {
                for ( VertexIdIterator iit = it->begin_vertex_ids ( ); iit != it->end_vertex_ids ( ); ++iit )
                {
                    os << vtk_id_map[*iit] << " ";

                    // keep largest and smallest Id
                    id_min = vtk_id_map[*iit] < id_min ? vtk_id_map[*iit] : id_min;
                    id_max = vtk_id_map[*iit] > id_max ? vtk_id_map[*iit] : id_max;
                }
            }
            TiXmlText* conns = new TiXmlText ( os.str ( ) );
            data_array->SetAttribute ( "RangeMin", id_min );
            data_array->SetAttribute ( "RangeMax", id_max );

            cells->LinkEndChild ( data_array );
            data_array->LinkEndChild ( conns );
            os.str ( "" );
            os.clear ( );

            // offsets
            data_array = new TiXmlElement ( "DataArray" );
            data_array->SetAttribute ( "type", "Int64" );
            data_array->SetAttribute ( "Name", "offsets" );
            data_array->SetAttribute ( "format", "ascii" );

            EntityCount offset = 0;
            EntityCount off_min = 0, off_max = 0;
            for ( EntityIterator it = mesh.begin ( tdim ); it != mesh.end ( tdim ); ++it )
            {
                offset += it->num_vertices ( );
                os << offset << " ";

                // keep largest and smallest offset
                off_min = offset < off_min ? offset : off_min;
                off_max = offset > off_max ? offset : off_max;
            }
            data_array->SetAttribute ( "RangeMin", off_min );
            data_array->SetAttribute ( "RangeMax", off_max );

            cells->LinkEndChild ( data_array );
            TiXmlText* offs = new TiXmlText ( os.str ( ) );
            data_array->LinkEndChild ( offs );
            os.str ( "" );
            os.clear ( );

            // types
            data_array = new TiXmlElement ( "DataArray" );
            data_array->SetAttribute ( "type", "UInt8" );
            data_array->SetAttribute ( "Name", "types" );
            data_array->SetAttribute ( "format", "ascii" );

            int typ_min = 0, typ_max = 0;
            for ( EntityIterator it = mesh.begin ( tdim ); it != mesh.end ( tdim ); ++it )
            {
                const CellType& type = it->cell_type ( );
                os << vtk_cell_types[static_cast < int > ( type.tag ( ) )] << " ";

                // keep largest and smallest type
                typ_min = vtk_cell_types[static_cast < int > ( type.tag ( ) )] < typ_min ? vtk_cell_types[static_cast < int > ( type.tag ( ) )] : typ_min;
                typ_max = vtk_cell_types[static_cast < int > ( type.tag ( ) )] > typ_max ? vtk_cell_types[static_cast < int > ( type.tag ( ) )] : typ_max;
            }
            data_array->SetAttribute ( "RangeMin", typ_min );
            data_array->SetAttribute ( "RangeMax", typ_max );

            cells->LinkEndChild ( data_array );
            TiXmlText* types = new TiXmlText ( os.str ( ) );

            data_array->LinkEndChild ( types );

            os.str ( "" );
            os.clear ( );

            FILE * pFile;
            pFile = fopen ( filename, "w" );
            if ( pFile != NULL )
            {
                doc.SaveFile ( pFile );
                fclose ( pFile );
            }
            else
            {
                std::stringstream err;
                err << "Path to write the files (" << filename << ") does not exist!";
                LOG_ERROR ( err.str ( ) );
                throw std::runtime_error ( err.str ( ) );
            }
        }

        //////// PVtkWriter ////////////////////////

        void PVtkWriter::write_file ( const char* filename, const Mesh & mesh ) const
        {
            // get MPI rank
            int rank = -1, num_procs = -1;
            MPI_Comm_rank ( mpi_comm_, &rank );
            MPI_Comm_size ( mpi_comm_, &num_procs );

            std::stringstream s;
            s << rank;

            // get the correct filename including the path
            std::istringstream filename_root_dir ( filename );
            std::size_t dir = filename_root_dir.str ( ).find_last_of ( "." );
            LOG_DEBUG ( 3, "Filename: " << filename );

            std::string filename_without_suffix = filename_root_dir.str ( ).substr ( 0, dir );
            assert ( !filename_without_suffix.empty ( ) );

            std::string str_src_filename = ( filename_without_suffix + "_" + s.str ( ) + ".vtu" );
            LOG_DEBUG ( 3, "Filename without suffix: " << filename_without_suffix );
            assert ( !str_src_filename.empty ( ) );

            // write all vtu files
            VtkWriter writer;
            writer.set_deformation ( deformation_ );
            for ( DataVectorPair::const_iterator it = attributes_to_write_.begin ( );
                  it != attributes_to_write_.end ( ); ++it )
            {
                // transfer attributes to write to local writer
                writer.add_attribute ( it->first, it->second );
            }
            writer.write ( str_src_filename.c_str ( ), mesh );

            // write pvtu file
            if ( rank == 0 )
            {
                TiXmlDocument doc;
                TiXmlDeclaration* decl = new TiXmlDeclaration ( "1.0", "", "" );
                doc.LinkEndChild ( decl );

                TiXmlElement* vtkFile = new TiXmlElement ( "VTKFile" );
                vtkFile->SetAttribute ( "type", "PUnstructuredGrid" );
                vtkFile->SetAttribute ( "version", "0.1" );
                vtkFile->SetAttribute ( "byte_order", "LittleEndian" );
                vtkFile->SetAttribute ( "compressor", "vtkZLibDataCompressor" );
                doc.LinkEndChild ( vtkFile );

                TiXmlElement* pumesh = new TiXmlElement ( "PUnstructuredGrid" );
                // TODO(Thomas): get the correct GhostLevel
                // GhostLevel in PUnstructuredGrid is always 0
                pumesh->SetAttribute ( "GhostLevel", 0 );
                vtkFile->LinkEndChild ( pumesh );

                TiXmlElement* p_point_data = new TiXmlElement ( "PPointData" );
                p_point_data->SetAttribute ( "Scalars", "Material Id" );
                pumesh->LinkEndChild ( p_point_data );

                TiXmlElement* p_cell_data = new TiXmlElement ( "PCellData" );
                p_cell_data->SetAttribute ( "Scalars", "Material Id" );
                pumesh->LinkEndChild ( p_cell_data );

                std::vector< std::string > attr_names = mesh.get_attribute_names ( mesh.tdim ( ) );
                LOG_INFO ( "writer", "Mesh cell attributes = " << string_from_range ( attr_names.begin ( ), attr_names.end ( ) ) );

                // write attributes
                TiXmlElement* data_array;
                for ( DataVectorPair::const_iterator it = attributes_to_write_.begin ( );
                      it != attributes_to_write_.end ( ); ++it )
                {

                    if ( !mesh.has_attribute ( it->first, it->second ) )
                    {
                        LOG_ERROR ( "No attribute " << it->first << " for dimension " << it->second << " found in mesh" );
                        continue;
                    }

                    if ( ( it->second != 0 ) && ( it->second != mesh.tdim ( ) ) )
                    {
                        LOG_ERROR ( "VTK does not support data for entities other than cells or vertices." );
                        continue;
                    }

                    data_array = new TiXmlElement ( "PDataArray" );
                    data_array->SetAttribute ( "Name", it->first );

                    AttributePtr attr_ptr = mesh.get_attribute ( it->first, it->second );
                    Attribute* attr = attr_ptr.get ( );
                    IntAttribute* int_attr;
                    DoubleAttribute* double_attr;
                    InheritedAttribute* inh_attr;

                    if ( ( inh_attr = dynamic_cast < InheritedAttribute* > ( attr ) ) != 0 )
                    {
                        attr = inh_attr->data_attribute ( ).get ( );
                    }

                    if ( ( int_attr = dynamic_cast < IntAttribute* > ( attr ) ) != 0 )
                    {
                        // int attribute
                        data_array->SetAttribute ( "type", "Int32" );

                    }
                    else if ( ( double_attr = dynamic_cast < DoubleAttribute* > ( attr ) ) != 0 )
                    {
                        // double attribute
                        data_array->SetAttribute ( "type", "Float64" );
                    }
                    else
                    {
                        // TODO(Thomas) Any attribute
                        LOG_ERROR ( "AnyAttributes are currently not supported" );
                        continue;
                    }

                    data_array->SetAttribute ( "format", "ascii" );

                    // set data for the right tdim
                    if ( it->second == 0 )
                    { // -> is point_data
                        p_point_data->LinkEndChild ( data_array );
                        // TODO(Thomas): Get type (Scalars, Vectors,
                        // Normals, ...)
                        p_point_data->SetAttribute ( "Scalars", it->first );
                    }
                    if ( it->second == mesh.tdim ( ) )
                    { // -> is cell_data
                        p_cell_data->LinkEndChild ( data_array );
                        // TODO(Thomas): Get type (Scalars, Vectors,
                        // Normals, ...)
                        p_cell_data->SetAttribute ( "Scalars", it->first );
                    }
                }

                // NB: This has to be AFTER the the other elements, since
                // the same order in the vtu and pvtu file is needed!
                TiXmlElement* p_cell_data_array = new TiXmlElement ( "PDataArray" );
                p_cell_data_array->SetAttribute ( "Name", "Material Id" );
                p_cell_data_array->SetAttribute ( "type", "Int32" );
                p_cell_data->LinkEndChild ( p_cell_data_array );

                TiXmlElement* p_points = new TiXmlElement ( "PPoints" );
                pumesh->LinkEndChild ( p_points );

                TiXmlElement* p_points_data_array = new TiXmlElement ( "PDataArray" );
                p_points_data_array->SetAttribute ( "type", "Float64" );
                p_points_data_array->SetAttribute ( "NumberOfComponents", "3" );
                p_points->LinkEndChild ( p_points_data_array );

                // get the correct filename without the path
                std::size_t pos = filename_root_dir.str ( ).find_last_of ( "/\\" );
                assert ( !filename_root_dir.str ( ).substr ( pos + 1, filename_root_dir.str ( ).length ( ) ).empty ( ) );

                std::stringstream str_proc_id;

                for ( int proc_id = 0; proc_id < num_procs; ++proc_id )
                {
                    TiXmlElement* piece = new TiXmlElement ( "Piece" ); // needs to be inside the loop!
                    str_proc_id << proc_id;
                    piece->SetAttribute ( "Source", filename_root_dir.str ( ).substr ( pos + 1, dir - pos - 1 ) + "_" + str_proc_id.str ( ) + ".vtu" );
                    pumesh->LinkEndChild ( piece );
                    str_proc_id.str ( "" );
                    str_proc_id.clear ( );
                }

                // TODO(Thomas): Look if file exists, if not, what to do?
                FILE * pFile;
                pFile = fopen ( filename, "w" );
                if ( pFile != NULL )
                {
                    doc.SaveFile ( pFile );
                    fclose ( pFile );
                }
                else
                {
                    std::stringstream err;
                    err << "Path to write the files (" << filename << ") does not exist!";
                    LOG_ERROR ( err.str ( ) );
                    throw std::runtime_error ( err.str ( ) );
                }
            }
        }
    }
} // namespace hiflow
