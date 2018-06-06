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

/// \author Teresa Beck

#ifndef _MY_CELL_VISUALIZATION_H_
#    define _MY_CELL_VISUALIZATION_H_

#    include <map>
#    include <vector>
#    include <string>

#    include "space/cell_visualization.h"
#    include "common/log.h"
#    include "fem/cell_transformation.h"
#    include "mesh/entity.h"
#    include "mesh/iterator.h"
#    include "space/vector_space.h"
#    include <cmath>
#    include <limits>
#    include <sstream>
#    include "../../contrib/tinyxml/tinyxml.h"

namespace hiflow
{
    namespace mesh
    {

        class MyCellVisualization : public ParallelCellVisualization<double>
        {
            /// \brief CellVisualization to be used for periodic domains
            /// \see CellVisualization
          private:
            //      const VectorSpace<double>& space_;
            //      CellVisualizationGrids<double> grids_;
            std::map< std::string, std::vector<double> > functions_;
            std::map< std::string, int > num_components_;
            double scale_factor_z_;
          public:

            MyCellVisualization ( const VectorSpace<double>& space,
                                  int num_intervals,
                                  const MPI_Comm& comm,
                                  const int master_rank,
                                  double scale_factor_z = 1. )
            : ParallelCellVisualization<double>( space, num_intervals, comm, master_rank )
            {
                scale_factor_z_ = scale_factor_z;
            };
            void write_sequential ( const std::string& filename ) const;
            void write ( const std::string& filename ) const;
            void visualize ( const EvalFunction& fun, const std::string& name, int num_components );

            void visualize_cell_data ( const std::vector<double>& cell_data, const std::string& name )
            {
                return CellVisualization::visualize_cell_data ( cell_data, name );
            };
        };
    };
};

inline void MyCellVisualization::write ( const std::string& filename ) const
{
    const mesh::Mesh& mesh = space_.mesh ( );
    const mesh::TDim tdim = mesh.tdim ( );
    const mesh::GDim gdim = mesh.gdim ( );
    // get MPI rank
    int rank = -1, num_procs = -1;
    MPI_Comm_rank ( mpi_comm_, &rank );
    MPI_Comm_size ( mpi_comm_, &num_procs );

    std::stringstream s;
    s << rank;

    // get the correct filename including the path
    std::istringstream filename_root_dir ( filename );

    std::size_t dir = filename_root_dir.str ( ).find_last_of ( "." );

    std::string filename_without_suffix = filename_root_dir.str ( ).substr ( 0, dir );

    assert ( !filename_without_suffix.empty ( ) );

    std::string str_src_filename = ( filename_without_suffix + "_" + s.str ( ) + ".vtu" );
    assert ( !str_src_filename.empty ( ) );

    // Each process writes its vtu file
    write_sequential ( str_src_filename );

    // Master writes pvtu file
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
        // GhostLevel in PUnstructuredGrid is always 0
        pumesh->SetAttribute ( "GhostLevel", 0 );
        vtkFile->LinkEndChild ( pumesh );

        TiXmlElement* p_point_data = new TiXmlElement ( "PPointData" );
        pumesh->LinkEndChild ( p_point_data );

        TiXmlElement* p_data_array;
        for ( std::map<std::string, std::vector<double> >::const_iterator it = functions_.begin ( ),
              end_it = functions_.end ( ); it != end_it; ++it )
        {

            p_data_array = new TiXmlElement ( "PDataArray" );
            p_data_array->SetAttribute ( "Name", it->first );
            p_data_array->SetAttribute ( "type", "Float64" );
            const int nc = num_components_.find ( it->first )->second;
            p_data_array->SetAttribute ( "NumberOfComponents", nc );
            p_data_array->SetAttribute ( "format", "ascii" );
            p_point_data->LinkEndChild ( p_data_array );
        }

        TiXmlElement* p_cell_data = new TiXmlElement ( "PCellData" );
        pumesh->LinkEndChild ( p_cell_data );
        int tdim = mesh.tdim ( );

        // write cell data
        // TODO currently only Float64 is supported
        TiXmlElement* data_array;
        for ( std::map< std::string, std::vector<double> >::const_iterator it = functions_cell_.begin ( );
              it != functions_cell_.end ( ); ++it )
        {
            data_array = new TiXmlElement ( "PDataArray" );
            data_array->SetAttribute ( "Name", it->first );
            data_array->SetAttribute ( "type", "Float64" );
            data_array->SetAttribute ( "format", "ascii" );
            p_cell_data->LinkEndChild ( data_array );
            p_cell_data->SetAttribute ( "Scalars", it->first );
        }

        // NB: This has to be AFTER the the other elements, since
        // the same order in the vtu and pvtu file is needed!

        TiXmlElement* p_points = new TiXmlElement ( "PPoints" );
        pumesh->LinkEndChild ( p_points );

        TiXmlElement* p_points_data_array = new TiXmlElement ( "PDataArray" );
        p_points_data_array->SetAttribute ( "type", "Float64" );
        p_points_data_array->SetAttribute ( "NumberOfComponents", 3 );
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

        std::string str_filename = ( filename_without_suffix + ".pvtu" );
        FILE * pFile;
        pFile = fopen ( str_filename.c_str ( ), "w" );
        if ( pFile != NULL )
        {
            doc.SaveFile ( pFile );
            fclose ( pFile );
        }
        else
        {
            std::stringstream err;
            err << "Path to write the files (" << str_filename << ") does not exist!";
            LOG_ERROR ( err.str ( ) );
            throw std::runtime_error ( err.str ( ) );
        }
    }
}

inline void MyCellVisualization::write_sequential ( const std::string& filename ) const
{
    const mesh::Mesh& mesh = space_.mesh ( );
    const mesh::TDim tdim = mesh.tdim ( );
    const mesh::GDim gdim = mesh.gdim ( );
    int num_mesh_cells = 0;
    for ( mesh::EntityIterator it = mesh.begin ( tdim ), end_it = mesh.end ( tdim );
          it != end_it; ++it )
    {

        if ( parallel_visualization_ )
        {
            int rem_ind = -100;
            it->get<int>( "_remote_index_", &rem_ind );
            if ( rem_ind != -1 ) continue;
        }
        ++num_mesh_cells;
    }

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
    piece->SetAttribute ( "NumberOfPoints", grids_.num_visu_points ( ) );
    piece->SetAttribute ( "NumberOfCells", grids_.num_visu_cells ( ) );
    umesh->LinkEndChild ( piece );

    // Scratch variables for data.
    TiXmlElement* data_array;
    std::stringstream os;

    //// Write points ////////////////////////////////////////////////////////////
    TiXmlElement* points = new TiXmlElement ( "Points" );
    piece->LinkEndChild ( points );
    data_array = new TiXmlElement ( "DataArray" );

    // Set correct length of float in dependence of double
    std::ostringstream type_float;
    type_float << "Float" << sizeof (double ) * 8;
    data_array->SetAttribute ( "type", type_float.str ( ) );

    data_array->SetAttribute ( "Name", "Array" );
    // always 3 comps, since vtk doesn:t handle 2D.
    data_array->SetAttribute ( "NumberOfComponents", "3" );
    data_array->SetAttribute ( "format", "ascii" );

    std::vector<double> ref_pt ( gdim, 0. ), mapped_pt ( 3, 0. );
    double range_min = std::numeric_limits<double>::max ( );
    double range_max = std::numeric_limits<double>::min ( );
    int cell_type_min = std::numeric_limits<int>::max ( );
    int cell_type_max = std::numeric_limits<int>::min ( );

    for ( mesh::EntityIterator it = mesh.begin ( tdim ), end_it = mesh.end ( tdim );
          it != end_it; ++it )
    {

        if ( parallel_visualization_ )
        {
            int rem_ind = -100;
            it->get<int>( "_remote_index_", &rem_ind );
            if ( rem_ind != -1 ) continue;
        }
        const doffem::CellTransformation<double>& cell_trans = space_.GetCellTransformation ( *it );

        // TODO: this can be done with a single call now
        //      for (int p = 0, p_end = grid_.num_points(); p != p_end; ++p) {
        for ( int p = 0, p_end = grids_.num_points ( it->cell_type ( ).tag ( ) ); p != p_end; ++p )
        {
            const int offset = gdim * p;
            for ( int c = 0; c < gdim; ++c )
            {
                ref_pt[c] = grids_.coords ( it->cell_type ( ).tag ( ) )[offset + c];
                mapped_pt[c] = 0.; // reset
            }

            // TODO: Would be great if we could do this with one call, as
            // with Hp-transformations.
            if ( gdim >= 1 ) mapped_pt[0] = cell_trans.x ( ref_pt );
            if ( gdim >= 2 ) mapped_pt[1] = cell_trans.y ( ref_pt );
            if ( gdim >= 3 ) mapped_pt[2] = cell_trans.z ( ref_pt );

            // z-scaling
            if ( gdim >= 3 ) mapped_pt[2] *= scale_factor_z_;

            for ( int c = 0; c < 3; ++c )
            {
                range_min = std::min ( range_min, mapped_pt[c] );
                range_max = std::max ( range_max, mapped_pt[c] );
                os << mapped_pt[c] << " ";
            }
        }

    }

    TiXmlText* coords = new TiXmlText ( os.str ( ) );
    data_array->SetAttribute ( "RangeMin", range_min );
    data_array->SetAttribute ( "RangeMax", range_max );

    points->LinkEndChild ( data_array );
    data_array->LinkEndChild ( coords );
    os.str ( "" );
    os.clear ( );
    //// End write points /////////////////////////////////////////////////////////

    //// Write cells //////////////////////////////////////////////////////////////
    TiXmlElement* cells = new TiXmlElement ( "Cells" );
    piece->LinkEndChild ( cells );

    int p_offset = 0, cell_offset = 0;

    // Connectivity, Offsets, and Types arrays
    std::ostringstream off_os, type_os;
    static const int vtk_cell_types[] = { 1, 3, 5, 9, 10, 12 };

    for ( mesh::EntityIterator it = mesh.begin ( tdim ); it != mesh.end ( tdim ); ++it )
    {
        if ( parallel_visualization_ )
        {
            int rem_ind = -100;
            it->get<int>( "_remote_index_", &rem_ind );
            if ( rem_ind != -1 ) continue;
        }

        for ( int c = 0, c_end = grids_.num_cells ( it->cell_type ( ).tag ( ) ); c != c_end; ++c )
        {
            const std::vector<int>& verts = grids_.vertices_of_cell ( it->cell_type ( ).tag ( ), c );
            for ( int v = 0, v_end = verts.size ( ); v != v_end; ++v )
            {
                os << verts[v] + p_offset << " ";
            }
            cell_offset += verts.size ( );
            off_os << cell_offset << " ";
            type_os << vtk_cell_types[static_cast < int > ( it->cell_type ( ).tag ( ) )] << " ";
            cell_type_min = std::min ( cell_type_min, vtk_cell_types[static_cast < int > ( it->cell_type ( ).tag ( ) )] );
            cell_type_max = std::max ( cell_type_max, vtk_cell_types[static_cast < int > ( it->cell_type ( ).tag ( ) )] );

        }
        p_offset += grids_.num_points ( it->cell_type ( ).tag ( ) );
    }

    data_array = new TiXmlElement ( "DataArray" );
    data_array->SetAttribute ( "type", "Int64" );
    data_array->SetAttribute ( "Name", "connectivity" );
    data_array->SetAttribute ( "format", "ascii" );
    data_array->SetAttribute ( "RangeMin", 0 );
    data_array->SetAttribute ( "RangeMax", grids_.num_visu_points ( ) );

    TiXmlText* conns = new TiXmlText ( os.str ( ) );
    data_array->LinkEndChild ( conns );
    cells->LinkEndChild ( data_array );
    os.str ( "" );
    os.clear ( );

    data_array = new TiXmlElement ( "DataArray" );
    data_array->SetAttribute ( "type", "Int64" );
    data_array->SetAttribute ( "Name", "offsets" );
    data_array->SetAttribute ( "format", "ascii" );
    data_array->SetAttribute ( "RangeMin", 0 );
    data_array->SetAttribute ( "RangeMax", cell_offset );

    TiXmlText* offs = new TiXmlText ( off_os.str ( ) );
    data_array->LinkEndChild ( offs );
    cells->LinkEndChild ( data_array );
    off_os.str ( "" );
    off_os.clear ( );

    data_array = new TiXmlElement ( "DataArray" );
    data_array->SetAttribute ( "type", "UInt8" );
    data_array->SetAttribute ( "Name", "types" );
    data_array->SetAttribute ( "format", "ascii" );
    data_array->SetAttribute ( "RangeMin", cell_type_min );
    data_array->SetAttribute ( "RangeMax", cell_type_max );

    TiXmlText* types = new TiXmlText ( type_os.str ( ) );
    data_array->LinkEndChild ( types );
    cells->LinkEndChild ( data_array );
    type_os.str ( "" );
    type_os.clear ( );

    //// End Write cells //////////////////////////////////////////////////////////

    //// Write point data /////////////////////////////////////////////////////////
    TiXmlElement* point_data = new TiXmlElement ( "PointData" );
    piece->LinkEndChild ( point_data );

    for ( typename std::map<std::string, std::vector<double> >::const_iterator it = functions_.begin ( ),
          end_it = functions_.end ( ); it != end_it; ++it )
    {
        data_array = new TiXmlElement ( "DataArray" );
        data_array->SetAttribute ( "Name", it->first );
        data_array->SetAttribute ( "type", type_float.str ( ) );
        data_array->SetAttribute ( "format", "ascii" );

        for ( int i = 0, end_i = it->second.size ( ); i != end_i; ++i )
        {
            os << ( it->second )[i] << " ";
        }

        TiXmlText* data = new TiXmlText ( os.str ( ) );
        data_array->LinkEndChild ( data );
        point_data->LinkEndChild ( data_array );
        os.str ( "" );
        os.clear ( );
    }

    TiXmlElement* cell_data = new TiXmlElement ( "CellData" );
    piece->LinkEndChild ( cell_data );

    for ( typename std::map<std::string, std::vector<double> >::const_iterator it = functions_cell_.begin ( ),
          end_it = functions_cell_.end ( ); it != end_it; ++it )
    {
        data_array = new TiXmlElement ( "DataArray" );
        data_array->SetAttribute ( "Name", it->first );
        data_array->SetAttribute ( "type", type_float.str ( ) );
        data_array->SetAttribute ( "format", "ascii" );

        for ( int i = 0, end_i = it->second.size ( ); i != end_i; ++i )
        {
            os << ( it->second )[i] << " ";
        }

        TiXmlText* data = new TiXmlText ( os.str ( ) );
        data_array->LinkEndChild ( data );
        cell_data->LinkEndChild ( data_array );
        os.str ( "" );
        os.clear ( );
    }

    doc.SaveFile ( filename );

    //   const mesh::Mesh& mesh = space_.mesh();
    //   const mesh::TDim tdim = mesh.tdim();
    //   const mesh::GDim gdim = mesh.gdim();
    //
    //
    //   int num_mesh_cells =0;
    //
    //   for (mesh::EntityIterator it = mesh.begin(tdim), end_it = mesh.end(tdim);
    //        it != end_it; ++it) {
    //
    //     if (parallel_visualization_)
    //       {
    //     int rem_ind = -100;
    //     it->get<int>("_remote_index_", &rem_ind);
    //     if (rem_ind != -1) continue;
    //     num_mesh_cells++;
    //       }
    //   }
    //
    //   TiXmlDocument doc;
    //   TiXmlDeclaration* decl = new TiXmlDeclaration("1.0", "", "");
    //   doc.LinkEndChild(decl);
    //
    //   TiXmlElement* vtkFile = new TiXmlElement("VTKFile");
    //   vtkFile->SetAttribute("type", "UnstructuredGrid");
    //   vtkFile->SetAttribute("version", "0.1");
    //   vtkFile->SetAttribute("byte_order", "LittleEndian");
    //   vtkFile->SetAttribute("compressor", "vtkZLibDataCompressor");
    //   doc.LinkEndChild(vtkFile);
    //
    //   TiXmlElement* umesh = new TiXmlElement("UnstructuredGrid");
    //   vtkFile->LinkEndChild(umesh);
    //
    //   TiXmlElement* piece = new TiXmlElement("Piece");
    //   piece->SetAttribute("NumberOfPoints", num_mesh_cells * grids_.num_points());
    //   piece->SetAttribute("NumberOfCells", num_mesh_cells * grids_.num_cells());
    //   umesh->LinkEndChild(piece);
    //
    //   // Scratch variables for data.
    //   TiXmlElement* data_array;
    //   std::stringstream os;
    //
    //   //// Write points ////////////////////////////////////////////////////////////
    //   TiXmlElement* points = new TiXmlElement("Points");
    //   piece->LinkEndChild(points);
    //   data_array = new TiXmlElement("DataArray");
    //   data_array->SetAttribute("type", "Float64");
    //   data_array->SetAttribute("Name", "Array");
    //   // always 3 comps, since vtk doesn:t handle 2D.
    //   data_array->SetAttribute("NumberOfComponents", "3");
    //   data_array->SetAttribute("format", "ascii");
    //
    //   std::vector<double> ref_pt(gdim, 0.), mapped_pt(3, 0.);
    //   double range_min = std::numeric_limits<double>::max();
    //   double range_max = std::numeric_limits<double>::min();
    //
    //   for (mesh::EntityIterator it = mesh.begin(tdim), end_it = mesh.end(tdim);
    //        it != end_it; ++it) {
    //     if (parallel_visualization_)
    //       {
    //     int rem_ind = -100;
    //     it->get<int>("_remote_index_", &rem_ind);
    //     if (rem_ind != -1) continue;
    //       }
    //     const doffem::CellTransformation& cell_trans = space_.GetCellTransformation(*it);
    //
    //     // TODO: this can be done with a single call now
    //     for (int p = 0, p_end = grids_.num_points(); p != p_end; ++p) {
    //       for (int c = 0; c < gdim; ++c) {
    //         ref_pt[c] = grids_.coords()[gdim * p + c];
    //         mapped_pt[c] = 0.;  // reset
    //       }
    //
    //       // TODO: Would be great if we could do this with one call, as
    //       // with Hp-transformations.
    //       if (gdim >= 1) mapped_pt[0] = cell_trans.x(ref_pt);
    //       if (gdim >= 2) mapped_pt[1] = cell_trans.y(ref_pt);
    //       if (gdim >= 3) mapped_pt[2] = cell_trans.z(ref_pt);
    //
    //       // z-scaling
    //       if (gdim >= 3) mapped_pt[2] *= scale_factor_z_;
    //
    //
    //       for (int c = 0; c < 3; ++c) {
    //         range_min = std::min(range_min, mapped_pt[c]);
    //         range_max = std::max(range_max, mapped_pt[c]);
    //         os << mapped_pt[c] << " ";
    //       }
    //     }
    //
    //   }
    //
    //   TiXmlText* coords = new TiXmlText(os.str());
    //   data_array->SetAttribute("RangeMin", range_min);
    //   data_array->SetAttribute("RangeMax", range_max);
    //
    //   points->LinkEndChild(data_array);
    //   data_array->LinkEndChild(coords);
    //   os.str("");
    //   os.clear();
    //   //// End write points /////////////////////////////////////////////////////////
    //
    //   //// Write cells //////////////////////////////////////////////////////////////
    //   TiXmlElement* cells = new TiXmlElement("Cells");
    //   piece->LinkEndChild(cells);
    //
    //
    //   int p_offset = 0, cell_offset = 0;
    //
    //   // Connectivity, Offsets, and Types arrays
    //   std::ostringstream off_os, type_os;
    //   assert(grids_.gdim() == 2 || grids_.gdim() == 3);
    //   const int cell_type = (grids_.gdim() == 2) ? 9 : 12;      // TODO: handle other types than quads
    //
    //   for (mesh::EntityIterator it = mesh.begin(tdim); it != mesh.end(tdim); ++it) {
    //
    //     if (parallel_visualization_)
    //       {
    //     int rem_ind = -100;
    //     it->get<int>("_remote_index_", &rem_ind);
    //     if (rem_ind != -1) continue;
    //       }
    //     for (int c = 0, c_end = grids_.num_cells(); c != c_end; ++c) {
    //       const std::vector<int>& verts = grids_.vertices_of_cell(c);
    //       for (int v = 0, v_end = verts.size(); v != v_end; ++v) {
    //         os << verts[v] + p_offset << " ";
    //       }
    //       cell_offset += verts.size();
    //       off_os << cell_offset << " ";
    //       type_os << cell_type << " ";
    //     }
    //     p_offset += grids_.num_points();
    //   }
    //
    //   data_array = new TiXmlElement("DataArray");
    //   data_array->SetAttribute("type", "Int64");
    //   data_array->SetAttribute("Name", "connectivity");
    //   data_array->SetAttribute("format", "ascii");
    //   data_array->SetAttribute("RangeMin", 0);
    //   data_array->SetAttribute("RangeMax", num_mesh_cells * grids_.num_points());
    //
    //   TiXmlText* conns = new TiXmlText(os.str());
    //   data_array->LinkEndChild(conns);
    //   cells->LinkEndChild(data_array);
    //   os.str("");
    //   os.clear();
    //
    //   data_array = new TiXmlElement("DataArray");
    //   data_array->SetAttribute("type", "Int64");
    //   data_array->SetAttribute("Name", "offsets");
    //   data_array->SetAttribute("format", "ascii");
    //   data_array->SetAttribute("RangeMin", 0);
    //   data_array->SetAttribute("RangeMax", cell_offset);
    //
    //   TiXmlText* offs = new TiXmlText(off_os.str());
    //   data_array->LinkEndChild(offs);
    //   cells->LinkEndChild(data_array);
    //   off_os.str("");
    //   off_os.clear();
    //
    //   data_array = new TiXmlElement("DataArray");
    //   data_array->SetAttribute("type", "UInt8");
    //   data_array->SetAttribute("Name", "types");
    //   data_array->SetAttribute("format", "ascii");
    //   data_array->SetAttribute("RangeMin", cell_type);
    //   data_array->SetAttribute("RangeMax", cell_type);
    //
    //   TiXmlText* types = new TiXmlText(type_os.str());
    //   data_array->LinkEndChild(types);
    //   cells->LinkEndChild(data_array);
    //   type_os.str("");
    //   type_os.clear();
    //
    //
    //   //// End Write cells //////////////////////////////////////////////////////////
    //
    //   //// Write point data /////////////////////////////////////////////////////////
    //
    //
    //   TiXmlElement* point_data = new TiXmlElement("PointData");
    //
    //   std::ostringstream scalar_sstr, vector_sstr;
    //   for (std::map< std::string, int >::const_iterator it = num_components_.begin(),
    //      end_it = num_components_.end(); it != end_it; ++it) {
    //     if (it->second == 1) {
    //       scalar_sstr << it->first << " ";
    //     } else {
    //       vector_sstr << it->first << " ";
    //     }
    //   }
    //   data_array->SetAttribute("Scalars", scalar_sstr.str());
    //   data_array->SetAttribute("Vectors", vector_sstr.str());
    //   piece->LinkEndChild(point_data);
    //
    //   for (std::map<std::string, std::vector<double> >::const_iterator it = functions_.begin(),
    //      end_it = functions_.end(); it != end_it; ++it) {
    //     data_array = new TiXmlElement("DataArray");
    //     data_array->SetAttribute("Name", it->first);
    //     data_array->SetAttribute("type", "Float64");
    //     data_array->SetAttribute("format", "ascii");
    //
    //     assert(num_components_.count(it->first) != 0);
    //
    //     const int nc = num_components_.find(it->first)->second;
    //     data_array->SetAttribute("NumberOfComponents", nc);
    //
    //     for (int i = 0, end_i = it->second.size(); i != end_i; ++i) {
    //       os << (it->second)[i] << " ";
    //     }
    //
    //     TiXmlText* data = new TiXmlText(os.str());
    //     data_array->LinkEndChild(data);
    //     point_data->LinkEndChild(data_array);
    //     os.str("");
    //     os.clear();
    //   }
    //
    //   TiXmlElement* cell_data = new TiXmlElement("CellData");
    //   piece->LinkEndChild(cell_data);
    //
    //   for (std::map<std::string, std::vector<double> >::const_iterator it = functions_cell_.begin(),
    //          end_it = functions_cell_.end(); it != end_it; ++it) {
    //     data_array = new TiXmlElement("DataArray");
    //     data_array->SetAttribute("Name", it->first);
    //     data_array->SetAttribute("type", "Float64");
    //     data_array->SetAttribute("format", "ascii");
    //
    //     for (int i = 0, end_i = it->second.size(); i != end_i; ++i) {
    //       os << (it->second)[i] << " ";
    //     }
    //
    //     TiXmlText* data = new TiXmlText(os.str());
    //     data_array->LinkEndChild(data);
    //     cell_data->LinkEndChild(data_array);
    //     os.str("");
    //     os.clear();
    //   }
    //  doc.SaveFile(filename);
}

inline void MyCellVisualization::visualize ( const EvalFunction& fun,
                                             const std::string& name,
                                             int num_components )
{
    const mesh::Mesh& mesh = space_.mesh ( );
    const mesh::TDim tdim = mesh.tdim ( );

    std::vector<double> values, cell_values;

    // before streamlining
    // cell_values.reserve(grids_.num_points() * num_components);
    // values.reserve(mesh.num_entities(tdim) * grids_.num_points() * num_components);

    // after streamlining
    values.reserve ( grids_.num_visu_points ( ) * num_components );

    for ( mesh::EntityIterator it = mesh.begin ( tdim ), end_it = mesh.end ( tdim );
          it != end_it; ++it )
    {

        if ( parallel_visualization_ )
        {
            int rem_ind = -100;
            it->get<int>( "_remote_index_", &rem_ind );
            if ( rem_ind != -1 ) continue;
        }

        cell_values.clear ( );
        cell_values.resize ( grids_.num_points ( it->cell_type ( ).tag ( ) ) * num_components, 1.e32 );

        fun ( *it, grids_.coords ( it->cell_type ( ).tag ( ) ), cell_values );
        values.insert ( values.end ( ), cell_values.begin ( ), cell_values.end ( ) );

    }

    functions_.insert ( std::make_pair ( name, values ) );
    num_components_.insert ( std::make_pair ( name, num_components ) );
}

#endif
