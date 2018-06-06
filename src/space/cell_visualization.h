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

#ifndef HIFLOW_SPACE_CELL_VISUALIZATION
#    define HIFLOW_SPACE_CELL_VISUALIZATION

/// \author Staffan Ronnas, Martin Baumann, Teresa Beck, Simon Gawlok, Jonas Kratzke
///
/// \brief Visualization of finite element functions.
///
/// Using this class a Vtk (http://www.vtk.org/) unstructured grid visualization
/// file can be created. Please find detailed information about Vtk's file
/// formats at http://www.vtk.org/VTK/img/file-formats.pdf.
/// This type of visualization writes out every cell and with function values
/// provided by a user-defined evaluation function.
///
/// Please note for simulations with multiple visualization calls, that this class
/// is NOT ment to be initialized once for several visualization calls. Please
/// construct a new instantiation of the CellVisualization every single time you want
/// to visualize your data.
///

#    include "mesh/mesh.h"

#    include <map>
#    include <string>
#    include <vector>

#    include <boost/function.hpp>
#    include <mesh/types.h>
#    include <mpi.h>

#    include "space/visualization_tools.h"
#    include "common/grid.h"
#    include "common/bbox.h"

namespace hiflow
{

    template<class DataType>
    class VectorSpace;

    /// \brief Description of a square 2d grid or cubic 3d grid.

    template<class DataType>
    class CellVisualizationGrids
    {
      public:
        CellVisualizationGrids ( const Mesh* mesh, int num_intervals, DataType origin, DataType side_length );
        int num_visu_points ( ) const;
        int num_visu_cells ( ) const;
        int num_points ( CellType::Tag cell_tag ) const;
        int num_cells ( CellType::Tag cell_tag ) const;
        const std::vector<int>& vertices_of_cell ( CellType::Tag cell_tag, int i ) const;
        const std::vector<DataType>& coords ( CellType::Tag cell_tag ) const;

      private:
        typename ScopedPtr< Grid<DataType> >::Type grids_[CellType::NUM_CELL_TYPES];
        const int tdim_;
        int num_visu_points_;
        int num_visu_cells_;
    };

    /// \brief Visualization of finite element solutions.

    template<class DataType>
    class CellVisualization
    {
      public:
        // Type of function for evaluation.
        typedef boost::function3<void,
        const mesh::Entity&, // cell
        const std::vector<DataType>&, // reference coordinates
        std::vector<DataType>& // values of function at the points
        > EvalFunction;

        explicit CellVisualization ( const VectorSpace<DataType>& space, int num_intervals );

        void visualize ( const EvalFunction& fun, const std::string& name );
        void visualize_cell_data ( const std::vector<DataType>& cell_data, const std::string& name );

        void write ( const std::string& filename ) const;
      protected:
        bool parallel_visualization_;

        const VectorSpace<DataType>& space_;
        CellVisualizationGrids<DataType> grids_;
        std::map< std::string, std::vector<DataType> > functions_;
        std::map< std::string, std::vector<DataType> > functions_cell_;
    };

    /// \brief Writer for Pvtk files.
    /// \details Write PVtk files and also the corresponding Vtk files.

    template<class DataType>
    class ParallelCellVisualization : public CellVisualization<DataType>
    {
      public:
        /// \brief Ctor for PVtkWriter.
        /// \param [in] mpi_comm MPI Communicator.

        explicit ParallelCellVisualization ( const VectorSpace<DataType>& space,
                                             int num_intervals,
                                             const MPI_Comm& mpi_comm,
                                             const int master_rank )
        : CellVisualization<DataType>( space, num_intervals ), mpi_comm_ ( mpi_comm ), master_rank_ ( master_rank )
        {
        }

        /// \brief Writes a parallel vtk unstructured grid.
        void write ( const std::string& filename, const std::string& path = "" ) const;

      protected:

        /// The MPI Communicator.
        MPI_Comm mpi_comm_;
        const int master_rank_;
    };
}

#endif
