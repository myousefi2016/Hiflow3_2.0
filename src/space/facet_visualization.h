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

#ifndef HIFLOW_SPACE_FACET_VISUALIZATION
#    define HIFLOW_SPACE_FACET_VISUALIZATION

/// Visualization of facet within each cell, w.r.t. required material number

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

    /// \brief Visualization of finite element solutions.

    template<class DataType>
    class FacetVisualization
    {
      public:
        // Type of function for evaluation.
        typedef boost::function3<void,
        const mesh::Entity&, // cell
        const std::vector<DataType>&, // reference coordinates
        std::vector<DataType>& // values of function at the points
        > EvalFunction;

        explicit FacetVisualization ( const VectorSpace<DataType>& space, std::vector<int>& material_numbers );

        void visualize ( const EvalFunction& fun, const std::string& name );
        void visualize_cell_data ( const std::vector<DataType>& cell_data, const std::string& name );

        void write ( const std::string& filename ) const;
      protected:
        bool parallel_visualization_;

        const VectorSpace<DataType>& space_;
        std::map< std::string, std::vector<DataType> > functions_;
        std::map< std::string, std::vector<DataType> > functions_cell_;
        int num_visu_points_;
        int num_visu_cells_;
        std::vector<int> material_numbers_;
    };

    /// \brief Writer for Pvtk files.
    /// \details Write PVtk files and also the corresponding Vtk files.

    template<class DataType>
    class ParallelFacetVisualization : public FacetVisualization<DataType>
    {
      public:
        /// \brief Ctor for PVtkWriter.
        /// \param [in] mpi_comm MPI Communicator.

        explicit ParallelFacetVisualization ( const VectorSpace<DataType>& space,
                                              std::vector<int>& material_numbers,
                                              const MPI_Comm& mpi_comm,
                                              const int master_rank )
        : FacetVisualization<DataType>( space, material_numbers ), mpi_comm_ ( mpi_comm ), master_rank_ ( master_rank )
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
