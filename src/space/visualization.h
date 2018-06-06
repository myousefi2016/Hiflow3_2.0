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

#ifndef __VISUALIZATION_H__
#    define __VISUALIZATION_H__

#    include <vector>

#    include "dof/degree_of_freedom.h"
#    include "dof/dof_partition.h"
#    include "mesh/entity.h"
#    include "fem/femanager.h"
#    include "mesh/iterator.h"
#    include "mesh/mesh.h"
#    include "mesh/types.h"
#    include "common/log.h"

namespace hiflow
{
    /// \deprecated This class is deprecated. Use CellVisualization
    /// if possible.
    ///
    /// @brief Visualization of finite element functions, f.e. solutions vectors.
    /// @author Martin Baumann
    ///
    /// Using this class a Vtk (http://www.vtk.org/) unstructured grid visualization
    /// file can be created. Please find detailed information about Vtk's file
    /// formats at http://www.vtk.org/VTK/img/file-formats.pdf.
    /// A finite element solution is given by a vector that represent the
    /// values for the degrees of freedom that describe the discrete solution
    /// function. Further the mesh, the description of the finite elements
    /// and the description of the degrees of freedom are needed for
    /// the visualization.
    ///
    /// There are two ways of visualization:
    ///
    /// -# visualization of a Solution instance \n
    ///    a Solution instance has information about its corresponding VectorSpace
    ///    -> mesh, dof and fe_manager are already known
    /// -# visualization of a std::vector<T> \n
    ///    a mesh, dof and fe_manager must be set before visualization
    /// \Example
    /// \code
    ///  Visualization vis;
    ///  vis.set_mesh       (mesh);
    ///  vis.set_dof        (dof);
    ///  vis.set_fe_manager (fe_manager);
    ///  vis.set_filename   ("solution.vtu");
    ///  vis.Visualize      (solution_vector);
    /// \endcode
    ///
    /// \todo
    /// - Visualization of Solution instance
    /// - use better names than 'val' (information in Solution)

    template<class DataType>
    class Visualization
    {
      public:

        typedef mesh::Mesh Mesh;
        typedef mesh::Entity MeshEntity;
        typedef mesh::EntityIterator MeshEntityIterator;
        typedef mesh::EntityNumber MeshEntityNumber;
        typedef mesh::IncidentEntityIterator MeshIncidentEntityIterator;
        typedef mesh::VertexIdIterator MeshVertexIdIterator;
        typedef mesh::Coordinate Coordinate;

        //typedef doffem::DegreeOfFreedom  DegreeOfFreedom;
        typedef doffem::DofPartition<DataType> DegreeOfFreedom;
        typedef doffem::FEManager<DataType> FEManager;

      protected:

        Mesh* mesh_;
        const DegreeOfFreedom* dof_;
        const FEManager* fe_manager_;
        std::string filename_;
        std::vector<std::string> var_names_;
        bool visualize_all_attributes_;

      public:

        Visualization ( );
        virtual ~Visualization ( );

        void set_filename ( std::string name )
        {
            filename_ = name;
            LOG_INFO ( "Filename", filename_ );
        }

        //void Visualize (Solution const& solution);

        virtual void Visualize ( const std::vector<double>& vec );

        void set_mesh ( Mesh* mesh )
        {
            mesh_ = mesh;
        }

        void set_dof ( const DegreeOfFreedom* dof )
        {
            dof_ = dof;
        }

        void set_fe_manager ( const FEManager* fem )
        {
            fe_manager_ = fem;
        }
        void set_var_names ( const std::vector<std::string>& names );
        void clear_var_names ( );
        void set_visualize_all_attributes ( bool flag );

    };

} // namespace hiflow

#endif
