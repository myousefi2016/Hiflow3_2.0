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

/// \author Thomas Gengenbach and Staffan Ronnas

#ifndef HIFLOW_MESH_WRITER_H
#    define HIFLOW_MESH_WRITER_H

#    include "mesh/mesh.h"
#    include "mesh/types.h"

#    include <mpi.h>

namespace hiflow
{
    namespace mesh
    {

        /// \brief Abstract interface for writing grids

        class Writer
        {
          public:
            typedef std::vector<std::pair<std::string, TDim> > DataVectorPair;
            /// \brief Ctor for writer class.

            Writer ( )
            {
            }
            /// \brief Virtual Dtor for writer class.

            virtual ~Writer ( )
            {
            }

            /// \brief Writes a grid to a file with the given name.
            void write ( const char* filename, const Mesh& mesh );

            /// \brief Adds data that is in the mesh and that the user
            /// want to write.
            void add_data_array ( std::string array_name, TDim tdim );

            /// \brief Adds attribute that is in the mesh and the user
            /// wants to write.
            void add_attribute ( const std::string& attribute_name, TDim tdim );

            /// \brief Add all attributes that are defined in the mesh and
            /// write them out.
            void add_all_attributes ( const Mesh& mesh, bool internal_mesh_attributes = false );

          protected:
            /// \brief Holds the name and the topological dimension of
            /// data arrays, that the user wants to write.
            DataVectorPair data_to_write_;
            DataVectorPair attributes_to_write_;
          private:
            /// \brief Virtual function to write the mesh.
            /// \details Pure virtual writer function implemented in the
            /// concrete writer subclass that writes the grid to the file
            /// with the given filename.
            /// \param [in] filename The name of the mesh file to be written.
            /// \param [in] mesh The view of the mesh to be written.
            virtual void write_file ( const char* filename, const Mesh& mesh ) const = 0;
        };

        /// \brief Concrete UCD writer class. Not yet implemented.

        class UcdWriter : public Writer
        {
          private:
            /// \brief Writes a ucd grid. Is not yet implemented.
            void write_file ( const char* filename, const Mesh& mesh ) const;
        };

        /// \brief Concrete VTK writer class.

        class VtkWriter : public Writer
        {
          public:
            VtkWriter ( );
            void set_deformation ( const std::vector<double>* deformation );
          private:
            /// \brief Writes a vtk unstructured grid.
            void write_file ( const char* filename, const Mesh& mesh ) const;

            const std::vector<double>* deformation_;
        };

        /// \brief Concrete PVTK writer class.
        /// \details This special parallel VTK writer calls VtkWriter to
        /// write actual data.

        class PVtkWriter : public Writer
        {
          public:
            /// \brief Ctor for PVtkWriter.
            /// \param [in] mpi_comm MPI Communicator.

            explicit PVtkWriter ( const MPI_Comm& mpi_comm )
            : mpi_comm_ ( mpi_comm ), deformation_ ( 0 )
            {
            }
            void set_deformation ( const std::vector<double>* deformation );
          private:
            /// The MPI Communicator.
            MPI_Comm mpi_comm_;

            /// \brief Writes a parallel vtk unstructured grid.
            void write_file ( const char* filename, const Mesh& mesh ) const;

            const std::vector<double>* deformation_;
        };
    }
} // namespace hiflow
#endif
