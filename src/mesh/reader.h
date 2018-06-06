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

#ifndef HIFLOW_MESH_READER_H
#    define HIFLOW_MESH_READER_H

#    include <exception>
#    include <string>

#    include <mpi.h>
#    include "mesh/mesh_builder.h"

namespace hiflow
{
    namespace mesh
    {

        class Mesh;

        /// \brief Abstract interface for reading in grids from files.

        class Reader
        {
          public:
            /// \brief Ctor for abstract reader.
            explicit Reader ( MeshBuilder* mb );

            /// \brief Dtor for abstract reader.

            virtual ~Reader ( )
            {
            }

            /// \brief Reads a mesh.
            void read ( const char* filename, MeshPtr& mesh );

          protected:
            /// Returns pointer to MeshBuilder object.
            MeshBuilder* mesh_builder ( );

          private:
            /// \brief Function to be implemented by concrete reader
            /// class.
            ///
            /// \details Pure virtual function implemented in concrete
            /// reader subclass which reads the vertices and entities in
            /// the given file, hands them to the MeshBuilder, calls the
            /// MeshBuilder to create a mesh object containing the
            /// entities.
            /// \param [in]  filename The name of the mesh file to be read.
            /// \param [out] mesh     The view of the read-in mesh.
            virtual void read_file ( const char* filename, MeshPtr& mesh ) = 0;

            /// Pointer to the MeshBuilder object that creates the mesh.
            MeshBuilder* mesh_builder_;
        };

        /// \brief Exception class that detects errors while reading a
        /// mesh.

        class ReadGridException : public std::exception
        {
          public:
            /// \brief Ctor for ReadGridException.
            /// \param [in] msg The message to the user if some exception
            /// occurs.
            /// \param [in] line The line in the mesh file where the
            /// exception occured.

            ReadGridException ( const std::string& msg,
                                const std::string& line ) throw ( )
            : _msg ( msg ), _line ( line )
            {
            }

            /// \brief Dtor for ReadGridException

            virtual ~ReadGridException ( ) throw ( )
            {
            }

            /// \brief Function to determine what kind of exception
            /// occured.

            virtual const char* what ( ) const throw ( )
            {
                std::string error ( _msg );
                error.append ( std::string ( "\nparse line: " ) );
                error.append ( _line );
                return error.c_str ( );
            }

          private:
            /// Message to give to the user.
            std::string _msg;
            /// Line in which the exception occured.
            std::string _line;
        };

        /// \brief Concrete grid reader, that reads a grid in the UCD
        /// format.

        class UcdReader : public Reader
        {
          public:
            explicit UcdReader ( MeshBuilder* mb );
          private:
            /// \brief Reads a UCD grid
            void read_file ( const char* filename, MeshPtr& mesh );
        };

        /// \brief Concrete grid reader, that reads a grid in the VTK
        /// Unstructured Grid format.

        class VtkReader : public Reader
        {
          public:
            explicit VtkReader ( MeshBuilder* mb );
          private:
            /// \brief Reads a VTK Unstructured grid.
            void read_file ( const char* filename, MeshPtr& mesh );
        };

        /// \brief Concrete grid reader, that reads a grid in the parallel
        /// VTK Unstructured Grid format.
        ///
        /// \details This special reader reads a description (*.pvtu)
        /// file that has the filenames of the actual data in it. These
        /// actual data files are then passed to the VtkReader according
        /// to the (MPI) process number and the datafile number.

        class PVtkReader : public Reader
        {
          public:
            explicit PVtkReader ( MeshBuilder* mb, const MPI_Comm& mpi_comm );
          private:
            /// The MPI Communicator.
            const MPI_Comm mpi_comm_;

            /// \brief Reads a PVTK Parallel Unstructured grid.
            void read_file ( const char* filename, MeshPtr& mesh );
        };

        /// \brief Concrete grid reader, that reads a grid in the UCD
        /// format by using p4est routines

        class UcdPXestReader : public Reader
        {
          public:
            explicit UcdPXestReader ( MeshBuilder* mb );
          private:
            /// \brief Reads a UCD grid
            void read_file ( const char* filename, MeshPtr& mesh );
            //        MeshP4estBuilder* p4est_mb_;
        };
    }
} // namespace hiflow

#endif /* _GRID_READER_H_ */
