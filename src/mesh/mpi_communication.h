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

/// \author Staffan Ronnas

#ifndef HIFLOW_MESH_MPI_COMMUNICATION_H
#    define HIFLOW_MESH_MPI_COMMUNICATION_H

#    include "communication.h"

#    include <mpi.h>

namespace hiflow
{
    namespace mesh
    {

        /// \brief Broadcasting communicator

        class MpiBroadcast : public Communicator
        {
          public:
            MpiBroadcast ( const MPI_Comm& mpi_communicator, int root );

            // Broadcast entities in sent_entities (non-zero on root proc only),
            // to recv_entities (defined on all procs).
            void communicate ( /*const*/ EntityPackage* sent_entities,
                               EntityPackage* received_entities ) const;

          private:
            bool is_root ( ) const;
            int rank ( ) const;

            const int root_;
            int rank_;
            MPI_Comm mpi_comm_;
        };

        /// \brief Scattering communicator

        class MpiScatter : public Communicator
        {
          public:
            MpiScatter ( const MPI_Comm& mpi_communicator, int sender,
                         const std::vector<EntityCount>& num_entities_on_proc );

            void communicate ( /*const*/ EntityPackage* sent_entities,
                               EntityPackage* received_entities ) const;

            void communicate ( std::vector<EntityPackage>& sent_entities,
                               EntityPackage& received_entities ) const;

          private:
            bool is_sender ( ) const;
            int rank ( ) const;

            const int sender_;
            MPI_Comm mpi_comm_;
            const std::vector<EntityCount>& num_entities_on_proc_;
        };

        /// \brief Non-blocking point-to-point communicator

        class MpiNonBlockingPointToPoint : public NonBlockingCommunicator
        {
          public:
            MpiNonBlockingPointToPoint ( const MPI_Comm& communicator, int sender, int receiver );

            void communicate ( /*const*/ EntityPackage* sent_entities,
                               EntityPackage* received_entities ) const;
            void wait ( ) const;

          private:
            int rank ( ) const;
            bool is_sender ( ) const;
            bool is_receiver ( ) const;

            const int sender_;
            const int receiver_;
            MPI_Comm comm_;

            mutable bool is_communicating_;

            mutable int local_buf_[5];

            mutable MPI_Request local_buf_request_;
            mutable MPI_Request coords_request_;
            mutable MPI_Request offsets_request_;
            mutable MPI_Request connections_request_;
            mutable MPI_Request materials_request_;
            mutable MPI_Status status_;
        };

    }
} // namespace hiflow

#endif /* _MPI_COMMUNICATION_H_ */
