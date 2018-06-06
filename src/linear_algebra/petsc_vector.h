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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-11-17

#ifndef HIFLOW_LINEARALGEBRA_PETSC_VECTOR_H_
#    define HIFLOW_LINEARALGEBRA_PETSC_VECTOR_H_

#    include "common/log.h"
#    include "config.h"
#    include "linear_algebra/vector.h"
#    include "linear_algebra/coupled_vector.h"
#    include "linear_algebra/la_couplings.h"
#    include "tools/mpi_tools.h"
#    include "mpi.h"

// TODO: Fix order dependency of next inclusion
#    include "common/smart_pointers.h"

namespace hiflow
{
    namespace la
    {

        /// Forwarding PETSc vector wrapper class
        namespace petsc
        {
            class Vec_wrapper;
        }

        /// Forwarding PETSc matrix class
        template <class DataType>
        class PETScMatrix;

        /// Forwarding PETSc linear solver
        template <class LAD>
        class PETScCG;

        /// @brief Wrapper class to PETSC vector

        template <class DataType>
        class PETScVector : public Vector<DataType>
        {
          public:
            /// Default constructor
            PETScVector ( );

            /// Destructor
            virtual ~PETScVector ( );

            virtual Vector<DataType> *Clone ( ) const;

            /// Initialize matrix
            /// @param[in] comm MPI communicator to be used by vector
            void Init ( MPI_Comm &comm );

            /// Initialize matrix
            /// @param[in] comm MPI communicator to be used by Vector
            /// @param[in] cp Linear algebra couplings describing ownerships
            void Init ( MPI_Comm &comm, const LaCouplings &cp );

            /// Clear all allocated data
            void Clear ( );

            /// Clone vector structure without content
            void CloneFromWithoutContent ( const PETScVector<DataType> &vec );

            /// Copy vector content to this Vector
            void CopyFrom ( const PETScVector<DataType> &vec );

            void CloneFrom ( const PETScVector<DataType> &vec )
            {
                this->CloneFromWithoutContent ( vec );
                this->CopyFrom ( vec );
            }

            /// Get laCouplings

            const LaCouplings &la_couplings ( ) const
            {
                return *cp_;
            }
            /// Get MPI communicator

            const MPI_Comm &comm ( ) const
            {
                return comm_;
            }
            /// Return mpi rank

            int my_rank ( ) const
            {
                return my_rank_;
            }

            /// Local size of Vector
            int size_local ( ) const;
            /// Global size of Vector
            int size_global ( ) const;

            /// Update operator, i.e. exchange values for distributed vectors
            void Update ( );
            /// Initiate update

            void begin_update ( )
            {
            }
            /// Finalize update

            void end_update ( )
            {
            }

            /// Send border values
            void SendBorder ( );
            /// Receive ghost values
            void ReceiveGhost ( );
            /// Wait for finishing of send operations
            void WaitForSend ( );
            /// Wait for finishing of receive operations
            void WaitForRecv ( );
            /// Update ghost values
            void UpdateGhost ( );
            /// Initialize post-processing data structures
            void InitializePostProcessing ( const Couplings<DataType> *couplings );
            /// Update post-processing values
            void UpdatePpValues ( );
            /// Set ghost values
            void SetGhostValues ( const DataType *values );
            /// Get all dofs and values (interior, ghost and pp data)
            void GetAllDofsAndValues ( std::vector<int> &id, std::vector<DataType> &val );

            /// @return Local size of ghost

            int size_local_ghost ( ) const
            {
                assert ( this->ghost_ != NULL );
                return this->ghost_->get_size ( );
            }

            /// Set all elements to zero
            void Zeros ( );

            /// Get value at a known index
            DataType GetValue ( int index ) const;
            /// Get values at given indices
            void GetValues ( const int *indices, int length, DataType *values ) const;

            /// Euclidean length of vector
            DataType Norm2 ( ) const;
            /// l1 norm of vector
            DataType Norm1 ( ) const;
            /// maximum absolute value of entries
            DataType NormMax ( ) const;

            /// Scalar product
            DataType Dot ( const Vector<DataType> &vec ) const;
            /// Scalar product
            DataType Dot ( const PETScVector<DataType> &vec ) const;

            /// Add value to a given index
            void Add ( int index, DataType scalar );
            /// Add values to given indices
            void Add ( const int *indices, int length, const DataType *values );
            /// Set given global index to value

            void SetValue ( int index, DataType value );
            /// Set given global indices to given values
            void SetValues ( const int *indices, int length,
                             const DataType *values );

            /// this <- this + alpha * vecx
            void Axpy ( const Vector<DataType> &vecx, DataType alpha );
            void Axpy ( const PETScVector<DataType> &vecx, DataType alpha );

            /// this <- alpha * this + vecx
            void ScaleAdd ( const Vector<DataType> &vecx, DataType alpha );
            void ScaleAdd ( const PETScVector<DataType> &vecx, DataType alpha );

            /// this <- alpha * this
            void Scale ( DataType alpha );

            const lVector<DataType> &ghost ( ) const
            {
                return *( this->ghost_ );
            }

            lVector<DataType> &ghost ( )
            {
                return *( this->ghost_ );
            }

            const PpData<DataType> &pp_data ( ) const
            {
                return *( this->pp_data_ );
            }

            bool HasPpData ( ) const
            {
                return (this->pp_data_ != NULL );
            }

          private:

            /// Friends
            friend class PETScMatrix<DataType>;
            template <class LAD> friend class PETScCG;
            template <class LAD> friend class PETScGeneralKSP;

            /// MPI communicator
            MPI_Comm comm_;
            /// Global number of first component owned by this process
            int ilower_;
            /// Global number of last component owned by this process
            int iupper_;

            /// Linear algebra couplings describing global dof distribution
            const LaCouplings *cp_;

            lVector<DataType> *ghost_;

            /// Number of processes.
            int nb_procs_;
            /// Rank of this process.
            int my_rank_;
            /// Number of send operations for border values
            int nb_sends_;
            /// Number of receive operations for ghost values
            int nb_recvs_;
            /// MPI requests
            std::vector<MPI_Request> mpi_req_;
            /// MPI status
            std::vector<MPI_Status> mpi_stat_;
            /// border values to be send
            std::vector<DataType> border_val_;
            /// ghost values to be received
            std::vector<DataType> ghost_val_;
            /// Data for post processing with the rule: pp_data_ != NULL iff
            /// pp_data_ initialized.
            PpData<DataType> *pp_data_;
            std::vector<int> border_indices_;

            /// Is true if we already checked for a dof partition in la_couplings
            bool checked_for_dof_partition_;

            /// Flag if vector is already initialized
            bool initialized_;

            /// Pointer to PETSc matrix object
            hiflow::scoped_ptr<petsc::Vec_wrapper> ptr_vec_wrapper_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_PETSC_VECTOR_H_
