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

/// @author Chandramowli Subramanian, Nico Trost, Dimitar Lukarski, Martin Wlotzka

#ifndef HIFLOW_LINEARALGEBRA_COUPLED_VECTOR_H_
#    define HIFLOW_LINEARALGEBRA_COUPLED_VECTOR_H_

#    include <iostream>
#    include "linear_algebra/vector.h"
#    include "linear_algebra/la_couplings.h"
#    include "linear_algebra/couplings.h"
#    include "linear_algebra/coupled_vector_factory.h"
#    include "linear_algebra/lmp/la_global.h"
#    include "linear_algebra/lmp/cuda/lvector_gpu.h"
#    include "common/property_tree.h"
#    include "linear_algebra/lmp/platform_management.h"

#    include "dof/dof_partition.h" //for PostProcessing

#    include "mpi.h"

namespace hiflow
{
    namespace la
    {

        template <class DataType> class lVector;
        template <class DataType> class PpData;

        /// @brief Distributed vector.
        ///
        /// This vector stores also couplings. After matrix-vector multiplication
        /// these couplings have to be updated.

        template<class DataType>
        class CoupledVector : public Vector<DataType>
        {
          public:
            /// Standard constructor
            CoupledVector ( );

            /// Destructor
            virtual ~CoupledVector ( );

            /// Inits empty vector.
            /// @param comm MPI communicator
            /// @param cp LaCouplings (see la_couplings.h)
            /// @param plat System platform (see lmp/la_global.h)
            /// @param impl Implementation (see lmp/la_global.h)
            void Init ( const MPI_Comm& comm, const LaCouplings& cp,
                        PLATFORM plat, IMPLEMENTATION impl );

            void Init ( const MPI_Comm& comm, const LaCouplings& cp,
                        PLATFORM plat, IMPLEMENTATION impl,
                        const SYSTEM& my_system );

            /// 2 new funcitons for Init
            void Init ( const MPI_Comm& comm, const LaCouplings& cp );
            void Init_la_system ( PLATFORM plat, IMPLEMENTATION impl );
            /// Initializes the structure of the vector.
            /// LaCouplings need to be initialized.
            void InitStructure ( );

            /// Initialize the Data needed for post processing (e.g. visualization)
            void InitializePostProcessing ( const Couplings<DataType>* couplings );

            /// Sets every element to zero.
            void Zeros ( );

            // Functions inherited from Vector class
            DataType Dot ( const Vector<DataType>& vec ) const;
            void Axpy ( const Vector<DataType>& vec, const DataType alpha );
            void ScaleAdd ( const Vector<DataType>& vec, const DataType alpha );
            void Add ( const int* indices, const int size_indices, const DataType* values );

            /// Update operator, i.e. exchange values for distributed vectors

            void Update ( )
            {
                this->begin_update ( );
                this->end_update ( );
            }
            void begin_update ( );
            void end_update ( );

            // Specializations of inherited member functions of Vector to CoupledVector
            DataType Dot ( const CoupledVector<DataType>& vec ) const;
            void Axpy ( const CoupledVector<DataType>& vec, const DataType alpha );
            void ScaleAdd ( const CoupledVector<DataType>& vec, const DataType alpha );

            /// Adds a value. Component must be owned by process.
            /// @param global_dof_id Global dof index
            /// @param val Value to add
            void Add ( int global_dof_id, DataType val );

            /// Scales this vector.
            /// @param alpha Scale factor
            void Scale ( const DataType alpha );

            /// Returns the one-norm of this vector.
            /// @return One-norm
            DataType Norm1 ( ) const;

            /// Returns the two-norm of this vector.
            /// @return Two-norm
            DataType Norm2 ( ) const;

            /// Returns the maximum-norm of this vector.
            /// @return Maximum-norm
            DataType NormMax ( ) const;

            /// Retrieve certain entries of the vector.
            /// The user has to provide the allocated array @em values.
            /// @param indices Global indices of the entries to retrieve
            ///        These indices have to be owned by the process invoking this function
            ///        or have to be included as couplings, i.e. in the ghost vector.
            ///        Otherwise, zero is returned (so that assembling with linearization
            ///        on a cell where no dof is owned works)
            /// @param size_indices Size of the array @em indices
            /// @param values Array of retrieved values
            virtual void GetValues ( const int* indices, const int size_indices, DataType* values ) const;

            /// Retrieves the whole vector (interior).
            /// @param values Array of values, must fit size of interior
            void GetValues ( DataType* values ) const;

            /// Get value at a known index
            ///@param index Global index of the entry to retrieve. Index has to be owned by
            ///             the process invoking this function or has to be included as couplings.
            ///             Otherwise, zero is returned.
            virtual DataType GetValue ( const int index ) const;

            /// @return All Dofs and values that are in interior, ghost and pp_data_.
            /// They are NOT sorted.
            void GetAllDofsAndValues ( std::vector<int>& id, std::vector<DataType>& val ) const;

            /// Sets value at global index.
            /// @param index Global index.
            /// @param value Value to be set.
            virtual void SetValue ( const int index, const DataType value );

            /// Sets values at given global indices.
            /// @param indices Global indices of the entries to be set.
            ///        These indices have to be owned by the process invoking this function.
            /// @param size_indices Size of the array @em indices
            /// @param values Array of values to be set
            virtual void SetValues ( const int* indices, const int size_indices, const DataType* values );

            /// Sets all the values in the vector (interior).
            /// @param values Array of values, must fit size of interior
            void SetValues ( const DataType* values );

            /// Sets all ghost values. Used inside UpdateCouplings.
            /// @param values Values to set of size @em size_ghost
            void SetGhostValues ( const DataType* values );

            /// Gathers the vector (only CPU vectors).
            /// @param recv_id Process id which gathers
            /// @param values Array of values to be stored (needs to be allocated: only on @em recv_id)
            void Gather ( int recv_id, DataType* values ) const;

            virtual Vector<DataType> * Clone ( ) const;
            /// Clones the whole vector (everything).
            /// @param vec Vector to be cloned
            void CloneFrom ( const CoupledVector<DataType>& vec );

            /// Copies interior and ghost (no structure, no platform).
            /// Vector already needs to be initialized.
            /// @param vec Vector to copy
            void CopyFrom ( const CoupledVector<DataType>& vec );

            /// Copies only the interior (no ghost, no structure, no platform).
            void CopyInteriorFrom ( const CoupledVector<DataType>& vec );

            /// Cast data of interior from another CoupledVector.
            void CastInteriorFrom ( const CoupledVector<double>& vec );
            void CastInteriorFrom ( const CoupledVector<float>& vec );

            /// Cast data to interior of another CoupledVector.
            void CastInteriorTo ( CoupledVector<double>& vec ) const;
            void CastInteriorTo ( CoupledVector<float>& vec ) const;

            void CopyTo ( CoupledVector<DataType>& vec ) const;

            /// Copies only the structure.
            void CopyStructureFrom ( const CoupledVector<DataType>& vec );

            /// Copies everything but the values.
            /// @param vec Vector whose strucutre is to be copied
            void CloneFromWithoutContent ( const CoupledVector<DataType>& vec );

            /// Clears allocated local vectors, i.e. interior and ghost.
            void Clear ( );

            /// Sends border values.
            void SendBorder ( );

            /// Receives ghost values.
            void ReceiveGhost ( );

            /// Waits for the completion of the send operations of the border values.
            void WaitForSend ( );

            /// Waits for the completion of the receive operations of the ghost values
            /// and stores them.
            void WaitForRecv ( );

            /// Updates couplings, i.e. exchange border and ghost values and
            /// values needed for post processing (if possible) between processors.
            void UpdateCouplings ( );

            /// Updates Ghost, i.e. exchange border and ghost values between processors.
            void UpdateGhost ( );

            /// Communicate the values needed for post processing and not owned
            /// by this process
            void UpdatePpValues ( );

            /// @return Global size of the vector
            virtual int size_global ( ) const;

            /// Prints information to stream \em out.
            /// @param out Stream for output
            void Print ( std::ostream &out = std::cout ) const;

            void WriteHDF5 ( const std::string& filename, const std::string& groupname,
                             const std::string& datasetname );

            void ReadHDF5 ( const std::string& filename, const std::string& groupname,
                            const std::string& datasetname );

            // inline functions
            /// @return Local size of the vector, i.e. size of the interior.

            virtual int size_local ( ) const
            {
                assert ( this->interior_ != NULL );
                return this->interior_->get_size ( );
            }

            /// @return Local size of ghost

            int size_local_ghost ( ) const
            {
                assert ( this->ghost_ != NULL );
                return this->ghost_->get_size ( );
            }

            /// @param begin Global index of the first local entry
            /// @param end One more than the global index of the last local entry

            void GetOwnershipRange ( int* begin, int* end ) const
            {
                assert ( this->ownership_begin_ > -1 );
                assert ( this->ownership_end_ > -1 );

                *begin = this->ownership_begin_;
                *end = this->ownership_end_;
            }

            /// @return Global index of the first local entry

            int ownership_begin ( ) const
            {
                return this->ownership_begin_;
            }

            /// @return One more than the global index of the last local entry

            int ownership_end ( ) const
            {
                return this->ownership_end_;
            }

            /// @return LaCouplings

            const LaCouplings& la_couplings ( ) const
            {
                return *( this->la_couplings_ );
            }

            /// @return Local vector holding the interior entries

            const lVector<DataType>& interior ( ) const
            {
                return *( this->interior_ );
            }

            lVector<DataType>& interior ( )
            {
                return *( this->interior_ );
            }

            /// @return Local vector holding the couplings

            const lVector<DataType>& ghost ( ) const
            {
                return *( this->ghost_ );
            }

            lVector<DataType>& ghost ( )
            {
                return *( this->ghost_ );
            }

            /// @return All the pp_data (including communication data)

            const PpData<DataType>& pp_data ( ) const
            {
                return *( this->pp_data_ );
            }

            /// @return MPI communicator

            const MPI_Comm& comm ( ) const
            {
                return this->comm_;
            }

            /// @return Number of processes

            int nb_procs ( ) const
            {
                return this->nb_procs_;
            }

            /// @return Rank of this process

            int my_rank ( ) const
            {
                return this->my_rank_;
            }

            /// @return Number of send operations for border values

            int nb_sends ( ) const
            {
                return this->nb_sends_;
            }

            /// @return Number of receive operations for ghost values

            int nb_recvs ( ) const
            {
                return this->nb_recvs_;
            }

            /// @return std::vector of MPI requests for communication

            std::vector<MPI_Request> mpi_req ( ) const
            {
                return this->mpi_req_;
            }

            /// @return std::vector of MPI status for communication

            std::vector<MPI_Status> mpi_stat ( ) const
            {
                return this->mpi_stat_;
            }

            /// @return std::vector of border values to be sent

            std::vector<DataType> border_val ( ) const
            {
                return this->border_val_;
            }

            /// @return std::vector for ghost values to be received

            std::vector<DataType> ghost_val ( ) const
            {
                return this->ghost_val_;
            }

            /// @return true if pp_data_ has been initialized
            bool HasPpData ( ) const;

          private:
            // no implementation of copy constructor or assignement operator
            CoupledVector ( const CoupledVector<DataType>& );
            CoupledVector<DataType>& operator= ( const CoupledVector<DataType>& );

            /// Computes ownership range via LaCouplings.
            void ComputeOwnershipRange ( );
            /// Global index of the first local entry.
            int ownership_begin_;
            /// One more than the global index of the last local entry.
            int ownership_end_;

            const LaCouplings* la_couplings_;
            lVector<DataType>* interior_; // entries of own domain
            lVector<DataType>* ghost_;

            /// MPI communicator.
            MPI_Comm comm_;
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
            PpData<DataType>* pp_data_;

            /// Is true if we already checked for a dof partition in la_couplings
            bool checked_for_dof_partition_;
        };

        /// This struct has all the data needed to do post processing within
        /// a CoupledVector

        template<class DataType>
        class PpData
        {
          public:

            PpData ( )
            {
            };
            /// Values and DoFs not owned by this process but which belongs to
            /// this' process mesh.
            std::vector<DataType> values;
            std::vector<int> dof_ids;
            SortedArray<int> sorted_dof_ids;
            std::vector<int> dof_ids_offsets;
            /// Number of send operations for postprocessing values
            std::vector<int> nb_sends;
            /// Number of receive operations for postprocessing values
            std::vector<int> nb_recvs;
            /// DoFs other processes need
            std::vector<int> send_dofs;
            std::vector<int> send_offsets;

            /// Creates a clone without values and returns a pointer pointing at it

            PpData* Clone ( ) const
            {
                PpData* pp_dat = new PpData;
                pp_dat->dof_ids = this->dof_ids;
                pp_dat->sorted_dof_ids = SortedArray<int>( this->dof_ids.begin ( ), this->dof_ids.end ( ) );
                pp_dat->dof_ids_offsets = this->dof_ids_offsets;
                pp_dat->nb_sends = this->nb_sends;
                pp_dat->nb_recvs = this->nb_recvs;
                pp_dat->send_dofs = this->send_dofs;
                pp_dat->send_offsets = this->send_offsets;
                pp_dat->values = this->values;
                return pp_dat;
            }
            /// Simply copies everything contained by the struct

            PpData ( const PpData& pp_dat )
            {
                this->dof_ids = pp_dat.dof_ids;
                this->sorted_dof_ids.resize ( pp_dat.sorted_dof_ids.size ( ) );
                for ( int i = 0; i < this->sorted_dof_ids.size ( ); ++i )
                {
                    this->sorted_dof_ids[i] = pp_dat.sorted_dof_ids[i];
                }
                this->dof_ids_offsets = pp_dat.dof_ids_offsets;
                this->nb_sends = pp_dat.nb_sends;
                this->nb_recvs = pp_dat.nb_recvs;
                this->send_dofs = pp_dat.send_dofs;
                this->send_offsets = pp_dat.send_offsets;
                this->values = pp_dat.values;
            }
        };

        template<class DataType>
        class CoupledVectorCreator
        {
          public:

            CoupledVector<DataType>* params ( const PropertyTree& c )
            {
                CoupledVector<DataType>* newCoupledVector = new CoupledVector<DataType>( );
                if ( c.contains ( "Platform" ) && c.contains ( "Implementation" ) &&
                     c.contains ( "MatrixFormat" ) )
                {
                    const std::string platform_str = c["Platform"].template get<std::string>( );
                    if ( platform_str == "CPU" )
                    {
                        la_sys_.Platform = CPU;
                    }
                    else if ( platform_str == "GPU" )
                    {
                        la_sys_.Platform = GPU;
                    }
                    else
                    {
                        LOG_ERROR ( "CoupledVectorCreator::params: No format of this name registered(platform)." );
                        return NULL;
                    }
                    init_platform ( la_sys_ );
                    const std::string impl_str = c["Implementation"].template get<std::string>( );
                    if ( impl_str == "Naive" )
                    {
                        la_impl_ = NAIVE;
                    }
                    else if ( impl_str == "BLAS" )
                    {
                        la_impl_ = BLAS;
                    }
                    else if ( impl_str == "MKL" )
                    {
                        la_impl_ = MKL;
                    }
                    else if ( impl_str == "OPENMP" )
                    {
                        la_impl_ = OPENMP;
                    }
                    else if ( impl_str == "SCALAR" )
                    {
                        la_impl_ = SCALAR;
                    }
                    else if ( impl_str == "SCALAR_TEX" )
                    {
                        la_impl_ = SCALAR_TEX;
                    }
                    else
                    {
                        LOG_ERROR ( "CoupledVectorCreator::params: No format of this name registered(implementation)." );
                        return NULL;
                    }
                    const std::string matrix_str = c["MatrixFormat"].template get<std::string>( );
                    if ( matrix_str == "CSR" )
                    {
                        la_matrix_format_ = CSR;
                    }
                    else if ( matrix_str == "COO" )
                    {
                        la_matrix_format_ = COO;
                    }
                    else
                    {
                        LOG_ERROR ( "CoupledVectorCreator::params: No format of this name registered(matrix format)." );
                        return NULL;
                    }
                    newCoupledVector->Init_la_system ( la_sys_.Platform, la_impl_ );
                    return newCoupledVector;
                }
                return 0;
            }
          private:
            SYSTEM la_sys_;
            IMPLEMENTATION la_impl_;
            MATRIX_FORMAT la_matrix_format_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_COUPLED_VECTOR_H_
