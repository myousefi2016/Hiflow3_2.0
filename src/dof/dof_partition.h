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

#ifndef __DOF_PARTITION_H_
#    define __DOF_PARTITION_H_

#    include <mpi.h>
#    include <vector>
#    include <map>
#    include "dof/degree_of_freedom.h"
#    include "fem/femanager.h"
#    include "dof/dof_fem_types.h"
#    include "common/sorted_array.h"

namespace hiflow
{
    namespace doffem
    {

        ///
        /// \class DofPartition dof_partition.h
        /// \brief This class manages the dof partitioning for parallel computation
        /// \author Michael Schick<br>Martin Baumann

        template<class DataType>
        class DofPartition : public DegreeOfFreedom<DataType>
        {
          public:

            /// Default Constructor
            DofPartition ( );
            /// Destructor
            ~DofPartition ( );

            /// Set MPI Communicator if needed, default value is MPI_COMM_WORLD
            void set_mpi_comm ( const MPI_Comm& comm );

            /// Get dof offset of this subdomain
            inline int get_my_dof_offset ( ) const;
            /// Get subdomain index (subdomain which owns dof)
            inline int get_my_subdomain ( ) const;
            /// Get MPI communicator
            inline const MPI_Comm& get_mpi_comm ( ) const;

            /// Returns the number of dofs on a specific subdomain owned by the 
            /// subdomain
            inline int ndofs_on_sd ( int subdomain ) const;
            /// Returns the number of dofs on the whole global computational domain
            /// (including the other subdomains)
            inline int ndofs_global ( ) const;

            /// Returns dofs on cell for a specific variable. The DofIDs are local 
            /// w.r.t. subdomain (Including the dofs which are not owned)
            void get_dofs_on_cell_sd ( int var, int cell_index,
                                       std::vector<DofID>& ids ) const;

            /// Create local numbering for dofs on the subdomain
            /// @param[in] order Ordering strategy for DoFs.
            void number ( DOF_ORDERING order = HIFLOW_CLASSIC );

            /// Permute the local DofIDs w.r.t. the subdomain by a given permutation
            /// (Including the dofs which are not owned)
            void apply_permutation_on_sd ( const std::vector<DofID>& permutation );

            void GetPostProcessingStructure ( std::vector<DofID>& dof, std::vector<int>& dof_offset ) const;

            /// Check if global dof is on subdomain (only if dof is owned)
            inline bool is_dof_on_sd ( DofID global_id ) const;
            /// Returns the lowest subdomain index, which shares the global dof ID
            inline int owner_of_dof ( DofID global_id ) const;

            /// For a given local DofId on a subdomain, this routine computes the global
            /// DofId (including the local dofs which are not owned have local number on
            /// subdomain)
            inline void local2global ( DofID local_id, DofID* global_id ) const;
            /// For a given global DofId, this routine computes the local DofId on the 
            /// current subdomain (including the dofs which are not owned have local 
            /// number on subdomain)
            inline void global2local ( DofID global_id, DofID* local_id ) const;

          private:

            /// Create ownerships of dofs (needed for identification of sceleton dofs)
            void create_ownerships ( );
            /// Renumber dofs according to subdomains to achive global unique numbering
            void renumber ( );
            /// Create local and global correspondences
            void consecutive_numbering ( );

            /// Subdomain index
            int my_subdomain_;
            /// Offset for this subdomain
            int my_dof_offset_;
            /// Number of dofs including ghost layer
            int ndofs_incl_ghost_;
            /// Ghost subdomain indices, which are needed in case of point to point 
            /// communication
            std::vector<bool> shared_subdomains_;

            /// MPI Communicator
            MPI_Comm comm_;

            /// Total number of subdomains
            int number_of_subdomains_;

            /// Number of dofs for a specific subdomain owned by the subdomain
            std::vector<int> ndofs_on_sd_;
            /// Number of dofs on the global computational domain
            /// (including the other subdomains)
            int ndofs_global_;

            /// Vector of ownership
            std::vector<int> ownership_;

            /// Vector storing information for the mapping local 2 global
            std::vector<DofID> local2global_;
            /// Map storing information for the mapping global 2 local
            std::map<DofID, DofID> global2local_;
            /// Vector storing local (w.r.t. subdomain) numer_ field
            std::vector<DofID> numer_sd_;
        };

        // INLINE FUNCTIONS FOR DOFPARTITION

        template<class DataType>
        const MPI_Comm& DofPartition<DataType>::get_mpi_comm ( ) const
        {
            return comm_;
        }

        template<class DataType>
        int DofPartition<DataType>::get_my_subdomain ( ) const
        {
            return my_subdomain_;
        }

        template<class DataType>
        int DofPartition<DataType>::get_my_dof_offset ( ) const
        {
            return my_dof_offset_;
        }

        template<class DataType>
        int DofPartition<DataType>::ndofs_on_sd ( int subdomain ) const
        {
            assert ( subdomain >= 0 && subdomain < number_of_subdomains_ );
            return ndofs_on_sd_[subdomain];
        }

        template<class DataType>
        int DofPartition<DataType>::ndofs_global ( ) const
        {
            return ndofs_global_;
        }

        template<class DataType>
        void DofPartition<DataType>::local2global ( DofID local_id, DofID* global_id ) const
        {
            assert ( local_id >= 0 && local_id < ndofs_on_sd_[my_subdomain_] );
            *global_id = local_id + my_dof_offset_;
            //  *global_id = local2global_[local_id];
        }

        template<class DataType>
        void DofPartition<DataType>::global2local ( DofID global_id, DofID* local_id ) const
        {
            assert ( global_id >= my_dof_offset_ );
            assert ( global_id < my_dof_offset_ + ndofs_on_sd_[my_subdomain_] );
            *local_id = global_id - my_dof_offset_;
            //   std::map<DofID,DofID>::const_iterator it = global2local_.find(global_id);
            //   *local_id = (*it).second;
        }

        template<class DataType>
        bool DofPartition<DataType>::is_dof_on_sd ( DofID global_id ) const
        {
            if ( global_id >= my_dof_offset_
                 && global_id < my_dof_offset_ + ndofs_on_sd_[my_subdomain_] )
                return true;

            return false;
        }

        template<class DataType>
        int DofPartition<DataType>::owner_of_dof ( DofID global_id ) const
        {
            int result = 0;
            int bound = 0;

            for ( size_t s = 0, e_s = number_of_subdomains_; s != e_s; ++s )
            {
                bound += ndofs_on_sd_[s];
                if ( global_id < bound )
                {
                    result = s;
                    break;
                }
            }

            return result;
        }

    }
} // namespace hiflow

#endif
