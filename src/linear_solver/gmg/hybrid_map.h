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

/// \author Aksel Alpay

#ifndef HYBRID_MAP_H
#    define HYBRID_MAP_H

#    include <cassert>
#    include "boost/unordered_set.hpp"

#    include "util.h"
#    include "dof/dof_partition.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {
            namespace util
            {

                /// Implements a divider for the hybrid_map to distinguish between local
                /// and non-local dofs using their id.

                template<class DataType>
                class dof_locality_index_space_divider
                {
                  public:

                    /// Construct object
                    /// @param dof_partition A pointer to a \c DofPartition object that shall
                    /// be used to query dof ids.

                    dof_locality_index_space_divider ( const doffem::DofPartition<DataType>* dof_partition )
                    : dof_partition_ ( dof_partition )
                    {
                        assert ( dof_partition != NULL );
                    }

                    /// Construct object

                    dof_locality_index_space_divider ( )
                    : dof_partition_ ( NULL )
                    {
                    }

                    /// Checks if a dof is non-local, i.e. cannot be represented by the continous
                    /// local dof id enumeration (i.e. the dof is "scattered")
                    /// @return whether the supplied dof id is scattered
                    /// @param global_dof the global dof id that shall be checked

                    bool is_scattered ( int global_dof ) const
                    {
                        assert ( dof_partition_ );
                        return !dof_partition_->is_dof_on_sd ( global_dof );
                    }

                    /// @return The index in the continous local enumeration of the supplied dof.
                    /// The supplied dof id must not be scattered.
                    /// @param global_dof The global dof id of the dof

                    int get_bitmap_index ( int global_dof ) const
                    {
                        assert ( dof_partition_ );
                        assert ( !is_scattered ( global_dof ) );

                        int local_dof = 0;
                        dof_partition_->global2local ( global_dof, &local_dof );

                        return local_dof;
                    }
                  private:

                    const doffem::DofPartition<DataType>* dof_partition_;
                };

                /// Implements a "hybrid bitmap". For indices from a continous range, a bitmap
                /// is used, for scattered indices a set is used. Whether an index is from
                /// the continous range is determined by the \c IndexSpaceDivider specidied
                /// (e.g. a dof_localility_index_space_divider)

                template<class IndexSpaceDivider>
                class hybrid_bitmap
                {
                    typedef boost::unordered_set<int> SetType;

                    SetType scattered_entries_;
                    std::vector<bool> continous_entries_;

                    IndexSpaceDivider divider_;
                  public:
                    /// Initialize object
                    /// @param expected_num_scattered_entries The expected number of non-contiguous entries
                    /// @param num_contiguous_entries The number of contiguous entries
                    /// @param space_divider A Index space divider object that specifies if an index
                    /// is contiguous (i.e. will be stored in the bitmap) or scattered (i.e. will be
                    /// stored in a set)

                    hybrid_bitmap ( std::size_t expected_num_scattered_entries,
                                    std::size_t num_contiguous_entries,
                                    const IndexSpaceDivider& space_divider )
                    : scattered_entries_ ( 1.2 * expected_num_scattered_entries ),
                    divider_ ( space_divider )
                    {
                        continous_entries_ = std::vector<bool>(
                                num_contiguous_entries,
                                false );
                    }

                    /// Lookup an entry in the bitmap
                    /// @param idx The index that shall be looked up
                    /// @return whether the index is marked in the hybrid bitmap

                    bool query_bitmap ( int idx ) const
                    {
                        if ( !divider_.is_scattered ( idx ) )
                        {
                            int continous_idx = divider_.get_bitmap_index ( idx );
                            assert ( 0 <= continous_idx && continous_idx < continous_entries_.size ( ) );
                            return continous_entries_[continous_idx];
                        }
                        else return scattered_entries_.find ( idx ) != scattered_entries_.end ( );
                    }

                    /// Mark an index in the hybrid bitmap
                    /// @param idx The index that shall be marked

                    void mark_dof ( int idx )
                    {
                        if ( !divider_.is_scattered ( idx ) )
                        {
                            int continous_idx = divider_.get_bitmap_index ( idx );
                            assert ( 0 <= continous_idx && continous_idx < continous_entries_.size ( ) );
                            continous_entries_[continous_idx] = true;
                        }
                        else scattered_entries_.insert ( idx );
                    }

                    /// Unmark index in hybrid bitmap
                    /// @param index The index that shall be unmarked

                    void unmark_dof ( int idx )
                    {
                        if ( !divider_.is_scattered ( ) )
                        {
                            int continous_idx = divider_.get_bitmap_index ( idx );
                            assert ( 0 <= continous_idx && continous_idx < continous_entries_.size ( ) );
                            continous_entries_[continous_idx] = false;
                        }
                        else
                        {
                            SetType::iterator it = scattered_entries_.find ( idx );
                            if ( it != scattered_entries_.end ( ) )
                                scattered_entries_.erase ( it );
                        }
                    }
                };

                /// Implements a hybrid bitmap of dof indices.

                template<class DataType>
                class HybridDofBitmap : public hybrid_bitmap<dof_locality_index_space_divider<DataType> >
                {
                  public:
                    /// Construct object
                    /// @param num_ghosts The estimated number of ghost dofs on the local subdomain.
                    /// @param dof_partition A pointer to a \c dof_partition object to look up dofs

                    HybridDofBitmap ( std::size_t num_ghosts,
                                      const doffem::DofPartition<DataType>* dof_partition )
                    : hybrid_bitmap<dof_locality_index_space_divider<DataType> >(
                    num_ghosts,
                    dof_partition->ndofs_on_sd ( dof_partition->get_my_subdomain ( ) ),
                    dof_locality_index_space_divider<DataType>( dof_partition ) )
                    {
                    }

                    /// Construct object. Objects created by this constructor will be unable
                    /// to perform any meaningful operations.

                    HybridDofBitmap ( )
                    : hybrid_bitmap<dof_locality_index_space_divider<DataType> >(
                    0,
                    0,
                    dof_locality_index_space_divider<DataType>( ) )
                    {
                    }

                };

            } // namespace util
        } // namespace gmg
    } // namespace la
} // namespace hiflow

#endif
