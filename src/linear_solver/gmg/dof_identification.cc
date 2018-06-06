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

/// \author Aksel Alpay, Martin Wlotzka

#include "dof_identification.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            namespace detail
            {

                bool operator== ( const DofIdentificationExchangedValues& lhs,
                        const DofIdentificationExchangedValues& rhs )
                {
                    double rel_epsilon = 100.0 * std::numeric_limits<double>::epsilon ( );

                    if ( lhs.coordinates.size ( ) != rhs.coordinates.size ( ) )
                        return false;

                    for ( std::size_t i = 0; i < lhs.coordinates.size ( ); ++i )
                        if ( std::abs ( lhs.coordinates[i] - rhs.coordinates[i] ) >
                             rel_epsilon * std::abs ( lhs.coordinates[i] ) )
                            return false;

                    if ( lhs.variable != rhs.variable )
                        return false;
                    //We do not check for the rank, because where it lies is not
                    //important for the equality of two dofs

                    return true;
                }

                std::size_t hash_value ( const DofIdentificationExchangedValues& data )
                {
                    boost::hash<std::vector<float> > coordinate_hasher;
                    std::vector<float> float_coordinates ( data.coordinates.size ( ) );
                    for ( std::size_t i = 0; i < data.coordinates.size ( ); ++i )
                        float_coordinates[i] = static_cast < float > ( data.coordinates[i] );

                    std::size_t result = coordinate_hasher ( float_coordinates );

                    boost::hash<int> int_hasher;

                    //We do not hash the rank
                    result ^= int_hasher ( data.variable );

                    return result;
                }

            } // namespace detail

            template<class LAD>
            DofIdentification<LAD>::DofIdentification ( BasicLevel<LAD>& fine_level,
                                                        BasicLevel<LAD>& coarse_level,
                                                        const boost::shared_ptr<const DataTransferInformation<LAD> >& info,
                                                        const boost::shared_ptr<DataTransferFineToCoarse<LAD> >& transfer )
            : fine_level_ ( fine_level ), coarse_level_ ( coarse_level ),
            coarse_rank_ ( -1 ), fine_rank_ ( -1 ), num_coarse_processes_ ( -1 ),
            num_fine_processes_ ( -1 ),
            coarse_dof_list_size_ ( 0 )
            {
                process_is_in_coarse_level_ = coarse_level_.is_scheduled_to_this_process ( );
                process_is_in_fine_level_ = fine_level_.is_scheduled_to_this_process ( );

                info_ = info;

                if ( process_is_in_coarse_level_ )
                {
                    // Processes in the coarse grid have to be a subset
                    // of the processes in the fine grid
                    assert ( process_is_in_fine_level_ );

                    assert ( fine_level.space ( )->get_dim ( ) ==
                             coarse_level.space ( )->get_dim ( ) );
                }

                if ( process_is_in_fine_level_ )
                {
                    dimensions_ = fine_level.space ( )->get_dim ( );

                    coarse_rank_ = info_->coarse_rank ( );
                    fine_rank_ = info_->fine_rank ( );
                    num_coarse_processes_ = info_->num_coarse_procs ( );
                    num_fine_processes_ = info_->num_fine_procs ( );

                }
                transfer_to_coarse_ = transfer;

                exchanged_values_size_ = exchanged_values::number_of_values ( dimensions_ );

                assert ( num_coarse_processes_ <= num_fine_processes_ );
            }

            template<class LAD>
            void DofIdentification<LAD>::identify_dofs ( )
            {
                if ( process_is_in_fine_level_ )
                {
                    // First, create lists of the DoFs residing on this process,
                    // both for the fine and coarse grids.
                    create_local_dof_lists ( );

                    // Next distribute the sizes of the DoF lists of the coarse processes
                    // to all fine processes
                    fetch_coarse_list_size ( );

                    // Distribute the DoF lists of all coarse processes to all fine
                    // processes. The result is stored as a hash table to speed
                    // up processing.
                    CoordinateMapType coarse_maps;
                    fetch_coarse_dof_list ( coarse_maps );

                    //       if (process_is_in_fine_level_)
                    // compute which fine DoFs equal which of the coarse DoFs
                    // that have just been fetched by comparing coordinates
                    // and variable ids. The result will be stored in
                    // fine_to_coarse_list_ (it maps fine DoF id -> (coarse DoF id,
                    //       rank of the coarse process in the fine communicator))
                    generate_coordinate_intersection_map ( fine_dof_list_,
                                                           fine_level_,
                                                           coarse_maps,
                                                           matched_dofs_,
                                                           fine_to_coarse_list_ );

                    // empty the local lists, they are not needed anymore.
                    fine_dof_list_.clear ( );
                    coarse_dof_list_.clear ( );
                    coarse_list_sizes_.clear ( );

                    std::vector<std::vector<int> > compiled_coarse_dof_information;

                    // Iterate the fine_to_coarse_map_, and create a lists of DoFs
                    // to transfer for every coarse process.
                    // This function call fills the fine_dofs_to_transfer_ object
                    // such, that the i-th entry contains a list of fine DoFs that have
                    // to be transferred to the coarse process with coarse rank i.
                    // Furthermore, after the call compiled_coarse_dof_information
                    // will be filled such that the i-th entry contains a list of
                    // coarse DoFs that have to be transferred to this process
                    // by the coarse process of rank i.
                    compile_dof_transfer_information ( compiled_coarse_dof_information );

                    /// This list is not needed anymore
                    fine_to_coarse_list_.clear ( );

                    // The final step: Tell the coarse processes
                    // which of their coarse DoFs they are supposed to transfer
                    // to this fine process. This is the information we have
                    // generated with compile_dof_transfer_information()
                    transfer_to_coarse_->transfer_vectors ( compiled_coarse_dof_information,
                                                            coarse_dofs_to_transfer_ );

                }
            }

            template<class LAD>
            void DofIdentification<LAD>::fetch_coarse_list_size ( void )
            {
                if ( process_is_in_fine_level_ )
                {
                    coarse_dof_list_size_ = 0;
                    if ( process_is_in_coarse_level_ )
                    {
                        assert ( fine_level_.partial_rank ( ) == 0 );
                        coarse_dof_list_size_ = coarse_dof_list_.size ( );
                    }

                    MPI_Bcast ( &coarse_dof_list_size_, 1, MPI_INT, 0, fine_level_.partial_comm ( ) );
                    if ( !process_is_in_coarse_level_ )
                    {
                        assert ( fine_level_.partial_rank ( ) != 0 );
                        coarse_dof_list_.resize ( coarse_dof_list_size_ );
                    }
                }
            }

            /// Transmits the local coarse DoF lists from all coarse processes
            /// to all fine processes. Requires a call to fetch_coarse_dof_list_sizes
            /// beforehand.
            /// Collective on the fine communicator.
            /// @param out A map where the result of this function will be stored.
            /// After the function call, this map will map the exchanged_values objects
            /// to its corresponding coarse DoF ids. Thus, it will be possible to
            /// obtain the coarse DoF id that is located at a certain coordinate by
            /// querying the out map.

            template<class LAD>
            void DofIdentification<LAD>::fetch_coarse_dof_list ( CoordinateMapType& out )
            {
                if ( process_is_in_fine_level_ )
                {
                    MPI_Bcast ( util::raw_array ( coarse_dof_list_ ), coarse_dof_list_size_, MPI_DOUBLE,
                                0, fine_level_.partial_comm ( ) );

                    load_dof_list_into_map ( coarse_dof_list_, out );
                }
            }

            /// Calculates which coarse processes are supposed to send which of their
            /// DoFs. It also fills the fine_dofs_to_transfer_ array, such that
            /// fine_dofs_to_transfer_[i] contains a vector of fine DoFs that are
            /// to be sent to the coarse process with rank i in the fine communicator.
            /// The fine_to_coarse_list_ object must have been filled before calling
            /// this function.
            /// @param ids The method will fill ids[i] with a vector of the coarse DoF ids
            /// that are to be transmitted by coarse process with fine rank i

            template<class LAD>
            void DofIdentification<LAD>::compile_dof_transfer_information ( std::vector<std::vector<int> >& ids )
            {
                if ( process_is_in_fine_level_ )
                {
                    //     ids = std::vector<std::vector<int> >(num_fine_processes_);
                    //     fine_dofs_to_transfer_ = std::vector<std::vector<int> >(
                    //             num_fine_processes_);
                    ids = std::vector<std::vector<int> >( num_coarse_processes_ );
                    fine_dofs_to_transfer_ = std::vector<std::vector<int> >(
                            num_coarse_processes_ );

                    for ( typename IdListType::const_iterator it = fine_to_coarse_list_.begin ( );
                          it != fine_to_coarse_list_.end ( ); ++it )
                    {

                        //insert coarse id
                        ids[it->second.rank].push_back ( it->second.id );
                        //insert fine id into the fine_dofs_to_transfer array
                        fine_dofs_to_transfer_[it->second.rank].push_back ( it->first );
                    }
                }
            }

            template<class LAD>
            void DofIdentification<LAD>::create_dof_list ( const BasicLevel<LAD>& level,
                                                           DofListType& out, const int rank,
                                                           bool include_ghosts ) const
            {
                const DofPartition<DataType>& dof = level.space ( )->dof ( );

                out.clear ( );
                const int rank_on_lvl = level.rank ( );
                const int num_dofs_on_lvl = dof.ndofs_on_sd ( rank_on_lvl );
                out.reserve ( num_dofs_on_lvl * ( 1 + exchanged_values_size_ ) ); // id, coords, variable, rank

                const int tdim = level.mesh ( )->tdim ( );
                const int ncells = level.mesh ( )->num_entities ( tdim );

                const int num_variables = level.space ( )->get_nb_var ( );

                // We assume 15% of all dofs are ghosts (which is very pessimistic),
                // but it's better to use somewhat more memory than to cause a rehash
                // of the hash table due to overloading it...
                util::HybridDofBitmap<DataType> processed_dofs ( 0.15 * num_dofs_on_lvl, &dof );

                for ( int cellid = 0; cellid < ncells; ++cellid )
                {
                    for ( int var = 0; var < num_variables; ++var )
                    {
                        std::vector<Coord> coords;

                        dof.get_coord_on_cell ( var, cellid, coords );

                        std::vector<DofID> ids;
                        dof.get_dofs_on_cell ( var, cellid, ids );
                        assert ( coords.size ( ) == ids.size ( ) );

                        for ( int i = 0; i < coords.size ( ); ++i )
                        {
                            if ( include_ghosts || dof.is_dof_on_sd ( ids[i] ) )
                            {
                                // only do work if we haven't processed the dof already
                                if ( !processed_dofs.query_bitmap ( ids[i] ) )
                                {
                                    exchanged_values values ( dimensions_ );

                                    std::vector<double> double_coordinates;
                                    for ( std::size_t pos = 0; pos < coords.at ( i ).size ( ); ++pos )
                                        double_coordinates.push_back ( static_cast < double > ( coords.at ( i )[pos] ) );

                                    values.coordinates = double_coordinates;

                                    //we always work with the fine rank
                                    values.rank = rank;
                                    values.variable = var;

                                    out.push_back ( static_cast < double > ( ids[i] ) );
                                    values.append_to ( out );

                                    processed_dofs.mark_dof ( ids[i] );
                                }
                            }
                        }
                    }
                }
            }

            /// Parses a dof list, extracts the DoF ids and exchanged_values objects,
            /// and puts them into a map.
            /// @param array the DoF list to be parsed.
            /// @param result The resulting map. After a call to this function, result
            /// will map the exchanged_values objects to their corresponding DoF ids.

            template<class LAD>
            void DofIdentification<LAD>::load_dof_list_into_map ( const DofListType& array,
                                                                  CoordinateMapType& result ) const
            {
                // make sure the array has a valid size - it is supposed to store
                // a sequence of dof ids with their exchanged_values object.
                // Thus, we need to make sure that the number of elements is a
                // multiple of exchanged_values_size_ + 1
                assert ( array.size ( ) % ( exchanged_values_size_ + 1 ) == 0 );

                unsigned num_dofs = array.size ( ) / ( exchanged_values_size_ + 1 );

                // We use 20% more buckets than we have elements to avoid collisions
                // in the hash table
                result = CoordinateMapType ( 1.2 * num_dofs );

                for ( std::vector<double>::const_iterator it = array.begin ( );
                      it != array.end ( ); )
                {
                    int id = static_cast < int > ( *it );
                    ++it;

                    exchanged_values values ( dimensions_ );
                    it = values.load_from ( it );

                    result.insert ( std::pair<exchanged_values, int>( values, id ) );
                }
            }

            /// Compares a local DoF list with received DoF lists and computes
            /// the DoFs that are present in both lists by comparing their coordinates
            /// and variable ids. This process is optimized by the choice of
            /// a hash table for one of the lists.
            /// @param local_list the lcoal DoF List to be compared
            /// @param received_dofs A map of dofs as created by load_dof_list_into_map()
            /// @param result After a call to this function, this will pair
            /// all fine DoF ids of this process with their corresponding coarse dof ids
            /// and the rank of the coarse process on which the coarse dofs are located.
            /// Ghost DoFs will not be included in \c result.
            /// @param local_grid The grid to which \c local_list belongs
            /// @param matched_dofs a vector<bool> bitmap that stores which DoFs had
            /// a match. If a DoF was found to exist in both lists,
            /// \c matched_dofs[local id of DoF] will be set to true. Ghost DoFs
            /// will be stored in a hash table of the \c dof_match_table

            template<class LAD>
            void DofIdentification<LAD>::generate_coordinate_intersection_map ( const DofListType& local_list,
                                                                                const BasicLevel<LAD>& local_grid,
                                                                                const CoordinateMapType& received_dofs,
                                                                                util::HybridDofBitmap<DataType>& matched_dofs,
                                                                                IdListType& result ) const
            {
                //result: vector<pair<fineid, (coarse id, rank)> >

                assert ( local_list.size ( ) % ( exchanged_values_size_ + 1 ) == 0 );
                unsigned num_elements = local_list.size ( ) / ( exchanged_values_size_ + 1 );

                result.clear ( );

                const DofPartition<DataType>& dof = local_grid.space ( )->dof ( );

                unsigned num_local_elements = dof.ndofs_on_sd ( local_grid.rank ( ) );

                result.reserve ( num_local_elements );

                matched_dofs = util::HybridDofBitmap<DataType>( num_elements - num_local_elements,
                        &( local_grid.space ( )->dof ( ) ) );

                for ( DofListType::const_iterator it = local_list.begin ( );
                      it != local_list.end ( ); )
                {
                    // Load dof id and exchanged_values from the local dof list
                    int dof_id = static_cast < int > ( *it );
                    ++it;
                    exchanged_values local_values ( dimensions_ );
                    it = local_values.load_from ( it );

                    typename CoordinateMapType::const_iterator search_iterator =
                            received_dofs.find ( local_values );

                    if ( search_iterator != received_dofs.end ( ) )
                    {
                        // Only process DoF Id if it is not a ghost
                        if ( dof.is_dof_on_sd ( dof_id ) )
                        {
                            id_map_element elem;

                            // The search iterator points to a pair<exchanged_values, dof id>
                            // Thus, the dof id is at search_iterator->second
                            elem.id = search_iterator->second;
                            // The rank can then be accessed via the the exchanged_values
                            // object at search_iterator->first
                            elem.rank = search_iterator->first.rank;

                            result.push_back ( std::pair<int, id_map_element>( dof_id, elem ) );
                        }

                        // mark dof as matched
                        matched_dofs.mark_dof ( dof_id );
                    }
                }
            }

            template class DofIdentification<LADescriptorCoupledD>;
            template class DofIdentification<LADescriptorCoupledS>;

        } // namespace gmg
    } // namespace la
} // namespace hiflow
