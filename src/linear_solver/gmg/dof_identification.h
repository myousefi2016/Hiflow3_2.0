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

#ifndef DOF_IDENTIFICATION_H
#    define DOF_IDENTIFICATION_H

#    include <vector>
#    include "hybrid_map.h"
#    include "data_transfer.h"

#    include "boost/unordered_map.hpp"
#    include <stdexcept>
#    include <cassert>
#    include <limits>

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            namespace detail
            {
                /// A data structure that stores the data that will be transmitted between
                /// the processes and the two grids for each DoF to do the identification.
                /// This is: the DoF coordinates, the variable id to which the DoF belongs,
                /// and the process rank on which the DoF resides.

                struct DofIdentificationExchangedValues
                {
                    std::vector<double> coordinates;
                    int variable;
                    int rank;

                    /// Initializes the object with a dimension of zero. This constructor's
                    /// purpose is mainly to enable the use in some STL containers/function
                    /// that require a standard constructor. Unless you really know what
                    /// you are doing, you should not use this constructor explicitly.

                    DofIdentificationExchangedValues ( )
                    : dimensions_ ( 0 )
                    {
                    }

                    /// Initializes the exchanged_values object.
                    /// @param dimensions The dimensions of the problem. This is neccessary,
                    /// because the exchanged_values needs to know how many coordinate
                    /// dimensions exist.

                    explicit DofIdentificationExchangedValues ( unsigned dimensions )
                    : dimensions_ ( dimensions )
                    {
                        coordinates.reserve ( dimensions );
                    }

                    /// Fills the object with data. The data is expected to be formatted
                    /// in the way the append_to does it:
                    /// [coordinate0, coordinate1, ..., coordinateN, variable, rank]
                    /// @param source an InputIterator that specifies from the data
                    /// shall be read. The data type of the elements of the container
                    /// to which the iterator belongs must be convertible to double as
                    /// well as int. e.g, a std::vector<double> will do fine.

                    template<class InputIterator>
                    InputIterator load_from ( InputIterator source )
                    {
                        for ( unsigned i = 0; i < dimensions_; ++i )
                        {
                            coordinates.push_back ( static_cast < double > ( *source ) );
                            ++source;
                        }

                        variable = static_cast < int > ( *source );
                        ++source;

                        rank = static_cast < int > ( *source );
                        ++source;

                        return source;
                    }

                    /// Appends the content of the object to a container, in the same
                    /// format that the load_from method expects it:
                    /// [coordinate0, coordinate1, ..., coordinateN, variable, rank]
                    /// @param out The container to which the content of this object
                    /// shall be appended. It must support the push_back() operation
                    /// and double values must be convertible to data the data type
                    /// of its elements. E.g., a std::vector<double> will do fine.

                    template<class Container>
                    void append_to ( Container& out ) const
                    {
                        for ( unsigned i = 0; i < coordinates.size ( ); ++i )
                            out.push_back ( coordinates[i] );

                        out.push_back ( static_cast < double > ( variable ) );
                        out.push_back ( static_cast < double > ( rank ) );
                    }

                    unsigned size ( )
                    {
                        return number_of_values ( dimensions_ );
                    }

                    static inline unsigned number_of_values ( unsigned dimension )
                    {
                        return dimension + 2;
                    }

                  private:
                    unsigned dimensions_;
                };

                /// Compares two exchanged_values objects,
                /// and enables the use of exchanged_values objects as key in boost::unordered_map
                /// hashtables.
                /// The rank stored in the exchanged_values objects is not compared, because
                /// two exchanged_values objects shall be equal if their coordinates and
                /// variable ids are equal. On which process the DoFs are located is
                /// not relevant, because we want to know which DoFs correspond to each other.
                /// @return whether the DoFs tow which the exchanged_values objects refer
                /// correspond to each other, i.e. whether the coordinates and variable ids
                /// of the exchanged_values objects are equal.
                /// @param lhs left hand side exchanged_values object
                /// @param rhs right hand side exchanged_values object

                bool operator== ( const DofIdentificationExchangedValues& lhs,
                        const DofIdentificationExchangedValues& rhs );

                /// Hashes an exchanged_values object.
                /// This enables the use of exchanged_values objects as key in boost::unordered_map
                /// hashtables, which optimizes the coordinate comparison of two DoF lists.
                /// (see DofIdentification::generate_coordinate_intersection_map(...))
                /// The rank stored in an exchanged_values object is not hashed, because
                /// two exchanged_values objects shall be equal if their coordinates and
                /// variable ids are equal. On which process the DoFs are located is
                /// not relevant, because we want to know which DoFs correspond to each other.
                /// @return The hashed value for an exchanged_values object.
                /// @param data The exchanged_values object to hash.
                std::size_t hash_value ( const DofIdentificationExchangedValues& data );

            } // namespace detail

            /// Identifies which DoFs on two grids are at the same coordinate positions.

            template<class LAD>
            class DofIdentification
            {
              public:

                TYPE_FROM_CLASS ( LAD, DataType );

                typedef std::vector<DataType> Coord;

                /// Initializes the DoF Identification. DoF identifcation must
                /// be initiated explicitly by calling the identify_dofs() member
                /// function.
                /// @param fine_level The finer grid
                /// @param coarse_level The coarser grid. The processes in the communicator
                /// of the coarser grid must be a subset of the processes in the finer grid.
                /// @param info A pointer to a DataTransferInformation object for the two levels
                /// specified in fine_level and coarse_level
                /// @param transfer A pointer to a DataTransferFineToCoarse object
                /// that has been constructed with the fine and coarse levels specified
                /// in the fine_level and coarse_level arguments.

                DofIdentification ( BasicLevel<LAD>& fine_level,
                                    BasicLevel<LAD>& coarse_level,
                                    const boost::shared_ptr<const DataTransferInformation<LAD> >& info,
                                    const boost::shared_ptr<DataTransferFineToCoarse<LAD> >& transfer );

                ~DofIdentification ( )
                {
                }

                /// A data structure to store the id of a DoF and on which process
                /// rank this DoF is stored.

                struct id_map_element
                {
                    int id;
                    int rank;
                };

                typedef detail::DofIdentificationExchangedValues exchanged_values;

                /// This type stores a list of dof ids (cast to a double), each
                /// followed by the double vector representation of their respective
                /// exchanged_values object.
                typedef std::vector<double> DofListType;

                /// This is the map type that is used to check if local exchanged_values
                /// are also present in a received DoF list.
                typedef boost::unordered_map<exchanged_values, int> CoordinateMapType;

                /// This is the type of one of the final results of the DoF identification.
                /// This type will be used to map DoF ids on the fine grid to DoF ids
                /// on the coarse grid and the rank of the process on which these coarse
                /// DoF ids lie.
                typedef std::vector<std::pair<int, id_map_element> > IdListType;

                /// Runs the DoF identification. Collective on the processes of the
                /// communicator of the fine grid that has been specified in the
                /// constructor.
                void identify_dofs ( void );

                /// Check if a fine dof of the calling subdomain has been identified
                /// on the coarse grid
                /// @param fine_global_dof_id A global dof id on the fine grid
                /// @return whether a dof exists at the same position for the same variable
                /// on the coarse grid as well

                inline bool dof_exists_on_coarse_level ( int fine_global_dof_id ) const
                {
                    return matched_dofs_.query_bitmap ( fine_global_dof_id );
                }

                /// Returns a vector of coarse dof ids which are to be sent to a given process
                /// on the fine grid
                /// @param rank The process rank that is to be queried for dofs to transmit to
                /// @return a vector of dofs which specify which dofs to transfer

                const std::vector<int>& get_coarse_dofs_to_transfer_to ( int rank ) const
                {
                    assert ( rank >= 0 );
                    assert ( rank < coarse_dofs_to_transfer_.size ( ) );
                    return coarse_dofs_to_transfer_[rank];
                }

                /// @return A vector of vectors containing the coarse DoF ids
                /// that are to be transferred, with the i-th vector being the DoFs
                /// that shall be sent to the process with rank i on the communicator
                /// of the fine grid.

                const std::vector<std::vector<int> >& get_coarse_dofs_to_transfer ( ) const
                {
                    return coarse_dofs_to_transfer_;
                }

                /// Returns a vector of fine dof ids which are to be sent to a given process
                /// on the coarse grid
                /// @param rank The process rank that is to be queried for dofs to transmit to
                /// @return a vector of dofs which specify which dofs to transfer

                const std::vector<int>& get_fine_dofs_to_transfer_to ( int rank ) const
                {
                    assert ( rank >= 0 );
                    assert ( rank < fine_dofs_to_transfer_.size ( ) );
                    return fine_dofs_to_transfer_[rank];
                }

                /// @return A vector of vectors containing the fine DoF ids
                /// that are to be transferred to the coarse processes, with the i-th
                /// vector being the DoFs that shall be sent to the process with rank i
                /// on the communicator of the fine grid.

                const std::vector<std::vector<int> >& get_fine_dofs_to_transfer ( ) const
                {
                    return fine_dofs_to_transfer_;
                }

                const BasicLevel<LAD>& fine_level ( ) const
                {
                    return fine_level_;
                }

                const BasicLevel<LAD>& coarse_level ( ) const
                {
                    return coarse_level_;
                }

              private:

                /// Creates a vector of the DoFs that reside on this process, along with
                /// the information stored in an exchanged_values object for each DoF.
                /// Any DoFs on the coarse grid will be saved in coarse_dof_list_,
                /// DoFs on the fine grid will be saved in fine_dof_list_.

                void create_local_dof_lists ( )
                {
                    create_coarse_dof_list ( );
                    create_fine_dof_list ( );
                }

                /// Transmits the sizes of the local DoF lists (coarse_dof_list_) of all
                /// coarse processes to all fine processes.
                /// They will be saved in coarse_list_sizes_, with the i-th entry being
                /// size of the coarse DoF list of the process with rank i in the fine
                /// communicator.
                /// Additionally, sets global_coarse_list_size_ to the total number of
                /// coarse DoFs, i.e. the sum of the sizes received from all coarse processes.
                /// Collective on the fine communicator.
                void fetch_coarse_list_size ( void );

                /// Transmits the local coarse DoF lists from all coarse processes
                /// to all fine processes. Requires a call to fetch_coarse_dof_list_sizes
                /// beforehand.
                /// Collective on the fine communicator.
                /// @param out A map where the result of this function will be stored.
                /// After the function call, this map will map the exchanged_values objects
                /// to its corresponding coarse DoF ids. Thus, it will be possible to
                /// obtain the coarse DoF id that is located at a certain coordinate by
                /// querying the out map.
                void fetch_coarse_dof_list ( CoordinateMapType& out );

                /// Calculates which coarse processes are supposed to send which of their
                /// DoFs. It also fills the fine_dofs_to_transfer_ array, such that
                /// fine_dofs_to_transfer_[i] contains a vector of fine DoFs that are
                /// to be sent to the coarse process with rank i in the fine communicator.
                /// The fine_to_coarse_list_ object must have been filled before calling
                /// this function.
                /// @param ids The method will fill ids[i] with a vector of the coarse DoF ids
                /// that are to be transmitted by coarse process with fine rank i
                void compile_dof_transfer_information ( std::vector<std::vector<int> >& ids );

                /// Fills the local coarse DoF list coarse_dof_list_ with the coarse DoF ids
                /// followed by their exchanged_values object in its vector<double> representation.

                void create_coarse_dof_list ( )
                {
                    if ( process_is_in_coarse_level_ )
                        create_dof_list ( coarse_level_, coarse_dof_list_, coarse_rank_, true );
                }

                /// Fills the local coarse DoF list fine_dof_list_ with the fine DoF ids
                /// followed by their exchanged_values object in its vector<double> representation.

                void create_fine_dof_list ( )
                {
                    if ( process_is_in_fine_level_ )
                        create_dof_list ( fine_level_, fine_dof_list_, fine_rank_, true );
                }

                /// Fills a DofListType object with the local DoF list of a BasicLevel
                /// object. Each DofId (cast to double) will be followed by the corresponding
                /// exchanged_values object, serialized to vector<double>.
                /// @param level The level of which to create a DoF list
                /// @param out where the resulting DoF list will be stored.
                /// @param rank The rank of this process in the fine communicator.
                /// This rank will be saved in the exchanged_values object that will be
                /// put into the output DoF list.
                /// @param include_ghosts whether to include ghost DoFs in the list

                void create_dof_list ( const BasicLevel<LAD>& level,
                                       DofListType& out, const int rank, bool include_ghosts ) const;

                /// Parses a dof list, extracts the DoF ids and exchanged_values objects,
                /// and puts them into a map.
                /// @param array the DoF list to be parsed.
                /// @param result The resulting map. After a call to this function, result
                /// will map the exchanged_values objects to their corresponding DoF ids.

                void load_dof_list_into_map ( const DofListType& array,
                                              CoordinateMapType& result ) const;

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

                void generate_coordinate_intersection_map ( const DofListType& local_list,
                                                            const BasicLevel<LAD>& local_grid,
                                                            const CoordinateMapType& received_dofs,
                                                            util::HybridDofBitmap<DataType>& matched_dofs,
                                                            IdListType& result ) const;

                boost::shared_ptr<const DataTransferInformation<LAD> > info_;
                boost::shared_ptr<DataTransferFineToCoarse<LAD> > transfer_to_coarse_;

                DofListType fine_dof_list_;
                DofListType coarse_dof_list_;

                //this map holds the DoF identification. It is a map of type
                //(fine DoF id, {coarse DoF id, process rank})
                IdListType fine_to_coarse_list_;

                /// A bitmap where the i-th entry tells whether a dof with local id i
                /// on the fine grid has been identified on the coarse grid
                util::HybridDofBitmap<DataType> matched_dofs_;

                BasicLevel<LAD>& fine_level_;
                BasicLevel<LAD>& coarse_level_;

                bool process_is_in_coarse_level_;
                bool process_is_in_fine_level_;

                int coarse_rank_;
                int num_coarse_processes_;

                int fine_rank_;
                int num_fine_processes_;

                std::vector<int> coarse_list_sizes_;

                //This Object tells coarse processes which coarse DoF ids to transfer
                //to rank i
                std::vector<std::vector<int> > coarse_dofs_to_transfer_;
                // which fine DoFs to transfer to coarse process of rank i
                std::vector<std::vector<int> > fine_dofs_to_transfer_;

                int coarse_dof_list_size_;

                unsigned dimensions_;
                unsigned exchanged_values_size_;
            };

        } // namespace gmg
    } // namespace la
} // namespace hiflow

#endif
