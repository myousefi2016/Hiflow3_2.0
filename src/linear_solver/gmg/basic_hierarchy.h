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

#ifndef BASIC_HIERARCHY_H
#    define BASIC_HIERARCHY_H

#    include <mpi.h>
#    include <cassert>
#    include <stdexcept>

#    include "mesh/mesh_tools.h"
#    include "util.h"
#    include "gmg_level.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            namespace communicator_hierarchy_generators
            {
                /// Base class for all classes that generate a communicator hierarchy,
                /// i.e. the generation of one communicator per level of the MultiLevelHierarchy

                class BasicGenerator
                {
                  public:
                    /// @param root_comm The communicator from which to generate the hierarchy.

                    BasicGenerator ( MPI_Comm root_comm )
                    : root_comm_ ( root_comm )
                    {
                        assert ( root_comm_ != MPI_COMM_NULL );
                    }

                    virtual ~BasicGenerator ( )
                    {
                    }

                    /// Generates the hierarchy. Unless noted otherwise, the implementations
                    /// of this function can be assumed to be collective on the root communicator
                    /// supplied in the constructor.
                    /// Implementations have to make sure that after a call to this function,
                    /// the number of elements of the out vector must equal num_levels.
                    /// @param num_levels How many global_comms shall be generated.
                    /// Generally, this corresponds to the number of levels in the
                    /// MultiLevelHierarchy.
                    /// @param out std::vector in which the generated global_comms will be stored.
                    /// The first element is the largest communicator, while the following elements
                    /// are global_comms of decreasing size.
                    virtual void generate ( int num_levels,
                                            std::vector<MPI_Comm>& global_comms,
                                            std::vector<MPI_Comm>& partial_comms ) = 0;

                    /// @return The root communicator that was supplied in the constructor

                    MPI_Comm root_comm ( void ) const
                    {
                        return root_comm_;
                    }

                  protected:

                    /// @return a communicator containing only the first processes of a
                    /// communicator (as ordered by their rank). The ranks of the processes
                    /// in the new communicator should remain unchanged compared to
                    /// their ranks in the old communicator.
                    /// On processes that are not part of the returned communicator,
                    /// this method will return MPI_COMM_NULL.
                    /// @param comm The communicator from which to generate the returned
                    /// communicator
                    /// @param size Specifies how many processes will be in the
                    /// returned commmunicator.

                    inline MPI_Comm shrink_communicator_to_size ( MPI_Comm comm,
                                                                  int size ) const
                    {
                        if ( comm == MPI_COMM_NULL )
                            return MPI_COMM_NULL;

                        int rank = 0;
                        MPI_Comm_rank ( comm, &rank );

                        if ( rank >= size )
                            return split_communicator ( comm, MPI_UNDEFINED, rank );

                        return split_communicator ( comm, 0, rank );
                    }

                    /// A wrapper for MPI_Comm_split
                    /// @return the splitted communicator
                    /// @param comm The communicator to be splitted
                    /// @param color the color along which to split
                    /// @param key the key

                    inline MPI_Comm split_communicator ( MPI_Comm comm, int color, int key ) const
                    {
                        MPI_Comm result;
                        if ( MPI_Comm_split ( comm, color, key, &result ) == MPI_SUCCESS )
                            return result;
                        else throw std::runtime_error ( "could not split communicator" );
                    }

                  private:
                    MPI_Comm root_comm_;
                };

                class IndividualGrowthFactorGenerator : public BasicGenerator
                {
                  public:

                    IndividualGrowthFactorGenerator ( MPI_Comm root_comm, const std::vector<int>& factors )
                    : BasicGenerator ( root_comm ), factors_ ( factors )
                    {
                    }

                    virtual ~IndividualGrowthFactorGenerator ( )
                    {
                    }

                    virtual void generate ( int num_levels,
                                            std::vector<MPI_Comm>& global_comms,
                                            std::vector<MPI_Comm>& partial_comms );

                  protected:

                    std::vector<int> factors_;
                };

                class ConstantGrowthFactorGenerator
                : public IndividualGrowthFactorGenerator
                {
                  public:

                    /// @param root_comm The root communicator to be used by this generator.
                    /// All other global_comms will contain a subset of the processes contained
                    /// by the root communicator.
                    /// @param factor Specifies by how much the number of processes in
                    /// subsequent levels shall increase. If, for example root_comm contains
                    /// 10 processes and factor = 2, the generated global_comms
                    /// will contain 10, 5, 3, 2 processes.

                    ConstantGrowthFactorGenerator ( MPI_Comm root_comm, const int factor )
                    : IndividualGrowthFactorGenerator ( root_comm, std::vector<int>( 1, factor ) ), factor_ ( factor )
                    {
                        assert ( factor_ >= 1 );
                    }

                    virtual ~ConstantGrowthFactorGenerator ( )
                    {
                    }

                    virtual void generate ( int num_levels,
                                            std::vector<MPI_Comm>& global_comms,
                                            std::vector<MPI_Comm>& partial_comms );

                    const int factor_;
                };

            } // namespace communicator_hierarchy_generators

            /// Implements a hierarchy of arbitrary types

            template<class LevelType>
            class AbstractHierarchy
            {
              public:
                typedef typename std::vector<LevelType>::iterator IteratorFromCoarsest;
                typedef typename std::vector<LevelType>::const_iterator ConstIteratorFromCoarsest;

                typedef typename std::vector<LevelType>::reverse_iterator IteratorFromFinest;
                typedef typename std::vector<LevelType>::const_reverse_iterator ConstIteratorFromFinest;

                virtual ~AbstractHierarchy ( )
                {
                }

                /// @return An iterator to the coarsest level. Incrementing the iterator
                /// will traverse the hierarchy to the finer levels.

                IteratorFromCoarsest begin_from_coarsest_level ( )
                {
                    return levels_.begin ( );
                }

                /// @return A const_iterator to the coarsest level. Incrementing the iterator
                /// will traverse the hierarchy to the finer levels.

                ConstIteratorFromCoarsest begin_from_coarsest_level ( ) const
                {
                    return levels_.begin ( );
                }

                /// @return An iterator pointing one element beyond the finest level.

                IteratorFromCoarsest end_at_finest_level ( )
                {
                    return levels_.end ( );
                }

                /// @return A const_iterator pointing one element beyond the finest level.

                ConstIteratorFromCoarsest end_at_finest_level ( ) const
                {
                    return levels_.end ( );
                }

                /// @return An iterator to the finest level. Incrementing the iterator
                /// will traverse the hierarchy to the coarser levels.

                IteratorFromFinest begin_from_finest_level ( )
                {
                    return levels_.rbegin ( );
                }

                /// @return A const_iterator to the finest level. Incrementing the iterator
                /// will traverse the hierarchy to the coarser levels.

                ConstIteratorFromFinest begin_from_finest_level ( ) const
                {
                    return levels_.rbegin ( );
                }

                /// @return An iterator pointing one element beyond the coarsest level.

                IteratorFromFinest end_at_coarsest_level ( )
                {
                    return levels_.rend ( );
                }

                /// @return A const_iterator pointing one element beyond the coarsest level.

                ConstIteratorFromFinest end_at_coarsest_level ( ) const
                {
                    return levels_.rend ( );
                }

                /// @return The number of levels in this hierarchy.

                unsigned get_levels ( ) const
                {
                    return levels_.size ( );
                }

                /// @return A reference to a level.
                /// @param index Index of the level to be returned. Must be
                /// smaller than the number of levels stored.

                LevelType& get_level ( unsigned index )
                {
                    assert ( index < levels_.size ( ) );
                    return levels_[index];
                }

                /// @return A reference to a const level.
                /// @param index Index of the level to be returned. Must be
                /// smaller than the number of levels stored.

                const LevelType& get_level ( unsigned index ) const
                {
                    assert ( index < levels_.size ( ) );
                    return levels_[index];
                }

                /// @return The index to use with get_level() to access the coarsest level.

                inline unsigned get_index_of_coarsest ( ) const
                {
                    return 0;
                }

                /// @return The index to use with get_level() to access the finest level.

                inline unsigned get_index_of_finest ( ) const
                {
                    return levels_.size ( ) - 1;
                }

                /// This function is synonymous to get_level() const.
                /// @return A reference to a const level.
                /// @param index Index of the level to be returned. Must be
                /// smaller than the number of levels stored.

                const LevelType& operator[] ( unsigned index ) const
                {
                    return get_level ( index );
                }

                /// This function is synonymous to get_level().
                /// @return A reference to a level.
                /// @param index Index of the level to be returned. Must be
                /// smaller than the number of levels stored.

                LevelType& operator[] ( unsigned index )
                {
                    return get_level ( index );
                }

              protected:
                std::vector<LevelType> levels_;
            };

            /// Implements the basic multilevel hierarchy. The levels are stored in an array
            /// with the finest grid at index n_levels - 1 and the coarsest grid
            /// at index 0

            template<class LevelType>
            class BasicHierarchy : public AbstractHierarchy<LevelType>
            {
              public:

                IMPORT_FROM_BASECLASS ( AbstractHierarchy<LevelType>, IteratorFromFinest );
                IMPORT_FROM_BASECLASS ( AbstractHierarchy<LevelType>, IteratorFromCoarsest );
                IMPORT_FROM_BASECLASS ( AbstractHierarchy<LevelType>, ConstIteratorFromFinest );
                IMPORT_FROM_BASECLASS ( AbstractHierarchy<LevelType>, ConstIteratorFromCoarsest );

                /// Initializes the hierarchy. Collective on the root communicator that has
                /// been supplied to the communicator hierarchy generator.
                /// @param master_rank The rank of the master process in the
                /// root communicator.
                /// @param num_levels How many levels the hierarchy is supposed to have.
                /// @param master_mesh The mesh that shall be used as coarsest mesh in the hierarchy.
                /// It will be uniformly refined in successive steps until the desired number of levels is reached.
                /// Each mesh will be individually distributed among the processes in the corresponding MPI communicator of the hierarchy.
                /// This parameter only needs to be a valid pointer on the master process.
                /// @param fe_degrees Specifies the finite element polynomial degrees. The i-th entry of the
                /// vector is the polynomial degree of the i-th variable of the problem.
                /// @param settings The settings that are to be used by the linear algebra
                /// routines.
                /// @param comm_generator A pointer to a communicator hierarchy generator.

                BasicHierarchy ( unsigned master_rank,
                                 unsigned num_levels,
                                 const MeshPtr& master_mesh,
                                 const std::vector<int>& fe_degrees,
                                 const Settings& settings,
                                 communicator_hierarchy_generators::BasicGenerator* comm_generator )
                : comm_generator_ ( comm_generator ),
                number_levels_ ( num_levels ),
                master_rank_ ( master_rank ),
                master_mesh_ ( master_mesh )
                {
                    assert ( comm_generator != NULL );
                    assert ( num_levels >= 1 );
                    init ( fe_degrees, settings );
                }

                virtual ~BasicHierarchy ( )
                {
                }

                /// Resets all levels to their state before the assembly of the vectors
                /// and matrices.
                /// Collective on the root communicator of the communicator hierarchy.

                void reset_levels ( )
                {
                    for ( IteratorFromCoarsest it = this->begin_from_coarsest_level ( );
                          it != this->end_at_finest_level ( ); ++it )
                        it->reset_system ( );
                }

                /// Sets all solutions to zero.
                /// Collective on the root communicator of the communicator hierarchy.

                void reset_solutions ( )
                {
                    for ( IteratorFromCoarsest it = this->begin_from_coarsest_level ( );
                          it != this->end_at_finest_level ( ); ++it )
                        it->reset_solution ( );
                }

                /// @return The global assembler used by the hierarchy

                const typename LevelType::GlobalAssemblerType& get_global_assembler ( ) const
                {
                    return global_asm_;
                }
              private:

                /// Initializes the hierarchy. Collective on the root communicator supplied
                /// to the communicator hierarchy generator.
                /// @param fe_degrees The finite element degrees
                /// @param settings The linear algebra settings

                void init ( const std::vector<int>& fe_degrees, const Settings& settings )
                {
                    LOG_INFO ( "BasicHierarchy", "creating hierarchy with " << number_levels_ << " levels." );
                    std::vector<MPI_Comm> global_comms;
                    std::vector<MPI_Comm> partial_comms;

                    comm_generator_->generate ( number_levels_, global_comms, partial_comms );
                    assert ( global_comms.size ( ) == number_levels_ );
                    assert ( partial_comms.size ( ) == number_levels_ );

                    int rank_in_root_comm = 0;
                    MPI_Comm_rank ( global_comms.front ( ), &rank_in_root_comm );

                    // coarsest level index
                    int lvl = number_levels_ - 1;

                    // initial partitioning of coarsest mesh
                    MeshPtr local_mesh;
                    if ( global_comms.back ( ) != MPI_COMM_NULL )
                    {
                        local_mesh = partition_and_distribute ( master_mesh_, 0, global_comms.back ( ) );
                    }
                    master_mesh_.reset ( );

                    // create coarsest level
                    this->levels_.push_back ( LevelType ( global_comms[lvl],
                                                          partial_comms[lvl],
                                                          local_mesh,
                                                          master_rank_,
                                                          fe_degrees,
                                                          &global_asm_,
                                                          settings ) );

                    // other levels
                    for ( unsigned i = 1; i < number_levels_; ++i )
                    {
                        // this level index
                        --lvl;

                        if ( global_comms[lvl + 1] != MPI_COMM_NULL )
                        {
                            // proc is active on coarser level
                            local_mesh = local_mesh->refine ( );
                        }

                        if ( global_comms[lvl] != MPI_COMM_NULL )
                        {
                            // proc is active on this level
                            local_mesh = partition_and_distribute ( local_mesh, 0, partial_comms[lvl] );
                        }

                        this->levels_.push_back ( LevelType ( global_comms[lvl],
                                                              partial_comms[lvl],
                                                              local_mesh,
                                                              master_rank_,
                                                              fe_degrees,
                                                              &global_asm_,
                                                              settings ) );
                    }
                    assert ( lvl == 0 );
                }

                const unsigned number_levels_;
                const unsigned master_rank_;
                MeshPtr master_mesh_;

                typename LevelType::GlobalAssemblerType global_asm_;

                communicator_hierarchy_generators::BasicGenerator* comm_generator_;

            };

            /// Visualizes the solution vectors of all levels in multilevel hierarchy.
            /// Collective on the root communicator of the hierachy.
            /// @param hierarchy The multilevel hierarchy whose solutions shall be visualized.
            /// @param base_filename where the visualizations will be saved.
            /// the i-th level's solution will be stored at base_filename concatenated with i,
            /// where i = 0 corresponds to the coarsest level.

            template<class LevelType>
            void visualize_multilevel_solutions ( const BasicHierarchy<LevelType>& hierarchy,
                                                  const std::string& base_filename )
            {
                unsigned refinement_level = 0;

                for ( typename BasicHierarchy<LevelType>::ConstIteratorFromCoarsest it
                      = hierarchy.begin_from_coarsest_level ( );
                      it != hierarchy.end_at_finest_level ( ); ++it )
                {
                    if ( it->is_scheduled_to_this_process ( ) )
                        it->visualize_solution ( base_filename + boost::lexical_cast<std::string>( refinement_level ) );

                    ++refinement_level;
                }
            }

        } // namespace gmg
    } // namespace la
} // namespace hiflow

#endif
