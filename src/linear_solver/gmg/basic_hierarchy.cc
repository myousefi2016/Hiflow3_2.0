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

#include "basic_hierarchy.h"

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
                /*
                void ConstantReductionFactorGenerator::generate_communicators_of_all_sizes(std::vector<MPI_Comm>& global_comms) const
                {
                  MPI_Comm comm = root_comm();

                  global_comms.clear();

                  int size = 0;
                  MPI_Comm_size(comm, &size);

                  // Insert the root communicator
                  global_comms.push_back(comm);

                  // Make sure the factor is smaller than 1.0 to prevent
                  // an infinite loop in the following while loop
                  if (reduction_factor_ < 1.0)
                  {
                    // repeat as long as we have more than one process in
                    // the communicator
                    while (static_cast<int> (reduction_factor_ * size) > 1)
                    {
                      comm = shrink_communicator_to_size(global_comms.back(),
                                                          static_cast<int> (reduction_factor_ * size));

                      global_comms.push_back(comm);

                      size *= reduction_factor_;
                    }

                    // Make sure we include a final communicator of size 1
                    global_comms.push_back(shrink_communicator_to_size(root_comm(), 1));
                  }
                }*/

                /*
                  std::vector<MPI_Comm> available_communicators;
                  generate_communicators_of_all_sizes(available_communicators);

                  // do the global comms
                  global_comms.resize(num_levels, MPI_COMM_NULL);

                  if (num_levels <= available_communicators.size())
                  {
                    // We require less communicators than we have available,
                    // so there's no need to fill and we just copy
                    // as many communicators as desired.
                    std::copy(available_communicators.begin(),
                              available_communicators.begin() + num_levels,
                              global_comms.begin());
                  }
                  else
                  {
                    int num_fills = num_levels -
                            available_communicators.size();

                    // Fill the first elements of the result vector
                    // with available_communicators.front() which is the
                    // root communicator.
                    std::fill(global_comms.begin(), global_comms.begin() + num_fills,
                              available_communicators.front());

                    // Copy the remaining elements from our available communicators.
                    std::copy(available_communicators.begin(),
                              available_communicators.end(),
                              global_comms.begin() + num_fills);
                  }

                  // do the partial comms

                }*/

                // void BackFilledConstantReductionFactorGenerator::generate(int num_levels,
                //                                                           std::vector<MPI_Comm>& global_comms,
                //                                                           std::vector<MPI_Comm>& partial_comms)
                // {
                //   std::vector<MPI_Comm> available_communicators;
                //   generate_communicators_of_all_sizes(available_communicators);
                //
                //   global_comms.resize(num_levels, MPI_COMM_NULL);
                //
                //   if (num_levels <= available_communicators.size())
                //   {
                //     std::copy(available_communicators.begin(),
                //               available_communicators.begin() + num_levels,
                //               global_comms.begin());
                //   }
                //   else
                //   {
                //     std::copy(available_communicators.begin(),
                //               available_communicators.end(),
                //               global_comms.begin());
                //
                //     std::fill(global_comms.begin() + available_communicators.size(),
                //               global_comms.end(),
                //               available_communicators.back());
                //   }
                // }

                void IndividualGrowthFactorGenerator::generate ( int num_levels,
                                                                 std::vector<MPI_Comm>& global_comms,
                                                                 std::vector<MPI_Comm>& partial_comms )
                {
                    assert ( num_levels >= 1 );
                    assert ( factors_.size ( ) >= num_levels - 1 );

                    global_comms.clear ( );
                    global_comms.reserve ( num_levels );

                    // finest global comm is root comm
                    global_comms.push_back ( this->root_comm ( ) );

                    // partial comms vector
                    partial_comms.clear ( );
                    partial_comms.reserve ( num_levels );

                    int last_rank, last_size;

                    for ( int l = 0; l < num_levels - 1; ++l )
                    {
                        assert ( factors_[l] > 0 );

                        if ( global_comms.back ( ) != MPI_COMM_NULL )
                        {
                            // this proc is active on finer level
                            MPI_Comm_rank ( global_comms.back ( ), &last_rank );
                            MPI_Comm_size ( global_comms.back ( ), &last_size );

                            // derive partial comm for finer level
                            MPI_Comm partial_comm;
                            int color = last_rank / factors_[l];
                            MPI_Comm_split ( global_comms.back ( ), color, last_rank, &partial_comm );
                            partial_comms.push_back ( partial_comm );

                            // derive global comm for next coarser level
                            MPI_Comm global_comm;
                            if ( last_rank % factors_[l] == 0 )
                            {
                                MPI_Comm_split ( global_comms.back ( ), 0, last_rank, &global_comm );
                            }
                            else
                            {
                                MPI_Comm_split ( global_comms.back ( ), MPI_UNDEFINED, last_rank, &global_comm );
                            }
                            global_comms.push_back ( global_comm );
                        }
                        else
                        {
                            // proc is not active on finer level
                            global_comms.push_back ( MPI_COMM_NULL );
                            partial_comms.push_back ( MPI_COMM_NULL );
                        }
                    }

                    // partial comm for coarsest level is the global comm
                    partial_comms.push_back ( global_comms.back ( ) );
                }

                void ConstantGrowthFactorGenerator::generate ( int num_levels,
                                                               std::vector<MPI_Comm>& global_comms,
                                                               std::vector<MPI_Comm>& partial_comms )
                {
                    assert ( num_levels >= 1 );

                    this->factors_.clear ( );
                    this->factors_.resize ( num_levels, factor_ );

                    IndividualGrowthFactorGenerator::generate ( num_levels, global_comms, partial_comms );
                }

            } // namespace communicator_hierarchy_generators

            template class BasicHierarchy<BasicLevel<LADescriptorCoupledD> >;
            template class BasicHierarchy<BasicLevel<LADescriptorCoupledS> >;

            template void visualize_multilevel_solutions<BasicLevel<LADescriptorCoupledD> >( const BasicHierarchy<BasicLevel<LADescriptorCoupledD> >&, const std::string& );
            template void visualize_multilevel_solutions<BasicLevel<LADescriptorCoupledS> >( const BasicHierarchy<BasicLevel<LADescriptorCoupledS> >&, const std::string& );

        } // namespace multigrid
    } // namespace la
} // namespace hiflow
