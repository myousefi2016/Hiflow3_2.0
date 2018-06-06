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

#ifndef GMG_RESTRICTION_H
#    define GMG_RESTRICTION_H

#    include "gmg_interpolation.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            /// Prepares a cell of a vector for restriction by adding contributions from dofs that
            /// are ignored during the vector transfer onto the transferred dofs.

            template<class LAD>
            class LinearCellRestrictor
            {
              public:

                TYPE_FROM_CLASS ( LAD, MatrixType );
                TYPE_FROM_CLASS ( LAD, VectorType );
                TYPE_FROM_CLASS ( LAD, DataType );

                /// Construct object
                /// @param dof_ident A pointer to a dof identification object between the level
                /// on which the vector lives that shall be restricted and the next coarser level
                /// @param vector A pointer to the vector that shall be prepared for restriction

                LinearCellRestrictor ( const boost::shared_ptr<const DofIdentification<LAD> >& dof_ident,
                                       VectorType* vector,
                                       const std::vector<int>& dirichlet_dofs
                                       )
                : dof_ident_ ( dof_ident ), vector_ ( vector ), dirichlet_dofs_ ( dirichlet_dofs )
                {
                    assert ( dof_ident );
                    assert ( vector );
                }

                /// Prepares the restriction matrix.
                /// @param element The cell that shall be used to derive the corresponding entries in the
                /// restriction matrix.
                /// @param matrix The matrix to be prepared.
                void operator() ( const Element<DataType>& element,
                        MatrixType* matrix,
                        std::set<int>& unchanged_dofs_set );

                /// Prepares a cell for restriction
                /// @param element The cell that shall be prepared. It must stem from the vector
                /// that has been supplied in the constructor.
                void operator() ( const Element<DataType>& element );

              private:

                /// Checks if a dof on an entity is a dof that must be modified, and if
                /// so, calculates and applies the contributions from the surrounding dofs.
                /// @param processed_dofs A bitmap indicating which dofs on the cell have already
                /// been processed, such that a dof with local id i is located at
                /// \c processed_dofs[i]. If this method modifies the supplied dof, its
                /// entry in the \c processed_dofs bitmap will be set to \c true.
                /// @param dof_partition The dof partition
                /// @param dof_converter A dof converter object for the cell in which the dof
                /// is located
                /// @param dof The global dof id of the dof to be processed
                /// @param var The variable id of the dof
                /// @param ent The (sub)entity on which the dof is located
                inline void process_dof ( std::vector<bool>& processed_dofs,
                                          const DofPartition<DataType>& dof_partition,
                                          const DofIdConverter<LAD>& dof_converter,
                                          const int dof,
                                          const int var,
                                          const Entity& ent,
                                          MatrixType* matrix );

                /// Checks if a dof on an entity is a dof that must be modified, and if
                /// so, calculates and applies the contributions from the surrounding dofs.
                /// @param processed_dofs A bitmap indicating which dofs on the cell have already
                /// been processed, such that a dof with local id i is located at
                /// \c processed_dofs[i]. If this method modifies the supplied dof, its
                /// entry in the \c processed_dofs bitmap will be set to \c true.
                /// @param dof_partition The dof partition
                /// @param dof_converter A dof converter object for the cell in which the dof
                /// is located
                /// @param dof The global dof id of the dof to be processed
                /// @param var The variable id of the dof
                /// @param ent The (sub)entity on which the dof is located
                inline void process_dof ( std::vector<bool>& processed_dofs,
                                          const DofPartition<DataType>& dof_partition,
                                          const DofIdConverter<LAD>& dof_converter,
                                          const int dof,
                                          const int var,
                                          const Entity& ent );

                /// Calculates all contributions to a dof from a given cell.
                /// @param dof_id The global dof id of the dof to be processed
                /// @param var The variable id of the dof
                /// @param cell The cell in which the dof is located
                /// @param dof_converter A dof converter object for this cell
                /// @param weight_sum After a successful call, will contain sum of the weights
                /// of all contributions
                /// @param contributing_dofs The method will append the global ids of the dofs that contributed
                /// to this vector
                /// @param contributing_weights The method will append the weights for contributing dofs.
                /// @return The sum of the weighted contributions from the surrounding dofs
                inline void get_contributions_from_cell ( int dof_id,
                                                          int var,
                                                          const Element<DataType>& cell,
                                                          const DofIdConverter<LAD>& dof_converter,
                                                          std::vector<int>& contributing_dofs,
                                                          std::vector<DataType>& contributing_weights ) const;

                /// Calculates all contributions to a dof from a given cell.
                /// @param dof_id The global dof id of the dof to be processed
                /// @param var The variable id of the dof
                /// @param cell The cell in which the dof is located
                /// @param dof_converter A dof converter object for this cell
                /// @param weight_sum After a successful call, will contain sum of the weights
                /// of all contributions
                /// @param contributing_dofs The method will append the global ids of the dofs that contributed
                /// to this vector
                /// @return The sum of the weighted contributions from the surrounding dofs
                inline DataType get_contributions_from_cell ( int dof_id,
                                                              int var,
                                                              const Element<DataType>& cell,
                                                              const DofIdConverter<LAD>& dof_converter,
                                                              std::vector<int>& contributing_dofs,
                                                              DataType& weight_sum ) const;

                /// Processes a list of ijk deltas, checks if the corresponding dofs are located
                /// in a given cell, and if so, adds their values and calculates their weight.
                /// @param delta_list An array of ijk deltas that will be applied onto the
                /// the origin ijk dof id to look for surrounding dofs
                /// @param origin The ijk dof id to which the ijk deltas refer, i.e. the dof
                /// that shall be prepared for restriction
                /// @param var The variable id of the dof
                /// @param dof_converter A dof converter object for this cell
                /// @param weight_of_dof The weight that shall be applied to dofs found via the
                /// delta list
                /// @param contributing_dofs The method will append the global ids of the dofs that contributed
                /// to this vector
                /// @param summed_result_weight After a successful call, will contain the sum
                /// of the weights of the dofs that have contributed
                /// @return The sum of the contributions of the dofs

                template<unsigned N>
                inline void get_contributions_of_delta_list ( const int (&delta_list )[N][3],
                                                              const ijk_dof_id& origin,
                                                              int var,
                                                              const DofIdConverter<LAD>& dof_converter,
                                                              const DataType weight_of_dof,
                                                              std::vector<int>& contributing_dofs,
                                                              std::vector<DataType>& contributing_weights ) const
                {
                    assert ( origin.size ( ) >= 2 && origin.size ( ) <= 3 );

                    ijk_dof_id current_id = origin;

                    for ( unsigned i = 0; i < N; ++i )
                    {
                        for ( unsigned j = 0; j < origin.size ( ); ++j )
                            current_id[j] = origin[j] + delta_list[i][j];

                        // check if dof exists and obtain global id
                        int local_id = -1;
                        if ( dof_converter.map_ijk2l ( current_id, var, local_id ) )
                        {
                            int global_id = -1;
                            if ( dof_converter.map_l2g ( local_id, var, global_id ) )
                            {
                                int lid;
                                dof_converter.map_ijk2l ( origin, var, lid );
                                int gid;
                                dof_converter.map_l2g ( lid, var, gid );
                                //           std::cout << "for dof " << origin[0] << "," << origin[1] << " with gid " << gid << " found contributing dof " << current_id[0] << "," << current_id[1] << " with gid " << global_id << " and weight " << weight_of_dof << std::endl;
                                // if dof has not yet contributed
                                if ( std::find ( contributing_dofs.begin ( ), contributing_dofs.end ( ), global_id )
                                     == contributing_dofs.end ( ) )
                                {
                                    //             std::cout << "for dof " << origin[0] << "," << origin[1] << " with gid " << gid
                                    //             << " found contributing dof " << current_id[0] << "," << current_id[1] << " with gid " << global_id
                                    //             << " and weight " << weight_of_dof << std::endl;
                                    contributing_dofs.push_back ( global_id );
                                    contributing_weights.push_back ( weight_of_dof );
                                }
                            }
                        }
                    }
                }

                /// Processes a list of ijk deltas, checks if the corresponding dofs are located
                /// in a given cell, and if so, adds their values and calculates their weight.
                /// @param delta_list An array of ijk deltas that will be applied onto the
                /// the origin ijk dof id to look for surrounding dofs
                /// @param origin The ijk dof id to which the ijk deltas refer, i.e. the dof
                /// that shall be prepared for restriction
                /// @param var The variable id of the dof
                /// @param dof_converter A dof converter object for this cell
                /// @param weight_of_dof The weight that shall be applied to dofs found via the
                /// delta list
                /// @param contributing_dofs The method will append the global ids of the dofs that contributed
                /// to this vector
                /// @param summed_result_weight After a successful call, will contain the sum
                /// of the weights of the dofs that have contributed
                /// @return The sum of the contributions of the dofs

                template<unsigned N>
                inline DataType get_contributions_of_delta_list ( const int (&delta_list )[N][3],
                                                                  const ijk_dof_id& origin,
                                                                  int var,
                                                                  const DofIdConverter<LAD>& dof_converter,
                                                                  const DataType weight_of_dof,
                                                                  std::vector<int>& contributing_dofs,
                                                                  DataType& summed_result_weight ) const
                {
                    //std::cout << "Investigating cell " << dof_converter.get_element()->get_cell_index() << std::endl;

                    assert ( origin.size ( ) >= 2 && origin.size ( ) <= 3 );

                    ijk_dof_id current_id = origin;
                    DataType sum_contributions = 0.0;
                    summed_result_weight = 0.0;

                    for ( unsigned i = 0; i < N; ++i )
                    {
                        for ( unsigned j = 0; j < origin.size ( ); ++j )
                            current_id[j] = origin[j] + delta_list[i][j];

                        // check if dof exists and obtain global id
                        int local_id = -1;
                        if ( dof_converter.map_ijk2l ( current_id, var, local_id ) )
                        {
                            int global_id = -1;
                            if ( dof_converter.map_l2g ( local_id, var, global_id ) )
                            {
                                // if dof has not yet contributed
                                if ( std::find ( contributing_dofs.begin ( ), contributing_dofs.end ( ), global_id )
                                     == contributing_dofs.end ( ) )
                                {
                                    // obtain dof value
                                    DataType value = 0.0;
                                    vector_->GetValues ( &global_id, 1, &value );

                                    sum_contributions += value * weight_of_dof;
                                    summed_result_weight += weight_of_dof;

                                    //std::cout << "Adding " << weight_of_dof << " from " << dof_converter.get_element()->get_cell_index() << std::endl;

                                    contributing_dofs.push_back ( global_id );
                                }
                            }
                        }

                    }

                    return sum_contributions;
                }

                /// @return whether a given cell has already been processed
                /// @param current_cell The cell that is currently processed
                /// @param other The cell that shall be checked

                inline bool is_cell_already_processed ( const Entity& current_cell, const Entity& other ) const
                {
                    return other.index ( ) < current_cell.index ( );
                }

                /// @return the number of expected contributions of a dof located in a given cell
                /// @param cell Specifies the cell type of which in turn is used to determine
                /// expected number of contributions.
                inline unsigned get_expected_num_contributions ( const Element<DataType>& cell ) const;

                boost::shared_ptr<const DofIdentification<LAD> > dof_ident_;
                VectorType* vector_;
                OnDemandDofIdConverter<LAD> on_demand_dof_converter_;
                const std::vector<int>& dirichlet_dofs_;
            };

            /// Prepares a vector for restriction. The vector will be modified in the process.

            template<class LAD>
            class LinearRestriction
            {
              public:

                TYPE_FROM_CLASS ( LAD, VectorType );

                LinearRestriction<LAD>( )
                : use_restriction_matrix_ ( false )
                {
                }

                LinearRestriction<LAD>( const bool use )
                : use_restriction_matrix_ ( use )
                {
                }

                void use_restriction_matrix ( const bool use )
                {
                    use_restriction_matrix_ = use;
                }

                template<class ConnectedLevelType>
                void use_transposed_of_interpolation_matrix ( ConnectedLevelType& lvl )
                {
                    if ( lvl.is_scheduled_to_this_process ( ) )
                    {
                        assert ( lvl.interpolation_matrix ( ) );

                        lvl.create_restriction_matrix ( );
                        lvl.restriction_matrix ( )->CreateTransposedFrom ( *( lvl.interpolation_matrix ( ).get ( ) ) );
                    }
                }

                template<class ConnectedLevelType>
                void build_restriction_matrix ( ConnectedLevelType& lvl )
                {
                    if ( lvl.is_scheduled_to_this_process ( ) )
                    {
                        lvl.create_restriction_matrix ( );

                        lvl.restriction_matrix ( )->Zeros ( );

                        ForEachCell<LAD> for_each_cell ( lvl.mesh ( ), lvl.space ( ) );

                        LinearCellRestrictor<LAD> cell_preparator (
                                                                    lvl.get_connection_to_next_coarser_grid ( )->get_dof_identification ( ),
                                                                    lvl.res ( ).get ( ),
                                                                    lvl.dirichlet_dofs ( ) );

                        std::set<int> unchanged_dofs_set ( lvl.dirichlet_dofs ( ).begin ( ),
                                                           lvl.dirichlet_dofs ( ).end ( ) );

                        for_each_cell ( cell_preparator, lvl.restriction_matrix ( ), unchanged_dofs_set );

                        std::vector<int> unchanged_dofs_vec;

                        for ( std::set<int>::const_iterator it = unchanged_dofs_set.begin ( );
                              it != unchanged_dofs_set.end ( ); ++it )
                        {
                            if ( lvl.space ( )->dof ( ).is_dof_on_sd ( *it ) )
                            {
                                unchanged_dofs_vec.push_back ( *it );
                            }
                        }

                        lvl.restriction_matrix ( )->diagonalize_rows ( vec2ptr ( unchanged_dofs_vec ),
                                                                       unchanged_dofs_vec.size ( ), 1.0 );

                        lvl.restriction_matrix ( )->Compress ( );
                    }
                }

                /// Run the preparation algorithm.
                /// @param lvl The level in the hierarchy on which the vector is located
                /// @param vector The vector that shall be prepared for restriction. It must
                /// originate from the level given by \c lvl

                template<class ConnectedLevelType>
                void operator() ( const ConnectedLevelType& lvl,
                        VectorType& vector,
                        const std::vector<int>& dirichlet_dofs
                        ) const
                {
                    if ( lvl.is_scheduled_to_this_process ( ) )
                    {
                        if ( use_restriction_matrix_ )
                        {
                            assert ( lvl.tmp_vector ( ) );
                            assert ( lvl.restriction_matrix ( ) );

                            lvl.restriction_matrix ( )->VectorMult ( vector, lvl.tmp_vector ( ).get ( ) );
                            vector.CopyFrom ( *( lvl.tmp_vector ( ) ) );
                        }
                        else
                        {
                            vector.Update ( );

                            ForEachCell<LAD> for_each_cell ( lvl.mesh ( ), lvl.space ( ) );

                            LinearCellRestrictor<LAD> cell_preparator (
                                                                        lvl.get_connection_to_next_coarser_grid ( )->get_dof_identification ( ),
                                                                        &vector,
                                                                        dirichlet_dofs );

                            for_each_cell ( cell_preparator );
                        }
                    }
                }

              private:
                bool use_restriction_matrix_;
            };

        } // namespace gmg
    } // namespace la
} // namespace hiflow

#endif
