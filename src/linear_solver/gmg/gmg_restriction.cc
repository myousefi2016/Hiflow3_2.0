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

#include "gmg_restriction.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            /// Prepares a cell of a vector for restriction by adding contributions from dofs that
            /// are ignored during the vector transfer onto the transferred dofs.

            template<class LAD>
            void LinearCellRestrictor<LAD>::operator() ( const Element<DataType>& element,
                    MatrixType* matrix,
                    std::set<int>& unchanged_dofs_set )
            {
                unsigned tdim = element.get_cell ( ).tdim ( );

                const DofIdConverter<LAD>& dof_converter = on_demand_dof_converter_.get_converter ( &element );
                const DofPartition<DataType>& dof_partition = element.space ( ).dof ( );

                std::vector<int> dofs;

                for ( int var = 0; var < element.get_num_variables ( ); ++var )
                {
                    std::vector<bool> is_dof_processed ( element.get_num_dofs ( var ), false );

                    // skip dofs on dirichlet boundary
                    const int n = dirichlet_dofs_.size ( );
                    for ( int i = 0; i < n; ++i )
                    {
                        int gid = dirichlet_dofs_[i];
                        int lid;
                        if ( dof_converter.map_g2l ( gid, lid ) )
                        {
                            assert ( lid >= 0 );
                            assert ( lid < is_dof_processed.size ( ) );
                            is_dof_processed[lid] = true;
                        }
                    }

                    // vertices, edges, faces
                    for ( unsigned dim = 0; dim < tdim; ++dim )
                    {
                        int sindex = 0;
                        for ( IncidentEntityIterator ent = element.get_cell ( ).begin_incident ( dim );
                              ent != element.get_cell ( ).end_incident ( dim );
                              ++ent, ++sindex )
                        {
                            // obtain dofs on the subentity
                            dof_partition.get_dofs_on_subentity ( var,
                                                                  element.get_cell_index ( ),
                                                                  dim,
                                                                  sindex,
                                                                  dofs );

                            // Check if the neighboring cells have already processed the dof
                            for ( IncidentEntityIterator neighbor_cell = ent->begin_incident ( tdim );
                                  neighbor_cell != ent->end_incident ( tdim );
                                  ++neighbor_cell )
                            {
                                if ( is_cell_already_processed ( element.get_cell ( ), *neighbor_cell ) )
                                {
                                    // Mark all dofs on this entity as processed
                                    for ( int i = 0; i < dofs.size ( ); ++i )
                                    {
                                        int current_local_id = -1;
                                        if ( dof_converter.map_g2l ( dofs[i], current_local_id ) )
                                        {
                                            assert ( current_local_id < is_dof_processed.size ( ) );
                                            is_dof_processed[current_local_id] = true;
                                        }
                                    }
                                }
                            }

                            for ( int i = 0; i < dofs.size ( ); ++i )
                            {
                                int current_local_id = -1;
                                if ( dof_converter.map_g2l ( dofs[i], current_local_id ) )
                                {
                                    if ( !is_dof_processed[current_local_id] )
                                    {
                                        process_dof ( is_dof_processed,
                                                      dof_partition,
                                                      dof_converter,
                                                      dofs[i],
                                                      var,
                                                      *ent,
                                                      matrix );
                                    }
                                }
                            }
                        }
                    }

                    // the cell itself
                    element.get_dof_indices ( var, dofs );

                    for ( int i = 0; i < dofs.size ( ); ++i )
                    {
                        // dofs which do not exist on the coarse grid stay unchanged
                        if ( !dof_ident_->dof_exists_on_coarse_level ( dofs[i] ) )
                        {
                            unchanged_dofs_set.insert ( dofs[i] );
                        }
                    }

                    for ( int i = 0; i < dofs.size ( ); ++i )
                        if ( !is_dof_processed[i] )
                            process_dof ( is_dof_processed,
                                          dof_partition,
                                          dof_converter,
                                          dofs[i],
                                          var,
                                          element.get_cell ( ),
                                          matrix );
                }
            }

            /// Prepares a cell for restriction
            /// @param element The cell that shall be prepared. It must stem from the vector
            /// that has been supplied in the constructor.

            template<class LAD>
            void LinearCellRestrictor<LAD>::operator() ( const Element<DataType>& element )
            {
                unsigned tdim = element.get_cell ( ).tdim ( );

                const DofIdConverter<LAD>& dof_converter = on_demand_dof_converter_.get_converter ( &element );
                const DofPartition<DataType>& dof_partition = element.space ( ).dof ( );

                std::vector<int> dofs;

                for ( int var = 0; var < element.get_num_variables ( ); ++var )
                {
                    std::vector<bool> is_dof_processed ( element.get_num_dofs ( var ), false );

                    // skip dofs on dirichlet boundary
                    const int n = dirichlet_dofs_.size ( );
                    for ( int i = 0; i < n; ++i )
                    {
                        int gid = dirichlet_dofs_[i];
                        int lid;
                        if ( dof_converter.map_g2l ( gid, lid ) )
                        {
                            assert ( lid >= 0 );
                            assert ( lid < is_dof_processed.size ( ) );
                            is_dof_processed[lid] = true;
                        }
                    }

                    for ( unsigned dim = 0; dim < tdim; ++dim )
                    {
                        int sindex = 0;
                        for ( IncidentEntityIterator ent = element.get_cell ( ).begin_incident ( dim );
                              ent != element.get_cell ( ).end_incident ( dim );
                              ++ent, ++sindex )
                        {
                            // obtain dofs
                            dof_partition.get_dofs_on_subentity ( var,
                                                                  element.get_cell_index ( ),
                                                                  dim,
                                                                  sindex,
                                                                  dofs );

                            // Check if the neighboring cells have already processed the dof
                            for ( IncidentEntityIterator neighbor_cell = ent->begin_incident ( tdim );
                                  neighbor_cell != ent->end_incident ( tdim );
                                  ++neighbor_cell )
                            {
                                if ( is_cell_already_processed ( element.get_cell ( ), *neighbor_cell ) )
                                {
                                    // Mark all dofs on this entity as processed
                                    for ( int i = 0; i < dofs.size ( ); ++i )
                                    {
                                        int current_local_id = -1;
                                        if ( dof_converter.map_g2l ( dofs[i], current_local_id ) )
                                        {
                                            assert ( current_local_id < is_dof_processed.size ( ) );
                                            is_dof_processed[current_local_id] = true;
                                        }
                                    }
                                }
                            }

                            for ( int i = 0; i < dofs.size ( ); ++i )
                            {
                                int current_local_id = -1;
                                if ( dof_converter.map_g2l ( dofs[i], current_local_id ) )
                                {
                                    if ( !is_dof_processed[current_local_id] )
                                    {
                                        process_dof ( is_dof_processed,
                                                      dof_partition,
                                                      dof_converter,
                                                      dofs[i],
                                                      var,
                                                      *ent );
                                    }
                                }
                            }
                        }
                    }
                    element.get_dof_indices ( var, dofs );

                    for ( int i = 0; i < dofs.size ( ); ++i )
                        if ( !is_dof_processed[i] )
                            process_dof ( is_dof_processed,
                                          dof_partition,
                                          dof_converter,
                                          dofs[i],
                                          var,
                                          element.get_cell ( ) );
                }
            }

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

            template<class LAD>
            void LinearCellRestrictor<LAD>::process_dof ( std::vector<bool>& processed_dofs,
                                                          const DofPartition<DataType>& dof_partition,
                                                          const DofIdConverter<LAD>& dof_converter,
                                                          const int dof,
                                                          const int var,
                                                          const Entity& ent,
                                                          MatrixType* matrix )
            {
                int cell_dim = dof_converter.get_element ( )->get_cell ( ).tdim ( );

                int local_id = -1;
                if ( dof_converter.map_g2l ( dof, local_id ) )
                {
                    if ( dof_ident_->dof_exists_on_coarse_level ( dof ) )
                    {
                        if ( dof_partition.is_dof_on_sd ( dof ) )
                        {
                            std::vector<int> contributing_dofs;
                            std::vector<DataType> contributing_weights;

                            contributing_dofs.push_back ( dof );
                            contributing_weights.push_back ( 1.0 );

                            // Get contributions and weights of surrounding dofs
                            get_contributions_from_cell ( dof, var, *dof_converter.get_element ( ),
                                                          dof_converter,
                                                          contributing_dofs,
                                                          contributing_weights );

                            // If we are at a cell border we need to consider the neighboring cells
                            if ( ent.tdim ( ) < cell_dim )
                            {
                                for ( IncidentEntityIterator cell = ent.begin_incident ( cell_dim );
                                      cell != ent.end_incident ( cell_dim ); ++cell )
                                {
                                    if ( cell->index ( ) != dof_converter.get_element ( )->get_cell_index ( ) )
                                    {
                                        Element<DataType> neighbor_cell ( dof_converter.get_element ( )->space ( ),
                                                                          cell->index ( ) );

                                        DofIdConverter<LAD> neighbor_converter ( &neighbor_cell );

                                        get_contributions_from_cell ( dof, var, neighbor_cell,
                                                                      neighbor_converter,
                                                                      contributing_dofs,
                                                                      contributing_weights );
                                    }
                                }
                            }

                            std::map<int, DataType> sorted;
                            for ( int i = 0; i < contributing_dofs.size ( ); ++i )
                                sorted.insert ( std::make_pair ( contributing_dofs[i], contributing_weights[i] ) );

                            contributing_dofs.clear ( );
                            contributing_weights.clear ( );

                            for ( typename std::map<int, DataType>::const_iterator it = sorted.begin ( );
                                  it != sorted.end ( ); ++it )
                            {
                                contributing_dofs.push_back ( it->first );
                                contributing_weights.push_back ( it->second );
                            }
                            matrix->Add ( &dof, 1,
                                          vec2ptr ( contributing_dofs ), contributing_dofs.size ( ),
                                          vec2ptr ( contributing_weights ) );
                        }
                    }
                    processed_dofs[local_id] = true;
                }
            }

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

            template<class LAD>
            void LinearCellRestrictor<LAD>::process_dof ( std::vector<bool>& processed_dofs,
                                                          const DofPartition<DataType>& dof_partition,
                                                          const DofIdConverter<LAD>& dof_converter,
                                                          const int dof,
                                                          const int var,
                                                          const Entity& ent )
            {
                int cell_dim = dof_converter.get_element ( )->get_cell ( ).tdim ( );

                int local_id = -1;
                if ( dof_converter.map_g2l ( dof, local_id ) )
                {
                    if ( dof_ident_->dof_exists_on_coarse_level ( dof ) )
                    {
                        if ( dof_partition.is_dof_on_sd ( dof ) )
                        {
                            std::vector<int> contributing_dofs;

                            DataType contributions = 0.0;
                            DataType weight_sum = 0.0;

                            // Get contributions and weights of surrounding dofs
                            contributions += get_contributions_from_cell ( dof, var, *dof_converter.get_element ( ),
                                                                           dof_converter,
                                                                           contributing_dofs,
                                                                           weight_sum );

                            // If we are at a cell border we need to consider the neighboring cells
                            if ( ent.tdim ( ) < cell_dim )
                            {
                                for ( IncidentEntityIterator cell = ent.begin_incident ( cell_dim );
                                      cell != ent.end_incident ( cell_dim ); ++cell )
                                {
                                    if ( cell->index ( ) != dof_converter.get_element ( )->get_cell_index ( ) )
                                    {
                                        Element<DataType> neighbor_cell ( dof_converter.get_element ( )->space ( ),
                                                                          cell->index ( ) );

                                        DofIdConverter<LAD> neighbor_converter ( &neighbor_cell );

                                        DataType weights = 0.0;
                                        contributions += get_contributions_from_cell ( dof, var, neighbor_cell,
                                                                                       neighbor_converter,
                                                                                       contributing_dofs,
                                                                                       weights );
                                        weight_sum += weights;
                                    }
                                }
                            }

                            // Add previous value
                            DataType previous_value = 0.0;
                            vector_->GetValues ( &dof, 1, &previous_value );

                            contributions += previous_value;
                            weight_sum += 1.0;

                            //std::cout << "dof = " << dof << " prev = " << previous_value << " contributions = " << contributions << " weight = " << weight_sum << std::endl;

                            DataType new_value = contributions;

                            vector_->SetValues ( &dof, 1, &new_value );
                            //std::cout << "---------------------\n";
                        }
                    }
                    processed_dofs[local_id] = true;
                }
            }

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

            template<class LAD>
            void LinearCellRestrictor<LAD>::get_contributions_from_cell ( int dof_id,
                                                                          int var,
                                                                          const Element<DataType>& cell,
                                                                          const DofIdConverter<LAD>& dof_converter,
                                                                          std::vector<int>& contributing_dofs,
                                                                          std::vector<DataType>& contributing_weights ) const
            {
                int local_id = -1;
                if ( dof_converter.map_g2l ( dof_id, local_id ) )
                {
                    ijk_dof_id ijk_id;
                    if ( dof_converter.map_l2ijk ( local_id, var, ijk_id ) )
                    {
                        CellType::Tag t = cell.get_cell ( ).cell_type ( ).tag ( );

                        switch ( t )
                        {
                            case CellType::QUADRILATERAL:
                                get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas2d_on_axis,
                                                                  ijk_id, var, dof_converter, 0.5, contributing_dofs, contributing_weights );

                                get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas2d_diagonals,
                                                                  ijk_id, var, dof_converter, 0.25, contributing_dofs, contributing_weights );
                                break;
                            case CellType::TRIANGLE:
                                get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::PDeltas::deltas2d,
                                                                  ijk_id, var, dof_converter, 0.5, contributing_dofs, contributing_weights );

                                break;
                            case CellType::HEXAHEDRON:
                                get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas3d_on_axis,
                                                                  ijk_id, var, dof_converter, 0.5, contributing_dofs, contributing_weights );

                                get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas3d_two_axis_diagonals,
                                                                  ijk_id, var, dof_converter, 0.25, contributing_dofs, contributing_weights );

                                get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas3d_three_axis_diagonals,
                                                                  ijk_id, var, dof_converter, 0.125, contributing_dofs, contributing_weights );
                                break;
                            case CellType::TETRAHEDRON:
                                //TODO
                                throw std::runtime_error ( "Restriction Preparator: Tetrahedrons not yet supported." );
                                break;
                            case CellType::PYRAMID:
                                //TODO
                                throw std::runtime_error ( "Restriction Preparator: Pyramids not yet supported." );
                                break;
                            case CellType::POINT:
                                //TODO
                                throw std::runtime_error ( "Restriction Preparator: Points not yet supported." );
                                break;
                            case CellType::LINE:
                                //TODO
                                throw std::runtime_error ( "Restriction Preparator: Lines not yet supported." );
                                break;
                            default:
                                throw std::runtime_error ( "Not a implemented element." );
                                break;
                        }
                    }
                }
            }

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

            template<class LAD>
            typename LinearCellRestrictor<LAD>::DataType LinearCellRestrictor<LAD>::get_contributions_from_cell ( int dof_id,
                                                                                                                  int var,
                                                                                                                  const Element<DataType>& cell,
                                                                                                                  const DofIdConverter<LAD>& dof_converter,
                                                                                                                  std::vector<int>& contributing_dofs,
                                                                                                                  DataType& weight_sum ) const
            {
                weight_sum = 0.0;

                DataType result = 0.0;

                int local_id = -1;
                if ( dof_converter.map_g2l ( dof_id, local_id ) )
                {
                    ijk_dof_id ijk_id;
                    if ( dof_converter.map_l2ijk ( local_id, var, ijk_id ) )
                    {
                        CellType::Tag t = cell.get_cell ( ).cell_type ( ).tag ( );

                        DataType weight_temp = 0.0;
                        switch ( t )
                        {
                            case CellType::QUADRILATERAL:
                                result += get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas2d_on_axis,
                                                                            ijk_id, var, dof_converter, 0.5, contributing_dofs, weight_temp );
                                weight_sum += weight_temp;

                                result += get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas2d_diagonals,
                                                                            ijk_id, var, dof_converter, 0.25, contributing_dofs, weight_temp );
                                weight_sum += weight_temp;

                                break;
                            case CellType::TRIANGLE:
                                result += get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::PDeltas::deltas2d,
                                                                            ijk_id, var, dof_converter, 0.5, contributing_dofs, weight_sum );

                                break;
                            case CellType::HEXAHEDRON:
                                result += get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas3d_on_axis,
                                                                            ijk_id, var, dof_converter, 0.5, contributing_dofs, weight_temp );
                                weight_sum += weight_temp;

                                result += get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas3d_two_axis_diagonals,
                                                                            ijk_id, var, dof_converter, 0.25, contributing_dofs, weight_temp );
                                weight_sum += weight_temp;

                                result += get_contributions_of_delta_list ( NearestKnownDofsFinder<LAD>::QDeltas::deltas3d_three_axis_diagonals,
                                                                            ijk_id, var, dof_converter, 0.125, contributing_dofs, weight_temp );
                                weight_sum += weight_temp;

                                break;
                            case CellType::TETRAHEDRON:
                                //TODO
                                throw std::runtime_error ( "Restriction Preparator: Tetrahedrons not yet supported." );
                                break;
                            case CellType::PYRAMID:
                                //TODO
                                throw std::runtime_error ( "Restriction Preparator: Pyramids not yet supported." );
                                break;
                            case CellType::POINT:
                                //TODO
                                throw std::runtime_error ( "Restriction Preparator: Points not yet supported." );
                                break;
                            case CellType::LINE:
                                //TODO
                                throw std::runtime_error ( "Restriction Preparator: Lines not yet supported." );
                                break;
                            default:
                                throw std::runtime_error ( "Not a implemented element." );
                                break;
                        }

                        return result;
                    }
                }

                return 0.0;
            }

            /// @return the number of expected contributions of a dof located in a given cell
            /// @param cell Specifies the cell type of which in turn is used to determine
            /// expected number of contributions.

            template<class LAD>
            unsigned LinearCellRestrictor<LAD>::get_expected_num_contributions ( const Element<DataType>& cell ) const
            {
                CellType::Tag t = cell.get_cell ( ).cell_type ( ).tag ( );

                switch ( t )
                {
                    case CellType::QUADRILATERAL:
                        return 8;
                        break;
                    case CellType::TRIANGLE:
                        return 6;
                        break;
                    case CellType::HEXAHEDRON:
                        return 26;
                        break;
                    case CellType::TETRAHEDRON:
                        return 12;
                        break;
                    default:
                        throw std::runtime_error ( "Elements are not supported yet." );
                        break;
                }
            }

            template class LinearCellRestrictor<LADescriptorCoupledD>;
            template class LinearCellRestrictor<LADescriptorCoupledS>;

            template class LinearRestriction<LADescriptorCoupledD>;
            template class LinearRestriction<LADescriptorCoupledS>;

        } // namespace gmg
    } // namespace la
} // namespace hiflow
