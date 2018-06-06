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

#include "gmg_interpolation.h"

namespace hiflow
{
    namespace la
    {
        namespace gmg
        {

            void QuadBasedIjkIterator::next ( )
            {
                if ( is_end_ )
                    return;

                if ( current_ijk_id_[0] < fe_degree_ )
                    ++( current_ijk_id_[0] );
                else
                {
                    if ( current_ijk_id_[1] < fe_degree_ )
                    {
                        ++( current_ijk_id_[1] );
                        current_ijk_id_[0] = 0;
                    }
                    else
                    {
                        if ( is_three_dimensional_ )
                        {
                            if ( current_ijk_id_[2] < fe_degree_ )
                            {
                                ++( current_ijk_id_[2] );
                                current_ijk_id_[0] = 0;
                                current_ijk_id_[1] = 0;
                            }
                            else
                            {
                                is_end_ = true;
                                return;
                            }
                        }
                        else
                        {
                            //We have reached the final dof, and cannot go further.
                            is_end_ = true;
                            return;
                        }
                    }
                }

                ++current_local_id_;
            }

            void TriangleIjkIterator::next ( )
            {
                if ( is_end_ )
                    return;

                if ( current_ijk_id_[0] < last_index_in_row ( current_ijk_id_[1] ) )
                    ++( current_ijk_id_[0] );
                else
                {
                    if ( current_ijk_id_[1] < fe_degree_ )
                    {
                        ++( current_ijk_id_[1] );
                        current_ijk_id_[0] = 0;
                    }
                    else
                    {
                        //We have reached the final dof, and cannot go further.
                        is_end_ = true;
                        return;
                    }
                }

                ++current_local_id_;
            }

            void TetrahedronIjkIterator::next ( )
            {
                if ( is_end_ )
                    return;

                if ( current_ijk_id_[0] < last_index_in_row ( current_ijk_id_[1],
                                                              current_ijk_id_[2] ) )
                    ++( current_ijk_id_[0] );
                else
                {
                    if ( current_ijk_id_[1] < last_index_in_layer ( current_ijk_id_[2] ) )
                    {
                        ++( current_ijk_id_[1] );
                        current_ijk_id_[0] = 0;
                    }
                    else
                    {
                        if ( current_ijk_id_[2] < fe_degree_ )
                        {
                            ++( current_ijk_id_[2] );
                            current_ijk_id_[0] = 0;
                            current_ijk_id_[1] = 0;
                        }
                        else
                        {
                            //We have reached the final dof, and cannot go further.
                            is_end_ = true;
                            return;
                        }
                    }
                }

                ++current_local_id_;
            }

            /// For a given cell and variable, returns an ijk iterator of corresponding type

            template<class LAD>
            IjkIterator* IjkIteratorFactory<LAD>::get ( int var ) const
            {
                CellType::Tag t = element_->get_cell ( ).cell_type ( ).tag ( );

                const FEType<DataType>* fetype = element_->get_fe_type ( var );

                int degree = fetype->get_fe_deg ( );

                switch ( t )
                {
                    case CellType::QUADRILATERAL:
                        return new QuadrilateralIjkIterator ( degree );
                    case CellType::TRIANGLE:
                        return new TriangleIjkIterator ( degree );
                    case CellType::HEXAHEDRON:
                        return new HexahedronIjkIterator ( degree );
                    case CellType::TETRAHEDRON:
                        return new TetrahedronIjkIterator ( degree );
                    default:
                        throw std::invalid_argument ( "invalid cell type in IjkIteratorFactory" );
                }
            }

            template class IjkIteratorFactory<LADescriptorCoupledD>;
            template class IjkIteratorFactory<LADescriptorCoupledS>;

            template<class LAD>
            DofIdConverter<LAD>::DofIdConverter ( const Element<DataType>* element )
            : element_ ( element )
            {
                assert ( element );

                tag_ = element->get_cell ( ).cell_type ( ).tag ( );

                // build global -> local map

                int total_num_dofs = 0;
                for ( int var = 0; var < element_->get_num_variables ( ); ++var )
                    total_num_dofs += element_->get_num_dofs ( var );

                g2lmap_ = boost::unordered_map<int, int>( 1.2 * total_num_dofs );

                init_g2l_map ( );

                //next, build ijk translation data structures

                l2ijk_map_.resize ( element->get_num_variables ( ),
                                    std::vector<ijk_dof_id>( ) );
                ijk2l_map_.resize ( element->get_num_variables ( ) );

                for ( int var = 0; var < element->get_num_variables ( ); ++var )
                {
                    IjkIteratorFactory<LAD> iter_fact ( element );
                    boost::shared_ptr<IjkIterator> ijk_iter ( iter_fact.get ( var ) );

                    int num_dofs = element->get_num_dofs ( var );

                    l2ijk_map_[var] =
                            std::vector<ijk_dof_id>( num_dofs,
                            ijk_dof_id ( ijk_iter->get_dimension ( ) ) );

                    unsigned fe_degree = ijk_iter->get_fe_degree ( );

                    unsigned k_size = 1;
                    if ( ijk_iter->get_dimension ( ) == 3 )
                    {
                        util::multi_array<int>::size_type array_sizes [] = { fe_degree + 1, fe_degree + 1, fe_degree + 1 };

                        ijk2l_map_[var] = IjkToLocalMapType ( array_sizes );
                    }
                    else
                    {
                        util::multi_array<int>::size_type array_sizes [] = { fe_degree + 1, fe_degree + 1 };

                        ijk2l_map_[var] = IjkToLocalMapType ( array_sizes );
                    }

                    // Unused values of the ijk->local map are supposed to be -1,
                    // thus we fill the map with -1 values.
                    // This is important because the ijk2l_map_ will, as a multidimensional
                    // array, always be a cuboid and tri based geometries such
                    // as a tetrahedron will be embedded within the cuboid of the array.
                    // Hence, there are elements that will never be used and must be
                    // marked accordingly
                    std::fill ( ijk2l_map_[var].begin ( ),
                                ijk2l_map_[var].end ( ),
                                -1 );

                    for (; !ijk_iter->is_at_end ( ); ++( *ijk_iter ) )
                    {
                        l2ijk_map_[var][ijk_iter->get_current_local_id ( )]
                                = ijk_iter->get_current_ijk_id ( );

                        set_ijk2l_map_element ( var, ijk_iter->get_current_ijk_id ( ),
                                                ijk_iter->get_current_local_id ( ) );

                    }

                    assert ( l2ijk_map_[var].size ( ) == element->get_num_dofs ( var ) );
                }
            }

            template class DofIdConverter<LADescriptorCoupledD>;
            template class DofIdConverter<LADescriptorCoupledS>;

            template class OnDemandDofIdConverter<LADescriptorCoupledD>;
            template class OnDemandDofIdConverter<LADescriptorCoupledS>;

            template class IsDofKnownLookup<LADescriptorCoupledD>;
            template class IsDofKnownLookup<LADescriptorCoupledS>;

            template class DirectIsDofKnownLookup<LADescriptorCoupledD>;
            template class DirectIsDofKnownLookup<LADescriptorCoupledS>;

            template class OverlayBasedIsDofKnownLookup<LADescriptorCoupledD>;
            template class OverlayBasedIsDofKnownLookup<LADescriptorCoupledS>;

            template<class LAD>
            void NearestKnownDofsFinder<LAD>::find_nearest_known_dofs ( int local_dof_id,
                                                                        int var,
                                                                        std::vector<int>& nearest_dofs ) const
            {
                nearest_dofs.clear ( );

                const Element<DataType>* element = dof_converter_->get_element ( );

                ijk_dof_id search_origin;
                if ( dof_converter_->map_l2ijk ( local_dof_id, var, search_origin ) )
                {
                    if ( element->get_cell ( ).cell_type ( ).tag ( ) == CellType::QUADRILATERAL ||
                         element->get_cell ( ).cell_type ( ).tag ( ) == CellType::HEXAHEDRON )
                        q_proximity_search ( search_origin, var, nearest_dofs );
                    else
                        p_proximity_search ( search_origin, var, nearest_dofs );
                }
            }

            /// Finds the nearest known dof on a subentity
            /// @return The global dof id of the nearest known dof or -1, if no
            /// known dof has been found.
            /// @param tdim The dimension of the subentity
            /// @param sindex The sindex of the subentity (unused if tdim == tdim of cell)
            /// @param global_dof The global dof of which the nearest dofs shall be found.
            /// It must also reside on the supplied subentity.

            template<class LAD>
            int NearestKnownDofsFinder<LAD>::find_nearest_known_dof_on_subent ( TDim tdim,
                                                                                int sindex,
                                                                                int global_dof,
                                                                                int var ) const
            {
                std::vector<Coord> coordinates;
                std::vector<int> dofs;

                const Element<DataType>* element = dof_converter_->get_element ( );

                if ( tdim < element->get_cell ( ).tdim ( ) )
                {
                    element->space ( ).dof ( ).get_coord_on_subentity ( var,
                                                                        element->get_cell_index ( ),
                                                                        tdim,
                                                                        sindex,
                                                                        coordinates );

                    element->space ( ).dof ( ).get_dofs_on_subentity ( var,
                                                                       element->get_cell_index ( ),
                                                                       tdim,
                                                                       sindex,
                                                                       dofs );
                }
                else
                {
                    element->get_dof_indices ( var, dofs );
                    element->space ( ).dof ( ).get_coord_on_cell ( var,
                                                                   element->get_cell_index ( ),
                                                                   coordinates );
                }

                double min_distance = std::numeric_limits<double>::max ( );
                int min_dof = -1;

                assert ( coordinates.size ( ) == dofs.size ( ) );

                Coord search_dof_coord;
                for ( std::size_t i = 0; i < dofs.size ( ) && search_dof_coord.empty ( ); ++i )
                {
                    if ( dofs[i] == global_dof )
                        search_dof_coord = coordinates[i];
                }

                if ( search_dof_coord.empty ( ) )
                    // We have not found the dof we are looking for in the coordinate list
                    return -1;

                for ( std::size_t i = 0; i < coordinates.size ( ); ++i )
                {
                    assert ( coordinates[i].size ( ) == search_dof_coord.size ( ) );

                    if ( dofs[i] != global_dof &&
                         ( *is_dof_known_ )( dofs[i] ) )
                    {
                        double distance = 0.0;

                        // calculate distance
                        for ( std::size_t axis = 0; axis < search_dof_coord.size ( ); ++axis )
                            distance += std::abs ( search_dof_coord[axis] - coordinates[i][axis] );

                        if ( distance < min_distance )
                        {
                            min_distance = distance;
                            min_dof = dofs[i];
                        }
                    }
                }

                return min_dof;
            }

            /// Performs the search for nearest known dofs on Q (quad based) cells
            /// @param origin A ijk dof id representing the origin of the search
            /// @param var The variable to which the ijk dof id \c origin belongs
            /// @param found_dofs A vector that will contain the nearest known dofs
            /// that have been found once a call to this function returns.

            template<class LAD>
            void NearestKnownDofsFinder<LAD>::q_proximity_search ( ijk_dof_id& search_origin,
                                                                   int var,
                                                                   std::vector<int>& nearest_dofs ) const
            {
                assert ( search_origin.size ( ) == 2 || search_origin.size ( ) == 3 );

                if ( search_origin.size ( ) == 2 )
                {
                    // 2d case

                    // on axis neighbors
                    process_delta_list ( QDeltas::deltas2d_on_axis, search_origin, var,
                                         nearest_dofs );

                    if ( nearest_dofs.empty ( ) )
                        // look at diagonals
                        process_delta_list ( QDeltas::deltas2d_diagonals, search_origin, var,
                                             nearest_dofs );
                }
                else
                {
                    // 3d case

                    // First, try the on axis neighbors
                    process_delta_list ( QDeltas::deltas3d_on_axis, search_origin, var,
                                         nearest_dofs );

                    if ( nearest_dofs.empty ( ) )
                    {
                        // Then look at two axis diagonals
                        process_delta_list ( QDeltas::deltas3d_two_axis_diagonals, search_origin,
                                             var, nearest_dofs );

                        if ( nearest_dofs.empty ( ) )
                            // Finally, look at three axis diagonals
                            process_delta_list ( QDeltas::deltas3d_three_axis_diagonals,
                                                 search_origin, var, nearest_dofs );
                    }
                }
            }

            template class NearestKnownDofsFinder<LADescriptorCoupledD>;
            template class NearestKnownDofsFinder<LADescriptorCoupledS>;

            template class ForEachCell<LADescriptorCoupledD>;
            template class ForEachCell<LADescriptorCoupledS>;

            /// Interpolates the unkown dofs on a cell.
            /// This only works if the mesh on which the interpolation takes place is
            /// not the coarsest mesh in the hierarchy. The algorithm used here can only
            /// interpolate vectors that have been transferred from the next coarser grid.

            template<class LAD>
            void LinearCellInterpolator<LAD>::operator() ( const Element<DataType>& element )
            {

                int cell_dimension = element.get_cell ( ).tdim ( );

                const DofPartition<DataType>& dof_partition = element.space ( ).dof ( );

                const DofIdConverter<LAD>& dof_converter = on_demand_dof_converter_.get_converter ( &element );

                OverlayBasedIsDofKnownLookup<LAD> is_dof_known ( dof_ident_, &dof_converter );

                //BestSurroundingDofSubsetCalculator subset_calculator(element);

                for ( int var = 0; var < element.get_num_variables ( ); ++var )
                {
                    // allows us to modify the overlay
                    is_dof_known.resize_overlay_to_var ( element, var );
                    // disable overlay by default
                    is_dof_known.disable_overlay ( );

                    if ( interpolate_from_interpolated_ )
                    {
                        // configure IsDofKnownLookup such that all dofs
                        // are known except the ones that are marked in
                        // the still unknown bitmap
                        for ( int local_dof_id = 0;
                              local_dof_id < element.get_num_dofs ( var );
                              ++local_dof_id )
                        {
                            int global_id = -1;
                            dof_converter.map_l2g ( local_dof_id, var, global_id );

                            if ( dof_partition.is_dof_on_sd ( global_id ) )
                            {
                                if ( !is_local_dof_still_unknown ( local_dof_id, var, element, dof_partition ) )
                                    is_dof_known.mark_dof_as_known ( local_dof_id );
                            }
                            else is_dof_known.mark_dof_as_known ( local_dof_id );
                        }
                    }

                    std::vector<int> dofs;
                    for ( int dim = 0; dim < cell_dimension; ++dim )
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
                            interpolate_entity ( var, element, *ent, dofs,
                                                 dof_converter, is_dof_known );
                        }
                    }
                    // for dim == cell_dimension
                    element.get_dof_indices ( var, dofs );
                    interpolate_entity ( var, element, element.get_cell ( ), dofs,
                                         dof_converter, is_dof_known );
                }
            }

            /// Interpolate an entity
            /// @param var The variable id
            /// @param element The cell on which the entity resides
            /// @param ent The entity that shall be interpolated
            /// @param dofs The dofs that are on this entity
            /// @param dof_converter A DofIdConverter for the cell on which the entity
            /// resides
            /// @param is_dof_known An \c OverlayBasedIsDofKnownLookup that will be used
            /// to query if a dof is unknown, ie needs to be interpolated.

            template<class LAD>
            void LinearCellInterpolator<LAD>::interpolate_entity ( int var,
                                                                   const Element<DataType>& element,
                                                                   const Entity& ent,
                                                                   const std::vector<int>& dofs,
                                                                   const DofIdConverter<LAD>& dof_converter,
                                                                   OverlayBasedIsDofKnownLookup<LAD>& is_dof_known )
            {
                const DofPartition<DataType>& dof_partition = element.space ( ).dof ( );
                unsigned cell_dimension = element.get_cell ( ).tdim ( );

                if ( ent.tdim ( ) < cell_dimension )
                {
                    // check if a neighboring cell has already been
                    // interpolated completely
                    if ( has_entity_interpolated_neighbor_cells ( ent,
                                                                  element ) )
                    {
                        // mark all dofs on this entity as known and return
                        for ( std::size_t i = 0; i < dofs.size ( ); ++i )
                        {
                            int local_id = -1;
                            dof_converter.map_g2l ( dofs[i], local_id );
                            assert ( local_id != -1 );

                            is_dof_known.mark_dof_as_known ( local_id );
                        }
                        return;
                    }
                }

                for ( std::size_t i = 0; i < dofs.size ( ); ++i )
                {
                    int local_id = -1;
                    dof_converter.map_g2l ( dofs[i], local_id );
                    assert ( local_id != -1 );

                    // overlay must be enabled to avoid interpolating the same
                    // dof several times
                    is_dof_known.enable_overlay ( );
                    if ( dof_partition.is_dof_on_sd ( dofs[i] ) &&
                         !is_dof_known ( dofs[i] ) )
                    {
                        is_dof_known.use_overlay ( interpolate_from_interpolated_ );

                        unsigned num_contributions = 0;
                        DataType sum_contributions = 0.0;

                        std::vector<int> found_nearest_dofs;
                        found_nearest_dofs.reserve ( 12 );

                        if ( ent.tdim ( ) < element.get_cell ( ).tdim ( ) )
                        {

                            // for all cells that contain dof
                            for ( IncidentEntityIterator cell = ent.begin_incident ( cell_dimension );
                                  cell != ent.end_incident ( cell_dimension );
                                  ++cell )
                            {
                                Element<DataType> neighboring_cell ( element.space ( ), cell->index ( ) );

                                DofIdConverter<LAD> neighboring_converter ( &neighboring_cell );

                                OverlayBasedIsDofKnownLookup<LAD> neighboring_lookup ( dof_ident_,
                                                                                       &neighboring_converter );

                                neighboring_lookup.use_overlay ( interpolate_from_interpolated_ );

                                OverlayBasedIsDofKnownLookup<LAD>* used_lookup = &neighboring_lookup;
                                if ( neighboring_cell.get_cell_index ( ) == element.get_cell_index ( ) )
                                {
                                    // we are looking at the cell that we are currently interpolating,
                                    // we can just use our old is dof known lookup
                                    used_lookup = &is_dof_known;
                                }
                                else
                                {
                                    // if we want to interpolate from interpolated dofs,
                                    // we need to adjust the is dof known lookup
                                    if ( interpolate_from_interpolated_ )
                                    {
                                        neighboring_lookup.resize_overlay_to_var ( neighboring_cell, var );

                                        // configure IsDofKnownLookup such that all dofs
                                        // are known except the ones that are marked in
                                        // the still unknown bitmap
                                        for ( int local_dof_id = 0;
                                              local_dof_id < neighboring_cell.get_num_dofs ( var );
                                              ++local_dof_id )
                                        {
                                            if ( !is_local_dof_still_unknown ( local_dof_id, var, neighboring_cell, dof_partition ) )
                                                neighboring_lookup.mark_dof_as_known ( local_dof_id );
                                        }

                                    }
                                }

                                NearestKnownDofsFinder<LAD> neighboring_finder ( &neighboring_converter,
                                                                                 used_lookup );

                                std::vector<int> nearest_dofs;

                                int neighbor_local_id = -1;
                                if ( neighboring_converter.map_g2l ( dofs[i], neighbor_local_id ) )
                                {
                                    neighboring_finder.find_nearest_known_dofs ( neighbor_local_id,
                                                                                 var,
                                                                                 nearest_dofs );
                                    for ( std::vector<int>::const_iterator it = nearest_dofs.begin ( );
                                          it != nearest_dofs.end ( ); ++it )
                                    {
                                        // make sure dofs only contribute once
                                        if ( std::find ( found_nearest_dofs.begin ( ), found_nearest_dofs.end ( ), *it )
                                             == found_nearest_dofs.end ( ) )
                                        {
                                            found_nearest_dofs.push_back ( *it );
                                        }
                                    }
                                }
                            }
                        }
                        else // dim == cell tdim
                        {
                            NearestKnownDofsFinder<LAD> nearest_dof_finder ( &dof_converter,
                                                                             &is_dof_known );

                            nearest_dof_finder.find_nearest_known_dofs ( local_id,
                                                                         var,
                                                                         found_nearest_dofs );
                        }

                        //std::vector<int> dofset_to_optimize = found_nearest_dofs;
                        //if(found_nearest_dofs.size() != 12)
                        //    subset_calculator(dofset_to_optimize, dofs[i], found_nearest_dofs);

                        num_contributions = found_nearest_dofs.size ( );

                        // 1 known Dof is not enough to interpolate..
                        if ( num_contributions > 1 )
                        {
                            std::vector<DataType> nearest_values ( num_contributions );

                            vector_->GetValues ( util::raw_array ( found_nearest_dofs ),
                                                 num_contributions,
                                                 util::raw_array ( nearest_values ) );

                            sum_contributions = std::accumulate ( nearest_values.begin ( ),
                                                                  nearest_values.end ( ),
                                                                  0.0 );

                            DataType new_value = sum_contributions / num_contributions;

                            vector_->SetValues ( &dofs[i], 1, &new_value );

                            is_dof_known.mark_dof_as_known ( local_id );
                            mark_dof_in_still_unknown_map ( dof_partition, dofs[i], false );
                        }
                        else
                        {
                            mark_dof_in_still_unknown_map ( dof_partition, dofs[i], true );
                        }

                    }

                }
            }

            template class LinearCellInterpolator<LADescriptorCoupledD>;
            template class LinearCellInterpolator<LADescriptorCoupledS>;

            template class LinearInterpolation<LADescriptorCoupledD>;
            template class LinearInterpolation<LADescriptorCoupledS>;

        } // namespace multigrid
    } // namespace la
} // namespace hiflow
