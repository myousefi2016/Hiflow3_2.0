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

#include "geometric_search.h"
#include "common/log.h"
#include "mesh_database.h"
#include "mesh_db_view.h"

const int DEBUG_LEVEL = 1;

namespace hiflow
{
    namespace mesh
    {

        GridGeometricSearch::GridGeometricSearch ( MeshPtr mesh )
        : GeometricSearch ( mesh )
        {
            GDim gdim = mesh_->gdim ( );
            assert ( gdim == 2 || gdim == 3 );

            // Initialise container to find cells efficiently
            // 1. Compute bounding box of the mesh
            BBox<Coordinate> bbox ( gdim );
            for ( int vert = 0; vert < mesh_->num_entities ( 0 ); ++vert )
            {
                bbox.add_point ( mesh_->get_coordinates ( 0, vert ) );
            }
            bbox.uniform_extension ( 10 * GEOM_TOL );

            // 2. Determine a suitable number of intervals for the underlying grid
            const int num_edges = mesh_->num_entities ( 1 );
            const int sqrt_num_edges = std::sqrt ( num_edges );
            std::vector<Coordinate> mean_edge_length ( gdim, 0. );
            int count = 0;
            for ( int index = 0; index < num_edges; index += sqrt_num_edges, ++count )
            {
                std::vector<Coordinate> coords = mesh_->get_coordinates ( 1, index );
                assert ( static_cast < int > ( coords.size ( ) ) == 2 * gdim );
                for ( int d = 0; d < gdim; ++d )
                {
                    mean_edge_length[d] += std::abs ( coords[d] - coords[d + gdim] );
                }
            }
            for ( int d = 0; d < gdim; ++d )
            {
                mean_edge_length[d] /= static_cast < Coordinate > ( count );
            }
            std::vector<int> num_intervals ( gdim );
            for ( int d = 0; d < gdim; ++d )
            {
                num_intervals[d] = static_cast < int > ( std::max ( 1., ( bbox.max ( d ) - bbox.min ( d ) ) / mean_edge_length[d] ) );
            }

            // 3. Construct grid
            grid_.reset ( new Grid<Coordinate>( num_intervals, bbox ) );
            assert ( grid_->get_num_points ( ) > 0 );
            assert ( grid_->get_num_cells ( ) > 0 );
            for ( int d = 0; d < gdim; ++d )
            {
                assert ( grid_->delta ( d ) > 0 );
            }

            // 4. Initialisation of the cell and boundary facet index map
            compute_cell_index_map ( );
            compute_facet_index_map ( );
        }

        void GridGeometricSearch::find_cell ( const std::vector<Coordinate>& point, std::vector<int>& cells ) const
        {
            assert ( !point.empty ( ) );
            assert ( static_cast < int > ( point.size ( ) ) == mesh_->gdim ( ) );
            TDim tdim = mesh_->tdim ( );

            cells.clear ( );
            cells.reserve ( 8 );
            // get the grid cell containing the point
            const int grid_index = grid_->cell_with_point ( point );
            // check, whether the grid cell is part of the grid and if so,
            // check, whether this grid cell has any points
            if ( grid_index == -1 || cell_index_map_[grid_index].empty ( ) )
            {
                return;
            }
            // only look for the point in the mesh cells, that were assigned to the grid cell
            for ( std::list<int>::const_iterator ind_it = cell_index_map_[grid_index].begin ( ),
                  e_ind_it = cell_index_map_[grid_index].end ( );
                  ind_it != e_ind_it; ++ind_it )
            {
                if ( point_inside_entity ( point, tdim, mesh_->get_coordinates ( tdim, *ind_it ) ) )
                {
                    cells.push_back ( *ind_it );
                }
            }
            return;
        }

        void GridGeometricSearch::find_cell ( const std::vector<Coordinate>& point, const std::vector<int>& trial_cells, std::vector<int>& cells ) const
        {
            assert ( !point.empty ( ) );
            assert ( static_cast < int > ( point.size ( ) ) == mesh_->gdim ( ) );
            TDim tdim = mesh_->tdim ( );

            cells.clear ( );
            cells.reserve ( 8 );
            // check wether point lies in one of the trial cells
            for ( std::vector<int>::const_iterator ind_it = trial_cells.begin ( ),
                  e_ind_it = trial_cells.end ( );
                  ind_it != e_ind_it; ++ind_it )
            {
                if ( point_inside_entity ( point, tdim, mesh_->get_coordinates ( tdim, *ind_it ) ) )
                {
                    cells.push_back ( *ind_it );
                }
            }

            // if point is not contained in trial cells, search whole grid
            if ( cells.size ( ) == 0 )
            {
                find_cell ( point, cells );
            }
            return;
        }

        bool GridGeometricSearch::point_inside_mesh ( const std::vector<Coordinate>& point ) const
        {
            assert ( !point.empty ( ) );
            assert ( static_cast < int > ( point.size ( ) ) == mesh_->gdim ( ) );
            std::vector<int> cells;
            this->find_cell ( point, cells );
            return (!cells.empty ( ) );
        }

        std::vector<Coordinate> GridGeometricSearch::find_closest_point ( const std::vector<Coordinate>& point,
                                                                          int& facet_index,
                                                                          int material_id ) const
        {
            // TODO: This function only works for 3d with tetrahedrons and in 2d.

            // Finding the closest point by searching through the facet that are
            // near the given point. Those facets are found by searching through
            // the closest grid cells and than using the facet_index_map to get
            // the facet cells that are "close" to the grid cell.
            GDim gdim = mesh_->gdim ( );
            TDim tdim = mesh_->tdim ( );
            assert ( !point.empty ( ) );
            assert ( static_cast < int > ( point.size ( ) ) == gdim );

            // Allocate return variables.
            std::vector<Coordinate> closest_point;
            facet_index = -1;

            // Compute closest point on the grid by a projection on the grid.
            BBox<Coordinate> grid_box = grid_->get_bbox ( );
            std::vector<Coordinate> nearest_grid_point ( point );
            for ( GDim i = 0; i < gdim; i++ )
            {
                if ( point[i] < grid_box.min ( i ) )
                {
                    nearest_grid_point[i] = grid_box.min ( i );
                }
                else if ( point[i] > grid_box.max ( i ) )
                {
                    nearest_grid_point[i] = grid_box.max ( i );
                }
                else
                {
                    nearest_grid_point[i] = point[i];
                }
            }

            // Compute mean grid spacing.
            Coordinate mean_delta = 0;
            for ( GDim d = 0; d < gdim; d++ )
            {
                mean_delta += grid_->delta ( d );
            }
            mean_delta /= ( Coordinate ) gdim;

            // Compute the initial search radius.
            Coordinate initial_radius = distance_point_point ( nearest_grid_point, point ) + mean_delta;

            // The maximal search radius is given by the distance
            // of the point to the outer vertices of the grid.
            Coordinate max_search_radius = 0;
            std::vector<Coordinate> grid_box_vertices = grid_box.get_vertices ( );
            int num_vert = grid_box_vertices.size ( ) / gdim;
            assert ( num_vert == ( int ) std::pow ( ( double ) 2, ( double ) gdim ) );
            for ( GDim n = 0; n < num_vert; n++ )
            {
                std::vector<Coordinate> box_vertex ( grid_box_vertices.begin ( ) + n*gdim, grid_box_vertices.begin ( ) + ( n + 1 ) * gdim );
                max_search_radius = std::max ( max_search_radius, distance_point_point ( point, box_vertex ) );
            }
            max_search_radius += mean_delta;

            // Initialize the distance variable with infinity.
            Coordinate dist = std::numeric_limits<Coordinate>::max ( );

            // Scanned_facets saves all mesh facets that were already scanned
            // and prevents from searching a facet twice.
            SortedArray<int> scanned_cells ( 0 );

            // Search for the closest boundary facet in successively
            // increasing the search radius if not any facet
            // could be found so far.
            for ( Coordinate search_radius = initial_radius;
                  search_radius < max_search_radius;
                  search_radius += mean_delta )
            {
                // The search area is represented by a sphere.
                BSphere<Coordinate> search_sphere ( point, search_radius );
                // Get the indices of the grid cells that intersect the sphere.
                std::vector<int> grid_cell_indices;
                grid_->intersect ( search_sphere, grid_cell_indices );
                // Iterate the grid cell selection.
                for ( std::vector<int>::iterator grid_it = grid_cell_indices.begin ( );
                      grid_it != grid_cell_indices.end ( );
                      ++grid_it )
                {
                    // Iterate the boundary facets assigned to the grid cell.
                    for ( std::list<int>::const_iterator facet_it = facet_index_map_[*grid_it].begin ( ),
                          e_facet_it = facet_index_map_[*grid_it].end ( );
                          facet_it != e_facet_it;
                          ++facet_it )
                    {
                        if ( !scanned_cells.find_insert ( *facet_it ) )
                        {
                            if ( ( material_id != -1 ) && ( boundary_mesh_->get_entity ( tdim - 1, *facet_it ).get_material_number ( ) != material_id ) ) continue;
                            // Calculate the distance between the point and the current facet.
                            // At this, retrieve the closest point on the facet.
                            std::vector<Coordinate> temp_point;
                            std::vector<Coordinate> facet_vertices;
                            boundary_mesh_->get_entity ( tdim - 1, *facet_it ).get_coordinates ( facet_vertices );
                            Coordinate temp_dist = distance_point_facet ( point, facet_vertices, temp_point );
                            // we have a closer point if the distance is smaller
                            if ( temp_dist < dist )
                            {
                                dist = temp_dist;
                                closest_point = temp_point;
                                facet_index = *facet_it;
                            }
                        }
                    }
                }
                // If the determined point does not lie within the search sphere,
                // there could be other boundary facets being closer to the point.
                if ( dist < search_radius )
                {
                    return closest_point;
                }
            }

            // If the grid does not contain any boundary facets, an empty vector is returned.
            return closest_point;
        }

        std::vector<Coordinate> GridGeometricSearch::find_closest_point_parallel ( const std::vector<Coordinate>& point,
                                                                                   int& facet_index,
                                                                                   const MPI_Comm& comm, int material_id ) const
        {
            GDim gdim = mesh_->gdim ( );
            assert ( !point.empty ( ) );
            assert ( static_cast < int > ( point.size ( ) ) == gdim );
            int rank, num_partitions;
            MPI_Comm_rank ( comm, &rank );
            MPI_Comm_size ( comm, &num_partitions );

            facet_index = -1;
            int facet_index_local = -1;
            // has to be double because MPI_DOUBLE is used
            std::vector<double> closest_point = find_closest_point ( point, facet_index_local, material_id );
            // has to be double because MPI_DOUBLE is used
            double local_distance = std::numeric_limits<double>::max ( );
            if ( !closest_point.empty ( ) )
            {
                local_distance = distance_point_point ( point, closest_point );
            }
            else
            {
                closest_point.resize ( gdim );
            }

            std::vector<MPI_Request> req ( num_partitions );
            std::vector<MPI_Status> status ( num_partitions );
            // has to be double because MPI_DOUBLE is used
            std::vector<double> all_distances ( num_partitions );
            // send all local_distances to the process with rank 0
            int rk_with_closest_point = -1;
            if ( rank == 0 )
            {
                all_distances[0] = local_distance;
                for ( int part = 1; part < num_partitions; part++ )
                {
                    MPI_Irecv ( &( all_distances[part] ), 1, MPI_DOUBLE, part, 0, comm, &req[part] );
                }
            }
            else
            {
                MPI_Send ( &local_distance, 1, MPI_DOUBLE, 0, 0, comm );
            }

            if ( rank == 0 )
            {
                double global_distance = local_distance;
                rk_with_closest_point = rank;
                for ( int part = 1; part < num_partitions; part++ )
                {
                    MPI_Wait ( &req[part], &status[part] );
                    if ( all_distances[part] < global_distance )
                    {
                        global_distance = all_distances[part];
                        rk_with_closest_point = part;
                    }
                }
            }
            // at this point rank 0 is the only process who knows which
            // process has the closest point
            // rank 0 Bcasts his knowledge
            MPI_Bcast ( &rk_with_closest_point, 1, MPI_INT, 0, comm );
            assert ( rk_with_closest_point != -1 );
            if ( rank == rk_with_closest_point )
            {
                //assert(facet_index_local != -1);
                facet_index = facet_index_local;
            }
            else
            {
                assert ( facet_index == -1 );
            }
            // now since every point knows who is the owner of the closest
            // point, they send/wait for the data
            assert ( static_cast < int > ( closest_point.size ( ) ) == gdim );
            MPI_Bcast ( &closest_point[0], gdim, MPI_DOUBLE, rk_with_closest_point, comm );
            // since the facet_id is a local index, only the owner of the point
            // has to know it. The other processes will return an index of -1
            // to indicate that they are not the owner of the point.

            return closest_point;
        }

        std::vector<Coordinate> GridGeometricSearch::intersect_boundary ( const std::vector<Coordinate>& point_a,
                                                                          const std::vector<Coordinate>& point_b ) const
        {
            GDim gdim = mesh_->gdim ( );
            TDim tdim = mesh_->tdim ( );
            assert ( !point_a.empty ( ) );
            assert ( !point_b.empty ( ) );
            assert ( static_cast < int > ( point_a.size ( ) ) == gdim );
            assert ( static_cast < int > ( point_b.size ( ) ) == gdim );

            std::vector<Coordinate> intersections;

            // compute bounding box of a and b
            BBox<Coordinate> line_bbox ( gdim );
            line_bbox.add_point ( point_a );
            line_bbox.add_point ( point_b );
            // small enlargement of the box for the case, that the line is parallel to a coordinate axis
            line_bbox.uniform_extension ( GEOM_TOL );

            // get grid cells
            std::vector<int> grid_cell_indices;
            grid_->intersect ( line_bbox, grid_cell_indices );

            // a container to make sure to check a facet only once
            SortedArray<int> scanned_cells;

            // iterate grid cell selection
            for ( std::vector<int>::iterator grid_it = grid_cell_indices.begin ( );
                  grid_it != grid_cell_indices.end ( );
                  ++grid_it )
            {
                // iterate boundary facets
                for ( std::list<int>::const_iterator facet_it = facet_index_map_[*grid_it].begin ( ),
                      e_facet_it = facet_index_map_[*grid_it].end ( );
                      facet_it != e_facet_it;
                      ++facet_it )
                {
                    if ( !scanned_cells.find_insert ( *facet_it ) )
                    {
                        // compute intersection
                        const std::vector<Coordinate> vertices = boundary_mesh_->get_coordinates ( tdim - 1, *facet_it );
                        const std::vector<Coordinate> current_intersection = intersect_facet ( point_a, point_b, vertices );
                        if ( !current_intersection.empty ( ) )
                        {
                            intersections.insert ( intersections.end ( ), current_intersection.begin ( ), current_intersection.end ( ) );
                        }
                    }
                }
            }
            return intersections;
        }

        bool GridGeometricSearch::crossed_boundary ( const std::vector<Coordinate>& point_a,
                                                     const std::vector<Coordinate>& point_b ) const
        {
            GDim gdim = mesh_->gdim ( );
            TDim tdim = mesh_->tdim ( );
            assert ( !point_a.empty ( ) );
            assert ( !point_b.empty ( ) );
            assert ( static_cast < int > ( point_a.size ( ) ) == gdim );
            assert ( static_cast < int > ( point_b.size ( ) ) == gdim );

            bool crossed = false;

            // compute bounding box of a and b
            BBox<Coordinate> line_bbox ( gdim );
            line_bbox.add_point ( point_a );
            line_bbox.add_point ( point_b );
            // small enlargement of the box for the case, that the line is parallel to a coordinate axis
            line_bbox.uniform_extension ( GEOM_TOL );

            // get grid cells
            std::vector<int> grid_cell_indices;
            grid_->intersect ( line_bbox, grid_cell_indices );

            // a container to make sure to check a facet only once
            SortedArray<int> scanned_cells;

            // iterate grid cell selection
            for ( std::vector<int>::iterator grid_it = grid_cell_indices.begin ( );
                  grid_it != grid_cell_indices.end ( );
                  ++grid_it )
            {
                // iterate boundary facets
                for ( std::list<int>::const_iterator facet_it = facet_index_map_[*grid_it].begin ( ),
                      e_facet_it = facet_index_map_[*grid_it].end ( );
                      facet_it != e_facet_it;
                      ++facet_it )
                {
                    if ( !scanned_cells.find_insert ( *facet_it ) )
                    {
                        // compute intersection
                        const std::vector<Coordinate> vertices = boundary_mesh_->get_coordinates ( tdim - 1, *facet_it );
                        if ( crossed_facet ( point_a, point_b, vertices ) )
                        {
                            return true;
                        }
                    }
                }
            }

            return crossed;
        }

        int GridGeometricSearch::max_mesh_cells ( ) const
        {
            int max_mesh_cells = 0;
            for ( std::vector<std::list<int> >::const_iterator ind_it = cell_index_map_.begin ( ), ind_it_end = cell_index_map_.end ( );
                  ind_it != ind_it_end;
                  ++ind_it )
            {
                max_mesh_cells = std::max ( ( int ) ind_it->size ( ), max_mesh_cells );
            }
            return max_mesh_cells;
        }

        double GridGeometricSearch::mean_mesh_cells ( ) const
        {
            double mean_mesh_cells = 0;
            for ( std::vector<std::list<int> >::const_iterator ind_it = cell_index_map_.begin ( ), ind_it_end = cell_index_map_.end ( );
                  ind_it != ind_it_end;
                  ++ind_it )
            {
                mean_mesh_cells += ind_it->size ( );
            }
            assert ( cell_index_map_.size ( ) > 0 );
            return mean_mesh_cells / ( double ) cell_index_map_.size ( );
        }

        void GridGeometricSearch::compute_cell_index_map ( )
        {
            cell_index_map_.clear ( );
            cell_index_map_.resize ( grid_->get_num_cells ( ) );
            TDim tdim = mesh_->tdim ( );
            GDim gdim = mesh_->gdim ( );

            // iterate mesh cells
            for ( mesh::EntityIterator it = mesh_->begin ( tdim ), end_it = mesh_->end ( tdim );
                  it != end_it;
                  ++it )
            {
                // skip ghost cells
                if ( mesh_->has_attribute ( "_remote_index_", tdim ) )
                {
                    int remote_index;
                    mesh_->get_attribute_value ( "_remote_index_", tdim, it->index ( ),
                                                 &remote_index );
                    if ( remote_index != -1 ) continue;
                }
                // create a bounding box of the current mesh cell
                BBox<Coordinate> cell_bbox ( gdim );
                std::vector<Coordinate> coords;
                it->get_coordinates ( coords );
                cell_bbox.add_points ( coords );
                // get list of grid cell indices that intersect the bounding box of the current mesh cell
                std::vector< int> grid_cell_indices;
                grid_->intersect ( cell_bbox, grid_cell_indices );
                // assign the cell indices to the current grid cell
                for ( std::vector<int>::iterator ind_it = grid_cell_indices.begin ( );
                      ind_it != grid_cell_indices.end ( );
                      ++ind_it )
                {
                    cell_index_map_[*ind_it].push_back ( it->index ( ) );
                }
            }
        }

        void GridGeometricSearch::compute_facet_index_map ( )
        {
            facet_index_map_.clear ( );
            facet_index_map_.resize ( grid_->get_num_cells ( ) );
            TDim tdim = mesh_->tdim ( );
            GDim gdim = mesh_->gdim ( );

            // iterate boundary facets
            for ( mesh::EntityIterator it = boundary_mesh_->begin ( tdim - 1 ), end_it = boundary_mesh_->end ( tdim - 1 );
                  it != end_it;
                  ++it )
            {
                // skip ghost facets
                if ( mesh_->has_attribute ( "_remote_index_", tdim ) )
                {
                    int mesh_facet_index;
                    boundary_mesh_->get_attribute_value ( "_mesh_facet_index_", tdim - 1, it->index ( ), &mesh_facet_index );
                    int remote_index;
                    const Entity mesh_facet = mesh_->get_entity ( tdim - 1, mesh_facet_index );
                    assert ( mesh_facet.num_incident_entities ( tdim ) == 1 );
                    const IncidentEntityIterator cell = mesh_facet.begin_incident ( tdim );
                    mesh_->get_attribute_value ( "_remote_index_", tdim, cell->index ( ), &remote_index );
                    if ( remote_index != -1 ) continue;
                }
                // only consider physical boundaries with a material number != -1
                int material_number = boundary_mesh_->get_material_number ( tdim - 1, it->index ( ) );
                if ( material_number == -1 ) continue;

                // create a bounding box of the current boundary facet
                BBox<double> facet_bbox ( gdim );
                std::vector<Coordinate> coords;
                it->get_coordinates ( coords );
                facet_bbox.add_points ( coords );
                // small enlargement of the box for the case, that the box lies on the grid boundary
                facet_bbox.uniform_extension ( GEOM_TOL );
                // get list of grid cell indices that intersect the bounding box of the current boundary facet
                std::vector<int> grid_cell_indices;
                grid_->intersect ( facet_bbox, grid_cell_indices );
                // assign the cell indices to the current grid cell
                for ( std::vector<int>::iterator ind_it = grid_cell_indices.begin ( );
                      ind_it != grid_cell_indices.end ( );
                      ++ind_it )
                {
                    facet_index_map_[*ind_it].push_back ( it->index ( ) );
                }
            }
        }
    }
}
