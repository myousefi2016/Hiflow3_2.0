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

#include "geometric_tools.h"
#include "../fem/trilinearhexahedrontransformation.h"
#include "../fem/linearpyramidtransformation.h"
#include <numeric>
#include "common/log.h"

const int DEBUG_LEVEL = 1;

namespace hiflow
{
    namespace mesh
    {

        std::vector<Coordinate> Axpy ( const std::vector<Coordinate>& vector_a,
                                       const std::vector<Coordinate>& vector_b,
                                       const Coordinate alpha )
        {
            assert ( !vector_a.empty ( ) );
            assert ( !vector_b.empty ( ) );
            assert ( vector_a.size ( ) == vector_b.size ( ) );
            std::vector<Coordinate> result ( vector_a );
            for ( size_t i = 0; i != result.size ( ); ++i )
            {
                result[i] += alpha * vector_b[i];
            }
            return result;
        }

        std::vector<Coordinate> Scale ( const std::vector<Coordinate>& vector,
                                        const Coordinate value )
        {
            assert ( !vector.empty ( ) );
            std::vector<Coordinate> result ( vector );
            for ( size_t i = 0; i != vector.size ( ); ++i )
            {
                result[i] *= value;
            }
            return result;
        }

        Coordinate Dot ( const std::vector<Coordinate>& vector_a,
                         const std::vector<Coordinate>& vector_b )
        {
            assert ( !vector_a.empty ( ) );
            assert ( !vector_b.empty ( ) );
            assert ( vector_a.size ( ) == vector_b.size ( ) );
            return std::inner_product ( vector_a.begin ( ), vector_a.end ( ), vector_b.begin ( ), 0. );
        }

        std::vector<Coordinate> cross ( const std::vector<Coordinate>& vector_a,
                                        const std::vector<Coordinate>& vector_b )
        {
            assert ( vector_a.size ( ) == 3 );
            assert ( vector_b.size ( ) == 3 );

            std::vector<Coordinate> result ( 3, 0 );
            result[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
            result[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];
            result[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
            return result;
        }

        std::vector<Coordinate> normal ( const std::vector<Coordinate>& directions )
        {
            assert ( !directions.empty ( ) );

            std::vector<Coordinate> result;
            GDim gdim;
            if ( directions.size ( ) == 2 )
            {
                gdim = 2;
                result.resize ( gdim );
                Coordinate dir_norm = Norm2 ( directions );
                assert ( dir_norm > 0 );
                result[0] = directions[1] / dir_norm;
                result[1] = -directions[0] / dir_norm;
            }
            else if ( directions.size ( ) == 6 )
            {
                gdim = 3;
                std::vector<Coordinate> dir_a ( directions.begin ( ), directions.begin ( ) + gdim );
                std::vector<Coordinate> dir_b ( directions.begin ( ) + gdim, directions.begin ( ) + 2 * gdim );
                result = cross ( dir_a, dir_b );
                Coordinate res_norm = Norm2 ( result );
                assert ( res_norm > 0 );
                result = Scale ( result, 1. / res_norm );
            }
            else
            {
                // directions has the wrong size
                assert ( 0 );
            }
            return result;
        }

        Coordinate Norm2 ( const std::vector<Coordinate>& vector )
        {
            assert ( !vector.empty ( ) );
            return std::sqrt ( std::inner_product ( vector.begin ( ), vector.end ( ), vector.begin ( ), 0. ) );
        }

        Coordinate distance_point_point ( const std::vector<Coordinate>& point_a,
                                          const std::vector<Coordinate>& point_b )
        {
            assert ( !point_a.empty ( ) );
            assert ( !point_b.empty ( ) );
            assert ( point_a.size ( ) == point_b.size ( ) );

            Coordinate dist = 0;
            for ( size_t i = 0; i != point_a.size ( ); i++ )
            {
                dist += ( point_a[i] - point_b[i] ) * ( point_a[i] - point_b[i] );
            }
            return std::sqrt ( dist );
        }

        Coordinate distance_point_hyperplane ( const std::vector<Coordinate>& point,
                                               const std::vector<Coordinate>& origin,
                                               const std::vector<Coordinate>& normal )
        {
            assert ( !point.empty ( ) );
            assert ( !origin.empty ( ) );
            assert ( !normal.empty ( ) );
            assert ( origin.size ( ) == point.size ( ) );
            assert ( normal.size ( ) == point.size ( ) );

            return Dot ( point, normal ) - Dot ( origin, normal );
        }

        Coordinate distance_point_line ( const std::vector<Coordinate>& point,
                                         const std::vector<Coordinate>& origin,
                                         const std::vector<Coordinate>& direction )
        {
            assert ( !point.empty ( ) );
            assert ( !origin.empty ( ) );
            assert ( !direction.empty ( ) );
            assert ( origin.size ( ) == point.size ( ) );
            assert ( direction.size ( ) == point.size ( ) );

            std::vector<Coordinate> foot = foot_point_line ( point, origin, direction );
            return distance_point_point ( point, foot );
        }

        std::vector<Coordinate> foot_point_hyperplane ( const std::vector<Coordinate>& point,
                                                        const std::vector<Coordinate>& origin,
                                                        const std::vector<Coordinate>& normal )
        {
            assert ( !point.empty ( ) );
            assert ( !origin.empty ( ) );
            assert ( !normal.empty ( ) );
            assert ( origin.size ( ) == point.size ( ) );
            assert ( normal.size ( ) == point.size ( ) );

            // orthogonal projection
            // f = p - n * (p-o) dot n / (n dot n)
            Coordinate factor = Dot ( Axpy ( point, origin, -1 ), normal ) / Dot ( normal, normal );
            return Axpy ( point, normal, -factor );
        }

        std::vector<Coordinate> foot_point_line ( const std::vector<Coordinate>& point,
                                                  const std::vector<Coordinate>& origin,
                                                  const std::vector<Coordinate>& direction )
        {
            assert ( !point.empty ( ) );
            assert ( !origin.empty ( ) );
            assert ( !direction.empty ( ) );
            assert ( origin.size ( ) == point.size ( ) );
            assert ( direction.size ( ) == point.size ( ) );

            // foot of point on line
            // foot = o + d * ((p-o) dot d) / (d dot d)
            Coordinate factor = Dot ( Axpy ( point, origin, -1 ), direction ) / Dot ( direction, direction );
            return Axpy ( origin, direction, factor );
        }

        bool gauss ( std::vector<Coordinate>& mat,
                     std::vector<Coordinate>& vec )
        {
            GDim gdim = vec.size ( );
            assert ( static_cast < int > ( mat.size ( ) ) == gdim * gdim );

            for ( GDim i = 0; i < gdim; ++i )
            {

                const int i_ind = i * gdim;

                // find pivot row
                int pivot = i;
                for ( GDim p = i; p < gdim; ++p )
                {
                    if ( std::abs ( mat[pivot * gdim + i] ) < std::abs ( mat[p * gdim + i] ) )
                    {
                        pivot = p;
                    }
                }
                // check if system is solvable
                if ( std::abs ( mat[pivot * gdim + i] ) < GEOM_TOL )
                {
                    return false;
                }
                // swap rows
                if ( pivot != i )
                {
                    for ( GDim n = i; n < gdim; ++n )
                    {
                        std::swap ( mat[pivot * gdim + n], mat[i_ind + n] );
                    }
                    std::swap ( vec[pivot], vec[i] );
                }
                // reduce
                const Coordinate pivot_elem = mat[i_ind + i];
                vec[i] /= pivot_elem;
                for ( GDim n = i; n < gdim; ++n )
                {
                    mat[i_ind + n] /= pivot_elem;
                }
                // elimination forwards
                for ( GDim m = i + 1; m < gdim; ++m )
                {
                    const int m_ind = m * gdim;
                    const Coordinate mat_elem = mat[m_ind + i];
                    vec[m] -= vec[i] * mat_elem;
                    for ( GDim n = i; n < gdim; ++n )
                    {
                        mat[m_ind + n] -= mat[i_ind + n] * mat_elem;
                    }
                }
            }

            // elimination backwards
            for ( GDim i = gdim - 1; i > 0; --i )
            {
                const int i_ind = i * gdim;
                for ( GDim m = i - 1; m >= 0; --m )
                {
                    const int m_ind = m * gdim;
                    const Coordinate mat_elem = mat[m_ind + i];
                    vec[m] -= vec[i] * mat_elem;
                    for ( GDim n = i; n < gdim; ++n )
                    {
                        mat[m_ind + n] -= mat[i_ind + n] * mat_elem;
                    }
                }
            }
            return true;
        }

        Coordinate triangle_area ( const std::vector<Coordinate>& vertices )
        {
            assert ( !vertices.empty ( ) );
            assert ( !( vertices.size ( ) % 3 ) );
            const GDim gdim = vertices.size ( ) / 3;
            assert ( gdim > 1 );

            // 0.5 * |A-B|*dist(C,AB)
            std::vector<Coordinate> p_a ( vertices.begin ( ), vertices.begin ( ) + gdim );
            std::vector<Coordinate> p_b ( vertices.begin ( ) + gdim, vertices.begin ( ) + gdim * 2 );
            std::vector<Coordinate> p_c ( vertices.begin ( ) + gdim * 2, vertices.begin ( ) + gdim * 3 );
            std::vector<Coordinate> dir_a_b ( p_a );
            Coordinate dir_a_b_norm = 0;
            for ( int i = 0; i < gdim; ++i )
            {
                dir_a_b[i] -= p_b[i];
                dir_a_b_norm += dir_a_b[i] * dir_a_b[i];
            }
            dir_a_b_norm = std::sqrt ( dir_a_b_norm );

            return 0.5 * dir_a_b_norm * distance_point_line ( p_c, p_a, dir_a_b );
        }

        Coordinate facet_area ( const std::vector<Coordinate>& vertices,
                                const GDim gdim )
        {
            assert ( !vertices.empty ( ) );
            assert ( gdim == 2 || gdim == 3 );
            assert ( !( vertices.size ( ) % gdim ) );
            Coordinate area = 0;
            switch ( gdim )
            {
                case (2 ):
                {
                    // in 2D: Facet = Line
                    assert ( vertices.size ( ) == 4 );
                    std::vector<Coordinate> vert_a ( vertices.begin ( ), vertices.begin ( ) + gdim );
                    std::vector<Coordinate> vert_b ( vertices.begin ( ) + gdim, vertices.begin ( ) + 2 * gdim );
                    area = distance_point_point ( vert_a, vert_b );
                }
                    break;
                case (3 ):
                {
                    // The vertices have to lie in one plane.
                    const int num_vertices = vertices.size ( ) / gdim;
                    assert ( num_vertices < 4 || points_inside_one_hyperplane ( vertices, gdim ) );
                    for ( int i = 1; i < num_vertices - 1; ++i )
                    {
                        std::vector<Coordinate> tri_verts ( vertices.begin ( ), vertices.begin ( ) + gdim );
                        tri_verts.insert ( tri_verts.end ( ), vertices.begin ( ) + i * gdim, vertices.begin ( ) + ( i + 2 ) * gdim );
                        area += triangle_area ( tri_verts );
                    }
                }
                    break;
                default:
                    NOT_YET_IMPLEMENTED;
                    break;
            }
            return area;
        }

        bool in_plane ( const std::vector<Coordinate>& point,
                        const std::vector<Coordinate>& origin,
                        const std::vector<Coordinate>& normal,
                        const Coordinate eps )
        {
            assert ( !point.empty ( ) );
            assert ( !origin.empty ( ) );
            assert ( !normal.empty ( ) );
            assert ( origin.size ( ) == point.size ( ) );
            assert ( normal.size ( ) == point.size ( ) );

            const Coordinate distance = distance_point_hyperplane ( point, origin, normal );
            return (std::abs ( distance ) < eps );
        }

        bool crossed_plane ( const std::vector<Coordinate>& point_a,
                             const std::vector<Coordinate>& point_b,
                             const std::vector<Coordinate>& origin,
                             const std::vector<Coordinate>& normal )
        {
            assert ( !point_a.empty ( ) );
            assert ( !point_b.empty ( ) );
            assert ( !origin.empty ( ) );
            assert ( !normal.empty ( ) );
            assert ( origin.size ( ) == point_a.size ( ) );
            assert ( normal.size ( ) == point_a.size ( ) );
            assert ( point_b.size ( ) == point_a.size ( ) );

            const Coordinate distance_a = distance_point_hyperplane ( point_a, origin, normal );
            const Coordinate distance_b = distance_point_hyperplane ( point_b, origin, normal );
            return (distance_a * distance_b < 0 );
        }

        bool crossed_facet ( const std::vector<Coordinate>& point_a,
                             const std::vector<Coordinate>& point_b,
                             const std::vector<Coordinate>& vertices )
        {
            assert ( !point_a.empty ( ) );
            assert ( !point_b.empty ( ) );
            assert ( !vertices.empty ( ) );
            assert ( point_b.size ( ) == point_a.size ( ) );
            assert ( !( vertices.size ( ) % point_a.size ( ) ) );

            return !intersect_facet ( point_a, point_b, vertices ).empty ( );
        }

        std::vector<Coordinate> intersect_facet ( const std::vector<Coordinate>& point_a,
                                                  const std::vector<Coordinate>& point_b,
                                                  const std::vector<Coordinate>& vertices )
        {
            assert ( !point_a.empty ( ) );
            assert ( !point_b.empty ( ) );
            assert ( !vertices.empty ( ) );
            assert ( point_b.size ( ) == point_a.size ( ) );
            assert ( !( vertices.size ( ) % point_a.size ( ) ) );

            std::vector<Coordinate> intersection;
            const GDim gdim = point_a.size ( );
            assert ( ( gdim == 2 ) || ( gdim == 3 ) );
            if ( static_cast < int > ( vertices.size ( ) ) != gdim * gdim )
            {
                // TODO: implementation for hexahedron facets
                NOT_YET_IMPLEMENTED;
                return intersection;
            }

            // implemented for a facet beeing a line in 2D and a triangle in 3D

            /*
             * 2D line:    3D triangle:    connection from point_a to point_b:
             *   B        C                       F = point_a
             *   *        *                      *
             *   |        |\                    /
             *   |        | \ E                / g
             *   |        |  \                *
             *   *        *---*              G = point_b
             *   A        A    B
             *
             *   E: x = A + x1(B-A) [+ x2(C-A)]  with xi >= 0 and sum xi <= 1
             *   g: x = F + x3(G-F)              with 0 <= x3 <= 1
             *   equating yields GLS:
             *   x1(B-A) [+ x2(C-A)] + x3(F-G) = F-A
             */

            // compute intersection
            std::vector<Coordinate> mat ( gdim * gdim, 0 );
            std::vector<Coordinate> vec ( gdim, 0 );
            for ( int d = 0; d < gdim; ++d )
            {
                for ( int dir = 0; dir < gdim - 1; ++dir )
                {
                    mat[d * gdim + dir] = vertices[( dir + 1 ) * gdim + d] - vertices[d];
                }
                mat[d * gdim + gdim - 1] = point_a[d] - point_b[d];
                vec[d] = point_a[d] - vertices[d];
            }

            // If the system is not solvable, line and facet are parallel
            const bool solved = gauss ( mat, vec );
            if ( !solved )
            {
                // check parallelism
                std::vector<Coordinate> directions ( gdim * ( gdim - 1 ) );
                for ( int d = 0; d < gdim; ++d )
                {
                    for ( int dir = 0; dir < gdim - 1; ++dir )
                    {
                        directions[dir * gdim + d] = vertices[( dir + 1 ) * gdim + d] - vertices[d];
                    }
                }
                std::vector<Coordinate> facet_normal = normal ( directions );
                std::vector<Coordinate> dir_a_b = Axpy ( point_b, point_a, -1 );
                if ( std::abs ( Dot ( dir_a_b, facet_normal ) ) < GEOM_TOL )
                {
                    // vectors are parallel
                    // TODO: check intersection in this case
                    NOT_YET_IMPLEMENTED;
                }
                return intersection;
            }

            // the facet is intersected if
            // 0 <= x3 <= 1
            // xi >= 0 and sum xi <= 1
            if ( vec[gdim - 1] < -GEOM_TOL ||
                 vec[gdim - 1] > 1 + GEOM_TOL )
            {
                return intersection;
            }
            Coordinate sum = 0;
            for ( int d = 0; d < gdim - 1; ++d )
            {
                if ( vec[d] < -GEOM_TOL )
                {
                    return intersection;
                }
                sum += vec[d];
            }
            if ( sum > 1 + GEOM_TOL )
            {
                return intersection;
            }

            // fill intersection coordinate vector
            intersection.resize ( gdim, 0.0 );
            for ( int d = 0; d < gdim; ++d )
            {
                intersection[d] = point_a[d] + vec[gdim - 1] * ( point_b[d] - point_a[d] );
            }

            return intersection;
        }

        Coordinate distance_point_facet ( const std::vector<Coordinate>& point,
                                          const std::vector<Coordinate>& vertices,
                                          std::vector<Coordinate>& closest_point )
        {
            assert ( !point.empty ( ) );
            assert ( !vertices.empty ( ) );
            assert ( !( vertices.size ( ) % point.size ( ) ) );

            // Working for 3d and 2d
            GDim gdim = point.size ( );
            closest_point.clear ( );

            int num_verts = vertices.size ( ) / gdim;

            switch ( gdim )
            {
                case 1:
                {
                    NOT_YET_IMPLEMENTED;
                }
                    break;
                case 2:
                {
                    // 2D: Distance to a finite line.
                    // Check first if the closest point is the orthogonal projection on the facet.
                    assert ( num_verts == 2 );
                    std::vector< Coordinate> origin ( vertices.begin ( ), vertices.begin ( ) + gdim );
                    std::vector< Coordinate> direction ( gdim );
                    for ( int i = 0; i < gdim; i++ )
                    {
                        direction[i] = vertices[i + gdim] - vertices[i];
                    }
                    std::vector<Coordinate> foot_point = foot_point_line ( point, origin, direction );
                    if ( point_inside_entity ( foot_point, 1, vertices ) )
                    {
                        closest_point = foot_point;
                        assert ( static_cast < int > ( closest_point.size ( ) ) == gdim );
                        return distance_point_point ( closest_point, point );
                    }
                }
                    break;
                case 3:
                {
                    // 3D: first we do a orthogonal projection onto the plane spanned
                    // by the facet. The projected point is called foot_point.
                    assert ( num_verts == 3 || num_verts == 4 );
                    std::vector< Coordinate> origin ( vertices.begin ( ), vertices.begin ( ) + gdim );
                    std::vector< Coordinate> directions ( vertices.begin ( ) + gdim, vertices.begin ( ) + 3 * gdim );
                    for ( int i = 0; i < 2 * gdim; ++i )
                    {
                        directions[i] -= vertices[i % gdim];
                    }
                    std::vector< Coordinate> facet_normal = normal ( directions );
                    // TODO: The vertices of the facet currently have to lie in one plane.
                    assert ( num_verts < 4 || points_inside_one_hyperplane ( vertices, gdim ) );
                    std::vector<Coordinate> foot_point = foot_point_hyperplane ( point, origin, facet_normal );
                    // if the foot_point is inside the entity we are done and return
                    // the foot_point as closest_point
                    if ( point_inside_entity ( foot_point, 2, vertices ) )
                    {
                        closest_point = foot_point;
                        assert ( static_cast < int > ( closest_point.size ( ) ) == gdim );
                        return distance_point_point ( closest_point, point );
                    }
                    // else we project our point onto the subentities (lines)
                    Coordinate dist = std::numeric_limits< Coordinate>::max ( );
                    bool found = false;
                    for ( int j = 0; j < num_verts; j++ )
                    {
                        std::vector<Coordinate> line ( 2 * gdim );
                        std::vector<Coordinate> direction ( gdim );
                        for ( int i = 0; i < gdim; i++ )
                        {
                            line[i] = vertices[i + j * gdim];
                            origin[i] = line[i];
                            line[i + gdim] = vertices[i + ( ( j + 1 ) % num_verts ) * gdim];
                            direction[i] = line[i + gdim] - line[i];
                        }
                        foot_point = foot_point_line ( point, origin, direction );
                        // if the projected point is inside the entity we are done!
                        if ( point_inside_entity ( foot_point, 1, line ) )
                        {
                            Coordinate tmp_dist = distance_point_point ( foot_point, point );
                            if ( tmp_dist < dist )
                            {
                                dist = tmp_dist;
                                closest_point = foot_point;
                                found = true;
                            }
                        }
                    }
                    if ( found )
                    {
                        assert ( static_cast < int > ( closest_point.size ( ) ) == gdim );
                        return dist;
                    }
                }
                    break;
            }
            // If the closest point is not an orthogonal projection of the point
            // on the facet (or its edges in 3D), one of the vertices is the
            // closest point
            Coordinate dist = std::numeric_limits< Coordinate>::max ( );
            for ( int n = 0; n < num_verts; n++ )
            {
                std::vector<Coordinate> current_vertex ( vertices.begin ( ) + n * gdim, vertices.begin ( ) + ( n + 1 ) * gdim );
                Coordinate temp_dist = distance_point_point ( current_vertex, point );
                if ( temp_dist < dist )
                {
                    dist = temp_dist;
                    closest_point = current_vertex;
                }
            }
            assert ( static_cast < int > ( closest_point.size ( ) ) == gdim );
            return dist;
        }

        bool point_inside_entity ( const std::vector<Coordinate>& point,
                                   const int tdim,
                                   const std::vector<Coordinate>& vertices )
        {
            // implemented for lines,
            //                 triangles,
            //                 quadrilaterals in one plane,
            //                 tetrahedrons and
            //                 hexahedrons
            // in the dimensions 2D and 3D
            assert ( !point.empty ( ) );
            assert ( tdim == 1 || tdim == 2 || tdim == 3 );
            assert ( !vertices.empty ( ) );
            assert ( !( vertices.size ( ) % point.size ( ) ) );

            GDim gdim = point.size ( );
            assert ( gdim == 2 || gdim == 3 );
            bool inside = false;

            switch ( tdim )
            {
                    // lines, rectangles and triangles can be handled via distance / volume computation:
                    // if vol(ABCD) == vol(PAB) + vol(PBC) + vol(PCD) + vol(PDA) the point lies in the entity
                case 1:
                {
                    assert ( gdim >= 1 );
                    assert ( static_cast < int > ( vertices.size ( ) ) == gdim * 2 );
                    std::vector<Coordinate> p_a ( vertices.begin ( ), vertices.begin ( ) + gdim );
                    std::vector<Coordinate> p_b ( vertices.begin ( ) + gdim, vertices.begin ( ) + gdim * 2 );
                    Coordinate point_to_a = distance_point_point ( point, p_a );
                    Coordinate point_to_b = distance_point_point ( point, p_b );
                    Coordinate a_to_b = distance_point_point ( p_a, p_b );
                    return (point_to_a + point_to_b < a_to_b + GEOM_TOL );
                }
                    break;
                case 2:
                {
                    assert ( gdim >= 2 );

                    Coordinate area_sum = 0;
                    const int num_vertices = vertices.size ( ) / gdim;
                    // Attention: Here we assume that all vertices lie in one plane!
                    assert ( gdim < 3 || num_vertices < 4 || points_inside_one_hyperplane ( vertices, gdim ) );
                    for ( int i = 0; i < num_vertices; ++i )
                    {
                        std::vector<Coordinate> tri_verts ( point );
                        if ( point[0] == vertices[i * gdim] && point[1] == vertices[i * gdim + 1] )
                        {
                            return true;
                        }
                        tri_verts.insert ( tri_verts.end ( ), vertices.begin ( ) + i * gdim, vertices.begin ( ) + ( i + 1 ) * gdim );
                        tri_verts.insert ( tri_verts.end ( ), vertices.begin ( ) + ( ( i + 1 ) % num_vertices ) * gdim, vertices.begin ( ) + ( ( i + 1 ) % num_vertices + 1 ) * gdim );
                        area_sum += triangle_area ( tri_verts );
                    }
                    Coordinate entitity_area = 0;
                    for ( int i = 1; i < num_vertices - 1; ++i )
                    {
                        std::vector<Coordinate> tri_verts ( vertices.begin ( ), vertices.begin ( ) + gdim );
                        tri_verts.insert ( tri_verts.end ( ), vertices.begin ( ) + i * gdim, vertices.begin ( ) + ( i + 2 ) * gdim );
                        entitity_area += triangle_area ( tri_verts );
                    }
                    return (area_sum < entitity_area + GEOM_TOL );
                }
                    break;
                case 3:
                {
                    assert ( gdim == 3 );
                    if ( static_cast < int > ( vertices.size ( ) ) == ( gdim * ( gdim + 1 ) ) )
                    {
                        // Tetrahedron
                        // Parametric equation to check, where the point is:
                        // P = A + x0(B-A) [ + x1(C-A) ( + x2(D-A) ) ]
                        // This algorithm could also handle lines in 1D and triangles in 2D
                        std::vector<Coordinate> mat ( gdim * gdim, 0 );
                        std::vector<Coordinate> vec ( gdim, 0 );
                        for ( int i = 0; i < gdim; ++i )
                        {
                            for ( int j = 0; j < gdim; ++j )
                            {
                                mat[i * gdim + j] = vertices[( j + 1 ) * gdim + i] - vertices[i];
                            }
                            vec[i] = point[i] - vertices[i];
                        }
                        // solve this linear system of equations
                        const bool solved = gauss ( mat, vec );
                        // if the system is not solvable, the cell is degenerated
                        assert ( solved );

                        // the point lies in the cell, if
                        // xi >= 0
                        // sum xi <= 1
                        Coordinate sum = 0;
                        for ( int d = 0; d < gdim; ++d )
                        {
                            if ( vec[d] < -GEOM_TOL )
                            {
                                return false;
                            }
                            sum += vec[d];
                        }
                        return (sum < 1 + GEOM_TOL );
                    }
                    else if ( static_cast < int > ( vertices.size ( ) ) == ( gdim * 8 ) )
                    {
                        // TODO: implementation of a more performant solutions for hexahedrons
                        // For an arbitrary hex, the four vertices of one facet may not lie in one plane.
                        hiflow::doffem::TriLinearHexahedronTransformation<double> ht ( gdim );
                        ht.reinit ( vertices );
                        std::vector<Coordinate> ref_pt;
                        return (ht.contains_physical_point ( point, &ref_pt ) );
                    }
                    else
                    {
                        //Pyramid
                        hiflow::doffem::LinearPyramidTransformation<double> ht ( gdim );
                        ht.reinit ( vertices );
                        std::vector<Coordinate> ref_pt;
                        return (ht.contains_physical_point ( point, &ref_pt ) );
                    }
                }
                    break;
                default:
                    assert ( 0 );
                    return false;
                    break;
            }
            // A return should have been called before.
            assert ( 0 );
            return inside;
        }

        bool points_inside_one_hyperplane ( const std::vector<Coordinate>& points,
                                            const GDim gdim )
        {
            assert ( !points.empty ( ) );
            assert ( !( points.size ( ) % gdim ) );
            int num_points = points.size ( ) / gdim;
            assert ( num_points > gdim );

            std::vector<Coordinate> directions;
            for ( int i = 0; i < gdim - 1; ++i )
            {
                directions.insert ( directions.end ( ), points.begin ( ), points.begin ( ) + gdim );
                for ( GDim d = 0; d < gdim; ++d )
                {
                    directions[i * gdim + d] -= points[( i + 1 ) * gdim + d];
                }
            }
            std::vector<Coordinate> origin ( points.begin ( ), points.begin ( ) + gdim );
            std::vector<Coordinate> plane_normal = normal ( directions );
            for ( int n = gdim; n < num_points; ++n )
            {
                std::vector<Coordinate> test_point ( points.begin ( ) + n * gdim, points.begin ( ) + ( n + 1 ) * gdim );
                if ( !in_plane ( test_point, origin, plane_normal, GEOM_TOL ) )
                {
                    return false;
                }
            }
            // at this point, it is proved, that all points lie in one plane.
            return true;
        }

        std::vector<Coordinate> project_point ( const BoundaryDomainDescriptor &bdd,
                                                const std::vector<Coordinate> &p,
                                                const MaterialNumber mat )
        {
            LOG_DEBUG ( 2, "Starting vector: " << string_from_range ( p.begin ( ), p.end ( ) ) );
            // starting vector is p
            std::vector<Coordinate> xk = p;
            const int DIM = xk.size ( );
            int num_it = 0;

            // step_length is used as a stopping criteria. Initially chosen 1.
            // to not fullfill the criteria
            Coordinate step_length = 1.;
            // if steplength < TOL the algorithm will stop and it is assumed that
            // a solution has been found.
            const Coordinate TOL = 1e-8;
            //maximal number of iteration
            const int max_it = 1e3;

            //leave early if initial point is already a zero of g.
            {
                const Coordinate init_g_xk = bdd.eval_func ( xk, mat );
                const Coordinate ABS_TOL = 1e-8;
                if ( std::abs ( init_g_xk ) < ABS_TOL )
                {
                    LOG_DEBUG ( 2, "Left early because point already on boundary" );
                    return xk;
                }
            }

            //Following does solve the problem:
            //  min_{x \in R^3} 1/2|| p - x ||^2
            //  s.t.  g(x) = 0
            // Where g is given through the BoundaryDomainDescriptor.
            // The algorithm is based on an SQP approach:
            //
            // x_0 starting point
            // for k=0 to ...
            //   min_{d_k \in R^3} 1/2|| p -(x_k + d_k)||^2
            //   s.t.  g(x_k) + d_k*grad(g)(x_k) = 0
            //
            //   x_{k+1} = x_k + d_k
            // endfor
            //
            // The solution of the minimizing problem is done via the Lagrange
            // function: L(d, l) = 1/2|| p -(x_k + d_k)||^2 + l*(g(x_k)+d*grad(g)(x_k))
            // Need to solve grad(L) = 0.
            // This can be done on an algebraic level to get an "exact" solution.

            while ( ( step_length > TOL ) && ( num_it < max_it ) )
            {
                const Coordinate g_xk = bdd.eval_func ( xk, mat );
                const std::vector<Coordinate> grad_xk = bdd.eval_grad ( xk, mat );
                std::vector<Coordinate> pmxk ( DIM ); // p - xk
                for ( int i = 0; i < DIM; ++i )
                    pmxk[i] = p[i] - xk[i];
                const Coordinate grad_g_dot_pmxk = Dot ( grad_xk, pmxk );

                // lambda_k = (grad(g)(x_k)*(p - x_k) + g(x_k))/||grad(g)(x_k)||^2
                const Coordinate lambdak = ( grad_g_dot_pmxk + g_xk ) / Dot ( grad_xk, grad_xk );

                // d_k = p - x_k - lambda_k*grad(x_k)
                std::vector<Coordinate> dk = Axpy ( pmxk, grad_xk, -lambdak );

                //damping?
                Coordinate damping_factor = 1.;
                // x_{k+1} = x_k + d_k
                xk = Axpy ( xk, dk, damping_factor );

                step_length = sqrt ( Norm2 ( dk ) );
                ++num_it;

                //Some high level debug information.
                LOG_DEBUG ( 99, "lambdak: " << lambdak );
                LOG_DEBUG ( 99, "dk: " << string_from_range ( dk.begin ( ), dk.end ( ) ) );
                LOG_DEBUG ( 99, "xk(updated): " << string_from_range ( xk.begin ( ), xk.end ( ) ) );
                LOG_DEBUG ( 99, "g(xk): " << g_xk );
                LOG_DEBUG ( 99, "grad_g(xk): " << string_from_range ( grad_xk.begin ( ), grad_xk.end ( ) ) );
                LOG_DEBUG ( 99, "steplength: " << step_length );
            }
            if ( num_it == max_it )
                LOG_DEBUG ( 2, "Stopped Iteration because max Iteration was reached" );
            LOG_DEBUG ( 2, "Final vector: " << string_from_range ( xk.begin ( ), xk.end ( ) ) );
            LOG_DEBUG ( 2, "Final defect: " << bdd.eval_func ( xk ) );
            LOG_DEBUG ( 2, "Number of Iterations needed: " << num_it );

            return xk;
        }

    }
}
