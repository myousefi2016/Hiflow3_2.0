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

/// \author Jonathan Schwegler and Jonas Kratzke

#ifndef HIFLOW_GEOMETRIC_TOOLS_H_
#    define HIFLOW_GEOMETRIC_TOOLS_H_

#    include <vector>
#    include <cmath>
#    include "mesh/types.h"
#    include "mesh/boundary_domain_descriptor.h"

namespace
{
    // Tolerance for geometrical comparisons.
    const double GEOM_TOL = 1.e-14;
}

namespace hiflow
{
    namespace mesh
    {
        /// \brief vector_a + alpha * vector_b
        std::vector<Coordinate> Axpy ( const std::vector<Coordinate>& vector_a,
                                       const std::vector<Coordinate>& vector_b,
                                       const Coordinate alpha );

        /// \brief value*vector
        std::vector<Coordinate> Scale ( const std::vector<Coordinate>& vector,
                                        const Coordinate value );

        /// \brief < vector_a , vector_b >
        Coordinate Dot ( const std::vector<Coordinate>& vector_a,
                         const std::vector<Coordinate>& vector_b );

        /// \brief Euclidian norm
        Coordinate Norm2 ( const std::vector<Coordinate>& vector );

        /// \brief 3d cross product: vector_a x vector_b.
        std::vector<Coordinate> cross ( const std::vector<Coordinate>& vector_a,
                                        const std::vector<Coordinate>& vector_b );

        /// \brief normal of a hyperplane given by directional vectors
        std::vector<Coordinate> normal ( const std::vector<Coordinate>& directions );

        /// \brief (Euclidian) Distance between two points
        Coordinate distance_point_point ( const std::vector<Coordinate>& point_a,
                                          const std::vector<Coordinate>& point_b );

        /// \brief Distance between point and hyperplane with orientation
        /// It should hold: ||normal|| = 1
        Coordinate distance_point_hyperplane ( const std::vector<Coordinate>& point,
                                               const std::vector<Coordinate>& origin,
                                               const std::vector<Coordinate>& normal );

        /// \brief Distance between point and line (without orientation)
        Coordinate distance_point_line ( const std::vector<Coordinate>& point,
                                         const std::vector<Coordinate>& origin,
                                         const std::vector<Coordinate>& direction );

        /// \brief returns the foot of a point on a hyperplane
        std::vector<Coordinate> foot_point_hyperplane ( const std::vector<Coordinate>& point,
                                                        const std::vector<Coordinate>& origin,
                                                        const std::vector<Coordinate>& normal );

        /// \brief returns the foot of a point on a line
        std::vector<Coordinate> foot_point_line ( const std::vector<Coordinate>& point,
                                                  const std::vector<Coordinate>& origin,
                                                  const std::vector<Coordinate>& direction );

        /// \brief Gaussian elimination with pivoting of a linear system of equations
        /// transforms the given Matrix and vector
        /// Matrix needs to be square and stored row wise -> [a b; c d] =  [a b c d].
        bool gauss ( std::vector<Coordinate>& mat,
                     std::vector<Coordinate>& vec );

        /// \brief returns the area of a triangle
        Coordinate triangle_area ( const std::vector<Coordinate>& vertices );

        /// \brief returns the area of a facet
        /// the area can be a convex simplex with sorted vertices
        /// TODO: For quadrilaterals, the vertices currently have to lie in one plane.
        Coordinate facet_area ( const std::vector<Coordinate>& vertices,
                                const GDim gdim );

        /// \brief checks if a point lies in a hyperplane
        /// with an precision of eps
        bool in_plane ( const std::vector<Coordinate>& point_a,
                        const std::vector<Coordinate>& origin,
                        const std::vector<Coordinate>& normal,
                        const Coordinate eps );

        /// \brief checks if the connection of two points crosses a plane
        bool crossed_plane ( const std::vector<Coordinate>& point_a,
                             const std::vector<Coordinate>& point_b,
                             const std::vector<Coordinate>& origin,
                             const std::vector<Coordinate>& normal );

        /// \brief checks if a facet is crossed by a line from a to b
        /// The input facet is given by the coordinates of its vertices
        bool crossed_facet ( const std::vector<Coordinate>& point_a,
                             const std::vector<Coordinate>& point_b,
                             const std::vector<Coordinate>& vertices );

        /// \brief Computes the intersection of a line with a facet.
        /// The input facet is given by the coordinates of its vertices
        // TODO: implementation for hexahedron facets
        std::vector<Coordinate> intersect_facet ( const std::vector<Coordinate>& point_a,
                                                  const std::vector<Coordinate>& point_b,
                                                  const std::vector<Coordinate>& vertices );

        /// \brief Calculates the distance from a point to a facet where the
        /// facet is of dimension DIM - 1. It also returns the closest point.
        /// The facet is given by it's vertices and the ordering of them is
        /// crucial for quadrilateral facets.
        /// TODO: For quadrilaterals, the vertices currently have to lie in one plane.
        Coordinate distance_point_facet ( const std::vector<Coordinate>& point,
                                          const std::vector<Coordinate>& vertices,
                                          std::vector<Coordinate>& closest_point );

        /// \brief Determines, whether a point lies in an entity of topologic dimension tdim.
        /// The entity is given by the coordinates of its vertices.
        /// TODO: For quadrilaterals, the vertices currently have to lie in one plane.
        bool point_inside_entity ( const std::vector<Coordinate>& point,
                                   const TDim tdim,
                                   const std::vector<Coordinate>& vertices );

        /// \brief Determines, whether a set of point lies in one and the same
        /// hyperplane and if true returns the normal of that hyperplane.
        bool points_inside_one_hyperplane ( const std::vector<Coordinate>& points,
                                            const GDim gdim );

        /// \brief Returns the point that is closest to p but still on the
        /// domain given by the BoundaryDomainDescriptor
        std::vector<Coordinate> project_point ( const BoundaryDomainDescriptor &bdd,
                                                const std::vector<Coordinate> &p,
                                                const MaterialNumber mat );

    }
}

#endif
