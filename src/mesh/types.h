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

#ifndef HIFLOW_MESH_TYPES_H
#    define HIFLOW_MESH_TYPES_H

/// This file contains definitions of simple types used in the Mesh
/// module. These types include integer types for e.g. entity id:s,
/// and imported types from tr1 and boost libraries.

#    include "common/permutation.h"
#    include "common/pointers.h"

#    include <boost/function.hpp>
#    include <boost/noncopyable.hpp>
#    include <boost/intrusive_ptr.hpp>

#    include <iostream>
#    include <vector>

namespace hiflow
{
    namespace mesh
    {
        // Forward declarations
        class Mesh;
        class GeometricSearch;
    }
}

// Support for intrusive reference counting
// NB: These functions are not thread safe.
// These functions are implemented in mesh.cc.
//namespace  {
//    void intrusive_ptr_add_ref(const hiflow::mesh::Mesh* ptr);
//    void intrusive_ptr_release(const hiflow::mesh::Mesh* ptr);
//}

namespace hiflow
{
    namespace mesh
    {
        void intrusive_ptr_add_ref ( const hiflow::mesh::Mesh* ptr );
        void intrusive_ptr_release ( const hiflow::mesh::Mesh* ptr );

        class Entity;
        class RefinedCell;

        /// Types for Mesh pointers
        typedef boost::intrusive_ptr<Mesh> MeshPtr;
        typedef boost::intrusive_ptr<const Mesh> ConstMeshPtr;

        /// Types for Geometric Search pointer
        typedef SharedPtr<GeometricSearch>::Type GeometricSearchPtr;

        /// Type for representing entity Id:s. Id:s are arbitrary, non-negative integer values.
        typedef int Id;

        /// Type for representing entity numbers. A sequence of
        /// EntityNumber:s always starts with 0 and is consecutive.
        typedef int EntityNumber;

        /// Type for representing sizes of entity collections
        typedef int EntityCount;

        /// Type for topological dimensions. This is the maximum dimension
        /// possible for entities in the mesh.
        typedef int TDim;

        /// Type for geometrical dimension. This is the dimension of the
        /// space in which a mesh is embedded.
        typedef int GDim;

        /// Type for the material number.
        typedef int MaterialNumber;

        /// Type for subdomain id
        typedef int SubDomainId;

        /// Type representing a geometrical coordinate.
        typedef double Coordinate;

        /// Iterator for sequences of vertex id:s
        typedef std::vector<Id>::const_iterator VertexIdIterator;

        typedef boost::function<void(const Entity&, std::vector<Coordinate>& ) > RefinedGeometryFunction;

        typedef boost::function<void(const Entity&, const RefinedCell&, std::vector<Coordinate>& ) > RecursiveRefinedGeometryFunction;

        /// Base class for classes without copy ctor and assignment
        /// operator
        typedef boost::noncopyable NonCopyable;

        /// Macros
#    define NOT_YET_IMPLEMENTED                                     \
        { std::clog << "\nError: Function \"" << __FUNCTION__   \
                    << "\" in file \""        << __FILE__       \
                    << "\" at line \""        << __LINE__       \
                    << "\" isn't implemented yet!"              \
                    << std::endl;                               \
            assert(0); }

        template<class T, class I>
        void range_check ( const T& container, I i )
        {
            assert ( i >= 0 );
            assert ( i < static_cast < I > ( container.size ( ) ) );
        }

    } // namespace mesh
} // namespace hiflow

#endif /* _TYPES_H_ */
