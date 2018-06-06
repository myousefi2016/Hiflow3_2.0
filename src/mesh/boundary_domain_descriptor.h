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

/// \author Jonathan Schwegler

#ifndef HIFLOW_MESH_BOUNDARY_DESCRIPTOR
#    define HIFLOW_MESH_BOUNDARY_DESCRIPTOR

#    include "mesh/types.h"

namespace hiflow
{
    namespace mesh
    {

        /// \brief Abstract class to describe the boundary of a mesh
        /// via the zero set of a function. The function itself may depend
        /// on the MaterialNumber to project different parts of the meshs boundary
        /// onto different domains.
        /// One has to describe both, the function as well as its gradient.

        class BoundaryDomainDescriptor
        {
          public:

            virtual Coordinate eval_func ( const std::vector< Coordinate >& x, MaterialNumber mat_num = -1 ) const = 0;

            virtual std::vector<Coordinate> eval_grad ( const std::vector< Coordinate >& x, MaterialNumber mat_num = -1 ) const = 0;
        };

    } //namespace mesh
} //namespace hiflow
#endif
