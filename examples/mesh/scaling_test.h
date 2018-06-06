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

/// \author Philipp Gerstner

// System includes.
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <mpi.h>
#include "hiflow.h"

// All names are imported for simplicity.
using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

// Shorten some datatypes with typedefs.
#ifdef WITH_HYPRE
typedef LADescriptorHypreD LAD;
#else
typedef LADescriptorCoupledD LAD;
#endif

typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;

typedef std::vector<double> Coord;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
const int DIMENSION = 3;

// Choose mesh implementation
#define USE_MESH_P4EST
#define nPARALLEL_READ

// check compatibility
#ifdef PARALLEL_READ
#    define USE_MESH_P4EST
#endif

#ifndef USE_MESH_P4EST
#    undef PARALLEL_READ
#endif

#ifndef WITH_P4EST
#    undef USE_MESH_P4EST
#    undef PARALLEL_READ
#endif

// Dirichlet boundary condition
// Functor used to impose u = 0 on the boundary.

struct DirichletZero
{

    std::vector<double> evaluate (
                                   const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face
                                   ) const
    {
        // return array with Dirichlet values for dof:s on boundary face
        if ( face.get_material_number ( ) == 11 )
        {
            return std::vector<double>( coords_on_face.size ( ), 0.0 );
        }
        // Neumann boundary values
        return std::vector<double>( 0, 0.0 );
    }
};
// Right hand side f
