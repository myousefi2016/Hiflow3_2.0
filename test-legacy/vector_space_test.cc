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

/// \author Martin Baumann

#include <cmath>
#include <string>
#include <vector>
#include <mpi.h>

#include "test.h"
#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

/// Vector Space test
///
/// \brief For a geometry of two neighbouring hexahedra, Q1 and Q2 FE ansatzs
///        are initialized, also on refined mesh, and correct number of dofs
///        is checked
///

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESH_DATADIR;

typedef std::vector<double> Coord;

int main ( int argc, char *argv[] )
{

    MPI_Init ( &argc, &argv );

    // mesh

    const string filename = string ( datadir ) + string ( "two_hexas_3d.inp" );

    MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
    ScopedPtr<Reader>::Type reader ( new UcdReader ( mb ) );

    MeshPtr mesh;
    reader->read ( filename.c_str ( ), mesh );

    MeshPtr refined_mesh;
    refined_mesh = mesh->refine ( );

    // space

    VectorSpace<double> space;

    // 4 tests

    bool success = true;
    bool verbosity = false;
    int n_dof = 0;

    // -> 1

    space.Init ( 1, *mesh );
    n_dof = space.dof ( ).get_nb_dofs ( );
    if ( n_dof != 12 )
        success = false;
    if ( verbosity == true )
        std::cout << "Level 0, Lagrange Q1, #DoFs: " << n_dof << std::endl;

    // -> 2

    space.Init ( 2, *mesh );
    n_dof = space.dof ( ).get_nb_dofs ( );
    if ( n_dof != 45 )
        success = false;
    if ( verbosity == true )
        std::cout << "Level 0, Lagrange Q2, #DoFs: " << n_dof << std::endl;

    // -> 3

    space.Init ( 1, *refined_mesh );
    n_dof = space.dof ( ).get_nb_dofs ( );
    if ( n_dof != 45 )
        success = false;
    if ( verbosity == true )
        std::cout << "Level 1, Lagrange Q1, #DoFs: " << n_dof << std::endl;

    // -> 4

    space.Init ( 2, *refined_mesh );
    n_dof = space.dof ( ).get_nb_dofs ( );
    if ( n_dof != 225 )
        success = false;
    if ( verbosity == true )
        std::cout << "Level 1, Lagrange Q2, #DoFs: " << n_dof << std::endl;

    // Check

    if ( success == true )
        std::cout << "TEST PASSED SUCCESSFULLY." << std::endl;

    if ( verbosity == true )
        std::cout << "TEST FAILED." << std::endl;

    delete mb;

    MPI_Finalize ( );

    return 0;
}
