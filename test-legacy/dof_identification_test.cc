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

/// \author Martin Baumann, Thomas Gengenbach, Michael Schick

#include <string>
#include <vector>
#include <mpi.h>

#include "test.h"
#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

/// DegreeOfFreedom Identification test
///
/// \brief For a list of files it is analysed, whether the number of DoFs
///        after identification is correct.
///
/// For a continuous FE ansatz (linear/quadratic) the number of DoFs after
/// identification of common DoFs are checked by the number of points of meshes.

static const char* datadir = MESH_DATADIR;

int main ( int argc, char *argv[] )
{
    int refinement_levels = 1;

    MPI_Init ( &argc, &argv );

    // Which files should be checked?

    std::vector<std::string> filenames;
    std::vector<TDim> tdims;
    std::vector<GDim> gdims;

    filenames.push_back ( std::string ( datadir ) + std::string ( "two_triangles_2d.inp" ) );
    tdims.push_back ( 2 );
    gdims.push_back ( 2 );

    filenames.push_back ( std::string ( datadir ) + std::string ( "two_quads_2d.inp" ) );
    tdims.push_back ( 2 );
    gdims.push_back ( 2 );

    filenames.push_back ( std::string ( datadir ) + std::string ( "2d_lung_4g.inp" ) );
    tdims.push_back ( 2 );
    gdims.push_back ( 2 );

    filenames.push_back ( std::string ( datadir ) + std::string ( "two_tetras_3d.inp" ) );
    tdims.push_back ( 3 );
    gdims.push_back ( 3 );

    filenames.push_back ( std::string ( datadir ) + std::string ( "two_hexas_3d.inp" ) );
    tdims.push_back ( 3 );
    gdims.push_back ( 3 );

    filenames.push_back ( std::string ( datadir ) + std::string ( "unit_cube_tetras_3d.inp" ) );
    tdims.push_back ( 3 );
    gdims.push_back ( 3 );

    for ( int test_number = 0; test_number < filenames.size ( ); ++test_number )
    {

        std::string filename = filenames.at ( test_number );
        TDim tdim = tdims.at ( test_number );
        GDim gdim = gdims.at ( test_number );

        /////////////////////////////////////
        // mesh

        MeshBuilder * mb ( new MeshDbViewBuilder ( tdim, gdim ) );
        ScopedPtr<Reader>::Type reader ( new UcdReader ( mb ) );

        std::vector<MeshPtr> mesh ( refinement_levels + 1 );
        reader->read ( filename.c_str ( ), mesh[0] );

        for ( int i = 1; i < refinement_levels + 1; ++i )
            mesh[i] = mesh[i - 1]->refine ( );

        for ( int i = 0; i < refinement_levels; ++i )
        {
            // space

            VectorSpace<double> space;

            // tests

            std::cout << "Testing " << filename << " on mesh level: " << i << std::endl;

            int n_dof = 0;

            // -> Linear Ansatz on Mesh_Level 0

            int q1_dofs = mesh[i]->num_entities ( 0 );
            space.Init ( 1, *mesh[i] );
            n_dof = space.dof ( ).get_nb_dofs ( );
            // space.dof().print_interface_modes();
            TEST_EQUAL ( q1_dofs, n_dof );

            // -> Quadratic Ansatz on Mesh_Level 0

            int q2_dofs = mesh[i + 1]->num_entities ( 0 );
            space.Init ( 2, *mesh[i] );
            n_dof = space.dof ( ).get_nb_dofs ( );
            TEST_EQUAL ( q2_dofs, n_dof );

            std::cout << "===============================" << std::endl;

        }
        delete mb;
    }

    MPI_Finalize ( );

    return 0;

}
