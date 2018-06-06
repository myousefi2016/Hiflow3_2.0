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

/// \author Staffan Ronnas

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "test.h"
#include "hiflow.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

const TDim tdim = 3;
const GDim gdim = 3;

static const char* datadir = MESH_DATADIR;

// Refinement level of unit cube.
const int REF_LEVEL = 1;

// Polynomial element degree. This should be chosen so that TestFunction can be
// represented exactly (to machine precision).
const int FE_DEGREE = 2;

// Number of points in each direction for evaluation (including end
// points). This should not correspond to dof nodes for the given FE_DEGREE.
const int N = 4;

// Absolute tolerance for error in FE function evaluation.
const double TOL = 1.e-13;

struct TestFunction
{

    double operator() ( const std::vector<double>& pt ) const
    {
        return 2. * pt.at ( 0 ) + pt.at ( 1 ) - 3. * pt.at ( 2 );
    }

};

int main ( int argc, char** argv )
{

    MPI_Init ( &argc, &argv );
    {
        std::string filename = std::string ( datadir ) + std::string ( "unit_cube.inp" );

        // Create a VectorSpace.
        MeshPtr mesh = read_mesh_from_file ( filename, tdim, gdim, 0 );

        for ( int i = 0; i < REF_LEVEL; ++i )
        {
            mesh = mesh->refine ( );
        }

        VectorSpace<double> space;
        space.Init ( FE_DEGREE, *mesh );

        const int ndofs = space.dof ( ).ndofs_global ( );

        // Project TestFunction into discrete space. This function should be
        // represented "exactly".
        TestFunction f;
        std::vector<double> dof_values ( ndofs, 0 );

        for ( EntityIterator it = mesh->begin ( tdim ), end = mesh->end ( tdim );
              it != end; ++it )
        {
            std::vector< std::vector<double> > coords;
            std::vector<int> dofs;

            space.dof ( ).get_coord_on_cell ( 0, it->index ( ), coords );
            space.dof ( ).get_dofs_on_cell ( 0, it->index ( ), dofs );

            for ( int i = 0; i < static_cast < int > ( dofs.size ( ) ); ++i )
            {
                dof_values.at ( dofs.at ( i ) ) = f ( coords.at ( i ) );
            }
        }

        // Evaluate the corresponding FE function on a grid inside each cell, to
        // verify that the functions
        // Element::evaluate_fe_solution() and VectorSpace::get_solution_value()
        // work correctly.

        // First, create a uniform grid on the reference cell.
        std::vector< std::vector<double> > grid;
        std::vector< double > pt ( 3, 0. );

        for ( int i = 0; i < N; ++i )
        {
            for ( int j = 0; j < N; ++j )
            {
                for ( int k = 0; k < N; ++k )
                {
                    pt[0] = i / double(N - 1 );
                    pt[1] = j / double(N - 1 );
                    pt[2] = k / double(N - 1 );
                    grid.push_back ( pt );
                }
            }
        }

        int k = 0;

        for ( EntityIterator it = mesh->begin ( tdim ), end = mesh->end ( tdim );
              it != end; ++it )
        {
            Element<double> elem ( space, it->index ( ) );
            const doffem::CellTransformation<double>* trans = elem.get_cell_transformation ( );

            for ( int i = 0; i < static_cast < int > ( grid.size ( ) ); ++i )
            {
                // Transform grid point from reference cell to physical cell.
                pt[0] = trans->x ( grid.at ( i ) );
                pt[1] = trans->y ( grid.at ( i ) );
                pt[2] = trans->z ( grid.at ( i ) );

                // Compute exact value of solution.
                const double exact_val = f ( pt );

                // Compute FE value of solution.
                const double fe_val = elem.evaluate_fe_solution ( 0, pt, dof_values );

                TEST_EQUAL_EPS ( fe_val, exact_val, TOL );
                ++k;
            }
        }

        std::clog << "Values correspond at all " << k << " points!\n";

    }

    MPI_Finalize ( );

    return 0;
}
