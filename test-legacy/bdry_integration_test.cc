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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "test.h"
#include "hiflow.h"

using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

// parameters
const TDim tdim = 3;
const GDim gdim = 3;
const int fe_degree = 2;

const int DEBUG_LEVEL = 1;

static const char* datadir = MESH_DATADIR;

class Integrator : private AssemblyAssistant<tdim, double>
{
  public:

    Integrator ( )
    {
    }

    void operator() ( const Element<double>& element, int facet_number, const Quadrature<double>& quadrature,
            std::vector<double>& value )
    {
        AssemblyAssistant<tdim, double>::initialize_for_facet ( element, quadrature, facet_number );
        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            // error
            const double one = 1.0;

            // multiply with weight and transformation factor
            value[facet_number] += w ( q ) * one * ds ( q );
        }

    }
};

void compute_surface_area ( const std::string filename, const int max_refinements, const double correct_surface_area )
{
    // read mesh
    MeshPtr mesh;
    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    if ( rank == 0 )
    {
        mesh = read_mesh_from_file ( filename, tdim, gdim, 0 );
        // refinements
        for ( int refinements = 0; refinements < max_refinements; ++refinements )
        {
            mesh = mesh->refine ( );
        }
        LOG_DEBUG ( 1, "mesh->num_entities(0) == " << mesh->num_entities ( 0 ) );
    }

    NaiveGraphPartitioner naive_partitioner;
    MeshPtr local_mesh = partition_and_distribute ( mesh, 0, MPI_COMM_WORLD, &naive_partitioner );
    assert ( local_mesh != 0 );

    SharedVertexTable shared_verts;
    mesh = compute_ghost_cells ( *local_mesh, MPI_COMM_WORLD, shared_verts );

    // create finite element space
    VectorSpace<double> space;
    space.Init ( fe_degree, *mesh );

    // compute integral over the whole domain
    Integrator integrator;
    StandardGlobalAssembler<double> global_asm;

    double surface_area;
    global_asm.integrate_scalar_boundary ( space, integrator, surface_area );

    // get local surface area from each process
    double global_surface_area;
    MPI_Reduce ( &surface_area, &global_surface_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    // test if result is correct
    if ( rank == 0 )
    {
        LOG_DEBUG ( 1, filename << " has the surface area == " << global_surface_area );
        std::cerr << "Surface area diff for = " << correct_surface_area - global_surface_area << "\n";
        TEST_EQUAL_EPS ( correct_surface_area, global_surface_area, 1e-6 * correct_surface_area );
    }
}

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    std::ofstream info_log ( "bdry_integration_test_info_output.log" );
    LogKeeper::get_log ( "info" ).set_target ( &info_log );
    std::ofstream debug_file ( "bdry_integration_test_output.log" );
    LogKeeper::get_log ( "debug" ).set_target ( &std::cout );

    // test geometries

    // TEST WITH SMALL TETRA GEOMETRY //////////////////////
    std::string filename = std::string ( datadir ) + std::string ( "one_tet.inp" );
    int refinements = 4;
    double correct_surface_area = 1.5 + sqrt ( 2 * 1.5 ) / 2;
    compute_surface_area ( filename, refinements, correct_surface_area );

    // TEST WITH SMALL HEXA GEOMETRY ///////////////////////
    filename = std::string ( datadir ) + std::string ( "unit_cube.inp" );
    refinements = 3;
    correct_surface_area = 6.0;
    compute_surface_area ( filename, refinements, correct_surface_area );

    LogKeeper::get_log ( "debug" ).flush ( );
    MPI_Finalize ( );
    return 0;
}
