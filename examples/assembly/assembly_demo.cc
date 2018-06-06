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

#include "hiflow.h"

#include <iostream>
#include <cmath>
#include <vector>

#include <mpi.h>

/// \brief This program demonstrates the use of the
/// AssemblyAssitant for performing
/// local assembly, and performs some simple timing measurements.

/// \author Staffan Ronnas

using namespace std;
using namespace hiflow;

using mesh::TDim;
using mesh::GDim;
using mesh::EntityIterator;
using mesh::MeshPtr;
using mesh::Entity;

static const char* datadir = MESHES_DATADIR;
const TDim tdim = 3;
const GDim gdim = 3;
const int DIM = 3;

void new_assemble ( const Element<double>& element, const Quadrature<double>& quadrature,
                    AssemblyAssistant<DIM, double>::LocalMatrix& lm, AssemblyAssistant<DIM, double>::LocalVector& lv );

int main ( int argc, char *argv[] )
{
    MPI_Init ( &argc, &argv );

    int num_partitions = -1;
    MPI_Comm_size ( MPI_COMM_WORLD, &num_partitions );

    int rank = -1;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    MPI_Comm comm = MPI_COMM_WORLD;

    if ( num_partitions != 1 )
    {
        if ( rank == 0 )
        {
            std::cerr << "assembly_demo only works in sequential mode: exiting!\n\n\n";
        }
        return 1;
    }

    mesh::MeshPtr mesh = read_mesh_from_file ( std::string ( datadir ) + std::string ( "unit_cube.inp" ), tdim, gdim, &comm );

    VectorSpace<double> space;
    std::vector<int> fe_degree ( 4, 2 );
    fe_degree[3] = 1;
    space.Init ( fe_degree, *mesh );

    EntityIterator cell_it = mesh->begin ( tdim );
    Element<double> element ( space, cell_it->index ( ) );
    std::vector<Quadrature<double> > quadrature ( 1, Quadrature<double>( ) );
    quadrature[0].set_quadrature ( "GaussHexahedron", 27 );
    quadrature[0].set_cell_type ( mesh::CellType::HEXAHEDRON );

    AssemblyAssistant<DIM, double>::LocalMatrix lm;
    AssemblyAssistant<DIM, double>::LocalVector lv;

    Timer t_new;
    new_assemble ( element, quadrature[0], lm, lv );
    t_new.stop ( );
    std::cout << "New quadrature took " << t_new << "\n";

    Timer t_new2;
    new_assemble ( element, quadrature[0], lm, lv );
    t_new2.stop ( );
    std::cout << "New quadrature (2) took " << t_new2 << "\n";

    std::vector<int> dof;
    element.get_dof_indices ( dof );

    MPI_Finalize ( );
    return 0;
}

struct ExactSol
{

    double operator() ( const Vec<DIM, double>& pt ) const
    {
        const double x = pt[0];
        const double y = pt[1];
        const double z = pt[2];
        const double pi = M_PI;
        return std::cos ( 2. * pi * x ) * std::cos ( 2. * pi * y ) * std::cos ( 2. * pi * z );
    }
};

class LaplaceLocalAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
    }

    void assemble_local_matrix ( const Element<double>& element, LocalMatrix& lm ) const
    {
        const int num_dofs = num_dofs_total ( );
        const int num_q = num_quadrature_points ( );

        for ( int q = 0; q < num_q; ++q )
        {
            for ( int i = 0; i < num_dofs; ++i )
            {
                for ( int j = 0; j < num_dofs; ++j )
                {
                    lm ( i, j ) += dot ( grad_phi ( i, q ), grad_phi ( j, q ) ) * w ( q ) * std::abs ( detJ ( q ) );
                }
            }
        }
    }

    void assemble_local_vector ( const Element<double>& element, LocalVector& lv ) const
    {
        const int num_dofs = num_dofs_total ( );
        const int num_q = num_quadrature_points ( );

        for ( int q = 0; q < num_q; ++q )
        {
            for ( int i = 0; i < num_dofs; ++i )
            {
                lv[i] += my_f ( x ( q ) ) * phi ( i, q ) * w ( q ) * std::abs ( detJ ( q ) );
            }
        }
    }

  private:

    double my_f ( Vec<DIM, double> pt ) const
    {
        ExactSol sol;
        const double pi = M_PI;
        return 12. * pi * pi * sol ( pt );
    }
};

class StokesAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    StokesAssembler ( )
    {
    }

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
    {
        AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
    }

    void assemble_local_matrix ( const Element<double>& element, LocalMatrix& lm ) const
    {
        const int num_q = num_quadrature_points ( );
        const int total_dofs = num_dofs_total ( );
        lm.Resize ( total_dofs, total_dofs );
        lm.Zeros ( );

        for ( int q = 0; q < num_q; ++q )
        {
            // assemble u0, u1 components (only a(u0, u0), a(u1, u1) and a(u2, v2))
            for ( int u_var = 0; u_var < 3; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                { // test functions
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    { // trial functions
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                dot ( grad_phi ( i, q, u_var ), grad_phi ( j, q, u_var ) ) * w ( q ) * std::abs ( detJ ( q ) );
                    }
                }
            }

            // assemble b(p, v0), b(p, v1), b(p, v2)
            const int p_var = 3;
            for ( int v_var = 0; v_var < 3; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                { // test functions
                    for ( int j = 0; j < num_dofs ( p_var ); ++j )
                    { // trial functions
                        lm ( dof_index ( i, v_var ), dof_index ( j, p_var ) ) +=
                                -phi ( j, q, p_var ) * grad_phi ( i, q, v_var )[v_var] * w ( q ) * std::abs ( detJ ( q ) );
                    }
                }
            }

            // assemble b(u0, q), b(u1, q) and b(u2, q)
            const int q_var = 3;
            for ( int u_var = 0; u_var < 3; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( q_var ); ++i )
                { // test functions
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    { // trial functions
                        lm ( dof_index ( i, q_var ), dof_index ( j, u_var ) ) +=
                                phi ( i, q, q_var ) * grad_phi ( j, q, u_var )[u_var] * w ( q ) * std::abs ( detJ ( q ) );
                    }
                }
            }
        }
    }

    void assemble_local_vector ( const Element<double>& element, LocalVector& lv ) const
    {
        const int num_q = num_quadrature_points ( );
        const int total_dofs = num_dofs_total ( );

        lv.resize ( total_dofs, 0. );

        for ( int q = 0; q < num_q; ++q )
        {
            // assemble l(v0)
            int var = 0;
            for ( int i = 0; i < num_dofs ( var ); ++i )
            {
                lv[dof_index ( i, var )] +=
                        rhs_u0 ( x ( q ) ) * phi ( i, q, var ) * w ( q ) * std::abs ( detJ ( q ) );
            }

            // assemble l(v1)
            var = 1;
            for ( int i = 0; i < num_dofs ( var ); ++i )
            {
                lv[dof_index ( i, var )] +=
                        rhs_u1 ( x ( q ) ) * phi ( i, q, var ) * w ( q ) * std::abs ( detJ ( q ) );
            }

            // assemble l(v2)
            var = 2;
            for ( int i = 0; i < num_dofs ( var ); ++i )
            {
                lv[dof_index ( i, var )] +=
                        rhs_u2 ( x ( q ) ) * phi ( i, q, var ) * w ( q ) * std::abs ( detJ ( q ) );
            }

            // l(p, q) == 0 - no assembly needed
        }
    }

  private:

    double rhs_u0 ( Vec<DIM, double> pt ) const
    {
        return std::sin ( pt[0] * pt[1] );
    }

    double rhs_u1 ( Vec<DIM, double> pt ) const
    {
        return std::cos ( pt[0] * pt[1] );
    }

    double rhs_u2 ( Vec<DIM, double> pt ) const
    {
        return pt[0] + pt[1];
    }

};

void new_assemble ( const Element<double>& element, const Quadrature<double>& quadrature,
                    AssemblyAssistant<DIM, double>::LocalMatrix& lm, AssemblyAssistant<DIM, double>::LocalVector& lv )
{
    //    LaplaceAssembler local_asm;
    StokesAssembler local_asm;
    // assemble locally
    local_asm.initialize_for_element ( element, quadrature );
    Timer matrix_timer;
    local_asm.assemble_local_matrix ( element, lm );
    matrix_timer.stop ( );
    Timer vector_timer;
    local_asm.assemble_local_vector ( element, lv );
    vector_timer.stop ( );
    std::cout << "LocalAssembler matrix assembly took " << matrix_timer << "\n";
    std::cout << "LocalAssembler vector assembly took " << vector_timer << "\n";
}
