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

/// \author Aksel Alpay, Martin Wlotzka

#ifndef GMG_TUTORIAL_H
#    define GMG_TUTORIAL_H

// System includes.
#    include <cmath>
#    include <fstream>
#    include <iostream>
#    include <string>
#    include <vector>
#    include <mpi.h>
#    include "hiflow.h"

using namespace hiflow::la::gmg;
// Shorten some datatypes with typedefs.
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef std::vector<double> Coord;

// Rank of the master process.
const int MASTER_RANK = 0;

#    define master_here(msg)                   \
{                                       \
    int rank = 0;                       \
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);\
    if(rank == MASTER_RANK)             \
        here(msg);                         \
}

// Dimension of the problem.
const int DIMENSION = 2;

const double M = 2.0;
const double N = 4.0;

// Functor used to impose u = 0 on the boundary.

struct DirichletZero
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        // return array with Dirichlet values for dof:s on boundary face
        return std::vector<double>( coords_on_face.size ( ), 0.0 );
    }
};

// Functor used for the local assembly of the stiffness matrix and load vector.

class LocalPoissonAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // compute local matrix
        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            {
                for ( int j = 0; j < n_dofs; ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                            wq * dot ( grad_phi ( j, q ), grad_phi ( i, q ) ) * std::abs ( detJ ( q ) );
                }
            }
        }
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            {
                lv[dof_index ( i, 0 )] += wq * f ( x ( q ) ) * phi ( i, q ) * std::abs ( detJ ( q ) );
            }
        }
    }

    double f ( Vec<DIMENSION, double> pt )
    {
        double rhs_sol;

        switch ( DIMENSION )
        {
            case 1:
            {
                rhs_sol = 4.;
                break;
            }
            case 2:
            {
                rhs_sol = 4 * M_PI * M_PI * ( M * M + N * N ) * std::sin ( 2 * M_PI * M * pt[0] ) * std::sin ( 2 * M_PI * N * pt[1] );
                break;
            }
            case 3:
            {
                rhs_sol = 4.;
                break;
            }
            default: assert ( 0 );
        }

        return rhs_sol;
    }
};

class MultiVariableLocalPoissonAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // compute local matrix
        for ( int var = 0; var < element.get_num_variables ( ); ++var )
        {
            const int num_q = num_quadrature_points ( );
            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const int n_dofs = num_dofs ( var );
                for ( int i = 0; i < n_dofs; ++i )
                {
                    for ( int j = 0; j < n_dofs; ++j )
                    {
                        lm ( dof_index ( i, var ), dof_index ( j, var ) ) +=
                                wq * dot ( grad_phi ( j, q, var ), grad_phi ( i, q, var ) ) * std::abs ( detJ ( q ) );
                    }
                }
            }
        }
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        for ( int var = 0; var < element.get_num_variables ( ); ++var )
        {
            const int num_q = num_quadrature_points ( );
            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const int n_dofs = num_dofs ( var );
                for ( int i = 0; i < n_dofs; ++i )
                {
                    lv[dof_index ( i, var )] += wq * f ( x ( q ), var ) * phi ( i, q, var ) * std::abs ( detJ ( q ) );
                }
            }
        }
    }

    double f ( Vec<DIMENSION, double> pt, int var )
    {
        double rhs_sol;

        switch ( DIMENSION )
        {
            case 1:
            {
                rhs_sol = 4.;
                break;
            }
            case 2:
            {
                rhs_sol = pt[0] * ( 1.0 - pt[0] ) * pt[1] * ( 1.0 - pt[1] );
                break;
            }
            case 3:
            {
                rhs_sol = 4.;
                break;
            }
            default: assert ( 0 );
        }

        return static_cast < double > ( DIMENSION ) * rhs_sol;
    }
};

#endif
