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

/// \author Staffan Ronnas<br>Julian Kraemer

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
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
const int DIMENSION = 3;

/*
Domain:
          (0,10,2)______________________________________________________________(20,10,2)
                 /                                                             /|
                / |                                                           / |
               /                                                             /  |
              /   |                                                         /   |
          (0,10,0) _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ / _ _|(20,10,0)
            /    /|                                                       /    /|
           /                                                             /    / |
          /    /  |                                                     /    /  |
         /____ ___ ____________________________________________________/(20,0,2)|
 (0,0,2)|    /    |                                                   |    /    |
        |                                                             |   /     |
        |  /      |                     AIR                           |  /      |
        |                                                             | /       |
 (0,0,0)|/________|___________________________________________________|/(20,0,0)|
        |                      _______________________________        |         |
        |(0,10,-6)|_ _ _ _ _ _/| _ _ _ _ _ _ _ /|_ _ _  _ _  /| _ _ _ |_ _ _ _ _|(20,10,-6)
        |        /           / |              / |           / |       |        /
        |                   /  |_ _ _ _ _ _ _/_ |_ _ _ _ _ /_ |       |       /
        |      /           /___/____________/__/__________/   /       |      /
        |                  |  /     OIL    |  /    OIL    |  /        |     /
        |    /             | /      +      | /      -     | /         |    /
        |                  |/______________|/_____________|/          |   /
        |  /                                                          |  /
        |                             GROUND                          | /
        |/____________________________________________________________|/
      (0,0,-6)                                                       (20,0,-6)

 */
// Functor used for the local assembly of the stiffness matrix and load vector.

class LocalDirectInverseTutorialAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const double sigma = eval_sigma ( element );

        // Computes the local matrix.
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
                            sigma * wq * dot ( grad_phi ( j, q ), grad_phi ( i, q ) ) * std::abs ( detJ ( q ) );
                }
            }
        }
    }

    void operator() ( const Element<double>& element, int facet_num, const Quadrature<double>& quad, LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_facet ( element, quad, facet_num );

        // Computes contribution to local matrix from boundary facet integral.
        // This terms has its origin in the Robin boundary condition.
        const int num_q = num_quadrature_points ( );

        const int n_dofs = num_dofs ( 0 );
        const double sigma = eval_sigma ( element );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double g = eval_g ( x ( q ) );

            for ( int i = 0; i < n_dofs; ++i )
            {
                for ( int j = 0; j < n_dofs; ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                            wq * sigma * g * phi ( i, q ) * phi ( j, q ) * std::abs ( ds ( q ) );
                }
            }
        }
    }

    // Computes the right hand side.

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

    // Computes f. This is 1 in the left half of the oilfield, -1 in the other half
    // and 0 elsewhere.

    double f ( Vec<DIMENSION, double> pt )
    {
        double rhs_sol;

        rhs_sol = 0;
        if ( pt[1] <= 6 )
        {
            if ( pt[1] >= 4 )
            {
                if ( pt[2] <= -2 )
                {
                    if ( pt[2] >= -4 )
                    {
                        if ( pt[0] >= 6 )
                        {
                            if ( pt[0] <= 10 )
                            {
                                rhs_sol = 1;
                            }
                            else if ( pt[0] <= 14 )
                            {
                                rhs_sol = -1;
                            }
                        }
                    }
                }
            }
        }

        return rhs_sol;
    }

    // Computes Sigma, which is constant in the air and in the ground, but has a
    // different value in the air than in the ground.

    double eval_sigma ( const Element<double>& element )
    {
        const int mat = element.get_cell ( ).get_material_number ( );
        if ( mat == 2 )
        { // The air has the material number 2.
            return 1;
        }
        else if ( mat == 1 )
        { // The ground has the material number 1.
            return 2;
        }
        throw "Unknown material number!";
        return 0.;
    }

    // Computes g. It punishes the distance to the center.

    double eval_g ( Vec<DIMENSION, double> pt )
    {
        std::vector<double> center ( 3, 0.0 );
        center[0] = 10.0;
        center[1] = 5.0;
        center[2] = -3.0;
        const Vec<3, double> center_vec ( center );
        return 1. / ( norm ( pt - center_vec ) );
    }
};

// Functor used for the local assembly of the mass matrix.

class LocalDirectInverseTutorialMassAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lmm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // Computes the local matrix.
        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            {
                for ( int j = 0; j < n_dofs; ++j )
                {
                    lmm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                            wq * phi ( j, q, 0 ) * phi ( i, q, 0 ) * std::abs ( detJ ( q ) );
                }
            }
        }
    }
};
