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

/// \author Martin Wlotzka

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
typedef LADescriptorCoupledS LAS;

// double precision types
typedef LAD::DataType ScalarD;
typedef LAD::VectorType VectorD;
typedef LAD::MatrixType MatrixD;

// single precision types
typedef LAS::DataType ScalarS;
typedef LAS::VectorType VectorS;
typedef LAS::MatrixType MatrixS;

typedef std::vector<double> CoordD;
typedef std::vector<float> CoordS;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
const int DIMENSION = 2;

// Parameters M, N and O in solution. These decide the period in x-, y- and
// z-direction respectively.
const int M = 2;
const int N = 4;
const int O = 1;

// Functor used to impose u = 0 on the boundary.

struct DirichletZeroD
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<CoordD>& coords_on_face ) const
    {
        // return array with Dirichlet values for dof:s on boundary face
        return std::vector<double>( coords_on_face.size ( ), 0.0 );
    }
};

// Functor used to impose u = 0 on the boundary.

struct DirichletZeroS
{

    std::vector<float> evaluate ( const mesh::Entity& face,
                                  const std::vector<CoordS>& coords_on_face ) const
    {
        // return array with Dirichlet values for dof:s on boundary face
        return std::vector<float>( coords_on_face.size ( ), 0.0 );
    }
};

// Functor used for the local assembly of the stiffness matrix and load vector.

class LocalPoissonAssemblerD : private AssemblyAssistant<DIMENSION, double>
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
        double f_val = static_cast < double > ( -2.0 );
        double tmpx, tmpy, tmpz;

        switch ( DIMENSION )
        {
            case 1:
                break;

            case 2:
                tmpx = pt[0] * ( 1.0 - pt[0] );
                tmpy = pt[1] * ( 1.0 - pt[1] );
                f_val *= ( tmpx + tmpy );
                break;

            case 3:
                tmpx = pt[0] * ( 1.0 - pt[0] );
                tmpy = pt[1] * ( 1.0 - pt[1] );
                tmpz = pt[2] * ( 1.0 - pt[2] );
                f_val *= ( ( tmpx * tmpy ) + ( tmpx * tmpz ) + ( tmpy * tmpz ) );
                break;

            default: assert ( 0 );
        }

        return f_val;
    }
};

// Functor used for the local assembly of the stiffness matrix and load vector.

class LocalPoissonAssemblerF : private AssemblyAssistant<DIMENSION, float>
{
  public:

    void operator() ( const Element<float>& element, const Quadrature<float>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, float>::initialize_for_element ( element, quadrature );

        // compute local matrix
        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const float wq = w ( q );
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

    void operator() ( const Element<float>& element, const Quadrature<float>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, float>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const float wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            {
                lv[dof_index ( i, 0 )] += wq * f ( x ( q ) ) * phi ( i, q ) * std::abs ( detJ ( q ) );
            }
        }
    }

    float f ( Vec<DIMENSION, float> pt )
    {
        float f_val = static_cast < float > ( -2.0 );
        float tmpx, tmpy, tmpz;

        switch ( DIMENSION )
        {
            case 1:
                break;

            case 2:
                tmpx = pt[0] * ( 1.0 - pt[0] );
                tmpy = pt[1] * ( 1.0 - pt[1] );
                f_val *= ( tmpx + tmpy );
                break;

            case 3:
                tmpx = pt[0] * ( 1.0 - pt[0] );
                tmpy = pt[1] * ( 1.0 - pt[1] );
                tmpz = pt[2] * ( 1.0 - pt[2] );
                f_val *= ( ( tmpx * tmpy ) + ( tmpx * tmpz ) + ( tmpy * tmpz ) );
                break;

            default: assert ( 0 );
        }

        return f_val;
    }
};
