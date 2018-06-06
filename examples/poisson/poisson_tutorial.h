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

typedef std::vector<double> Coord;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
const int DIMENSION = 2;

// Choose mesh implementation
#define nUSE_MESH_P4EST

#ifndef WITH_P4EST
#    undef USE_MESH_P4EST
#endif

// Parameters M, N and O in solution. These decide the period in x-, y- and
// z-direction respectively.
const int M = 1;
const int N = 1;
const int O = 1;

// Functor to evaluate the exact solution u of the Poisson problem
// with Dirichlet BC, and its gradient \grad u.

struct ExactSol
{

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];
        const double y = ( DIMENSION > 1 ) ? pt[1] : 0;
        const double z = ( DIMENSION > 2 ) ? pt[2] : 0;
        const double pi = M_PI;
        double solution;

        switch ( DIMENSION )
        {
            case 1:
            {
                solution = 10.0 * std::sin ( 2. * M * pi * x );
                break;
            }
            case 2:
            {
                solution = 10.0 * std::sin ( 2. * M * pi * x ) * std::sin ( 2. * N * pi * y );
                break;
            }
            case 3:
            {
                solution = 10.0 * std::sin ( 2. * M * pi * x ) * std::sin ( 2. * N * pi * y ) * std::sin ( 2. * O * pi * z );
                break;
            }

            default: assert ( 0 );
        }
        return solution;
    }

    Vec<DIMENSION, double> eval_grad ( const Vec<DIMENSION, double>& pt ) const
    {
        Vec<DIMENSION, double> grad;
        const double x = pt[0];
        const double y = ( DIMENSION > 1 ) ? pt[1] : 0;
        const double z = ( DIMENSION > 2 ) ? pt[2] : 0;
        const double pi = M_PI;

        switch ( DIMENSION )
        {
            case 1:
            {
                grad[0] = 20. * M * pi * std::cos ( 2. * M * pi * x );
                break;
            }
            case 2:
            {
                grad[0] = 20. * M * pi * std::cos ( 2. * M * pi * x ) *
                        std::sin ( 2. * N * pi * y );
                grad[1] = 20. * N * pi * std::sin ( 2. * M * pi * x ) *
                        std::cos ( 2. * N * pi * y );
                break;
            }
            case 3:
            {
                grad[0] = 20. * M * pi * std::cos ( 2. * M * pi * x ) *
                        std::sin ( 2. * N * pi * y ) * std::sin ( 2. * O * pi * z );
                grad[1] = 20. * N * pi * std::sin ( 2. * M * pi * x ) *
                        std::cos ( 2. * N * pi * y ) * std::sin ( 2. * O * pi * z );
                grad[2] = 20. * O * pi * std::sin ( 2. * M * pi * x ) *
                        std::sin ( 2. * N * pi * y ) * std::cos ( 2. * O * pi * z );
                break;
            }
            default: assert ( 0 );
        }

        return grad;
    }
};

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

template<class ExactSol>
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
        ExactSol sol;
        double rhs_sol;

        switch ( DIMENSION )
        {
            case 1:
            {
                rhs_sol = 4. * M_PI * M_PI * ( M * M ) * sol ( pt );
                break;
            }
            case 2:
            {
                rhs_sol = 4. * M_PI * M_PI * ( M * M + N * N ) * sol ( pt );
                break;
            }
            case 3:
            {
                rhs_sol = 4. * M_PI * M_PI * ( M * M + N * N + O * O ) * sol ( pt );
                break;
            }
            default: assert ( 0 );
        }

        return rhs_sol;
    }
};

// Functor used for the local evaluation of the square of the L2-norm of the
// error on each element.

template<class ExactSol>
class L2ErrorIntegrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    L2ErrorIntegrator ( CoupledVector<Scalar>& pp_sol )
    : pp_sol_ ( pp_sol )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // Evaluate the computed solution at all quadrature points.
        evaluate_fe_function ( pp_sol_, 0, approx_sol_ );

        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double delta = sol_ ( x ( q ) ) - approx_sol_[q];
            value += wq * delta * delta * std::abs ( detJ ( q ) );
        }
    }

  private:
    // coefficients of the computed solution
    const CoupledVector<Scalar>& pp_sol_;
    // functor to evaluate exact solution
    ExactSol sol_;
    // vector with values of computed solution evaluated at each quadrature point
    FunctionValues< double > approx_sol_;
};

// Functor used for the local evaluation of the square of the H1-norm of the
// error on each element.

template<class ExactSol>
class H1ErrorIntegrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    H1ErrorIntegrator ( CoupledVector<Scalar>& pp_sol )
    : pp_sol_ ( pp_sol )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // Evaluate the gradient of the computed solution at all quadrature points.
        evaluate_fe_function_gradients ( pp_sol_, 0, approx_grad_u_ );

        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const Vec<DIMENSION, double> grad_u = sol_.eval_grad ( x ( q ) );
            value += w ( q ) * ( dot ( grad_u, grad_u )
                    - 2. * dot ( grad_u, approx_grad_u_[q] )
                    + dot ( approx_grad_u_[q], approx_grad_u_[q] ) )
                    * std::abs ( detJ ( q ) );
        }
    }

  private:
    // coefficients of the computed solution
    const CoupledVector<Scalar>& pp_sol_;
    // functor to evaluate exact solution
    ExactSol sol_;
    // gradient of computed solution evaluated at each quadrature point
    FunctionValues< Vec<DIMENSION, double> > approx_grad_u_;
};
