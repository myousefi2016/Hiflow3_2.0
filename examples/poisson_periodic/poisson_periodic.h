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

/// \author Teresa Beck, Staffan Ronnas, Martin Baumann

// System includes.
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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

// Functor to evaluate the exact solution u of the Poisson problem
// with Dirichlet BC, and its gradient \grad u.

struct ExactSolPeriodX
{

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        const double y = pt[1];
        return -0.5 * ( y - 0.5 )*( y - 0.5 ) + 0.125;
    }

    Vec<DIMENSION, double> eval_grad ( const Vec<DIMENSION, double>& pt ) const
    {

        Vec<DIMENSION, double> grad;
        const double y = pt[1];

        grad[0] = 0.0;
        grad[1] = -0.5 * ( 2. * y - 1. );
        return grad;
    }
};

// Functor to evaluate the exact solution u of the Poisson problem
// with Dirichlet BC, and its gradient \grad u.

struct ExactSolPeriodY
{

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];
        return -0.5 * ( x - 0.5 )*( x - 0.5 ) + 0.125;
    }

    Vec<DIMENSION, double> eval_grad ( const Vec<DIMENSION, double>& pt ) const
    {

        Vec<DIMENSION, double> grad;
        const double x = pt[0];

        grad[0] = -0.5 * ( 2. * x - 1. );
        grad[1] = 0;
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
        return std::vector<double>( coords_on_face.size ( ), 0. );
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
        const int total_dofs = this->num_dofs_total ( );

        lm.Clear ( );
        lm.Resize ( total_dofs, total_dofs );

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
        const int total_dofs = this->num_dofs_total ( );

        lv.clear ( );
        lv.resize ( total_dofs, 0. );

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
        return 1.0;
    }

};

// Functor used for the local evaluation of the square of the L2-norm of the
// error on each element.

class L2ErrorIntegrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    L2ErrorIntegrator ( const CoupledVector<Scalar>& pp_sol, const int periodic_boundaries )
    : pp_sol_ ( pp_sol ), periodic_boundaries_ ( periodic_boundaries )
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
            double sol = 0.;
            switch ( periodic_boundaries_ )
            {
                case 0:
                {
                    // Unique solution but no analytical expression
                    sol = approx_sol_[q];
                    break;
                }
                case 1:
                {
                    sol = sol_period_x_ ( x ( q ) );
                    break;
                }
                case 2:
                {
                    sol = sol_period_y_ ( x ( q ) );
                    break;
                }
                case 3:
                {
                    // Solution not unique!!!
                    sol = 1.0 / 0.0;
                    break;
                }
                default:
                {
                    exit ( -1 );
                    break;
                }
            }
            const double wq = w ( q );
            const double delta = sol - approx_sol_[q];
            value += wq * delta * delta * std::abs ( detJ ( q ) );
        }
    }

  private:
    // coefficients of the computed solution
    const CoupledVector<Scalar>& pp_sol_;
    // Periodicity
    const int periodic_boundaries_;
    // functor to evaluate exact solution
    ExactSolPeriodX sol_period_x_;
    // functor to evaluate exact solution
    ExactSolPeriodY sol_period_y_;
    // vector with values of computed solution evaluated at each quadrature point
    FunctionValues< double > approx_sol_;
};

// Functor used for the local evaluation of the square of the H1-norm of the
// error on each element.

class H1ErrorIntegrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    H1ErrorIntegrator ( const CoupledVector<Scalar>& pp_sol, const int periodic_boundaries )
    : pp_sol_ ( pp_sol ), periodic_boundaries_ ( periodic_boundaries )
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
            Vec<DIMENSION, double> grad_u;
            switch ( periodic_boundaries_ )
            {
                case 0:
                {
                    // Unique solution but no analytical expression
                    grad_u[0] = approx_grad_u_[q][0];
                    grad_u[1] = approx_grad_u_[q][1];
                    break;
                }
                case 1:
                {
                    grad_u = sol_period_x_.eval_grad ( x ( q ) );
                    break;
                }
                case 2:
                {
                    grad_u = sol_period_y_.eval_grad ( x ( q ) );
                    break;
                }
                case 3:
                {
                    // Solution not unique!!!
                    grad_u[0] = 1.0 / 0.0;
                    grad_u[1] = 1.0 / 0.0;
                    break;
                }
                default:
                {
                    exit ( -1 );
                    break;
                }
            }
            value += w ( q ) * ( dot ( grad_u, grad_u )
                    - 2. * dot ( grad_u, approx_grad_u_[q] )
                    + dot ( approx_grad_u_[q], approx_grad_u_[q] ) )
                    * std::abs ( detJ ( q ) );
        }
    }

  private:
    // coefficients of the computed solution
    const CoupledVector<Scalar>& pp_sol_;
    // Periodicity
    const int periodic_boundaries_;
    // functor to evaluate exact solution
    ExactSolPeriodX sol_period_x_;
    // functor to evaluate exact solution
    ExactSolPeriodY sol_period_y_;
    // gradient of computed solution evaluated at each quadrature point
    FunctionValues< Vec<DIMENSION, double> > approx_grad_u_;
};
