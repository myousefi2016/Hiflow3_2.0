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

/// \author Julian Kraemer

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

// Desired Dimension of the problem. (Implemented dimensions: 1, 2, 3)
const int DIMENSION = 2;

// Functors to evaluate the exact solution y, p and u of the distributed
// control poisson problem with Dirichlet BC

class ExactSol_p
{
  public:

    ExactSol_p ( double lambda ) : lambda_ ( lambda )
    {
    }

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];

        double y;
        if ( DIMENSION > 1 )
        {
            y = pt[1];
        }

        double z;
        if ( DIMENSION > 2 )
        {
            z = pt[2];
        }

        const double pi = M_PI;
        double p_solution;

        switch ( DIMENSION )
        {
            case 1:
            {
                p_solution = ( -1 ) * lambda_ * pi * pi * std::sin ( pi * x );
                break;
            }
            case 2:
            {
                p_solution = ( -2 ) * lambda_ * pi * pi * std::sin ( pi * x ) * std::sin ( pi * y );
                break;
            }
            case 3:
            {
                p_solution = ( -3 ) * lambda_ * pi * pi * std::sin ( pi * x ) * std::sin ( pi * y )
                        * std::sin ( pi * z );
                break;
            }
            default: assert ( 0 );
        }

        return p_solution;
    }
  private:
    const double lambda_;
};

struct ExactSol_y
{

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];

        double y;
        if ( DIMENSION > 1 )
        {
            y = pt[1];
        }

        double z;
        if ( DIMENSION > 2 )
        {
            z = pt[2];
        }

        const double pi = M_PI;
        double y_solution;

        switch ( DIMENSION )
        {
            case 1:
            {
                y_solution = std::sin ( pi * x );
                break;
            }
            case 2:
            {
                y_solution = std::sin ( pi * x ) * std::sin ( pi * y );
                break;
            }
            case 3:
            {
                y_solution = std::sin ( pi * x ) * std::sin ( pi * y ) * std::sin ( pi * z );
                break;
            }
            default: assert ( 0 );
        }
        return y_solution;
    }
};

struct ExactSol_u
{

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];

        double y;
        if ( DIMENSION > 1 )
        {
            y = pt[1];
        }

        double z;
        if ( DIMENSION > 2 )
        {
            z = pt[2];
        }

        const double pi = M_PI;
        double u_solution;

        switch ( DIMENSION )
        {
            case 1:
            {
                u_solution = pi * pi * std::sin ( pi * x );
                break;
            }
            case 2:
            {
                u_solution = 2 * pi * pi * std::sin ( pi * x ) * std::sin ( pi * y );
                break;
            }
            case 3:
            {
                u_solution = 3 * pi * pi * std::sin ( pi * x ) * std::sin ( pi * y )
                        * std::sin ( pi * z );
                break;
            }
            default: assert ( 0 );
        }

        return u_solution;
    }
};

// Functor used to impose y = 0 and p = 0 on the boundary.

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

// The following system of equations is assembled:
//         |  0  A  -M  | |p|    | 0 |
//         |  A -M   0  | |y| =  |-Mŷ|
//         | -M   0 -LM | |u|    | 0 |
//
// with p being the adjoint state, y being the state, u being the optimal
// solution, ŷ being the desired state, L being the regularisation parameter
// lambda, A being the (n x n)-stiffness matrix and M being the
// (n x n)-mass matrix.

class LocalDistributedControlPoissonAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    LocalDistributedControlPoissonAssembler ( double lambda ) : lambda_ ( lambda )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // compute local matrix
        lm.Resize ( num_dofs_total ( ), num_dofs_total ( ) );
        lm.Zeros ( );

        const int num_q = num_quadrature_points ( );
        const int n_dofs_trial_p = num_dofs ( 0 );
        const int n_dofs_trial_y = num_dofs ( 1 );
        const int n_dofs_trial_u = num_dofs ( 2 );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double det_j = std::abs ( detJ ( q ) );

            for ( int i = 0; i < n_dofs_trial_p; ++i )
            {
                for ( int j = 0; j < n_dofs_trial_y; ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 1 ) ) +=
                            wq * dot ( grad_phi ( i, q, 0 ), grad_phi ( j, q, 1 ) ) * det_j;
                }
                for ( int j = 0; j < n_dofs_trial_u; ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 2 ) ) +=
                            ( -1 ) * wq * phi ( i, q, 0 ) * phi ( j, q, 2 ) * det_j;
                }
            }
            for ( int i = 0; i < n_dofs_trial_y; ++i )
            {
                for ( int j = 0; j < n_dofs_trial_p; ++j )
                {
                    lm ( dof_index ( i, 1 ), dof_index ( j, 0 ) ) +=
                            wq * dot ( grad_phi ( i, q, 1 ), grad_phi ( j, q, 0 ) ) * det_j;
                }
                for ( int j = 0; j < n_dofs_trial_y; ++j )
                {
                    lm ( dof_index ( i, 1 ), dof_index ( j, 1 ) ) +=
                            ( -1 ) * wq * phi ( i, q, 1 ) * phi ( j, q, 1 ) * det_j;
                }
            }
            for ( int i = 0; i < n_dofs_trial_u; ++i )
            {
                for ( int j = 0; j < n_dofs_trial_p; ++j )
                {
                    lm ( dof_index ( i, 2 ), dof_index ( j, 0 ) ) +=
                            ( -1 ) * wq * phi ( i, q, 2 ) * phi ( j, q, 0 ) * det_j;
                }
                for ( int j = 0; j < n_dofs_trial_u; ++j )
                {
                    lm ( dof_index ( i, 2 ), dof_index ( j, 2 ) ) +=
                            ( -1 ) * lambda_ * wq * phi ( i, q, 2 ) * phi ( j, q, 2 ) * det_j;
                }
            }
        }
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        // Sets the vector to 0.
        for ( int i = 0; i < static_cast < int > ( lv.size ( ) ); ++i )
        {
            lv[i] = 0;
        }

        const int num_q = num_quadrature_points ( );
        const int n_dofs_trial_y = num_dofs ( 1 );
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            for ( int i = 0; i < n_dofs_trial_y; ++i )
            {
                lv[dof_index ( i, 1 )] += ( -1 ) * wq * y_omega ( x ( q ) )
                        * phi ( i, q, 1 ) * std::abs ( detJ ( q ) );
            }
        }
    }

    // computes y_omega outside the loop.

    double y_omega ( Vec<DIMENSION, double> pt )
    {
        double x = pt[0];
        double y;
        if ( DIMENSION > 1 )
        {
            y = pt[1];
        }

        double z;
        if ( DIMENSION > 2 )
        {
            z = pt[2];
        }

        const double pi = M_PI;
        double y_omega_solution;
        switch ( DIMENSION )
        {
            case 1:
            {
                y_omega_solution = ( 1 + lambda_ * pi * pi * pi * pi ) * std::sin ( pi * x );
                break;
            }
            case 2:
            {
                y_omega_solution = ( 1 + 4 * lambda_ * pi * pi * pi * pi )
                        * std::sin ( pi * x ) * std::sin ( pi * y );
                break;
            }
            case 3:
            {
                y_omega_solution = ( 1 + 9 * lambda_ * pi * pi * pi * pi )
                        * std::sin ( pi * x ) * std::sin ( pi * y ) * std::sin ( pi * z );
                break;
            }
            default: assert ( 0 );
        }
        return y_omega_solution;
    }

  private:
    const double lambda_;
};

// Functor used for the local evaluation of the square of the L2-norm of the
// error on each element.

template<class ExactSol>
class L2ErrorIntegrator_p : private AssemblyAssistant<DIMENSION, double>
{
  public:

    L2ErrorIntegrator_p ( const CoupledVector<Scalar>& pp_sol, ExactSol& sol )
    : pp_sol_ ( pp_sol ), sol_ ( sol )
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
    // functor to evaluate exact solutions
    ExactSol sol_;

    // vector with values of computed solution evaluated at each quadrature point
    FunctionValues< double > approx_sol_;
};

template<class ExactSol>
class L2ErrorIntegrator_y : private AssemblyAssistant<DIMENSION, double>
{
  public:

    L2ErrorIntegrator_y ( const CoupledVector<Scalar>& pp_sol, ExactSol& sol )
    : pp_sol_ ( pp_sol ), sol_ ( sol )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // Evaluate the computed solution at all quadrature points.
        evaluate_fe_function ( pp_sol_, 1, approx_sol_ );

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
    // functor to evaluate exact solutions
    ExactSol sol_;

    // vector with values of computed solution evaluated at each quadrature point
    FunctionValues< double > approx_sol_;
};

template<class ExactSol>
class L2ErrorIntegrator_u : private AssemblyAssistant<DIMENSION, double>
{
  public:

    L2ErrorIntegrator_u ( const CoupledVector<Scalar>& pp_sol, ExactSol& sol )
    : pp_sol_ ( pp_sol ), sol_ ( sol )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // Evaluate the computed solution at all quadrature points.
        evaluate_fe_function ( pp_sol_, 2, approx_sol_ );

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
    // functor to evaluate exact solutions
    ExactSol sol_;

    // vector with values of computed solution evaluated at each quadrature point
    FunctionValues< double > approx_sol_;
};
