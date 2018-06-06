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

/// \author Katrin Mang<br>Simon Gawlok

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
#ifdef WITH_HYPRE
typedef LADescriptorHypreD LAD;
#else
typedef LADescriptorCoupledD LAD;
#endif

typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;

typedef std::vector<double> Coord;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
const int DIMENSION = 2;

// Choose mesh implementation
#define USE_MESH_P4EST

// Set to true if boundary should be adapted to ellipsoid
#define nELLIPSOID_BOUNDARY

// If set, rhs f=const=1. Otherwise, f is defined as below. In this case (#ifndef CONST_F) comparison with the exact solution is possible for
// unit_squre and unit_cubbe
#define CONST_F

#ifndef WITH_P4EST
#    undef USE_MESH_P4EST
#endif

struct QuadratureSelection
{
    /// Constructor
    /// \param[in] order Desired order of quadrature rule

    QuadratureSelection ( const int order ) : order_ ( order )
    {
    }

    /// Operator to obtain quadrature rule on given element
    /// \param[in] elem Element on which quadrature should be done
    /// \param[out] quadrature Quadrature rule on given Element elem

    void operator() ( const Element<double>& elem,
            Quadrature<double>& quadrature )
    {
        // Get ID of FE type
        const FEType<double>::FiniteElement fe_id =
                elem.get_fe_type ( 0 )->get_my_id ( );

        // Switch by FE type for quadrature selection
        switch ( fe_id )
        {
            case FEType<double>::LAGRANGE_TRI:
            {
                quadrature.set_cell_type ( 2 );
                quadrature.set_quadrature_by_order
                        ( "GaussTriangle", order_ );
                break;
            }
            case FEType<double>::LAGRANGE_QUAD:
            {
                quadrature.set_cell_type ( 3 );
                quadrature.set_quadrature_by_order
                        ( "EconomicalGaussQuadrilateral", order_ );
                break;
            }
            case FEType<double>::LAGRANGE_TET:
            {
                quadrature.set_cell_type ( 4 );
                quadrature.set_quadrature_by_order
                        ( "EconomicalGaussTetrahedron", order_ );
                break;
            }
            case FEType<double>::LAGRANGE_HEX:
            {
                quadrature.set_cell_type ( 5 );
                quadrature.set_quadrature_by_order
                        ( "GaussHexahedron", order_ );
                break;
            }
            default:
            {
                assert ( false );
            }
        };
    }
    // Order of quadrature
    const int order_;
};

class Ellipsoid : public BoundaryDomainDescriptor
{
  public:

    //Constructor for a circle or a sphere

    Ellipsoid ( const Coordinate radius )
    : a_ ( radius ),
    b_ ( radius ),
    c_ ( radius )
    {
    };

    //Constructor for a 2D ellipse

    Ellipsoid ( const Coordinate a, const Coordinate b )
    : a_ ( a ),
    b_ ( b ),
    c_ ( 0. )
    {
    };

    //Implementation of the ellipsoid formula:
    // ((x-0.5)/a)^2 + ((y-0.5)/b)^2 + ((z-0.5)/c)^2 = 1.

    Coordinate eval_func ( const std::vector<Coordinate> &x,
                           MaterialNumber mat_num ) const
    {
        return 1. - x[0] * x[0]
                / ( a_ * a_ )
                - x[1] * x[1]
                / ( b_ * b_ );
    }

    // Computation of the gradiant of the ellipsoid formula

    std::vector<Coordinate> eval_grad ( const std::vector<Coordinate> &x,
                                        MaterialNumber mat_num ) const
    {
        std::vector<Coordinate> grad ( DIMENSION );
        grad[0] = -2. * ( x[0] ) / ( a_ * a_ );
        grad[1] = -2. * ( x[1] ) / ( b_ * b_ );
        return grad;
    }

    // Parameter describing the ellipsoid
    const Coordinate a_;
    const Coordinate b_;
    const Coordinate c_;
};

// Class for creating ellipsoids with holes.

class EllipsoidWithHole : public BoundaryDomainDescriptor
{
  public:

    //Constructor for a circle (or sphere) with a circle (or sphere) hole

    EllipsoidWithHole ( const Coordinate outer_radius,
                        const Coordinate inner_radius )
    : outer_a_ ( outer_radius ),
    inner_a_ ( inner_radius ),
    outer_b_ ( outer_radius ),
    inner_b_ ( inner_radius ),
    outer_c_ ( outer_radius ),
    inner_c_ ( inner_radius )
    {
    };

    //Constructor for a 2D ellipse with a ellipse hole

    EllipsoidWithHole ( const Coordinate outer_a, const Coordinate outer_b,
                        const Coordinate inner_a, const Coordinate inner_b )
    : outer_a_ ( outer_a ),
    inner_a_ ( inner_a ),
    outer_b_ ( outer_b ),
    inner_b_ ( inner_b ),
    outer_c_ ( 0. ),
    inner_c_ ( 0. )
    {
    };

    //Implementation of the ellipsoid formulas.
    // Everything with MaterialNumber 11 will be mapped to the inner ellipsoid
    // Everything with MaterialNumber 12 to the outer one.

    Coordinate eval_func ( const std::vector<Coordinate> &x,
                           MaterialNumber mat_num ) const
    {
        Coordinate a;
        Coordinate b;
        Coordinate c;
        if ( mat_num == 12 )
        {
            a = inner_a_;
            b = inner_b_;
            c = inner_c_;
        }
        else if ( mat_num == 11 )
        {
            a = outer_a_;
            b = outer_b_;
            c = outer_c_;
        }
        else
        {
            return 0.;
        }
        return 1. - ( x[0] )*( x[0] ) / ( a * a )-( x[1] )*( x[1] ) / ( b * b );
    }

    // Computation of the corresponding gradients.

    std::vector<Coordinate> eval_grad ( const std::vector<Coordinate> &x,
                                        MaterialNumber mat_num ) const
    {
        Coordinate a;
        Coordinate b;
        Coordinate c;
        if ( mat_num == 12 )
        {
            a = inner_a_;
            b = inner_b_;
            c = inner_c_;
        }
        else if ( mat_num == 11 )
        {
            a = outer_a_;
            b = outer_b_;
            c = outer_c_;
        }
        else
        {
            return std::vector<Coordinate>( DIMENSION, 0. );
        }
        std::vector<Coordinate> grad ( DIMENSION );
        grad[0] = -2. * ( x[0] ) / ( a * a );
        grad[1] = -2. * ( x[1] ) / ( b * b );
        return grad;
    }

    //Parameters to describe the ellipsoids.
    const Coordinate outer_a_;
    const Coordinate outer_b_;
    const Coordinate outer_c_;

    const Coordinate inner_a_;
    const Coordinate inner_b_;
    const Coordinate inner_c_;
};

// Functor to evaluate the exact solution u of the Poisson problem
// with Dirichlet BC, and its gradient \grad u.
// Parameters M, N and O in solution. These decide the period in x-, y- and
// z-direction respectively.
const double M = 2.;
const double N = 2.;
const double O = 0.5;

// Exact solution for Ellipsoid domain and square domain

struct ExactSol
{

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];
        const double y = ( DIMENSION > 1 ) ? pt[1] : 0;
        const double z = ( DIMENSION > 2 ) ? pt[2] : 0;
        double solution;
#ifdef ELLIPSOID_BOUNDARY
        // a and b are the axes of the ellipse geometry
        // unit circle a = b = 1
        const double a = 1.;
        const double b = 1.;
        const double k = 10.0;
        solution = -std::pow ( ( ( x / a )*( x / a )+( y / b )*( y / b ) ),
                               4. * k )
                + 1.;
#else
        const double pi = M_PI;
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
        }
#endif
        return solution;
    }

    // Partial derivations of the exact solution

    Vec<DIMENSION, double> eval_grad ( const Vec<DIMENSION, double>& pt ) const
    {
        Vec<DIMENSION, double> grad;
        const double x = pt[0];
        const double y = ( DIMENSION > 1 ) ? pt[1] : 0;

#ifdef ELLIPSOID_BOUNDARY
        const double a = 1.;
        const double b = 1.;
        const double k = 10.0;
        grad[0] = -4. * k * std::pow ( ( x / a )*( x / a )+( y / b )*( y / b ),
                                       ( 4. * k - 1. ) )
                * ( 2. / ( a * a ) )
                * x;

        grad[1] = -4. * k * std::pow ( ( x / a )*( x / a )+( y / b )*( y / b ),
                                       ( 4. * k - 1. ) )
                * ( 2. / ( b * b ) )
                * y;
#else
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
        }
#endif
        return grad;
    }
};

// Dirichlet boundary condition
// Functor used to impose u = 0 on the boundary.

struct DirichletZero
{

    std::vector<double> evaluate (
                                   const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face
                                   ) const
    {
        // return array with Dirichlet values for dof:s on boundary face
        if ( face.get_material_number ( ) == 11 )
        {
            return std::vector<double>( coords_on_face.size ( ), 0.0 );
        }
        // Neumann boundary values
        return std::vector<double>( 0, 0.0 );
    }
};
// Right hand side f

double f ( Vec<DIMENSION, double> pt )
{
    ExactSol sol;
    double rhs_sol;
#ifdef ELLIPSOID_BOUNDARY
    const double x = pt[0];
    const double y = ( DIMENSION > 1 ) ? pt[1] : 0;
    const double pi = M_PI;
    const double a = 1.;
    const double b = 1.;
    const double k = 10.0;
    rhs_sol = 8.
            * k
            * std::pow ( ( x * x ) / ( a * a )+( y * y ) / ( b * b ), 4. * k )
            * (
            a * a * a * a * ( 8. * k - 1 ) * y * y
            + a * a * b * b * ( x * x + y * y )
            + b * b * b * b * ( 8. * k - 1 ) * x * x
            )
            / (
            std::pow ( a * a * y * y + b * b * x*x, 2 )
            );
#else
#    ifdef CONST_F
    rhs_sol = 1;
#    else
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
    }
#    endif
#endif
    return rhs_sol;
}

// Assembling of the linear system
// Functor used for the local assembly of the stiffness matrix and load vector.

template<class ExactSol>
class LocalPoissonAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element
                ( element, quadrature );

        // Local stiffness matrix
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
                            wq
                            * dot ( grad_phi ( j, q ), grad_phi ( i, q ) )
                            * std::abs ( detJ ( q ) );
                }
            }
        }
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element
                ( element, quadrature );

        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const int n_dofs = num_dofs ( 0 );
            for ( int i = 0; i < n_dofs; ++i )
            {
                lv[dof_index ( i, 0 )] += wq
                        * f ( x ( q ) ) * phi ( i, q )
                        * std::abs ( detJ ( q ) );
            }
        }
    }

};

// Functor used for the local evaluation of the square of the L2-norm of the
// error on each element.

template<class ExactSol>
class L2ErrorIntegrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    L2ErrorIntegrator ( CVector& pp_sol )
    : pp_sol_ ( pp_sol )
    {
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element
                ( element, quadrature );

        // Evaluate the computed solution at all quadrature points.
        approx_sol_.clear ( );
        evaluate_fe_function ( pp_sol_, 0., approx_sol_ );

        const int num_q = num_quadrature_points ( );
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            // Delta is the quadrature point wise error.
            const double delta = sol_ ( x ( q ) ) - approx_sol_[q];
            value += wq * delta * delta * std::abs ( detJ ( q ) );
        }
    }

  private:
    // Coefficients of the computed solution
    const CVector& pp_sol_;
    // functor to evaluate exact solution
    ExactSol sol_;
    // Vector with values of computed solution evaluated at
    // each quadrature point
    FunctionValues< double > approx_sol_;
};

// Functor used for the local evaluation of the square of the H1-norm of the
// error on each element.

template<class ExactSol>
class H1ErrorIntegrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    H1ErrorIntegrator ( CVector& pp_sol )
    : pp_sol_ ( pp_sol )
    {
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element
                ( element, quadrature );
        // Evaluate the gradient of the computed solution at all
        // quadrature points.
        approx_grad_u_.clear ( );
        evaluate_fe_function_gradients ( pp_sol_, 0., approx_grad_u_ );

        const int num_q = num_quadrature_points ( );
        for ( int q = 0.; q < num_q; ++q )
        {
            const Vec<DIMENSION, double> grad_u = sol_.eval_grad ( x ( q ) );
            value += w ( q )
                    * (
                    dot ( grad_u, grad_u )
                    - 2. * dot ( grad_u, approx_grad_u_[q] )
                    + dot ( approx_grad_u_[q], approx_grad_u_[q] )
                    )
                    * std::abs ( detJ ( q ) );
        }
    }

  private:
    // Coefficients of the computed solution
    const CVector& pp_sol_;
    // Functor to evaluate exact solution
    ExactSol sol_;
    // Gradient of computed solution evaluated at each quadrature point
    FunctionValues< Vec<DIMENSION, double> > approx_grad_u_;
};

// Cell Term Assembler for inner cells

class CellTermAssembler : private AssemblyAssistant <DIMENSION, double>
{
  public:

    CellTermAssembler ( CVector& u_h )
    : u_h_ ( u_h )
    {
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            double& value )
    {

        AssemblyAssistant<DIMENSION, double>::initialize_for_element
                ( element, quadrature );

        approx_hess_u_h.clear ( );
        evaluate_fe_function_hessians ( u_h_, 0, approx_hess_u_h );

        const int num_q = num_quadrature_points ( );

        for ( int q = 0.; q < num_q; ++q )
        {
            const double laplace_u_h = trace ( approx_hess_u_h[q] );
            value += this->h ( ) * this->h ( )
                    * w ( q )
                    * ( f ( x ( q ) ) + laplace_u_h )
                    * ( f ( x ( q ) ) + laplace_u_h )
                    * std::abs ( detJ ( q ) );
        }
    }
    const CVector& u_h_;
    FunctionValues< Mat<DIMENSION, DIMENSION, double> > approx_hess_u_h;
};

// Cell Term Assembler for the semicircular boundary cells

class BoundCellTermAssembler : private DGAssemblyAssistant <DIMENSION, double>
{
  public:

    BoundCellTermAssembler ( CVector& u_h )
    : u_h_ ( u_h )
    {
    }

    void operator() ( const Element<double>& left_elem,
            const Element<double>& right_elem,
            const Quadrature<double>& left_quad,
            const Quadrature<double>& right_quad,
            int left_facet_number,
            int right_facet_number,
            InterfaceSide left_if_side,
            InterfaceSide right_if_side,
            double& local_val )
    {

        const bool is_boundary =
                (
                right_if_side == DGGlobalAssembler<double>::INTERFACE_BOUNDARY
                );

        if ( !is_boundary )
        {
            local_val = 0.;
            return;
        }

        initialize_for_interface ( left_elem, right_elem,
                                   left_quad, right_quad,
                                   left_facet_number, right_facet_number,
                                   left_if_side, right_if_side );

        const int num_q = num_quadrature_points ( );
        const int num_dofs_trial = trial ( ).num_dofs ( 0 );

        double vol = 0.;
        for ( int q = 0.; q < num_q; ++q )
        {
            vol += w ( q ) * std::abs ( ds ( q ) );
        }
        double h_A = std::pow ( vol, 1. / DIMENSION );

        /*
        double h_A = 0.;
        for (int i = 0; i < num_q; ++i) {
            for (int j = 0; j < num_q; ++j) {
                h_A = std::max(h_A, norm(x(i) - x(j)));
            }
        }
         */

        const double alpha = std::acos ( ( 1. - ( h_A * h_A ) ) / 2. );
        const double dA = 0.5 * ( alpha - std::sin ( alpha ) );
        const double h_S = h_A / 2. * std::tan ( 0.25 * alpha );
        Vec<DIMENSION, double> q_first = x ( 0 );
        Vec<DIMENSION, double> q_last = x ( num_q - 1 );
        Vec<DIMENSION, double> lot_point;
        lot_point[0] = 0.5 * ( q_first[0] + q_last[0] );
        lot_point[1] = 0.5 * ( q_first[1] + q_last[1] );
        Vec<DIMENSION, double> middle_point;
        middle_point = lot_point + 0.5 * h_S * trial ( ).n ( num_q / 2 );

        local_val += h_A * h_A * h_A * h_A
                * std::abs ( f ( middle_point ) * f ( middle_point ) )
                * dA * dA;
    }
    const CVector& u_h_;
};

// Jump Term Assembler for the jumps over the inner edges

class JumpTermAssembler : private DGAssemblyAssistant <DIMENSION, double>
{
  public:

    JumpTermAssembler ( CVector& u_h )
    : u_h_ ( u_h )
    {
    }

    void operator() ( const Element<double>& left_elem,
            const Element<double>& right_elem,
            const Quadrature<double>& left_quad,
            const Quadrature<double>& right_quad,
            int left_facet_number,
            int right_facet_number,
            InterfaceSide left_if_side,
            InterfaceSide right_if_side,
            double& local_val )
    {
        const bool is_boundary =
                (
                right_if_side == DGGlobalAssembler<double>::INTERFACE_BOUNDARY
                );

        if ( is_boundary
             || ( left_if_side == right_if_side ) )
        {
            local_val = 0.;
            return;
        }

        this->initialize_for_interface ( left_elem, right_elem,
                                         left_quad, right_quad,
                                         left_facet_number, right_facet_number,
                                         left_if_side, right_if_side );

        left_grad_u_h.clear ( );
        right_grad_u_h.clear ( );

        this->trial ( ).evaluate_fe_function_gradients ( u_h_, 0, left_grad_u_h );
        this->test ( ).evaluate_fe_function_gradients ( u_h_, 0, right_grad_u_h );

        const int num_q = num_quadrature_points ( );
        double vol = 0.;
        for ( int q = 0.; q < num_q; ++q )
        {
            vol += w ( q ) * std::abs ( ds ( q ) );
        }
        double h_E = std::pow ( vol, 1. / DIMENSION );

        /*
        double h_E = 0.;
        for (int i = 0; i < num_q; ++i) {
            for (int j = 0; j < num_q; ++j) {
                h_E = std::max(h_E, norm(this->x(i) - this->x(j)));
            }
        }
         */

        // Loop over quadrature points on each edge
        for ( int q = 0.; q < num_q; ++q )
        {
            const double dS = std::abs ( this->ds ( q ) );
            local_val += 0.25
                    * h_E
                    * this->w ( q )
                    * std::abs (
                                 dot ( this->trial ( ).n ( q ), left_grad_u_h[q] )
                                 + dot ( this->test ( ).n ( q ), right_grad_u_h[q] ) )
                    * std::abs (
                                 dot ( this->trial ( ).n ( q ), left_grad_u_h[q] )
                                 + dot ( this->test ( ).n ( q ), right_grad_u_h[q] ) )
                    * dS;
        }
    }
    const CVector& u_h_;
    FunctionValues< Vec<DIMENSION, double> > left_grad_u_h;
    FunctionValues< Vec<DIMENSION, double> > right_grad_u_h;
};

// Jump Term Assembler for the jumps over the outer edges

class BoundaryTermAssembler : private DGAssemblyAssistant <DIMENSION, double>
{
  public:

    BoundaryTermAssembler ( CVector& u_h )
    : u_h_ ( u_h )
    {
    }

    void operator() ( const Element<double>& left_elem,
            const Element<double>& right_elem,
            const Quadrature<double>& left_quad,
            const Quadrature<double>& right_quad,
            int left_facet_number,
            int right_facet_number,
            InterfaceSide left_if_side,
            InterfaceSide right_if_side,
            double& local_val )
    {

        const bool is_boundary =
                (
                right_if_side == DGGlobalAssembler<double>::INTERFACE_BOUNDARY
                );

        // To sort the inner and outer edges
        if ( !is_boundary )
        {
            local_val = 0.;
            return;
        }

        initialize_for_interface ( left_elem, right_elem,
                                   left_quad, right_quad,
                                   left_facet_number, right_facet_number,
                                   left_if_side, right_if_side );
        left_grad_u_h.clear ( );
        trial ( ).evaluate_fe_function_gradients ( u_h_, 0, left_grad_u_h );
        const int num_q = this->num_quadrature_points ( );
        const double C = is_boundary ? 1 : 0;

        double vol = 0.;
        for ( int q = 0.; q < num_q; ++q )
        {
            vol += w ( q ) * std::abs ( ds ( q ) );
        }
        double h_A = std::pow ( vol, 1. / DIMENSION );
        /*
        for (int i = 0; i < num_q; ++i) {
            for (int j = 0; j < num_q; ++j) {
                h_A = std::max(h_A, norm(x(i) - x(j)));
            }
        }
         * */
        // Loop over quadrature points on the edge
        for ( int q = 0.; q < num_q; ++q )
        {
            const double dS = std::abs ( ds ( q ) );
            local_val += C
                    * h_A * h_A
                    * w ( q )
                    * std::abs ( dot ( trial ( ).n ( q ), left_grad_u_h[q] ) )
                    * std::abs ( dot ( trial ( ).n ( q ), left_grad_u_h[q] ) )
                    * dS;
        }
    }
    const CVector& u_h_;
    FunctionValues< Vec<DIMENSION, double> > left_grad_u_h;
};
