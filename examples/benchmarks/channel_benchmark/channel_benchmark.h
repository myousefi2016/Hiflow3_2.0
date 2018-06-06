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

#ifndef HIFLOW_CHANNEL_BENCHMARK_H
#    define HIFLOW_CHANNEL_BENCHMARK_H

#    include <mpi.h>

#    include <exception>
#    include <fstream>
#    include <string>
#    include <stdexcept>
#    include <vector>

#    include "hiflow.h"

/// \author Staffan Ronnas

using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

#    define DIMENSION 3
// debugging example with constructed right hand side and exact solution,
// see explication of ChannelBenchmark::eval_exact_sol() in source code
// EXACTSOL == 1 else == 0
#    define EXACTSOL 0

// Linear Algebra type renaming.
typedef LADescriptorCoupledD LAD;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;
typedef LAD::DataType Scalar;

typedef std::vector<double> Coord;

// Exception types

struct UnexpectedParameterValue : public std::runtime_error
{

    UnexpectedParameterValue ( const std::string& name,
                               const std::string& value )
    : std::runtime_error (
                           "Unexpected value '" + value + "' for parameter " + name )
    {
    }
};

/// The default quadrature selection chooses a quadrature rule that is accurate to 2 * max(fe_degree).

struct QuadratureSelection
{

    QuadratureSelection ( int order ) : order_ ( order )
    {
    }

    void operator() ( const Element<double>& elem, Quadrature<double>& quadrature )
    {
        const FEType<double>::FiniteElement fe_id = elem.get_fe_type ( 0 )->get_my_id ( );

        switch ( fe_id )
        {
            case FEType<double>::LAGRANGE_TRI:
                quadrature.set_cell_type ( 2 );
                quadrature.set_quadrature_by_order ( "GaussTriangle", order_ );
                break;
            case FEType<double>::LAGRANGE_QUAD:
                quadrature.set_cell_type ( 3 );
                quadrature.set_quadrature_by_order ( "GaussQuadrilateral", order_ );
                break;
            case FEType<double>::LAGRANGE_TET:
                quadrature.set_cell_type ( 4 );
                quadrature.set_quadrature_by_order ( "GaussTetrahedron", order_ );
                break;
            case FEType<double>::LAGRANGE_HEX:
                quadrature.set_cell_type ( 5 );
                quadrature.set_quadrature_by_order ( "GaussHexahedron", order_ );
                break;
            default:
                assert ( false );
        };
    }

    int order_;
};

struct ChannelFlowBC2d
{
    // Parameters:
    // var - variable
    // H - channel height, W - channel width
    // Um - maximum inflow
    // inflow_bdy - material number of inflow boundary
    // outflow_bdy - material number of outflow boundary

    ChannelFlowBC2d ( int var, double H, double Um, int inflow_bdy, int outflow_bdy )
    : var_ ( var ), H_ ( H ), Um_ ( Um ), inflow_bdy_ ( inflow_bdy ), outflow_bdy_ ( outflow_bdy )
    {
        assert ( var_ == 0 || var_ == 1 );
        assert ( DIMENSION == 2 );
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );

        const bool outflow = ( material_num == outflow_bdy_ );
        const bool inflow = ( material_num == inflow_bdy_ );

        if ( !outflow )
        {
            // All boundaries except outflow have Dirichlet BC.
            values.resize ( coords_on_face.size ( ) );

            // loop over dof points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                // evaluate dirichlet function at each point
                const Coord& pt = coords_on_face[i];

                if ( inflow )
                {
                    if ( var_ == 0 )
                    { // x-component
                        values[i] = 4. * Um_ * pt[1] * ( H_ - pt[1] ) / ( H_ * H_ );
                    }
                    else if ( var_ == 1 )
                    { // y-components
                        values[i] = 0.;
                    }
                    else
                    {
                        assert ( false );
                    }
                }
                else
                {
                    // not inflow: u = 0
                    values[i] = 0.;
                }
            }
        }
        return values;
    }

    const int var_;
    const double H_; // size in y- direction, respectively.
    const double Um_; // max inflow velocity
    const int inflow_bdy_, outflow_bdy_;
};

class CylinderDescriptor : public BoundaryDomainDescriptor
{
  public:

    CylinderDescriptor ( Coordinate radius, std::vector<Coordinate> center ) :
    radius_ ( radius ), center_ ( center )
    {
        assert ( center_.size ( ) > 1 );
    }

    Coordinate eval_func ( const std::vector<Coordinate>& p, MaterialNumber mat_num ) const
    {
        if ( mat_num == 13 )
        {
            const Coordinate x = p[0] - center_[0];
            const Coordinate y = p[1] - center_[1];
            return radius_ * radius_ - x * x - y*y;
        }
        else
        {
            return 0.;
        }
    }

    std::vector<Coordinate> eval_grad ( const std::vector<Coordinate>& p, MaterialNumber mat_num ) const
    {
        if ( mat_num == 13 )
        {
            const Coordinate x = p[0] - center_[0];
            const Coordinate y = p[1] - center_[1];
            std::vector<Coordinate> grad ( DIMENSION, 0. );
            grad[0] = -2 * x;
            grad[1] = -2 * y;
            return grad;
        }
        else
        {
            return std::vector<Coordinate>( DIMENSION, 0. );
        }
    }

  private:
    const Coordinate radius_;
    const std::vector<Coordinate> center_;
};

struct ChannelFlowBC3d
{
    // Parameters:
    // var - variable
    // H - channel height, W - channel width
    // Um - maximum inflow
    // inflow_bdy - material number of inflow boundary
    // outflow_bdy - material number of outflow boundary

    ChannelFlowBC3d ( int var, double W, double H, double Um, int inflow_bdy, int outflow_bdy )
    : var_ ( var ), W_ ( W ), H_ ( H ), Um_ ( Um ), inflow_bdy_ ( inflow_bdy ), outflow_bdy_ ( outflow_bdy )
    {
        assert ( var_ == 0 || var_ == 1 || var_ == 2 );
        assert ( DIMENSION == 3 );
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );

        const bool outflow = ( material_num == outflow_bdy_ );
        const bool inflow = ( material_num == inflow_bdy_ );

        // Dirichlet boundary conditions. Check whether
        // the face lies on an inflow boundary, and if so set
        // u_x = 4*Um * y * (1-y) / H^2 and u_y = 0. Otherwise, set u_x = u_y = 0 .

        if ( !outflow )
        {
            // All boundaries except outflow have Dirichlet BC.
            values.resize ( coords_on_face.size ( ) );

            // loop over dof points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                // evaluate dirichlet function at each point
                const Coord& pt = coords_on_face[i];

                if ( inflow )
                {
                    if ( var_ == 0 )
                    { // x-component
                        values[i] = 16. * Um_ * pt[1] * pt[2] * ( W_ - pt[1] ) * ( H_ - pt[2] ) / ( W_ * W_ * H_ * H_ );
                    }
                    else if ( var_ == 1 || var_ == 2 )
                    { // y- and z-components
                        values[i] = 0.;
                    }
                    else
                    {
                        assert ( false );
                    }
                }
                else
                {
                    // not inflow: u = 0
                    values[i] = 0.;
                }
            }
        }
        return values;
    }

    const int var_;
    const double W_, H_; // size in y- and z- direction, respectively.
    const double Um_; // max inflow velocity
    const int inflow_bdy_, outflow_bdy_;
};

// boundary conditions for example with constructed exact solution
// u = (z-y, x-z, y-x) and p = x + 2y + 3z (divergence free, do nothing boundary
// condition does not hold hence complete Dirichlet boundary must be used.
// Used for debugging purpose to exclude errors in the mesh and dof numbering
// Use this as start solution for stationary case, if no bug, residual of newton
// method should fulfill stop criteria from the beginning.

struct ExactSolChannelFlowBC3d
{
    // Parameters:
    // var - variable

    ExactSolChannelFlowBC3d ( int var )
    : var_ ( var )
    {
        assert ( var_ == 0 || var_ == 1 || var_ == 2 );
        assert ( DIMENSION == 3 );
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        // Dirichlet boundary conditions everywhere.
        // Set u_x = z-y, u_y = x-z, u_z = y-x

        // All boundaries except outflow have Dirichlet BC.
        values.resize ( coords_on_face.size ( ) );

        // loop over dof points on the face
        for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
        {
            // evaluate dirichlet function at each point
            const Coord& pt = coords_on_face[i];

            if ( var_ == 0 )
            { // x-component
                values[i] = pt[2] - pt[1];
            }
            else if ( var_ == 1 )
            {
                values[i] = pt[0] - pt[2];
            }
            else if ( var_ == 2 )
            { // y- and z-components
                values[i] = pt[1] - pt[0];
            }
            else
            {
                assert ( false );
            }
        }
        return values;
    }

    const int var_;
};
//////////////// End boundary conditions /////////////////////////////

// Functor to evaluate the exact solution, constructed example for debugging
// purposes, see explication in source code of function
// ChannelBenchmark::eval_exact_sol(),
// u = (z-y, x-z, y-x) and p = x + 2y + 3z of the Navier-Stokes problem
// corresponding to the boundary values of the struct ExactSolChannelFlowBC3d

struct ExactSol3D
{

    double operator() ( const Vec<DIMENSION, double>& pt, int var ) const
    {
        if ( var == 0 )
        {
            return pt[2] - pt[1];
        }
        else if ( var == 1 )
        {
            return pt[0] - pt[2];
        }
        else if ( var == 2 )
        {
            return pt[1] - pt[0];
        }
        else if ( var == 3 )
        {
            return pt[0] + 2 * pt[1] + 3 * pt[2];
        }
        else
        {
            assert ( false );
            return -1.e99;
        }
    }
};

//////////////// Stationary assembler ////////////////////////////////

class StationaryFlowAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    StationaryFlowAssembler ( const CVector& solution, double nu, double rho )
    : solution_ ( solution ), nu_ ( nu ), inv_rho_ ( 1. / rho )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // recompute previous solution values
        for ( int v = 0; v < DIMENSION; ++v )
        {
            prev_vel_[v].clear ( );
            grad_prev_vel_[v].clear ( );
            evaluate_fe_function ( solution_, v, prev_vel_[v] );
            evaluate_fe_function_gradients ( solution_, v, grad_prev_vel_[v] );
        }

        const int num_q = num_quadrature_points ( );

        // loop q
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // get previous solution in vector form
            Vec<DIMENSION, double> vel_k;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_k[var] = prev_vel_[var][q];
            }

            // assemble a1(u,v) = \int {\grad(u) : \grad(v)}
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * ( nu_ * dot ( grad_phi ( j, q, u_var ), grad_phi ( i, q, u_var ) ) ) * dJ;
                    }
                }
            }

            // assemble a2(u,v) = \int { (vel_k*\grad{u})*v }
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * ( dot ( vel_k, grad_phi ( j, q, u_var ) ) * phi ( i, q, u_var ) ) * dJ;
                    }
                }
            }

            // assemble a3(u,v) = \int { (u\grad{u_k}*v }
            for ( int test_var = 0; test_var < DIMENSION; ++test_var )
            {
                for ( int trial_var = 0; trial_var < DIMENSION; ++trial_var )
                {
                    for ( int i = 0; i < num_dofs ( test_var ); ++i )
                    {
                        for ( int j = 0; j < num_dofs ( trial_var ); ++j )
                        {
                            lm ( dof_index ( i, test_var ), dof_index ( j, trial_var ) ) +=
                                    wq * ( grad_prev_vel_[test_var][q][trial_var] *
                                    phi ( j, q, trial_var ) *
                                    phi ( i, q, test_var ) ) * dJ;
                        }
                    }
                }
            }

            // assemble b(p, v) = - \int{p div{v}}
            const int p_var = DIMENSION;
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( p_var ); ++j )
                    {
                        lm ( dof_index ( i, v_var ), dof_index ( j, p_var ) ) +=
                                -wq * ( inv_rho_ * phi ( j, q, p_var ) *
                                grad_phi ( i, q, v_var )[v_var] ) * dJ;
                    }
                }
            }

            // assemble bT(u, q) = \int{q div(u)}
            const int q_var = DIMENSION;
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( q_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, q_var ), dof_index ( j, u_var ) ) +=
                                wq * ( inv_rho_ * phi ( i, q, q_var ) *
                                grad_phi ( j, q, u_var )[u_var] ) * dJ;
                    }
                }
            }
        }
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // recompute previous solution values
        for ( int v = 0; v < DIMENSION; ++v )
        {
            prev_vel_[v].clear ( );
            grad_prev_vel_[v].clear ( );
            evaluate_fe_function ( solution_, v, prev_vel_[v] );
            evaluate_fe_function_gradients ( solution_, v, grad_prev_vel_[v] );
        }
        pressure_k_.clear ( );
        evaluate_fe_function ( solution_, DIMENSION, pressure_k_ );

        const int num_q = num_quadrature_points ( );

        // loop over quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // get previous solution in vector form
            Vec<DIMENSION, double> vel_k;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_k[var] = prev_vel_[var][q];
            }

            // l1(v) = \nu * \int( \grad{u_k} : \grad{v} )
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            wq * ( nu_ * dot ( grad_phi ( i, q, v_var ), grad_prev_vel_[v_var][q] ) ) * dJ;
                }
            }

            // l2(v) = \int(u_k*\grad{u_k}*v)
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            wq * ( dot ( grad_prev_vel_[v_var][q], vel_k ) * phi ( i, q, v_var ) ) * dJ;
                }
            }

            // l3(v) = -1/rho*\int(p_k*div(v))
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            -wq * ( inv_rho_ * pressure_k_[q] * grad_phi ( i, q, v_var )[v_var] ) * dJ;
                }
            }

            // l4(q) = \int(q * div(u_k))
            const int q_var = DIMENSION;
            double div_u_k = 0.;
            for ( int d = 0; d < DIMENSION; ++d )
            {
                div_u_k += grad_prev_vel_[d][q][d];
            }

            for ( int i = 0; i < num_dofs ( q_var ); ++i )
            {
                lv[dof_index ( i, q_var )] +=
                        wq * ( inv_rho_ * div_u_k * phi ( i, q, q_var ) ) * dJ;
            }
            // debugging example with constructed right hand side and exact solution,
            // see explication of ChannelBenchmark::eval_exact_sol() in source code
#    if EXACTSOL == 1
            // l5(v) = f(x,y,z) * v
            // f = (  inv_rho - 2x +  y +  z,
            //      2*inv_rho +  x - 2y +  z,
            //      3*inv_rho +  x +  y - 2z)
            Vec<DIMENSION, double> point = x ( q );
            int v_var = 0;
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] -=
                        wq * ( inv_rho_ - 2 * point[0] + point[1] + point[2] )
                        * phi ( i, q, v_var ) * dJ;
            }
            v_var = 1;
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] -=
                        wq * ( inv_rho_ * 2 + point[0] - 2 * point[1] + point[2] )
                        * phi ( i, q, v_var ) * dJ;
            }
            v_var = 2;
            for ( int i = 0; i < num_dofs ( v_var ); ++i )
            {
                lv[dof_index ( i, v_var )] -=
                        wq * ( inv_rho_ * 3 + point[0] + point[1] - 2 * point[2] )
                        * phi ( i, q, v_var ) * dJ;
            }
#    endif
        }
    }

  private:
    const CVector& solution_;
    double nu_, inv_rho_;
    FunctionValues<double> prev_vel_[DIMENSION];
    FunctionValues<double> pressure_k_;
    FunctionValues< Vec<DIMENSION, double> > grad_prev_vel_[DIMENSION];
};

class StationaryFlowMassVelocityAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // assembel mass matrix for velocity components
            // m(vel,phi) = int { vel * phi}
            // loop over variables u,v,w belonging to velocity
            for ( int vel_var = 0; vel_var < DIMENSION; ++vel_var )
            {
                for ( int i = 0; i < num_dofs ( vel_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( vel_var ); ++j )
                    {
                        lm ( dof_index ( i, vel_var ), dof_index ( j, vel_var ) ) +=
                                wq * phi ( j, q, vel_var ) * phi ( i, q, vel_var ) * dJ;
                    }
                }
            }
        }
    }
};

class StationaryFlowStiffnessVelocityAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // assembel stiffness matrix for velocity components
            // a(vel,phi) = int { grad(vel) * grad(phi)}
            // loop over variables u,v,w belonging to velocity
            for ( int vel_var = 0; vel_var < DIMENSION; ++vel_var )
            {
                for ( int i = 0; i < num_dofs ( vel_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( vel_var ); ++j )
                    {
                        lm ( dof_index ( i, vel_var ), dof_index ( j, vel_var ) ) +=
                                wq * dot ( grad_phi ( j, q, vel_var ), grad_phi ( i, q, vel_var ) ) * dJ;
                    }
                }
            }
        }
    }
};

//////////////// InstationaryAssembler ////////////////////////////////

class InstationaryFlowAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    InstationaryFlowAssembler ( double nu, double rho )
    : nu_ ( nu ), inv_rho_ ( 1. / rho )
    {
    }

    void set_newton_solution ( const CVector* newton_sol )
    {
        prev_newton_sol_ = newton_sol;
    }

    void set_time_solution ( const CVector* time_sol )
    {
        prev_time_sol_ = time_sol;
    }

    void set_time_stepping_weights ( double alpha1, double alpha2, double alpha3 )
    {
        alpha1_ = alpha1;
        alpha2_ = alpha2;
        alpha3_ = alpha3;
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        // recompute previous solution values
        for ( int d = 0; d < DIMENSION; ++d )
        {
            vel_ns_[d].clear ( );
            grad_vel_ns_[d].clear ( );

            evaluate_fe_function ( *prev_newton_sol_, d, vel_ns_[d] );
            evaluate_fe_function_gradients ( *prev_newton_sol_, d, grad_vel_ns_[d] );
        }

        // indices j -> trial variable, i -> test variable
        // basis functions \phi -> velocity components, \eta -> pressure

        const int num_q = num_quadrature_points ( );

        // loop quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // Vector with velocity of previous newton point
            Vec<DIMENSION, double> vel_ns;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_ns[var] = vel_ns_[var][q];
            }

            // Velocity terms symmetric in test and trial variables
            for ( int var = 0; var < DIMENSION; ++var )
            {
                const int n_dofs = num_dofs ( var );
                for ( int i = 0; i < n_dofs; ++i )
                {
                    for ( int j = 0; j < n_dofs; ++j )
                    {
                        lm ( dof_index ( i, var ), dof_index ( j, var ) ) +=
                                wq *
                                // a0(\phi_i,\phi_j) = \int{ dot(\phi_j,\phi_i) }
                                ( phi ( j, q, var ) * phi ( i, q, var ) +

                                // a1(\phi_i, \phi_j) = alpha1 * \int{ \nu \grad{phi_j} : \grad{\phi_i} }
                                alpha1_ * nu_ * dot ( grad_phi ( j, q, var ), grad_phi ( i, q, var ) ) +

                                // c1(\phi_i,\phi_j) = \alpha1 * \int { (vel_ns*\grad{\phi_j}) * \phi_i }
                                alpha1_ * dot ( vel_ns, grad_phi ( j, q, var ) ) * phi ( i, q, var ) )
                                * dJ;
                    }
                }
            }

            // c2(\phi_i,\phi_j) = \int \alpha1 * { (\phi_j*\grad{vel_ns}*\phi_i }
            for ( int test_var = 0; test_var < DIMENSION; ++test_var )
            {
                for ( int trial_var = 0; trial_var < DIMENSION; ++trial_var )
                {
                    for ( int i = 0; i < num_dofs ( test_var ); ++i )
                    {
                        for ( int j = 0; j < num_dofs ( trial_var ); ++j )
                        {
                            lm ( dof_index ( i, test_var ), dof_index ( j, trial_var ) ) +=
                                    wq * alpha1_ * (
                                    grad_vel_ns_[test_var][q][trial_var] *
                                    phi ( j, q, trial_var ) *
                                    phi ( i, q, test_var ) ) * dJ;
                        }
                    }
                }
            }

            // b(\phi_i, \eta_j) = -\alpha2 / rho * \int{\eta_j div{\phi_i}}
            const int p_var = DIMENSION;
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( p_var ); ++j )
                    {
                        lm ( dof_index ( i, v_var ), dof_index ( j, p_var ) ) +=
                                -wq * alpha2_ * inv_rho_ * ( phi ( j, q, p_var ) * grad_phi ( i, q, v_var )[v_var] ) * dJ;
                    }
                }
            }

            // bT(\eta_i, \phi_j) = 1/rho \int{ \eta_i * div{\phi_j}}
            const int q_var = DIMENSION;
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( q_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, q_var ), dof_index ( j, u_var ) ) +=
                                wq * inv_rho_ * ( grad_phi ( j, q, u_var )[u_var] * phi ( i, q, q_var ) ) * dJ;
                    }
                }
            }
        } // end loop q
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        // recompute previous solution values
        for ( int d = 0; d < DIMENSION; ++d )
        {
            vel_ns_[d].clear ( );
            vel_ts_[d].clear ( );
            grad_vel_ns_[d].clear ( );
            grad_vel_ts_[d].clear ( );

            evaluate_fe_function ( *prev_newton_sol_, d, vel_ns_[d] );
            evaluate_fe_function ( *prev_time_sol_, d, vel_ts_[d] );
            evaluate_fe_function_gradients ( *prev_newton_sol_, d, grad_vel_ns_[d] );
            evaluate_fe_function_gradients ( *prev_time_sol_, d, grad_vel_ts_[d] );
        }

        // recompute pressure
        p_ns_.clear ( );
        evaluate_fe_function ( *prev_newton_sol_, DIMENSION, p_ns_ );

        // indices j -> trial variable, i -> test variable
        // basis functions \phi -> velocity components, \eta -> pressure

        const int num_q = num_quadrature_points ( );

        // loop quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // get previous newton and time step solutions in vector form
            Vec<DIMENSION, double> vel_ts, vel_ns;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_ns[var] = vel_ns_[var][q];
                vel_ts[var] = vel_ts_[var][q];
            }

            // Residual without incompressibility
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            wq *
                            // l0(\phi_i) = \int {dot(vel_ns - vel_ts, \phi_i)}
                            ( ( vel_ns[v_var] - vel_ts[v_var] ) * phi ( i, q, v_var ) +

                            // l1(\phi_i) =
                            // \nu * \int{(\alpha1\grad{vel_ns} + \alpha3\grad{vel_ts}) * \grad{\phi_i}}
                            nu_ * ( dot ( ( alpha1_ * grad_vel_ns_[v_var][q] ) +
                                          ( alpha3_ * grad_vel_ts_[v_var][q] ), grad_phi ( i, q, v_var ) ) ) +

                            // l2(\phi_i) = \int{ (\alpha1*dot(vel_ns, \grad{\vel_ns}
                            //                   + \alpha3*dot(vel_ts, \grad{\vel_ts}) * \phi_i}
                            ( alpha1_ * dot ( grad_vel_ns_[v_var][q], vel_ns ) +
                            alpha3_ * dot ( grad_vel_ts_[v_var][q], vel_ts ) ) * phi ( i, q, v_var ) +

                            // l3(\phi_i) = -\alpha2/rho*\int{p_ns * div(\phi_i)}
                            -( alpha2_ * inv_rho_ * p_ns_[q] * grad_phi ( i, q, v_var )[v_var] ) )
                            * dJ;
                }
            }

            // Incompressibility term
            const int q_var = DIMENSION;
            double div_u_k = 0.;
            for ( int d = 0; d < DIMENSION; ++d )
            {
                div_u_k += grad_vel_ns_[d][q][d];
            }

            // l4(\eta_i) = 1/rho*\int{\eta_i * div(vel_ns)}
            for ( int i = 0; i < num_dofs ( q_var ); ++i )
            {
                lv[dof_index ( i, q_var )] +=
                        wq * ( inv_rho_ * div_u_k * phi ( i, q, q_var ) ) * dJ;
            }
        } // end quadrature loop
    }

  private:
    const CVector* prev_time_sol_;
    const CVector* prev_newton_sol_;

    double alpha1_, alpha2_, alpha3_;

    double nu_, inv_rho_;
    FunctionValues<double> vel_ns_[DIMENSION]; // velocity at previous newton step
    FunctionValues<double> vel_ts_[DIMENSION]; // velocity at previous timestep
    FunctionValues<double> p_ns_; // pressure at previous newton step
    FunctionValues< Vec<DIMENSION, double> > grad_vel_ns_[DIMENSION]; // gradient of velocity at previous newton step
    FunctionValues< Vec<DIMENSION, double> > grad_vel_ts_[DIMENSION]; // gradient of velocity at previous timestep
};

// Functor used for the local evaluation of the square of the L2-norm
// of a set of variables on each element.

class L2NormIntegratorPp : private AssemblyAssistant<DIMENSION, double>
{
  public:

    L2NormIntegratorPp ( CoupledVector<Scalar>& pp_sol, const std::vector<int>& vars )
    : pp_sol_ ( pp_sol ), vars_ ( vars )
    {
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );
        for ( int v = 0, end_v = vars_.size ( ); v != end_v; ++v )
        {
            const int var = vars_[v];
            evaluate_fe_function ( pp_sol_, var, approx_sol_ );

            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const double dJ = std::abs ( detJ ( q ) );

                value += wq * approx_sol_[q] * approx_sol_[q] * dJ;
            }
        }
    }

  private:
    // coefficients of the computed solution
    CoupledVector<Scalar>& pp_sol_;
    // vector with values of computed solution evaluated at each quadrature point
    FunctionValues< double > approx_sol_;
    // variables for which to compute the norm
    std::vector<int> vars_;
};

// Functor used for the local evaluation of the square of the H1-seminorm
// of a set of variables on each element.

class H1semiNormIntegratorPp : private AssemblyAssistant<DIMENSION, double>
{
  public:

    H1semiNormIntegratorPp ( CoupledVector<Scalar>& pp_sol, const std::vector<int>& vars )
    : pp_sol_ ( pp_sol ), vars_ ( vars )
    {
    }

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );
        for ( int v = 0, end_v = vars_.size ( ); v != end_v; ++v )
        {
            const int var = vars_[v];
            evaluate_fe_function_gradients ( pp_sol_, var, approx_sol_grad_ );

            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const double dJ = std::abs ( detJ ( q ) );

                value += wq * dot ( approx_sol_grad_[q], approx_sol_grad_[q] ) * dJ;
            }
        }
    }

  private:
    // coefficients of the computed solution
    CoupledVector<Scalar>& pp_sol_;
    // vector with values of computed solution evaluated at each quadrature point
    FunctionValues< Vec<DIMENSION, double> > approx_sol_grad_;
    // variables for which to compute the norm
    std::vector<int> vars_;
};

// Functor used for the local evaluation of the square of the L2-norm of the
// difference between the solution of last end penultimate timestep.

class InstationaryL2ErrorIntegrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    InstationaryL2ErrorIntegrator ( const LAD::VectorType& coeff, const LAD::VectorType& coeff_penult )
    : coeff_ ( coeff ), coeff_penult_ ( coeff_penult )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );
        for ( int var = 0; var < DIMENSION; ++var )
        {
            evaluate_fe_function ( coeff_, var, approx_sol_ );
            evaluate_fe_function ( coeff_penult_, var, approx_sol_penult_ );
            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const double dJ = std::abs ( detJ ( q ) );
                const double delta = approx_sol_penult_[q] - approx_sol_[q];
                value += wq * delta * delta * dJ;
            }
        }
    }

  private:
    // coefficients of soluition of last timestep
    const LAD::VectorType& coeff_;
    // coefficients of solution of penultimate timestep
    const LAD::VectorType& coeff_penult_;
    // vector with values of computed solution evaluated at each quadrature point
    // for last and penultimate timestep
    FunctionValues< double > approx_sol_, approx_sol_penult_;
};

// Functor used for the local evaluation of the square of the L2-norm of the
// solution in one timestep on each element.

class InstationaryL2Integrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    InstationaryL2Integrator ( const LAD::VectorType& coeff )
    : coeff_ ( coeff )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            double& value )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );
        for ( int var = 0; var < DIMENSION; ++var )
        {
            evaluate_fe_function ( coeff_, var, approx_sol_ );
            for ( int q = 0; q < num_q; ++q )
            {
                const double wq = w ( q );
                const double dJ = std::abs ( detJ ( q ) );
                value += wq * approx_sol_[q] * approx_sol_[q] * dJ;
            }
        }
    }

  private:
    // coefficients of soluition of a timestep
    const LAD::VectorType& coeff_;
    // vector with values of computed solution evaluated at each quadrature point
    // for a timestep
    FunctionValues< double > approx_sol_;
};

#endif
