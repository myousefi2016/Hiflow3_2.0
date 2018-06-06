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

/// \author Ulrich Blunck, Tobias Hahn, Julian Kraemer

#ifndef FLOW_TUTORIAL_H
#    define FLOW_TUTORIAL_H

#    include <algorithm>
#    include <mpi.h>
#    include <time.h>
#    include <stdio.h>
#    include <exception>
#    include <fstream>
#    include <string>
#    include <stdexcept>
#    include <vector>

#    include "hiflow.h"

using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

#    define PI M_PI;
#    define DIMENSION 2    //Define the spatial dimensions the problem is solved for
#    define debug 1

// Linear Algebra type renaming.
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;

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

//////////////// Boundary conditions ////////////////////////////////
/// \brief class defining boundary conditions for 2-dimensional case

struct ChannelFlowBC
{
    /// \brief constructor
    ///
    /// \param [in] var  index of variable
    /// \param [in] D  channel diameter
    /// \param [in] Um  maximum inflow
    /// \param [in] inflow_bdy  material number of inflow boundary
    /// \param [in] outflow_bdy  material number of outflow boundary

    ChannelFlowBC ( int var, double D, double Um, int inflow_bdy,
                    int outflow_bdy )
    : var_ ( var ), R_ ( D / 2 ), Um_ ( Um ), inflow_bdy_ ( inflow_bdy ),
    outflow_bdy_ ( outflow_bdy )
    {
        assert ( var_ == 0 || var_ == 1 );
    }

    /// \brief evaluates values of all Dirichlet degrees of freedom which
    /// 				lie on face.
    ///
    /// \param [in] face face on which values of dirichlet dofs will be
    /// 							evaluated
    /// \param [in] coords_on_face coordinates of dofs on face

    std::vector<double> evaluate ( const Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        // stores values of Dirichlet dofs
        std::vector<double> values;

        // material number of face needed for evaluation
        const int material_num = face.get_material_number ( );

        // variables to distinguish type of boundary
        const bool outflow = ( material_num == outflow_bdy_ );
        const bool inflow = ( material_num == inflow_bdy_ );

        // Following Dirichlet boundary conditions are set: if the face lies on an
        // inflow boundary. On outflow, no Dirichlet boundaries are applied. In case
        // of no outflow we distinguish and if so set
        // u_x = 4*Um * y * (1-y) / D^2 and u_y = 0. Otherwise, set u_x = u_y = 0 .

        // check if face lies on the outflow boundary or not
        if ( !outflow )
        {
            // depending on number of dofs on face, set size of values
            values.resize ( coords_on_face.size ( ) );

            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                // evaluate dirichlet function at each point and store coordinates of
                // each dof
                const Coord& pt = coords_on_face[i];

                if ( inflow )
                {
                    if ( var_ == 0 )
                    { // x-component
                        // values[i] = 4. * Um_ * (pt[1] - 0.35) * (0.65 - pt[1])/(D_ * D_);
                        values[i] = -Um_ * ( ( pt[1] / R_ - 1 ) * ( pt[1] / R_ - 1 ) - 1 );
                    }
                    else if ( var_ == 1 )
                    {
                        // y-component
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
    // index of variable
    const int var_;
    // radius of channel
    const double R_;
    // maximum inflow velocity
    const double Um_;
    // material number of inflow and outflow boundary
    const int inflow_bdy_, outflow_bdy_;
};

/// \brief class defining boundary conditons of 3 dimensional case

struct ChannelFlowBC3d
{
    /// \brief constructor
    ///
    /// \param [in] var  index of variable
    /// \param [in] D  channel diameter
    /// \param [in] Um  maximum inflow
    /// \param [in] inflow_bdy  material number of inflow boundary
    /// \param [in] outflow_bdy  material number of outflow boundary

    ChannelFlowBC3d ( int var, double D, double Um, int inflow_bdy,
                      int outflow_bdy )
    : var_ ( var ), R_ ( D / 2 ), Um_ ( Um ), inflow_bdy_ ( inflow_bdy ),
    outflow_bdy_ ( outflow_bdy )
    {
        assert ( var_ == 0 || var_ == 1 || var_ == 2 );
        assert ( DIMENSION == 3 );
    }

    /// \brief evaluates values of all Dirichlet degrees of freedom which
    ///					lie on face.
    ///
    /// \param [in] face face on which values of dirichlet dofs will be
    ///								evaluated
    /// \param [in] coords_on_face coordinates of dofs on face

    std::vector<double> evaluate ( const Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        // stores values of Dirichlet dofs
        std::vector<double> values;

        // material number of face needed for evaluation
        const int material_num = face.get_material_number ( );

        // variables to distinguish type of boundary
        const bool outflow = ( material_num == outflow_bdy_ );
        const bool inflow = ( material_num == inflow_bdy_ );

        // Dirichlet boundary conditions. Check whether
        // the face lies on an inflow boundary, and if so set
        // u_x = -4 * Um *((y-.5)^2 + (z-.5)^2 - .25*D_^2)/(D_^2) and u_y = u_z = 0.

        // check if face lies on the outflow boundary or not
        if ( !outflow )
        {
            // depending on number of dofs on face, set size of values
            values.resize ( coords_on_face.size ( ) );

            // loop over points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                // evaluate dirichlet function at each point and store coordinates of
                // each dof
                const Coord& pt = coords_on_face[i];

                if ( inflow )
                {
                    if ( var_ == 0 )
                    {
                        // x-component
                        values[i] = -Um_ * ( pow ( ( pt[1] - .5 ), 2 ) + pow ( ( pt[2] - .5 ), 2 ) - R_ * R_ )
                                / ( R_ * R_ );
                    }
                    else if ( var_ == 1 || var_ == 2 )
                    {
                        // y / z -component
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

    // index of variable
    const int var_;
    // radius of channel
    const double R_;
    // maximum inflow velocity
    const double Um_;
    // material number of inflow and outflow boundary
    const int inflow_bdy_, outflow_bdy_;
};
////////////////// End boundary conditions /////////////////////////////

////////////////// Assembler ////////////////////////////////

class PorousMediaAssembler : private AssemblyAssistant<DIMENSION, double>
{
    /// \brief constructor
    ///
    /// \param [in] solution vector with values of the current solution
    /// \param [in] nu viscosity
    /// \param [in] rho density
    /// \param [in] u_av inflow speed
    /// \param [in] l_av reference height
    /// \param [in] conti contribution weight of the continuum equation
    ///							 and the GPME
    /// \param [in] eps_free porosity for the part of the domain
    ///							 where no porous media is
    /// \param [in] eps_por porosity of the porous media
    /// \param [in] kap_free peremeability for the part of the
    ///							 domain where no porous media is
    /// \param [in] kap_por peremeability of the porous media
    /// \param [in] geometry used geometry

  public:

    PorousMediaAssembler ( const CVector& solution, double nu, double rho,
                           double u_av, double l_av, double conti,
                           double eps_free, double eps_por,
                           double kap_free, double kap_por,
                           const std::string geometry
                           )
    :
    solution_ ( solution ),
    nu_ ( nu / rho ),
    inv_rho_ ( 1. / ( rho ) ),
    ref_u ( u_av ),
    ref_l ( l_av ),
    ref_p ( rho*u_av*u_av ),
    inv_reynolds ( nu / ( u_av*l_av ) ),
    eps_free_ ( eps_free ),
    eps_por_ ( eps_por ),
    kap_free_ ( kap_free ),
    kap_por_ ( kap_por ),
    conti_ ( conti ),
    geometry_ ( geometry )
    {
    }

    /// \brief operator defining the local system matrix on each element
    ///
    /// \param [in] element element for which the system matrix is
    ///								computed
    /// \param [in]	quadrature used quadrature for integrations
    /// \param [out] lm variable in which the matrix is stored in for
    ///								 the element

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // interpolate previous values of solution and the gradient of the solution
        // to quadrature points
        for ( int v = 0; v < DIMENSION; ++v )
        {
            prev_vel_[v].clear ( );
            grad_prev_vel_[v].clear ( );
            evaluate_fe_function ( solution_, v, prev_vel_[v] );
            evaluate_fe_function_gradients ( solution_, v, grad_prev_vel_[v] );
        }

        const int num_q = num_quadrature_points ( );

        double eps_local, kap_local;

        // enforce local constant porosity (epsilon) and permeability (kappa)
        // on each element
        double eps = 0;
        double kap = 0;
        for ( int q = 0; q < num_q; ++q )
        {
            eps_local = evaluate_epsilon ( x ( q ) );
            kap_local = evaluate_kappa ( x ( q ) );
            if ( eps_local > eps )
            {
                eps = eps_local;
            }
            if ( kap_local > kap )
            {
                kap = kap_local;
            }
        }

        // compute the needed constants for computation of local element matrix
        double inveps = 1. / eps;
        double poweps = 1. / pow ( eps, ( 3 / 2 ) );
        inv_darcy = ref_l * ref_l / kap;
        inv_darcy_root = sqrt ( inv_darcy );

        // loop over quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            // compute weight for quadrature points
            const double wq = w ( q );
            // compute determinant of Jacobi matrix for each quadrature point
            const double dJ = std::abs ( detJ ( q ) );

            // get previous solution in vector form
            Vec<DIMENSION, double> vel_k;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_k[var] = prev_vel_[var][q];
            }

            // compute norm of previous solution
            norm_ = sqrt ( dot ( vel_k, vel_k ) );

            // assemble L1(p, v) = - \int{p div{v}}
            // index of pressure variable
            const int p_var = DIMENSION;
            // loop over variables of test functions
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    // loop over ansatz functions
                    for ( int j = 0; j < num_dofs ( p_var ); ++j )
                    {
                        // compute integral
                        lm ( dof_index ( i, v_var ), dof_index ( j, p_var ) ) +=
                                -wq * ( phi ( j, q, p_var ) *
                                grad_phi ( i, q, v_var )[v_var] ) * dJ;
                    }
                }
            }

            // assemble L2(u,v) = nu / eps * \int {\grad(u) : \grad(v)}
            // loop over variables of test functions and ansatz functions
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    // loop over ansatz functions
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        // compute integral
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=

                                wq * ( inv_reynolds * inveps
                                * dot ( grad_phi ( j, q, u_var ), grad_phi ( i, q, u_var ) ) )
                                * dJ;
                    }
                }
            }

            // assemble L3(u,v) = \int { nu/kappa * v }
            // loop over variables of test functions and ansatz functions
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    // loop over ansatz functions
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        // compute integral
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * ( inv_reynolds * inv_darcy
                                * phi ( i, q, u_var ) ) * dJ;
                    }
                }
            }

            // assemble L4(u, q) = \int{q div(u)}
            // index of pressure variable
            const int q_var = DIMENSION;
            // loop over variables of ansatz functions
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( q_var ); ++i )
                {
                    // loop over ansatz functions
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        // compute integral
                        lm ( dof_index ( i, q_var ), dof_index ( j, u_var ) ) +=
                                wq / ref_u * conti_ * ( phi ( i, q, q_var )
                                * grad_phi ( j, q, u_var )[u_var] ) * dJ;
                    }
                }
            }

            // assemble N1_1(u,v) = \int { (vel_k*\grad{u})*v } / eps
            // loop over variables of test functions and ansatz functions
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    // loop over ansatz functions
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        // compute integral
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * ( inveps * inveps * dot ( vel_k, grad_phi ( j, q, u_var ) )
                                * phi ( i, q, u_var ) ) * dJ;
                    }
                }
            }

            // assemble N1_2(u,v) = \int { (u\grad{u_k}*v } / eps^2
            // loop over variables of test functions
            for ( int test_var = 0; test_var < DIMENSION; ++test_var )
            {
                // loop over variables of ansatz functions
                for ( int trial_var = 0; trial_var < DIMENSION; ++trial_var )
                {
                    // loop over test functions
                    for ( int i = 0; i < num_dofs ( test_var ); ++i )
                    {
                        // loop over ansatz functions
                        for ( int j = 0; j < num_dofs ( trial_var ); ++j )
                        {
                            // compute integral
                            lm ( dof_index ( i, test_var ), dof_index ( j, trial_var ) ) +=
                                    wq * ( inveps * inveps * grad_prev_vel_[test_var][q][trial_var]
                                    * phi ( j, q, trial_var )
                                    * phi ( i, q, test_var ) ) * dJ;
                        }
                    }
                }
            }

            // assemble N2(u,v) =
            // \int{C * norm_ * v} * 1.75/(sqrt(150) * sqrt(kappa) * eps^(3/2))
            // loop over variables of test functions and ansatz functions
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    // loop over ansatz functions
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        // compute integral
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * ( inv_darcy_root * 1.75 / ( sqrt ( 150 ) ) * poweps * norm_
                                * phi ( i, q, u_var ) ) * dJ;
                    }
                }
            }
        }
    }

    /// \brief operator defining the local residual on each element
    ///
    /// \param [in] element element for which the system matrix is
    ///								computed
    /// \param [in]	quadrature used quadrature for integrations
    /// \param [out] lm variable in which the residual is stored in for
    ///								 the element

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // interpolate previous values of solution and the gradient of the solution
        // to quadrature points
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

        double eps_local, kap_local;

        // enforce local constant porosity (epsilon) and permeability (kappa)
        // on each element
        double eps = 0;
        double kap = 0;
        for ( int q = 0; q < num_q; ++q )
        {
            eps_local = evaluate_epsilon ( x ( q ) );
            kap_local = evaluate_kappa ( x ( q ) );
            if ( eps_local > eps )
            {
                eps = eps_local;
            }
            if ( kap_local > kap )
            {
                kap = kap_local;
            }
        }

        // compute the needed constants for computation of local element residual
        double inveps = 1. / eps;
        double poweps = 1. / pow ( eps, ( 3 / 2 ) );
        inv_darcy = ref_l * ref_l / kap;
        inv_darcy_root = sqrt ( inv_darcy );

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

            // compute norm of previous solution
            norm_ = sqrt ( dot ( vel_k, vel_k ) );
            // if (norm_ < 0.05) norm_=.05;

            // L1(v) = -1/rho*\int(p_k*div(v))
            // loop over variables of test functions
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    // compute integral
                    lv[dof_index ( i, v_var )] +=
                            -wq * ( pressure_k_[q] * grad_phi ( i, q, v_var )[v_var] ) * dJ;
                }
            }

            // L2(v) = nu / eps * \int( \grad{u_k} : \grad{v} )
            // loop over variables of test functions
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    // compute integral
                    lv[dof_index ( i, v_var )] +=
                            wq * ( inv_reynolds * inveps * dot ( grad_phi ( i, q, v_var ),
                                                                 grad_prev_vel_[v_var][q] ) ) * dJ;
                }
            }

            // L3(v) = nu/kappa*\int(v)
            // loop over variables of test functions
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    // compute integral
                    lv[dof_index ( i, v_var )] +=
                            wq * ( inv_reynolds * inv_darcy * phi ( i, q, v_var ) ) * dJ;
                }
            }

            // L4(q) = \int(q * div(u_k))
            // index of pressure variable
            const int q_var = DIMENSION;
            // variable for divergence
            double div_u_k = 0.;
            // computing the divergence
            for ( int d = 0; d < DIMENSION; ++d )
            {
                div_u_k += grad_prev_vel_[d][q][d];
            }
            // loop over test functions
            for ( int i = 0; i < num_dofs ( q_var ); ++i )
            {
                // compute integral
                lv[dof_index ( i, q_var )] +=
                        wq / ref_u * conti_ * ( div_u_k * phi ( i, q, q_var ) ) * dJ;
            }

            // N1(v) = \int(u_k*\grad{u_k}*v)
            // loop over variables of test functions
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    // compute integral
                    lv[dof_index ( i, v_var )] +=
                            wq * ( inveps * inveps * dot ( grad_prev_vel_[v_var][q], vel_k )
                            * phi ( i, q, v_var ) ) * dJ;
                }
            }

            // N2(v) = 1.75/(sqrt(150)*kappa*epsilon^(3/2)) * norm(u_k) * v
            // loop over variables of test functions
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                // loop over test functions
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    // compute integral
                    lv[dof_index ( i, v_var )] +=
                            wq * ( inv_darcy_root * 1.75 / ( sqrt ( 150 ) ) * poweps * norm_
                            * phi ( i, q, v_var ) ) * dJ;
                }
            }
        }
    }

    /// \brief evaluates epsilon depending on geometry
    ///
    /// \param [in] pt coordinates of point which will be evaluated

    double evaluate_epsilon ( const Vec<DIMENSION, double> pt )
    {
        // return value
        double ret;
        // coordinates
        double x = pt[0];
        double y = pt[1];
        double z = 0;
        if ( DIMENSION == 3 )
        {
            z = pt[2];
        }

        // set ret to the default value eps_por_, which is te porosity of the porous
        // media. Correct it later to eps_free_ if pt is in a free-flow zone.
        ret = eps_por_;

        // the geometry channel has no free flow zone
        if ( geometry_ == "Channel" )
        {
            ;
        }

        if ( geometry_ == "Column" )
        {
            // Inflow/Outflow are free flow zones, correct ret
            if ( ( x < 1 ) or (x > 5 ) )
            {
                ret = eps_free_;
            }
            // the free flow zone of the inflow is longer then the inflow tube
            if ( DIMENSION == 2 )
            {
                if ( ( x >= 1 ) and (x <= 1.35 ) and (y < .65 ) and (y > .35 ) )
                {
                    ret = eps_free_;
                }
            }
            else
            {
                // in the 3-dimensional case, the free flow zone of the inflow is longer
                // then the inflow tube as well.
                if ( ( x >= 1 ) and (x <= 1.25 ) and (pow ( ( y - .5 ), 2 ) + pow ( ( z - .5 ), 2 )
                     < pow ( .15, 2 ) ) )
                {
                    ret = eps_free_;
                }
            }
            ;
        }
        return ret;
    }

    /// \brief evaluates kappa depending on geometry
    ///
    /// \param [in] pt coordinates of point which will be evaluated

    double evaluate_kappa ( Vec<DIMENSION, double> pt )
    {
        // return value
        double ret;
        // coordinates
        double x = pt[0];
        double y = pt[1];
        double z = 0;
        if ( DIMENSION == 3 )
        {
            z = pt[2];
        }
        // set ret to the default value kap_por_, which is te permeability of the
        // porous media. Correct it later to kap_free_ if pt is in a free-flow zone.
        ret = kap_por_;

        // the geometry channel has no free flow zone
        if ( geometry_ == "Channel" )
        {
            ;
        }

        if ( geometry_ == "Column" )
        {
            // Inflow/Outflow are free flow zones, correct ret
            if ( ( x < 1 ) or (x > 5 ) )
            {
                ret = eps_free_;
            }
            // the free flow zone of the inflow is longer then the inflow tube
            if ( DIMENSION == 2 )
            {
                if ( ( x >= 1 ) and (x <= 1.35 ) and (y < .65 ) and (y > .35 ) )
                {
                    ret = eps_free_;
                }
            }
            else
            {
                // in the 3-dimensional case, the free flow zone of the inflow is longer
                // then the inflow tube as well.
                if ( ( x >= 1 ) and (x <= 1.25 ) and (pow ( ( y - .5 ), 2 ) + pow ( ( z - .5 ), 2 )
                     < pow ( .15, 2 ) ) )
                {
                    ret = eps_free_;
                }
            }
            ;
        }
        return ret;
    }

  private:
    // vector containing solution of last Newton step
    const CVector& solution_;
    // constants
    double nu_, inv_rho_, norm_, ref_u, ref_l, ref_p, inv_reynolds, inv_darcy,
    eps_free_, eps_por_, kap_free_, kap_por_, inv_darcy_root;
    // contribution weight
    double conti_;
    // interpolation of previous solution for velocity
    FunctionValues<double> prev_vel_[DIMENSION];
    // interpolation of previous solution for pressure
    FunctionValues<double> pressure_k_;
    // interpolation of gradient of previous solution
    FunctionValues< Vec<DIMENSION, double> > grad_prev_vel_[DIMENSION];
    // string to distinguish the different geometries
    const std::string geometry_;
};

#endif
