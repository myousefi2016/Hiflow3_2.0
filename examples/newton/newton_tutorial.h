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

#ifndef HIFLOW_NEWTON_TUTORIAL_H
#    define HIFLOW_NEWTON_TUTORIAL_H

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

#    define DIMENSION 2

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

            Vec<DIMENSION, double> vel_k;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_k[var] = prev_vel_[var][q];
            }

            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            wq * ( nu_ * dot ( grad_phi ( i, q, v_var ), grad_prev_vel_[v_var][q] ) ) * dJ;
                }
            }

            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            wq * ( dot ( grad_prev_vel_[v_var][q], vel_k ) * phi ( i, q, v_var ) ) * dJ;
                }
            }

            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            -wq * ( inv_rho_ * pressure_k_[q] * grad_phi ( i, q, v_var )[v_var] ) * dJ;
                }
            }

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
        }
    }

  private:
    const CVector& solution_;
    double nu_, inv_rho_;
    FunctionValues<double> prev_vel_[DIMENSION];
    FunctionValues<double> pressure_k_;
    FunctionValues< Vec<DIMENSION, double> > grad_prev_vel_[DIMENSION];
};

#endif
