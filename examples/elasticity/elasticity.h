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

/// \author Nicolai Schoch

#ifndef ELASTICITY_H
#    define ELASTICITY_H

#    include <mpi.h>

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

#    define DIMENSION 3

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

struct StationaryElasticity_DirichletBC_3D
{
    // Parameters:
    // bdy - material number of boundary

    StationaryElasticity_DirichletBC_3D ( int var, int bdy1, int bdy2, int bdy3 )
    : var_ ( var ), bdy1_ ( bdy1 ), bdy2_ ( bdy2 ), bdy3_ ( bdy3 )
    {
        assert ( var_ == 0 || var_ == 1 || var_ == 2 );
        assert ( DIMENSION == 3 );
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        //Return array with Dirichlet values for dof:s on boundary face.
        std::vector<double> values;

        const int material_num = face.get_material_number ( );

        if ( material_num == bdy1_ )
        { // NOTE: this is the fixed Dirichlet BC part, and holds during the whole simulation.
            values.resize ( coords_on_face.size ( ) );

            // loop over dof points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                // evaluate dirichlet function at each point
                values[i] = 0.;
            }
        }
        else if ( material_num == bdy2_ )
        { // NOTE: this is a displacement Dirichlet BC, and holds ONLY for the time dt*timestep <= 1.0.
            values.resize ( coords_on_face.size ( ) );

            // loop over dof points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                // evaluate dirichlet function at each point
                values[i] = 0.;
                if ( var_ == 0 )
                {
                    values[i] = 0.;
                }
                else if ( var_ == 1 )
                {
                    values[i] = 0.04;
                }
                else if ( var_ == 2 )
                {
                    values[i] = 0.02;
                }
            }

        }
        else if ( material_num == bdy3_ )
        { // NOTE: this is a displacement Dirichlet BC, and holds ONLY for the time dt*timestep <= 1.0.
            values.resize ( coords_on_face.size ( ) );

            // loop over dof points on the face
            for ( int i = 0; i < static_cast < int > ( coords_on_face.size ( ) ); ++i )
            {
                // evaluate dirichlet function at each point
                values[i] = 0.;
                if ( var_ == 0 )
                {
                    values[i] = 0.;
                }
                else if ( var_ == 1 )
                {
                    values[i] = 0.04;
                }
                else if ( var_ == 2 )
                {
                    values[i] = 0.;
                }
            }

        }
        else
        {
            assert ( 0 );
        }

        return values;
    }

    const int var_;
    const int bdy1_, bdy2_, bdy3_;
};

//////////////// Stationary assembler (SystemMatrix and RHS) ////////////////////////////////

class StationaryElasticityAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    StationaryElasticityAssembler ( double lambda, double mu, double rho, double gravity )
    : lambda_ ( lambda ), mu_ ( mu ), rho_ ( rho ), gravity_ ( gravity )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );

        // loop over quadrature points.
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // assemble \int {lambda \div(u) \div(phi)}
            for ( int test_var = 0; test_var < DIMENSION; ++test_var )
            {
                for ( int trial_var = 0; trial_var < DIMENSION; ++trial_var )
                {
                    for ( int i = 0; i < num_dofs ( test_var ); ++i )
                    {
                        for ( int j = 0; j < num_dofs ( trial_var ); ++j )
                        {
                            lm ( dof_index ( i, test_var ), dof_index ( j, trial_var ) ) +=
                                    wq * lambda_ * grad_phi ( j, q, trial_var )[trial_var] * grad_phi ( i, q, test_var )[test_var] * dJ;
                        }
                    }
                }
            }

            // assemble \int {mu [ \frob(\nabla(u)\nabla(\phi)) + \frob(\nabla(u)^T\nabla(\phi)) ] }
            for ( int var = 0; var < DIMENSION; ++var )
            {
                for ( int i = 0; i < num_dofs ( var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( var ); ++j )
                    {
                        for ( int var_frob = 0; var_frob < DIMENSION; ++var_frob )
                        {
                            lm ( dof_index ( i, var ), dof_index ( j, var ) ) +=
                                    wq * mu_ * ( grad_phi ( j, q, var )[var_frob] ) * grad_phi ( i, q, var )[var_frob] * dJ;
                            lm ( dof_index ( i, var ), dof_index ( j, var_frob ) ) +=
                                    wq * mu_ * ( grad_phi ( j, q, var_frob )[var] ) * grad_phi ( i, q, var )[var_frob] * dJ;
                        }
                    }
                }
            }

        }
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        const int num_q = num_quadrature_points ( );

        const double source[3] = { 0., gravity_ * rho_, 0. };

        // loop over quadrature points.
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // assemble l(v) = \rho * \int( f_ext \phi )
            for ( int test_var = 0; test_var < DIMENSION; ++test_var )
            {
                for ( int i = 0; i < num_dofs ( test_var ); ++i )
                {
                    lv[dof_index ( i, test_var )] +=
                            wq * source[test_var] * phi ( i, q, test_var ) * dJ;
                }
            }
        }
    }

  private:
    double lambda_, mu_, rho_, gravity_;
};

#endif
