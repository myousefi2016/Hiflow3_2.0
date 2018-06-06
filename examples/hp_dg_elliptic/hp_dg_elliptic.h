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

/// \author Staffan Ronnas

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
typedef LAD::VectorType Vector;
typedef LAD::MatrixType Matrix;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
#define DIMENSION 3

// Options
#define USE_DG
#define STRONG_BC

#define TRIG_SOL
//#define SING_SOL

// Parameters M and N in solution. These decide the period in x- and
// y-direction respectively.
const int M = 2;
const int N = 2;
const int L = 2;

// Functor to evaluate the exact solution u of the Poisson problem
// with Dirichlet BC, and its gradient \grad u.

#if DIMENSION == 2
static const double Ac[] = {
                            1., 0.,
                            0., 1.
};
#endif

#if DIMENSION == 3
// Coefficient operator -- should be spd.
static const double Ac[] = {
                            1., 0., 0.,
                            0., 1., 0.,
                            0., 0., 1.
};
#endif

#if DIMENSION == 2

struct ExactSol
{
    const Mat<DIMENSION, DIMENSION, double> A; // coefficient matrix

    const double c;

    ExactSol ( ) : A ( &Ac[0] ), c ( 0. )
    {
    }

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        return eval_sol ( pt );
    }

    double eval_sol ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];
        const double y = pt[1];
        const double pi = M_PI;

        return std::sin ( 2. * M * pi * x ) * std::sin ( 2. * N * pi * y );
    }

    Vec<DIMENSION> eval_grad ( const Vec<DIMENSION, double>& pt ) const
    {
        const double pi = M_PI;

        Vec<DIMENSION, double> grad;
        const double x = pt[0];
        const double y = pt[1];

        grad[0] = 2. * M * pi * std::cos ( 2. * M * pi * x ) *
                std::sin ( 2. * N * pi * y );
        grad[1] = 2. * N * pi * std::sin ( 2. * M * pi * x ) *
                std::cos ( 2. * N * pi * y );

        return grad;
    }

    double eval_f ( const Vec<DIMENSION, double>& pt ) const
    {
        return 4. * M_PI * M_PI * ( M * M + N * N ) * eval_sol ( pt );
    };

};

#endif

#if DIMENSION == 3

#    ifdef TRIG_SOL

struct ExactSol
{
    const Mat<DIMENSION, DIMENSION, double> A; // coefficient matrix

    const double c;

    const double F;

    ExactSol ( )
    : A ( &Ac[0] ), c ( 0. ), F ( 1 )
    {
    }

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        return eval_sol ( pt );
    }

    double eval_sol ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];
        const double y = pt[1];
        const double z = pt[2];
        const double pi = M_PI;

        return F * std::sin ( 2. * M * pi * x ) * std::sin ( 2. * N * pi * y ) * std::sin ( 2. * L * pi * z );
    }

    Vec<DIMENSION, double> eval_grad ( const Vec<DIMENSION, double>& pt ) const
    {
        const double pi = M_PI;

        Vec<DIMENSION, double> grad;
        const double x = pt[0];
        const double y = pt[1];
        const double z = pt[2];

        grad[0] = F * 2. * M * pi * std::cos ( 2. * M * pi * x ) *
                std::sin ( 2. * N * pi * y ) * std::sin ( 2. * L * pi * z );
        grad[1] = F * 2. * N * pi * std::sin ( 2. * M * pi * x ) *
                std::cos ( 2. * N * pi * y ) * std::sin ( 2. * L * pi * z );
        grad[2] = F * 2. * L * pi * std::sin ( 2. * M * pi * x ) *
                std::sin ( 2. * N * pi * y ) * std::cos ( 2. * L * pi * z );

        return grad;
    }

    Mat<DIMENSION, DIMENSION, double> eval_hessian ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];
        const double y = pt[1];
        const double z = pt[2];
        const double pi = M_PI;
        const double r = F * 4. * M_PI * M_PI;

        Mat<DIMENSION, DIMENSION, double> H;

        const double C[DIMENSION] = { M, N, L };
        const double si[DIMENSION] = { std::sin ( 2. * M * pi * x ),
                                      std::sin ( 2. * N * pi * y ),
                                      std::sin ( 2. * L * pi * z ) };
        const double co[DIMENSION] = { std::cos ( 2. * M * pi * x ),
                                      std::cos ( 2. * N * pi * y ),
                                      std::cos ( 2. * L * pi * z ) };

        H ( 0, 0 ) = -r * C[0] * C[0] * si[0] * si[1] * si[2];
        H ( 0, 1 ) = r * C[0] * C[1] * co[0] * co[1] * si[2];
        H ( 0, 2 ) = r * C[0] * C[2] * co[0] * si[1] * co[2];

        H ( 1, 0 ) = r * C[1] * C[0] * co[0] * co[1] * si[2];
        H ( 1, 1 ) = -r * C[1] * C[1] * si[0] * si[1] * si[2];
        H ( 1, 2 ) = r * C[1] * C[2] * si[0] * co[1] * co[2];

        H ( 2, 0 ) = r * C[2] * C[0] * co[0] * si[1] * co[2];
        H ( 2, 1 ) = r * C[2] * C[1] * si[0] * co[1] * co[2];
        H ( 2, 2 ) = -r * C[2] * C[2] * si[0] * si[1] * si[2];

        return H;
    }

    double eval_divAgrad ( const Vec<DIMENSION, double>& pt ) const
    {
        double val = 0.;
        const Mat<DIMENSION, DIMENSION, double> H = eval_hessian ( pt );
        for ( int i = 0; i < DIMENSION; ++i )
        {
            for ( int j = 0; j < DIMENSION; ++j )
            {
                val += A ( i, j ) * H ( i, j );
            }
        }
        return val;
    }

    double eval_f ( const Vec<DIMENSION, double>& pt ) const
    {
        const double u = eval_sol ( pt );
        const double divAgrad_u = eval_divAgrad ( pt );
        return -divAgrad_u + c*u;
    };

};

#    endif // TRIG_SOL

#    ifdef SING_SOL

struct ExactSol
{
    const Mat<DIMENSION, DIMENSION, double> A; // coefficient matrix

    const double ax, ay, az;

    const double c;

    ExactSol ( )
    : A ( &Ac[0] ), ax ( 0.8 ), ay ( 2.5 ), az ( 0. ), c ( 0. )
    {
    }

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        return eval_sol ( pt );
    }

    double eval_sol ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];
        const double y = pt[1];
        const double z = pt[2];

        return std::pow ( x, ax ) * std::pow ( y, ay ) * std::pow ( z, az );
    }

    Vec<DIMENSION, double> eval_grad ( const Vec<DIMENSION, double>& pt ) const
    {
        Vec<DIMENSION, double> grad;

        const double x = pt[0];
        const double y = pt[1];
        const double z = pt[2];

        const double u = eval_sol ( pt );

        if ( std::abs ( u ) > 1.e-15 )
        {
            if ( ax > 0. )
            {
                grad[0] = ax * u / x;
            }
            if ( ay > 0. )
            {
                grad[1] = ay * u / y;
            }
            if ( az > 0. )
            {
                grad[2] = az * u / z;
            }
        }

        return grad;
    }

    Mat<DIMENSION, DIMENSION, double> eval_hessian ( const Vec<DIMENSION, double>& pt ) const
    {
        const double x = pt[0];
        const double y = pt[1];
        const double z = pt[2];

        const double u = eval_sol ( pt );

        Mat<DIMENSION, DIMENSION, double> H;

        H ( 0, 0 ) = ax * ( ax - 1. ) * u / ( x * x );
        H ( 0, 1 ) = ax * ay * u / ( x * y );
        H ( 0, 2 ) = ax * az * u / ( x * z );

        H ( 1, 0 ) = H ( 0, 1 );
        H ( 1, 1 ) = ay * ( ay - 1. ) * u / ( y * y );
        H ( 1, 2 ) = ay * az * u / ( y * z );

        H ( 2, 0 ) = H ( 0, 2 );
        H ( 2, 1 ) = H ( 1, 2 );
        H ( 2, 2 ) = az * ( az - 1. ) * u / ( z * z );

        return H;
    }

    double eval_divAgrad ( const Vec<DIMENSION, double>& pt ) const
    {
        double val = 0.;
        const Mat<DIMENSION, DIMENSION, double> H = eval_hessian ( pt );
        for ( int i = 0; i < DIMENSION; ++i )
        {
            for ( int j = 0; j < DIMENSION; ++j )
            {
                val += A ( i, j ) * H ( i, j );
            }
        }
        return val;
    }

    double eval_f ( const Vec<DIMENSION, double>& pt ) const
    {
        const double u = eval_sol ( pt );
        const double divAgrad_u = eval_divAgrad ( pt );
        return -divAgrad_u + c*u;
    };

};

#    endif

#endif

// Functor used to impose u = 0 on the boundary.

struct DirichletZero
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector< std::vector<double> >& coords_on_face ) const
    {
        // return array with Dirichlet values for dof:s on boundary face
        return std::vector<double>( coords_on_face.size ( ), 0.0 );
    }
};

struct DirichletExact
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector< std::vector<double> >& coords_on_face ) const
    {
        std::vector<double> values ( coords_on_face.size ( ), 0. );
        Vec<DIMENSION, double> pt;

        for ( int i = 0; i < coords_on_face.size ( ); ++i )
        {
            for ( int c = 0; c < DIMENSION; ++c )
            {
                pt[c] = coords_on_face[i][c];
            }
            values[i] = sol_ ( pt );
        }
        return values;
    }

    ExactSol sol_;
};

// Functor used for the local evaluation of the square of the L2-norm of the
// error on each element.

template<class ExactSol>
class L2ErrorIntegrator : private AssemblyAssistant<DIMENSION, double>
{
  public:

    L2ErrorIntegrator ( const CoupledVector<Scalar>& pp_sol )
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
            const double delta = sol_.eval_sol ( x ( q ) ) - approx_sol_[q];
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

    H1ErrorIntegrator ( const CoupledVector<Scalar>& pp_sol )
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
            const Vec<DIMENSION, double> diff = grad_u - approx_grad_u_[q];

            for ( int c = 0; c < DIMENSION; ++c )
            {
                if ( std::isnan ( grad_u[c] ) )
                {
                    std::cerr << "exact grad_u is NaN\n";
                }
                if ( std::isnan ( approx_grad_u_[q][c] ) )
                {
                    std::cerr << "approx grad_u is NaN\n";
                }
            }

            value += w ( q ) * dot ( diff, diff ) * std::abs ( detJ ( q ) );

        }
        //    if (std::isnan(value))  { std::cerr << "value is NaN\n"; }
        if ( value < 0. )
        {
            if ( std::abs ( value ) < 1.e-16 )
            {
                value = 0.;
            }
            else
            {
                std::cerr << "Funny negative H1 error: " << value << "\n";
            }
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

class DGJumpMatrixAssembler : DGAssemblyAssistant<DIMENSION, double>
{
  public:

    DGJumpMatrixAssembler ( double theta, double gamma ) : theta_ ( theta ), gamma_ ( gamma )
    {
    }

    void operator() ( const Element<double>& left_elem, const Element<double>& right_elem,
            const Quadrature<double>& left_quad, const Quadrature<double>& right_quad,
            int left_facet_number, int right_facet_number,
            InterfaceSide left_if_side, InterfaceSide right_if_side,
            LocalMatrix& lm )
    {
        const bool is_boundary = ( right_if_side == DGGlobalAssembler<double>::INTERFACE_BOUNDARY );

#ifdef STRONG_BC
        if ( is_boundary )
        {
            lm.Zeros ( );
            return; // do nothing.
        }
#endif

        initialize_for_interface ( left_elem, right_elem,
                                   left_quad, right_quad,
                                   left_facet_number, right_facet_number,
                                   left_if_side, right_if_side );

        const int num_q = num_quadrature_points ( );
        const int num_dofs_trial = trial ( ).num_dofs ( 0 );
        const int num_dofs_test = test ( ).num_dofs ( 0 );

        assert ( lm.nrows ( ) >= num_dofs_trial );
        assert ( lm.ncols ( ) >= num_dofs_test );

        const double C = is_boundary ? 1 : 0.5; // TODO: not needed?

        const double h = std::pow ( std::abs ( trial ( ).detJ ( 0 ) ), 1. / double(DIMENSION ) );
        const double p = left_elem.get_fe_type ( 0 )->get_fe_deg ( );
        const double alpha = p * p / ( std::pow ( h, 2. ) );

        for ( int q = 0; q < num_q; ++q )
        {
            const double dS = std::abs ( ds ( q ) );
            const double W = w ( q ) * dS;

            const double dot_n = dot ( test ( ).n ( q ), trial ( ).n ( q ) );

            for ( int i = 0; i < num_dofs_test; ++i )
            {
                for ( int j = 0; j < num_dofs_trial; ++j )
                {
                    lm ( test ( ).dof_index ( i, 0 ), trial ( ).dof_index ( j, 0 ) ) +=
                            W * (
                            -C * dot ( sol_.A * test ( ).grad_phi ( i, q ), trial ( ).n ( q ) ) * trial ( ).phi ( j, q )
                            + C * theta_ * dot ( sol_.A * trial ( ).grad_phi ( j, q ), test ( ).n ( q ) ) * test ( ).phi ( i, q )
                            + gamma_ * alpha * test ( ).phi ( i, q ) * trial ( ).phi ( j, q ) * dot_n
                            );
                }
            }
        }
    }

  private:
    ExactSol sol_;
    double theta_;
    double gamma_;
};

class DGBoundaryVectorAssembler : DGAssemblyAssistant<DIMENSION, double>
{
  public:
    typedef DGGlobalAssembler<double>::InterfaceSide InterfaceSide;

    DGBoundaryVectorAssembler ( double theta, double gamma ) : theta_ ( theta ), gamma_ ( gamma )
    {
    }

    void operator() ( const Element<double>& left_elem, const Element<double>& right_elem,
            const Quadrature<double>& left_quad, const Quadrature<double>& right_quad,
            int left_facet_number, int right_facet_number,
            InterfaceSide left_if_side, InterfaceSide right_if_side,
            LocalVector& lv )
    {
        const bool is_boundary = ( right_if_side == DGGlobalAssembler<double>::INTERFACE_BOUNDARY );

        // Only add contributions to boundary facets.
        if ( !is_boundary ) return;

        initialize_for_interface ( left_elem, right_elem,
                                   left_quad, right_quad,
                                   left_facet_number, right_facet_number,
                                   left_if_side, right_if_side );

        const int num_q = num_quadrature_points ( );
        const int ndofs = test ( ).num_dofs ( 0 );

        const double h = std::pow ( std::abs ( trial ( ).detJ ( 0 ) ), 1. / double(DIMENSION ) );
        const double p = left_elem.get_fe_type ( 0 )->get_fe_deg ( );
        const double alpha = p * p / ( std::pow ( h, 2. ) );

        for ( int q = 0; q < num_q; ++q )
        {
            const double dS = std::abs ( ds ( q ) );
            const double W = w ( q ) * dS;
            const Vec<DIMENSION, double> pt = x ( q );

            const double sol_val = sol_.eval_sol ( pt );

            for ( int i = 0; i < ndofs; ++i )
            {

                lv[test ( ).dof_index ( i, 0 )] += W * (
                        // TODO: sign of these terms?
                        +theta_ * dot ( sol_.A * test ( ).grad_phi ( i, q ), test ( ).n ( q ) ) * sol_val
                        + gamma_ * alpha * sol_val * test ( ).phi ( i, q )
                        );
            }
        }

    }

  private:
    ExactSol sol_;
    double theta_;
    double gamma_;
};

// Functor used for the local assembly of the stiffness matrix and load vector.

class DGCellAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    void operator() ( const Element<double>& element,
            const Quadrature<double>& quadrature,
            LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // compute local matrix
        const int num_q = num_quadrature_points ( );
        const int n_dofs = num_dofs ( 0 );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            for ( int i = 0; i < n_dofs; ++i )
            {
                for ( int j = 0; j < n_dofs; ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                            wq * (
                            dot ( sol_.A * grad_phi ( j, q ), grad_phi ( i, q ) )
                            //+ sol_.c * phi(j, q) * phi(i, q)
                            ) * std::abs ( detJ ( q ) );
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
                lv[dof_index ( i, 0 )] += wq * sol_.eval_f ( x ( q ) ) * phi ( i, q ) * std::abs ( detJ ( q ) );
            }
        }
    }

  private:
    ExactSol sol_;
};

class VisualizeCellValues
{
  public:

    VisualizeCellValues ( const std::vector<double>& cell_values )
    : cell_values_ ( cell_values )
    {
    }

    void operator() ( const Entity& cell,
            const std::vector<double>& ref_coords,
            std::vector<double>& values ) const
    {
        values.clear ( );
        const int num_points = ref_coords.size ( ) / cell.gdim ( );
        values.resize ( num_points, cell_values_.at ( cell.index ( ) ) );
    }

  private:
    const std::vector<double>& cell_values_;

};

template<class F>
class EvalAnalyticFunction
{
  public:

    EvalAnalyticFunction ( const VectorSpace<double>& space, const F fun ) : space_ ( space ), fun_ ( fun )
    {
    }

    void operator() ( const Entity& cell,
            const std::vector<double>& ref_coords,
            std::vector<double>& values ) const
    {
        const int gdim = cell.gdim ( );
        const int num_points = ref_coords.size ( ) / gdim;

        const CellTransformation<double>& cell_trans = space_.GetCellTransformation ( cell );

        std::vector<double> refp ( DIMENSION, 0. );
        Vec<DIMENSION, double> physp;
        values.resize ( num_points, 0. );
        for ( int i = 0; i < num_points; ++i )
        {
            for ( int c = 0; c < DIMENSION; ++c )
            {
                refp[c] = ref_coords.at ( gdim * i + c );
            }
            physp[0] = cell_trans.x ( refp );
            physp[1] = cell_trans.y ( refp );
            if ( DIMENSION == 3 )
            {
                physp[2] = cell_trans.z ( refp );
            }
            values[i] = fun_ ( physp );
        }
    }

    const VectorSpace<double>& space_;
    const F fun_;
};
