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

#ifndef jumpterm_assembler_h
#    define jumpterm_assembler_h

#    include <cmath>
#    include <fstream>
#    include <iostream>
#    include <string>

#    include "hiflow.h"
//#include "met_flow.h"

#    ifndef DIM_VAR
#        define DIM_VAR

#        ifdef APPDIM
const static int DIM = APPDIM;
#        else
const static int DIM = 2;
#        endif

#    endif

///
/// \file jumpterm_assembler.h
/// \brief Assembler for jumps of the solution over element interfaces.
///
/// \author Martin Baumann
///

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;
using namespace hiflow::la;

#    ifdef USE_HYPRE
typedef LADescriptorHypreD LAD;
#    else
typedef LADescriptorCoupledD LAD;
#    endif

///
/// \brief Auxiliary construction for dual JumpTermAssembler with left-hand and right-hand-sided AssemblyAssistants. 
///

template<int DIM, int MY_N>
class IndexedAssemblyAssistant : protected AssemblyAssistant<DIM, double>
{

    enum Index
    {
        N = MY_N
    };
};

typedef IndexedAssemblyAssistant<DIM, 0> LeftAA;
typedef IndexedAssemblyAssistant<DIM, 1> RightAA;

///
/// \brief Assembles jumps of discrete solution at element interfaces.
/// \details JumpTermAssembler is a dual AssemblyAssistant, consisting of a left-hand-sided and a right-hand-sided AssemblyAssistant.
///

class JumpTermAssembler : private LeftAA, private RightAA
{
  public:
    typedef StandardGlobalAssemblerWithJumps<double>::InterfaceSide InterfaceSide;

    JumpTermAssembler ( )
    {
        ;
    }

    void set_sol ( const LAD::VectorType* sol )
    {
        sol_ = sol;
    }

    void operator() ( const Element<double>& left_elem, const Element<double>& right_elem,
            const Quadrature<double>& left_quad, const Quadrature<double>& right_quad,
            int left_facet_number, int right_facet_number,
            InterfaceSide left_if_side, InterfaceSide right_if_side,
            double& value )
    {
        int debug_level = 0;

        LeftAA::initialize_for_facet ( left_elem, left_quad, left_facet_number );
        RightAA::initialize_for_facet ( right_elem, right_quad, right_facet_number );

        //    evaluate_fe_function_gradients();

        assert ( LeftAA::num_quadrature_points ( ) == RightAA::num_quadrature_points ( ) );

        const int num_q = LeftAA::num_quadrature_points ( );
        const int num_dofs_left = LeftAA::num_dofs ( 0 );
        const int num_dofs_right = RightAA::num_dofs ( 0 );

        const bool is_boundary = ( right_if_side == StandardGlobalAssemblerWithJumps<double>::INTERFACE_BOUNDARY );
        //const double C = is_boundary ? 1. : 0.5;   // TODO: not needed?

        FunctionValues<Vec<DIM, double> > left_gradients_u;
        LeftAA::evaluate_fe_function_gradients ( *sol_, 0, left_gradients_u );

        FunctionValues<Vec<DIM, double> > right_gradients_u;
        RightAA::evaluate_fe_function_gradients ( *sol_, 0, right_gradients_u );

        FunctionValues<Vec<DIM, double> > left_gradients_v;
        LeftAA::evaluate_fe_function_gradients ( *sol_, 1, left_gradients_v );

        FunctionValues<Vec<DIM, double> > right_gradients_v;
        RightAA::evaluate_fe_function_gradients ( *sol_, 1, right_gradients_v );

        // NB: h can only be known for a cell -- never for an edge
        //double h = std::pow(std::abs(LeftAA::detJ(0)), 1./double(DIMENSION));

        for ( int q = 0; q < num_q; ++q )
        {
            const double dS =
                    ( is_boundary || left_if_side == StandardGlobalAssemblerWithJumps<double>::INTERFACE_SLAVE ) ?
                    std::abs ( LeftAA::ds ( q ) ) : std::abs ( RightAA::ds ( q ) );
            const double W = LeftAA::w ( q ) * dS;

            if ( std::abs ( LeftAA::w ( q ) - RightAA::w ( q ) ) > 1.e-12 )
            {
                if ( debug_level > 0 )
                    std::cerr << "check w1 = " << LeftAA::w ( q ) << ", w2 = " << RightAA::w ( q ) << "\n";
            }

#    if 1
            if ( std::abs ( LeftAA::ds ( q ) - RightAA::ds ( q ) ) > 1.e-12 )
            {
                if ( debug_level > 0 )
                    std::cerr << "check ds1 = " << LeftAA::ds ( q ) << ", ds2 = " << RightAA::ds ( q ) << "\n";
            }
#    endif
            if ( debug_level > 0 )
                std::cout << "dS = " << dS << std::endl;

            if ( !is_boundary )
            {

                if ( left_if_side != right_if_side )
                {
                    if ( norm ( LeftAA::n ( q ) + RightAA::n ( q ) ) > 1.e-9 )
                    {
                        if ( debug_level > 0 )
                            std::cerr << "check opposite:  n1 = " << LeftAA::n ( q ) << ", n2 = " << RightAA::n ( q ) << "\n";
                    }
                }
                else
                {
                    if ( norm ( LeftAA::n ( q ) - RightAA::n ( q ) ) > 1.e-9 )
                    {
                        if ( debug_level > 0 )
                            std::cerr << "check same:  n1 = " << LeftAA::n ( q ) << ", n2 = " << RightAA::n ( q ) << "\n";
                    }
                }

                if ( norm ( LeftAA::x ( q ) - RightAA::x ( q ) ) > 1.e-9 )
                {
                    if ( debug_level > 0 )
                        std::cerr << "check:  x1 = " << LeftAA::x ( q ) << ", x2 = " << RightAA::x ( q ) << "\n";
                }
            }

            //            std::cerr << "nl*nr = " << dot(LeftAA::n(q), RightAA::n(q)) << "\n";

            //      for (int i = 0; i < num_dofs_right; ++i)
            //      {
            //        for (int j = 0; j < num_dofs_left; ++j)
            //        {
            value += W * ( dot ( left_gradients_u[q], LeftAA::n ( q ) )
                    * dot ( right_gradients_u[q], RightAA::n ( q ) )
                    + dot ( left_gradients_v[q], LeftAA::n ( q ) )
                    * dot ( right_gradients_v[q], RightAA::n ( q ) )
                    );
            if ( debug_level > 0 )
                std::cout << "value = " << value << std::endl;
            //        }
            //      }

            /*
                for (int i = 0; i < num_dofs_right; ++i) {
                  for (int j = 0; j < num_dofs_left; ++j) {
                    lv(RightAA::dof_index(i, 0), LeftAA::dof_index(j, 0)) +=
                      W * (
                        - C * dot(exact_sol_.A * LeftAA::grad_phi(j,q), RightAA::n(q)) * RightAA::phi(i,q)
                        + C * theta_ * dot(exact_sol_.A * RightAA::grad_phi(i,q), LeftAA::n(q)) * LeftAA::phi(j,q)
                        + gamma_ * alpha * LeftAA::phi(j,q) * RightAA::phi(i,q) * dot(LeftAA::n(q), RightAA::n(q))
                        );
                  }
                }
             */
        }

        // if (is_boundary) {
        //   std::cerr << "Boundary facet: (cell " << left_elem.get_cell().id()
        //             << ", facet = " << left_facet_number << "), x(0) = " << LeftAA::x(0) << ", left_n = "
        //             << LeftAA::n(0) << ", right_n = " << RightAA::n(0) << "\n";
        // }

        //    std::cerr << "IF1 = " << left_if_side << ", " << "IF2 = " << right_if_side << "\n";

        //    std::cerr << "A(l=" << left_elem.get_cell().index() << ", r=" << right_elem.get_cell().index() << ") = \n" << lv << "\n";
    }

  private:
    const LAD::VectorType* sol_;
    double theta_;
    double gamma_;
};

///
/// \brief Assembles jumps of discrete solution at element interfaces.
/// \details WeightedJumpTermAssembler is a AssemblyAssistant, consisting of \
///          a left-hand-sided and a right-hand-sided AssemblyAssistant. \
///          Three assembly assistants are used for the calculation:
///          i)   LeftAA:  Discrete function of one cell
///          ii)  RightAA: Discrete function of neighbouring cell
///          iii) HOAA:    Higher-Order function on one cell
///

class WeightedJumpTermAssembler : private LeftAA, private RightAA
{
  public:
    typedef StandardGlobalAssemblerWithJumps<double>::InterfaceSide InterfaceSide;

    WeightedJumpTermAssembler ( const VectorSpace<double>& space,
                                const double nu,
                                const double nu_vertical )
    : space_ho_ ( space ),
    nu_ ( nu ),
    nu_vertical_ ( nu_vertical )
    {
        // Crank Nicolson
        theta1_ = 0.5;
        theta2_ = 0.5;

        dT_ = -9999999999999.9;
    }

    void set_solP ( const LAD::VectorType* sol )
    {
        solP_ = sol;
    }

    void set_solP_prev ( const LAD::VectorType* sol )
    {
        solP_prev_ = sol;
    }

    void set_solD ( const LAD::VectorType* sol )
    {
        solD_ = sol;
    }

    void set_dT ( const double dT )
    {
        dT_ = dT;
    }

  protected:

    // For the CCI-3D application within MetStroem project

    double nu ( int quadrant, int derivative_direction ) const
    {
        static double visc [5][3] = {
            { this->nu_, this->nu_, this->nu_vertical_, }, // u-component
            { this->nu_, this->nu_, this->nu_vertical_, }, // v-component
            { this->nu_, this->nu_, this->nu_vertical_, }, // w-component
            { 0., 0., 0., }, // pressure -> not used
            { this->nu_, this->nu_, this->nu_vertical_, } // theta
        };

        assert ( quadrant >= 0 );
        assert ( quadrant < 5 );
        assert ( quadrant != 3 ); // no viscosity for pressure
        assert ( derivative_direction >= 0 );
        assert ( derivative_direction < 3 );

        return visc[quadrant][derivative_direction];
    }

    const double nu_;
    const double nu_vertical_;

    double theta1_; // linear part at t_i    -> implicit
    double theta2_; // linear part at t_i-1  -> explicit

    double dT_;

  public:

    void operator() ( const Element<double>& left_elem,
            const Element<double>& right_elem,
            const Quadrature<double>& left_quad,
            const Quadrature<double>& right_quad,
            int left_facet_number,
            int right_facet_number,
            InterfaceSide left_if_side,
            InterfaceSide right_if_side,
            double& value )
    {
        value = 0.;
        int debug_level = 0;
        const int t_var = DIM + 1;

        // Initialize third assembly assistant for higher-order test function
        // TODO: Should be placed in outer instance for performance
        AssemblyAssistant<DIM, double> HOAA;

        // Initialize the Assembly Assistants
        HOAA. initialize_for_facet ( Element<double>( space_ho_, left_elem.get_cell_index ( ) ),
                                     left_quad,
                                     left_facet_number );
        LeftAA::initialize_for_facet ( left_elem,
                                       left_quad,
                                       left_facet_number );
        RightAA::initialize_for_facet ( right_elem,
                                        right_quad,
                                        right_facet_number );

        // Set up local left and right gradients of discrete function and 
        // values of weighting

        assert ( left_quad.size ( ) == right_quad.size ( ) );

        assert ( LeftAA::num_vars ( ) == RightAA::num_vars ( ) );
        assert ( LeftAA::num_vars ( ) == HOAA.num_vars ( ) );

        const int num_vars = LeftAA::num_vars ( );

        std::vector<FunctionValues<Vec<DIM, double> > > left_gradients;
        std::vector<FunctionValues<Vec<DIM, double> > > right_gradients;

        std::vector<FunctionValues<Vec<DIM, double> > > left_gradients_prev;
        std::vector<FunctionValues<Vec<DIM, double> > > right_gradients_prev;

        std::vector<FunctionValues<double> > ho_sol;

        left_gradients.resize ( num_vars );
        right_gradients.resize ( num_vars );

        left_gradients_prev.resize ( num_vars );
        right_gradients_prev.resize ( num_vars );
        ho_sol.resize ( num_vars );

        for ( int var = 0; var < num_vars; ++var )
        {
            // set up gradients corresponding to solP
            LeftAA::evaluate_fe_function_gradients ( *solP_, var, left_gradients.at ( var ) );
            RightAA::evaluate_fe_function_gradients ( *solP_, var, right_gradients.at ( var ) );

            // set up gradients corresponding to solP_prev
            LeftAA::evaluate_fe_function_gradients ( *solP_prev_, var, left_gradients_prev.at ( var ) );
            RightAA::evaluate_fe_function_gradients ( *solP_prev_, var, right_gradients_prev.at ( var ) );

            // set up solution values corresponding to solD
            HOAA.evaluate_fe_function ( *solD_, var, ho_sol.at ( var ) );
        }

        const bool is_boundary = ( right_if_side == StandardGlobalAssemblerWithJumps<double>::INTERFACE_BOUNDARY );

        // Initialize Quadrature

        assert ( LeftAA::num_quadrature_points ( ) == RightAA::num_quadrature_points ( ) );
        assert ( LeftAA::num_quadrature_points ( ) == HOAA.num_quadrature_points ( ) );

        const int num_q = LeftAA::num_quadrature_points ( );

        // Quadrature loop

        for ( int q = 0; q < num_q; ++q )
        {
            const double dS =
                    ( is_boundary || left_if_side == StandardGlobalAssemblerWithJumps<double>::INTERFACE_SLAVE ) ?
                    std::abs ( LeftAA::ds ( q ) ) : std::abs ( RightAA::ds ( q ) );

            const double W = LeftAA::w ( q ) * dS;

            // Contribution -> velocity and temperature, but not the pressure

            int ignore_var = DIM; // not the pressure quadrant

            for ( int var = 0; var < num_vars; ++var )
            {
                if ( var != ignore_var )
                {
                    for ( int dim = 0; dim < DIM; ++dim )
                    {
                        // contribution of solP
                        value += W * dT_
                                * theta1_
                                * nu ( var, dim ) * ( left_gradients.at ( var )[q][dim] * LeftAA::n ( q )[dim]
                                + right_gradients.at ( var )[q][dim] * RightAA::n ( q )[dim] )
                                * ho_sol.at ( var )[q];

                        // contribution of solP_prev
                        value += W * dT_
                                * theta2_
                                * nu ( var, dim ) * ( left_gradients_prev.at ( var )[q][dim] * LeftAA::n ( q )[dim]
                                + right_gradients_prev.at ( var )[q][dim] * RightAA::n ( q )[dim] )
                                * ho_sol.at ( var )[q];
                    }
                }
            }
        }
        //     std::cout << "Left:  " << LeftAA::n(0)[0] << "\t" << LeftAA::n(0)[1] << "\t" << LeftAA::n(0)[2] << std::endl;
        //     std::cout << "Right: " << RightAA::n(0)[0] << "\t" << RightAA::n(0)[1] << "\t" << RightAA::n(0)[2] << std::endl;
        //     std::cout << std::endl;
    }

  private:
    const LAD::VectorType* solP_;
    const LAD::VectorType* solP_prev_;
    const LAD::VectorType* solD_;

    const VectorSpace<double>& space_ho_;

    double theta_;
    double gamma_;
};

#endif
