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

/// \brief This program benchmarks the assembly times for three
/// different Local Assemblers: Laplace equation and two variants of an
/// instationary Navier-Stokes equation. It measures the execution
/// times and computes the mean value and the standard deviation of the
/// execution times.

/// \author Staffan Ronnas

#ifndef HIFLOW_ASSEMBLY_BENCH_H
#    define HIFLOW_ASSEMBLY_BENCH_H

#    include "hiflow.h"

#    include <fstream>
#    include <iostream>
#    include <unistd.h>

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::la;
using namespace hiflow::doffem;

typedef LADescriptorCoupledD LAD;

// Input parameters.
const char* g_mesh_filename; // name of mesh                                (required input argument)
int g_fe_degree; // FE degree to use                            (option -p)
int g_refine_level; // mesh refinement level                       (option -l)
int g_dim; // mesh dimension                              (option -d)
int g_num_runs; // num repetitions of global assembly          (option -r)
int g_verbosity; // verbosity level                             (option -v)
int g_first_test, g_last_test; // range of tests to perform                   (options -b and -e)
bool g_write_matrix; // whether or not to write out test matrices   (option -m)

// HiFlow objects
SYSTEM g_sys;
MeshPtr g_mesh;
VectorSpace<double> g_space;
Couplings<double> g_couplings;
LAD::MatrixType g_matrix;
LAD::VectorType g_vec;
std::vector<int> g_diag_rows, g_diag_cols, g_off_diag_rows, g_off_diag_cols;

// Timing objects
std::vector<double> g_global_times;
std::vector< std::vector<double> > g_local_times;
int g_run;

// Global variables common to all tests
const PLATFORM APP_PLATFORM = CPU;
const IMPLEMENTATION APP_LINALG_IMPLEMENTATION = NAIVE;
const MATRIX_FORMAT APP_MATRIX_FORMAT = CSR;
const MPI_Comm WORLD_COMM = MPI_COMM_WORLD;

// Function declarations
void test_laplace_aa ( );
void test_navier_stokes_aa ( );
void test_navier_stokes_aa_alt ( );
void setup ( );
void compute_stats ( const std::vector<double>& values, double& mean, double& variance );
void print_stats ( );
void write_matrix ( );

//////////////// AssemblyAssistant Laplace ////////////////

template<int DIMENSION>
struct ExactSol
{
    double operator() ( const Vec<DIMENSION, double>& pt ) const;
    Vec<DIMENSION, double> eval_grad ( const Vec<DIMENSION, double>& pt ) const;
};

template<>
double ExactSol<2>::operator() ( const Vec<2, double>& pt ) const
{
    const double x = pt[0];
    const double y = pt[1];
    const double pi = M_PI;
    return std::cos ( 2. * pi * x ) * std::cos ( 2. * pi * y );
}

template<>
Vec<2, double> ExactSol<2>::eval_grad ( const Vec<2, double>& pt ) const
{
    const double pi = M_PI;
    Vec<2, double> grad;
    const double x = pt[0];
    const double y = pt[1];
    grad[0] = -2. * pi * std::sin ( 2. * pi * x ) * std::cos ( 2. * pi * y );
    grad[1] = -2. * pi * std::cos ( 2. * pi * x ) * std::sin ( 2. * pi * y );
    return grad;
}

template<>
double ExactSol<3>::operator() ( const Vec<3, double>& pt ) const
{
    const double x = pt[0];
    const double y = pt[1];
    const double z = pt[2];
    const double pi = M_PI;
    return std::cos ( 2. * pi * x ) * std::cos ( 2. * pi * y ) * std::cos ( 2. * pi * z );
}

template<>
Vec<3, double> ExactSol<3>::eval_grad ( const Vec<3, double>& pt ) const
{
    const double pi = M_PI;
    Vec<3, double> grad;
    const double x = pt[0];
    const double y = pt[1];
    const double z = pt[2];
    grad[0] = -2. * pi * std::sin ( 2. * pi * x ) * std::cos ( 2. * pi * y ) * std::cos ( 2. * pi * z );
    grad[1] = -2. * pi * std::cos ( 2. * pi * x ) * std::sin ( 2. * pi * y ) * std::cos ( 2. * pi * z );
    grad[2] = -2. * pi * std::cos ( 2. * pi * x ) * std::cos ( 2. * pi * y ) * std::sin ( 2. * pi * z );
    return grad;
}

template<int DIMENSION>
class LocalLaplaceAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:
    // necessary since we have template parameter DIMENSION
    using AssemblyAssistant<DIMENSION, double>::num_quadrature_points;
    using AssemblyAssistant<DIMENSION, double>::num_dofs;
    using AssemblyAssistant<DIMENSION, double>::detJ;
    using AssemblyAssistant<DIMENSION, double>::w;
    using AssemblyAssistant<DIMENSION, double>::phi;
    using AssemblyAssistant<DIMENSION, double>::grad_phi;
    using AssemblyAssistant<DIMENSION, double>::dof_index;
    using AssemblyAssistant<DIMENSION, double>::x;

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            typename AssemblyAssistant<DIMENSION, double>::LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        Timer local_timer;

        const int num_q = num_quadrature_points ( );
        const int total_dofs = this->num_dofs_total ( );
        const int num_local_dofs = num_dofs ( 0 );

        lm.Clear ( );
        lm.Resize ( total_dofs, total_dofs );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            for ( int i = 0; i < num_local_dofs; ++i )
            {
                for ( int j = 0; j < num_local_dofs; ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                            wq * dot ( grad_phi ( i, q ), grad_phi ( j, q ) ) * dJ;
                }
            }
        }

        local_timer.stop ( );
        g_local_times.at ( g_run ).at ( element.get_cell_index ( ) ) = local_timer.get_duration ( );
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            typename AssemblyAssistant<DIMENSION, double>::LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );
        const int num_q = num_quadrature_points ( );
        const int num_local_dofs = num_dofs ( 0 );
        const int total_dofs = this->num_dofs_total ( );

        lv.clear ( );
        lv.resize ( total_dofs, 0.0 );

        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            for ( int i = 0; i < num_local_dofs; ++i )
            {
                lv[dof_index ( i, 0 )] +=
                        wq * ( my_f ( x ( q ) ) * phi ( i, q ) ) * dJ;
            }
        }
    }

    // driving force

    double my_f ( Vec<DIMENSION, double> pt ) const
    {
        ExactSol<DIMENSION> sol;
        const double pi = M_PI;
        return 8. * pi * pi * sol ( pt );
    }
};

//////////////// AssemblyAssistant InstationaryNavierStokes ////////////////

template<int DIMENSION>
class InstationaryFlowAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:
    using AssemblyAssistant<DIMENSION, double>::num_quadrature_points;
    using AssemblyAssistant<DIMENSION, double>::num_dofs;
    using AssemblyAssistant<DIMENSION, double>::detJ;
    using AssemblyAssistant<DIMENSION, double>::w;
    using AssemblyAssistant<DIMENSION, double>::phi;
    using AssemblyAssistant<DIMENSION, double>::grad_phi;
    using AssemblyAssistant<DIMENSION, double>::dof_index;
    using AssemblyAssistant<DIMENSION, double>::x;
    using AssemblyAssistant<DIMENSION, double>::evaluate_fe_function;
    using AssemblyAssistant<DIMENSION, double>::evaluate_fe_function_gradients;

    InstationaryFlowAssembler ( double nu, double rho )
    : nu_ ( nu ), inv_rho_ ( 1. / rho )
    {
    }

    void set_newton_solution ( const LAD::VectorType* newton_sol )
    {
        newton_sol_ = newton_sol;
    }

    void set_prev_solution ( const LAD::VectorType* prev_sol )
    {
        prev_sol_ = prev_sol;
    }

    void set_timestep_parameters ( const std::vector<double>& alphas )
    {
        assert ( alphas.size ( ) == 5 );
        alpha1_ = alphas[0];
        alpha2_ = alphas[1];
        alpha3_ = alphas[2];
        alpha4_ = alphas[3];
        alpha5_ = alphas[4];
    }

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
    {
        AssemblyAssistant<DIMENSION, double>::initialize_for_element ( element, quadrature );

        // recompute previous solution values
        for ( int v = 0; v < DIMENSION; ++v )
        {
            prev_ns_vel_[v].clear ( );
            prev_ts_vel_[v].clear ( );
            grad_prev_ns_vel_[v].clear ( );
            grad_prev_ts_vel_[v].clear ( );

            evaluate_fe_function ( *newton_sol_, v, prev_ns_vel_[v] );
            evaluate_fe_function ( *prev_sol_, v, prev_ts_vel_[v] );
            evaluate_fe_function_gradients ( *newton_sol_, v, grad_prev_ns_vel_[v] );
            evaluate_fe_function_gradients ( *prev_sol_, v, grad_prev_ts_vel_[v] );
        }
        pressure_k_.clear ( );
        evaluate_fe_function ( *newton_sol_, DIMENSION, pressure_k_ );
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, typename AssemblyAssistant<DIMENSION, double>::LocalMatrix& lm )
    {
        initialize_for_element ( element, quadrature );

        Timer local_timer;

        const int num_q = num_quadrature_points ( );
        const int total_dofs = this->num_dofs_total ( );

        lm.Clear ( );
        lm.Resize ( total_dofs, total_dofs );

        // loop q
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // get previous solution in vector form
            Vec<DIMENSION, double> vel_k;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_k[var] = prev_ns_vel_[var][q];
            }

            // assemble a0(u,v) = \int {dot(u,v)}
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * ( phi ( j, q, u_var ) * phi ( i, q, u_var ) ) * dJ;
                    }
                }
            }

            // assemble a1(u,v) = \int \alpha_1 * {\grad(u) : \grad(v)}
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * alpha1_ * ( nu_ * dot ( grad_phi ( j, q, u_var ), grad_phi ( i, q, u_var ) ) ) * dJ;
                    }
                }
            }

            // assemble a2(u,v) = \int \alpha_1 * { (vel_k*\grad{u})*v }
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( u_var ); ++j )
                    {
                        lm ( dof_index ( i, u_var ), dof_index ( j, u_var ) ) +=
                                wq * alpha1_ * ( dot ( vel_k, grad_phi ( j, q, u_var ) ) * phi ( i, q, u_var ) ) * dJ;
                    }
                }
            }

            // assemble a3(u,v) = \int \alpha_1 * { (u\grad{u_k}*v }
            for ( int test_var = 0; test_var < DIMENSION; ++test_var )
            {
                for ( int trial_var = 0; trial_var < DIMENSION; ++trial_var )
                {
                    for ( int i = 0; i < num_dofs ( test_var ); ++i )
                    {
                        for ( int j = 0; j < num_dofs ( trial_var ); ++j )
                        {
                            lm ( dof_index ( i, test_var ), dof_index ( j, trial_var ) ) +=
                                    wq * alpha1_ * ( grad_prev_ns_vel_[test_var][q][trial_var] *
                                    phi ( j, q, trial_var ) *
                                    phi ( i, q, test_var ) ) * dJ;
                        }
                    }
                }
            }

            // assemble b(p, v) = - \alpha_2 * \int{p div{v}}
            const int p_var = DIMENSION;
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( p_var ); ++j )
                    {
                        lm ( dof_index ( i, v_var ), dof_index ( j, p_var ) ) +=
                                -wq * alpha2_ * ( inv_rho_ * phi ( j, q, p_var ) *
                                grad_phi ( i, q, v_var )[v_var] ) * dJ;
                    }
                }
            }

            // assemble bT(u, q) = - \int{q div(u)}
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
        local_timer.stop ( );
        g_local_times.at ( g_run ).at ( element.get_cell_index ( ) ) = local_timer.get_duration ( );
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, typename AssemblyAssistant<DIMENSION, double>::LocalVector& lv )
    {
        initialize_for_element ( element, quadrature );
        const int num_q = num_quadrature_points ( );
        const int total_dofs = this->num_dofs_total ( );

        lv.clear ( );
        lv.resize ( total_dofs, total_dofs );

        // loop over quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::abs ( detJ ( q ) );

            // get previous newton step solution in vector form
            Vec<DIMENSION, double> vel_k;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_k[var] = prev_ns_vel_[var][q];
            }

            // get previous time step solution in vector form
            Vec<DIMENSION, double> vel_n;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_n[var] = prev_ts_vel_[var][q];
            }

            // l0(v) = \int( dot(u_n - u_k, v))
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            wq * ( ( vel_n[v_var] - vel_k[v_var] ) * phi ( i, q, v_var ) ) * dJ;
                }
            }

            // l1n(v) = -\alpha_3 * \nu * \int( \grad{u_n} : \grad{v} )
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            -wq * ( nu_ * dot ( grad_phi ( i, q, v_var ), grad_prev_ts_vel_[v_var][q] ) ) * dJ;
                }
            }

            // l1k(v) = -\alpha_1 * \nu * \int( \grad{u_k} : \grad{v} )
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            -wq * alpha1_ * ( nu_ * dot ( grad_phi ( i, q, v_var ), grad_prev_ns_vel_[v_var][q] ) ) * dJ;
                }
            }

            // l2n(v) = -\alpha_3 * \int(u_n*\grad{u_n}*v)
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            -wq * ( dot ( grad_prev_ts_vel_[v_var][q], vel_n ) * phi ( i, q, v_var ) ) * dJ;
                }
            }

            // l2k(v) = -\alpha_1 * \int(u_k*\grad{u_k}*v)
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            -wq * alpha1_ * ( dot ( grad_prev_ns_vel_[v_var][q], vel_k ) * phi ( i, q, v_var ) ) * dJ;
                }
            }

            // l3(v) = 1/rho*\alpha_2*\int(p_k*div(v))
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            wq * alpha2_ * ( inv_rho_ * pressure_k_[q] * grad_phi ( i, q, v_var )[v_var] ) * dJ;
                }
            }

            // l4(q) = -\int(q * div(u_k))
            const int q_var = DIMENSION;
            double div_u_k = 0.;
            for ( int d = 0; d < DIMENSION; ++d )
            {
                div_u_k += grad_prev_ns_vel_[d][q][d];
            }

            for ( int i = 0; i < num_dofs ( q_var ); ++i )
            {
                lv[dof_index ( i, q_var )] +=
                        -wq * ( inv_rho_ * div_u_k * phi ( i, q, q_var ) ) * dJ;
            }
        }
    }

  private:
    const LAD::VectorType* newton_sol_;
    const LAD::VectorType* prev_sol_;

    double nu_, inv_rho_;
    double alpha1_, alpha2_, alpha3_, alpha4_, alpha5_;
    FunctionValues<double> prev_ns_vel_[DIMENSION];
    FunctionValues<double> prev_ts_vel_[DIMENSION];
    FunctionValues<double> pressure_k_;
    FunctionValues< Vec<DIMENSION, double> > grad_prev_ns_vel_[DIMENSION];
    FunctionValues< Vec<DIMENSION, double> > grad_prev_ts_vel_[DIMENSION];
};

//////////////// AssemblyAssistant InstationaryNavierStokes Alt. ////////////////

template<int DIMENSION>
class AltInstationaryFlowAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:
    using AssemblyAssistant<DIMENSION, double>::num_quadrature_points;
    using AssemblyAssistant<DIMENSION, double>::num_dofs;
    using AssemblyAssistant<DIMENSION, double>::detJ;
    using AssemblyAssistant<DIMENSION, double>::w;
    using AssemblyAssistant<DIMENSION, double>::phi;
    using AssemblyAssistant<DIMENSION, double>::grad_phi;
    using AssemblyAssistant<DIMENSION, double>::dof_index;
    using AssemblyAssistant<DIMENSION, double>::x;
    using AssemblyAssistant<DIMENSION, double>::evaluate_fe_function;
    using AssemblyAssistant<DIMENSION, double>::evaluate_fe_function_gradients;

    AltInstationaryFlowAssembler ( double nu, double rho )
    : nu_ ( nu ), inv_rho_ ( 1. / rho )
    {
    }

    void set_newton_solution ( const LAD::VectorType* newton_sol )
    {
        prev_newton_sol_ = newton_sol;
    }

    void set_time_solution ( const LAD::VectorType* time_sol )
    {
        prev_time_sol_ = time_sol;
    }

    void set_time_stepping_weights ( double alpha1, double alpha2, double alpha3 )
    {
        alpha1_ = alpha1;
        alpha2_ = alpha2;
        alpha3_ = alpha3;
    }

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
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
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, typename AssemblyAssistant<DIMENSION, double>::LocalMatrix& lm )
    {
        initialize_for_element ( element, quadrature );
        // indices j -> trial variable, i -> test variable
        // basis functions \phi -> velocity components, \eta -> pressure
        Timer local_timer;

        const int num_q = num_quadrature_points ( );
        const int total_dofs = this->num_dofs_total ( );

        lm.Clear ( );
        lm.Resize ( total_dofs, total_dofs );

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
                                wq * (
                                // a0(\phi_i,\phi_j) = \int{ dot(\phi_j,\phi_i) }
                                phi ( j, q, var ) * phi ( i, q, var ) +

                                // a1(\phi_i, \phi_j) = alpha1 * \int{ \nu \grad{phi_j} : \grad{\phi_i} }
                                alpha1_ * nu_ * dot ( grad_phi ( j, q, var ), grad_phi ( i, q, var ) ) +

                                // c1(\phi_i,\phi_j) = \alpha1 * \int { (vel_ns*\grad{\phi_j}) * \phi_i }
                                alpha1_ * dot ( vel_ns, grad_phi ( j, q, var ) ) * phi ( i, q, var )
                                ) * dJ;
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
        local_timer.stop ( );
        g_local_times.at ( g_run ).at ( element.get_cell_index ( ) ) = local_timer.get_duration ( );
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, typename AssemblyAssistant<DIMENSION, double>::LocalVector& lv )
    {
        initialize_for_element ( element, quadrature );
        // indices j -> trial variable, i -> test variable
        // basis functions \phi -> velocity components, \eta -> pressure

        const int num_q = num_quadrature_points ( );
        const int total_dofs = this->num_dofs_total ( );

        lv.clear ( );
        lv.resize ( total_dofs, 0. );

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
                            wq * (
                            // l0(\phi_i) = \int {dot(vel_ns - vel_ts, \phi_i)}
                            ( vel_ns[v_var] - vel_ts[v_var] ) * phi ( i, q, v_var ) +

                            // l1(\phi_i) =
                            // \nu * \int{(\alpha1\grad{vel_ns} + \alpha3\grad{vel_ts}) * \grad{\phi_i}}
                            nu_ * ( dot ( ( alpha1_ * grad_vel_ns_[v_var][q] ) +
                                          ( alpha3_ * grad_vel_ts_[v_var][q] ), grad_phi ( i, q, v_var ) ) ) +

                            // l2(\phi_i) = \int{ (\alpha1*dot(vel_ns, \grad{\vel_ns}
                            //                   + \alpha3*dot(vel_ts, \grad{\vel_ts}) * \phi_i}
                            ( alpha1_ * dot ( grad_vel_ns_[v_var][q], vel_ns ) +
                            alpha3_ * dot ( grad_vel_ts_[v_var][q], vel_ts ) ) * phi ( i, q, v_var ) +

                            // l3(\phi_i) = -\alpha2/rho*\int{p_ns * div(\phi_i)}
                            -( alpha2_ * inv_rho_ * p_ns_[q] * grad_phi ( i, q, v_var )[v_var] )
                            ) * dJ;
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
    const LAD::VectorType* prev_time_sol_;
    const LAD::VectorType* prev_newton_sol_;

    double alpha1_, alpha2_, alpha3_;

    double nu_, inv_rho_;
    FunctionValues<double> vel_ns_[DIMENSION]; // velocity at previous newton step
    FunctionValues<double> vel_ts_[DIMENSION]; // velocity at previous timestep
    FunctionValues<double> p_ns_; // pressure at previous newton step
    FunctionValues< Vec<DIMENSION, double> > grad_vel_ns_[DIMENSION]; // gradient of velocity at previous newton step
    FunctionValues< Vec<DIMENSION, double> > grad_vel_ts_[DIMENSION]; // gradient of velocity at previous timestep
};
#endif
