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

/// \author Jonas Kratzke, Katrin Mang

#ifndef BLOOD_FLOW_TUTORIAL_H
#    define BLOOD_FLOW_TUTORIAL_H

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

// Linear Algebra type renaming
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;

typedef std::vector<Scalar> Coord;

int main ( int argc, char** argv );

class BloodFlowTutorial : public NonlinearProblem<LAD>
{
  public:
    BloodFlowTutorial ( const std::string& param_filename );
    ~BloodFlowTutorial ( );

    virtual void run ( );

    // Helper functions for nonlinear solver
    // Updates res_ with the residual
    virtual void EvalFunc ( const CVector& newton_sol, CVector* res );

    // Updates matrix_ with the jacobian matrix
    virtual void EvalGrad ( const CVector& newton_sol, CMatrix* jac );

  private:
    // Read, refine and partition mesh.
    void read_mesh ( );

    // Set up datastructures and read in some parameters
    void prepare_parameters ( );

    // Prepare the FEM spaces for velocity and pressure
    void prepare_space ( );

    // Prepare the linear and the non-linear solver
    void prepare_solver ( );

    // Set up boundary conditions
    void prepare_bc ( );

    // Prepare interpolation of boundary flow rates
    void calculate_spline ( const int bdy_id );

    // Interpolation of boundary flow rates in time
    Scalar evaluate_spline ( const int bdy_id, Scalar current_time );

    // Compute instationary solution by a time-stepping method.
    void run_time_loop ( );

    // Solving step of the non-linear problem
    void solve ( );

    // Evaluate results
    void evaluate ( );
    Scalar evaluate_bdy_flowrate ( const int bdy_material );
    Scalar evaluate_bdy_pressure ( const int bdy_material );

    // Visualize the solution of the current time step
    void visualize ( );

    // MPI objects
    MPI_Comm comm_;
    int rank_, num_partitions_;

    // Parameter data from xml file.
    PropertyTree params_;

    // Simulation name used as prefix for output files
    std::string simul_name_;

    // FEM-mesh
    MeshPtr mesh_;

    // FEM-space
    VectorSpace<Scalar> space_;

    // Linear algebra objects
    Couplings<Scalar> couplings_;
    CMatrix* matrix_;
    CVector *sol_, *prev_sol_, *res_, *rhs_;

    // Assembler for linear system of equations
    StandardGlobalAssembler<Scalar> global_asm_;

    // Linear solver objects
    LinearSolver<LAD>* linsolver_;
#    ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> ilupp_;
    bool use_ilupp_;
#    endif

    // Non-linear solver objects
    NonlinearSolver<LAD>* nls_;
    EWForcing<LAD>* nls_forcing_;
    ArmijoDamping<LAD>* nls_damping_;

    // Time-stepping variables
    int ts_;
    Scalar dt_;
    Scalar alpha1_, alpha2_, stab_fac_;
    int visu_interval_;

    // Flow model variables
    Scalar rho_, nu_;

    // Dirichlet/Poiseuille boundary Information
    int wall_bdy_;

    std::vector<int> poiseuille_materials_;
    std::vector<Scalar> radius_, exponent_;
    std::vector<Coord> normal_, center_;
    Scalar period_;
    Scalar smooth_start_up_time_;
    std::vector<std::vector<Scalar> > timestamps_, flow_rate_, b_, c_, d_;
    std::vector<Scalar> smooth_initial_c_, smooth_initial_d_;

    std::vector<int> dirichlet_dofs_;
    std::vector<Scalar> dirichlet_values_;

    // Neumann/Windkessel boundary Information
    std::vector<int> windkessel_materials_;
    std::vector<Scalar> resistance_, compliance_, decay_;

    // Evaluation variables for Poiseuille boundary
    std::vector<Scalar> eval_poiseuille_flowrates_;
    std::vector<Scalar> eval_poiseuille_pressures_;

    // Evaluation variables for Windkessel boundary
    std::vector<Scalar> eval_windkessel_flowrates_;
    std::vector<Scalar> eval_windkessel_pressures_;

    // CSV Writer for evaluated results
    CSVWriter<Scalar> results_writer_;
};

/// Outer wall zero Dirichlet boundary condition

struct WallBC
{

    WallBC ( int wall_bdy )
    : wall_bdy_ ( wall_bdy )
    {
    }

    std::vector<Scalar> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<Scalar> values;
        if ( wall_bdy_ == face.get_material_number ( ) )
        {
            // Set zero values
            values.resize ( coords_on_face.size ( ), 0. );
        }
        return values;
    }

    const int wall_bdy_;
};

/// Dirichlet boundary condition for Poiseuille in- and outflow profiles

struct PoiseuilleBC
{

    PoiseuilleBC ( int var,
                   std::vector<int> material,
                   std::vector<Scalar> exponent,
                   std::vector<Scalar> radius,
                   std::vector<Coord> center,
                   std::vector<Coord> normal,
                   std::vector<Scalar> flow_rate )
    : var_ ( var ), material_ ( material ),
    exponent_ ( exponent ),
    radius_ ( radius ), center_ ( center ), normal_ ( normal ),
    flow_rate_ ( flow_rate )
    {
        assert ( var_ == 0 || var_ == 1 || var_ == 2 );
        assert ( material_.size ( ) == radius_.size ( ) );
        assert ( material_.size ( ) == normal_.size ( ) );
        assert ( material_.size ( ) == exponent_.size ( ) );
        assert ( material_.size ( ) == center_.size ( ) );
        assert ( material_.size ( ) == flow_rate_.size ( ) );
    }

    std::vector<Scalar> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {

        std::vector<Scalar> values;
        for ( int m = 0; m < material_.size ( ); ++m )
        {
            if ( material_[m] == face.get_material_number ( ) )
            {
                values.resize ( coords_on_face.size ( ) );
                // Loop over dof points on the face
                for ( int i = 0; i < coords_on_face.size ( ); ++i )
                {
                    // Get point
                    const Coord& pt = coords_on_face[i];
                    // Compute distance
                    Scalar r = 0;
                    for ( int d = 0; d < DIMENSION; ++d )
                    {
                        r += ( pt[d] - center_[m][d] )*( pt[d] - center_[m][d] );
                    }
                    r = std::sqrt ( r );
                    // Compute values
                    if ( r < radius_[m] )
                    {
                        values[i] = ( std::pow ( r / radius_[m], ( Scalar ) exponent_[m] ) - 1. );
                        assert ( values[i] >= -1 && values[i] <= 1 );
                        values[i] *= normal_[m][var_] * flow_rate_[m];
                    }
                    else
                    {
                        values[i] = 0;
                    }
                }
            }
        }
        return values;
    }

    const int var_;
    std::vector<int> material_;
    std::vector<Scalar> exponent_, radius_, flow_rate_;
    std::vector<Coord> normal_, center_;
};

/// Neumann BC for pressure: 2-element Windkessel model

class WindkesselBoundaryAssembler : private AssemblyAssistant<DIMENSION, Scalar>
{
  public:

    WindkesselBoundaryAssembler ( const std::vector<int>& material,
                                  const Scalar rho, const Scalar dt,
                                  const std::vector<Scalar>& prev_flowrate,
                                  const std::vector<Scalar>& prev_pressure,
                                  const std::vector<Scalar>& decay,
                                  const std::vector<Scalar>& resistance,
                                  const std::vector<Scalar>& compliance )
    : material_ ( material )
    {
        pressure_timestep_rhoinv_.resize ( material_.size ( ) );

        for ( int m = 0; m < material_.size ( ); ++m )
        {
            // P_n = R*dt*Q_(n-1) / (dt + C*R) + C*R*P_(n-1) / (dt + C*R)
            pressure_timestep_rhoinv_[m] = resistance[m] * dt * prev_flowrate[m] / ( dt + compliance[m] * resistance[m] )
                    + compliance[m] * resistance[m] * prev_pressure[m] / ( dt + compliance[m] * resistance[m] );
            // dt * P_n / rho = dt * ( (1 - decay) * P_n + decay * P_n-1 ) / rho
            pressure_timestep_rhoinv_[m] = dt * ( ( 1. - decay[m] ) * pressure_timestep_rhoinv_[m] + decay[m] * prev_pressure[m] ) / rho;
        }
    }

    void operator() ( const Element<Scalar>& element,
            const int facet_number,
            const Quadrature<Scalar>& quadrature,
            LocalVector& lv )
    {
        // Procedure to get the facet entity
        mesh::IncidentEntityIterator facet = element.get_cell ( ).begin_incident ( DIMENSION - 1 );
        for ( int i = 0; i < facet_number; ++i, ++facet )
        {
        }

        // Iterate Windkessel materials (boundaries)
        for ( int m = 0; m < material_.size ( ); ++m )
        {
            // Check if the facet has the current material number
            if ( material_[m] == facet->get_material_number ( ) )
            {

                // Initialize the quadrature for integration over the facet
                AssemblyAssistant<DIMENSION, Scalar>::initialize_for_facet ( element, quadrature, facet_number );

                // Loop over quadrature points
                const int num_q = num_quadrature_points ( );
                for ( int q = 0; q < num_q; ++q )
                {
                    const Scalar wq_dsurf = w ( q ) * ds ( q );
                    const Vec<DIMENSION, Scalar> normal = n ( q );

                    for ( int test_var = 0; test_var < DIMENSION; ++test_var )
                    {
                        for ( int i = 0; i < num_dofs ( test_var ); ++i )
                        {
                            lv[dof_index ( i, test_var )] +=
                                    wq_dsurf * pressure_timestep_rhoinv_[m] * phi ( i, q, test_var ) * normal[test_var];
                        }
                    }
                }
                return;
            }
        }
    }
  private:
    const std::vector<int>& material_;
    std::vector<Scalar> pressure_timestep_rhoinv_;
};

/// Instationary incompressible Navier-Stokes Assembler with partial convection stabilization

class InstationaryFlowAssembler : private AssemblyAssistant<DIMENSION, Scalar>
{
  public:

    InstationaryFlowAssembler ( const Scalar nu, const Scalar rho,
                                const CVector& newton_sol, const CVector& time_sol,
                                const Scalar alpha1, const Scalar dt, const Scalar alpha2,
                                const Scalar stab_fac )
    : nu_ ( nu ), inv_rho_ ( 1. / rho ),
    prev_newton_sol_ ( newton_sol ),
    prev_time_sol_ ( time_sol ),
    alpha1_ ( alpha1 ), dt_ ( dt ), alpha2_ ( alpha2 ),
    stab_fac_ ( stab_fac )
    {
    }

    void operator() ( const Element<Scalar>& element, const Quadrature<Scalar>& quadrature, LocalMatrix& lm )
    {
        AssemblyAssistant<DIMENSION, Scalar>::initialize_for_element ( element, quadrature );

        // Recompute previous solution values
        for ( int d = 0; d < DIMENSION; ++d )
        {
            vel_ns_[d].clear ( );
            grad_vel_ns_[d].clear ( );

            evaluate_fe_function ( prev_newton_sol_, d, vel_ns_[d] );
            evaluate_fe_function_gradients ( prev_newton_sol_, d, grad_vel_ns_[d] );
        }

        // Indices j -> trial variable, i -> test variable
        // Basis functions \phi -> velocity components, \eta -> pressure

        const int num_q = num_quadrature_points ( );

        // Loop quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const Scalar w_dJ_q = w ( q ) * std::abs ( detJ ( q ) );
            const Scalar h = std::pow ( std::abs ( detJ ( q ) ), 1. / ( Scalar ) DIMENSION );

            // Vector with velocity of previous Newton point
            Vec<DIMENSION, Scalar> vel_ns;
            std::vector< Vec<DIMENSION, Scalar> > grad_vel_ns;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_ns[var] = vel_ns_[var][q];
                grad_vel_ns.push_back ( grad_vel_ns_[var][q] );
            }
            const Scalar delta_stab_ns = stab_fac_ * h * h / ( 6. * nu_ + h * norm ( vel_ns ) );

            // Velocity terms symmetric in test and trial variables
            for ( int var = 0; var < DIMENSION; ++var )
            {
                const int n_dofs = num_dofs ( var );
                for ( int i = 0; i < n_dofs; ++i )
                {
                    for ( int j = 0; j < n_dofs; ++j )
                    {
                        lm ( dof_index ( i, var ), dof_index ( j, var ) ) +=
                                w_dJ_q *
                                // a0(\phi_i,\phi_j) = \int{ dot(\phi_j,\phi_i) }
                                ( ( phi ( j, q, var ) +

                                // c1(\phi_i,\phi_j) = \alpha1 * \int { (vel_ns*\grad{\phi_j}) * \phi_i }
                                alpha1_ * dot ( vel_ns, grad_phi ( j, q, var ) )
                                ) * phi ( i, q, var ) +

                                // c1_stab(\phi_i,\phi_j) = delta_stab * \alpha1 * \int { (vel_ns*\grad{\phi_j}) * \vel_ns\dot\grad{\phi_i} }
                                delta_stab_ns * alpha1_ * dot ( vel_ns, grad_phi ( j, q, var ) ) * dot ( vel_ns, grad_phi ( i, q, var ) ) +

                                // a1(\phi_i, \phi_j) = alpha1 * \int{ \nu \grad{phi_j} : \grad{\phi_i} }
                                alpha1_ * nu_ * dot ( grad_phi ( j, q, var ), grad_phi ( i, q, var ) ) );
                    }
                }
            }

            // c2(\phi_i,\phi_j) = \int \alpha1 * { (\phi_j*\grad{vel_ns}*\phi_i }
            // c2_stab(\phi_i,\phi_j) = delta_stab * \int \alpha1 * { (\phi_j*\grad{vel_ns}*\vel_ns\dot\grad{\phi_i} }
            // c3_stab(\phi_i,\phi_j) = der_delta_stab * \int \alpha1 * { (vel_ns*\grad{vel_ns}*\vel_ns\dot\grad{\phi_i} }
            // c4_stab(\phi_i,\phi_j) = delta_stab * \int \alpha1 * { (vel_ns*\grad{vel_ns}*\phi_j\dot\grad{\phi_i} }
            const Scalar der_stab_factor = -stab_fac_ * h * h * h / ( ( 6. * nu_ + h * norm ( vel_ns ) )*( 6. * nu_ + h * norm ( vel_ns ) ) );
            for ( int test_var = 0; test_var < DIMENSION; ++test_var )
            {
                for ( int trial_var = 0; trial_var < DIMENSION; ++trial_var )
                {
                    for ( int i = 0; i < num_dofs ( test_var ); ++i )
                    {
                        for ( int j = 0; j < num_dofs ( trial_var ); ++j )
                        {
                            lm ( dof_index ( i, test_var ), dof_index ( j, trial_var ) ) +=
                                    w_dJ_q * alpha1_ * (
                                    grad_vel_ns[test_var][trial_var] * phi ( j, q, trial_var ) *
                                    phi ( i, q, test_var ) +
                                    delta_stab_ns *
                                    grad_vel_ns[test_var][trial_var] * phi ( j, q, trial_var ) *
                                    dot ( vel_ns, grad_phi ( i, q, test_var ) ) +
                                    der_stab_factor * std::abs ( phi ( j, q, trial_var ) ) *
                                    dot ( vel_ns, grad_vel_ns[test_var] ) *
                                    dot ( vel_ns, grad_phi ( i, q, test_var ) ) +
                                    delta_stab_ns *
                                    dot ( vel_ns, grad_vel_ns[test_var] ) *
                                    phi ( j, q, trial_var ) * grad_phi ( i, q, test_var )[trial_var]
                                    );
                        }
                    }
                }
            }

            const int p_var = DIMENSION;
            for ( int u_var = 0; u_var < DIMENSION; ++u_var )
            {
                for ( int i = 0; i < num_dofs ( u_var ); ++i )
                {
                    for ( int j = 0; j < num_dofs ( p_var ); ++j )
                    {
                        // b(\phi_i, \eta_j) = -dt / rho * \int{\eta_j div{\phi_i}}
                        lm ( dof_index ( i, u_var ), dof_index ( j, p_var ) ) +=
                                -w_dJ_q * dt_ * inv_rho_ * ( phi ( j, q, p_var ) * grad_phi ( i, q, u_var )[u_var] );
                        // bT(\eta_j, \phi_i) = 1/rho \int{ \eta_j * div{\phi_i}}
                        lm ( dof_index ( j, p_var ), dof_index ( i, u_var ) ) +=
                                w_dJ_q * inv_rho_ * ( grad_phi ( i, q, u_var )[u_var] * phi ( j, q, p_var ) );
                    }
                }
            }
        } // End loop over quadrature points
    }

    void operator() ( const Element<Scalar>& element, const Quadrature<Scalar>& quadrature, LocalVector& lv )
    {
        AssemblyAssistant<DIMENSION, Scalar>::initialize_for_element ( element, quadrature );

        // Recompute previous solution values
        for ( int d = 0; d < DIMENSION; ++d )
        {
            vel_ns_[d].clear ( );
            vel_ts_[d].clear ( );
            grad_vel_ns_[d].clear ( );
            grad_vel_ts_[d].clear ( );

            evaluate_fe_function ( prev_newton_sol_, d, vel_ns_[d] );
            evaluate_fe_function ( prev_time_sol_, d, vel_ts_[d] );
            evaluate_fe_function_gradients ( prev_newton_sol_, d, grad_vel_ns_[d] );
            evaluate_fe_function_gradients ( prev_time_sol_, d, grad_vel_ts_[d] );
        }

        // Recompute pressure
        p_ns_.clear ( );
        evaluate_fe_function ( prev_newton_sol_, DIMENSION, p_ns_ );

        // Indices j -> trial variable, i -> test variable
        // Basis functions \phi -> velocity components, \eta -> pressure

        const int num_q = num_quadrature_points ( );

        // Loop quadrature points
        for ( int q = 0; q < num_q; ++q )
        {
            const Scalar w_dJ_q = w ( q ) * std::abs ( detJ ( q ) );
            const Scalar h = std::pow ( std::abs ( detJ ( q ) ), 1. / ( Scalar ) DIMENSION );

            // Get previous Newton and time step solutions in vector form
            Scalar p_ns = p_ns_[q];
            Vec<DIMENSION, Scalar> vel_ts, vel_ns;
            std::vector< Vec<DIMENSION, Scalar> > grad_vel_ts, grad_vel_ns;
            for ( int var = 0; var < DIMENSION; ++var )
            {
                vel_ns[var] = vel_ns_[var][q];
                vel_ts[var] = vel_ts_[var][q];
                grad_vel_ts.push_back ( grad_vel_ts_[var][q] );
                grad_vel_ns.push_back ( grad_vel_ns_[var][q] );
            }

            // Stabilization parameters
            const Scalar delta_stab_ns = stab_fac_ * h * h / ( 6. * nu_ + h * norm ( vel_ns ) );
            const Scalar delta_stab_ts = stab_fac_ * h * h / ( 6. * nu_ + h * norm ( vel_ts ) );

            // Residual without incompressibility
            for ( int v_var = 0; v_var < DIMENSION; ++v_var )
            {
                for ( int i = 0; i < num_dofs ( v_var ); ++i )
                {
                    lv[dof_index ( i, v_var )] +=
                            w_dJ_q *
                            // l0(\phi_i) = \int {dot(vel_ns - vel_ts, \phi_i)}
                            ( ( vel_ns[v_var] - vel_ts[v_var] ) * phi ( i, q, v_var ) +

                            // l1(\phi_i) =
                            // \nu * \int{(\alpha1\grad{vel_ns} + \alpha2\grad{vel_ts}) * \grad{\phi_i}}
                            nu_ * ( dot ( ( alpha1_ * grad_vel_ns[v_var] ) +
                                          ( alpha2_ * grad_vel_ts[v_var] ), grad_phi ( i, q, v_var ) ) ) +

                            // l2(\phi_i) = \int{ (\alpha1*dot(vel_ns, \grad{\vel_ns}
                            //                   + \alpha2*dot(vel_ts, \grad{\vel_ts}) * \phi_i}
                            ( alpha1_ * dot ( grad_vel_ns[v_var], vel_ns ) +
                            alpha2_ * dot ( grad_vel_ts[v_var], vel_ts ) ) * phi ( i, q, v_var ) +

                            // l2_stab(\phi_i) = \int{ delta_stab_ns * (\alpha1*dot(vel_ns, \grad{\vel_ns} * \vel_ns\dot\grad{\phi_i}
                            //                        + delta_stab_ts * \alpha2*dot(vel_ts, \grad{\vel_ts} * \vel_ts\dot\grad{\phi_i}}
                            delta_stab_ns * alpha1_ * dot ( grad_vel_ns[v_var], vel_ns ) * dot ( vel_ns, grad_phi ( i, q, v_var ) ) +
                            delta_stab_ts * alpha2_ * dot ( grad_vel_ts[v_var], vel_ts ) * dot ( vel_ts, grad_phi ( i, q, v_var ) ) +

                            // l3(\phi_i) = -dt/rho*\int{p_ns * div(\phi_i)}
                            -( dt_ * inv_rho_ * p_ns * grad_phi ( i, q, v_var )[v_var] ) );
                }
            }

            // Incompressibility term
            const int q_var = DIMENSION;
            Scalar div_u_k = 0.;
            for ( int d = 0; d < DIMENSION; ++d )
            {
                div_u_k += grad_vel_ns_[d][q][d];
            }

            // l4(\eta_i) = 1/rho*\int{\eta_i * div(vel_ns)}
            for ( int i = 0; i < num_dofs ( q_var ); ++i )
            {
                lv[dof_index ( i, q_var )] +=
                        w_dJ_q * ( inv_rho_ * div_u_k * phi ( i, q, q_var ) );
            }
        } // End quadrature loop
    }

  private:
    const CVector& prev_time_sol_;
    const CVector& prev_newton_sol_;

    const Scalar alpha1_, dt_, alpha2_;

    const Scalar nu_, inv_rho_, stab_fac_;
    // Velocity at previous Newton step
    FunctionValues<Scalar> vel_ns_[DIMENSION];
    // Velocity at previous timestep
    FunctionValues<Scalar> vel_ts_[DIMENSION];
    // Pressure at previous Newton step
    FunctionValues<Scalar> p_ns_;
    // Gradient of velocity at previous Newton step
    FunctionValues< Vec<DIMENSION, Scalar> > grad_vel_ns_[DIMENSION];
    // Gradient of velocity at previous timestep
    FunctionValues< Vec<DIMENSION, Scalar> > grad_vel_ts_[DIMENSION];
};

/// Flowrate evaluation at in- and outflow boundaries

class FlowrateIntegrator : private AssemblyAssistant<DIMENSION, Scalar>
{
  public:

    FlowrateIntegrator ( const CVector& sol, const int material_num )
    : sol_ ( sol ), mat_num_ ( material_num )
    {
    }

    void operator() ( const Element<Scalar>& element, int facet_num, const Quadrature<Scalar>& quadrature, std::vector<Scalar>& flux )
    {
        // Check the material number of the facet
        mesh::IncidentEntityIterator fac_it = element.get_cell ( ).begin_incident ( DIMENSION - 1 );
        for ( int i = 0; i < facet_num; ++i, ++fac_it )
        {
        }
        if ( fac_it->get_material_number ( ) != mat_num_ ) return;
        AssemblyAssistant<DIMENSION, Scalar>::initialize_for_facet ( element, quadrature, facet_num );

        for ( int d = 0; d < DIMENSION; ++d )
        {
            vel_[d].clear ( );
            evaluate_fe_function ( sol_, d, vel_[d] );
        }

        for ( int q = 0; q < num_quadrature_points ( ); ++q )
        {
            const Vec<DIMENSION, Scalar> normal = n ( q );
            const Scalar w_ds_q = w ( q ) * ds ( q );
            for ( int d = 0; d < DIMENSION; ++d )
            {
                flux[facet_num] += w_ds_q * normal[d] * vel_[d][q];
            }
        }
    }

  private:
    // Solution vector
    const CVector& sol_;
    // Velocity part of the solution vector
    FunctionValues<Scalar> vel_[DIMENSION];
    const int mat_num_;
};

/// Area evaluation at in- and outflow boundaries

class BoundaryAreaIntegrator : private AssemblyAssistant<DIMENSION, Scalar>
{
  public:

    BoundaryAreaIntegrator ( const int material_num )
    : mat_num_ ( material_num )
    {
    }

    void operator() ( const Element<Scalar>& element, int facet_num, const Quadrature<Scalar>& quadrature, std::vector<Scalar>& value )
    {
        mesh::IncidentEntityIterator fac_it = element.get_cell ( ).begin_incident ( DIMENSION - 1 );
        for ( int i = 0; i < facet_num; ++i, ++fac_it )
        {
        }
        if ( fac_it->get_material_number ( ) != mat_num_ ) return;
        AssemblyAssistant<DIMENSION, Scalar>::initialize_for_facet ( element, quadrature, facet_num );

        for ( int q = 0; q < num_quadrature_points ( ); ++q )
        {
            value[facet_num] += w ( q ) * ds ( q );
        }
    }

  private:
    // Material number of boundary
    const int mat_num_;
};

/// Pressure evaluation at in- and outflow boundaries

class BoundaryValueIntegrator : private AssemblyAssistant<DIMENSION, Scalar>
{
  public:

    BoundaryValueIntegrator ( const int var, const CVector& sol, const int material_num )
    : var_ ( var ), sol_ ( sol ), mat_num_ ( material_num )
    {
    }

    void operator() ( const Element<Scalar>& element, int facet_num, const Quadrature<Scalar>& quadrature, std::vector<Scalar>& value )
    {

        mesh::IncidentEntityIterator fac_it = element.get_cell ( ).begin_incident ( DIMENSION - 1 );
        for ( int i = 0; i < facet_num; ++i, ++fac_it )
        {
        }
        if ( fac_it->get_material_number ( ) != mat_num_ ) return;
        AssemblyAssistant<DIMENSION, Scalar>::initialize_for_facet ( element, quadrature, facet_num );

        sol_loc_.clear ( );
        evaluate_fe_function ( sol_, var_, sol_loc_ );

        for ( int q = 0; q < num_quadrature_points ( ); ++q )
        {
            value[facet_num] += w ( q ) * sol_loc_[q] * ds ( q );
        }
    }

  private:
    const int var_, mat_num_;
    // Solution vector
    const CVector& sol_;
    // Local part of the solution vector
    FunctionValues<Scalar> sol_loc_;
};

/// Visualization function for the wall shear stress magnitude

template<class LAD>
class EvalWSSmagn
{
    /*
      // Type of function for evaluation.
      typedef boost::function3<void,
          // Cell
          const mesh::Entity&,
          // Reference coordinates
          const std::vector<Scalar>&,
          // Values of function at the points
          std::vector<Scalar>&
          > EvalFunction;
     */
  public:

    EvalWSSmagn ( const VectorSpace<Scalar>& space,
                  const CVector& fun,
                  const int material,
                  const Scalar dyn_visc )
    : space_ ( space ), fun_ ( fun ), material_ ( material ), dyn_visc_ ( dyn_visc )
    {
    }

    std::vector<Scalar> compute_facet_normal ( const int dim,
                                               const std::vector<Scalar>& f_coords,
                                               const CellTransformation<Scalar>& ct ) const
    {
        // Compute outer the normal of the face
        std::vector<Scalar> dir ( dim * ( dim - 1 ), 0. );
        Scalar facet_size = 0;
        for ( int i = 0; i < dim - 1; ++i )
        {
            Scalar dir_size = 0;
            for ( int d = 0; d < dim; ++d )
            {
                dir[i * dim + d] = f_coords[( i + 1 ) * dim + d] - f_coords[d];
                dir_size += dir[i * dim + d] * dir[i * dim + d];
            }
            facet_size += std::sqrt ( dir_size );
        }
        facet_size /= dim - 1;
        std::vector<Scalar> n = normal ( dir );
        std::vector<Scalar> mid ( dim, 0. );
        std::vector<Scalar> mid_ref ( dim, 0. );
        for ( int d1 = 0; d1 < dim; ++d1 )
        {
            for ( int d2 = 0; d2 < dim; ++d2 )
            {
                mid[d1] += f_coords[d2 * dim + d1];
            }
            mid[d1] = mid[d1] / dim - n[d1] * facet_size * 1.e-5;
        }
        bool point_inside = ct.contains_physical_point ( mid, &mid_ref );
        if ( !point_inside )
        {
            for ( int d = 0; d < dim; d++ )
            {
                n[d] *= -1;
            }
        }
        return n;
    }

    void operator() ( const Entity& cell,
            const std::vector<Scalar>& ref_coords,
            std::vector<Scalar>& values ) const
    {
        const int gdim = space_.mesh ( ).gdim ( );
        const int num_points = ref_coords.size ( ) / gdim;

        values.clear ( );
        values.resize ( num_points, 0. );

        // Check if it is a boundary cell
        bool is_bdy_cell = false;
        for ( IncidentEntityIterator f_it = cell.begin_incident ( gdim - 1 );
              f_it != cell.end_incident ( gdim - 1 ); ++f_it )
        {
            if ( material_ == f_it->get_material_number ( ) )
            {
                is_bdy_cell = true;
            }
        }
        if ( !is_bdy_cell ) return;

        // Compute NSEstressTensor
        std::vector<Mat<DIMENSION, DIMENSION, Scalar> > stress_tensor ( num_points );

        for ( int row = 0; row < DIMENSION; ++row )
        {
            // Get global dof ids
            std::vector<int> global_dof_ids;
            space_.GetDofIndices ( row, cell, &global_dof_ids );
            const int num_dofs = global_dof_ids.size ( );

            // Vector of shape function gradient values
            std::vector< std::vector< std::vector<Scalar> > > shape_fun_grad ( gdim, std::vector< std::vector<Scalar> >( num_points, std::vector<Scalar>( num_dofs, 1.e13 ) ) );

            // Compute the shape function gradient values
            std::vector<Scalar> pt ( gdim, 0. );
            // Collect points to compute the Jacobian of the cell transformation later on
            // Quadrature points on evaluation element
            std::vector< Vec<DIMENSION, Scalar> > points;
            int k = 0;
            for ( int i = 0; i < num_points; ++i )
            {
                for ( int c = 0; c < gdim; ++c )
                {
                    pt[c] = ref_coords[k++];
                }
                space_.fe_manager ( ).get_fe_on_cell ( cell.index ( ), row )->N_x ( pt, shape_fun_grad[0][i] );
                space_.fe_manager ( ).get_fe_on_cell ( cell.index ( ), row )->N_y ( pt, shape_fun_grad[1][i] );
                space_.fe_manager ( ).get_fe_on_cell ( cell.index ( ), row )->N_z ( pt, shape_fun_grad[2][i] );
                Vec<DIMENSION, Scalar> Vpt ( pt );
                points.push_back ( Vpt );
            }
            // Get dof values
            std::vector<Scalar> dof_values ( global_dof_ids.size ( ), 1.e25 );
            fun_.GetValues ( &global_dof_ids[0], num_dofs, &dof_values[0] );

            // Compute inverse of the transposed Jacobian
            // Jacobian matrix at evaluation points
            FunctionValues< Mat<DIMENSION, DIMENSION, Scalar> > J, JinvT;
            const doffem::CellTransformation<Scalar>& ct = space_.GetCellTransformation ( cell );
            J.compute ( points, EvalPhysicalJacobian<DIMENSION, Scalar>( ct ) );
            JinvT.compute ( J, EvalInvTranspose<DIMENSION, Scalar>( ) );
            // Compute the velocity gradients
            for ( int i = 0; i < num_points; ++i )
            {
                // Compute the velocity gradients on the reference cell
                // by multiplying the shape function gradients with the dof values
                std::vector<Scalar> ref_velo_grad ( gdim, 0. );
                for ( int d = 0; d < gdim; ++d )
                {
                    for ( int j = 0; j < num_dofs; ++j )
                    {
                        ref_velo_grad[d] += shape_fun_grad.at ( d ).at ( i ).at ( j ) * dof_values.at ( j );
                    }
                }
                // Map the gradients from the reference to the physical cell
                for ( int col = 0; col < DIMENSION; ++col )
                {
                    for ( int d = 0; d < gdim; ++d )
                    {
                        stress_tensor[i]( row, col ) += JinvT[i]( col, d ) * ref_velo_grad[d];
                        stress_tensor[i]( col, row ) += JinvT[i]( col, d ) * ref_velo_grad[d];
                    }
                }
            }
        }
        for ( int i = 0; i < num_points; ++i )
        {
            // Consider only shear stresses
            for ( int d = 0; d < gdim; ++d )
            {
                stress_tensor[i]( d, d ) = 0.;
            }
            stress_tensor[i] *= dyn_visc_;
        }

        // For each point determine the facet it lies in
        // and compute the normal of the facet
        for ( int i = 0; i < num_points; ++i )
        {
            const CellTransformation<Scalar>& ct = space_.GetCellTransformation ( cell );
            const std::vector<Scalar> pt_ref ( ref_coords.begin ( ) + i*gdim, ref_coords.begin ( )+( i + 1 ) * gdim );
            std::vector<Scalar> pt_phys ( gdim );
            pt_phys[0] = ct.x ( pt_ref );
            pt_phys[1] = ct.y ( pt_ref );
            if ( gdim == 3 )
                pt_phys[2] = ct.z ( pt_ref );
            // Check if point is on the boundary
            int facet_count = 0;
            for ( IncidentEntityIterator f_it = cell.begin_incident ( gdim - 1 );
                  f_it != cell.end_incident ( gdim - 1 ); ++f_it )
            {
                if ( material_ == f_it->get_material_number ( ) )
                {
                    std::vector<Scalar> f_coords;
                    f_it->get_coordinates ( f_coords );
                    const std::vector<Scalar> n = compute_facet_normal ( gdim, f_coords, ct );
                    const std::vector<Scalar> p_a ( pt_phys );
                    const std::vector<Scalar> origin ( f_coords.begin ( ), f_coords.begin ( ) + gdim );
                    const std::vector<Scalar> o ( origin );
                    Scalar distance = std::abs ( Dot ( p_a, n ) - Dot ( o, n ) );
                    if ( distance < 1.e-14 )
                    {
                        Scalar norm = 0;
                        for ( int row = 0; row < gdim; ++row )
                        {
                            Scalar force = 0;
                            for ( int col = 0; col < gdim; ++col )
                            {
                                force += stress_tensor[i]( row, col ) * n[col];
                            }
                            norm += force*force;
                        }
                        values[i] += std::sqrt ( norm );
                        ++facet_count;
                    }
                }
            }
            if ( facet_count > 1 ) values[i] /= ( Scalar ) facet_count;
        }
    }

  private:
    const VectorSpace<Scalar>& space_;
    const CVector& fun_;
    const int material_;
    const Scalar dyn_visc_;
};

#endif
