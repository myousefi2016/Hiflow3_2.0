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

#ifndef FLOW_TUTORIAL_H
#    define FLOW_TUTORIAL_H

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

int main ( int argc, char** argv );

class FlowTutorial : public NonlinearProblem<LAD>
{
  public:
    FlowTutorial ( const std::string& param_filename, const std::string& path_mesh );
    ~FlowTutorial ( );

    virtual void run ( );
    void ApplyFilter ( LAD::VectorType& u );
    virtual void EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F );
    virtual void EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF );

  private:

    const MPI_Comm& communicator ( )
    {
        return comm_;
    }

    int rank ( )
    {
        return rank_;
    }

    int num_partitions ( )
    {
        return num_partitions_;
    }

    // Read, refine and partition mesh.
    std::string path_mesh;
    void read_mesh ( );

    // Set up datastructures and read in some parameters.
    void prepare ( );

    // Set up boundary conditions
    void prepare_bc ( );

    // Solve nonlinear problem
    void solve ( );

    // Computate instationary solution by time-stepping method.
    void run_time_loop ( );

    // Visualize the solution in a file. In stationary mode, the
    // filename contains 'stationary', in instationary mode, it contains the
    // current time-step ts_.
    void visualize ( );

    // Adapt the mesh, before restarting computation.
    void adapt ( );

    // Helper functions for nonlinear solver
    // updates res_ with the residual
    void compute_residual ( const LAD::VectorType& u, LAD::VectorType* F );
    // residual computation in stationary mode
    void compute_stationary_residual ( const LAD::VectorType& u, LAD::VectorType* F );
    // residual computation in instationary mode
    void compute_instationary_residual ( const LAD::VectorType& u, LAD::VectorType* F );

    // updates matrix_ with the jacobian matrix
    void compute_jacobian ( const LAD::VectorType& u, LAD::MatrixType* DF );
    // jacobi matrix computation in stationary mode
    void compute_stationary_matrix ( const LAD::VectorType& u, LAD::MatrixType* DF );
    // jacobi matrix computation in instationary mode
    void compute_instationary_matrix ( const LAD::VectorType& u, LAD::MatrixType* DF );

    // Pressure filter: substracts the mean of the pressure from each
    // pressure dof in sol_ .
    void filter_pressure ( );

    // Linear algebra set up
    // void setup_linear_algebra();

    // MPI stuff
    MPI_Comm comm_;
    int rank_, num_partitions_;

    // Parameter data read in from file.
    PropertyTree params_;

    // parameter 'OutputPrefix': prefix for output files
    std::string simul_name_;
    // parameter 'FlowModel.Type': either "Channel" or "Cavity" .
    std::string flow_model_;

    // Time-stepping variables
    int ts_;
    double dt_;
    double alpha1_, alpha2_, alpha3_;

    // Flow model variables
    double Um_, H_, W_, rho_, nu_;

    // Flag for pressure filter -- parameter 'UsePressureFilter'
    bool use_pressure_filter_;

    // Meshes
    MeshPtr master_mesh_, mesh_;
    int refinement_level_;

    VectorSpace<double> space_;

    // Linear algebra objects
    Couplings<double> couplings_;
    CoupledMatrix<Scalar>* matrix_;
    CoupledVector<Scalar>* sol_, *prev_sol_, *res_, *rhs_;
    CVector pressure_correction_;

    DGGlobalAssembler<double> global_asm_;

    bool is_done_, solve_instationary_;

    std::vector<int> dirichlet_dofs_;
    std::vector<Scalar> dirichlet_values_;

    NonlinearSolver<LAD>* nls_;

    PreconditionerVanka<LAD> precond_;
};

//////////////// Boundary conditions ////////////////////////////////

struct ChannelFlowBC2d
{
    // Parameters:
    // var - variable
    // H - channel height
    // Um - maximum inflow
    // inflow_bdy - material number of inflow boundary
    // outflow_bdy - material number of outflow boundary

    ChannelFlowBC2d ( int var, double H, double Um, int inflow_bdy, int outflow_bdy )
    : var_ ( var ), H_ ( H ), Um_ ( Um ), inflow_bdy_ ( inflow_bdy ),
    outflow_bdy_ ( outflow_bdy )
    {
        assert ( var_ == 0 || var_ == 1 );
    }

    std::vector<double> evaluate ( const Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
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
            values.resize ( coords_on_face.size ( ) );

            // loop over points on the face
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
                    { // y-component
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
    const double H_;
    const double Um_;
    const int inflow_bdy_, outflow_bdy_;
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

struct LidCavityBC
{

    LidCavityBC ( int var, double U, int lid_bdy )
    : var_ ( var ), U_ ( U ), lid_bdy_ ( lid_bdy )
    {
    }

    std::vector<double> evaluate ( const Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        const int material_num = face.get_material_number ( );

        if ( material_num == lid_bdy_ && var_ == 0 )
        {
            return std::vector<double>( coords_on_face.size ( ), U_ );
        }
        else
        {
            return std::vector<double>( coords_on_face.size ( ), 0. );
        }
    }
    const int var_;
    const double U_;
    const int lid_bdy_;
};

//////////////// End boundary conditions /////////////////////////////

//////////////// Stationary assembler ////////////////////////////////

class StationaryFlowAssembler : private AssemblyAssistant<DIMENSION, double>
{
  public:

    StationaryFlowAssembler ( const CVector& solution, double nu, double rho )
    : solution_ ( solution ), nu_ ( nu ), inv_rho_ ( 1. / rho )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
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

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv )
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
                            wq * ( nu_ * dot ( grad_phi ( i, q, v_var ), grad_prev_vel_[v_var][q] ) )
                            * dJ;
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
        }
    }

  private:
    const CVector& solution_;
    double nu_, inv_rho_;
    FunctionValues<double> prev_vel_[DIMENSION];
    FunctionValues<double> pressure_k_;
    FunctionValues< Vec<DIMENSION, double> > grad_prev_vel_[DIMENSION];
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

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm )
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

        // indices i -> trial variable, j -> test variable
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
                                wq * (
                                // a0(\phi_i,\phi_j) = \int{ dot(\phi_j,\phi_i) }
                                phi ( j, q, var ) * phi ( i, q, var ) +

                                // a1(\phi_i, \phi_j) = alpha1 * \int{ \nu \grad{phi_j} :
                                // \grad{\phi_i} }
                                alpha1_ * nu_ * dot ( grad_phi ( j, q, var ), grad_phi ( i, q, var ) ) +

                                // c1(\phi_i,\phi_j) = \alpha1 * \int { (vel_ns*\grad{\phi_j})
                                // * \phi_i }
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
                                -wq * alpha2_ * inv_rho_ *
                                ( phi ( j, q, p_var ) * grad_phi ( i, q, v_var )[v_var] ) * dJ;
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
                                wq * inv_rho_ *
                                ( grad_phi ( j, q, u_var )[u_var] * phi ( i, q, q_var ) ) * dJ;
                    }
                }
            }
        } // end loop q
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv )
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

        // indices i -> trial variable, j -> test variable
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
                            wq * (
                            // l0(\phi_i) = \int {dot(vel_ns - vel_ts, \phi_i)}
                            ( vel_ns[v_var] - vel_ts[v_var] ) * phi ( i, q, v_var ) +

                            // l1(\phi_i) =
                            // \nu * \int{(\alpha1\grad{vel_ns} +
                            // \alpha3\grad{vel_ts}) * \grad{\phi_i}}
                            nu_ * ( dot ( ( alpha1_ * grad_vel_ns_[v_var][q] ) +
                                          ( alpha3_ * grad_vel_ts_[v_var][q] ), grad_phi ( i, q, v_var ) ) ) +

                            // l2(\phi_i) = \int{ (\alpha1*dot(vel_ns, \grad{\vel_ns}
                            //             + \alpha3*dot(vel_ts, \grad{\vel_ts}) * \phi_i}
                            ( alpha1_ * dot ( grad_vel_ns_[v_var][q], vel_ns ) +
                            alpha3_ * dot ( grad_vel_ts_[v_var][q], vel_ts ) ) *
                            phi ( i, q, v_var ) +

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
    // velocity at previous newton step
    FunctionValues<double> vel_ns_[DIMENSION];
    // velocity at previous timestep
    FunctionValues<double> vel_ts_[DIMENSION];
    // pressure at previous newton step
    FunctionValues<double> p_ns_;
    // gradient of velocity at previous newton step
    FunctionValues< Vec<DIMENSION, double> > grad_vel_ns_[DIMENSION];
    // gradient of velocity at previous timestep
    FunctionValues< Vec<DIMENSION, double> > grad_vel_ts_[DIMENSION];
};

#endif
