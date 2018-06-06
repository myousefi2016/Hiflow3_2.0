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

/// \author Simon Gawlok, Jonathan Schwegler, Carmen Straub

#ifndef _INSTATIONARY_CONVDIFF_TUTORIAL_H_
#    define _INSTATIONARY_CONVDIFF_TUTORIAL_H_

#    define COMPUTE_HESSIAN

//#define RUNGE_KUTTA

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>
#    include <sstream>
#    include <limits>
#    include <algorithm>
#    include <iostream>
#    include <fstream>
#    include <boost/function.hpp>

#    include "hiflow.h"

// All names are imported for simplicity
using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

// Dimension of the problem
const int DIM = 2;

// Shorten some datatypes with typedefs
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;

typedef std::vector<double> Coord;

// Linear Algebra Parameters.
const PLATFORM la_platform = CPU;
const IMPLEMENTATION la_impl = NAIVE;
const MATRIX_FORMAT la_matrix_format = CSR;

int main ( int argc, char** argv );

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Stiffness-Matrix Assembler.

class StiffnessAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    StiffnessAssembler ( )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );
};

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Mass-Matrix Assembler.

class MassAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    MassAssembler ( )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );
};

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Convection-Diffusion Assembler.

class ConvectionAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    ConvectionAssembler ( )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );

    std::vector<double> get_convective_velocity ( int q );
};

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Source Assembler.

class SourceAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    SourceAssembler ( )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv );

    double get_source_term ( int q );
};

// Source function for visualization

class EvalSourceFunction
{
  public:

    EvalSourceFunction ( const VectorSpace<double>& space, double time )
    : space_ ( space ), time_ ( time )
    {
    }

    void operator() ( const Entity& cell,
            const std::vector<double>& ref_coords,
            std::vector<double>& values ) const
    {
        double s = 0.0;
        const int num_points = ref_coords.size ( ) / DIM;

        const doffem::CellTransformation<double>& trans = space_.GetCellTransformation ( cell );

        double phys_x;
        double phys_y;
        std::vector<double> pt ( DIM, 0.0 );
        int k = 0;
        for ( int i = 0; i < num_points; i++ )
        {
            for ( int j = 0; j < DIM; j++ )
            {
                pt[j] = ref_coords[k];
                k++;
            }
            phys_x = trans.x ( pt );
            phys_y = trans.y ( pt );
            double fac = std::sqrt ( ( phys_x - 0.25 ) * ( phys_x - 0.25 ) + ( phys_y - 0.25 ) * ( phys_y - 0.25 ) );
            s = 0;
            if ( fac < 0.25 )
            {
                s = std::exp ( -pow ( time_, 10.0 ) ) * std::cos ( 2 * M_PI * fac );
            }
            values[i] = s;
        }
    }
  private:
    const VectorSpace<double>& space_;
    double time_;
};

//eval function of the convection velocity vector for visualization

class EvalConvectionFunction
{
  public:

    EvalConvectionFunction ( const VectorSpace<double>& space, int comp )
    : space_ ( space ), comp_ ( comp )
    {
    }

    void operator() ( const Entity& cell,
            const std::vector<double>& ref_coords,
            std::vector<double>& values ) const
    {
        const int num_points = ref_coords.size ( ) / DIM;

        const doffem::CellTransformation<double>& trans = space_.GetCellTransformation ( cell );

        double phys_x;
        double phys_y;
        std::vector<double> pt ( DIM, 0.0 );
        int k = 0;
        for ( int i = 0; i < num_points; i++ )
        {
            for ( int j = 0; j < DIM; j++ )
            {
                pt[j] = ref_coords[k];
                k++;
            }
            phys_x = trans.x ( pt );
            phys_y = trans.y ( pt );
            double x_q = phys_x - 0.5;
            double y_q = phys_y - 0.5;

            double fac = std::sqrt ( x_q * x_q + y_q * y_q );
            double a_x = fac * ( -y_q );
            double a_y = fac * x_q;
            if ( comp_ == 0 )
            {
                values[i] = a_x;
            }
            else
            {
                values[i] = a_y;
            }
        }
    }
  private:
    const VectorSpace<double>& space_;
    int comp_;
};

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Stiffness-Matrix Assembler.

void StiffnessAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    // Compute local matrix.
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );

    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const int n_dofs = num_dofs ( 0 );
        for ( int i = 0; i < n_dofs; ++i )
        { //loop over test functions
            for ( int j = 0; j < n_dofs; ++j )
            { //loop over ansatz functions
                lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * dot ( grad_phi ( j, q ), grad_phi ( i, q ) ) * std::abs ( detJ ( q ) );
            }
        }
    }
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Mass-Matrix Assembler.

void MassAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    // Compute local matrix.
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );

    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const int n_dofs = num_dofs ( 0 );
        for ( int i = 0; i < n_dofs; ++i )
        { //loop over test functions
            for ( int j = 0; j < n_dofs; ++j )
            { //loop over ansatz functions
                lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * phi ( j, q ) * phi ( i, q ) * std::abs ( detJ ( q ) );
            }
        }
    }
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Convection-Matrix Assembler.

void ConvectionAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    // Compute local matrix.
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );

    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    // Loop over all quadrature points q.
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const int n_dofs = num_dofs ( 0 );
        for ( int i = 0; i < n_dofs; ++i )
        { //loop over test functions
            for ( int j = 0; j < n_dofs; ++j )
            { //loop over ansatz functions
                std::vector<double> vel = get_convective_velocity ( q );
                const double prod = vel[0] * grad_phi ( j, q )[0] + vel[1] * grad_phi ( j, q )[1];
                lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * prod * phi ( i, q ) * std::abs ( detJ ( q ) );
            }
        }
    }
}

// -------------------------------------------------------------------

std::vector<double> ConvectionAssembler::get_convective_velocity ( int q )
{
    // Vector field for rotating pulse.
    std::vector<double> a ( DIM, 0.0 );

    // Implement the vector field for the rotating pulse!
    double x_q = x ( q )[0] - 0.5;
    double y_q = x ( q )[1] - 0.5;
    double fac = std::sqrt ( x_q * x_q + y_q * y_q );
    a[0] = -fac * y_q;
    a[1] = fac * x_q;

    return a;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Source Assembler

void SourceAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );

    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const int n_dofs = num_dofs ( 0 );
        for ( int i = 0; i < n_dofs; ++i )
        {
            lv[dof_index ( i, 0 )] += wq * get_source_term ( q ) * phi ( i, q ) * std::abs ( detJ ( q ) );
        }
    }
}

// -------------------------------------------------------------------

double SourceAssembler::get_source_term ( int q )
{
    double s = 0.0;
    // Implement the source term!
    double x_q = x ( q )[0] - 0.25;
    double y_q = x ( q )[1] - 0.25;

    double fac = std::sqrt ( x_q * x_q + y_q * y_q );
    if ( fac <= 0.25 )
    {
        s = std::cos ( 2 * M_PI * fac );
    }
    return s;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Convection-Diffusion Assembler

class ConvectionDiffusionAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    ConvectionDiffusionAssembler ( double nu )
    : nu_ ( nu ),
    time_ ( 0.0 )
    {
    }

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalMatrix& lm );
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature,
            LocalVector& lv );

    void set_timestep_parameters ( double theta, double delta_t );
    void set_time ( double time );
    void set_prev_solution ( const CVector* sol_prev_ );

    double get_source_term ( int q, double t );
    std::vector<double> get_convective_velocity ( int q );

  private:
    const double nu_;

    // For timestepping

    // These parameters specify which time discretization method you choose
    // See function set_timestep_parameters(double theta, double delta_t)
    double alpha1_;
    double alpha2_;
    double alpha3_;

    double time_;
    const CVector* prev_time_sol_;
    FunctionValues<double> prev_ts_c_;
    FunctionValues< Vec<DIM, double> > grad_prev_ts_c_;
};

// -------------------------------------------------------------------

// Convection-Diffusion Assembler.

void ConvectionDiffusionAssembler::operator() ( const Element<double>& element,
        const Quadrature<double>& quadrature,
        LocalMatrix& lm )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
    // For timestepping.
    prev_ts_c_.clear ( );
    grad_prev_ts_c_.clear ( );
    evaluate_fe_function ( *prev_time_sol_, 0, prev_ts_c_ );
    evaluate_fe_function_gradients ( *prev_time_sol_, 0, grad_prev_ts_c_ );

    // Compute local matrix.
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );

    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    // Loop over all quadrature points q.
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const int n_dofs = num_dofs ( 0 );
        for ( int i = 0; i < n_dofs; ++i )
        { // Loop over test functions
            for ( int j = 0; j < n_dofs; ++j )
            { // Loop over ansatz functions
                std::vector<double> vel = get_convective_velocity ( q );
                const double prod = vel[0] * grad_phi ( j, q )[0] + vel[1] *
                        grad_phi ( j, q )[1];
                lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) += wq * ( phi ( j, q ) * phi ( i, q ) +
                        alpha1_ * ( prod * phi ( i, q ) +
                        nu_ * dot ( grad_phi ( j, q ),
                                    grad_phi ( i, q ) ) ) ) *
                        std::abs ( detJ ( q ) );
            }
        }
    }
}

// -------------------------------------------------------------------

void ConvectionDiffusionAssembler::operator() ( const Element<double>& element,
        const Quadrature<double>& quadrature,
        LocalVector& lv )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
    // For timestepping
    prev_ts_c_.clear ( );
    grad_prev_ts_c_.clear ( );
    evaluate_fe_function ( *prev_time_sol_, 0, prev_ts_c_ );
    evaluate_fe_function_gradients ( *prev_time_sol_, 0, grad_prev_ts_c_ );

    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );

    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    // Compute local vector
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const int n_dofs = num_dofs ( 0 );
        for ( int i = 0; i < n_dofs; ++i )
        { // Loop over test functions
            for ( int j = 0; j < n_dofs; ++j )
            { // Loop over ansatz functions
                std::vector<double> vel = get_convective_velocity ( q );
                const double prod = vel[0] * grad_phi ( j, q )[0] + vel[1] *
                        grad_phi ( j, q )[1];
                const double prev_q = prev_ts_c_[q];
                lv[dof_index ( i, 0 )] += wq * ( prev_ts_c_[q] * phi ( j, q ) * phi ( i, q ) +
                        alpha1_ * get_source_term ( q, time_ ) * phi ( j, q ) *
                        phi ( i, q ) + alpha3_ *
                        ( get_source_term ( q, time_ - alpha2_ ) *
                        phi ( j, q ) * phi ( i, q ) - prev_q * prod *
                        phi ( i, q ) - nu_ * prev_q *
                        dot ( grad_phi ( i, q ), grad_phi ( j, q ) ) ) ) *
                        std::abs ( detJ ( q ) );
            }
        }
    }
}

// -------------------------------------------------------------------

double ConvectionDiffusionAssembler::get_source_term ( int q, double t )
{
    double s = 0.0;
    // Source term
    double x_q = x ( q )[0] - 0.25;
    double y_q = x ( q )[1] - 0.25;

    double fac = std::sqrt ( x_q * x_q + y_q * y_q );
    if ( fac <= 0.25 )
    {
        s = std::exp ( -pow ( t, 10.0 ) ) * std::cos ( 2 * M_PI * fac );
    }
    return s;
}

// -------------------------------------------------------------------

std::vector<double> ConvectionDiffusionAssembler::get_convective_velocity ( int q )
{
    // Vector field for rotating pulse
    std::vector<double> a ( DIM, 0.0 );

    double x_q = x ( q )[0] - 0.5;
    double y_q = x ( q )[1] - 0.5;
    double fac = std::sqrt ( x_q * x_q + y_q * y_q );
    a[0] = -fac * y_q;
    a[1] = fac * x_q;

    return a;
}

// -------------------------------------------------------------------

void ConvectionDiffusionAssembler::set_prev_solution ( const CVector* prev_sol )
{
    prev_time_sol_ = prev_sol;
}

// -------------------------------------------------------------------

// Specify which time discretization method you choose
// theta = 0: ExplicitEuler, theta=0.5: CrankNicolson, theta=1: ImplicitEuler

void ConvectionDiffusionAssembler::set_timestep_parameters ( double theta,
                                                             double delta_t )
{
    alpha1_ = theta * delta_t;
    alpha2_ = delta_t;
    alpha3_ = ( 1. - theta ) * delta_t;
}

// -------------------------------------------------------------------

void ConvectionDiffusionAssembler::set_time ( double time )
{
    time_ = time;
}

// -------------------------------------------------------------------
// -------------------------------------------------------------------

// Convection-Diffusion Application

class ConvectionDiffusionApp
{
  public:
    ConvectionDiffusionApp ( const std::string& param_filename,
                             const std::string& path_mesh );
    ~ConvectionDiffusionApp ( ); //  Destructor.
    void run ( ); //  Run the application.

  private:
    // Member functions.

    // Read the paramaters from provided XML file.
    void read_parameters ( );

    // Read and distribute mesh.
    void read_and_distribute_mesh ( );

    // Setup space, linear algebra, and compute Dirichlet values.
    void prepare_system ( );

    // Compute the matrix and rhs.
    void assemble_system ( );

    // Compute solution.
    void solve_system ( );

    // Visualize the results at time t.
    void visualize ( int time );

    // Member variables.
    PropertyTree params_; // XML tree with parameters.

    MPI_Comm comm_; // MPI communicator.
    int rank_; // Local process rank.
    int num_partitions_; // Number of processes.
    const int master_rank_; // Master rank.

    int refinement_level_; // Refinement level.
    int fe_degree_; // Degree of local shape functions-
    int maxits_; // Maximum number of iterations in linear solver-
    //  Absolute tolerance for residual: if |b - Au| < abs_tol, the iteration is
    //  finished.
    double abstol_;
    //  Relative tolerance for residual: if |b - Au|/|u| < rel_tol, the iteration
    //  is finished.
    double reltol_;
    // Tolerance for stopping the iteration in case of divergence.
    double divtol_;
    double basis_size_; // Size of the Krylov basis for GMRES.

    std::string path_mesh; // Path to geometry data.
    MeshPtr mesh_; // Local mesh.
    VectorSpace<double> space_; // Solution space.
    Couplings<double> couplings_; // Linear algebra couplings helper object.
    CoupledMatrix<Scalar>* matrix_; // System matrix.
    // Vectors for solution, for previous time solution and load vector.
#    ifndef RUNGE_KUTTA
    CoupledVector<Scalar>* sol_, *sol_prev_, *rhs_;
#    else
    CVector sol_, sol_prev_, rhs_;
#    endif

    CVector source_; // Vector for source term.

    CVector k1; // Vector for time-stepping
    CVector k2; // Vector for time-stepping
    CVector k3; // Vector for time-stepping
    CVector k4; // Vector for time-stepping
    CVector u_; // Vector for time-stepping
    CVector help1; // Vector for time-stepping
    CVector help2; // Vector for time-stepping
    CVector help3; // Vector for time-stepping
    CVector help4; // Vector for time-stepping
    CVector s; // Source Vector

    CMatrix mass_; // Mass matrix.
    CMatrix stiff_; // Stiffness matrix.
    CMatrix conv_; // Convection matrix.

    // Dof id:s for Dirichlet boundary conditions.
    std::vector<int> dirichlet_dofs_;
    // Dof values for Dirichlet boundary conditions.
    std::vector<Scalar> dirichlet_values_;

    StandardGlobalAssembler<double> global_asm_; // Global assembler.

#    ifndef RUNGE_KUTTA
    LinearSolver<LAD>* solver_;
    LinearSolverFactory<LAD> SolFact; // Solver for the linear system.
#    else
    CG<LAD> solver_; // Solver for the linear system.
#    endif

    double nu_; // Coefficient of diffusivity.

    std::string method; // Defines which method is used for time discretization.

    double theta_; // Timestepping parameter to choose method.
    double delta_t_; // Size of time step.
    int Tmax_; // Number of time iterations.
    double time_; // Actual simulation time.
};

#endif /* _INSTATIONARY_CONVDIFF_TUTORIAL_H_ */
