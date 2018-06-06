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

#ifndef _CHORIN_H_
#    define _CHORIN_H_

/// \brief This application solves the instationary incompressible Navier Stokes
/// equation with Chorin:s method, allows for pressure boundary conditions.
/// \author Thomas Gengenbach

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>
#    include <sstream>
#    include <fstream>
#    include <iostream>
#    include <iomanip>

#    include "hiflow.h"

using namespace hiflow;
using namespace hiflow::la;
using namespace hiflow::mesh;
using namespace hiflow::doffem;

#    define NEW_ASSEMBLY

static const int DIM = 2;
//static const int DIM = 3;
typedef LADescriptorCoupledD LAD;

typedef std::vector<double> Coord;

int main ( int argc, char** argv );

class IntermediateVelocityAssembler;
class PressureCorrectionAssembler;
class VelocityCorrectionAssembler;

class InstationaryFlowApp : public NonlinearProblem<LAD>
{
  public:

    static const PLATFORM APP_PLATFORM = CPU;
    static const IMPLEMENTATION APP_LINALG_IMPLEMENTATION = NAIVE; //MKL;
    static const MATRIX_FORMAT APP_MATRIX_FORMAT = CSR;

    InstationaryFlowApp ( );
    ~InstationaryFlowApp ( );

    void run ( );

    // Callback functions for nonlinear solver
    virtual void EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F );
    virtual void EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF );

    void compute_residual ( const LAD::VectorType& u, LAD::VectorType* F );
    void compute_matrix ( const LAD::VectorType& u, LAD::MatrixType* F );

  private:

    void read_and_distribute_mesh ( const std::string& filename );

    void prepare_space ( );
    void prepare_linear_solver ( );
    void prepare_lin_alg_structures ( );
    void prepare_bc ( double time = 0.0 );
    void prepare_nls ( );

    // Chorin:s method
    void solve_nlp ( );
    void solve_pressure_correction ( );
    void solve_velocity_correction ( );

    void visualize ( LAD::VectorType& sol, VectorSpace<double>& space, int time_step, std::string name, std::vector<int> vars );

    void compute_alphas ( int sub_step, double delta_t, std::vector<double>* alphas );

    MPI_Comm comm_; // must use concrete object here -> hence C interface
    int rank_;
    int num_partitions_;
    const int master_rank_;

    MeshPtr mesh_;
    VectorSpace<double> intermediate_space_, pressure_space_, velocity_space_;
    int refinement_level_;

    // Three assemblers needes for Chorin:s method
    IntermediateVelocityAssembler* intermediate_asm_;
    PressureCorrectionAssembler* pressure_asm_;
    VelocityCorrectionAssembler* velocity_asm_;

    /// Time-stepping method
    std::string method_;

    /// Variable names for visualization
    std::vector<std::string> visualization_names_;

    SYSTEM la_sys_;
    Couplings<double> intermediate_couplings_, pressure_couplings_, velocity_couplings_;

    /// Density rho
    double rho_;

    /// Viscosity nu
    double nu_;

    /// Height of the flow channel
    double H_;

    LAD::MatrixType intermediate_matrix_, pressure_matrix_, velocity_matrix_;
    LAD::VectorType intermediate_sol_, pressure_sol_, velocity_sol_, sol_prev_,
    intermediate_rhs_, pressure_rhs_, velocity_rhs_, res_;
    std::vector<int> intermediate_dirichlet_dofs_, pressure_dirichlet_dofs_, velocity_dirichlet_dofs_;
    std::vector<LAD::DataType> intermediate_dirichlet_values_, pressure_dirichlet_values_, velocity_dirichlet_values_;
    // Global assembler.
    StandardGlobalAssembler<double> global_asm_;
    // would be better to use pointer to base class here, but not
    // possible since these objects do not have sensible constructors...
    ScopedPtr< Newton<LAD> >::Type nls_;
#    ifdef WITH_ILUPP
    // ilupp preconditioner
    ScopedPtr< PreconditionerIlupp<LAD> >::Type ilupp_;
#    endif
    ScopedPtr< GMRES<LAD> >::Type intermediate_linear_solver_;
    ScopedPtr< CG<LAD> >::Type pressure_linear_solver_;
    ScopedPtr< GMRES<LAD> >::Type velocity_linear_solver_;

};

// IntermediateVelocityAssembler ////////////////

class IntermediateVelocityAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    IntermediateVelocityAssembler ( double nu, double rho )
    : nu_ ( nu ), inv_rho_ ( 1. / rho )
    {
    }

    void set_newton_solution ( const LAD::VectorType& newton_sol );
    void set_prev_solution ( const LAD::VectorType& prev_sol );
    void set_timestep_parameters ( const std::vector<double>& alphas );

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature );

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv );

  private:
    LAD::VectorType const* newton_sol_;
    LAD::VectorType const* prev_sol_;

    double nu_, inv_rho_;
    double alpha1_, alpha2_, alpha3_, alpha4_, alpha5_;
    FunctionValues<double> prev_ns_vel_[DIM];
    FunctionValues<double> prev_ts_vel_[DIM];
    FunctionValues< Vec<DIM, double> > grad_prev_ns_vel_[DIM];
    FunctionValues< Vec<DIM, double> > grad_prev_ts_vel_[DIM];
};

// PressureCorrectionAssembler ////////////////

class PressureCorrectionAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    PressureCorrectionAssembler ( double rho )
    : rho_ ( rho )
    {
    }

    void set_intermediate_velocity_solution ( const LAD::VectorType& intermediate_sol );
    void set_timestep ( double delta_t );

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature );

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv );

  private:
    LAD::VectorType const* intermediate_sol_;

    double rho_;
    double delta_t_;
    FunctionValues< Vec<DIM, double> > intermediate_velocity_gradient_[DIM];
};

// VelocityCorrectionAssembler ////////////////

class VelocityCorrectionAssembler : private AssemblyAssistant<DIM, double>
{
  public:

    VelocityCorrectionAssembler ( double rho )
    : rho_ ( rho )
    {
    }

    void set_pressure_solution ( const LAD::VectorType& pressure_sol );
    void set_intermediate_velocity_solution ( const LAD::VectorType& intermediate_sol );
    void set_timestep ( double delta_t );

    void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature );

    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm );
    void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv );

  private:
    double rho_;
    double delta_t_;

    LAD::VectorType const* pressure_sol_;
    FunctionValues< Vec<DIM, double> > pressure_gradient_;

    LAD::VectorType const* intermediate_sol_;
    FunctionValues<double> intermediate_velocity_[DIM];
};

#endif /* _CHORIN_H_ */
