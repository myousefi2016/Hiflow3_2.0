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

#ifndef HIFLOW_FSI_BENCHMARK_H_
#    define HIFLOW_FSI_BENCHMARK_H_

/// \author Jonas Kratzke

#    include <mpi.h>
#    include "hiflow.h"

using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

int main ( int argc, char** argv );

typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType CVector;
typedef LAD::MatrixType CMatrix;
typedef std::vector<Scalar> Coord;

/// Geometric dimension
static const int DIM = 2;

/// \brief Instationary assembler for monolithic FSI in the
///        Arbitrary Lagrangian Eulerian frame

class FSIBenchmarkAssembler : private AssemblyAssistant<DIM, Scalar>
{
  public:
    // parameters:
    /// \param rho        density of both the fluid and solid material
    /// \param nu         kinematic viscosity
    /// \param lambda     first Lamee cofficient
    /// \param mu         second Lamee coefficient
    /// \param mesh_diff  fluid mesh regularity coefficient
    /// \param linear     switch for linear/non-linear elasticity

    FSIBenchmarkAssembler (
                            Scalar rho,
                            Scalar nu,
                            Scalar lambda,
                            Scalar mu,
                            Scalar mesh_diff,
                            bool linear )
    : nu_ ( nu ), rho_ ( rho ),
    lambda_ ( lambda ), mu_ ( mu ),
    mesh_diff_ ( mesh_diff ),
    linear_ ( linear )
    {
    }

    void SetMaterials ( int fluid_mat, int solid_mat, int outflow_mat )
    {
        fluid_mat_ = fluid_mat;
        solid_mat_ = solid_mat;
        outflow_mat_ = outflow_mat;
    }

    void SetTimeStepping ( Scalar theta, Scalar dt )
    {
        theta_ = theta;
        dt_ = dt;
    }

    void SetNewtonSolution ( const CVector* u )
    {
        u_ = u;
    }

    void SetTimeSolution ( const CVector* prev_sol )
    {
        prev_sol_ = prev_sol;
    }

    void operator() ( const Element<Scalar>& element,
            const Quadrature<Scalar>& quadrature, LocalMatrix& lm );
    void operator() ( const Element<Scalar>& element,
            const Quadrature<Scalar>& quadrature, LocalVector& lv );
    void operator() ( const Element<Scalar>& element, const int facet_number,
            const Quadrature<Scalar>& quadrature, LocalMatrix& lm );
    void operator() ( const Element<Scalar>& element, const int facet_number,
            const Quadrature<Scalar>& quadrature, LocalVector& lv );

  protected:
    int fluid_mat_, solid_mat_, outflow_mat_;
    Scalar nu_, rho_, lambda_, mu_, mesh_diff_;
    Scalar theta_, dt_;
    bool linear_;

    const CVector *u_, *prev_sol_;
    FunctionValues<Scalar> vel_ns_[DIM], vel_ts_[DIM];
    FunctionValues<Scalar> disp_ns_[DIM], disp_ts_[DIM];
    FunctionValues<Scalar> p_ns_;
    FunctionValues< Vec<DIM, Scalar> > grad_vel_ns_[DIM], grad_vel_ts_[DIM];
    FunctionValues< Vec<DIM, Scalar> > grad_disp_ns_[DIM], grad_disp_ts_[DIM];
};

/// \brief DG-Assembler for the integration of the fluid flow forces
///        acting on the obstacle and elastic bar

class FSIForceIntegrator : private DGAssemblyAssistant<DIM, Scalar>
{
  public:
    /// \param sol             FEM solution vector
    /// \param material_num    interface material number
    /// \param fluid_material  material number of fluid cells
    /// \param dyn_visc        dynamic viscosity (\$rho*\nu$)
    /// \param var             direction of the force to be calculated

    FSIForceIntegrator (
                         const CVector& sol,
                         const int material_num,
                         const int fluid_material,
                         const Scalar dyn_visc,
                         const int var )
    : sol_ ( sol ),
    mat_num_ ( material_num ),
    fluid_material_ ( fluid_material ),
    dyn_visc_ ( dyn_visc ),
    var_ ( var )
    {
    }

    void operator() ( const Element<Scalar>& left_elem, const Element<Scalar>& right_elem,
            const Quadrature<Scalar>& left_quad, const Quadrature<Scalar>& right_quad,
            int left_facet_number, int right_facet_number,
            InterfaceSide left_if_side, InterfaceSide right_if_side,
            Scalar& force );

  private:
    const CVector& sol_;
    FunctionValues< Scalar > pressure_;
    FunctionValues< Vec<DIM, Scalar> > grad_vel_[DIM], grad_disp_[DIM];
    const int mat_num_, fluid_material_, var_;
    const Scalar dyn_visc_;
};

/// \brief DG-Assembler for the integration of the fluid flow magnitude forces
///        acting on fluid-structure interface of the Couette scenario

class FSIForceMagnIntegrator : private DGAssemblyAssistant<DIM, Scalar>
{
  public:
    /// \param sol             FEM solution vector
    /// \param material_num    interface material number
    /// \param fluid_material  material number of fluid cells
    /// \param dyn_visc        dynamic viscosity (\$rho*\nu$)

    FSIForceMagnIntegrator (
                             const CVector& sol,
                             const int material_num,
                             const int fluid_material,
                             const Scalar dyn_visc )
    : sol_ ( sol ),
    mat_num_ ( material_num ),
    fluid_material_ ( fluid_material ),
    dyn_visc_ ( dyn_visc )
    {
    }

    void operator() ( const Element<Scalar>& left_elem, const Element<Scalar>& right_elem,
            const Quadrature<Scalar>& left_quad, const Quadrature<Scalar>& right_quad,
            int left_facet_number, int right_facet_number,
            InterfaceSide left_if_side, InterfaceSide right_if_side,
            Scalar& force );

  private:
    const CVector& sol_;
    FunctionValues< Vec<DIM, Scalar> > grad_vel_[DIM], grad_disp_[DIM];
    const int mat_num_, fluid_material_;
    const Scalar dyn_visc_;
};

/// \brief Assembler for the integration of the L2 error of the Couette fluid velocity field

class flowL2error : private AssemblyAssistant<DIM, Scalar>
{
  public:
    /// \param sol             FEM solution vector
    /// \param fluid_material  material number of fluid cells
    /// \param vel0            radial velocity at the inner boundary
    /// \param radii           radii of the geometry

    flowL2error (
                  const CVector& sol,
                  const int fluid_material,
                  const Scalar vel0,
                  const std::vector<Scalar>& radii )
    : sol_ ( sol ),
    fluid_material_ ( fluid_material ),
    vel0_ ( vel0 ),
    radii_ ( radii )
    {
    }

    Scalar exact_vel ( int var, const Vec<DIM, Scalar>& x );

    void operator() ( const Element<Scalar>& element,
            const Quadrature<Scalar>& quadrature, Scalar& lm );

  private:
    const CVector& sol_;
    FunctionValues< Scalar > vel_[DIM], disp_[DIM];
    const int fluid_material_;
    const Scalar vel0_;
    const std::vector<Scalar>& radii_;
};

/// \brief Assembler for the integration of the L2 error of the Couette pressure field

class pressureL2error : private AssemblyAssistant<DIM, Scalar>
{
  public:
    /// \param sol             FEM solution vector
    /// \param fluid_material  material number of fluid cells
    /// \param vel0            radial velocity at the inner boundary
    /// \param radii           radii of the geometry

    pressureL2error (
                      const CVector& sol,
                      const int fluid_material,
                      const Scalar vel0,
                      const Scalar rho,
                      const std::vector<Scalar>& radii )
    : sol_ ( sol ),
    fluid_material_ ( fluid_material ),
    vel0_ ( vel0 ),
    rho_ ( rho ),
    radii_ ( radii )
    {
    }

    Scalar exact_press ( const Vec<DIM, Scalar>& x );

    void operator() ( const Element<Scalar>& element,
            const Quadrature<Scalar>& quadrature, Scalar& lm );

  private:
    const CVector& sol_;
    FunctionValues< Scalar > press_, disp_[DIM];
    const int fluid_material_;
    const Scalar vel0_, rho_;
    const std::vector<Scalar>& radii_;
};

/// \brief Assembler for the integration of the L2 error of the Couette displacement field

class dispL2error : private AssemblyAssistant<DIM, Scalar>
{
  public:
    /// \param sol             FEM solution vector
    /// \param solid_material  material number of solid cells
    /// \param disp1           radial displacement at the interface
    /// \param radii           radii of the geometry

    dispL2error (
                  const CVector& sol,
                  const int solid_material,
                  const Scalar disp1,
                  const std::vector<Scalar>& radii )
    : sol_ ( sol ),
    solid_material_ ( solid_material ),
    disp1_ ( disp1 ),
    radii_ ( radii )
    {
    }

    Scalar exact_disp ( const int var, const Vec<DIM, Scalar>& x );

    void operator() ( const Element<Scalar>& element,
            const Quadrature<Scalar>& quadrature, Scalar& lm );

  private:
    const CVector& sol_;
    FunctionValues< Scalar > disp_[DIM];
    const int solid_material_;
    const Scalar disp1_;
    const std::vector<Scalar>& radii_;
};

/// \brief Assembler for the integration of the L2 error of the Couette displacement diffusion field

class diffusionL2error : private AssemblyAssistant<DIM, Scalar>
{
  public:
    /// \param sol             FEM solution vector
    /// \param solid_material  material number of solid cells
    /// \param disp1           radial displacement at the interface
    /// \param radii           radii of the geometry

    diffusionL2error ( const CVector& sol,
                       const int fluid_material,
                       const Scalar disp1,
                       const std::vector<Scalar>& radii )
    : sol_ ( sol ),
    fluid_material_ ( fluid_material ),
    disp1_ ( disp1 ),
    radii_ ( radii )
    {
    }

    Scalar exact_diff ( const int var, const Vec<DIM, Scalar>& x );

    void operator() ( const Element<Scalar>& element,
            const Quadrature<Scalar>& quadrature, Scalar& lm );

  private:
    const CVector& sol_;
    FunctionValues< Scalar > disp_[DIM];
    const int fluid_material_;
    const Scalar disp1_;
    const std::vector<Scalar>& radii_;
};

/// Main class

class FSIBenchmark : public NonlinearProblem<LAD>
{
  public:
    /// Constructor using property tree
    FSIBenchmark ( PropertyTree& config );

    /// Destructor
    ~FSIBenchmark ( );

    /// Run application
    void run ( );

    /// Newton method functions for the RHS and the Jacobian
    virtual void EvalFunc ( const CVector& u, CVector* F );
    virtual void EvalGrad ( const CVector& u, CMatrix* DF );

  protected:

    enum BenchmarkType
    {
        NOT_SET = 0,
        Channel,
        Couette
    };

    /// Initialization
    void initialize_parameters ( );
    void initialize_mesh ( );
    void initialize_space ( );
    void initialize_la_structures ( );
    void initialize_linear_solver ( );
    void initialize_nonlinear_solver ( );
    void initialize_postprocessing ( );

    /// Assembly routines
    void assemble_residual ( const CVector& u, CVector* F );
    void assemble_matrix ( const CVector& u, CMatrix* DF );

    /// Boundary conditions
    void compute_bc ( );

    /// Handling of unused pressure DoFs in structure domain
    void fix_solid_pressure_dofs ( std::vector<int> &dirichlet_dofs,
                                   std::vector<Scalar> &dirichlet_values );

    /// Solve routine
    void solve_nonlinear_system ( );

    /// Visualization of computation of results
    void postprocessing ( );
    void visualize_solution ( CVector& u, std::string const& filename ) const;

    /// MPI stuff
    const MPI_Comm comm_;
    int rank_, num_partitions_;

    /// Parameter input tree
    PropertyTree& config_;

    /// Prefix for output files
    std::string simul_name_;

    /// Application parameters
    BenchmarkType type_;
    Scalar Vm_, Um_, nu_, rho_, H_, lambda_, mu_, mesh_diff_;
    bool linear_;
    std::vector<Scalar> radii_;
    int inflow_mat_, outflow_mat_;
    int fluid_mat_, solid_mat_;
    int obstacle_mat_, interface_mat_, side_mat_;

    /// Timestepping parameters
    int ts_, visu_intervall_;
    Scalar theta_, t_, Dt_, dt_, start_up_time_;

    /// Mesh parameters
    int refinement_level_;
    MeshPtr mesh_;

    /// FEM parameters
    VectorSpace<Scalar> space_;
    Couplings<Scalar> couplings_;

    std::vector<int> dirichlet_dofs_;
    std::vector<Scalar> dirichlet_values_;

    /// Linear Algebra objects
    CVector *sol_, *rhs_, *res_, *prev_sol_;
    CMatrix* matrix_;

    /// Assembler object
    StandardGlobalAssembler<Scalar> global_asm_;

    /// Solver
    NonlinearSolver<LAD>* nls_;
#    ifdef WITH_MUMPS
    MumpsSolver<LAD, MumpsStructureD>* linear_solver_;
#    else
    LinearSolver<LAD>* linear_solver_;
#    endif
    /// Preconditioning
#    ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> ilupp_;
    bool use_ilupp_;
#    endif

    /// Forcing parameters
    std::string forcing_strategy_;
    Scalar eta_initial_, eta_max_, gamma_EW2_, alpha_EW2_;
    EWForcing<LAD>* ew_forcing_;
    bool newton_forcing_;

    /// Damping parameters
    std::string damping_strategy_;
    Scalar theta_initial_, theta_min_, armijo_dec_, suff_dec_;
    int max_armijo_ite_;
    ArmijoDamping<LAD>* armijo_damping_;
    bool newton_damping_;

    /// CSV Writer for evaluated results
    CSVWriter<Scalar> results_writer_;

    /// location of evaluation point
    Coord eval_point_;

    /// L2 error vectors
    std::vector<Scalar> L2_err_vel_, L2_err_disp_, L2_err_press_, L2_err_diff_;
};

#endif
