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

#ifndef MET_FLOW_CONVDIFF_H
#    define MET_FLOW_CONVDIFF_H

///
/// \file conv_diff.h
/// \brief Fully implicit solver for thermo electro dynamic Boussinesq approximation in cylindrical coordinates
///
/// \author Philipp Gerstner
///

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>

#    include "hiflow.h"
#    include "metflow.h"

const static int XDIM = 1;
const static int COSYSTEM = 0;

///
/// Derived application class for coupled Boussinesq solver in cartesian coordinates
///

class MetFlowConvDiffApp : public MetFlowApp
{
  public:

    /// Constructor
    MetFlowConvDiffApp ( const std::string& root_path, const std::string& param_filename, const std::string& base_param_filename, bool resumed_run );

    /// Destructor
    virtual ~MetFlowConvDiffApp ( );

    /// ***************************************************************************
    /// MAIN LOOPS
    /// ***************************************************************************
    void pp_run ( );

    /// ***************************************************************************
    /// INIT
    /// ***************************************************************************
    void initial_prepare ( );
    void dwr_loop_prepare ( );
    void prepare_goal_functional ( );
    void prepare_space ( VectorSpace<double>& space, std::vector< std::vector< bool> >& coupling_vars, MeshPtr mesh, int mode );
    void prepare_space_convection ( );
    void prepare_lin_alg_structures_convection ( );

    // Primal problem
    void prepare_assembler ( );
    void prepare_parameters ( );
    void prepare_ic ( );
    void prepare_bc ( );
    void prepare_bc ( double rotation, double time );
    void perturb_ic ( );

    // Dual Problem
    void prepare_assembler_dual ( );
    void prepare_ic_dual ( );
    void prepare_bc_dual ( );

    /// Error estimation
    void prepare_assembler_est ( );

    void prepare_parameters_est ( )
    {
        /*TODO*/;
    }

    void clear_problem ( );

    /// ***************************************************************************
    /// IO
    /// ***************************************************************************
    void write_ic ( );
    void visualize_solution ( int step );
    void visualize_solution ( LAD::VectorType& sol, const VectorSpace<double>* space, std::string const& prefix, int time_step );
    void visualize_solution_dual ( int step );
    void create_log_files ( int offset );
    void create_log_files_dual ( int offset );
    void create_ic_hdf5_filename ( );

    /// ***************************************************************************
    /// POSTPROCESSING
    /// ***************************************************************************
    void compute_flow_quantities ( std::string filename, LAD::VectorType* sol, double time );
    void compute_solution_norm ( std::string filename, LAD::VectorType* sol, LAD::VectorType* sol_prev, double time );
    void compute_quant_scalar ( );
    void compute_quant_vector ( std::string filename, double time );
    void compute_char_quant ( );

    /// ***************************************************************************
    /// SOLVER
    /// ***************************************************************************
    virtual void update_preconditioner ( const LAD::VectorType& u, LAD::MatrixType* DF );
    virtual void update_preconditioner_dual ( LAD::MatrixType* DF );
    virtual void prepare_linear_solver ( bool ic );
    virtual void prepare_linear_solver_dual ( );

    /// ***************************************************************************
    /// IN-TIME LOOP FUNCTIONS
    /// ***************************************************************************
    virtual void post_processing ( );
    virtual void post_processing_dual ( );
    virtual void update_assembler ( );
    virtual void update_assembler_dual ( );
    virtual void filter_solution ( );
    virtual void update_convection ( double time, bool backward );

  protected:
    // Linear algebra structures
    bool computed_T_mass_;
    bool computed_T_mass_dual_;

    // Other stuff
    LAD::MatrixType T_mass_matrix_;
    LAD::MatrixType T_mass_matrix_dual_;

    // convection field
    LAD::VectorType V_;
    LAD::VectorType V_prev_;
    LAD::VectorType V_next_;
    std::vector<std::vector<bool> > V_coupling_vars_;
    Couplings<double> V_couplings_;

    // Linear solver
    // primal
    FGMRES<LAD>* T_linear_solver_;

    // dual
    FGMRES<LAD>* T_linear_solver_dual_;

    // general

#    ifdef USE_HYPRE
    CG<LAD> T_mass_solver_;
    HypreBoomerAMG<LAD> T_mass_precond_;
    HypreBoomerAMG<LAD> T_precond_;
    HypreBoomerAMG<LAD> T_precond_dual_;
#    else
    GMRES<LAD> T_mass_solver_;
#        ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> T_mass_precond_;
    PreconditionerIlupp<LAD> T_precond_;
    PreconditionerIlupp<LAD> T_precond_dual_;
#        endif
#    endif

    VectorSpace<double> V_space_;

    // Local assembler
    MetFlowConvDiffCartAssembler <DIM, DATATYPE> T_local_asm_;
    MetFlowConvDiffCartDualAssembler <DIM, DATATYPE> T_local_asm_dual_;
    MetFlowConvDiffCartEstimatorAssembler<DIM, DATATYPE> T_local_asm_est_;
    MassMatrixCylAssembler <DIM, VARDIM, DATATYPE> T_local_asm_mass_;
    ScalarProdCylAssembler <DIM, VARDIM, DATATYPE> T_local_asm_prod_;

    // Goal functional
    GoalFunctionalVariableSubDomain<DIM, DATATYPE> j_var_on_sub_;

    // variable indices

    int t_var_;

    // Fluid parameters
    double kappa_; // kappa/(rho*cv)                       [m^2/s]
    double lambda_;
    double gamma_;

    // Experimental Setup
    double start_T_; // initial temperature                  [°C]
    std::vector<int> dirichlet_bdy_;
    std::vector<double> dirichlet_val_;
};

#endif
