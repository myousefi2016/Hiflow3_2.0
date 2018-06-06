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

#ifndef MET_FLOW_BOUSSINESQ_3D_H
#    define MET_FLOW_BOUSSINESQ_3D_H

///
/// \file boussinesq.h
/// \brief Fully implicit version of Boussinesq approximation version of the MetFlow solver.
///
/// \author Martin Baumann, Philipp Gerstner, Jonathan Schwegler
///

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>

#    include "hiflow.h"
#    include "metflow.h"

///
/// Application base class
///

class MetFlowBousCyl3dApp : public MetFlowApp
{
  public:

    MetFlowBousCyl3dApp ( const std::string& root_path, const std::string& param_filename, const std::string base_param_filename, bool resume_run );

    virtual ~MetFlowBousCyl3dApp ( );

    void read ( );

  protected:

    // Linear algebra structures
    std::vector<LAD::DataType> div2_;

    // Linear solver
#    ifdef SCHUR_SOLVER
    FGMRES<LAD>* PVT_linear_solver_;
#    else
    GMRES<LAD>* PVT_linear_solver_;
#        ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> PVT_precond_;
#        endif
#    endif

    // Local assembler
    MetFlowBousCylAssembler<DIM, DATATYPE> PVT_local_asm_;
    ScalarProdCylAssembler<DIM, DIM + 2, DATATYPE> PVT_local_asm_prod_;

    /// ***************************************************************************
    /// INIT
    /// ***************************************************************************
    virtual void initial_prepare ( );
    virtual void dwr_loop_prepare ( );

    virtual void prepare_goal_functional ( )
    {
        ;
    }
    virtual void prepare_space ( VectorSpace<double>& space, std::vector< std::vector< bool> >& coupling_vars, MeshPtr mesh, int mode );

    virtual void prepare_assembler ( );
    virtual void prepare_parameters ( );
    virtual void prepare_ic ( );
    virtual void prepare_bc ( );
    virtual void prepare_bc ( double rotation, double time );
    virtual void perturb_ic ( );

    /// ***************************************************************************
    /// IO
    /// ***************************************************************************
    virtual void write_ic ( );
    virtual void visualize_solution ( int step );
    virtual void visualize_solution ( LAD::VectorType& sol, const VectorSpace<double>* space, std::string const& prefix, int time_step );
    virtual void create_log_files ( int offset );

    /// ***************************************************************************
    /// POSTPROCESSING
    /// ***************************************************************************
    virtual void compute_flow_quantities ( std::string filename, LAD::VectorType* sol, double time );
    virtual void compute_solution_norm ( std::string filename, LAD::VectorType* sol, LAD::VectorType* sol_prev, double time );
    virtual void compute_quant_vector ( std::string filename, double time );

    virtual void compute_char_quant ( )
    {
    };

    /// ***************************************************************************
    /// SOLVER
    /// ***************************************************************************
    virtual void update_preconditioner ( const LAD::VectorType& u, LAD::MatrixType* DF );
    virtual void prepare_linear_solver ( bool ic );

    /// ***************************************************************************
    /// IN-TIME LOOP FUNCTIONS
    /// ***************************************************************************
    virtual void post_processing ( );
    virtual void update_assembler ( );
    virtual void filter_solution ( );

    // ************************************************************************
    // Wavetank-specific variables and functions

    double rad_i_;
    double rad_o_;
    double height_;
    // Physical parameters
    double nu_; // kinematic viscosity                  [m^2/s]
    double rho_; // density                              [kg/m^3]
    double inv_rho_; // inverse density                      [m^3/kg]

    double kappa_; // kappa/(rho*cv)                       [m^2/s]
    double alpha_g_; // thermal expansion coefficient        [1/K]

    // Temperature
    double ref_T_; // reference temperature                [°C]
    double start_T_; // initial temperature                  [°C]
    double warm_T_; // temperature outer boundary (warm)    [°C]
    double cold_T_; // temperature inner temperature (cold) [°C]
    int warm_material_; // material number
    int cold_material_; // material number
    std::vector<double> gravity_;

    // Rotation
    double omega_; // rotation speed                       [U/s]

    // Structure
    const int T_var_;

    // Stabilization
    int graddiv_mode_; // VMS or SUPG -> different choice of parameter tau
    double graddiv_gamma_; // parameter for grad-div stabilization: gamma * tau * (div , div)
    double graddiv_C_; // parameter occuring in tau of VMS (usually approx 9)

    int skew_mode_;

    int temp_supg_mode_;
    double temp_supg_gamma_;

    double Nu_; // Nusselt number
    double Nu0_; // initial nusselt number
    double nusselt_r_min_;
    double nusselt_r_max_;
    int nusselt_surface_id_;
};

#endif
