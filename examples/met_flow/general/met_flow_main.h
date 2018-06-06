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

#ifndef MET_FLOW_MAIN_H
#    define MET_FLOW_MAIN_H

///
/// \file met_flow.h
/// \brief Assembler and application base class for meteorologic application.
///
/// \author Martin Baumann, Teresa Beck, Volker Lange, Simon Gawlok, Philipp Gerstner
///

// BUGS:
// fehler in letzten (t=0) zeitschritt in dualer lösung für galerkin CD
// mesh balance

// RESTRUKTURIERUNG:
// DynamicMeshHandler
// FeInterpolator
// Ordner Adaptivity

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>
#    include <sstream>
#    include <algorithm>

#    include "hiflow.h"
#    include "goal_functional.h"
#    include "../assembler/scalar_prod_assembler.h"
#    include "../assembler/met_flow_assembler.h"
#    include "../assembler/met_flow_est_assembler.h"
#    include "my_newton.h"
#    include "write_csv.h"
#    include "../tmp_config/met_flow_vars.h"

/// The default quadrature selection chooses a quadrature rule that is accurate to 2 * max(fe_degree).

struct QuadratureSelection
{

    QuadratureSelection ( int order ) : order_ ( order )
    {
    }

    void operator() ( const Element<double>& elem, Quadrature<double>& quadrature )
    {
        const FEType<double>::FiniteElement fe_id = elem.get_fe_type ( 0 )->get_my_id ( );

        switch ( fe_id )
        {
            case FEType<double>::LAGRANGE_LINE:
                quadrature.set_cell_type ( 1 );
                quadrature.set_quadrature_by_order ( "GaussLine", order_ );
                break;
            case FEType<double>::LAGRANGE_TRI:
                quadrature.set_cell_type ( 2 );
                quadrature.set_quadrature_by_order ( "GaussTriangle", order_ );
                break;
            case FEType<double>::LAGRANGE_QUAD:
                quadrature.set_cell_type ( 3 );
                quadrature.set_quadrature_by_order ( "GaussQuadrilateral", order_ );
                break;
            case FEType<double>::LAGRANGE_TET:
                quadrature.set_cell_type ( 4 );
                quadrature.set_quadrature_by_order ( "GaussTetrahedron", order_ );
                break;
            case FEType<double>::LAGRANGE_HEX:
                quadrature.set_cell_type ( 5 );
                quadrature.set_quadrature_by_order ( "GaussHexahedron", order_ );
                break;
            default:
                assert ( false );
        };
    }

    int order_;
};

///
/// Application base class
///

class MetFlowApp : public NonlinearProblem<LAD>, public DynamicMeshProblem<LAD, MESHIMPL, DIM>
{
  public:

    static const PLATFORM APP_PLATFORM = CPU;
    static const IMPLEMENTATION APP_LINALG_IMPLEMENTATION = NAIVE;
    static const MATRIX_FORMAT APP_MATRIX_FORMAT = CSR;

    MetFlowApp ( const std::string& root_path, const std::string& param_filename, const std::string base_param_filename, bool resumed_run );
    virtual ~MetFlowApp ( );

    /// *****************************************
    /// MAIN LOOPS
    /// *****************************************
    virtual void run ( );
    virtual void dwr_run ( );

    virtual void pp_run ( )
    {
    };
    virtual void qoi_run ( );
    virtual void perturb_test_run ( );

    /// *****************************************
    /// Virtual functions of class NonLinearProblem
    /// *****************************************
    virtual void EvalFunc ( const LAD::VectorType& u, LAD::VectorType* F );
    virtual void EvalGrad ( const LAD::VectorType& u, LAD::MatrixType* DF );

    void compute_residual ( const LAD::VectorType* u, LAD::VectorType* F, int mode );
    void compute_matrix ( const LAD::VectorType* u, LAD::MatrixType* DF, int mode );

    /// *****************************************
    /// Virtual functions of DynamicMeshProblem
    /// *****************************************
    virtual void setup_LA_primal ( VectorSpace<DATATYPE>& space, VectorSpace<DATATYPE>& space_dual,
                                   std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                   Couplings<DATATYPE>& couplings, Couplings<DATATYPE>& couplings_dual );

    virtual void setup_LA_dual ( VectorSpace<DATATYPE>& space, VectorSpace<DATATYPE>& space_dual,
                                 std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                 Couplings<DATATYPE>& couplings, Couplings<DATATYPE>& couplings_dual );

    virtual void setup_LA_est ( VectorSpace<DATATYPE>& space, VectorSpace<DATATYPE>& space_dual,
                                std::vector< std::vector< bool> >& coupling_vars, std::vector< std::vector< bool> >& coupling_vars_dual,
                                Couplings<DATATYPE>& couplings, Couplings<DATATYPE>& couplings_dual );

    virtual void setup_space ( VectorSpace<DATATYPE>& space, std::vector< std::vector< bool> >& coupling_vars, MeshPtr mesh, int mode );

    virtual void set_active_space ( VectorSpace<typename LAD::DataType>* space, int mode );

    virtual void set_active_mesh ( MeshPtr mesh );

    virtual void init_mesh_change_list ( );

    virtual void set_update_vectors ( );

  protected:
    /// *****************************************
    /// INITIALIZATION
    /// *****************************************
    virtual void initial_prepare ( ) = 0;
    virtual void dwr_loop_prepare ( );
    virtual void periodify_space ( VectorSpace<DATATYPE>& space, MeshPtr mesh );
    virtual void prepare_space ( VectorSpace<DATATYPE>& space, std::vector< std::vector< bool> >& coupling_vars, MeshPtr mesh, int mode ) = 0;
    virtual void build_initial_mesh ( int adapt_counter );
    virtual void build_initial_mesh_parallel ( int adapt_counter );
    virtual void build_initial_mesh_p4est ( int adapt_counter );
    virtual void build_initial_time_mesh ( int adapt_counter );
    virtual void prepare_periodicity ( );
    virtual void set_bc ( int mode );

    /// Primal problem
    virtual void prepare_assembler ( );
    virtual void prepare_parameters ( );
    virtual void prepare_time_method ( );
    virtual void prepare_bc ( ) = 0;
    virtual void prepare_ic ( ) = 0;
    virtual void prepare_nls ( bool stationary );
    virtual void prepare_lin_alg_structures ( Couplings<DATATYPE>& couplings );
    virtual int init_run ( bool run_mode0 );

    /// Dual Problem

    virtual void prepare_goal_functional ( )
    {
    };
    virtual void prepare_assembler_dual ( );
    virtual void prepare_time_method_dual ( );

    virtual void prepare_ic_dual ( )
    {
    };

    virtual void prepare_bc_dual ( )
    {
    };
    virtual void prepare_lin_alg_structures_dual ( Couplings<DATATYPE>& couplings );
    virtual int init_run_dual ( );

    /// Error estimation
    virtual void prepare_assembler_est ( );

    /// clear data structures
    virtual void clear_problem ( );

    /// *****************************************
    /// SOLVER
    /// *****************************************

    /// time loop
    virtual int solve_primal ( int adapt, int time_step );
    virtual int solve_dual ( int adapt, int time_step );

    /// algebraic solver
    virtual bool solve_nlp ( int mode );
    virtual bool solve_lp ( int mode );

    /// linear solver
    virtual void prepare_linear_solver ( bool ic ) = 0;

    virtual void prepare_linear_solver_dual ( )
    {
    };

    /// preconditioner
    virtual void update_preconditioner ( const LAD::VectorType& u, LAD::MatrixType* DF ) = 0;

    virtual void update_preconditioner_dual ( LAD::MatrixType* DF )
    {
    };

    /// ***************************************************************************
    /// IN-TIME LOOP FUNCTIONS
    /// ***************************************************************************
    /// functions that are called for each time step

    virtual void post_processing ( )
    {
    };

    virtual void post_processing_dual ( )
    {
    };

    virtual void update_assembler ( )
    {
    };

    virtual void update_assembler_dual ( )
    {
    };

    virtual void filter_solution ( )
    {
    };
    virtual void perturb_solution ( int time_step );

    /// *****************************************
    /// IO
    /// *****************************************
    /// Visualization

    virtual void visualize_solution ( int time_step )
    {
        ;
    }

    virtual void visualize_solution_dual ( int time_step )
    {
        ;
    };
    virtual void visualize_patch_interpolation ( int time_step );
    virtual void visualize_function ( const LAD::VectorType& sol, const VectorSpace<double>* space, std::string const& prefix, int time_step, std::vector<std::string>& var_names );

    /// HDF5 IO
    virtual void read_file ( LAD::VectorType& sol, std::string const& filename, std::string const& groupname = "NoHDF5", std::string const& datasetname = "NoHDF5" );
    virtual void read_file ( LAD::VectorType& sol, int time_step, std::string const& filename, std::string const& groupname = "NoHDF5", std::string const& datasetname = "NoHDF5", int mode = 1 );
    virtual void write_file ( LAD::VectorType& sol, std::string const& filename, std::string const& groupname = "NoHDF5", std::string const& datasetname = "NoHDF5" );
    virtual void write_ic ( ) = 0;

    virtual void create_ic_hdf5_filename ( )
    {
    };

    /// Log data

    virtual void create_log_files ( int offset )
    {
    };

    virtual void create_log_files_dual ( int offset )
    {
    };
    virtual void compute_comm_pattern ( std::string const& filename );
    void create_log_file ( int offset, std::string filename, std::string suffix = "txt", std::string seperator = " ", bool print_counter = true );

    /// *****************************************
    /// POSTPROCESSING
    /// *****************************************

    virtual void prepare_postprocessing ( int adaptation_counter = -1 )
    {
        ;
    }

    virtual void postprocessing ( int time_step, int adaptation_counter = -1 )
    {
        ;
    }
    virtual void evaluate_inner_prod ( std::string filename, LAD::VectorType* solA, LAD::VectorType* solB, double time );
    virtual void evaluate_qoi_fin ( std::string filename, LAD::VectorType* sol, double time );
    virtual void evaluate_qoi_int ( std::string filename, LAD::VectorType* sol, double time );

    /// *****************************************
    /// Adaptivity
    /// *****************************************
    /// Time mesh
    virtual double get_delta_t ( int time_step );
    virtual int get_num_intervals ( );
    virtual double get_time ( int time_step );

    /// computation of error indicators
    virtual void prepare_higher_order_space ( );
    virtual void compute_cell_diameters ( std::vector<double>& hK );

    virtual void estimate_error ( int adaption_counter, int num_time_steps );

    virtual void assemble_intra_time_indicators ( int rel_time,
                                                  std::vector< std::vector<double> >& est_cell_residual,
                                                  std::vector< std::vector<double> >& est_interface );

    virtual void accumulate_residuals_and_weights ( int est_type,
                                                    const std::vector< std::vector< std::vector<double> > >& est_cell_residual,
                                                    const std::vector< std::vector< std::vector<double> > >& est_cell_timejump,
                                                    const std::vector< std::vector< std::vector<double> > >& est_interface,
                                                    std::vector< std::vector<double> >& element_residual,
                                                    std::vector< std::vector<double> >& element_weight );

    virtual void compute_estimators ( const std::vector< std::vector<double> >& element_residual_h,
                                      const std::vector< std::vector<double> >& element_residual_tau,
                                      const std::vector< std::vector<double> >& element_weight_h,
                                      const std::vector< std::vector<double> >& element_weight_tau );

    virtual void compute_reduced_estimators ( );

    /// Mesh adaption
    virtual void adapt_spatial_mesh ( );
    virtual void adapt_temporal_mesh ( );
    virtual void adapt_mesh_change_list ( );

    /// Mesh change
    virtual void read_mesh_change_list ( int counter, std::vector<DATATYPE>& mesh_change_times );
    virtual void write_mesh_change_list ( int counter, std::vector<DATATYPE>& mesh_change_times );

    /// Visualization of error indicators
    virtual void visualize_error_indicators ( int timestep );
    virtual void visualize_reduced_error_indicators ( int timestep, int mesh_index );
    virtual void visualize_reduced_space_indicators ( );
    virtual void visualize_reduced_time_indicators ( );

    /// *****************************************
    /// TOOLS
    /// *****************************************

    void set_problem_mode_to_primal ( )
    {
        problem_mode_ = 1;
    }

    void set_problem_mode_to_dual ( )
    {
        problem_mode_ = -1;
    }
    virtual void ApplyFilter ( VectorType& u );
    virtual void ApplyFilter ( VectorType& u, int var );
    virtual double compute_pressure_int ( VectorType& u, int var );
    virtual double compute_volume_int ( );
    int rank ( ) const;
    int master_rank ( ) const;
    int num_partitions ( ) const;
    bool use_pressure_filter ( );

    /// MPI stuff
    MPI_Comm comm_;
    int rank_;
    int num_partitions_;
    const int master_rank_;

    /// IO stuff
    std::string root_;
    std::string filename_nonlinear_solver_;
    std::string filename_linear_solver_;
    std::string filename_start_;
    std::string filename_base_;

    PropertyTree params_;
    PropertyTree base_params_;
    bool use_hdf5_;
    bool resumed_run_;

    /// Mesh stuff
    MeshPtr mesh_, master_mesh_;
    int refinement_level_;
    std::vector<MasterSlave> period_;

    /// Other stuff
    bool test_res_;
    int problem_mode_; // 1: primal, -1: dual
    int pp_step_;
    bool use_pressure_filter_;
    bool use_pressure_filter_dual_;
    bool cyl_coord_;
    int aug_p_var_;
    SYSTEM la_sys_;
    int num_vars_;
    int num_eq_;
    int update_every_newton_step_;
    int update_time_step_;
    int test_step_;
    double lp_div_tol_;
    int stationary_;
    int print_level_;
    std::vector<int> vel_var_;

    /// Primal problem
    VectorSpace<DATATYPE>* space_;
    MetFlowAssembler<DIM, DATATYPE>* local_asm_primal_;
    ScalarProdAssembler<DIM, VARDIM, DATATYPE>* local_asm_prod_;

    LAD::MatrixType* matrix_;
    LAD::VectorType* perturb_;
    LAD::VectorType* perturb_prev_;
    LAD::VectorType* solP_;
    LAD::VectorType* solP_prev_;
    LAD::VectorType* solP_next_;
    LAD::VectorType* rhs_;
    LAD::VectorType* res_;
    LAD::VectorType* base_;
    std::vector<int> dirichlet_dofs_;
    std::vector<LAD::DataType> dirichlet_values_;
    std::vector<std::vector<bool> > coupling_vars_;
    bool is_linear_problem_;
    bool matrix_assembled_;

    /// Dual problem
    VectorSpace<DATATYPE>* space_dual_;
    MetFlowAssembler<DIM, DATATYPE>* local_asm_dual_;

    LAD::MatrixType* matrix_dual_;
    LAD::VectorType* solD_;
    LAD::VectorType* solD_prev_;
    LAD::VectorType* solD_next_;
    LAD::VectorType* rhs_dual_;
    LAD::VectorType* res_dual_;
    std::vector<int> dirichlet_dofs_dual_;
    std::vector<LAD::DataType> dirichlet_values_dual_;
    std::vector<std::vector<bool> > coupling_vars_dual_;
    bool dual_matrix_assembled_;

    HpFemAssembler<double> global_asm_;
    StandardGlobalAssembler<double> pp_asm_;

    /// Error estimation
    DGGlobalAssembler<double> jump_term_asm_;
    MetFlowEstimatorAssembler<DIM, DATATYPE>* local_asm_est_;
    GoalFunctional<DIM, DATATYPE>* j_;

    std::vector< std::vector<double> > time_indicator_; // [time_step][cell_index]
    std::vector< std::vector<double> > space_indicator_; // [time_step][cell_index]
    std::vector< std::vector<double> > reduced_space_indicator_; // [mesh_index][cell_index]
    std::vector<double> reduced_time_indicator_; // [time_step]

    std::vector<double> local_equation_estimator_;
    double local_time_estimator_;
    double local_space_estimator_;
    std::vector<double> global_equation_estimator_;
    double global_time_estimator_;
    double global_space_estimator_;
    double global_estimator_;

    VectorSpace<DATATYPE>* fine_space_;
    VectorSpace<DATATYPE>* fine_space_dual_;
    LAD::VectorType fineP_;
    LAD::VectorType fineP_prev_;
    LAD::VectorType fineP_next_;
    LAD::VectorType fineD_;
    LAD::VectorType fineD_prev_;
    LAD::VectorType fineD_next_;

    std::vector<int> indicator_mesh_indices_;

#    ifdef USE_HYPRE
    SpacePatchInterpolation<LADescriptorHypreD, mesh::IMPL_P4EST> patch_interpolation_;
    SpacePatchInterpolation<LADescriptorHypreD, mesh::IMPL_P4EST> patch_interpolation_dual_;
#    else
    SpacePatchInterpolation<LADescriptorCoupledD, mesh::IMPL_P4EST> patch_interpolation_;
    SpacePatchInterpolation<LADescriptorCoupledD, mesh::IMPL_P4EST> patch_interpolation_dual_;
#    endif

    /// Linear solver
    int ilupp_time_step_;
    int ilupp_time_step_dual_;

    LinearSolver<LAD>* linear_solver_;
    LinearSolver<LAD>* linear_solver_dual_;

    /// Time discretization
    double cur_time_;
    double delta_t_;
    double delta_t_next_;
    double duration_;
    int time_step_; // Note: Theta-scheme and cGdG -> time_step refers to point in time, dGcG -> time_step refers to time interval
    double theta_;
    std::string method_;
    bool mod_galerkin_;
    int backup_step_;
    int visual_step_;
    int num_time_steps_;
    int use_azimuthal_filter_;
    int min_step_dual_;
    int max_step_dual_;
    int min_step_primal_;
    int max_step_primal_;
    double theta_dual_;
    std::string method_dual_;
    bool is_stationary_;

    /// Nonlinear solver
    MyNewton<LAD>* nls_;
    MyNewton<LAD>* nls_dual_;
    DampingStrategy<LAD>* damping_strategy_;
    ForcingStrategy<LAD>* forcing_strategy_;
    ForcingStrategy<LAD>* forcing_strategy_dual_;

    int adapt_counter_;

};
#endif
