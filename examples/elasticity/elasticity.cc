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

#include "elasticity.h"

#include <iomanip>
#include <iterator>
#include <algorithm>
#include <vector>

#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif

namespace
{
    static const char* DATADIR = MESHES_DATADIR;
    static const int MASTER_RANK = 0;
    static const char* PARAM_FILENAME = "elasticity_Bunny_Scenario_DirBC.xml"; // replaced if function main is called with argument "elasticity.xml".
    static bool CONSOLE_OUTPUT_ACTIVE = true;
    static const int CONSOLE_THRESHOLD_LEVEL = 3;
}
#define CONSOLE_OUTPUT(lvl, x) { if (CONSOLE_OUTPUT_ACTIVE &&            \
                                    lvl <= CONSOLE_THRESHOLD_LEVEL) {   \
            for (int i = 0; i < lvl; ++i) { std::cout << "  "; }        \
            std::cout << x << "\n"; }}

struct TimingData
{
    double time_elapsed;
};

class TimingScope
{
  public:

    TimingScope ( const std::string& name )
    {
        if ( report_ )
        {
            report_->begin_section ( name );
        }
    }

    TimingScope ( int iteration )
    {
        if ( report_ )
        {
            std::stringstream sstr;
            sstr << "Iteration " << iteration;
            report_->begin_section ( sstr.str ( ) );
            timer_.reset ( );
            timer_.start ( );
        }
    }

    ~TimingScope ( )
    {
        timer_.stop ( );
        if ( report_ )
        {
            TimingData* data = report_->end_section ( );
            data->time_elapsed = timer_.get_duration ( );
        }
    }

    static void set_report ( HierarchicalReport<TimingData>* report )
    {
        report_ = report;
    }

  private:
    static HierarchicalReport<TimingData>* report_;
    Timer timer_;
};

HierarchicalReport<TimingData>* TimingScope::report_ = 0;

class TimingReportOutputVisitor
{
  public:

    TimingReportOutputVisitor ( std::ostream& os )
    : os_ ( os ), level_ ( 0 )
    {
    }

    void enter ( const std::string& name, TimingData* data )
    {
        if ( name == "root" )
        {
            os_ << "+++ Timing Report +++\n\n";
        }
        else
        {
            for ( int l = 0; l < level_; ++l )
            {
                os_ << "  ";
            }
            os_ << name << " took " << data->time_elapsed << " s.\n";
            ++level_;
        }
    }

    void exit ( const std::string& name, TimingData* data )
    {
        if ( name == "root" )
        {
            os_ << "\n+++ End Timing Report +++\n\n";
        }
        else
        {
            --level_;
        }
    }

  private:
    std::ostream& os_;
    int level_;
};

class Elasticity
{
  public:

    Elasticity ( const std::string& param_filename )
    : comm_ ( MPI_COMM_WORLD ),
    params_ ( param_filename.c_str ( ), MASTER_RANK, MPI_COMM_WORLD ),
    refinement_level_ ( 0 )
    {
    }

    virtual ~Elasticity ( )
    {
    }

    virtual void run ( )
    {
        simul_name_ = params_["OutputPathAndPrefix"].get<std::string>( );
        u_deg = params_["FiniteElements"]["DisplacementDegree"].get<int>( );

        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );

        // Turn off INFO log except on master proc.
        if ( rank_ != MASTER_RANK )
        {
            INFO = false;
            CONSOLE_OUTPUT_ACTIVE = false;
        }

        std::ofstream info_log ( ( simul_name_ + "_info_log" ).c_str ( ) );
        LogKeeper::get_log ( "info" ).set_target ( &info_log );
        std::ofstream debug_log ( ( simul_name_ + "_debug_log" ).c_str ( ) );
        LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

        CONSOLE_OUTPUT ( 0, "=========================================================" );
        CONSOLE_OUTPUT ( 0, "==== Elasticity Simulation                            ===" );
        CONSOLE_OUTPUT ( 0, "====       built using HiFlow3.                       ===" );
        CONSOLE_OUTPUT ( 0, "====                                                  ===" );
        CONSOLE_OUTPUT ( 0, "==== Engineering Mathematics and Computing Lab (EMCL) ===" );
        CONSOLE_OUTPUT ( 0, "=========================================================" );
        CONSOLE_OUTPUT ( 0, "" );

        // output parameters for debugging
        LOG_INFO ( "parameters", params_ );

        // setup timing report
        TimingScope::set_report ( &time_report_ );

        {
            TimingScope tscope ( "Setup" );

            setup_linear_algebra ( );

            read_mesh ( );

            // The simulation has two modes: stationary and
            // instationary. Which one is used, depends on the parameter
            // Instationary.SolveInstationary. In stationary mode, solve()
            // is only called once, whereas in instationary mode, it is
            // called several times via the time-stepping method implemented in run_time_loop() .
            solve_instationary_ = params_["Instationary"]["SolveInstationary"].get<bool>( );

            prepare ( );

        }

        if ( solve_instationary_ )
        {
            TimingScope tscope ( "Complete_instationary_simulation_loop" );
            // See private source code by Nicolai Schoch.

        }
        else
        {
            TimingScope tscope ( "Complete_stationary_simulation_loop" );

            LOG_INFO ( "simulation", "Solving stationary problem" );

            assemble_system ( );

            solve_system ( );

            visualize ( );
        }

        CONSOLE_OUTPUT ( 0, "" );

        if ( rank_ == MASTER_RANK )
        {
            // Output time report
            TimingReportOutputVisitor visitor ( std::cout );
            time_report_.traverse_depth_first ( visitor );
        }

        LogKeeper::get_log ( "info" ).flush ( );
        LogKeeper::get_log ( "debug" ).flush ( );
        LogKeeper::get_log ( "info" ).set_target ( 0 );
        LogKeeper::get_log ( "debug" ).set_target ( 0 );

        CONSOLE_OUTPUT ( 0, "============================================================\n" );
    }

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

    PLATFORM la_platform ( ) const
    {
        return la_sys_.Platform;
    }

    IMPLEMENTATION la_implementation ( ) const
    {
        return la_impl_;
    }

    MATRIX_FORMAT la_matrix_format ( ) const
    {
        return la_matrix_format_;
    }

    // Read, refine and partition mesh.
    void read_mesh ( );

    // Set up datastructures and read in some parameters.
    void prepare ( );

    // Set up boundary conditions
    void prepare_bc_stationary ( );

    // Assemble system (includes solve_system() for the instationary case).
    void assemble_system ( );

    // Solve System (excludes solver procedure for the instationary case).
    void solve_system ( );

    // Visualize the solution in a file.
    // In stationary mode, the filename contains 'stationary',
    // in instationary mode, it contains the current time-step ts_.
    void visualize ( /*ts_*/ );

    // Linear algebra set up
    void setup_linear_algebra ( );

    // MPI stuff
    MPI_Comm comm_;
    int rank_, num_partitions_;

    // Linear algebra stuff
    SYSTEM la_sys_;
    IMPLEMENTATION la_impl_;
    MATRIX_FORMAT la_matrix_format_;

    // Parameter data read in from file.
    PropertyTree params_;

    std::string simul_name_; // parameter 'OutputPathAndPrefix': prefix for output files
    int u_deg; // Finite Element Degree (of variable u)

    // Elasticity model variables
    double lambda_, mu_, rho_, gravity_;

    // Meshes
    MeshPtr mesh_;
    int refinement_level_;

    VectorSpace<double> space_;

    // linear algebra objects
    Couplings<double> couplings_;
    CMatrix matrix_; // Moreover: stiffM_, dampM_, massM_, iterM_;
    CVector sol_, rhs_, nbc_, cbc_; // solution-vector u, rhs-vector f, nbc-vector s, cbc-vector c;

    // linear solver parameters
    std::string solver_name_;
    int lin_max_iter;
    double lin_abs_tol;
    double lin_rel_tol;
    double lin_div_tol;
    int basis_size;

    // preconditioner parameters
    MATRIX_FREE_PRECOND matrix_precond_; // preconditioner name
    double omega_; // for SOR and SSOR
    int ilu_p_; // for ILUp

    // preconditioner
    PreconditionerBlockJacobiStand<LAD> preconditioner_; // Base class for all block Jacobi preconditioners
    // PreconditionerBlockJacobi: base class for all block Jacobi preconditioners (Jacobi, GaussSeidel, S-GaussSeidel, SOR, SSOR, ...).
#ifdef WITH_ILUPP
    PreconditionerIlupp<LAD> ilupp_;
#endif

    // linear solver
    // instead of putting up two solvers, alternatively put up an abstract LinearSolver-Object which is filled during runtime via XML-Input.
    GMRES<LAD> gmres_; // GMRES works best with ILUpp.
    CG<LAD> cg_; // CG works best with SymmetricGaussSeidel.

    StandardGlobalAssembler<double> global_asm_;

    bool solve_instationary_, use_ilupp_;

    // Dirichlet BC vectors (DoFs and Values).
    std::vector<int> dirichlet_dofs_;
    std::vector<Scalar> dirichlet_values_;
    std::vector<int> fixed_dirichlet_dofs_;
    std::vector<Scalar> fixed_dirichlet_values_;
    std::vector<int> displacement_dirichlet_dofs_;
    std::vector<Scalar> displacement_dirichlet_values_;

    // Visualization vector.
    std::vector<std::string> visu_names_;

    HierarchicalReport<TimingData> time_report_;
};

// ------------------------------------------------------------------

// program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    std::string param_filename ( PARAM_FILENAME );
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }

    try
    {
        Elasticity app ( param_filename );
        app.run ( );
    }
    catch ( std::exception& e )
    {
        std::cerr << "\nProgram ended with uncaught exception.\n";
        std::cerr << e.what ( ) << "\n";
        return -1;
    }
    MPI_Finalize ( );

    return 0;
}

// ------------------------------------------------------------------

void Elasticity::read_mesh ( )
{
    TimingScope tscope ( "read_mesh" );

    MeshPtr master_mesh;
    MeshPtr local_mesh;

    // get mesh_filename:
    const std::string mesh_name = params_["Mesh"]["Filename"].get<std::string>( );
    std::string mesh_filename = std::string ( DATADIR ) + mesh_name;

    // find position of last '.' in mesh_filename:
    int dot_pos = mesh_filename.find_last_of ( "." );
    assert ( dot_pos != std::string::npos );

    // get suffix of mesh_filename:
    std::string suffix = mesh_filename.substr ( dot_pos );

    // switch cases according to the suffix of mesh_filename --> pvtu or vtu/inp:
    // this tutorial allows for ".vtu" or ".inp" only.

    if ( rank ( ) == MASTER_RANK )
    {

        master_mesh = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, &comm_ );
        CONSOLE_OUTPUT ( 1, "Read mesh with " << master_mesh->num_entities ( DIMENSION ) << " cells." );

        if ( suffix == std::string ( ".vtu" ) )
        {
            // in case the input mesh does not specify material_IDs for the boundary facets,
            // the function set_default_material_number_on_bdy(MeshPtr*, int) sets default bdy_facet_IDs
            // in order for (among others) the GridGeometricSearch to work properly
            // especially w.r.t. find_closest_point() -> see also: documentation.
            set_default_material_number_on_bdy ( master_mesh, 1 );
        }
        refinement_level_ = 0;
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );

        for ( int r = 0; r < initial_ref_lvl; ++r )
        {
            master_mesh = master_mesh->refine ( );
            ++refinement_level_;
        }
        LOG_INFO ( "mesh", "Initial refinement level = " << refinement_level_ );
        CONSOLE_OUTPUT ( 1, "Refined mesh (level " << refinement_level_ << ") has "
                         << master_mesh->num_entities ( DIMENSION ) << " cells." );

    }

    if ( num_partitions ( ) > 1 )
    {
#ifdef WITH_METIS
        MetisGraphPartitioner partitioner;
#else
        NaiveGraphPartitioner partitioner;
#endif

        const GraphPartitioner* p = &partitioner;

        local_mesh = partition_and_distribute ( master_mesh, MASTER_RANK, comm_, p );
    }
    else
    {
        NaiveGraphPartitioner partitioner;
        local_mesh = partition_and_distribute ( master_mesh, MASTER_RANK, comm_, &partitioner );
    }
    assert ( local_mesh != 0 );

    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );

    std::ostringstream rank_str;
    rank_str << rank ( );

    PVtkWriter writer ( comm_ );
    std::string output_file = simul_name_ + std::string ( "_initial_mesh_local.pvtu" );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

void Elasticity::prepare ( )
{
    TimingScope tscope ( "prepare" );

    // prepare modelling problem parameters
    rho_ = params_["ElasticityModel"]["density"].get<double>( );
    mu_ = params_["ElasticityModel"]["mu"].get<double>( );
    lambda_ = params_["ElasticityModel"]["lambda"].get<double>( );
    gravity_ = params_["ElasticityModel"]["gravity"].get<double>( );

    // prepare space
    std::vector< int > degrees ( DIMENSION );
    for ( int c = 0; c < DIMENSION; ++c )
    {
        degrees.at ( c ) = u_deg;
    }

    // Initialize the VectorSpace object
    space_.Init ( degrees, *mesh_ );

    CONSOLE_OUTPUT ( 1, "Total number of dofs = " << space_.dof ( ).ndofs_global ( ) );

    for ( int p = 0; p < num_partitions ( ); ++p )
    {
        CONSOLE_OUTPUT ( 2, "Num dofs on process " << p << " = " << space_.dof ( ).ndofs_on_sd ( p ) );
    }

    // prepare visualization structures
    visu_names_.push_back ( "u1" );
    visu_names_.push_back ( "u2" );
    visu_names_.push_back ( "u3" );

    // Setup couplings object and prepare linear algebra structures
    couplings_.Clear ( );
    couplings_.Init ( communicator ( ), space_.dof ( ) );

    // prepare global assembler
    QuadratureSelection q_sel ( params_["QuadratureOrder"].get<int>( ) );
    global_asm_.set_quadrature_selection_function ( q_sel );

    // compute matrix graph
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows, sparsity.off_diagonal_cols );

    // Initialize system matrix (K_eff), solution vector and rhs vector (R_eff) [and possibly nbc vector].
    matrix_.Init ( communicator ( ), couplings_, la_platform ( ), la_implementation ( ), la_matrix_format ( ) ); // System (Stiffness) Matrix.
    sol_.Init ( communicator ( ), couplings_, la_platform ( ), la_implementation ( ) ); // Solution Vector.
    rhs_.Init ( communicator ( ), couplings_, la_platform ( ), la_implementation ( ) ); // RHS Vector.

    matrix_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                            vec2ptr ( sparsity.diagonal_cols ),
                            sparsity.diagonal_rows.size ( ),
                            vec2ptr ( sparsity.off_diagonal_rows ),
                            vec2ptr ( sparsity.off_diagonal_cols ),
                            sparsity.off_diagonal_rows.size ( ) );
    matrix_.Zeros ( );

    sol_.InitStructure ( );
    sol_.Zeros ( );

    rhs_.InitStructure ( );
    rhs_.Zeros ( );

    // setup linear solver parameters
    solver_name_ = params_["LinearSolver"]["SolverName"].get<std::string>( );
    lin_max_iter = params_["LinearSolver"]["MaximumIterations"].get<int>( );
    lin_abs_tol = params_["LinearSolver"]["AbsoluteTolerance"].get<double>( );
    lin_rel_tol = params_["LinearSolver"]["RelativeTolerance"].get<double>( );
    lin_div_tol = params_["LinearSolver"]["DivergenceLimit"].get<double>( );
    basis_size = params_["LinearSolver"]["BasisSize"].get<int>( );

    // Setup Solver and Preconditioner //////////////////////////////////////////////////////////
    if ( solver_name_ == "CG" )
    {

        // Setup CG Preconditioner //////////////////////////////////////////
        if ( matrix_precond_ != NOPRECOND )
        { // With Preconditioning.

            if ( matrix_precond_ == JACOBI )
            {
                preconditioner_.Init_Jacobi ( rhs_ );
            }
            if ( matrix_precond_ == GAUSS_SEIDEL )
            {
                preconditioner_.Init_GaussSeidel ( );
            }
            if ( matrix_precond_ == SGAUSS_SEIDEL )
            {
                preconditioner_.Init_SymmetricGaussSeidel ( );
            }
            if ( matrix_precond_ == SOR )
            {
                omega_ = params_["LinearSolver"]["Omega"].get<double>( );
                preconditioner_.Init_SOR ( omega_ );
            }
            if ( matrix_precond_ == SSOR )
            {
                omega_ = params_["LinearSolver"]["Omega"].get<double>( );
                preconditioner_.Init_SSOR ( omega_ );
            }
            if ( matrix_precond_ == ILU )
            {
                std::cout << "By definition: CG is not supposed to pair with ILU as solver-preconditioner-combination." << std::endl;
            }
            if ( matrix_precond_ == ILU2 )
            { // Note that: ILU2 = ILUpp; handled for GMRES only.
                std::cout << "By definition: CG is not supposed to pair with ILU2 (a.k.a. ILUpp) as solver-preconditioner-combination." << std::endl;
            }

            cg_.InitParameter ( "Preconditioning" );
            cg_.SetupPreconditioner ( preconditioner_ );

        }
        else
        { // Without Preconditioning.

            cg_.InitParameter ( "NoPreconditioning" );

        }

        // Setup CG Solver ////////////////////////////////////////////////
        cg_.InitControl ( lin_max_iter, lin_abs_tol, lin_rel_tol, lin_div_tol );
        cg_.SetupOperator ( matrix_ );

    } // end of CG Setup.
    else /*if (solver_name_ == "GMRES")*/
    { // if not CG then GMRES by default.

        // Setup GMRES Preconditioner //////////////////////////////////////////
        use_ilupp_ = 0; // default.
        if ( matrix_precond_ == ILU2 )
        { // Note that: ILU2 = ILU_pp.
            use_ilupp_ = 1; // otherwise: no preconditioner for GMRES.
        }
#ifdef WITH_ILUPP
        if ( use_ilupp_ )
        {
            ilupp_.InitParameter ( params_["ILUPP"]["PreprocessingType"].get<int>( ),
                                   params_["ILUPP"]["PreconditionerNumber"].get<int>( ),
                                   params_["ILUPP"]["MaxMultilevels"].get<int>( ),
                                   params_["ILUPP"]["MemFactor"].get<double>( ),
                                   params_["ILUPP"]["PivotThreshold"].get<double>( ),
                                   params_["ILUPP"]["MinPivot"].get<double>( ) );

            gmres_.SetupPreconditioner ( ilupp_ );
            gmres_.InitParameter ( basis_size, "RightPreconditioning" );
        }
        else
        {
            gmres_.InitParameter ( basis_size, "NoPreconditioning" );
        }
#else
        gmres_.InitParameter ( basis_size, "NoPreconditioning" );
#endif

        // Setup GMRES Solver ////////////////////////////////////////////////
        gmres_.InitControl ( lin_max_iter, lin_abs_tol, lin_rel_tol, lin_div_tol );
        gmres_.SetupOperator ( matrix_ );

    } // end of GMRES Setup.

    // prepare dirichlet BC
    if ( !solve_instationary_ )
    {
        prepare_bc_stationary ( );
    }

}

void Elasticity::prepare_bc_stationary ( )
{
    TimingScope tscope ( "prepare_bc_stationary" );

    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    // >-----------------------------------------------------------------------------------------------------
    // Set DirichletBCs for prescribed facets (i.e. material_IDs given in hiflow3_scene.xml-File):
    const int dir_bdy1 = params_["Boundary"]["DirichletMaterial1"].get<int>( );
    const int dir_bdy2 = params_["Boundary"]["DirichletMaterial2"].get<int>( );
    const int dir_bdy3 = params_["Boundary"]["DirichletMaterial3"].get<int>( );

    // create InstationaryElasticity_DirichletBC_3D-Object.
    StationaryElasticity_DirichletBC_3D bc[3] = { StationaryElasticity_DirichletBC_3D ( 0, dir_bdy1, dir_bdy2, dir_bdy3 ),
                                                 StationaryElasticity_DirichletBC_3D ( 1, dir_bdy1, dir_bdy2, dir_bdy3 ),
                                                 StationaryElasticity_DirichletBC_3D ( 2, dir_bdy1, dir_bdy2, dir_bdy3 ) };

    // and compute Dirichlet values for pre-set dofs.
    for ( int var = 0; var < DIMENSION; ++var )
    {
        compute_dirichlet_dofs_and_values ( bc[var], space_, var,
                                            dirichlet_dofs_, dirichlet_values_ );
    }
}

void Elasticity::assemble_system ( )
{

    if ( !solve_instationary_ )
    {
        TimingScope tscope ( "Assemble_system_stationary" );

        StationaryElasticityAssembler local_asm ( lambda_, mu_, rho_, gravity_ );

        global_asm_.assemble_matrix ( space_, local_asm, matrix_ );

        global_asm_.assemble_vector ( space_, local_asm, rhs_ );

        rhs_.UpdateCouplings ( );

        if ( !dirichlet_dofs_.empty ( ) )
        {
            // Correct Dirichlet dofs.
            matrix_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
            rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                             vec2ptr ( dirichlet_values_ ) );
            sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                             vec2ptr ( dirichlet_values_ ) );
        }

        rhs_.UpdateCouplings ( );
        sol_.UpdateCouplings ( );

        // Setup Operator for Preconditioner/Solver:
        // (needs to be done AFTER DirichletBC-Setup)
        if ( solver_name_ == "CG" && matrix_precond_ != NOPRECOND )
        {
            preconditioner_.SetupOperator ( matrix_ );
        }
        else if ( solver_name_ == "CG" && matrix_precond_ == NOPRECOND )
        {
            cg_.SetupOperator ( matrix_ );
        }
        else
        {
            // GMRES (and ILUPP) by default.
            if ( use_ilupp_ )
            {
#ifdef WITH_ILUPP
                ilupp_.SetupOperator ( matrix_ );
#endif
            }
        }

    }

}

void Elasticity::solve_system ( )
{

    // Solver setup is done in prepare_system() and assemble_system().

    if ( !solve_instationary_ )
    {
        TimingScope tscope ( "Solve_system_stationary" );

        // Solve linear system.
        if ( solver_name_ == "CG" )
        {
            cg_.Solve ( rhs_, &sol_ );
        }
        else /*if (solver_name_ == "GMRES")*/
        {
            gmres_.Solve ( rhs_, &sol_ );
        }

        sol_.UpdateCouplings ( );

    }

    if ( solver_name_ == "CG" )
    {
        if ( matrix_precond_ == 0 )
        {
            CONSOLE_OUTPUT ( 3, "Linear solver (CG) not using preconditioner." );
        }
        else
        {
            CONSOLE_OUTPUT ( 3, "Linear solver (CG) using the following preconditioner: " << matrix_precond_ << "." );
            CONSOLE_OUTPUT ( 3, "Note: Jacobi = 1, GaussSeidel = 2, SymmetricGaussSeidel = 3, SOR = 4, SSOR = 5, ILU = 6, ILU2 = ILU_pp = 7, ILU_P = 8." );
        }
        CONSOLE_OUTPUT ( 3, "Linear solver (CG) computed solution in " << cg_.iter ( ) << " iterations." );
        CONSOLE_OUTPUT ( 3, "Residual norm for solution = " << cg_.res ( ) );
    }
    else
    {
        if ( matrix_precond_ == 0 )
        {
            CONSOLE_OUTPUT ( 3, "Linear solver (GMRES) not using preconditioner." );
        }
        else
        {
            CONSOLE_OUTPUT ( 3, "Linear solver (GMRES) using the following preconditioner: " << matrix_precond_ << "." );
            CONSOLE_OUTPUT ( 3, "Note: Jacobi = 1, GaussSeidel = 2, SymmetricGaussSeidel = 3, SOR = 4, SSOR = 5, ILU = 6, ILU2 = ILU_pp = 7, ILU_P = 8." );
        }
        CONSOLE_OUTPUT ( 3, "Linear solver (GMRES) computed solution in " << gmres_.iter ( ) << " iterations." );
        CONSOLE_OUTPUT ( 3, "Residual norm for solution = " << gmres_.res ( ) );
    }
}

void Elasticity::visualize ( /* int ts_*/ )
{
    TimingScope tscope ( "Visualization" );

    // Setup visualization object.
    int num_intervals = u_deg; // num_intervals for CellVisualization is reasonably less or equal the FE degree.
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    // Generate filename for mesh which includes solution-array.
    std::stringstream input;
    input << simul_name_ << "_solution";
    const int xml_init_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );

    input << "_np" << num_partitions ( ) << "_RefLvl" << xml_init_ref_lvl << "_Tstep";
    input << "_stationary";

    if ( num_partitions ( ) > 1 )
        input << ".pvtu";
    else
        input << ".vtu";

    std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> material_number ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) ); it != mesh_->end ( mesh_->tdim ( ) ); ++it )
    {
        int temp1, temp2;
        mesh_->get_attribute_value ( "_remote_index_", mesh_->tdim ( ),
                                     it->index ( ),
                                     &temp1 );
        mesh_->get_attribute_value ( "_sub_domain_", mesh_->tdim ( ),
                                     it->index ( ),
                                     &temp2 );
        remote_index.at ( it->index ( ) ) = temp1;
        sub_domain.at ( it->index ( ) ) = temp2;
        material_number.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
    }

    sol_.UpdateCouplings ( );

    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 0 ), "u0" );
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 1 ), "u1" );
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_, 2 ), "u2" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );

    // Create visualization helper vector from (post-processing/solution) vector (pp_sol_/sol_).
    std::vector<LAD::DataType> visu_vec ( sol_.size_global ( ), 1.e20 );
    std::vector<int> dof_ids;
    std::vector<double> values;
    sol_.GetAllDofsAndValues ( dof_ids, values ); // solution values (u1,u2,u3) from "sol_" are given into the "visu" object.
    for ( int i = 0; i < static_cast < int > ( values.size ( ) ); ++i )
    {
        visu_vec.at ( dof_ids[i] ) = values.at ( i ); // values.at(i) correspond to displacement values (as in (pp_)sol_).
        // Note: visu_vec does not contain the body's deformed state
        // (the deformed state would correspond to orig values + displacement values)
    }

    // ---------------------------------------------------------
    // START THE ADDITIONAL DEFORMATION VISUALIZATION PART HERE:
    // Compute "deformed state"
    // as sum of initial coords (= mesh input) and displacement values (= solution vector u).

    std::vector<double> deformation_vector ( 3 * mesh_->num_entities ( 0 ) ); // This vector contains the mesh/geometry in deformed state.

    std::vector<int> local_dof_id;

    AttributePtr sub_domain_attr;
    // if attribute _sub_domain_ exists, there are several subdomains, and hence ghost cells, too.
    const bool has_ghost_cells = mesh_->has_attribute ( "_sub_domain_", mesh_->tdim ( ) );
    if ( has_ghost_cells )
    {
        sub_domain_attr = mesh_->get_attribute ( "_sub_domain_", mesh_->tdim ( ) );
    }

    for ( EntityIterator cell = mesh_->begin ( mesh_->tdim ( ) ), end_c = mesh_->end ( mesh_->tdim ( ) ); cell != end_c; ++cell )
    {

        int vertex_num = 0;

        if ( has_ghost_cells && sub_domain_attr->get_int_value ( cell->index ( ) ) != rank_ )
        {
            continue; // in order not to write values twice/thrice/... for parallel processes
        }

        for ( IncidentEntityIterator vertex_it = cell->begin_incident ( 0 ); vertex_it != cell->end_incident ( 0 ); ++vertex_it )
        {
            //for (int vertex = 0; vertex < cell->num_vertices(); ++vertex) {
            for ( int var = 0; var < 3; ++var )
            {
                space_.dof ( ).get_dofs_on_subentity ( var, cell->index ( ), 0, vertex_num, local_dof_id );
                deformation_vector.at ( 3 * vertex_it->index ( ) + var ) = visu_vec.at ( local_dof_id[0] );
            }
            ++vertex_num;
        }
    }

    // Generate output filename and write output mesh/geometry of deformed state.
    std::stringstream input2;
    input2 << simul_name_ << "_deformedSolution";
    input2 << "_np" << num_partitions ( ) << "_RefLvl" << xml_init_ref_lvl << "_Tstep";
    input2 << "_stationary";

    if ( num_partitions ( ) == 1 )
    {
        input2 << ".pvtu";
        PVtkWriter deformation_writer ( comm_ );
        deformation_writer.set_deformation ( &deformation_vector );
        deformation_writer.add_all_attributes ( *mesh_, true );
        deformation_writer.write ( input2.str ( ).c_str ( ), *mesh_ );

    }
    else if ( num_partitions ( ) > 1 )
    {

        mesh::MeshPtr mesh_ptr;
        mesh::TDim tdim = mesh_->tdim ( );

        if ( mesh_->has_attribute ( "_sub_domain_", tdim ) )
        {
            // create new mesh without ghost layers
            mesh::MeshDbViewBuilder builder ( ( static_cast < mesh::MeshDbView* > ( mesh_.get ( ) ) )->get_db ( ) );

            // cell index map
            std::vector<int> index_map;

            // get __remote_indices__ attribute
            mesh::AttributePtr subdomain_attr = mesh_->get_attribute ( "_sub_domain_", tdim );

            // loop over cells that are not ghosts
            for ( mesh::EntityIterator cell_it = mesh_->begin ( tdim ); cell_it != mesh_->end ( tdim ); ++cell_it )
            {
                if ( subdomain_attr->get_int_value ( cell_it->index ( ) ) == rank_ )
                {
                    // loop over vertices of cell
                    std::vector<mesh::MeshBuilder::VertexHandle> vertex_handle;
                    for ( mesh::IncidentEntityIterator inc_vert_it = cell_it->begin_incident ( 0 ); inc_vert_it != cell_it->end_incident ( 0 ); ++inc_vert_it )
                    {
                        std::vector<double> coord ( 3, 0.0 );
                        inc_vert_it->get_coordinates ( coord );
                        vertex_handle.push_back ( builder.add_vertex ( coord ) );
                    }
                    builder.add_entity ( tdim, vertex_handle );
                    index_map.push_back ( cell_it->index ( ) );
                }
            }

            mesh_ptr = builder.build ( );
            mesh::AttributePtr index_map_attr ( new mesh::IntAttribute ( index_map ) );
            mesh_ptr->add_attribute ( "__index_map__", mesh_ptr->tdim ( ), index_map_attr );

            // transfer cell attributes
            std::vector< std::string > cell_attr_names = mesh_->get_attribute_names ( mesh_ptr->tdim ( ) );
            for ( std::vector< std::string >::const_iterator it = cell_attr_names.begin ( ),
                  end_it = cell_attr_names.end ( ); it != end_it; ++it )
            {
                mesh::AttributePtr mapped_attr (
                                                 new mesh::InheritedAttribute ( mesh_->get_attribute ( *it, mesh_->tdim ( ) ), index_map_attr ) );
                mesh_ptr->add_attribute ( *it, mesh_->tdim ( ), mapped_attr );
            }
        }
        else
        {
            std::cout << "Error in Visualization Routine! Missing subdomain attribute!\n";
        }

        // Since deformation_vector still contains some vertices twice/thrice/... (according to num processes, and ghost cells),
        // every vertex is now (by means of an EntityIterator) set/mapped (once only) into a mapped_deformation_vector.
        // Thus, finally, the mapped_deformation_vector has the size 3*mesh_ptr->num_entities(0).
        std::vector<double> mapped_deformation_vector ( 3 * mesh_ptr->num_entities ( 0 ) ); // dim = number of nodes in mesh * dimension
        for ( EntityIterator v_it = mesh_ptr->begin ( 0 ), end_v = mesh_ptr->end ( 0 ); v_it != end_v; ++v_it )
        {
            mesh::Id v_id = v_it->id ( );
            mesh::EntityNumber v_index = -1;
            bool found = mesh_->find_entity ( 0, v_id, &v_index );
            if ( !found )
            {
                std::cout << "Error in Visualization Routine! Missing vertex!\n";
            }
            for ( int c = 0; c < 3; ++c )
            {
                mapped_deformation_vector[3 * v_it->index ( ) + c] = deformation_vector[3 * v_index + c];
            }
        }

        input2 << ".pvtu";
        PVtkWriter p_deformation_writer ( comm_ );
        p_deformation_writer.set_deformation ( &mapped_deformation_vector );
        p_deformation_writer.add_all_attributes ( *mesh_ptr, true );
        p_deformation_writer.write ( input2.str ( ).c_str ( ), *mesh_ptr );

    }
    else
    {
        std::cout << "Error in Visualization Routine! Wrong counting and managing the number of partitions: num_partitions_ <= 0 is wrong." << std::endl;
    }
}

void Elasticity::setup_linear_algebra ( )
{

    TimingScope tscope ( "setup_linear_algebra" );

    const std::string platform_str = params_["LinearAlgebra"]["Platform"].get<std::string>( );
    if ( platform_str == "CPU" )
    {
        la_sys_.Platform = CPU;
    }
    else if ( platform_str == "GPU" )
    {
        la_sys_.Platform = GPU;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.Platform", platform_str );
    }
    init_platform ( la_sys_ );

    const std::string impl_str = params_["LinearAlgebra"]["Implementation"].get<std::string>( );
    if ( impl_str == "Naive" )
    {
        la_impl_ = NAIVE;
    }
    else if ( impl_str == "BLAS" )
    {
        la_impl_ = BLAS;
    }
    else if ( impl_str == "MKL" )
    {
        la_impl_ = MKL;
    }
    else if ( impl_str == "OPENMP" )
    {
        la_impl_ = OPENMP;
    }
    else if ( impl_str == "SCALAR" )
    {
        la_impl_ = SCALAR;
    }
    else if ( impl_str == "SCALAR_TEX" )
    {
        la_impl_ = SCALAR_TEX;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.Implementation", impl_str );
    }

    const std::string matrix_str = params_["LinearAlgebra"]["MatrixFormat"].get<std::string>( );
    if ( matrix_str == "CSR" )
    {
        la_matrix_format_ = CSR;
    }
    else if ( matrix_str == "COO" )
    {
        la_matrix_format_ = COO;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.MatrixFormat", impl_str );
    }

    // the following part is needed to initialize class member "matrix_precond_".
    const std::string precond_str = params_["LinearSolver"]["PreconditionerName"].get<std::string>( );
    if ( precond_str == "NOPRECOND" )
    {
        matrix_precond_ = NOPRECOND;
    }
    else if ( precond_str == "JACOBI" )
    {
        matrix_precond_ = JACOBI;
    }
    else if ( precond_str == "GAUSS_SEIDEL" )
    {
        matrix_precond_ = GAUSS_SEIDEL;
    }
    else if ( precond_str == "SGAUSS_SEIDEL" )
    {
        matrix_precond_ = SGAUSS_SEIDEL;
    }
    else if ( precond_str == "SOR" )
    {
        matrix_precond_ = SOR;
    }
    else if ( precond_str == "SSOR" )
    {
        matrix_precond_ = SSOR;
    }
    else if ( precond_str == "ILU" )
    {
        matrix_precond_ = ILU;
    }
    else if ( precond_str == "ILU2" )
    {
        matrix_precond_ = ILU2; // ILU2 = ILUpp;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearSolver.PreconditionerName", precond_str );
    }
}
