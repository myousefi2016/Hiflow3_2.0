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

/// \author Katrin Mang<br>Simon Gawlok<br>Philipp Gerstner

#include "poisson_adaptive.h"

const int DEBUG_LEVEL = 1;
const int PRINT_PROC = 1;
static const char* PARAM_FILENAME = "poisson_adaptive.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class PoissonAdaptive
{
  public:

    PoissonAdaptive ( const std::string& param_filename,
                      const std::string& path_mesh )
    : path_mesh ( path_mesh ),
    comm_ ( MPI_COMM_WORLD ),
    rank_ ( -1 ),
    num_partitions_ ( -1 ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    is_done_ ( false ),
    refinement_level_ ( 0 )
#ifndef WITH_HYPRE
    ,
    matrix_ ( new CoupledMatrix<Scalar> ),
    rhs_ ( new CoupledVector<Scalar> ),
    sol_ ( new CoupledVector<Scalar> ),
    exact_sol_ ( new CoupledVector<Scalar> )
#endif
    {
        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );

#ifdef ELLIPSOID_BOUNDARY
        if ( num_partitions_ > 1 )
        {
            std::cerr << "Adaptive poisson example with ellipsoid boundary can only be run sequentially!"
                    << std::endl;
            exit ( -1 );
        }
        if ( DIMENSION != 2 )
        {
            std::cerr << "Adaptive poisson example with ellipsoid boundary can only be run in 2D!"
                    << std::endl;
            exit ( -1 );
        }
#endif
#ifndef USE_MESH_P4EST
        if ( num_partitions_ > 1 )
        {
            if ( params_["Mesh"]["Refinement"].get<int>( ) != 1 )
            {
                std::cerr << "Adaptive poisson example with MeshDbView implementation can only be run sequentially!"
                        << std::endl;
                exit ( -1 );
            }
        }
#endif
    }

    // Main algorithm

    void run ( )
    {
        int patch_mode_level = params_["Mesh"]["PatchModeLevel"].get<int>( -1 );

        // Prepare headers of CSV file
        this->csv_names_.clear ( );
        this->csv_names_.push_back ( "Level" );
        this->csv_names_.push_back ( "Global number of DoFs " );
        this->csv_names_.push_back ( "Setup system time " );
        this->csv_names_.push_back ( "Assembly time " );
        this->csv_names_.push_back ( "Solve time " );
        this->csv_names_.push_back ( "Compute Error time " );
        this->csv_names_.push_back ( "L2 error " );
        this->csv_names_.push_back ( "H1 error " );
        this->csv_names_.push_back ( "Error estimator time " );
        this->csv_names_.push_back ( "Std estimator " );
        this->csv_names_.push_back ( "Patch interpolation time " );
        this->csv_names_.push_back ( "Mesh Adapt time " );
        this->csv_names_.push_back ( "Mesh Ghost time " );

        if ( this->rank_ == MASTER_RANK )
        {
            std::stringstream csv_file_name;
            csv_file_name << this->num_partitions_ << "_statistics.csv";
            this->csv_writer_.InitFilename ( csv_file_name.str ( ) );
            this->csv_writer_.Init ( this->csv_names_ );
        }

        // Construct / read in the initial mesh.
        build_initial_mesh ( );

        // Main adaptation loop.
        while ( !is_done_ )
        {
            // Prepare data for statistics in CSV file
            this->csv_quantities_.clear ( );

            // Add current refinement level
            this->csv_quantities_.push_back ( this->refinement_level_ );

            LOG_INFO ( "Refinement Level", " ===== " << refinement_level_ << "\n=====" );

            // Initialize space and linear algebra.
            prepare_system ( );

            // Compute the stiffness matrix and right-hand side.
            assemble_system ( );

            // Solve the linear system.
            solve_system ( );

            // Compute the error to the exact solution.
            // Includes the A posteriori error estimator.
            compute_error ( );
            error_estimator ( );

            // Visualize the solution and the errors.
            visualize ( );

#ifdef USE_MESH_P4EST
            // higher order patch interpolation
            if ( refinement_level_ > patch_mode_level && patch_mode_level >= 0 )
            {
                LOG_INFO ( "Higher order interpolation", true );
                higher_order_interpolation ( );
            }
            else
            {
                this->csv_quantities_.push_back ( 0.0 );
            }
#endif
            LOG_INFO ( "Visualize", true );
            // Uniform or adaptive refinement
            int refinement_ = params_["Mesh"]["Refinement"].get<int>( 3 );
            if ( refinement_ == 1 )
            {
                adapt_uniform ( );
                LOG_INFO ( "Adapt uniform", true );
            }
            else
            {
                adapt ( );
                LOG_INFO ( "Adapt locally", true );
                LOG_INFO ( "Is done", is_done_ );
            }

            // Write CSV output
            if ( this->rank_ == MASTER_RANK )
            {
                this->csv_writer_.write ( this->csv_quantities_ );
            }
        }
    }

    ~PoissonAdaptive ( )
    {
#ifndef WITH_HYPRE
        delete matrix_;
        delete sol_;
        delete rhs_;
#endif
    }

  private:
    // Member functions

    // Read and distribute mesh.
    std::string path_mesh;
    void build_initial_mesh ( );
    // Setup space, linear algebra, and compute Dirichlet values.
    void prepare_system ( );

    // Compute the matrix and rhs.
    void assemble_system ( );
    // Compute solution x.
    void solve_system ( );

    // Compute errors compared to exact solution.
    // Includes error estimators
    void compute_error ( );
    void error_estimator ( );
    // Visualize the results.
    void visualize ( );
    // Adapat the mesh global uniform
    void adapt_uniform ( );
    // Adapt the space (mesh and/or degree).
    void adapt ( );

    void higher_order_interpolation ( );

    // method of adaptive refinement
    int refinement_;
    // kind of estimator: standard or new constructed
    int estimator_;

    // member variables
    // MPI communicator.
    MPI_Comm comm_;
    // Local process rank and number of processes.
    int rank_, num_partitions_;

    // Parameter data read in from file.
    PropertyTree params_;

    // Local mesh and mesh on master process.
    MeshPtr mesh_, master_mesh_, mesh_without_ghost_;
    // Solution space.
    VectorSpace<double> space_;

    // Linear algebra couplings helper object.
    Couplings<double> couplings_;
#ifdef WITH_HYPRE
    HypreMatrix<double> matrix_;
    HypreVector<double> rhs_, sol_, exact_sol_;
#else
    // Vectors for solution and load vector.
    CoupledVector<Scalar>* rhs_, *sol_, *exact_sol_;
    // System matrix.
    CoupledMatrix<Scalar>* matrix_;
#endif

    // Global assembler.
    // HpFem instead of StandardAssembler to allow different polynomial degrees.
    HpFemAssembler<double> global_asm_;
    DGGlobalAssembler<double> jump_term_asm_;

    // Flag for stopping adaptive loop.
    bool is_done_;
    // Current refinement level.
    int refinement_level_;

    // Vectors for error norms and error estimators
    std::vector<Scalar> L2_err_, H1_err_;
    std::vector<Scalar> rho_cell_, rho_bcell_, rho_jump_, rho_boundary_;
    std::vector<Scalar> std_estimator_, new_estimator_;
    double global_std_estimator_, global_new_estimator_;
    int global_size_std_, global_size_new_;

    double relation;
    // Dof id:s for Dirichlet boundary conditions.
    std::vector<int> dirichlet_dofs_;
    // Dof values for Dirichlet boundary conditions.
    std::vector<Scalar> dirichlet_values_;

    /// CSV writer
    CSVWriter<double> csv_writer_;
    /// Headers of columns in CSV file
    std::vector<std::string> csv_names_;
    /// Data  in CSV file (filled and written to file after
    /// every level)
    std::vector<double> csv_quantities_;

}; // end class PoissonAdaptive

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    // Set default parameter file
    std::string param_filename ( PARAM_FILENAME );
    std::string path_mesh;
    // If set take parameter file specified on console
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }
    // If set take mesh following path specified on console
    if ( argc > 2 )
    {
        path_mesh = std::string ( argv[2] );
    }
    try
    {
        // Create log files for INFO and DEBUG output
        //    std::ofstream info_log("poisson_tutorial_info_log");
        LogKeeper::get_log ( "info" ).set_target ( &( std::cout ) );
        //std::ofstream debug_log("poisson_tutorial_debug_log");
        LogKeeper::get_log ( "debug" ).set_target ( &( std::cout ) );
        //     std::ofstream error_log ( "poisson_tutorial_error_log" );
        LogKeeper::get_log ( "error" ).set_target ( &( std::cout ) );

        // Create application object and run it
        PoissonAdaptive app ( param_filename, path_mesh );

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

void PoissonAdaptive::build_initial_mesh ( )
{
    Timer timer;
    timer.start ( );

    // Read in the mesh.
    //The mesh is chosen according to the dimension of the problem.
    std::string mesh_name;

    if ( DIMENSION == 2 )
    {
        mesh_name = params_["Mesh"]["Filename2"].get<std::string>( "unit_square_inner_square.inp" );
#ifdef ELLIPSOID_BOUNDARY
        mesh_name = "unit_square_inner_square.inp";
#endif
    }
    if ( DIMENSION == 3 )
    {
        mesh_name = params_["Mesh"]["Filename3"].get<std::string>( "unit_cube.inp" );
    }

    std::string mesh_filename;
    if ( path_mesh.empty ( ) )
    {
        mesh_filename = std::string ( DATADIR ) + mesh_name;
    }
    else
    {
        mesh_filename = path_mesh + mesh_name;
    }

    std::vector<MasterSlave> period ( 0, MasterSlave ( 0., 0., 0., 0 ) );

    SharedVertexTable shared_verts;
    if ( rank_ == MASTER_RANK )
    {
#ifdef USE_MESH_P4EST
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period, mesh::IMPL_P4EST );
#else
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period, mesh::IMPL_DBVIEW );
#endif
        // Refine the mesh until the initial refinement level is reached.
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( 3 );
        if ( initial_ref_lvl > 0 )
        {
            LOG_DEBUG ( 1, "Initial refinement level " << initial_ref_lvl );
            master_mesh_ = master_mesh_->refine_uniform_seq ( initial_ref_lvl );
            //master_mesh_->do_tests(0);
        }
        refinement_level_ = initial_ref_lvl;
    }
    // Broadcast information from master to slaves.
    MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
    int uniform_ref_steps;

#ifdef USE_MESH_P4EST
    mesh_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &uniform_ref_steps, mesh::IMPL_P4EST );
    mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_P4EST, 2 );
    if ( rank_ == 0 )
    {
        //mesh_->do_tests(1);
    }
#else
    mesh_without_ghost_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &uniform_ref_steps, mesh::IMPL_DBVIEW );
    mesh_ = compute_ghost_cells ( *mesh_without_ghost_, comm_, shared_verts, mesh::IMPL_DBVIEW, 1 );
#endif

    refinement_level_ += uniform_ref_steps;
    // Add subdomain and remote index information needed by the library
    if ( num_partitions_ == 1 )
    {
        std::vector<int> remote_index
                (
                  mesh_->num_entities ( mesh_->tdim ( ) ),
                  -1
                  );

        AttributePtr remote_index_attr ( new IntAttribute ( remote_index ) );
        mesh_->add_attribute ( "_remote_index_",
                               DIMENSION,
                               remote_index_attr );

        std::vector<int> subdomain
                ( mesh_->num_entities ( mesh_->tdim ( ) ),
                  0
                  );
        AttributePtr subdomain_attr ( new IntAttribute ( subdomain ) );
        mesh_->add_attribute ( "_sub_domain_",
                               DIMENSION,
                               subdomain_attr );
    }

#ifdef ELLIPSOID_BOUNDARY
    // circle with radius = 1
    Coordinate radius = 1.;
    Ellipsoid circle ( radius, radius );
    adapt_boundary_to_function ( mesh_, circle );

    // Refine the mesh until the initial refinement level is reached.
    const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( 3 );

    for ( int r = 0; r < initial_ref_lvl; ++r )
    {
        adapt_boundary_to_function ( mesh_, circle );
        ++refinement_level_;
    }
#endif

    timer.stop ( );
    LOG_INFO ( "Build initial mesh time ", timer.get_duration ( ) );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );

    PVtkWriter writer ( comm_ );
    std::ostringstream name;
    name << "poisson_adaptive_initial_mesh.pvtu";
    std::string output_file = name.str ( );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

void PoissonAdaptive::prepare_system ( )
{
    Timer timer;
    timer.start ( );

#ifndef WITH_HYPRE
    if ( matrix_ != NULL )
    {
        delete matrix_;
    }
    if ( rhs_ != NULL )
    {
        delete rhs_;
    }
    if ( sol_ != NULL )
    {
        delete sol_;
    }
    if ( exact_sol_ != NULL )
    {
        delete exact_sol_;
    }

    matrix_ = new CoupledMatrix<Scalar>;
    rhs_ = new CoupledVector<Scalar>;
    sol_ = new CoupledVector<Scalar>;
    exact_sol_ = new CoupledVector<Scalar>;
#endif

    // Assign degrees to each element.
    const int fe_degree = params_["Mesh"]["FeDegree"].get<int>( 1 );
    std::vector< int > degrees ( 1, fe_degree );

    // Initialize the VectorSpace object.
    space_.Clear ( );
    space_.Init ( degrees, *mesh_ );

    this->csv_quantities_.push_back ( space_.dof ( ).ndofs_global ( ) );

    // Setup couplings object.
    couplings_.Clear ( );
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

#ifdef WITH_HYPRE
    matrix_.Init ( comm_, couplings_ );
    matrix_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                            vec2ptr ( sparsity.diagonal_cols ),
                            sparsity.diagonal_rows.size ( ),
                            vec2ptr ( sparsity.off_diagonal_rows ),
                            vec2ptr ( sparsity.off_diagonal_cols ),
                            sparsity.off_diagonal_rows.size ( ) );
    //matrix_hypre_.print_statistics();

    rhs_.Init ( comm_, couplings_ );
    sol_.Init ( comm_, couplings_ );
    exact_sol_.Init ( comm_, couplings_ );
    //rhs_hypre_.print_statistics();
    //sol_hypre_.print_statistics();
#else

    matrix_->Init ( comm_, couplings_, CPU, NAIVE, CSR );
    rhs_->Init ( comm_, couplings_, CPU, NAIVE );
    sol_->Init ( comm_, couplings_, CPU, NAIVE );
    exact_sol_->Init ( comm_, couplings_, CPU, NAIVE );

    // Initialize structure of LA objects.
    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );

    rhs_->InitStructure ( );
    sol_->InitStructure ( );
    exact_sol_->InitStructure ( );

    // Zero all linear algebra objects.
    matrix_->Zeros ( );
    rhs_->Zeros ( );
    sol_->Zeros ( );
    exact_sol_->Zeros ( );
#endif

    // Compute Dirichlet BC dofs and values using known exact solution.
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    DirichletZero zero;

    // The function compute_dirichlet_dofs_and_values des not yet work for 1D.
    if ( DIMENSION == 1 )
    {
        dirichlet_values_.resize ( 2, 0.0 );

        // Loop over all cells.
        for ( EntityIterator facet_it = mesh_->begin ( DIMENSION - 1 ),
              facet_end = mesh_->end ( DIMENSION - 1 );
              facet_it != facet_end;
              ++facet_it )
        {

            // Returns the number of neighbors for each cell to check
            // if it is on the facet.
            const EntityCount num_cell_neighbors = facet_it
                    ->num_incident_entities
                    ( DIMENSION );

            // If it lies on the facet, the corresponding DOF is a Dirichlet
            // DOF and is added to dirichlet_dofs_.
            if ( num_cell_neighbors == 1 )
            {
                std::vector<int> dof_number_;
                space_.dof ( ).get_dofs_on_subentity
                        (
                          0,
                          facet_it->begin_incident ( DIMENSION )->index ( ),
                          0,
                          facet_it->index ( ),
                          dof_number_
                          );
                dirichlet_dofs_.push_back ( dof_number_[0] );
            }
        }
    }
        //homogeneous Dirichlet boundaries
    else
    {
        compute_dirichlet_dofs_and_values ( zero, space_, 0, dirichlet_dofs_,
                                            dirichlet_values_ );
    }
    timer.stop ( );
    LOG_INFO ( "Space setup time ", timer.get_duration ( ) );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );
}

void PoissonAdaptive::higher_order_interpolation ( )
{
    Timer timer;
    timer.start ( );

#ifdef USE_MESH_P4EST
    LOG_INFO ( "Higher order interpolation", " on level " << refinement_level_ );

    assert ( this->mesh_->is_uniformly_coarsenable ( ) );
#    ifdef WITH_HYPRE
    SpacePatchInterpolation<LADescriptorHypreD, mesh::IMPL_P4EST> patch_interpolation;
    HypreVector<double> inter_sol;
#    else
    SpacePatchInterpolation<LADescriptorCoupledD, mesh::IMPL_P4EST> patch_interpolation;
    CoupledVector<Scalar> inter_sol;
#    endif

    // init interpolating space and interpolation map
    patch_interpolation.init ( &this->space_ );

    // get interpolating space
    const VectorSpace<double>* inter_space = patch_interpolation.get_space ( );
    const Couplings<double>& inter_couplings = patch_interpolation.get_couplings ( );

    // interpolate solution
#    ifdef WITH_HYPRE
    inter_sol.Init ( comm_, inter_couplings );
    inter_sol.Zeros ( );
    sol_.Update ( );
    patch_interpolation.interpolate ( sol_, inter_sol );
#    else
    inter_sol.Init ( comm_, inter_couplings, CPU, NAIVE );
    inter_sol.InitStructure ( );
    inter_sol.Zeros ( );
    sol_->UpdateCouplings ( );
    patch_interpolation.interpolate ( *sol_, inter_sol );
#    endif

    timer.stop ( );
    LOG_INFO ( "Patch interpolation time ", timer.get_duration ( ) );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );

    // visualize
    // Setup visualization object.
    LOG_INFO ( "Visualize", " interpolation mesh " << refinement_level_ );
    const int num_intervals = params_["Mesh"]["FeDegree"].get<int>( 1 ) /** 2*/;
    ParallelCellVisualization<double> visu
            ( *inter_space,
              num_intervals,
              comm_,
              MASTER_RANK );

    // Generate filename.
    std::stringstream name;
    name << "high_solution" << refinement_level_;

#    ifdef WITH_HYPRE
    visu.visualize ( EvalFeFunction<LAD>( *inter_space, inter_sol ), "high_u" );
#    else
    visu.visualize ( EvalFeFunction<LAD>( *inter_space, inter_sol ), "high_u" );
#    endif

    visu.write ( name.str ( ) );
#endif
}

void PoissonAdaptive::assemble_system ( )
{
    Timer assembly_hf_la;
    assembly_hf_la.start ( );

    // Assemble matrix and right-hand-side vector.
    LocalPoissonAssembler<ExactSol> local_asm;
#ifdef WITH_HYPRE
    global_asm_.assemble_matrix ( space_, local_asm, matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, rhs_ );
#else
    global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, *rhs_ );
#endif

    assembly_hf_la.stop ( );
    LOG_INFO ( "Assembly time with HiFlow LA", assembly_hf_la.get_duration ( ) );
    this->csv_quantities_.push_back ( assembly_hf_la.get_duration ( ) );

#ifdef WITH_HYPRE
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        matrix_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ),
                                   dirichlet_dofs_.size ( ), 1.0 );
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                         vec2ptr ( dirichlet_values_ ) );
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                         vec2ptr ( dirichlet_values_ ) );
    }
    sol_.Update ( );
    rhs_.Update ( );
#else
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ),
                                    dirichlet_dofs_.size ( ), 1.0 );
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
    sol_->UpdateCouplings ( );
    rhs_->UpdateCouplings ( );
#endif

    ExactSol exact_sol_fn;
    for ( mesh::EntityIterator it = mesh_->begin ( DIMENSION ), end_it = mesh_->end ( DIMENSION ); it != end_it; ++it )
    {
        std::vector<int> global_dof_ids;
        space_.GetDofIndices ( 0, *it, &global_dof_ids );
        int num_dofs = global_dof_ids.size ( );
        std::vector<double> values;
        values.resize ( num_dofs, 0. );

        std::vector< Coord > coords;
        space_.dof ( ).get_coord_on_cell ( 0, it->index ( ), coords );
        for ( int i = 0; i < num_dofs; i++ )
        {
            if ( rank_ == space_.dof ( ).owner_of_dof ( global_dof_ids.at ( i ) ) )
            {
                Vec<DIMENSION, double> pt;
                for ( int l = 0; l < DIMENSION; ++l )
                {
                    pt[l] = coords[i][l];
                }
                double val = exact_sol_fn ( pt );
#ifdef WITH_HYPRE
                exact_sol_.SetValues ( &global_dof_ids.at ( i ), 1, &val );
#else
                exact_sol_->SetValues ( &global_dof_ids.at ( i ), 1, &val );
#endif
            }

            if ( DEBUG_LEVEL >= 3 )
            {
                if ( rank_ == space_.dof ( ).owner_of_dof ( global_dof_ids.at ( i ) ) )
                {
                    double val;
#ifdef WITH_HYPRE
                    matrix_.GetValues ( &global_dof_ids.at ( i ), 1, &global_dof_ids.at ( i ), 1, &val );
#else
                    matrix_->GetValues ( &global_dof_ids.at ( i ), 1, &global_dof_ids.at ( i ), 1, &val );
#endif

                    if ( val == 1. )
                    {
                        LOG_DEBUG ( 3, "[" << rank_ << "] Dof id " << global_dof_ids.at ( i ) << " has diagonal entry " << val << " and coord " << coords[i][0] << ", " << coords[i][1] );
                    }
                }
            }
        }
    }

#ifdef WITH_HYPRE
    exact_sol_.Update ( );
    interpolate_constrained_vector ( space_, exact_sol_ );
    exact_sol_.Update ( );
#else
    exact_sol_->UpdateCouplings ( );
    interpolate_constrained_vector ( space_, *exact_sol_ );
    exact_sol_->UpdateCouplings ( );
#endif
}

void PoissonAdaptive::solve_system ( )
{
    Timer timer;
    timer.start ( );

    const int max_iter = params_["LinearSolver"]["MaxIterations"].get<int>( 1000 );
    const double abs_tol = params_["LinearSolver"]["AbsTolerance"].get<double>( 1e-12 );
    const double rel_tol = params_["LinearSolver"]["RelTolerance"].get<double>( 1e-6 );
    const double div_tol = params_["LinearSolver"]["DivTolerance"].get<double>( 1e6 );
    const int basis_size = params_["LinearSolver"]["SizeBasis"].get<int>( 100 );
    const std::string solver_type = params_["LinearSolver"]["Name"].get<std::string>( "GMRES" );
    const std::string method = params_["LinearSolver"]["Method"].get<std::string>( "NoPreconditioning" );

#ifdef WITH_HYPRE
    HypreCG<LAD> solver_hypre;

    solver_hypre.SetupOperator ( matrix_ );
    solver_hypre.InitControl ( max_iter, abs_tol, rel_tol );
    solver_hypre.SetPrintLevel ( 0 );

    HypreBoomerAMG<LAD> precond_boomer;
    precond_boomer.InitControl ( 1, 0., 0. );
    precond_boomer.SetCycleType ( 2 );
    precond_boomer.SetRelaxType ( 6 );
    //precond_boomer.SetCycleRelaxType ( 6, 3 );
    precond_boomer.SetCycleNumSweeps ( 3, 1 );
    precond_boomer.SetCycleNumSweeps ( 3, 2 );
    precond_boomer.SetRelaxWt ( 0.5 );
    precond_boomer.SetCoarsenType ( 10 );
    precond_boomer.SetInterpType ( 6 );
#    if DIMENSION == 2
    precond_boomer.SetStrongThreshold ( 0.25 );
#    else
    precond_boomer.SetStrongThreshold ( 0.6 );
#    endif
    precond_boomer.SetAggNumLevels ( 25 );
    precond_boomer.SetSmoothType ( 0 );
    precond_boomer.SetupOperator ( matrix_ );

    precond_boomer.SetPrintLevel ( 0 );

    precond_boomer.Build ( );

    solver_hypre.SetupPreconditioner ( precond_boomer );

    solver_hypre.Solve ( rhs_, &sol_ );
#else
    LinearSolver<LAD>* solver;
    GMRES<LAD> solver_gmres;
    CG<LAD> solver_cg;

    if ( solver_type == "GMRES" )
    {
        solver = &solver_gmres;
        solver_gmres.InitParameter ( basis_size, method );
    }
    else if ( solver_type == "CG" )
    {
        solver = &solver_cg;
    }
    else
    {
        std::cout << " Solver Type not supported! Exit now " << std::endl;
        exit ( -1 );
    }

    solver->InitControl ( max_iter, abs_tol, rel_tol, div_tol );
    solver->SetPrintLevel ( 0 );

    PreconditionerBlockJacobiStand<LAD> precond;
    precond.InitParameter ( );
    precond.Init_SSOR ( 1.3 );

    precond.SetupOperator ( *matrix_ );
    if ( method != "NoPreconditioning" )
    {
        solver->SetupPreconditioner ( precond );
    }

    solver->SetupOperator ( *matrix_ );
    solver->Solve ( *rhs_, sol_ );
#endif

#ifdef WITH_HYPRE
    sol_.Update ( );
    interpolate_constrained_vector ( space_, sol_ );
    sol_.Update ( );
#else
    sol_->UpdateCouplings ( );
    interpolate_constrained_vector ( space_, *sol_ );
    sol_->UpdateCouplings ( );
#endif

    timer.stop ( );
#ifdef WITH_HYPRE
    LOG_INFO ( "Linear solver", " iter " << solver_hypre.iter ( ) << " CG res " << solver_hypre.res ( ) );
#else
    LOG_INFO ( "Linear solver", " iter " << solver->iter ( ) << " CG res " << solver->res ( ) );
#endif
    LOG_INFO ( "Solving time ", timer.get_duration ( ) );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );

}

void PoissonAdaptive::compute_error ( )
{
    Timer timer;
    timer.start ( );

    // Compute square of the L2 error on each element, putting the
    // values into L2_err_.
    L2_err_.clear ( );
#ifdef WITH_HYPRE
    L2ErrorIntegrator<ExactSol> L2_int ( sol_ );
#else
    L2ErrorIntegrator<ExactSol> L2_int ( *( sol_ ) );
#endif
    global_asm_.assemble_scalar ( space_, L2_int, L2_err_ );

    // Create attribute with L2 error for output.
    AttributePtr L2_err_attr ( new DoubleAttribute ( L2_err_ ) );

    //LOG_INFO("error", "[" << rank_ << " ] L2 error per cell: " << string_from_range(L2_err_.begin(), L2_err_.end()));

    mesh_->add_attribute ( "L2 error", DIMENSION, L2_err_attr );
    double total_L2_err = std::accumulate
            (
              L2_err_.begin ( ),
              L2_err_.end ( ),
              0.
              );
    double global_L2_err = 0.;
    MPI_Allreduce ( &total_L2_err,
                    &global_L2_err,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_ );
    LOG_INFO ( "error", "Local L2 error on partition " << rank_ << " = " << std::sqrt ( total_L2_err ) );

    // Compute square of the H1 error on each element, putting the
    // values into H1_err_.
    H1_err_.clear ( );
#ifdef WITH_HYPRE
    H1ErrorIntegrator<ExactSol> H1_int ( sol_ );
#else
    H1ErrorIntegrator<ExactSol> H1_int ( *( sol_ ) );
#endif
    global_asm_.assemble_scalar ( space_, H1_int, H1_err_ );

    // Create attribute with H1 error for output.
    //LOG_DEBUG(2,"[" << rank_ << " ] L2 error per cell: " << string_from_range(H1_err_.begin(), H1_err_.end()));

    AttributePtr H1_err_attr ( new DoubleAttribute ( H1_err_ ) );
    mesh_->add_attribute ( "H1 error", DIMENSION, H1_err_attr );
    double total_H1_err = std::accumulate
            (
              H1_err_.begin ( ),
              H1_err_.end ( ),
              0.
              );
    double global_H1_err = 0.;
    MPI_Allreduce ( &total_H1_err,
                    &global_H1_err,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_ );
    timer.stop ( );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );
    this->csv_quantities_.push_back ( std::sqrt ( global_L2_err ) );
    this->csv_quantities_.push_back ( std::sqrt ( global_H1_err ) );

    LOG_INFO ( "error", "Local H1 error on partition " << rank_ << " = " << std::sqrt ( total_H1_err ) );

    // Output on consol: number of elements and global H^1 error
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "error", "H1size: " << L2_err_.size ( ) << ", H1 error: " << std::sqrt ( global_H1_err ) );
    }
}

void PoissonAdaptive::error_estimator ( )
{
    Timer timer;
    timer.start ( );

    // ***************************************************************
    // get number of interior facets and boundary facets
    InterfaceList if_list = InterfaceList::create ( mesh_ );

    int Number_Boundary = 0;
    for ( InterfaceList::const_iterator it = if_list.begin ( ),
          end_it = if_list.end ( );
          it != end_it;
          ++it )
    {

        int remote_index_master = -10;
        mesh_->get_attribute_value
                ( "_remote_index_",
                  mesh_->tdim ( ),
                  it->master_index ( ),
                  &remote_index_master
                  );

        const int num_slaves = it->num_slaves ( );
        if ( remote_index_master == -1 )
        {
            if ( num_slaves == 0 )
            {
                Number_Boundary += 1;
            }
        }
    }

    // ***************************************************************
    // Create attribute with H1 error for output.
    AttributePtr H1_err_attr ( new DoubleAttribute ( H1_err_ ) );
    mesh_->add_attribute ( "H1 error", DIMENSION, H1_err_attr );
    double total_H1_err = std::accumulate
            (
              H1_err_.begin ( ),
              H1_err_.end ( ),
              0.
              );
    double global_H1_err = 0.;
    MPI_Allreduce ( &total_H1_err,
                    &global_H1_err,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_ );

    // ***************************************************************
    // rho_T : cellwise residual
    rho_cell_.clear ( );
#ifdef WITH_HYPRE
    CellTermAssembler rho_cell_int ( sol_ );
#else
    CellTermAssembler rho_cell_int ( *sol_ );
#endif

    global_asm_.assemble_scalar ( space_, rho_cell_int, rho_cell_ );

    // Create attribute with term on the cells for output.
    AttributePtr rho_cell_attr ( new DoubleAttribute ( rho_cell_ ) );
    mesh_->add_attribute ( "rho_T", DIMENSION, rho_cell_attr );
    double total_rho_cell = std::accumulate
            (
              rho_cell_.begin ( ),
              rho_cell_.end ( ),
              0.
              );
    double global_rho_cell = 0.;
    MPI_Allreduce ( &total_rho_cell,
                    &global_rho_cell,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_
                    );
    int local_size_T = rho_cell_.size ( );
    int global_size_T = 0;
    MPI_Allreduce ( &local_size_T,
                    &global_size_T,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    comm_
                    );

    LOG_INFO ( "rho_T", "partition " << rank_ << ": size = " << local_size_T << ", local value = " << std::sqrt ( total_rho_cell ) );
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "rho_T", "global size = " << global_size_T << ", global value = " << std::sqrt ( global_rho_cell ) );
    }

    // ***************************************************************
    // rho_S
    rho_bcell_.clear ( );
#ifdef WITH_HYPRE
    BoundCellTermAssembler rho_bcell_int ( sol_ );
#else
    BoundCellTermAssembler rho_bcell_int ( *sol_ );
#endif

    jump_term_asm_.assemble_interface_scalar_cells
            (
              space_,
              rho_bcell_int,
              rho_bcell_
              );

    // Create attribute with boundary cell term for output.
    AttributePtr rho_bcell_attr ( new DoubleAttribute ( rho_bcell_ ) );
    mesh_->add_attribute ( "rho_S", DIMENSION, rho_bcell_attr );
    double total_rho_bcell = std::accumulate
            (
              rho_bcell_.begin ( ),
              rho_bcell_.end ( ),
              0. );
    double global_rho_bcell = 0.;
    MPI_Allreduce ( &total_rho_bcell,
                    &global_rho_bcell,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_
                    );
    int local_size_S = rho_bcell_.size ( );
    int global_size_S = 0;
    MPI_Allreduce ( &local_size_S,
                    &global_size_S,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    comm_
                    );

    LOG_INFO ( "rho_S", "partition " << rank_ << ": size = " << local_size_S << ", local value = " << std::sqrt ( total_rho_bcell ) );
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "rho_S", "global size = " << global_size_S << ", global value = " << std::sqrt ( global_rho_bcell ) );
    }

    // ***************************************************************
    // rho_E : jump terms
    rho_jump_.clear ( );
#ifdef WITH_HYPRE
    JumpTermAssembler rho_jump_int ( sol_ );
#else
    JumpTermAssembler rho_jump_int ( *sol_ );
#endif
    jump_term_asm_.assemble_interface_scalar_cells
            (
              space_,
              rho_jump_int,
              rho_jump_
              );
    // Create attribute with inner edge jump term for output.
    AttributePtr rho_jump_attr ( new DoubleAttribute ( rho_jump_ ) );
    mesh_->add_attribute ( "rho_E", DIMENSION, rho_jump_attr );
    double total_rho_jump = std::accumulate
            (
              rho_jump_.begin ( ),
              rho_jump_.end ( ),
              0.
              );
    double global_rho_jump = 0.;
    MPI_Allreduce ( &total_rho_jump,
                    &global_rho_jump,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_
                    );
    int local_size_E = if_list.size ( ) - Number_Boundary;
    int global_size_E = 0;
    MPI_Allreduce ( &local_size_E,
                    &global_size_E,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    comm_
                    );

    LOG_INFO ( "rho_E", "partition " << rank_ << ": size = " << local_size_E << ", local value = " << std::sqrt ( total_rho_jump ) );
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "rho_E", "global size = " << global_size_E << ", global value = " << std::sqrt ( global_rho_jump ) );
    }

    // ***************************************************************
    // rho_A
    rho_boundary_.clear ( );
#ifdef WITH_HYPRE
    BoundaryTermAssembler rho_boundary_int ( sol_ );
#else
    BoundaryTermAssembler rho_boundary_int ( *sol_ );
#endif

    jump_term_asm_.assemble_interface_scalar_cells
            (
              space_,
              rho_boundary_int,
              rho_boundary_
              );

    // Create attribute with inner edge jump term for output.
    AttributePtr rho_boundary_attr ( new DoubleAttribute ( rho_boundary_ ) );
    mesh_->add_attribute ( "rho_A", DIMENSION, rho_boundary_attr );
    double total_rho_boundary = std::accumulate
            (
              rho_boundary_.begin ( ),
              rho_boundary_.end ( ),
              0.
              );
    double global_rho_boundary = 0.;
    MPI_Allreduce ( &total_rho_boundary,
                    &global_rho_boundary,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_
                    );
    int local_size_A = Number_Boundary;
    int global_size_A = 0;
    MPI_Allreduce ( &local_size_A,
                    &global_size_A,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    comm_
                    );

    LOG_INFO ( "rho_A", "partition " << rank_ << ": size = " << local_size_A << ", local value = " << std::sqrt ( total_rho_boundary ) );
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "rho_A", "global size = " << global_size_A << ", global value = " << std::sqrt ( global_rho_boundary ) );
    }

    // ***************************************************************
    // Standard estimator : rho_T + rho_E
    std_estimator_.clear ( );
    for ( int i = 0; i < L2_err_.size ( ); ++i )
    {
        std_estimator_.push_back ( rho_cell_[i] + rho_jump_[i] );
    }
    // Create attribute with term on the cells for output.
    AttributePtr std_estimator_attr ( new DoubleAttribute ( std_estimator_ ) );
    mesh_->add_attribute
            ( "std_estimator",
              DIMENSION,
              std_estimator_attr );
    double total_std_estimator = std::accumulate
            (
              std_estimator_.begin ( ),
              std_estimator_.end ( ),
              0.
              );
    this->global_std_estimator_ = 0.;
    MPI_Allreduce ( &total_std_estimator,
                    &global_std_estimator_,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_
                    );
    int local_size_std = std_estimator_.size ( );
    global_size_std_ = 0;
    MPI_Allreduce ( &local_size_std,
                    &global_size_std_,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    comm_
                    );

    LOG_INFO ( "std_estimator", "partition " << rank_ << ": size = " << local_size_std << ", local value = " << std::sqrt ( total_std_estimator ) );
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "std_estimator", "global size = " << global_size_std_ << ", global value = " << std::sqrt ( global_std_estimator_ ) );
    }

    // ***************************************************************
    // New estimator : rho_T + rho_E + rho_S + rho_A
    new_estimator_.clear ( );
    for ( int i = 0; i < L2_err_.size ( ); ++i )
    {
        new_estimator_.push_back
                (
                  rho_cell_[i]
                  + rho_jump_[i]
                  + rho_bcell_[i]
                  + rho_boundary_[i]
                  );
    }
    // Create attribute with term on the cells for output.
    AttributePtr new_estimator_attr ( new DoubleAttribute ( new_estimator_ ) );
    mesh_->add_attribute
            ( "new_estimator",
              DIMENSION,
              new_estimator_attr );
    double total_new_estimator = std::accumulate
            (
              new_estimator_.begin ( ),
              new_estimator_.end ( ),
              0.
              );
    this->global_new_estimator_ = 0.;
    MPI_Allreduce ( &total_new_estimator,
                    &global_new_estimator_,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    comm_ );
    int local_size_new = new_estimator_.size ( );
    global_size_new_ = 0;
    MPI_Allreduce ( &local_size_new,
                    &global_size_new_,
                    1,
                    MPI_INT,
                    MPI_SUM,
                    comm_
                    );

    LOG_INFO ( "new_estimator", "partition " << rank_ << ": size = " << local_size_new << ", local value = " << std::sqrt ( total_new_estimator ) );
    if ( rank_ == MASTER_RANK )
    {
        LOG_INFO ( "new_estimator", "global size = " << global_size_new_ << ", global value = " << std::sqrt ( global_new_estimator_ ) );
    }

    timer.stop ( );
    this->csv_quantities_.push_back ( timer.get_duration ( ) );
    this->csv_quantities_.push_back ( std::sqrt ( global_std_estimator_ ) );

    // ***************************************************************
    // Division new error estimator or std estimator / real H^1 error

    int estimator_ = params_["Mesh"]["Estimator"].get<int>( 3 );
    if ( estimator_ == 1 )
    {
        double relation = std::abs ( sqrt ( this->global_std_estimator_ )
                                     / sqrt ( global_H1_err ) );

        if ( rank_ == MASTER_RANK )
        {
            LOG_INFO ( "std_estimator", " relation: " << relation );
        }
    }
    else
    {
        double relation = std::abs ( sqrt ( this->global_new_estimator_ )
                                     / sqrt ( global_H1_err ) );

        if ( rank_ == MASTER_RANK )
        {
            LOG_INFO ( "new_estimator", "relation: " << relation );
        }
    }
}

void PoissonAdaptive::visualize ( )
{
    // Setup visualization object.
    LOG_INFO ( "Visualize", " solution on level " << refinement_level_ );
    const int num_intervals = 1; //params_["Mesh"]["FeDegree"].get<int>( 1 );
    ParallelCellVisualization<double> visu
            ( space_,
              num_intervals,
              comm_,
              MASTER_RANK );

    // Generate filename.
    std::stringstream name;
    name << "solution" << refinement_level_;

    std::vector<double> remote_index
            ( mesh_->num_entities ( mesh_->tdim ( ) ),
              0. );
    std::vector<double> sub_domain
            ( mesh_->num_entities ( mesh_->tdim ( ) ),
              0. );
    std::vector<double> material_number
            ( mesh_->num_entities ( mesh_->tdim ( ) ),
              0. );

    for ( mesh::EntityIterator it = mesh_->begin
          ( mesh_->tdim ( ) );
          it != mesh_->end ( mesh_->tdim ( ) );
          ++it )
    {
        int temp1, temp2, temp3;
        mesh_->get_attribute_value ( "_remote_index_",
                                     mesh_->tdim ( ),
                                     it->index ( ),
                                     &temp1 );
        mesh_->get_attribute_value ( "_sub_domain_",
                                     mesh_->tdim ( ),
                                     it->index ( ),
                                     &temp2 );

        remote_index.at ( it->index ( ) ) = temp1;
        sub_domain.at ( it->index ( ) ) = temp2;

        material_number.at ( it->index ( ) ) = mesh_->get_material_number
                ( mesh_->tdim ( ),
                  it->index ( ) );
    }

#ifdef WITH_HYPRE
    sol_.Update ( );
    exact_sol_.Update ( );
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_ ), "u" );
    visu.visualize ( EvalFeFunction<LAD>( space_, exact_sol_ ), "exact" );

    visu.visualize ( EvalDerivativeFeFunction<LAD, DIMENSION>( space_, sol_, 0, 0 ), "du_dx" );
    visu.visualize ( EvalDerivativeFeFunction<LAD, DIMENSION>( space_, sol_, 0, 1 ), "du_dy" );
    if ( DIMENSION > 2 )
        visu.visualize ( EvalDerivativeFeFunction<LAD, DIMENSION>( space_, sol_, 0, 2 ), "du_dz" );
#else
    sol_->UpdateCouplings ( );
    exact_sol_->UpdateCouplings ( );
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ) ), "u" );
    visu.visualize ( EvalFeFunction<LAD>( space_, *( exact_sol_ ) ), "exact" );

    visu.visualize ( EvalDerivativeFeFunction<LAD, DIMENSION>( space_, *( sol_ ), 0, 0 ), "du_dx" );
    visu.visualize ( EvalDerivativeFeFunction<LAD, DIMENSION>( space_, *( sol_ ), 0, 1 ), "du_dy" );

    if ( DIMENSION > 2 )
        visu.visualize ( EvalDerivativeFeFunction<LAD, DIMENSION>( space_, *( sol_ ), 0, 2 ), "du_dz" );
#endif

    // visualize error measures
    visu.visualize_cell_data ( L2_err_, "L2" );
    visu.visualize_cell_data ( H1_err_, "H1" );
    visu.visualize_cell_data ( rho_cell_, "rho_T" );
    visu.visualize_cell_data ( rho_bcell_, "rho_S" );
    visu.visualize_cell_data ( rho_jump_, "rho_E" );
    visu.visualize_cell_data ( rho_boundary_, "rho_A" );
    visu.visualize_cell_data ( std_estimator_, "std_estimator" );
    visu.visualize_cell_data ( new_estimator_, "new_estimator" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.write ( name.str ( ) );
}

void PoissonAdaptive::adapt_uniform ( )
{
    // Refine mesh on master process. 1D parallel execution is not
    // yet implemented.
    if ( DIMENSION == 1 )
    {
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"]
                    .get<int>( 6 );
            if ( refinement_level_ >= final_ref_level )
            {
                is_done_ = true;
            }
            else
            {
                master_mesh_ = master_mesh_->refine ( );
                ++refinement_level_;
            }
        }
    }
    else
    {
        const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
        if ( refinement_level_ >= final_ref_level )
        {
            is_done_ = true;
        }
        else
        {
            SharedVertexTable shared_verts;
            Timer timer;
            timer.start ( );
#ifdef USE_MESH_P4EST
            mesh_ = mesh_->refine ( );

#else
            mesh_without_ghost_ = mesh_without_ghost_->refine ( );
#endif
            timer.stop ( );
            this->csv_quantities_.push_back ( timer.get_duration ( ) );

            timer.reset ( );
            timer.start ( );
#ifdef USE_MESH_P4EST
            mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_P4EST, 2 );
#else
            mesh_ = compute_ghost_cells ( *mesh_without_ghost_, comm_, shared_verts, mesh::IMPL_DBVIEW, 1 );
#endif
            timer.stop ( );
            this->csv_quantities_.push_back ( timer.get_duration ( ) );

#ifdef ELLIPSOID_BOUNDARY
            Coordinate radius = 1.;
            Ellipsoid circle ( radius, radius );
            adapt_boundary_to_function ( mesh_, circle );
#endif
            ++refinement_level_;
        }
        if ( mesh_->num_entities ( mesh_->tdim ( ) ) > 100000 )
        {
            is_done_ = true;
        }

        PVtkWriter writer ( comm_ );
        std::ostringstream name;
        name << "poisson_adaptive_mesh_" << refinement_level_ << ".pvtu";
        std::string output_file = name.str ( );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );
    }
}

void PoissonAdaptive::adapt ( )
{
    // Refine mesh on master process. 1D parallel execution is not
    // yet implemented.
    if ( DIMENSION == 1 )
    {
        const int final_ref_level = params_["Mesh"]["FinalRefLevel"]
                .get<int>( 6 );
        if ( refinement_level_ >= final_ref_level )
        {
            is_done_ = true;
        }
        else
        {
            ++refinement_level_;
        }
    }
    else
    {
        const int final_ref_level = params_["Mesh"]["FinalRefLevel"]
                .get<int>( 6 );
        if ( refinement_level_ >= final_ref_level )
        {
            is_done_ = true;
        }
        else
        {
            is_done_ = false;
            int estimator_ = params_["Mesh"]["Estimator"].get<int>( 3 );
            int coarsen_threshold = params_["Mesh"]["CoarsenThreshold"].get<int>( 10000 );
            int coarsen_flag = params_["Mesh"]["CoarsenFlag"].get<int>( 1 );
            double refine_frac = params_["Mesh"]["FractionToRefine"].get<double>( 0.2 );
            double coarsen_frac = params_["Mesh"]["FractionToCoarsen"].get<double>( 0.1 );
            int strategy = params_["Mesh"]["Strategy"].get<int>( 1 );
            int connection_mode = params_["Mesh"]["BalanceConnectionMode"].get<int>( 0 );
            int patch_mode_level = params_["Mesh"]["PatchModeLevel"].get<int>( -1 );

            int num_local_cells = 0;

            for ( mesh::EntityIterator it = this->mesh_->begin ( this->mesh_->tdim ( ) ),
                  e_it = this->mesh_->end ( this->mesh_->tdim ( ) );
                  it != e_it;
                  ++it )
            {

                int subdomain_index = -10;
                this->mesh_->get_attribute_value ( "_sub_domain_", this->mesh_->tdim ( ), it->index ( ), &subdomain_index );
                if ( subdomain_index == this->rank_ )
                {
                    ++num_local_cells;
                }
            }

            LOG_INFO ( "Number of cells on process " << this->rank_, num_local_cells );
            LOG_INFO ( "Number of cells on process " << this->rank_, mesh_->num_local_cells ( ) );
            LOG_INFO ( "Number of ghost cells on process " << this->rank_, mesh_->num_ghost_cells ( ) );

            int num_global_cells = 0;
            MPI_Allreduce ( &num_local_cells,
                            &num_global_cells,
                            1,
                            MPI_INT,
                            MPI_SUM,
                            comm_
                            );
            if ( this->rank_ == MASTER_RANK )
            {
                LOG_INFO ( "Global number of cells", num_global_cells );
            }

            std::vector<double> indicators;
            if ( estimator_ == 1 )
            {
                indicators = std_estimator_;
            }
            if ( estimator_ == 2 )
            {
                indicators = new_estimator_;
            }

            std::vector<int> refinemnts ( indicators.size ( ), 0 );
            std::vector<int> sort_err ( indicators.size ( ), 0 );
            for ( int i = 0; i < std_estimator_.size ( ); ++i )
            {
                sort_err[i] = i;
            }
            sortingPermutation ( std_estimator_, sort_err );

            // Refinement of all cells over mean value of the std estimator
            if ( strategy == 1 )
            {
                int cell_index = sort_err.size ( ) - 1;
                while (
                        (
                        std_estimator_[sort_err[cell_index]]
                        >
                        ( this->global_std_estimator_ / num_global_cells )
                        )
                        && ( cell_index >= 0 )
                        )

                {
                    refinemnts[sort_err[cell_index]] = 1;
                    is_done_ = false;
                    cell_index--;
                }
            }
            if ( strategy == 2 )
            {
                // 1. Mark cells for refinement
                int first_cell = std::floor ( ( 1. - refine_frac ) * sort_err.size ( ) );

                for ( int l = first_cell; l < sort_err.size ( ); ++l )
                {
                    EntityNumber cell_index = sort_err[l];
                    refinemnts[cell_index] = 1;
                }

                // 2.Mark cells for coarsening
                if ( num_global_cells >= coarsen_threshold )
                {
                    int last_cell = std::ceil ( coarsen_frac * sort_err.size ( ) );
                    for ( int l = 0; l < last_cell; ++l )
                    {
                        EntityNumber cell_index = sort_err[l];
                        refinemnts[cell_index] = -coarsen_flag;
                    }
                }
                is_done_ = false;
                LOG_DEBUG ( 3, " ref : " << string_from_range ( refinemnts.begin ( ), refinemnts.end ( ) ) );
            }

            SharedVertexTable shared_verts;
            Timer timer;
            timer.start ( );
#ifdef USE_MESH_P4EST
            boost::intrusive_ptr<MeshPXest> mesh_pXest = boost::static_pointer_cast<MeshPXest> ( mesh_ );
            mesh_pXest->set_connection_mode ( connection_mode );

            if ( this->refinement_level_ >= patch_mode_level )
            {
                if ( patch_mode_level >= 0 )
                {
                    mesh_pXest->set_patch_mode ( true );
                    mesh_pXest->set_connection_mode ( 2 );
                    LOG_INFO ( "Adaption: ", " enable patch mode, curretn level: " << refinement_level_ );
                }
            }
#endif
            mesh_ = mesh_->refine ( refinemnts );
            timer.stop ( );
            this->csv_quantities_.push_back ( timer.get_duration ( ) );

            timer.reset ( );
            timer.start ( );
#ifdef USE_MESH_P4EST
            mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_P4EST, 2 );
#else
            mesh_ = compute_ghost_cells ( *mesh_, comm_, shared_verts, mesh::IMPL_DBVIEW, 1 );
#endif
            timer.stop ( );
            this->csv_quantities_.push_back ( timer.get_duration ( ) );

            LOG_INFO ( "Cells: ", "Actual new number " << mesh_->num_global_cells ( comm_ ) );
            LOG_INFO ( "Number of cells on process " << this->rank_, mesh_->num_local_cells ( ) );
            LOG_INFO ( "Number of ghost cells on process " << this->rank_, mesh_->num_ghost_cells ( ) );

#ifdef ELLIPSOID_BOUNDARY
            // Unit square -> Ellipsoid implementation a = b = radius = 1.
            Coordinate radius = 1.;
            Ellipsoid circle ( radius, radius );
            adapt_boundary_to_function ( mesh_, circle );
#endif

            // Regularize mesh
            const TDim td = DIMENSION;

#ifdef USE_MESH_P4EST
            bool mesh_is_regular = true;
#else
            bool mesh_is_regular = false;
#endif
            while ( !mesh_is_regular )
            {
                mesh_is_regular = true;

                // Vector to hold the cells that need to be refined
                const EntityCount num_cells = mesh_->num_entities ( td );
                std::vector<int> refinements;
                refinements.resize ( num_cells, -5 );

                // Create interface list
                InterfaceList if_list = InterfaceList::create ( mesh_ );

                // Loop over interfaces
                for ( InterfaceList::const_iterator it = if_list.begin ( ),
                      end_it = if_list.end ( );
                      it != end_it; ++it )
                {

                    // Do hanging nodes exist?
                    // Not more than 1 hanging node on each edge allowed.
                    const int num_slaves = it->num_slaves ( );
                    if ( num_slaves > 2 )
                    {
                        int default_num_children = 1;
                        if ( td == 3 )
                        {
                            Entity const& master_cell = mesh_
                                    ->get_entity ( mesh_->tdim ( ),
                                                   it->master_index ( ) );
                            if ( master_cell.cell_type ( ).tag ( )
                                 == CellType::HEXAHEDRON )
                                default_num_children = 4;
                            else if ( master_cell.cell_type ( ).tag ( )
                                      == CellType::TETRAHEDRON )
                                default_num_children = 3;
                            else
                            {
                                std::cout
                                        << "This cell type is not treated ..."
                                        << " implement me."
                                        << std::endl;
                                interminable_assert ( 0 );
                            }
                        }

                        // Check regularity of this interface
                        if ( num_slaves > default_num_children )
                        {
                            mesh_is_regular = false;

                            // Add master cell to refinement list
                            refinements[it->master_index ( )] = 0;
                        }
                    }
                }

                // Refine if needed.
                // If mesh is not regular -> build new mesh
                if ( !mesh_is_regular )
                {
                    mesh_ = mesh_->refine ( refinements );

#ifdef ELLIPSOID_BOUNDARY
                    // Unit square -> Ellipsoid implementation
                    // a = b = radius = 1.
                    Coordinate radius = 1.;
                    Ellipsoid circle ( radius, radius );
                    adapt_boundary_to_function ( mesh_, circle );
#endif
                }

            }
            ++refinement_level_;
        }

        if ( !is_done_ )
        {
            if ( num_partitions_ == 1 )
            {
                // Add subdomain and remote index information needed by the library.
                std::vector<int> remote_index
                        (
                          mesh_->num_entities ( mesh_->tdim ( ) ),
                          -1
                          );
                AttributePtr remote_index_attr ( new IntAttribute ( remote_index ) );
                mesh_->add_attribute
                        ( "_remote_index_",
                          DIMENSION,
                          remote_index_attr );

                std::vector<int> subdomain
                        (
                          mesh_->num_entities ( mesh_->tdim ( ) ),
                          0
                          );
                AttributePtr subdomain_attr ( new IntAttribute ( subdomain ) );
                mesh_->add_attribute
                        ( "_sub_domain_",
                          DIMENSION,
                          subdomain_attr );
            }
        }
    }

    int done_local = static_cast < int > ( is_done_ );
    int done_global = 0;
    MPI_Allreduce ( &done_local,
                    &done_global,
                    1,
                    MPI_INT,
                    MPI_LAND,
                    comm_
                    );
    is_done_ = static_cast < bool > ( done_global );

    LOG_INFO ( "Visualize", " adapted mesh " << refinement_level_ );
    PVtkWriter writer ( comm_ );
    std::ostringstream name;
    name << "poisson_adaptive_mesh_" << refinement_level_ << ".pvtu";
    std::string output_file = name.str ( );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

