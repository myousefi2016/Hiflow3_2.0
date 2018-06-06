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

/// \author Simon Gawlok

#include "poisson_precond.h"

static const char* PARAM_FILENAME = "poisson_precond.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class PoissonPrecond
{
  public:

    PoissonPrecond ( const std::string& param_filename, const std::string& path_mesh )
    : path_mesh ( path_mesh ),
    comm_ ( MPI_COMM_WORLD ),
    rank_ ( -1 ),
    num_partitions_ ( -1 ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    rhs_ ( 0 ), sol_ ( 0 ), matrix_ ( 0 ),
    is_done_ ( false ),
    refinement_level_ ( 0 )
    {
        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );
    }

    // Main algorithm

    void run ( )
    {
        // Prepare headers of CSV file
        this->csv_names_.clear ( );
        this->csv_names_.push_back ( "Level" );
        this->csv_names_.push_back ( "Global number of DoFs" );
        this->csv_names_.push_back ( "Assembly time HiFlow LA" );
#ifdef WITH_HYPRE
        this->csv_names_.push_back ( "Assembly time Hypre LA" );
#endif
#ifdef WITH_ILUPP
        this->csv_names_.push_back ( "Time setup of ILU++ preconditioner" );
        this->csv_names_.push_back ( "Time solution CG with ILU++ preconditioner" );
        this->csv_names_.push_back ( "Iterations CG with ILU++ preconditioner" );
#endif
        this->csv_names_.push_back ( "Time setup of Vanka preconditioner" );
        this->csv_names_.push_back ( "Time solution CG with Vanka preconditioner" );
        this->csv_names_.push_back ( "Iterations CG with Vanka preconditioner" );

        this->csv_names_.push_back ( "Time setup of SSOR preconditioner" );
        this->csv_names_.push_back ( "Time solution CG with SSOR preconditioner" );
        this->csv_names_.push_back ( "Iterations CG with SSOR preconditioner" );

#ifdef WITH_ILUPP
        this->csv_names_.push_back ( "Time setup of Naive Parallel (ILU++) preconditioner" );
        this->csv_names_.push_back ( "Time solution CG with Naive Parallel (ILU++) preconditioner" );
        this->csv_names_.push_back ( "Iterations CG with Naive Parallel (ILU++) preconditioner" );
#else
        this->csv_names_.push_back ( "Time setup of Naive Parallel (SSOR) preconditioner" );
        this->csv_names_.push_back ( "Time solution CG with Naive Parallel (SSOR) preconditioner" );
        this->csv_names_.push_back ( "Iterations CG with Naive Parallel (SSOR) preconditioner" );
#endif

#ifdef WITH_HYPRE
        this->csv_names_.push_back ( "Time solution CG with BoomerAMG preconditioner" );
        this->csv_names_.push_back ( "Iterations CG with BoomerAMG preconditioner" );

        this->csv_names_.push_back ( "Time solution BoomerAMG" );
        this->csv_names_.push_back ( "Iterations BoomerAMG" );

        /*this->csv_names_.push_back ( "Time solution BiCGStab with BoomerAMG preconditioner" );
        this->csv_names_.push_back ( "Iterations BiCGStab with BoomerAMG preconditioner" );

        this->csv_names_.push_back ( "Time solution CG with Euclid (ILUp) preconditioner" );
        this->csv_names_.push_back ( "Iterations CG with Euclid (ILUp) preconditioner" );

        this->csv_names_.push_back ( "Time solution GMRES with PILUT preconditioner" );
        this->csv_names_.push_back ( "Iterations GMRES with PILUT preconditioner" );*/
#endif

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

            // Initialize space and linear algebra.
            prepare_system ( );
            // Compute the stiffness matrix and right-hand side.
            assemble_system ( );
            // Solve the linear system.
            solve_system ( );

            // Write CSV output
            if ( this->rank_ == MASTER_RANK )
            {
                this->csv_writer_.write ( this->csv_quantities_ );
            }
            // Visualize the solution.
            visualize ( );
            adapt ( );
        }
    }

    ~PoissonPrecond ( )
    {
        delete matrix_;
        delete sol_;
        delete rhs_;
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
    // Visualize the results.
    void visualize ( );
    // Adapt the space (mesh and/or degree).
    void adapt ( );

    // member variables

    // MPI communicator.
    MPI_Comm comm_;
    // Local process rank and number of processes.
    int rank_, num_partitions_;

    // Parameter data read in from file.
    PropertyTree params_;

    // Local mesh and mesh on master process.
    MeshPtr mesh_, mesh_without_ghost_, master_mesh_;
    // Solution space.
    VectorSpace<double> space_;

    // Linear algebra couplings helper object.
    Couplings<double> couplings_;
    // Vectors for solution and load vector.
    CoupledVector<Scalar>* rhs_, *sol_;
    // System matrix.
    CoupledMatrix<Scalar>* matrix_;

#ifdef WITH_HYPRE
    HypreMatrix<double> matrix_hypre_;
    HypreVector<double> rhs_hypre_, sol_hypre_;
#endif

    // Global assembler.
    StandardGlobalAssembler<double> global_asm_;

    // Flag for stopping adaptive loop.
    bool is_done_;
    // Current refinement level.
    int refinement_level_;

    // Vectors for error norms
    std::vector<Scalar> L2_err_, H1_err_;

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
}; // end class PoissonPrecond

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    // set default parameter file
    std::string param_filename ( PARAM_FILENAME );
    std::string path_mesh;
    // if set take parameter file specified on console
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }
    // if set take mesh following path specified on console
    if ( argc > 2 )
    {
        path_mesh = std::string ( argv[2] );
    }
    try
    {
        // Create log files for INFO and DEBUG output
        //std::ofstream info_log("poisson_tutorial_info_log");
        LogKeeper::get_log ( "info" ).set_target ( &( std::cout ) );
        //std::ofstream debug_log("poisson_tutorial_debug_log");
        LogKeeper::get_log ( "debug" ).set_target ( &( std::cout ) );

        // Create application object and run it
        PoissonPrecond app ( param_filename, path_mesh );
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

//////////////// PoissonPrecond implementation //////////////

void PoissonPrecond::build_initial_mesh ( )
{
#ifdef WITH_PARMETIS
    if ( rank_ == MASTER_RANK )
    {
        std::string mesh_name;

        switch ( DIMENSION )
        {
            case 1:
            {
                mesh_name =
                        params_["Mesh"]["Filename1"].get<std::string>( "unit_line.inp" );
                break;
            }
            case 2:
            {
                mesh_name =
                        params_["Mesh"]["Filename2"].get<std::string>( "unit_square.inp" );
                break;
            }
            case 3:
            {
                mesh_name =
                        params_["Mesh"]["Filename3"].get<std::string>( "unit_cupe.inp" );
                break;
            }

            default: assert ( 0 );
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
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION,
                                             DIMENSION, 0 );

        // Refine the mesh until each process has at least 8 cells on average.
        while ( master_mesh_->num_entities ( DIMENSION ) < 8 * this->num_partitions_
                && refinement_level_ < params_["Mesh"]["FinalRefLevel"].get<int>( 6 ) )
        {
            master_mesh_ = master_mesh_->refine ( );
            ++refinement_level_;
        }
    }

    // 1D parallel execution is not yet implemented.
    if ( DIMENSION == 1 )
    {
        mesh_ = master_mesh_;
    }
    else
    {
        MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );

        mesh_without_ghost_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_ );
        assert ( mesh_without_ghost_ != 0 );
        SharedVertexTable shared_verts;
        mesh_ = compute_ghost_cells ( *mesh_without_ghost_, comm_, shared_verts );

        // Write out mesh of initial refinement level
        /*PVtkWriter writer ( comm_ );
        std::string output_file = std::string ( "poisson_tutorial_mesh.pvtu" );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );*/
    }
#else
    // Read in the mesh on the master process. The mesh is chosen according to the dimension of the problem.
    if ( rank_ == MASTER_RANK )
    {
        std::string mesh_name;

        switch ( DIMENSION )
        {
            case 1:
            {
                mesh_name =
                        params_["Mesh"]["Filename1"].get<std::string>( "unit_line.inp" );
                break;
            }
            case 2:
            {
                mesh_name =
                        params_["Mesh"]["Filename2"].get<std::string>( "unit_square.inp" );
                break;
            }
            case 3:
            {
                mesh_name =
                        params_["Mesh"]["Filename3"].get<std::string>( "unit_cupe.inp" );
                break;
            }

            default: assert ( 0 );
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
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION,
                                             DIMENSION, 0 );

        // Refine the mesh until the initial refinement level is reached.
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( 3 );
        for ( int r = 0; r < initial_ref_lvl; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
            ++refinement_level_;
        }
    }

    // 1D parallel execution is not yet implemented.
    if ( DIMENSION == 1 )
    {
        mesh_ = master_mesh_;
    }
    else
    {
        MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );

        mesh_without_ghost_ = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_ );
        assert ( mesh_without_ghost_ != 0 );
        SharedVertexTable shared_verts;
        mesh_ = compute_ghost_cells ( *mesh_without_ghost_, comm_, shared_verts );

        // Write out mesh of initial refinement level
        PVtkWriter writer ( comm_ );
        std::string output_file = std::string ( "poisson_tutorial_mesh.pvtu" );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );
    }
#endif
}

void PoissonPrecond::prepare_system ( )
{
    // Assign degrees to each element.
    const int fe_degree = params_["Mesh"]["FeDegree"].get<int>( 1 );
    std::vector< int > degrees ( 1, fe_degree );
    std::vector< bool > is_cg ( 1, true );

    // Initialize the VectorSpace object.
    space_.Init ( degrees, *mesh_, is_cg, KING );

    this->csv_quantities_.push_back ( space_.dof ( ).ndofs_global ( ) );

    // Setup couplings object.
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    // Setup linear algebra objects.
    delete matrix_;
    delete rhs_;
    delete sol_;

    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( "CoupledMatrix" ) )->
            params ( params_["LinearAlgebra"] );
    matrix_->Init ( comm_, couplings_ );
    CoupledVectorFactory<Scalar> CoupVecFact;
    rhs_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( "CoupledVector" ) )->
            params ( params_["LinearAlgebra"] );
    sol_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( "CoupledVector" ) )->
            params ( params_["LinearAlgebra"] );
    rhs_->Init ( comm_, couplings_ );
    sol_->Init ( comm_, couplings_ );

    // Initialize structure of LA objects.
    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );

    rhs_->InitStructure ( );
    sol_->InitStructure ( );

    // Zero all linear algebra objects.
    matrix_->Zeros ( );
    rhs_->Zeros ( );
    sol_->Zeros ( );

#ifdef WITH_HYPRE
    LOG_INFO ( "Also testing Hypre interface", true );
    matrix_hypre_.Init ( comm_, couplings_ );
    matrix_hypre_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                                  vec2ptr ( sparsity.diagonal_cols ),
                                  sparsity.diagonal_rows.size ( ),
                                  vec2ptr ( sparsity.off_diagonal_rows ),
                                  vec2ptr ( sparsity.off_diagonal_cols ),
                                  sparsity.off_diagonal_rows.size ( ) );
    //matrix_hypre_.print_statistics();

    rhs_hypre_.Init ( comm_, couplings_ );
    sol_hypre_.Init ( comm_, couplings_ );
    //rhs_hypre_.print_statistics();
    //sol_hypre_.print_statistics();
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
        for ( EntityIterator facet_it = mesh_->begin ( DIMENSION - 1 ), facet_end = mesh_->end ( DIMENSION - 1 );
              facet_it != facet_end; ++facet_it )
        {

            // Returns the number of neighbors for each cell, to check if it is on the facet.
            const EntityCount num_cell_neighbors = facet_it->num_incident_entities ( DIMENSION );

            // If it lies on the facet, the corresponding DOF is a Dirichlet DOF and is added to dirichlet_dofs_.
            if ( num_cell_neighbors == 1 )
            {
                std::vector<int> dof_number_;
                space_.dof ( ).get_dofs_on_subentity ( 0, facet_it->begin_incident ( DIMENSION )->index ( ), 0, facet_it->index ( ), dof_number_ );
                dirichlet_dofs_.push_back ( dof_number_[0] );
            }
        }
    }
    else
    {
        compute_dirichlet_dofs_and_values ( zero, space_, 0, dirichlet_dofs_,
                                            dirichlet_values_ );
    }
}

void PoissonPrecond::assemble_system ( )
{
    // Assemble matrix and right-hand-side vector.
    LocalPoissonAssembler local_asm;
    Timer assembly_hf_la;
    assembly_hf_la.start ( );
    global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, *rhs_ );
    assembly_hf_la.stop ( );
    LOG_INFO ( "Assembly time with HiFlow LA", assembly_hf_la.get_duration ( ) );
    this->csv_quantities_.push_back ( assembly_hf_la.get_duration ( ) );

#ifdef WITH_HYPRE
    Timer assembly_hypre_la;
    assembly_hypre_la.start ( );
    global_asm_.assemble_matrix ( space_, local_asm, matrix_hypre_ );
    global_asm_.assemble_vector ( space_, local_asm, rhs_hypre_ );
    assembly_hypre_la.stop ( );
    LOG_INFO ( "Assembly time with Hypre LA", assembly_hypre_la.get_duration ( ) );
    this->csv_quantities_.push_back ( assembly_hypre_la.get_duration ( ) );
#endif

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
#ifdef WITH_HYPRE
        matrix_hypre_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
        rhs_hypre_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                               vec2ptr ( dirichlet_values_ ) );
        sol_hypre_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                               vec2ptr ( dirichlet_values_ ) );
#endif
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }

    rhs_->UpdateCouplings ( );
    sol_->UpdateCouplings ( );
}

void PoissonPrecond::solve_system ( )
{
    LinearSolver<LAD>* solver_;
    LinearSolverFactory<LAD> SolFact;
    solver_ = SolFact.Get (
                            params_["LinearSolver"]["Name"].get<std::string>( "CG" ) )->
            params ( params_["LinearSolver"] );
    solver_->SetupOperator ( *matrix_ );

    //********************************
    // Test different preconditioners
    //********************************

#ifdef WITH_ILUPP
    if ( rank_ == MASTER_RANK )
    {
        std::cout << "********************************" << std::endl;
        std::cout << "ILU++" << std::endl;
        std::cout << "********************************" << std::endl;
    }
    sol_->Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }

    sol_->UpdateCouplings ( );
    PreconditionerIlupp<LAD> precond_ilupp;
    precond_ilupp.InitParameter (
                                  params_["Preconditioner"]["ILUPP"]["PreprocessingType"].get<int>( ),
                                  params_["Preconditioner"]["ILUPP"]["PreconditionerNumber"].get<int>( ),
                                  params_["Preconditioner"]["ILUPP"]["MaxMultilevels"].get<int>( ),
                                  params_["Preconditioner"]["ILUPP"]["MemFactor"].get<Scalar>( ),
                                  params_["Preconditioner"]["ILUPP"]["PivotThreshold"].get<Scalar>( ),
                                  params_["Preconditioner"]["ILUPP"]["MinPivot"].get<Scalar>( )
                                  );

    Timer setup_ilupp;
    setup_ilupp.start ( );

    precond_ilupp.SetupOperator ( *matrix_ );
    precond_ilupp.Build ( );

    setup_ilupp.stop ( );
    LOG_INFO ( "Time setup of ILU++ preconditioner", setup_ilupp.get_duration ( ) );
    this->csv_quantities_.push_back ( setup_ilupp.get_duration ( ) );

    solver_->SetupPreconditioner ( precond_ilupp );

    Timer solve_ilupp;
    solve_ilupp.start ( );

    solver_->Solve ( *rhs_, sol_ );

    solve_ilupp.stop ( );
    LOG_INFO ( "Time solution with ILU++ preconditioner", solve_ilupp.get_duration ( ) );
    this->csv_quantities_.push_back ( solve_ilupp.get_duration ( ) );
    this->csv_quantities_.push_back ( solver_->iter ( ) );

    /*	
            precond_ilupp.SetupOperator ( *matrix_ );
            solve_ilupp.reset();
            solve_ilupp.start ( );
            sol_->Zeros ( );
        if ( !dirichlet_dofs_.empty ( ) )
        {
            // Correct Dirichlet dofs.
            sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                              vec2ptr ( dirichlet_values_ ) );
        }

        sol_->UpdateCouplings ( );

        solver_->Solve ( *rhs_, sol_ );

        solve_ilupp.stop ( );
        LOG_INFO ( "Time solution with ILU++ preconditioner", solve_ilupp.get_duration ( ) );
     */
    precond_ilupp.Clear ( );
#endif
    if ( rank_ == MASTER_RANK )
    {
        std::cout << "********************************" << std::endl;
        std::cout << "Vanka" << std::endl;
        std::cout << "********************************" << std::endl;
    }
    sol_->Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }

    sol_->UpdateCouplings ( );
    PreconditionerVanka<LAD> precond_vanka;
    precond_vanka.InitParameter (
                                  space_,
                                  params_["Preconditioner"]["Vanka"]["DampingParam"].get<Scalar>( ),
                                  params_["Preconditioner"]["Vanka"]["NumIter"].get<int>( )
                                  );

    Timer setup_vanka;
    setup_vanka.start ( );

    precond_vanka.SetupOperator ( *matrix_ );
    precond_vanka.Build ( );

    setup_vanka.stop ( );
    LOG_INFO ( "Time setup of Vanka preconditioner", setup_vanka.get_duration ( ) );
    this->csv_quantities_.push_back ( setup_vanka.get_duration ( ) );

    solver_->SetupPreconditioner ( precond_vanka );

    Timer solve_vanka;
    solve_vanka.start ( );

    solver_->Solve ( *rhs_, sol_ );

    solve_vanka.stop ( );
    LOG_INFO ( "Time solution with Vanka preconditioner", solve_vanka.get_duration ( ) );
    this->csv_quantities_.push_back ( solve_vanka.get_duration ( ) );
    this->csv_quantities_.push_back ( solver_->iter ( ) );

    precond_vanka.Clear ( );

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "********************************" << std::endl;
        std::cout << "SSOR" << std::endl;
        std::cout << "********************************" << std::endl;
    }
    sol_->Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }

    sol_->UpdateCouplings ( );
    PreconditionerBlockJacobiStand<LAD> precond_ssor;
    precond_ssor.Init_SSOR ( params_["Preconditioner"]["SSOR"].get<Scalar>( ) );

    Timer setup_ssor;
    setup_ssor.start ( );

    precond_ssor.SetupOperator ( *matrix_ );
    precond_ssor.Build ( );

    setup_ssor.stop ( );
    LOG_INFO ( "Time setup of SSOR preconditioner", setup_ssor.get_duration ( ) );
    this->csv_quantities_.push_back ( setup_ssor.get_duration ( ) );

    solver_->SetupPreconditioner ( precond_ssor );

    Timer solve_ssor;
    solve_ssor.start ( );

    solver_->Solve ( *rhs_, sol_ );

    solve_ssor.stop ( );
    LOG_INFO ( "Time solution with SSOR preconditioner", solve_ssor.get_duration ( ) );
    this->csv_quantities_.push_back ( solve_ssor.get_duration ( ) );
    this->csv_quantities_.push_back ( solver_->iter ( ) );

    precond_ssor.Clear ( );

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "********************************" << std::endl;
        std::cout << "Naive parallel preconditioner" << std::endl;
        std::cout << "********************************" << std::endl;
    }
    sol_->Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }

    sol_->UpdateCouplings ( );
    PreconditionerParallelNaive<LAD> precond_naive;
    precond_naive.InitParameter ( params_["Preconditioner"]["Naive"]["NumIter"].get<int>( ) );

    Timer setup_naive;
    setup_naive.start ( );

#ifdef WITH_ILUPP
    precond_naive.SetPreconditioner ( precond_ilupp );
#else
    precond_ssor.Init_SSOR ( params_["Preconditioner"]["SSOR"].get<Scalar>( ) );
    precond_naive.SetPreconditioner ( precond_ssor );
#endif
    precond_naive.SetupOperator ( *matrix_ );
    precond_naive.Build ( );

    setup_naive.stop ( );
    LOG_INFO ( "Time setup of Naive parallel preconditioner", setup_naive.get_duration ( ) );
    this->csv_quantities_.push_back ( setup_naive.get_duration ( ) );

    solver_->SetupPreconditioner ( precond_naive );

    Timer solve_naive;
    solve_naive.start ( );

    solver_->Solve ( *rhs_, sol_ );

    solve_naive.stop ( );
    LOG_INFO ( "Time solution with Naive parallel preconditioner", solve_naive.get_duration ( ) );
    this->csv_quantities_.push_back ( solve_naive.get_duration ( ) );
    this->csv_quantities_.push_back ( solver_->iter ( ) );

    precond_naive.Clear ( );

#ifdef WITH_HYPRE
    HypreCG<LAD_H> solver_hypre;
    //	CG<LAD_H> solver_hypre;

    solver_hypre.SetupOperator ( matrix_hypre_ );
    solver_hypre.InitControl ( params_["LinearSolver"]["MaxIterations"].get<int>( ),
                               params_["LinearSolver"]["AbsTolerance"].get<Scalar>( ),
                               params_["LinearSolver"]["RelTolerance"].get<Scalar>( ) );
    //							   params_["LinearSolver"]["RelTolerance"].get<Scalar>( ), 1e6);

    Timer solve_hypre;
    Timer setup_hypre;

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "********************************" << std::endl;
        std::cout << "HypreCG with BoomerAMG" << std::endl;
        std::cout << "********************************" << std::endl;
    }

    sol_hypre_.Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        sol_hypre_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                               vec2ptr ( dirichlet_values_ ) );
    }
    sol_hypre_.Update ( );

    solver_hypre.Clear ( );
    solver_hypre.SetupOperator ( matrix_hypre_ );
    solver_hypre.InitControl ( params_["LinearSolver"]["MaxIterations"].get<int>( ),
                               params_["LinearSolver"]["AbsTolerance"].get<Scalar>( ),
                               params_["LinearSolver"]["RelTolerance"].get<Scalar>( ) );
    //								params_["LinearSolver"]["RelTolerance"].get<Scalar>( ), 1e6);

    HypreBoomerAMG<LAD_H> precond_boomer;
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
    precond_boomer.SetupOperator ( matrix_hypre_ );
    precond_boomer.SetReuse ( true );

    Timer solve_hypre_boomer;
    solve_hypre_boomer.start ( );

    precond_boomer.Build ( );

    solver_hypre.SetupPreconditioner ( precond_boomer );

    solver_hypre.Solve ( rhs_hypre_, &sol_hypre_ );

    solve_hypre_boomer.stop ( );

    LOG_INFO ( "Iterations with Hypre CG (BoomerAMG preconditioner)", solver_hypre.iter ( ) );
    LOG_INFO ( "Relative final residual with Hypre CG (BoomerAMG preconditioner)", solver_hypre.res ( ) );
    LOG_INFO ( "Time solution with Hypre CG (BoomerAMG preconditioner)", solve_hypre_boomer.get_duration ( ) );
    this->csv_quantities_.push_back ( solve_hypre_boomer.get_duration ( ) );
    this->csv_quantities_.push_back ( solver_hypre.iter ( ) );
    precond_boomer.Clear ( );

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "********************************" << std::endl;
        std::cout << "BoomerAMG (stand-alone)" << std::endl;
        std::cout << "********************************" << std::endl;
    }

    sol_hypre_.Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        sol_hypre_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                               vec2ptr ( dirichlet_values_ ) );
    }

    HypreBoomerAMG<LAD_H> hypre_boomer;
    hypre_boomer.InitControl ( params_["LinearSolver"]["MaxIterations"].get<int>( ),
                               params_["LinearSolver"]["AbsTolerance"].get<Scalar>( ),
                               params_["LinearSolver"]["RelTolerance"].get<Scalar>( ) );
    hypre_boomer.SetCycleType ( 2 );
    hypre_boomer.SetRelaxType ( 6 );
    hypre_boomer.SetCycleNumSweeps ( 3, 1 );
    hypre_boomer.SetCycleNumSweeps ( 3, 2 );
    hypre_boomer.SetRelaxWt ( 0.5 );
    hypre_boomer.SetCoarsenType ( 10 );
    hypre_boomer.SetInterpType ( 6 );
#    if DIMENSION == 2
    hypre_boomer.SetStrongThreshold ( 0.25 );
#    else
    hypre_boomer.SetStrongThreshold ( 0.6 );
#    endif
    hypre_boomer.SetAggNumLevels ( 25 );
    hypre_boomer.SetSmoothType ( 0 );

    hypre_boomer.SetupOperator ( matrix_hypre_ );

    Timer hypre_boomer_timer;
    hypre_boomer_timer.start ( );

    hypre_boomer.Solve ( rhs_hypre_, &sol_hypre_ );

    hypre_boomer_timer.stop ( );
    LOG_INFO ( "Iterations with Hypre BoomerAMG", hypre_boomer.iter ( ) );
    LOG_INFO ( "Relative final residual with Hypre BoomerAMG", hypre_boomer.res ( ) );
    LOG_INFO ( "Time solution with Hypre BoomerAMG", hypre_boomer_timer.get_duration ( ) );
    this->csv_quantities_.push_back ( hypre_boomer_timer.get_duration ( ) );
    this->csv_quantities_.push_back ( hypre_boomer.iter ( ) );

    /*if ( rank_ == MASTER_RANK )
    {
        std::cout << "********************************" << std::endl;
        std::cout << "HypreCG with Euclid" << std::endl;
        std::cout << "********************************" << std::endl;
    }

    sol_hypre_.Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        sol_hypre_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                               vec2ptr ( dirichlet_values_ ) );
    }

    solver_hypre.Clear ( );
    solver_hypre.SetupOperator ( matrix_hypre_ );
    solver_hypre.InitControl ( params_["LinearSolver"]["MaxIterations"].get<int>( ),
                               params_["LinearSolver"]["AbsTolerance"].get<Scalar>( ),
                               params_["LinearSolver"]["RelTolerance"].get<Scalar>( ) );

    HyprePreconditionerEuclid<LAD_H> precond_euclid ( comm_ );
    //precond_euclid.SetLevel(0);
    precond_euclid.SetBJ ( 0 );
    precond_euclid.SetSparseA ( 1.e-3 );
    precond_euclid.SetRowScale ( 1 );
    //precond_euclid.SetILUT(1.e-3);
    solver_hypre.SetupPreconditioner ( precond_euclid );

    solve_hypre_boomer.start ( );

    solver_hypre.Solve ( rhs_hypre_, &sol_hypre_ );

    solve_hypre_boomer.stop ( );
    LOG_INFO ( "Iterations with Hypre CG (Euclid preconditioner)", solver_hypre.iter ( ) );
    LOG_INFO ( "Relative final residual with Hypre CG (Euclid preconditioner)", solver_hypre.res ( ) );
    LOG_INFO ( "Time solution with Hypre CG (Euclid preconditioner)", solve_hypre_boomer.get_duration ( ) );
    this->csv_quantities_.push_back ( solve_hypre_boomer.get_duration ( ) );
    this->csv_quantities_.push_back ( solver_hypre.iter ( ) );

    precond_euclid.Clear ( );

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "********************************" << std::endl;
        std::cout << "HypreGMRES with PILUT" << std::endl;
        std::cout << "********************************" << std::endl;
    }

    sol_hypre_.Zeros ( );
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        sol_hypre_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                               vec2ptr ( dirichlet_values_ ) );
    }

    HypreGMRES<LAD_H> solver_hypre_pilut;
    solver_hypre_pilut.SetupOperator ( matrix_hypre_ );
    solver_hypre_pilut.InitControl ( params_["LinearSolver"]["MaxIterations"].get<int>( ),
                                     params_["LinearSolver"]["AbsTolerance"].get<Scalar>( ),
                                     params_["LinearSolver"]["RelTolerance"].get<Scalar>( ) );
    solver_hypre_pilut.InitParameter ( 250 );

    HyprePreconditionerPILUT<LAD_H> precond_pilut ( comm_ );
    precond_pilut.SetMaxIter ( 1000 );
    precond_pilut.SetDropTol ( 1.e-4 );
    precond_pilut.SetFactorRowSize ( 20 );

    solver_hypre_pilut.SetupPreconditioner ( precond_pilut );

    Timer solve_hypre_pilut;
    solve_hypre_pilut.start ( );

    solver_hypre_pilut.Solve ( rhs_hypre_, &sol_hypre_ );

    solve_hypre_pilut.stop ( );
    LOG_INFO ( "Iterations with Hypre GMRES (PILUT preconditioner)", solver_hypre_pilut.iter ( ) );
    LOG_INFO ( "Relative final residual with Hypre GMRES (PILUT preconditioner)", solver_hypre_pilut.res ( ) );
    LOG_INFO ( "Time solution with Hypre GMRES (PILUT preconditioner)", solve_hypre_pilut.get_duration ( ) );
    this->csv_quantities_.push_back ( solve_hypre_pilut.get_duration ( ) );
    this->csv_quantities_.push_back ( solver_hypre_pilut.iter ( ) );

    precond_pilut.Clear ( );*/
#endif

    delete solver_;
}

void PoissonPrecond::visualize ( )
{
    // Setup visualization object.
    const int num_intervals = params_["Mesh"]["FeDegree"].get<int>( 1 );
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    // Generate filename.
    std::stringstream name;
    name << "solution" << refinement_level_;

    std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> material_number ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

    for ( mesh::EntityIterator it = mesh_->begin ( mesh_->tdim ( ) );
          it != mesh_->end ( mesh_->tdim ( ) );
          ++it )
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

#ifdef WITH_HYPRE
    sol_hypre_.Update ( );
    visu.visualize ( EvalFeFunction<LAD_H>( space_, sol_hypre_ ), "u" );
#else
    sol_->Update ( );
    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ) ), "u" );
#endif

    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.write ( name.str ( ) );
}

void PoissonPrecond::adapt ( )
{
    // Refine mesh on master process. 1D parallel execution is not yet implemented.
    if ( DIMENSION == 1 )
    {
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
            if ( refinement_level_ >= final_ref_level )
            {
                is_done_ = true;
            }
            else
            {
                mesh_ = mesh_->refine ( );
                ++refinement_level_;
            }
        }
    }
    else
    {
#ifdef WITH_PARMETIS
        const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
        if ( refinement_level_ >= final_ref_level )
        {
            is_done_ = true;
        }
        else
        {
            mesh_without_ghost_ = mesh_without_ghost_->refine ( );
            ++refinement_level_;
        }

        if ( !is_done_ )
        {
            // Repartition the new mesh.
            ParMetisGraphPartitioner parmetis_partitioner;
            MeshPtr local_mesh = repartition_mesh ( mesh_without_ghost_, comm_, &parmetis_partitioner );
            assert ( local_mesh != 0 );
            mesh_.reset ( );
            SharedVertexTable shared_verts;
            mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
            local_mesh.reset ( );
        }
#else
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( 6 );
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
        // Broadcast information from master to slaves.
        MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
        MPI_Bcast ( &is_done_, 1, MPI_CHAR, MASTER_RANK, comm_ );

        if ( !is_done_ )
        {
            // Distribute the new mesh.
            MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_ );
            assert ( local_mesh != 0 );
            SharedVertexTable shared_verts;
            mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
        }
#endif
    }
}
