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

#include "poisson_eigenvalue.h"
#include "eigen_value/slepc_environment.h"

static const char* PARAM_FILENAME = "poisson_eigenvalue.xml";
#define MESHES_DATADIR "../data/"

static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class PoissonEigenValue
{
  public:

    PoissonEigenValue ( const std::string& param_filename, const std::string& path_mesh )
    : path_mesh ( path_mesh ),
    comm_ ( MPI_COMM_WORLD ),
    rank_ ( -1 ),
    num_partitions_ ( -1 ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    A_ ( 0 ),
    B_ ( 0 ),
    is_done_ ( false ),
    refinement_level_ ( 0 )
    {
        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );
    }

    // Main algorithm

    void run ( )
    {
        // Construct / read in the initial mesh.
        build_initial_mesh ( );

        // Main adaptation loop.
        while ( !is_done_ )
        {
            // Initialize space and linear algebra.
            prepare_system ( );

            // Compute the stiffness matrix and right-hand side.
            assemble_system ( );

            // Solve the linear system.
            compute_eigenvalues ( );

            // Modify the space through refinement. Set is_done_ = true when finished.
            adapt ( );

        }
    }

    ~PoissonEigenValue ( )
    {
        delete A_;
        delete B_;
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

    // Compute Eigenvalues
    void compute_eigenvalues ( );

    // Compute errors compared to exact solution.
    void compute_error ( );

    // Visualize the results.
    void visualize ( VectorType* vec, std::string vec_name );

    // Helper functions for residual computation
    void SetVectorsToZero ( );
    void UpdatetVectorCouplings ( );

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
    MeshPtr mesh_, master_mesh_;

    // Solution space.
    VectorSpace<double> space_;

    // Linear algebra couplings helper object.
    Couplings<double> couplings_;
    SparsityStructure sparsity_;

    VectorType* vec_ref_;

    // Real and imaginary part of computed eigenvalues
    std::vector< Scalar > k_real_;
    std::vector< Scalar > k_imag_;

    // Real and imaginary part of computed eigenvectors
    std::vector< VectorType* > v_real_;
    std::vector< VectorType* > v_imag_;

    VectorType Arxr_;
    VectorType Aixi_;
    VectorType Arxi_;
    VectorType Aixr_;
    VectorType Bxr_;
    VectorType Bxi_;

    VectorType res_r_;
    VectorType res_i_;

    // System matrix.
    MatrixType* A_;
    MatrixType* Ai_;
    MatrixType* B_;

    // Global assembler.
    StandardGlobalAssembler<double> global_asm_;

    // Flag for stopping adaptive loop.
    bool is_done_;
    // Current refinement level.
    int refinement_level_;

    // Dof id:s for Dirichlet boundary conditions.
    std::vector<int> dirichlet_dofs_;
    // Dof values for Dirichlet boundary conditions.
    std::vector<Scalar> dirichlet_values_;
}; // end class PoissonEigenValue

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );
    SLEPcEnvironment::initialize ( argc, argv );

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
        std::ofstream info_log ( "poisson_eigenvalue_info_log" );
        LogKeeper::get_log ( "info" ).set_target ( &info_log );
        std::ofstream debug_log ( "poisson_eigenvalue_error_log" );
        LogKeeper::get_log ( "error" ).set_target ( &debug_log );

        // Create application object and run it
        PoissonEigenValue app ( param_filename, path_mesh );
        app.run ( );

    }
    catch ( std::exception& e )
    {
        std::cerr << "\nProgram ended with uncaught exception.\n";
        std::cerr << e.what ( ) << "\n";
        return -1;
    }
    SLEPcEnvironment::finalize ( );
    MPI_Finalize ( );
    return 0;
}

//////////////// PoissonEigenValue implementation //////////////

void PoissonEigenValue::build_initial_mesh ( )
{
    // Read in the mesh on the master process. The mesh is chosen according to the dimension of the problem.

    if ( rank_ == MASTER_RANK )
    {

        std::string mesh_name;

        mesh_name = params_["Mesh"]["Filename"].get<std::string>( );

        std::string mesh_filename;
        if ( path_mesh.empty ( ) )
        {
            mesh_filename = std::string ( DATADIR ) + mesh_name;
        }
        else
        {
            mesh_filename = path_mesh + mesh_name;
        }

        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0 );

        // Refine the mesh until the initial refinement level is reached.
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );
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

        MeshPtr local_mesh;

        // Distribute mesh over all processes, and compute ghost cells
        if ( num_partitions_ > 1 )
        {
#ifdef WITH_METIS
            MetisGraphPartitioner partitioner;
#else
            NaiveGraphPartitioner partitioner;
#endif
            local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &partitioner );
        }
        else
        {
            NaiveGraphPartitioner partitioner;
            local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &partitioner );
        }

        assert ( local_mesh != 0 );
        SharedVertexTable shared_verts;

        mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );

        PVtkWriter writer ( comm_ );
        std::string output_file = std::string ( "poisson_eigenvalue_mesh.pvtu" );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );
    }
}

void PoissonEigenValue::prepare_system ( )
{

    A_ = new MatrixType;
    Ai_ = new MatrixType;
    B_ = new MatrixType;
    vec_ref_ = new VectorType;

    // Assign degrees to each element.
    const int fe_degree = params_["Mesh"]["FeDegree"].get<int>( );
    std::vector< int > degrees ( 1, fe_degree );

    // Initialize the VectorSpace object.
    space_.Init ( degrees, *mesh_ );

    // Setup couplings object.
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph.
    global_asm_.compute_sparsity_structure ( space_, sparsity_ );

    couplings_.InitializeCouplings ( sparsity_.off_diagonal_rows, sparsity_.off_diagonal_cols );

    // Setup linear algebra objects.
#ifdef USE_PETSC
    A_->Init ( comm_, couplings_ );
    Ai_->Init ( comm_, couplings_ );
    B_->Init ( comm_, couplings_ );
#else
#    ifdef USE_HYPRE
    A_->Init ( comm_, couplings_ );
    Ai_->Init ( comm_, couplings_ );
    B_->Init ( comm_, couplings_ );
#    else
    A_->Init ( comm_, couplings_, CPU, NAIVE, CSR );
    Ai_->Init ( comm_, couplings_, CPU, NAIVE, CSR );
    B_->Init ( comm_, couplings_, CPU, NAIVE, CSR );
#    endif
#endif

    A_->InitStructure ( vec2ptr ( sparsity_.diagonal_rows ),
                        vec2ptr ( sparsity_.diagonal_cols ),
                        sparsity_.diagonal_rows.size ( ),
                        vec2ptr ( sparsity_.off_diagonal_rows ),
                        vec2ptr ( sparsity_.off_diagonal_cols ),
                        sparsity_.off_diagonal_rows.size ( ) );

    Ai_->InitStructure ( vec2ptr ( sparsity_.diagonal_rows ),
                         vec2ptr ( sparsity_.diagonal_cols ),
                         sparsity_.diagonal_rows.size ( ),
                         vec2ptr ( sparsity_.off_diagonal_rows ),
                         vec2ptr ( sparsity_.off_diagonal_cols ),
                         sparsity_.off_diagonal_rows.size ( ) );

    B_->InitStructure ( vec2ptr ( sparsity_.diagonal_rows ),
                        vec2ptr ( sparsity_.diagonal_cols ),
                        sparsity_.diagonal_rows.size ( ),
                        vec2ptr ( sparsity_.off_diagonal_rows ),
                        vec2ptr ( sparsity_.off_diagonal_cols ),
                        sparsity_.off_diagonal_rows.size ( ) );

#ifdef USE_PETSC
    vec_ref_->Init ( comm_, couplings_ );
#else
#    ifdef USE_HYPRE
    vec_ref_->Init ( comm_, couplings_ );
#    else
    vec_ref_->Init ( comm_, couplings_, CPU, NAIVE );
    vec_ref_->InitStructure ( );
#    endif
#endif

    // Zero all linear algebra objects.
    A_->Zeros ( );
    Ai_->Zeros ( );
    B_->Zeros ( );
    vec_ref_->Zeros ( );

    // Compute Dirichlet BC dofs and values using known exact solution.
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    DirichletZero zero;

    // The function compute_dirichlet_dofs_and_values des not yet work for 1D.
    if ( DIMENSION == 1 )
    {
        dirichlet_values_.resize ( 2, 0.0 );

        // Loop over all cells.
        for ( EntityIterator facet_it = mesh_->begin ( DIMENSION - 1 ), facet_end = mesh_->end ( DIMENSION - 1 ); facet_it != facet_end; ++facet_it )
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
        compute_dirichlet_dofs_and_values ( zero, space_, 0, dirichlet_dofs_, dirichlet_values_ );
    }

    Arxr_.CloneFrom ( *vec_ref_ );
    Aixi_.CloneFrom ( *vec_ref_ );
    Arxi_.CloneFrom ( *vec_ref_ );
    Aixr_.CloneFrom ( *vec_ref_ );
    Bxr_ .CloneFrom ( *vec_ref_ );
    Bxi_.CloneFrom ( *vec_ref_ );

    res_r_.CloneFrom ( *vec_ref_ );
    res_i_.CloneFrom ( *vec_ref_ );
}

void PoissonEigenValue::assemble_system ( )
{
    // Assemble matrix and right-hand-side vector.
    LocalPoissonAssembler local_asm;
    local_asm.set_mode_to_real ( );
    global_asm_.assemble_matrix ( space_, local_asm, *A_ );

#ifdef COMPLEX_OPERATOR
#    ifdef WITH_COMPLEX_PETSC
    local_asm.set_mode_to_imag ( );
    global_asm_.assemble_matrix ( space_, local_asm, *Ai_ );
#    endif
#endif

    LocalMassAssembler local_mass_asm;
    global_asm_.assemble_matrix ( space_, local_mass_asm, *B_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        A_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
        Ai_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
        B_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
    }
}

void PoissonEigenValue::compute_eigenvalues ( )
{

#ifdef COMPLEX_OPERATOR
#    ifndef WITH_COMPLEX_PETSC
    std::cout << "CAUTION: Need complex version of PETSc for complex eigenvalue system -> Ignore imaginary part of operator " << std::endl;
#    endif
#endif

    // Type of method: ARNOLDI, KRYLOVSCHUR, LANCZOS, RQCG, GD, JD, POWER
    slepc::EpsType solver_type = slepc::KRYLOVSCHUR;

    // Type of eigenvalue problem:
    // HEP (Hermitian), NHEP (Non-Hermitian), GHEP (generalized Hermitian), GHIEP (Generalized indefinite Hermitian)
    // GNHEP (Generalized Non-Hermitian) , PGNHEP (GNHEP with positive (semi-) definite B  )
#ifdef ADD_CONVECTION
#    ifdef WITH_COMPLEX_PETSC
    slepc::EpsProblemType problem_type = slepc::GNHEP; // TODO: in principle, PGNHEP should be possible as well. However, this choice would lead to a SLEPc error
#    else
    slepc::EpsProblemType problem_type = slepc::PGNHEP;
#    endif
#else
#    ifdef WITH_COMPLEX_PETSC
    slepc::EpsProblemType problem_type = slepc::GNHEP; // TODO: in principle, GHEP should be possible as well. However, this choice would lead to a SLEPc error
#    else
    slepc::EpsProblemType problem_type = slepc::GHEP;
#    endif
#endif

    // Type of Desired Extremal Eigenvalue: L = Largest, S= Smallest, Mag = Magnitude, REAL = real part, IMAG = imaginary part
    slepc::EpsWhich eigenvalue_type = slepc::S_MAG;

    // Measure of convergence and error: REL, ABS, NORM
    slepc::EpsConv conv_type = slepc::REL;
    slepc::EpsConv error_type = slepc::REL;

    // Distance measure for Eigenvalue closest to target value: MAG, REAL, IMAG
    slepc::EpsMeasure meas_type = slepc::MAG;

    // Speactral Transformation type: SHIFT, SINVERT, CAYLEY, PRECOND
    slepc::StType st_type = slepc::SHIFT;

    // Linear solver used in Spectraltransformation ST, applied to B^{-1} (st_type = SHift) or (A-sig B)^{-1} (st_type = sinvert)
    slepc::KSPType st_ksp_type = slepc::PREONLY;
    slepc::PcType st_pc_type = slepc::LU;
    slepc::SolverPackage st_package = slepc::NOPACKAGE; // Use PETSc internal LU decomposition
    int st_ksp_max_iter = 1;
    int st_ksp_size = 1;
    double st_ksp_abs_tol = 1e-10;
    double st_ksp_rel_tol = 1e-6;

    Scalar target = params_["EigenSolver"]["TargetValue"].get<double>( );
    Scalar a = params_["EigenSolver"]["LowerBorder"].get<double>( );
    Scalar b = params_["EigenSolver"]["UpperBorder"].get<double>( );
    Scalar atol = params_["EigenSolver"]["AbsTolerance"].get<double>( );
    Scalar rtol = params_["EigenSolver"]["RelTolerance"].get<double>( );
    int maxits = params_["EigenSolver"]["MaxIterations"].get<int>( );
    int nev = params_["EigenSolver"]["NumberOfEigenvalues"].get<int>( );
    int ncv = params_["EigenSolver"]["MaxSubSpaceDimension"].get<int>( );
    int mpd = params_["EigenSolver"]["MaxProjectedDimension"].get<int>( );

    Timer timer;
#ifdef MATRIXFREE
    SLEPcGeneralEPS<LAD> solver ( solver_type );
#else
#    ifdef MATMIXED
    SLEPcGeneralEPS<LAD> solver ( solver_type );
#    else
    SLEPcGeneralEPS<LAD> solver ( solver_type );
#    endif
#endif

    timer.start ( );
    solver.Init ( comm_ );
    timer.stop ( );
    if ( rank_ == 0 )
        std::cout << "EigenSolver initialization took " << timer.get_duration ( ) << " sec " << std::endl;
    timer.reset ( );
    timer.start ( );

#ifdef MATRIXFREE
    // Matrix for EV problem A*x = k*x
    //	solver.SetupOperator(A_, A_);

    // Matrices for generalized EV problem A*x=k*M*x
    // NOTE: Some eigensolvers require multiplication with the transposed matrix ->  need to provide transposed operators in 3rd and 4th argument
    // NOTE: Its beneficial to provide B explicitly
    // NOTE: 3rd and 4th argument: imaginary part of operator and transposed operator
#    ifdef WITH_COMPLEX_PETSC
    solver.SetupOperatorA ( A_, A_, Ai_, Ai_, couplings_ );
#    else
    solver.SetupOperatorA ( A_, A_, NULL, NULL, couplings_ );
#    endif
    solver.SetupOperatorB ( B_, B_, NULL, NULL, couplings_ );

#else
#    ifdef MATMIXED

    // Matrices for generalized EV problem A*x=k*M*x
    // NOTE: 3rd and 4th argument: imaginary part of operator and transposed operator
#        ifdef WITH_COMPLEX_PETSC
    solver.SetupOperatorA ( A_, A_, Ai_, Ai_, couplings_ );
#        else
    solver.SetupOperatorA ( A_, A_, NULL, NULL, couplings_ );
#        endif

    // NOTE: 2nd argument: imaginary part of operator
    solver.SetupOperatorB ( B_, NULL, couplings_, sparsity_.diagonal_rows, sparsity_.diagonal_cols, sparsity_.off_diagonal_rows, sparsity_.off_diagonal_cols );
#    else
    // Matrices for EV problem A*x = k*x
    //	 solver.SetupOperator(A_);

    // Matrices for generalized EV problem A*x=k*M*x
    // NOTE: 2nd argument: imaginary part of operator
#        ifdef WITH_COMPLEX_PETSC
    solver.SetupOperatorA ( A_, Ai_, couplings_, sparsity_.diagonal_rows, sparsity_.diagonal_cols, sparsity_.off_diagonal_rows, sparsity_.off_diagonal_cols );
#        else
    solver.SetupOperatorA ( A_, NULL, couplings_, sparsity_.diagonal_rows, sparsity_.diagonal_cols, sparsity_.off_diagonal_rows, sparsity_.off_diagonal_cols );
#        endif

    solver.SetupOperatorB ( B_, NULL, couplings_, sparsity_.diagonal_rows, sparsity_.diagonal_cols, sparsity_.off_diagonal_rows, sparsity_.off_diagonal_cols );
#    endif
#endif

    // Tell eigenvalue solver that all reauired operators are set
    solver.PassOperatorsToSolver ( );
    timer.stop ( );
    if ( rank_ == 0 )
        std::cout << "EigenSolver operator setup took " << timer.get_duration ( ) << " sec " << std::endl;

    // Solver control
    solver.InitControl ( nev, ncv, mpd, maxits, rtol, conv_type, error_type, true );

    // Set type of problem and matrices
    solver.SetProblemType ( problem_type );

    // Set type of spectral transformation
    solver.SetSTType ( st_type );
    solver.SetSTKSP ( st_ksp_type, st_ksp_max_iter, st_ksp_size, st_ksp_abs_tol, st_ksp_rel_tol );
    solver.SetSTPC ( st_pc_type, st_package );

    // Set desired eigenvalues
    solver.ComputeExtremalEigenValue ( eigenvalue_type );
    //	solver.ComputeTargetEigenValue(target, meas_type;
    //	solver.ComputeIntervalEigenValue(a, b);

    // Solve problem

    timer.reset ( );
    timer.start ( );
    if ( rank_ == 0 )
        std::cout << "Start EigenSolver " << solver_type << std::endl;

    int state = solver.Solve ( );
    timer.stop ( );

    // Get number of converged eigenvalues and iterations
    int n = solver.ConvEigValues ( );
    int iter = solver.Iter ( );

    // Get eigenvalues
    k_real_.resize ( n, 0. );
    k_imag_.resize ( n, 0. );
    std::vector<Scalar> res ( n, 0. );

    for ( int l = 0; l < n; l++ )
    {
        res[l] = solver.Res ( l );
        solver.GetEigenValue ( l, k_real_[l], k_imag_[l] );
    }

    if ( rank_ == 0 )
    {
        std::cout << "EigenSolver ended after " << iter << " iterations and " << timer.get_duration ( ) << " sec with state " << state << std::endl << std::endl;
        std::cout << "The following " << n << " eigenvalues with respective relative (SLEPC) residual norms could be computed: " << std::endl;
        std::cout.precision ( 3 );
        std::cout << std::scientific;
        for ( int l = 0; l < n; l++ )
        {
            std::cout << "[" << l << "] ";
            if ( l < 10 )
                std::cout << " ";
            std::cout << k_real_[l] << " + i" << k_imag_[l] << " | rel_res: " << res[l] << std::endl;
        }
    }

    //	if (nev < n)
    //		n = nev;

    // Get eigenvectors
    v_real_.resize ( n );
    v_imag_.resize ( n );
    for ( int l = 0; l < n; l++ )
    {
        v_real_[l] = new VectorType;
        v_imag_[l] = new VectorType;
        v_real_[l]->CloneFromWithoutContent ( *vec_ref_ );
        v_imag_[l]->CloneFromWithoutContent ( *vec_ref_ );

        solver.GetEigenVector ( l, *v_real_[l], *v_imag_[l] );
    }

    for ( int l = 0; l < n; l++ )
    {
        std::stringstream rname;
        rname << "real-" << refinement_level_ << "-" << l;
        std::stringstream iname;
        iname << "imag-" << refinement_level_ << "-" << l;

        visualize ( v_real_[l], rname.str ( ) );
        visualize ( v_imag_[l], iname.str ( ) );
    }

    // recompute residuals
    if ( rank_ == 0 )
    {
        std::cout << std::endl << "Compute residuals explicitly " << std::endl;
    }

    for ( int l = 0; l < n; ++l )
    {
        SetVectorsToZero ( );
        UpdatetVectorCouplings ( );

        A_->VectorMult ( *v_real_[l], &Arxr_ );
        A_->VectorMult ( *v_imag_[l], &Arxi_ );
        Ai_->VectorMult ( *v_real_[l], &Aixr_ );
        Ai_->VectorMult ( *v_imag_[l], &Aixi_ );
        B_->VectorMult ( *v_real_[l], &Bxr_ );
        B_->VectorMult ( *v_imag_[l], &Bxi_ );

        UpdatetVectorCouplings ( );

        res_r_.Axpy ( Arxr_, 1. );
        res_r_.Axpy ( Bxr_, -k_real_[l] );
        res_r_.Axpy ( Bxr_, k_imag_[l] );
        res_r_.Axpy ( Aixi_, -1. );

        res_i_.Axpy ( Arxi_, 1. );
        res_i_.Axpy ( Bxr_, -k_imag_[l] );
        res_i_.Axpy ( Bxi_, -k_real_[l] );
        res_i_.Axpy ( Aixr_, 1. );

        UpdatetVectorCouplings ( );
        double norm_res_r = res_r_.Norm2 ( );
        double norm_res_i = res_i_.Norm2 ( );

        double norm_v_r = Bxr_.Norm2 ( );
        double norm_v_i = Bxi_.Norm2 ( );

        double abs_res = std::sqrt ( norm_res_r * norm_res_r + norm_res_i * norm_res_i );
        double abs_norm = std::sqrt ( norm_v_r * norm_v_r + norm_v_i * norm_v_i );
        double k_norm = std::sqrt ( k_real_[l] * k_real_[l] + k_imag_[l] * k_imag_[l] );

        if ( rank_ == 0 )
        {
            std::cout << "[" << l << "] ";
            if ( l < 10 )
                std::cout << " ";
            std::cout << k_real_[l] << " + i" << k_imag_[l] << " | real_abs_res: " << norm_res_r << " , imag_abs_res: " << norm_res_i << " , abs_res: " << abs_res << " , rel_res: " << abs_res / k_norm;
            if ( abs_res / k_norm > 0.01 )
                std::cout << " ! " << std::endl;
            else
                std::cout << std::endl;
        }
    }
}

void PoissonEigenValue::SetVectorsToZero ( )
{
    Arxr_.Zeros ( );
    Aixi_.Zeros ( );
    Arxi_.Zeros ( );
    Aixr_.Zeros ( );
    Bxr_.Zeros ( );
    Bxi_.Zeros ( );

    res_r_.Zeros ( );
    res_i_.Zeros ( );
}

void PoissonEigenValue::UpdatetVectorCouplings ( )
{
    Arxr_.Update ( );
    Aixi_.Update ( );
    Arxi_.Update ( );
    Aixr_.Update ( );
    Bxr_.Update ( );
    Bxi_.Update ( );

    res_r_.Update ( );
    res_i_.Update ( );
}

void PoissonEigenValue::visualize ( VectorType* vec, std::string vec_name )
{

    ParallelCellVisualization<double> visu ( space_, 1, comm_, MASTER_RANK );

    // Generate filename.
    std::stringstream name;
    name << "out/" << vec_name;

    std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

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
    }

#ifdef USE_PETSC
    //	vec->Update();
#else
#    ifdef USE_HYPRE
    vec->Update ( );
#    else
    vec->UpdateCouplings ( );
#    endif
#endif

#ifdef USE_PETSC
    //	visu.visualize(EvalFeFunction<LAD>(space_, vec), "u");
#else
    visu.visualize ( EvalFeFunction<LAD>( space_, *vec ), "u" );
#endif

    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( name.str ( ) );
}

void PoissonEigenValue::adapt ( )
{
    // Refine mesh on master process. 1D parallel execution is not yet implemented.
    if ( DIMENSION == 1 )
    {
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( );
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
        if ( rank_ == MASTER_RANK )
        {
            const int final_ref_level = params_["Mesh"]["FinalRefLevel"].get<int>( );
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

        // Distribute mesh over all processes, and compute ghost cells
        if ( !is_done_ )
        {
            // Distribute the new mesh.
            MeshPtr local_mesh;

            // Distribute mesh over all processes, and compute ghost cells
            if ( num_partitions_ > 1 )
            {
#ifdef WITH_METIS
                MetisGraphPartitioner partitioner;
#else
                NaiveGraphPartitioner partitioner;
#endif
                local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &partitioner );
            }
            else
            {
                NaiveGraphPartitioner partitioner;
                local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_, &partitioner );
            }

            assert ( local_mesh != 0 );
            SharedVertexTable shared_verts;

            mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
        }
    }
}
