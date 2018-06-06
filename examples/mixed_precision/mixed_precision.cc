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

/// \author Martin Wlotzka

#include "mixed_precision.h"

static const char* PARAM_FILENAME = "mixed_precision.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class MixedPrecisionExample
{
  public:

    MixedPrecisionExample ( const std::string& param_filename, const std::string& path_mesh )
    : is_done_ ( false ),
    refinement_level_ ( 0 ),
    comm_ ( MPI_COMM_WORLD ),
    rank_ ( -1 ),
    num_partitions_ ( -1 ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    Ad_ ( 0 ), xd_ ( 0 ), bd_ ( 0 ), resd_ ( 0 ), cord_ ( 0 ),
    Af_ ( 0 ), resf_ ( 0 ), corf_ ( 0 ),
    path_mesh ( path_mesh )
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
            solve_system ( );
            // Compute the error to the exact solution.
            visualize ( );
            adapt ( );
        }
    }

    ~MixedPrecisionExample ( )
    {
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
    void compute_error ( );
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
    MeshPtr mesh_, master_mesh_;

    // Solution space.
    VectorSpace<double> spd_;
    VectorSpace<float> spf_;

    // Linear algebra couplings helper object.
    Couplings<double> cpd_;
    Couplings<float> cpf_;

    // Vectors for solution and load vector.
    CoupledVector<double> *bd_, *xd_, *resd_, *cord_;
    CoupledVector<float> *resf_, *corf_;
    // System matrix.
    CoupledMatrix<double> *Ad_;
    CoupledMatrix<float> *Af_;

    // Global assembler.
    StandardGlobalAssembler<double> gad_;
    StandardGlobalAssembler<float> gaf_;

    // Flag for stopping adaptive loop.
    bool is_done_;
    // Current refinement level.
    int refinement_level_;

    int use_float_assembly_;

    // Dof id:s for Dirichlet boundary conditions.
    std::vector<int> ddd_;
    std::vector<int> ddf_;
    // Dof values for Dirichlet boundary conditions.
    std::vector<double> dvd_;
    std::vector<float> dvf_;
}; // end class MixedPrecisionExample

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

    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    try
    {
        // Create log files for INFO and DEBUG output
        std::ofstream info_log ( "info.log" );
        if ( !rank ) LogKeeper::get_log ( "info" ).set_target ( &info_log );
        std::ofstream debug_log ( "debug.log" );
        if ( !rank ) LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

        INFO = ( rank == 0 );

        // Create application object and run it
        MixedPrecisionExample app ( param_filename, path_mesh );
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

//////////////// MixedPrecisionExample implementation //////////////

void MixedPrecisionExample::build_initial_mesh ( )
{
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

        // Distribute mesh over all processes, and compute ghost cells
        MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_ );
        assert ( local_mesh != 0 );
        SharedVertexTable shared_verts;
        mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );

        // Write out mesh of initial refinement level
        PVtkWriter writer ( comm_ );
        std::string output_file = std::string ( "mpir_mesh.pvtu" );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );
    }
}

void MixedPrecisionExample::prepare_system ( )
{
    // Assign degrees to each element.
    const int fe_degree = params_["Mesh"]["FeDegree"].get<int>( 1 );
    std::vector< int > degrees ( 1, fe_degree );

    // Initialize the VectorSpace object.
    spd_.Init ( degrees, *mesh_ );
    spf_.Init ( degrees, *mesh_ );

    // Setup couplings object.
    cpd_.Init ( comm_, spd_.dof ( ) );
    cpf_.Init ( comm_, spf_.dof ( ) );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    gad_.compute_sparsity_structure ( spd_, sparsity );

    cpd_.InitializeCouplings ( sparsity.off_diagonal_rows,
                               sparsity.off_diagonal_cols );

    cpf_.InitializeCouplings ( sparsity.off_diagonal_rows,
                               sparsity.off_diagonal_cols );

    PLATFORM mplat, vplat;
    IMPLEMENTATION mimpl, vimpl;
    MATRIX_FORMAT mform;

    std::string s_mplat = params_["LinearAlgebra"]["CoupledMatrix"]["Platform"].get<std::string>( "CPU" );
    std::string s_mimpl = params_["LinearAlgebra"]["CoupledMatrix"]["Implementation"].get<std::string>( "NAIVE" );
    std::string s_mform = params_["LinearAlgebra"]["CoupledMatrix"]["Format"].get<std::string>( "CSR" );

    std::string s_vplat = params_["LinearAlgebra"]["CoupledVector"]["Platform"].get<std::string>( "CPU" );
    std::string s_vimpl = params_["LinearAlgebra"]["CoupledVector"]["Implementation"].get<std::string>( "NAIVE" );

    if ( s_mplat == "CPU" ) mplat = CPU;
    else if ( s_mplat == "GPU" ) mplat = GPU;
    else if ( s_mplat == "OPENCL" ) mplat = OPENCL;
    else quit_program ( );

    if ( s_mimpl == "NAIVE" ) mimpl = NAIVE;
    else if ( s_mimpl == "OPENMP" ) mimpl = OPENMP;
    else if ( s_mimpl == "MKL" ) mimpl = MKL;
    else if ( s_mimpl == "SCALAR" ) mimpl = SCALAR;
    else if ( s_mimpl == "SCALAR_TEX" ) mimpl = SCALAR_TEX;
    else quit_program ( );

    if ( s_mform == "CSR" ) mform = CSR;
    else if ( s_mform == "COO" ) mform = COO;
    else if ( s_mform == "ELL" ) mform = ELL;
    else if ( s_mform == "DENSE" ) mform = DENSE;
    else quit_program ( );

    if ( s_vplat == "CPU" ) vplat = CPU;
    else if ( s_vplat == "GPU" ) vplat = GPU;
    else if ( s_vplat == "OPENCL" ) vplat = OPENCL;
    else quit_program ( );

    if ( s_vimpl == "NAIVE" ) vimpl = NAIVE;
    else if ( s_vimpl == "OPENMP" ) vimpl = OPENMP;
    else if ( s_vimpl == "MKL" ) vimpl = MKL;
    else if ( s_vimpl == "BLAS" ) vimpl = BLAS;
    else quit_program ( );

    if ( Ad_ != 0 ) delete Ad_;
    if ( Af_ != 0 ) delete Af_;
    Ad_ = new CoupledMatrix<double>( );
    Af_ = new CoupledMatrix<float>( );

    Ad_->Init ( comm_, cpd_, mplat, mimpl, mform );
    Af_->Init ( comm_, cpf_, mplat, mimpl, mform );

    if ( xd_ != 0 ) delete xd_;
    if ( bd_ != 0 ) delete bd_;
    if ( resd_ != 0 ) delete resd_;
    if ( cord_ != 0 ) delete cord_;
    if ( resf_ != 0 ) delete resf_;
    if ( corf_ != 0 ) delete corf_;

    bd_ = new CoupledVector<double>( );
    xd_ = new CoupledVector<double>( );
    resd_ = new CoupledVector<double>( );
    cord_ = new CoupledVector<double>( );
    resf_ = new CoupledVector<float>( );
    corf_ = new CoupledVector<float>( );

    bd_->Init ( comm_, cpd_, vplat, vimpl );
    xd_->Init ( comm_, cpd_, vplat, vimpl );
    resd_->Init ( comm_, cpf_, vplat, vimpl );
    cord_->Init ( comm_, cpf_, vplat, vimpl );
    resf_->Init ( comm_, cpf_, vplat, vimpl );
    corf_->Init ( comm_, cpf_, vplat, vimpl );

    Ad_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                         vec2ptr ( sparsity.diagonal_cols ),
                         sparsity.diagonal_rows.size ( ),
                         vec2ptr ( sparsity.off_diagonal_rows ),
                         vec2ptr ( sparsity.off_diagonal_cols ),
                         sparsity.off_diagonal_rows.size ( ) );

    Af_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                         vec2ptr ( sparsity.diagonal_cols ),
                         sparsity.diagonal_rows.size ( ),
                         vec2ptr ( sparsity.off_diagonal_rows ),
                         vec2ptr ( sparsity.off_diagonal_cols ),
                         sparsity.off_diagonal_rows.size ( ) );

    bd_->InitStructure ( );
    xd_->InitStructure ( );
    resd_->InitStructure ( );
    cord_->InitStructure ( );
    resf_->InitStructure ( );
    corf_->InitStructure ( );

    // Zero all linear algebra objects.
    Ad_->Zeros ( );
    Af_->Zeros ( );
    bd_->Zeros ( );
    xd_->Zeros ( );

    // Compute Dirichlet BC dofs and values using known exact solution.
    ddd_.clear ( );
    dvd_.clear ( );
    dvf_.clear ( );

    DirichletZeroD zerod;
    DirichletZeroS zerof;

    // The function compute_dirichlet_dofs_and_values des not yet work for 1D.
    if ( DIMENSION == 1 )
    {
        dvd_.resize ( 2, 0.0 );

        // Loop over all cells.
        for ( EntityIterator facet_it = mesh_->begin ( DIMENSION - 1 ), facet_end = mesh_->end ( DIMENSION - 1 );
              facet_it != facet_end; ++facet_it )
        {

            // Returns the number of neighbors for each cell, to check if it is on the facet.
            const EntityCount num_cell_neighbors = facet_it->num_incident_entities ( DIMENSION );

            // If it lies on the facet, the corresponding DOF is a Dirichlet DOF and is added to dirichlet_dofs_.
            if ( num_cell_neighbors == 1 )
            {
                std::vector<int> dof_number_d;
                spd_.dof ( ).get_dofs_on_subentity ( 0, facet_it->begin_incident ( DIMENSION )->index ( ), 0, facet_it->index ( ), dof_number_d );
                ddd_.push_back ( dof_number_d[0] );

                std::vector<int> dof_number_f;
                spf_.dof ( ).get_dofs_on_subentity ( 0, facet_it->begin_incident ( DIMENSION )->index ( ), 0, facet_it->index ( ), dof_number_f );
                ddf_.push_back ( dof_number_f[0] );
            }
        }
    }
    else
    {
        compute_dirichlet_dofs_and_values ( zerod, spd_, 0, ddd_, dvd_ );
        compute_dirichlet_dofs_and_values ( zerof, spf_, 0, ddf_, dvf_ );
    }

    use_float_assembly_ = params_["MPIR"]["UseFloatAssembly"].get<int>( 0 );
}

void MixedPrecisionExample::assemble_system ( )
{
    // Assemble matrix and right-hand-side vector.
    LocalPoissonAssemblerD lad;
    gad_.assemble_matrix ( spd_, lad, *Ad_ );
    gad_.assemble_vector ( spd_, lad, *bd_ );

    if ( use_float_assembly_ )
    {
        LocalPoissonAssemblerF laf;
        gaf_.assemble_matrix ( spf_, laf, *Af_ );
    }
    else
    {
        Af_->CastFrom ( *Ad_ );
    }

    if ( !ddd_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        Ad_->diagonalize_rows ( vec2ptr ( ddd_ ), ddd_.size ( ), 1.0 );
        bd_->SetValues ( vec2ptr ( ddd_ ), ddd_.size ( ),
                         vec2ptr ( dvd_ ) );
        xd_->SetValues ( vec2ptr ( ddd_ ), ddd_.size ( ),
                         vec2ptr ( dvd_ ) );
    }

    if ( !ddf_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        Af_->diagonalize_rows ( vec2ptr ( ddf_ ), ddf_.size ( ), 1.0 );
    }
}

void MixedPrecisionExample::solve_system ( )
{

    CG<LAD> ecd;
    CG<LAS> ecf;

    Jacobi<LAD> pre_basync;
    Jacobi<LAD> post_basync;

    MPIR<LAD, LAS> mpir;

    mpir.SetHighPrecMatrix ( *Ad_ );
    mpir.SetLowPrecMatrix ( *Af_ );
    mpir.SetHighPrecResVector ( *resd_ );
    mpir.SetHighPrecCorVector ( *cord_ );
    mpir.SetLowPrecResVector ( *resf_ );
    mpir.SetLowPrecCorVector ( *corf_ );

    int max_ecd = params_["MPIR"]["MaxHighPrecErrCorLoops"].get<int>( 10 );
    int max_ecf = params_["MPIR"]["MaxLowPrecErrCorLoops"].get<int>( 10 );
    ScalarD abs_tol_mpir = params_["MPIR"]["AbsTol"].get<ScalarD>( 1.0e-12 );
    ScalarD rel_tol_mpir = params_["MPIR"]["RelTol"].get<ScalarD>( 1.0e-12 );
    ScalarD div_tol_mpir = params_["MPIR"]["DivTol"].get<ScalarD>( 1.0e6 );

    mpir.SetMpirParameters ( max_ecd, max_ecf, abs_tol_mpir, rel_tol_mpir, div_tol_mpir );

    int use_pre_jacd = params_["MPIR"]["HighPrecPreSmoother"]["Use"].get<int>( 0 );
    int max_iter_jacd = params_["MPIR"]["HighPrecPreSmoother"]["MaxIter"].get<int>( 3 );
    int n_iter = params_["MPIR"]["HighPrecPreSmoother"]["InnerIter"].get<int>( 3 );
    ScalarD rel_tol_jacd = params_["MPIR"]["HighPrecPreSmoother"]["RelTol"].get<ScalarD>( 1.0e-12 );

    pre_basync.Prepare ( *Ad_, *bd_ );
    pre_basync.SetSolveMode ( 3 );
    pre_basync.SetNumIter ( n_iter );
    pre_basync.SetInnerIter ( 1 );
    pre_basync.SetDampingParameter ( 0.5 );
    mpir.SetHighPrecPreSmoother ( pre_basync );
    mpir.SetHighPrecPreSmootherParameters ( max_iter_jacd, abs_tol_mpir, rel_tol_jacd, div_tol_mpir );
    mpir.UseHighPrecPreSmoother ( use_pre_jacd > 0 );

    int use_post_jacd = params_["MPIR"]["HighPrecPostSmoother"]["Use"].get<int>( 0 );
    max_iter_jacd = params_["MPIR"]["HighPrecPostSmoother"]["MaxIter"].get<int>( 3 );
    n_iter = params_["MPIR"]["HighPrecPostSmoother"]["InnerIter"].get<int>( 3 );
    rel_tol_jacd = params_["MPIR"]["HighPrecPostSmoother"]["RelTol"].get<ScalarD>( 1.0e-12 );

    post_basync.Prepare ( *Ad_, *bd_ );
    post_basync.SetSolveMode ( 3 );
    post_basync.SetNumIter ( n_iter );
    post_basync.SetInnerIter ( 1 );
    post_basync.SetDampingParameter ( 0.5 );
    mpir.SetHighPrecPostSmoother ( post_basync );
    mpir.SetHighPrecPostSmootherParameters ( max_iter_jacd, abs_tol_mpir, rel_tol_jacd, div_tol_mpir );
    mpir.UseHighPrecPostSmoother ( use_post_jacd > 0 );

    int max_iter_ecd = params_["MPIR"]["HighPrecErrCorSolver"]["MaxIter"].get<int>( 1000 );
    ScalarD rel_tol_ecd = params_["MPIR"]["HighPrecErrCorSolver"]["RelTol"].get<ScalarD>( 1.0e-12 );

    mpir.SetHighPrecErrCorSolver ( ecd );
    mpir.SetHighPrecErrCorSolverParameters ( max_iter_ecd, abs_tol_mpir, rel_tol_ecd, div_tol_mpir );

    int max_iter_ecf = params_["MPIR"]["LowPrecErrCorSolver"]["MaxIter"].get<int>( 1000 );
    ScalarS abs_tol_ecf = params_["MPIR"]["LowPrecErrCorSolver"]["AbsTol"].get<ScalarS>( 1.0e-6 );
    ScalarS rel_tol_ecf = params_["MPIR"]["LowPrecErrCorSolver"]["RelTol"].get<ScalarS>( 1.0e-6 );
    ScalarS div_tol_ecf = params_["MPIR"]["LowPrecErrCorSolver"]["DivTol"].get<ScalarS>( 1.0e6 );

    mpir.SetLowPrecErrCorSolver ( ecf );
    mpir.SetLowPrecErrCorSolverParameters ( max_iter_ecf, abs_tol_ecf, rel_tol_ecf, div_tol_ecf );

    mpir.Solve ( *bd_, xd_ );
}

void MixedPrecisionExample::visualize ( )
{
    // Setup visualization object.
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( spd_, num_intervals, comm_, MASTER_RANK );

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

    visu.visualize ( EvalFeFunction<LAD>( spd_, *( xd_ ) ), "u" );

    // visualize error measures
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.write ( name.str ( ) );
}

void MixedPrecisionExample::adapt ( )
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
            MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK,
                                                            comm_ );
            assert ( local_mesh != 0 );
            SharedVertexTable shared_verts;
            mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
        }
    }
}
