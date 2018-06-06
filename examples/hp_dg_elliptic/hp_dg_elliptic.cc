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

/// \author Staffan Ronnas

#include "hp_dg_elliptic.h"

#include <boost/bind.hpp>

static const char* PARAM_FILENAME = "hp_dg_elliptic_param.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class HpDgElliptic
{
  public:

    HpDgElliptic ( const std::string& param_filename, const std::string& path_mesh )
    : comm_ ( MPI_COMM_WORLD ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    path_mesh ( path_mesh ),
    space_ ( MPI_COMM_WORLD ),
    rhs_ ( 0 ), sol_ ( 0 ), matrix_ ( 0 ),
    is_done_ ( false ),
    refinement_level_ ( 0 )
    {
        MPI_Comm_rank ( comm_, &rank_ );
        MPI_Comm_size ( comm_, &num_partitions_ );
    }

    ~HpDgElliptic ( )
    {
        delete sol_;
        delete rhs_;
        delete matrix_;
    }

    // Main algorithm

    void run ( )
    {
        // Construct / read in the initial mesh.
        bool is_master = ( rank_ == MASTER_RANK );
        if ( is_master ) std::cout << "Begin building mesh ...";
        build_initial_mesh ( );
        if ( is_master ) std::cout << "... done." << std::endl;
        // Main adaptation loop.
        while ( !is_done_ )
        {
            // Initialize space and linear algebra.
            if ( is_master ) std::cout << "Preparing system ..." << std::endl;
            prepare_system ( );
            if ( is_master ) std::cout << "... done." << std::endl;
            // Compute the stiffness matrix and right-hand side.
            if ( is_master ) std::cout << "Assembling system..." << std::endl;
            assemble_system ( );
            if ( is_master ) std::cout << "... done." << std::endl;
            // Solve the linear system.
            if ( is_master ) std::cout << "Solving System ..." << std::endl;
            solve_system ( );
            if ( is_master ) std::cout << "... done." << std::endl;

            // Compute the error to the exact solution.
            if ( is_master ) std::cout << "Begin Postprocessing ..." << std::endl;
            compute_error ( );
            // Visualize the solution and the errors.
            visualize ( );
            if ( is_master ) std::cout << "... done." << std::endl;
            // Modify the space through refinement. Set is_done_ = true when finished.
            if ( is_master ) std::cout << "Adapting mesh ...";
            adapt ( );
            if ( is_master ) std::cout << "... done." << std::endl;
            if ( is_master ) std::cout << " ==================================== \n";
        }
    }

  private:
    // Member functions
    // Read and distribute mesh.
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

    void debug_mesh ( );

    // member variables

    // MPI communicator.
    const MPI_Comm comm_;
    // Local process rank and number of processes.
    int rank_, num_partitions_;

    // Parameter data read in from file.
    PropertyTree params_;

    // Local mesh and mesh on master process.
    MeshPtr mesh_, master_mesh_;
    std::string path_mesh;

    // Solution space.
    VectorSpace<double> space_;

    // Linear algebra couplings helper object.
    Couplings<double> la_couplings_; //TODO: rename to couplings
    // Vectors for solution and load vector.
    CoupledVector<Scalar>* rhs_, *sol_;
    // System matrix.
    CoupledMatrix<Scalar>* matrix_;

    // Global assembler.
#ifdef USE_DG
    DGGlobalAssembler<double> global_asm_;
#else
    StandardGlobalAssembler<double> global_asm_;
#endif

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
}; // end class HpDgElliptic

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
        std::ofstream info_log ( "hp_dg_elliptic_info_log" );
        LogKeeper::get_log ( "info" ).set_target ( &info_log );
        std::ofstream debug_log ( "hp_dg_elliptic_debug_log" );
        LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

        // Create application object and run it
        HpDgElliptic app ( param_filename, path_mesh );
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

//////////////// HpDgElliptic implementation //////////////

void HpDgElliptic::build_initial_mesh ( )
{

    // Read in the mesh on the master process.
    if ( rank_ == MASTER_RANK )
    {
        std::string mesh_name =
                params_["Mesh"]["Filename"].get<std::string>( );
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
        const int initial_ref_lvl = params_["Mesh"]["InitialRefLevel"].get<int>( );
        for ( int r = 0; r < initial_ref_lvl; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
            ++refinement_level_;
        }
    }

    MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );

    // Distribute mesh over all processes, and compute ghost cells
    MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK,
                                                    MPI_COMM_WORLD );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, MPI_COMM_WORLD, shared_verts );

    // InterfaceList if_list = InterfaceList::create(master_mesh_);
    // std::cout << "If-list = \n" << if_list << "\n";
}

void HpDgElliptic::prepare_system ( )
{

    int degree = params_["Mesh"]["Degree"].get<int>( );
    // Initialize the VectorSpace object.
    std::vector< int> deg_vec ( 1 );
#ifdef USE_DG
    std::vector<bool> is_cg ( 1, false );
#else
    std::vector<bool> is_cg ( 1, true );
#endif
    deg_vec[0] = degree;

    space_.Init ( deg_vec, *mesh_, is_cg );

    std::cout << "Num cells (" << rank_ << ") = " << mesh_->num_entities ( DIMENSION ) << "\n";
    std::cout << "Num dof (" << rank_ << ") = " << space_.dof ( ).get_nb_dofs ( ) << "\n";

    // Setup couplings object.
    la_couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    // extract input for LaCouplings
    int nb_procs = 0;
    int info = MPI_Comm_size ( comm_, &nb_procs );
    assert ( info == MPI_SUCCESS );
    assert ( nb_procs > 0 );

    // determine global offsets
    la_couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                        sparsity.off_diagonal_cols );

    // Setup linear algebra objects.
    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    matrix_->Init ( comm_, la_couplings_ );
    CoupledVectorFactory<Scalar> CoupVecFact;
    rhs_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    sol_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    rhs_->Init ( comm_, la_couplings_ );
    sol_->Init ( comm_, la_couplings_ );

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

    // Compute Dirichlet BC dofs and values using known exact solution.
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

#ifdef STRONG_BC

#    ifdef SING_SOL
    DirichletExact dir_bc;
#    endif

#    ifdef TRIG_SOL
    DirichletZero dir_bc;
#    endif

    compute_dirichlet_dofs_and_values ( dir_bc, space_, 0, dirichlet_dofs_,
                                        dirichlet_values_ );
#endif
}

void HpDgElliptic::assemble_system ( )
{
    // Assemble matrix and right-hand-side vector.
    DGCellAssembler local_cell_asm;

    matrix_->Zeros ( );
    rhs_->Zeros ( );

    global_asm_.should_reset_assembly_target ( true );
    global_asm_.assemble_matrix ( space_, local_cell_asm, *matrix_ );
    global_asm_.assemble_vector ( space_, local_cell_asm, *rhs_ );

#ifdef USE_DG
    global_asm_.should_reset_assembly_target ( false );

    const double theta = params_["DGMethod"]["Theta"].get<double>( );
    const double gamma = params_["DGMethod"]["Gamma"].get<double>( );
    DGJumpMatrixAssembler local_jump_asm ( theta, gamma );
    global_asm_.assemble_interface_matrix ( space_, local_jump_asm, *matrix_ );
    //    matrix_->diagonal().WriteFile("matrix_after.mat");

#    ifndef STRONG_BC
    DGBoundaryVectorAssembler local_bdy_asm ( theta, gamma );
    global_asm_.assemble_interface_vector ( space_, local_bdy_asm, *rhs_ );
#    endif

    global_asm_.should_reset_assembly_target ( true );
#endif

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        //matrix_->ZeroRows(vec2ptr(dirichlet_dofs_), dirichlet_dofs_.size(), 1.0);
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );

    }
    // Communicate modified values outside if-block
    rhs_->UpdateCouplings ( );
    sol_->UpdateCouplings ( );
}

void HpDgElliptic::solve_system ( )
{

    LinearSolver<LAD>* solver_;
#ifndef WITH_UMFPACK
    std::cout << "Solving with GMRES\n";
    LinearSolverFactory<LAD> SolFact;
    solver_ = SolFact.Get (
                            params_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( params_["LinearSolver"] );

#    ifdef WITH_ILUPP
    std::cout << "Preconditioning with ILU++\n";
    PreconditionerIlupp<LAD> precond;
    precond.InitParameter ( 1,
                            1010,
                            150,
                            0.8,
                            5.,
                            0.05 );
#    else
    std::cout << "Preconditioning with Gauss-Seidel\n";
    PreconditionerBlockJacobiStand<LAD> precond;
    precond.Init_GaussSeidel ( );
#    endif
    precond.SetupOperator ( *matrix_ );

    solver_->SetupOperator ( *matrix_ );
    solver_->SetupPreconditioner ( precond );
#else
    std::cout << "Solving with UMFPACK\n";
    solver_ = new UmfpackSolver<LAD>( );
#endif
    solver_->SetupOperator ( *matrix_ );
    const LinearSolverState state = solver_->Solve ( *rhs_, sol_ );

    delete solver_;

    if ( state == kSolverSuccess )
    {
        std::cout << "Linear solver finished successfully.\n";
    }
    else
    {
        std::cout << "Linear solver failed to converge.\n";
        throw "linear solver failed";
    }

    // Communicate new solution values.
    sol_->UpdateCouplings ( );
}

void HpDgElliptic::compute_error ( )
{
    L2_err_.clear ( );
    H1_err_.clear ( );

    // Compute square of the L2 error on each element, putting the
    // values into L2_err_.

    //global_asm_.set_quadrature_selection_function(qfun_err);

    L2ErrorIntegrator<ExactSol> L2_int ( *sol_ );
    global_asm_.assemble_scalar ( space_, L2_int, L2_err_ );
    // Create attribute with L2 error for output.
    AttributePtr L2_err_attr ( new DoubleAttribute ( L2_err_ ) );
    mesh_->add_attribute ( "L2 error", DIMENSION, L2_err_attr );

    double total_L2_err = std::accumulate ( L2_err_.begin ( ), L2_err_.end ( ), 0. );
    double global_L2_err = 0.;
    MPI_Reduce ( &total_L2_err, &global_L2_err, 1, MPI_DOUBLE, MPI_SUM, 0,
                 comm_ );
    LOG_INFO ( "error", "Local L2 error on partition " << rank_ << " = "
               << std::sqrt ( total_L2_err ) );

    // Compute square of the H1 error on each element, putting the
    // values into H1_err_.

    H1ErrorIntegrator<ExactSol> H1_int ( *sol_ );
    global_asm_.assemble_scalar ( space_, H1_int, H1_err_ );

    // Create attribute with H1 error for output.
    AttributePtr H1_err_attr ( new DoubleAttribute ( H1_err_ ) );
    mesh_->add_attribute ( "H1 error", DIMENSION, H1_err_attr );
    double total_H1_err = std::accumulate ( H1_err_.begin ( ), H1_err_.end ( ), 0. );
    double global_H1_err = 0.;
    MPI_Reduce ( &total_H1_err, &global_H1_err, 1, MPI_DOUBLE, MPI_SUM, 0,
                 comm_ );
    LOG_INFO ( "error", "Local H1 error on partition " << rank_ << " = "
               << std::sqrt ( total_H1_err ) );

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "Global L2 error = " << std::sqrt ( global_L2_err ) << "\n"
                << "Global H1 error = " << std::sqrt ( global_H1_err ) << "\n";
    }
}

struct EvalF
{

    double operator() ( const Vec<DIMENSION, double>& pt ) const
    {
        return sol_.eval_f ( pt );
    }

    ExactSol sol_;
};

void HpDgElliptic::visualize ( )
{

    ParallelCellVisualization<double> visu ( space_, params_["Mesh"]["Degree"].get<int>( ), comm_, MASTER_RANK );
    visu.visualize ( EvalFeFunction<LAD>( space_, *sol_ ), "u" );
    ExactSol s;
    visu.visualize ( EvalAnalyticFunction<ExactSol>( space_, s ), "u_exact" );

    std::vector<double> L2_err_sqr ( L2_err_.size ( ) ), H1_err_sqr ( H1_err_.size ( ) );

    std::transform ( L2_err_.begin ( ), L2_err_.end ( ), L2_err_sqr.begin ( ), &sqrtf );
    std::transform ( H1_err_.begin ( ), H1_err_.end ( ), H1_err_sqr.begin ( ), &sqrtf );

    visu.visualize ( VisualizeCellValues ( L2_err_sqr ), "L2 err" );
    visu.visualize ( VisualizeCellValues ( H1_err_sqr ), "H1 err" );

    // Generate filename.
    std::stringstream name;
    name << "solution" << refinement_level_;
    visu.write ( name.str ( ) );

}

void HpDgElliptic::adapt ( )
{
    // Refine mesh on master process.
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

    if ( !is_done_ )
    {
        // Distribute the new mesh.
        MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK,
                                                        MPI_COMM_WORLD );
        assert ( local_mesh != 0 );
        SharedVertexTable shared_verts;
        mesh_ = compute_ghost_cells ( *local_mesh, MPI_COMM_WORLD, shared_verts );
    }

    std::ostringstream sstr;
    sstr << "mesh" << refinement_level_ << ".vtu";

    VtkWriter writer;
    writer.write ( sstr.str ( ).c_str ( ), *mesh_ );

}

void HpDgElliptic::debug_mesh ( )
{
    std::cerr << "Vertices:\n";

    for ( EntityIterator it = mesh_->begin ( 0 ), end_it = mesh_->end ( 0 );
          it != end_it; ++it )
    {
        std::vector<double> vertex_coords;
        it->get_coordinates ( vertex_coords );
        std::cerr << "\tind = " << it->index ( ) << ", id = " << it->id ( )
                << ", coords = "
                << string_from_range ( vertex_coords.begin ( ),
                                       vertex_coords.end ( ) ) << "\n";
    }
    std::cerr << "\n\n";
    std::cerr << "Facets:\n";

    for ( EntityIterator it = mesh_->begin ( DIMENSION - 1 ), end_it = mesh_->end ( DIMENSION - 1 );
          it != end_it; ++it )
    {
        std::cerr << "\tind = " << it->index ( ) << ", id = " << it->id ( )
                << ", verts = ";
        for ( IncidentEntityIterator iit = it->begin_incident ( 0 ), end_iit = it->end_incident ( 0 ); iit != end_iit; ++iit )
        {
            std::cerr << iit->id ( ) << ", ";
        }
        std::cerr << "\n";
    }

    std::cerr << "\n\n";
    std::cerr << "Cells:\n";
    for ( EntityIterator it = mesh_->begin ( DIMENSION ), end_it = mesh_->end ( DIMENSION );
          it != end_it; ++it )
    {
        std::cerr << "\tind = " << it->index ( ) << ", id = " << it->id ( )
                << ", verts = ";
        for ( IncidentEntityIterator iit = it->begin_incident ( 0 ), end_iit = it->end_incident ( 0 ); iit != end_iit; ++iit )
        {
            std::cerr << iit->id ( ) << ", ";
        }
        std::cerr << "facets = ";
        for ( IncidentEntityIterator iit = it->begin_incident ( DIMENSION - 1 ), end_iit = it->end_incident ( DIMENSION - 1 );
              iit != end_iit; ++iit )
        {
            std::cerr << iit->id ( ) << ", ";
        }
        std::cerr << "\n";
    }
    std::cerr << "\n\n";
}
