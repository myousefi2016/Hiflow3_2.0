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

/// \author Staffan Ronnas<br>Julian Kraemer

#include "direct_inverse_tutorial.h"

static const char* PARAM_FILENAME = "direct_inverse_tutorial.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class DirectInverseTutorial
{
  public:

    DirectInverseTutorial ( const std::string& param_filename, const std::string& path_mesh )
    : path_mesh ( path_mesh ),
    comm_ ( MPI_COMM_WORLD ),
    rank_ ( -1 ),
    num_partitions_ ( -1 ),
    params_ ( param_filename, MASTER_RANK, MPI_COMM_WORLD ),
    rhs_ ( 0 ), sol_ ( 0 ), row_ ( 0 ), coupled_sol_ ( 0 ), matrix_ ( 0 ), mass_matrix_ ( 0 ),
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
            solve_system ( );
            // Visualize the solution.
            visualize ( );
            // Modify the space through refinement. Set is_done_ = true when finished.
            adapt ( );
        }
    }

    ~DirectInverseTutorial ( )
    {
        delete matrix_;
        delete sol_;
        delete rhs_;
        delete row_;
        delete coupled_sol_;
        delete mass_matrix_;
    }

  private:
    // Member functions

    // Read and distribute mesh.
    std::string path_mesh;
    void build_initial_mesh ( );
    // Setup space, linear algebra
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
    MeshPtr mesh_, master_mesh_;
    // Solution space.
    VectorSpace<double> space_;

    // Solution vector to visualize.
    std::vector<double> solution_vec;

    // Vectors which contains the rows and columns of the entries of the mass
    // matrix which or not zero.
    std::vector<int> row_entries_;
    std::vector<int> col_entries_;

    // Linear algebra couplings helper object.
    LaCouplings la_couplings_;

    // Vectors for solution and load vector.
    CoupledVector<Scalar>* rhs_, *sol_, *row_, *coupled_sol_;

    // System matrices.
    CoupledMatrix<Scalar>* matrix_;
    CoupledMatrix<Scalar>* mass_matrix_;

    // Trajectory Vectors
    std::vector<double> u_trajectory;

    // Global assembler.
    StandardGlobalAssembler<double> global_asm_;

    // Flag for stopping adaptive loop.
    bool is_done_;
    // Current refinement level.
    int refinement_level_;
}; // end class DirectInverseTutorial

// Program entry point

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    int num_partitions;
    MPI_Comm_size ( MPI_COMM_WORLD, &num_partitions );

    if ( num_partitions > 1 )
    {
        std::cout << "NO parallel execution supported! Works only sequentially!" << std::endl;
        exit ( -1 );
    }

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
        std::ofstream info_log ( "DirectInverseTutorial_info_log" );
        LogKeeper::create_log ( "info", &info_log );
        std::ofstream debug_log ( "DirectInverseTutorial_debug_log" );
        LogKeeper::create_log ( "debug", &debug_log );

        // Create application object and run it
        DirectInverseTutorial app ( param_filename, path_mesh );
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

//////////////// DirectInverseTutorial implementation //////////////

void DirectInverseTutorial::build_initial_mesh ( )
{
    // Read in the mesh on the master process.
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
                                                    comm_ );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );

    // Write out mesh of initial refinement level
    PVtkWriter writer ( comm_ );
    std::string output_file = std::string ( "DirectInverseTutorial_mesh.pvtu" );
    writer.add_all_attributes ( *mesh_, true );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

void DirectInverseTutorial::prepare_system ( )
{
    // Assign degrees to each element.
    const int fe_degree = params_["Mesh"]["FeDegree"].get<int>( );
    std::vector< int > degrees ( 1, fe_degree );

    // Initialize the VectorSpace object.
    space_.Init ( degrees, *mesh_ );

    // Setup couplings object.
    la_couplings_.Init ( comm_ );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    // extract input for LaCouplings
    int nb_procs = 0;
    int info = MPI_Comm_size ( comm_, &nb_procs );
    assert ( info == MPI_SUCCESS );
    assert ( nb_procs > 0 );

    // determine global offsets
    std::vector<int> global_offsets ( nb_procs + 1, 0 );
    for ( int id = 0; id < nb_procs; id++ )
    {
        global_offsets[id + 1] = global_offsets[id] +
                space_.dof ( ).ndofs_on_sd ( id );
    }
    assert ( global_offsets.back ( ) == space_.dof ( ).ndofs_global ( ) );

    std::vector<int> sorted_offdiag_cols = sparsity.off_diagonal_cols;
    std::sort ( sorted_offdiag_cols.begin ( ), sorted_offdiag_cols.end ( ) );

    // determine offdiagonal offsets
    std::vector<int> offdiag_offsets ( nb_procs + 1, 0 );
    std::vector<int> nb_dofs_per_proc ( nb_procs, 0 );
    int ind = 0;
    for ( int i = 0; i < static_cast < int > ( sparsity.off_diagonal_cols.size ( ) ); i++ )
    {
        ind = sorted_offdiag_cols[i];
        nb_dofs_per_proc[space_.dof ( ).owner_of_dof ( ind )] += 1;
    }
    assert ( std::accumulate ( nb_dofs_per_proc.begin ( ),
                               nb_dofs_per_proc.end ( ), 0 )
             == static_cast < int > ( sparsity.off_diagonal_cols.size ( ) ) );

    for ( int id = 0; id < nb_procs; id++ )
    {
        offdiag_offsets[id + 1] = offdiag_offsets[id] + nb_dofs_per_proc[id];
    }
    assert ( offdiag_offsets[nb_procs] == static_cast < int > ( sparsity.off_diagonal_cols.size ( ) ) );

    la_couplings_.InitializeCouplings ( global_offsets,
                                        sorted_offdiag_cols,
                                        offdiag_offsets );

    // Setup linear algebra objects.
    delete matrix_;
    delete rhs_;
    delete sol_;
    delete mass_matrix_;
    delete row_;
    delete coupled_sol_;

    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    matrix_->Init ( comm_, la_couplings_ );
    mass_matrix_ = CoupMaFact.Get (
                                    params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    mass_matrix_->Init ( comm_, la_couplings_ );

    CoupledVectorFactory<Scalar> CoupVecFact;
    rhs_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    sol_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    row_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    coupled_sol_ = CoupVecFact.Get (
                                     params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    rhs_->Init ( comm_, la_couplings_ );
    sol_->Init ( comm_, la_couplings_ );
    row_->Init ( comm_, la_couplings_ );
    coupled_sol_->Init ( comm_, la_couplings_ );

    // Initialize structure of LA objects.
    matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                             vec2ptr ( sparsity.diagonal_cols ),
                             sparsity.diagonal_rows.size ( ),
                             vec2ptr ( sparsity.off_diagonal_rows ),
                             vec2ptr ( sparsity.off_diagonal_cols ),
                             sparsity.off_diagonal_rows.size ( ) );

    mass_matrix_->InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                                  vec2ptr ( sparsity.diagonal_cols ),
                                  sparsity.diagonal_rows.size ( ),
                                  vec2ptr ( sparsity.off_diagonal_rows ),
                                  vec2ptr ( sparsity.off_diagonal_cols ),
                                  sparsity.off_diagonal_rows.size ( ) );
    rhs_->InitStructure ( );
    sol_->InitStructure ( );
    row_->InitStructure ( );
    coupled_sol_->InitStructure ( );

    // Gets the sparcity structure of the matrices. Only works sequential so far!
    row_entries_ = sparsity.diagonal_rows;
    col_entries_ = sparsity.diagonal_cols;

    // Zero all linear algebra objects.
    matrix_->Zeros ( );
    mass_matrix_->Zeros ( );
    rhs_->Zeros ( );
    sol_->Zeros ( );
    row_->Zeros ( );
    coupled_sol_->Zeros ( );
}

void DirectInverseTutorial::assemble_system ( )
{
    // Assemble matrix and right-hand-side vector.
    LocalDirectInverseTutorialAssembler local_asm;
    LocalDirectInverseTutorialMassAssembler local_mass_asm;
    global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );
    global_asm_.assemble_matrix ( space_, local_mass_asm, *mass_matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, *rhs_ );

    global_asm_.should_reset_assembly_target ( false );
    global_asm_.assemble_matrix_boundary ( space_, local_asm, *matrix_ );
    global_asm_.should_reset_assembly_target ( true );
}

void DirectInverseTutorial::solve_system ( )
{
    LinearSolver<LAD>* solver_;
    LinearSolverFactory<LAD> SolFact;
    solver_ = SolFact.Get (
                            params_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( params_["LinearSolver"] );

    solver_->SetupOperator ( *matrix_ );
    solver_->Solve ( *rhs_, sol_ );

    //////////////////////////// INVERSE PROBLEM //////////////////////////////////
    // Here the solution which was computed above is reduced to a small set of
    // solution points. This refers to some discrete measure point in the reality.
    // After that f is inversely computed by
    //
    // f_alpha = (lambda' * T' * T * lambda + alpha * penalty) lambda' * T' * u_d
    //
    // with lambda being A^(-1)M, with A and M from the direct problem, u_d being
    // the reduced solution, T being a map-matrix from u to u_d and penalty being
    // a diagonal matrix with entries, where the dofs are on the boundary or in the
    // air.

    // Dimension of the (square) mass matrix m and the (square) stiffness matrix a.
    const int n = matrix_->nrows_global ( );

    // Number of cells for the iteration.
    int number_of_cell = master_mesh_->num_entities ( DIMENSION );

    std::vector<int> dof_count_trajectory ( n, -1 );
    std::vector<int> penalty_dof_count ( n, -1 );
    std::vector<std::vector<double> > coord_vec;

    // Computes the values of u_trajectory, which is the array of the sample data.
    // Iteration over all cells.
    for ( int k = 0; k < number_of_cell; k++ )
    {
        space_.dof ( ).get_coord_on_cell ( 0, k, coord_vec );

        int number_dof_on_cells = space_.dof ( ).get_nb_dofs_on_cell ( k );

        for ( int p = 0; p < number_dof_on_cells; p++ )
        {
            // Checks if the dof lies on the trajectory.
            if ( coord_vec[p][2] == 0 || coord_vec[p][2] == 1 )
            {
                // Computes the global dof number of the local dof on the cell.
                int place = space_.dof ( ).mapl2g ( 0, k, p );
                dof_count_trajectory[place] = place;
            }

            // If the dof lies on the facet or if it is in the air, a penalty is later
            // added, to improve the solution.
            if ( ( coord_vec[p][2] >= 0 ) || ( coord_vec[p][0] == 0 ) ||
                 ( coord_vec[p][0] == 20 ) || ( coord_vec[p][1] == 0 ) ||
                 ( coord_vec[p][1] == 10 ) || ( coord_vec[p][2] == -6 ) )
            {
                int place2 = space_.dof ( ).mapl2g ( 0, k, p );
                penalty_dof_count[place2] = place2;
            }
        }
    }

    // Fills the vector u_trajectory with the values at the dofs on the trajectory.
    // The vector has full size (n) and is zero if the dof is not on the trajectory
    u_trajectory.resize ( n );
    for ( int i = 0; i < n; i++ )
    {
        u_trajectory[i] = 0;
        if ( dof_count_trajectory[i] != -1 )
        {
            sol_->GetValues ( &dof_count_trajectory.at ( i ), 1, &u_trajectory.at ( i ) );
        }
    }

    // Initialises the matrices needed in the following.
    SeqDenseMatrix<Scalar> lambda;
    lambda.Resize ( n, n );
    lambda.Zeros ( );

    SeqDenseMatrix<Scalar> penalty_mat;
    penalty_mat.Resize ( n, n );
    penalty_mat.Zeros ( );

    SeqDenseMatrix<Scalar> identity;
    identity.Resize ( n, n );
    lambda.Zeros ( );

    SeqDenseMatrix<Scalar> lambda_trans;
    lambda_trans.Resize ( n, n );
    lambda_trans.Zeros ( );

    SeqDenseMatrix<Scalar> product;
    product.Resize ( n, n );
    product.Zeros ( );

    // Creates lambda. Since A^(-1) is needed, it is retrieved by a series of
    // linear equations.
    // col_entries and row_entries contain the structure of the (sparse) mass
    // matrix, i.e. they contain the indices of all dofs which are not zero.
    for ( int j = 0; j < n; j++ )
    {
        for ( int i = 0; i < static_cast < int > ( col_entries_.size ( ) ); i++ )
        {
            if ( col_entries_[i] == j )
            {
                double value_;
                mass_matrix_->diagonal ( ).get_value ( row_entries_[i], j, &value_ );
                row_->SetValues ( &row_entries_[i], 1, &value_ );
            }
        }

        // Solves A*L = M with L being lambda. Every column of M is used as one
        // RHS to get one column of L.
        solver_->Solve ( *row_, coupled_sol_ );

        for ( int i = 0; i < n; i++ )
        {
            double val;
            coupled_sol_->GetValues ( &i, 1, &val );
            lambda ( j, i ) = val;
        }
        row_->Zeros ( );
        coupled_sol_->Zeros ( );
    }

    // Sets the unnecessary rows of lambda to 0.
    for ( int i = 0; i < n; i++ )
    {
        if ( dof_count_trajectory[i] == -1 )
        {
            for ( int j = 0; j < n; j++ )
            {
                lambda ( i, j ) = 0;
            }
        }
    }

    std::vector<double> right_hand_side ( n, 0 );
    solution_vec.resize ( n );
    lambda.transpose_me ( lambda_trans );
    lambda_trans.MatrixMult ( lambda, product );

    // Computes alpha. Alpha = 0.1 * max(product_ii)
    double alpha = 0;
    for ( int i = 0; i < n; i++ )
    {
        if ( alpha < product ( i, i ) )
        {
            alpha = product ( i, i );
        }
    }
    alpha = 0.1 * alpha;

    double penalty_val = 100;

    for ( int i = 0; i < n; i++ )
    {
        if ( penalty_dof_count[i] != -1 )
        {
            penalty_mat ( i, i ) = penalty_val * alpha;
        }
        else
        {
            penalty_mat ( i, i ) = alpha;
        }
    }

    product.Add ( penalty_mat );
    lambda_trans.VectorMult ( u_trajectory, right_hand_side );
    product.Solve ( right_hand_side, solution_vec );
    delete solver_;
}

void DirectInverseTutorial::visualize ( )
{

    sol_->UpdateCouplings ( );

    // Generate filename.
    std::stringstream input;
    input << "solution" << refinement_level_;
    if ( num_partitions_ > 1 )
    {
        input << ".pvtu";
    }
    else
    {
        input << ".vtu";
    }

    ParallelCellVisualization<double> visu ( space_, params_["Mesh"]["FeDegree"].get<int>( ), comm_, MASTER_RANK );

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
        material_number.at ( it->index ( ) ) = mesh_->get_material_number ( mesh_->tdim ( ), it->index ( ) );
        remote_index.at ( it->index ( ) ) = temp1;
        sub_domain.at ( it->index ( ) ) = temp2;
    }

    visu.visualize ( EvalFeFunction<LAD>( space_, *sol_, 0 ), "f" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );

    if ( rank_ == MASTER_RANK )
    {
        std::cout << "Solution for refinement level " << refinement_level_ <<
                " is computed and can be viewed." << "\n";
    }
}

void DirectInverseTutorial::adapt ( )
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
                                                        comm_ );
        assert ( local_mesh != 0 );
        SharedVertexTable shared_verts;
        mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
    }
}
