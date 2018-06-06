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

#include "poisson_periodic.h"

static const char* PARAM_FILENAME = "poisson_periodic.xml";
#ifndef MESHES_DATADIR
#    define MESHES_DATADIR "./"
#endif
static const char* DATADIR = MESHES_DATADIR;

// Main application class ///////////////////////////////////

class PoissonPeriodic
{
  public:

    PoissonPeriodic ( const std::string& param_filename, const std::string& path_mesh )
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
            compute_error ( );

            // Visualize the solution and the errors.
            visualize ( );

            // Modify the space through refinement. Set is_done_ = true when finished.
            adapt ( );

        }
    }

    ~PoissonPeriodic ( )
    {
        delete matrix_;
        delete sol_;
        delete rhs_;
    }

  private:
    // Member functions

    // Read and distribute mesh.
    std::string path_mesh;
    // vector that stores information about periodicity
    std::vector<MasterSlave> period;

    void build_initial_mesh ( );

    // set correct geometry in transformations
    void periodify_space ( );

    void prepare_periodicity ( );

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
    VectorSpace<double> space_;

    // Linear algebra couplings helper object.
    Couplings<double> couplings_;
    // Vectors for solution and load vector.
    CoupledVector<Scalar>* rhs_, *sol_;
    // System matrix.
    CoupledMatrix<Scalar>* matrix_;

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
}; // end class PoissonPeriodic

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
        std::ofstream info_log ( "poisson_tutorial_info_log" );
        LogKeeper::get_log ( "info" ).set_target ( &info_log );
        std::ofstream debug_log ( "poisson_tutorial_debug_log" );
        LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

        // Create application object and run it
        PoissonPeriodic app ( param_filename, path_mesh );
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

void PoissonPeriodic::prepare_periodicity ( )
{
    const int periodic_boundaries = params_["Mesh"]["PeriodicBoundaries"].get<int>( );
    period.clear ( );

    switch ( periodic_boundaries )
    {
        case 0: // no periodic boundaries
        {
            break;
        }
        case 1: // periodic boundaries in x-direction
        {
            period.push_back ( MasterSlave ( 1., 0., .5, 0 ) );
            break;
        }
        case 2: // periodic boundaries in y-direction
        {
            period.push_back ( MasterSlave ( 1., 0., .5, 1 ) );
            break;
        }
        case 3: // periodic boundaries in x- and y-direction
        {
            period.push_back ( MasterSlave ( 1., 0., .5, 0 ) );
            period.push_back ( MasterSlave ( 1., 0., .5, 1 ) );
            break;
        }
        default:
        {
            exit ( -1 );
            break;
        }
    }
}

//////////////// PoissonPeriodic implementation //////////////

void PoissonPeriodic::build_initial_mesh ( )
{
    // Read in the mesh on the master process. The mesh is chosen according to the dimension of the problem.

    prepare_periodicity ( );

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
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIMENSION, DIMENSION, 0, period );

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
        std::string output_file = std::string ( "poisson_periodic_mesh.pvtu" );
        writer.add_all_attributes ( *mesh_, true );
        writer.write ( output_file.c_str ( ), *mesh_ );

    }
}

void PoissonPeriodic::prepare_system ( )
{

    // Assign degrees to each element.
    const int fe_degree = params_["Mesh"]["FeDegree"].get<int>( );
    std::vector< int > degrees ( 1, fe_degree );

    // Initialize the VectorSpace object.
    space_.Init ( degrees, *mesh_ );

    // Modify space to account for periodic boundaries
    periodify_space ( );

    // Setup couplings object.
    couplings_.Init ( comm_, space_.dof ( ) );

    // Compute the matrix graph.
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );

    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows,
                                     sparsity.off_diagonal_cols );

    // Setup linear algebra objects.
    CoupledMatrixFactory<Scalar> CoupMaFact;
    matrix_ = CoupMaFact.Get (
                               params_["LinearAlgebra"]["NameMatrix"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    matrix_->Init ( comm_, couplings_ );
    CoupledVectorFactory<Scalar> CoupVecFact;
    rhs_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
            params ( params_["LinearAlgebra"] );
    sol_ = CoupVecFact.Get (
                             params_["LinearAlgebra"]["NameVector"].get<std::string>( ) )->
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

        if ( period.size ( ) < DIMENSION )
        {
            compute_dirichlet_dofs_and_values ( zero, space_, 0, dirichlet_dofs_,
                                                dirichlet_values_ );
        }
    }

}

void PoissonPeriodic::assemble_system ( )
{
    // Assemble matrix and right-hand-side vector.
    LocalPoissonAssembler local_asm;
    global_asm_.assemble_matrix ( space_, local_asm, *matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, *rhs_ );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // Correct Dirichlet dofs.
        matrix_->diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1.0 );
        rhs_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
        sol_->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                          vec2ptr ( dirichlet_values_ ) );
    }
    rhs_->UpdateCouplings ( );
    sol_->UpdateCouplings ( );

}

void PoissonPeriodic::solve_system ( )
{
    LinearSolver<LAD>* solver_;
    LinearSolverFactory<LAD> SolFact;
    solver_ = SolFact.Get (
                            params_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( params_["LinearSolver"] );
    solver_->SetupOperator ( *matrix_ );
    solver_->Solve ( *rhs_, sol_ );

    sol_->UpdateCouplings ( );

    delete solver_;
}

void PoissonPeriodic::compute_error ( )
{
    L2_err_.clear ( );
    H1_err_.clear ( );
    const int periodic_boundaries = params_["Mesh"]["PeriodicBoundaries"].get<int>( );

    if ( periodic_boundaries == 0 )
    {
        std::cout << "Unique solution exists, but no analytical expression available for error computation!" << std::endl;
    }
    else if ( periodic_boundaries == 3 )
    {
        std::cout << "The solution is not unique in this case!" << std::endl;
    }
    else
    {

        sol_->UpdateCouplings ( );
        // Compute square of the L2 error on each element, putting the
        // values into L2_err_.
        L2ErrorIntegrator L2_int ( *( sol_ ), periodic_boundaries );

        global_asm_.assemble_scalar ( space_, L2_int, L2_err_ );

        // Create attribute with L2 error for output.
        AttributePtr L2_err_attr ( new DoubleAttribute ( L2_err_ ) );

        mesh_->add_attribute ( "L2 error", DIMENSION, L2_err_attr );

        double total_L2_err = std::accumulate ( L2_err_.begin ( ), L2_err_.end ( ), 0. );
        double global_L2_err = 0.;
        MPI_Reduce ( &total_L2_err, &global_L2_err, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK,
                     comm_ );
        LOG_INFO ( "error", "Local L2 error on partition " << rank_ << " = "
                   << std::sqrt ( total_L2_err ) );

        // Compute square of the H1 error on each element, putting the
        // values into H1_err_.

        H1ErrorIntegrator H1_int ( *( sol_ ), periodic_boundaries );
        global_asm_.assemble_scalar ( space_, H1_int, H1_err_ );

        // Create attribute with H1 error for output.
        AttributePtr H1_err_attr ( new DoubleAttribute ( H1_err_ ) );
        mesh_->add_attribute ( "H1 error", DIMENSION, H1_err_attr );
        double total_H1_err = std::accumulate ( H1_err_.begin ( ), H1_err_.end ( ), 0. );
        double global_H1_err = 0.;
        MPI_Reduce ( &total_H1_err, &global_H1_err, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK,
                     comm_ );
        LOG_INFO ( "error", "Local H1 error on partition " << rank_ << " = "
                   << std::sqrt ( total_H1_err ) );

        if ( rank_ == MASTER_RANK )
        {
            std::cout << "Global L2 error = " << std::sqrt ( global_L2_err ) << "\n"
                    << "Global H1 error = " << std::sqrt ( global_H1_err ) << "\n";
        }
    }
}

void PoissonPeriodic::visualize ( )
{

    ParallelCellVisualization<double> visu ( space_, 1, comm_, MASTER_RANK );

    // Generate filename.
    std::stringstream name;
    name << "solution" << refinement_level_;

    std::vector<double> remote_index ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );
    std::vector<double> sub_domain ( mesh_->num_entities ( mesh_->tdim ( ) ), 0 );

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
    }
    // Update not nescessary since we update in compute_error
    //sol_->UpdateCouplings();

    visu.visualize ( EvalFeFunction<LAD>( space_, *( sol_ ) ), "u" );

    // visualize error measures
    const int periodic_boundaries = params_["Mesh"]["PeriodicBoundaries"].get<int>( );

    if ( periodic_boundaries == 1 || periodic_boundaries == 2 )
    {
        visu.visualize_cell_data ( L2_err_, "L2" );
        visu.visualize_cell_data ( H1_err_, "H1" );
    }
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( name.str ( ) );
}

void PoissonPeriodic::adapt ( )
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

void PoissonPeriodic::periodify_space ( )
{
    assert ( mesh_ != 0 );

    const int tdim = mesh_->tdim ( );
    // loop over all cells and init the cell transformations

    for ( mesh::EntityIterator it = mesh_->begin ( tdim );
          it != mesh_->end ( tdim );
          ++it )
    {
        Coord coord_vtx;
        it->get_coordinates ( coord_vtx );
        coord_vtx = unperiodify ( coord_vtx, it->gdim ( ), period );
        space_.fe_manager ( ).get_cell_transformation ( it->index ( ) )->reinit ( coord_vtx );
    }
}
