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

#include "poisson_uncertainty.h"
#include "common/macros.h"

static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );
    std::string fName;
    if ( argc > 1 )
    {
        fName = argv[1];
    }
    else
    {
        std::cout << "Pass XMl parameter file as argument!" << std::endl;
        std::cout << "Usage: ./poisson_uncertainty <xml-filename>" << std::endl;
        exit ( -1 );
    }

    // Setup timer.
    Timer main_timer;

    int num_partitions_space = 1;
    if ( argc > 2 )
        num_partitions_space = std::atoi ( argv[2] );
    // Read parameters
    PropertyTree config ( fName, MASTER_RANK, MPI_COMM_WORLD );

    // Start Logging
    std::ofstream info_file ( config["Output"]["InfoFilename"].get<std::string>( ).c_str ( ) );
    LogKeeper::get_log ( "info" ).set_target ( &info_file );
    // LogKeeper::get_log("info").set_target(&std::cout);
    int size;
    MPI_Comm_size ( MPI_COMM_WORLD, &size );
    LOG_INFO ( "MPI Processes", size );

    // Compute MPI subgroup for sequential domain treatment
    MPI_Comm comm_space;
    MPI_Comm comm_uq;

    MPI_Group orig_group, new_group_space, new_group_uq;
    MPI_Comm_group ( MPI_COMM_WORLD, &orig_group );

    int num_processors = size;
    int num_partitions_uq = num_processors / num_partitions_space;

    int my_rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );

    int* ranks = new int[num_partitions_space];
    for ( int i = 0; i < num_partitions_space; ++i )
        ranks[i] = my_rank % num_partitions_uq + i * num_partitions_uq;

    MPI_Group_incl ( orig_group, num_partitions_space, ranks, &new_group_space );
    MPI_Comm_create ( MPI_COMM_WORLD, new_group_space, &comm_space );

    delete[] ranks;
    ranks = new int[num_partitions_uq];
    for ( int i = 0; i < num_partitions_uq; ++i )
        ranks[i] = my_rank - my_rank % num_partitions_uq + i;

    MPI_Group_incl ( orig_group, num_partitions_uq, ranks, &new_group_uq );
    MPI_Comm_create ( MPI_COMM_WORLD, new_group_uq, &comm_uq );
    delete[] ranks;

    // Run application
    PoissonMPI application ( config, comm_space, comm_uq );
    application.run ( );

    // flush log here to avoid problems
    LogKeeper::get_log ( "info" ).flush ( );

    // Stop timer.
    main_timer.stop ( );
    LOG_INFO ( "Total program run time [without MPI Init and Finalize]",
               main_timer );

    MPI_Finalize ( );
    return 0;
}

PoissonMPI::PoissonMPI ( PropertyTree &config, const MPI_Comm& comm_space, const MPI_Comm& comm_uq )
: comm_space_ ( comm_space ),
comm_uq_ ( comm_uq ),
config_ ( &config ),
master_mesh_ ( 0 )
{
    MPI_Comm_size ( comm_space_, &num_partitions_ );
    MPI_Comm_rank ( comm_space_, &rank_ );

    setup_mesh ( );
    setup_application ( );
    setup_space ( );
    setup_system ( );
    setup_la_structures ( );
    setup_bc ( );
}

PoissonMPI::~PoissonMPI ( )
{
    stop_platform ( la_sys_ );

    for ( int i = 0; i < pcbasis_.N ( ) + 1; ++i )
        delete galerkin_matrix_[i];
    delete system_matrix_;

    delete galerkin_sol_;
    delete galerkin_res_;
    delete galerkin_cor_;

    delete space_;
}

void PoissonMPI::setup_mesh ( )
{
    assert ( DIM == 2 || DIM == 3 );
    LOG_INFO ( "Problem dimension", DIM );
    ( *config_ )["Mesh"]["RefinementLevel"].read<int>( refinement_level_ );
    LOG_INFO ( "Refinement level", refinement_level_ );

    if ( rank_ == MASTER_RANK )
    {
        // read mesh (sequential, dimension DIM
        std::string mesh_filename = std::string ( DATADIR ) + ( *config_ )["Mesh"]["Filename"].get<std::string>( );
        master_mesh_ = read_mesh_from_file ( mesh_filename, DIM, DIM, 0 );

        for ( int r = 0; r < refinement_level_; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
        }
    }

    MeshPtr local_mesh = partition_and_distribute ( master_mesh_, MASTER_RANK, comm_space_ );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_space_, shared_verts );

    std::ostringstream rank_str;
    rank_str << rank_;
    PVtkWriter writer ( comm_space_ );
    std::string output_file = ( *config_ )["Output"]["MeshFilename"].get<std::string>( );
    writer.write ( output_file.c_str ( ), *mesh_ );
}

void PoissonMPI::setup_application ( )
{
    // setup stochastic parameters
    ( *config_ )["PolynomialChaos"]["No"].read<int>( No_ );
    LOG_INFO ( "Polynomial Chaos order", No_ );
    ( *config_ )["PolynomialChaos"]["N"].read<int>( N_ );
    LOG_INFO ( "Dimension of random space", N_ );
    // Polynomial degree of random input -> here: linear random input
    No_input_ = 1;
    LOG_INFO ( "Polynomial degree of random input", No_input_ );
    std::vector<PCBasis::Distribution> dist ( N_, PCBasis::UNIFORM );
    LOG_INFO ( "Probability distribution", "UNIFORM" );

    pcbasis_.Init ( N_, No_, dist );
    pctensor_.SetMPIComm ( comm_uq_ );
    pctensor_.ComputeTensor ( No_input_, pcbasis_ );

    double nu0, sigma, q;
    ( *config_ )["Application"]["MeanViscosity"].read<double>( nu0 );
    ( *config_ )["Application"]["Variability"].read<double>( sigma );
    ( *config_ )["Application"]["Decay"].read<double>( q );

    nu_.resize ( N_ + 1, 0.0 );
    nu_[0] = nu0;
    for ( int i = 1; i < N_ + 1; ++i )
        nu_[i] = sigma * nu0 * pow ( q, i - 1.0 );
}

void PoissonMPI::setup_space ( )
{
    // setup space Q1 elements
    space_ = new VectorSpace<double>( comm_space_ );
    int p_degree = 0;
    ( *config_ )["PolynomialChaos"]["q"].read<int>( p_degree );
    std::vector<int> degrees ( 1, p_degree );
    space_->Init ( degrees, *mesh_ );

    std::cout << "Rank: " << pctensor_.MyRank ( ) << " using NDofs = " <<
            space_->dof ( ).get_nb_dofs ( )*( pctensor_.SizeLocal ( ) + pctensor_.SizeOffModes ( ) ) << "\n";

    if ( pctensor_.MyRank ( ) == 0 )
        std::cout << "Number of Total Dofs: " << space_->dof ( ).get_nb_dofs ( ) << "x"
        << pctensor_.Size ( ) << " = " << ( long ) space_->dof ( ).get_nb_dofs ( ) * pctensor_.Size ( ) << "\n";
}

void PoissonMPI::setup_system ( )
{
    la_sys_.Platform = CPU;
    la_sys_.rank = 0;
    init_platform ( la_sys_ );
    matrix_impl_ = NAIVE;
    vector_impl_ = NAIVE;
    la_matrix_format_ = CSR;
    matrix_precond_ = NOPRECOND;
    la_sys_.GPU_CUBLAS = false;
}

void PoissonMPI::setup_la_structures ( )
{
    // Initialize linear algebra structures
    couplings_.Clear ( );
    couplings_.Init ( comm_space_, space_->dof ( ) );

    // compute matrix graph to build matrix structure
    std::vector<int> diagonal_rows, diagonal_cols, off_diagonal_rows, off_diagonal_cols;
    std::vector < std::vector<bool> > var_coupling ( 1 );
    var_coupling[0].resize ( 1, true );
    InitStructure ( *space_, &diagonal_rows, &diagonal_cols,
                    &off_diagonal_rows, &off_diagonal_cols, &var_coupling );

    couplings_.InitializeCouplings ( off_diagonal_rows, off_diagonal_cols );

    // Initialize matrices and vectors
    matrix_.Init ( comm_space_, couplings_, la_sys_.Platform,
                   matrix_impl_, la_matrix_format_ );

    matrix_.InitStructure ( vec2ptr ( diagonal_rows ),
                            vec2ptr ( diagonal_cols ),
                            diagonal_rows.size ( ),
                            vec2ptr ( off_diagonal_rows ),
                            vec2ptr ( off_diagonal_cols ),
                            off_diagonal_rows.size ( ) );

    sol_.Init ( comm_space_, couplings_, la_sys_.Platform, vector_impl_ );
    sol_.InitStructure ( );
    sol_.Zeros ( );

    galerkin_matrix_.resize ( pcbasis_.N ( ) + 1 );
    for ( int i = 0; i < pcbasis_.N ( ) + 1; ++i )
        galerkin_matrix_[i] = new CoupledMatrix<double>;

    galerkin_sol_tmp_.resize ( pctensor_.SizeLocal ( ) + pctensor_.SizeOffModes ( ) );

    for ( int i = 0; i<static_cast < int > ( galerkin_sol_tmp_.size ( ) ); ++i )
    {
        galerkin_sol_tmp_[i] = new CoupledVector<double>;
        galerkin_sol_tmp_[i]->CloneFrom ( sol_ );
    }

    galerkin_sol_ = new PCGalerkinVector<double>;
    galerkin_res_ = new PCGalerkinVector<double>;
    galerkin_cor_ = new PCGalerkinVector<double>;

    galerkin_sol_->SetMPIComm ( comm_uq_ );
    galerkin_res_->SetMPIComm ( comm_uq_ );
    galerkin_cor_->SetMPIComm ( comm_uq_ );

    galerkin_sol_->SetModes ( galerkin_sol_tmp_, pctensor_.SizeLocal ( ), pctensor_.SizeOffModes ( ) );
    galerkin_sol_->Zeros ( );
    galerkin_res_->CloneFrom ( *galerkin_sol_ );
    galerkin_cor_->CloneFrom ( *galerkin_sol_ );

    for ( int mode = 0; mode<static_cast < int > ( galerkin_sol_tmp_.size ( ) ); ++mode )
    {
        delete galerkin_sol_tmp_[mode];
    }

    system_matrix_ = new PCGalerkinMatrix<double>( sol_ );
}

struct PoissonMPIDirichletBC2d
{

    PoissonMPIDirichletBC2d ( )
    {
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );
        const bool inflow = ( material_num == 11 );
        if ( inflow )
        {
            values.resize ( coords_on_face.size ( ), 1.0 );
        }
        return values;
    }
};

struct PoissonMPIDirichletBC3d
{

    PoissonMPIDirichletBC3d ( )
    {
    }

    std::vector<double> evaluate ( const Entity& face, const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values;

        const int material_num = face.get_material_number ( );
        const bool bdy = ( material_num == 10 );
        if ( bdy )
        {
            values.resize ( coords_on_face.size ( ), 0.0 );
        }
        return values;
    }
};

void PoissonMPI::setup_bc ( )
{
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    if ( DIM == 2 )
    {
        PoissonMPIDirichletBC2d bc;
        compute_dirichlet_dofs_and_values ( bc, *space_, 0,
                                            dirichlet_dofs_, dirichlet_values_ );

    }
    else
    {
        PoissonMPIDirichletBC3d bc;
        compute_dirichlet_dofs_and_values ( bc, *space_, 0,
                                            dirichlet_dofs_, dirichlet_values_ );
    }

    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        galerkin_sol_->Mode ( 0 )->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }
}

void PoissonMPI::setup_linear_solver ( )
{
    std::string solver = ( *config_ )["GalerkinLinearSolver"]["LinearSolver"].get<std::string>( );
    if ( solver == "CG" )
    {
        cg_solver_.InitControl ( 1000, 1.e-12, 1.e-12, 1.e6 );
        cg_solver_.SetupOperator ( *system_matrix_ );
        linear_solver_ = &cg_solver_;
    }
    else if ( solver == "Multilevel" )
    {
        linear_solver_ = new PCMultilevelSolver<LADescriptorPolynomialChaosD>( *config_ );
        linear_solver_->InitControl ( 100, 1.e-12, 1.e-12, 1.e6 );
        linear_solver_->SetupOperator ( *system_matrix_ );
        linear_solver_->Build ( );
    }
    else
    {
        std::cout << "ERROR: No solver defined!\n";
        interminable_assert ( 0 );
    }

    std::string preconditioner = ( *config_ )["GalerkinLinearSolver"]["Preconditioner"].get<std::string>( );
    if ( solver == "CG" && preconditioner == "Mean" )
    {
#ifdef WITH_UMFPACK
        mean_precond_ = new MeanbasedPreconditioner<LADescriptorPolynomialChaosD>;
        mean_precond_->SetupOperator ( *system_matrix_ );
        cg_solver_.InitParameter ( "Preconditioning" );
        cg_solver_.SetupPreconditioner ( *mean_precond_ );
#else
        interminable_assert ( 0 );
#endif
    }
    else if ( solver == "CG" && preconditioner == "ML" )
    {
        ml_precond_ = new PCMultilevelSolver<LADescriptorPolynomialChaosD>( *config_ );
        ml_precond_->SetupOperator ( *system_matrix_ );
        ml_precond_->Build ( );

        cg_solver_.InitParameter ( "Preconditioning" );
        cg_solver_.SetupPreconditioner ( *ml_precond_ );
    }
}

void PoissonMPI::visualize_solution ( LAD::VectorType& u, std::string const& filename ) const
{
    u.UpdateCouplings ( );
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( *( space_ ), num_intervals, comm_space_, 0 );

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

    std::stringstream input;
    input << filename;

    if ( num_partitions_ > 1 )
        input << ".pvtu";
    else
        input << ".vtu";

    visu.visualize ( EvalFeFunction<LAD>( *( space_ ), u ), "val" );

    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );
}

void PoissonMPI::compute_matrix ( )
{
    StandardGlobalAssembler<double> global_asm;

    for ( int mode = 0; mode < pcbasis_.N ( ) + 1; ++mode )
    {
        local_mode_asm_ = new PoissonMPIModeAssembler ( mode, nu_[mode] );
        global_asm.assemble_matrix ( *space_, *local_mode_asm_, matrix_ );

        // correct BC -- set Dirichlet rows
        if ( !dirichlet_dofs_.empty ( ) )
        {
            if ( mode == 0 )
            {
                matrix_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
            }
            else
            {
                matrix_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 0. );
            }
        }

        galerkin_matrix_[mode]->CloneFrom ( matrix_ );
        delete local_mode_asm_;
    }

    system_matrix_->SetNumThreads ( 1 );
    system_matrix_->SetMatrix ( galerkin_matrix_ );
    system_matrix_->SetTensor ( &pctensor_ );
}

void PoissonMPI::compute_residual ( )
{
    system_matrix_->VectorMult ( *galerkin_sol_, galerkin_res_ );
    galerkin_res_->Scale ( -1.0 );

    if ( pctensor_.MyRank ( ) == 0 )
    {
        StandardGlobalAssembler<double> global_asm;

        CoupledVector<double> mean_rhs;
        mean_rhs.CloneFrom ( sol_ );
        mean_rhs.Zeros ( );

        local_mode_asm_ = new PoissonMPIModeAssembler ( 0, nu_[0] );
        // Assemble mean rhs
        global_asm.assemble_vector ( *space_, *local_mode_asm_, mean_rhs );
        delete local_mode_asm_;

        galerkin_res_->ModeAxpy ( 0, mean_rhs, 1.0 );
    }

    if ( !dirichlet_dofs_.empty ( ) )
    {
        std::vector<LAD::DataType> zeros ( dirichlet_dofs_.size ( ), 0. );
        for ( int mode = 0; mode < pctensor_.SizeLocal ( ); ++mode )
        {
            galerkin_res_->Mode ( mode )->SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( zeros ) );
        }
    }
}

void PoissonMPI::run ( )
{
    compute_matrix ( );
    compute_residual ( );

    setup_linear_solver ( );

    Timer timer;
    timer.start ( );
    solve_linear_system ( );
    timer.stop ( );
    if ( pctensor_.MyRank ( ) == 0 )
        std::cout << "Solving linear system = " << timer.get_duration ( ) << "\n";

    std::string visu_folder = ( *config_ )["Output"]["VisuFolder"].get<std::string>( );

    for ( int mode = 0; mode < pctensor_.SizeLocal ( ); ++mode )
    {
        std::stringstream input;
        input << visu_folder;
        int glo_mode = pctensor_.L2G ( mode );
        input << "uq_poisson_mode_" << glo_mode;
        visualize_solution ( *galerkin_sol_->Mode ( mode ), input.str ( ) );
    }
}

void PoissonMPI::solve_linear_system ( )
{
    linear_solver_->Solve ( *galerkin_res_, galerkin_cor_ );
    galerkin_sol_->Axpy ( *galerkin_cor_, 1.0 );
}

void PoissonMPIModeAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    if ( mode_ == 0 )
    {
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::fabs ( detJ ( q ) );
            // assemble a1(u,v) = \int \alpha_1 * {\grad(u) : \grad(v)}
            for ( int i = 0; i < num_dofs ( 0 ); ++i )
            {
                for ( int j = 0; j < num_dofs ( 0 ); ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                            wq * nu_ * dJ * dot ( grad_phi ( j, q, 0 ), grad_phi ( i, q, 0 ) );
                }
            }
        }
    }
    else
    {
        double pi = acos ( 0.0 )*2.0;
        for ( int q = 0; q < num_q; ++q )
        {
            const double wq = w ( q );
            const double dJ = std::fabs ( detJ ( q ) );
            // assemble a1(u,v) = \int \alpha_1 * {\grad(u) : \grad(v)}
            for ( int i = 0; i < num_dofs ( 0 ); ++i )
            {
                for ( int j = 0; j < num_dofs ( 0 ); ++j )
                {
                    lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                            wq * nu_ * sin ( 2.0 * pi * x ( q )[0] * mode_ ) * sin ( 2.0 * pi * x ( q )[1] * mode_ ) * dJ * dot ( grad_phi ( j, q, 0 ), grad_phi ( i, q, 0 ) );
                }
            }
        }
    }
}

void PoissonMPIModeAssembler::operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );
    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::fabs ( detJ ( q ) );
        // assemble a1(u,v) = \int  * {\grad(u) : \grad(v)}
        for ( int i = 0; i < num_dofs ( 0 ); ++i )
        {
            lv[dof_index ( i, 0 )] +=
                    wq * dJ * phi ( i, q, 0 );
        }
    }
}
