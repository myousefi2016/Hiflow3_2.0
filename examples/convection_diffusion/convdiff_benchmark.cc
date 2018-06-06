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

#include "convdiff_benchmark.h"

#include <sys/time.h>

#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "hiflow.h"

static const char* DATADIR = MESHES_DATADIR;
static const int MASTER_RANK = 0;
static const int DEBUG_LEVEL = 0;
const TDim tdim = DIM;
const GDim gdim = DIM;

/// Program main loop: setup MPI, read parameters and start the application

int main ( int argc, char** argv )
{
    MPI_Init ( &argc, &argv );

    Timer main_timer;

    // Read parameters
    std::string param_filename ( "convdiff_benchmark.xml" );
    if ( argc > 1 )
    {
        param_filename = std::string ( argv[1] );
    }
    try
    {
        ConvDiff app ( param_filename );
        app.run ( );
    }
    catch ( const std::exception& e )
    {
        std::cerr << "\nProgram ended with uncaught exception.\n";
        std::cerr << e.what ( ) << "\n";
        return -1;
    }

    main_timer.stop ( );
    std::cout << "Total program run time [without MPI Init and Finalize]" <<
            main_timer << std::endl;

    MPI_Finalize ( );
    return 0;
}

// local assembling of system matrix

void ConvDiffAssembler::operator() ( const Element<double>& element,
        const Quadrature<double>& quadrature,
        LocalMatrix& lm )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );

    // warum???
    Vec<DIM, double> beta;
    for ( int var = 0; var < DIM; ++var )
    {
        beta[var] = beta_[var];
    }
    // warum??
    lm.Resize ( total_dofs, total_dofs );
    lm.Zeros ( );

    // loop q
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );

        // assemble a1(u,v) = \int \nu * {\grad(u) : \grad(v)}
        for ( int i = 0; i < num_dofs ( 0 ); ++i )
        {
            for ( int j = 0; j < num_dofs ( 0 ); ++j )
            {
                lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                        wq * ( nu_ * dot ( grad_phi ( j, q, 0 ), grad_phi ( i, q, 0 ) ) ) * dJ;
            }
        }

        // assemble a2(u,v) = \int { (beta*\grad{u})*v }
        for ( int i = 0; i < num_dofs ( 0 ); ++i )
        {
            for ( int j = 0; j < num_dofs ( 0 ); ++j )
            {
                lm ( dof_index ( i, 0 ), dof_index ( j, 0 ) ) +=
                        wq * ( dot ( beta, grad_phi ( j, q, 0 ) ) * phi ( i, q, 0 ) ) * dJ;
            }
        }
    }
}

// local assembling of rhs vector

void ConvDiffAssembler::operator() ( const Element<double>& element,
        const Quadrature<double>& quadrature,
        LocalVector& lv )
{
    AssemblyAssistant<DIM, double>::initialize_for_element ( element, quadrature );

    const int num_q = num_quadrature_points ( );
    const int total_dofs = num_dofs_total ( );
    lv.clear ( );
    lv.resize ( total_dofs, 0. );

    // loop over quadrature points
    for ( int q = 0; q < num_q; ++q )
    {
        const double wq = w ( q );
        const double dJ = std::abs ( detJ ( q ) );
        const double gamma = 1.0 / nu_;

        const double x1 = x ( q )[0];
        const double x2 = x ( q )[1];
        // l0(v) = \int( dot(u_n - u_k, v))
        for ( int i = 0; i < num_dofs ( 0 ); ++i )
        {
            lv[dof_index ( i, 0 )] +=
                    wq * ( gamma * ( -nu_ * ( gamma * exp ( ( x1 - 1 ) * gamma )
                    * ( x1 * gamma + 2 )
                    * sin ( M_PI * x2 ) + M_PI * M_PI * x1
                    * ( -1.0 * exp ( -1.0 * gamma ) )
                    * ( exp ( x1 * gamma ) - exp ( gamma ) )
                    * sin ( M_PI * x2 ) )
                    + ( beta_[0] ) * ( M_PI * x1 * exp ( -gamma )
                    * ( exp ( x1 * gamma ) - exp ( gamma ) )
                    * cos ( M_PI * x2 ) )
                    + ( beta_[1] ) * ( sin ( M_PI * x2 )
                    * ( exp ( ( x1 - 1 ) * gamma )
                    * ( x1 * gamma + 1 ) - 1 ) ) )
                    * phi ( i, q, 0 ) ) * dJ;
        }
    }
}

///////////// ConvectionDiffusionApplication /////////////

ConvDiff::ConvDiff ( const std::string& param_filename )
: param_ ( param_filename.c_str ( ), MASTER_RANK, MPI_COMM_WORLD ),
comm_ ( MPI_COMM_WORLD ),
rank_ ( -1 ),
num_partitions_ ( -1 ),
master_rank_ ( 0 ),
init_platform_ ( true )
{

    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    nu_ = param_["Model"]["nu"].get<Scalar>( );
    beta_.resize ( DIM );
    beta_[0] = param_["Model"]["beta1"].get<Scalar>( );
    beta_[1] = param_["Model"]["beta2"].get<Scalar>( );
    if ( DIM == 3 )
    {
        beta_[2] = param_["Model"]["beta3"].get<Scalar>( );
    }
}

ConvDiff::~ConvDiff ( )
{
    // Delete Platform Mat/Vec
    if ( la_sys_.Platform != CPU )
    {
        // matrix
        delete dev_matrix_;

        // vector
        delete dev_sol_;
        delete dev_rhs_;
    }
}

void ConvDiff::run ( )
{
    // Start Logging
    std::ofstream info_file ( param_["Output"]["LogFilename"].get<std::string>( ).c_str ( ) );
    LogKeeper::get_log ( "info" ).set_target ( &info_file );

    LOG_INFO ( "MPI Processes", num_partitions_ );
    // flush log here to avoid problems
    LogKeeper::get_log ( "info" ).flush ( );

    // setup timing report
    TimingScope::set_report ( &time_report_ );

    read_and_distribute_mesh ( );
    setup_system ( );
    prepare_space ( );
    prepare_bc ( );
    prepare_lin_alg_structures ( );
    assemble_system ( );
    prepare_linear_solver ( );
    solve_system ( );
    visualize_solution ( );

    if ( rank_ == MASTER_RANK )
    {
        // Output time report
        TimingReportOutputVisitor visitor ( std::cout );
        time_report_.traverse_depth_first ( visitor );
    }
}

void ConvDiff::read_and_distribute_mesh ( )
{
    TimingScope tscope ( "Read, refine and distribute mesh." );
    refinement_level_ = param_["Mesh"]["RefinementLevel"].get<int>( );
    std::string mesh_filename =
            std::string ( DATADIR ) + param_["Mesh"]["Filename"].get<std::string>( );
    std::cerr << "Filename " << mesh_filename << "\n";

    MeshPtr master_mesh ( 0 );

    // read in mesh and refine globally
    if ( rank_ == master_rank_ )
    {
        master_mesh = read_mesh_from_file ( mesh_filename, DIM, DIM, 0 );
        assert ( master_mesh != 0 );

        for ( int r = 0; r < refinement_level_; ++r )
        {
            master_mesh = master_mesh->refine ( );
        }
    }

    // Distribute mesh
    MeshPtr local_mesh = partition_and_distribute ( master_mesh, master_rank_, MPI_COMM_WORLD );
    assert ( local_mesh != 0 );

    // Compute ghost cells
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, MPI_COMM_WORLD, shared_verts );
}

void ConvDiff::setup_system ( )
{
    TimingScope tscope ( "Setup parameters for linear system." );
    const std::string platform_str =
            param_["LinearAlgebra"]["Platform"].get<std::string>( );
    if ( platform_str == "CPU" )
    {
        la_sys_.Platform = CPU;
    }
    else if ( platform_str == "GPU" )
    {
        la_sys_.Platform = GPU;
    }
    else if ( platform_str == "OpenCL" )
    {
        la_sys_.Platform = OPENCL;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.Platform", platform_str );
    }
    la_sys_.rank = rank_;
    init_platform ( la_sys_ );

    const std::string matrix_impl_str =
            param_["LinearAlgebra"]["MatrixImplementation"].get<std::string>( );
    if ( matrix_impl_str == "Naive" )
    {
        matrix_impl_ = NAIVE;
    }
    else if ( matrix_impl_str == "BLAS" )
    {
        matrix_impl_ = BLAS;
    }
    else if ( matrix_impl_str == "MKL" )
    {
        matrix_impl_ = MKL;
    }
    else if ( matrix_impl_str == "OPENMP" )
    {
        matrix_impl_ = OPENMP;
    }
    else if ( matrix_impl_str == "SCALAR" )
    {
        matrix_impl_ = SCALAR;
    }
    else if ( matrix_impl_str == "SCALAR_TEX" )
    {
        matrix_impl_ = SCALAR_TEX;
    }
    else if ( matrix_impl_str == "OpenCL" )
    {
        matrix_impl_ = OPEN_CL;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.MatrixImplementation",
                                         matrix_impl_str );
    }

    const std::string vector_impl_str =
            param_["LinearAlgebra"]["VectorImplementation"].get<std::string>( );
    if ( vector_impl_str == "Naive" )
    {
        vector_impl_ = NAIVE;
    }
    else if ( vector_impl_str == "BLAS" )
    {
        vector_impl_ = BLAS;
    }
    else if ( vector_impl_str == "MKL" )
    {
        vector_impl_ = MKL;
    }
    else if ( vector_impl_str == "OPENMP" )
    {
        vector_impl_ = OPENMP;
    }
    else if ( vector_impl_str == "SCALAR" )
    {
        vector_impl_ = SCALAR;
    }
    else if ( vector_impl_str == "SCALAR_TEX" )
    {
        vector_impl_ = SCALAR_TEX;
    }
    else if ( vector_impl_str == "OpenCL" )
    {
        vector_impl_ = OPEN_CL;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.VectorImplementation",
                                         vector_impl_str );
    }

    const std::string matrix_str =
            param_["LinearAlgebra"]["MatrixFormat"].get<std::string>( );
    if ( matrix_str == "DENSE" )
    {
        la_matrix_format_ = DENSE;
    }
    else if ( matrix_str == "CSR" )
    {
        la_matrix_format_ = CSR;
    }
    else if ( matrix_str == "COO" )
    {
        la_matrix_format_ = COO;
    }
    else if ( matrix_str == "ELL" )
    {
        la_matrix_format_ = ELL;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.MatrixFormat", matrix_str );
    }

    const std::string precond_str =
            param_["LinearAlgebra"]["MatrixFreePrecond"].get<std::string>( );
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
        matrix_precond_ = ILU2;
    }
    else
    {
        throw UnexpectedParameterValue ( "LinearAlgebra.MatrixFreePrecond",
                                         precond_str );
    }

    if ( la_sys_.Platform == GPU )
    {
        la_sys_.GPU_CUBLAS = true;
    }
    else
    {
        la_sys_.GPU_CUBLAS = false;
    }
}

void ConvDiff::prepare_space ( )
{
    TimingScope tscope ( "Prepare space." );
    int degree = param_["FiniteElements"]["Degree"].get<int>( );
    std::vector<int> degrees ( 1, degree );
    space_.Init ( degrees, *mesh_ );
}

void ConvDiff::prepare_bc ( )
{
    TimingScope tscope ( "Prepare boundary conditions." );
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    Dirichlet_0_BC bc;
    compute_dirichlet_dofs_and_values ( bc, space_, 0, dirichlet_dofs_,
                                        dirichlet_values_ );
}

// initialize couplings, set up structure of matrices and vectors

void ConvDiff::prepare_lin_alg_structures ( )
{
    TimingScope tscope ( "Prepare linear algebra structure." );
    // Initialize linear algebra structures
    if ( init_platform_ )
    {
        couplings_.Init ( comm_, space_.dof ( ) );
    }

    // compute matrix graph to build matrix structure
    std::vector<int> diagonal_rows, diagonal_cols, off_diagonal_rows,
            off_diagonal_cols;

    std::vector < std::vector<bool> > coupling_vars ( 1, std::vector<bool>( 1, true ) );

    InitStructure ( space_, &diagonal_rows, &diagonal_cols,
                    &off_diagonal_rows, &off_diagonal_cols, &coupling_vars );
    couplings_.InitializeCouplings ( off_diagonal_rows, off_diagonal_cols );

    // Initialize matrices and vectors
    if ( init_platform_ )
    {
        initialize_platform ( );
    }
    matrix_.InitStructure ( vec2ptr ( diagonal_rows ),
                            vec2ptr ( diagonal_cols ),
                            diagonal_rows.size ( ),
                            vec2ptr ( off_diagonal_rows ),
                            vec2ptr ( off_diagonal_cols ),
                            off_diagonal_rows.size ( ) );

    sol_.InitStructure ( );
    sol_.Zeros ( );

    rhs_.InitStructure ( );
    rhs_.Zeros ( );
}

void ConvDiff::initialize_platform ( )
{
    // init platform for matrix
    init_platform_mat<CMatrix>( la_sys_.Platform, matrix_impl_, la_matrix_format_,
            matrix_precond_, comm_, couplings_, &matrix_,
            &dev_matrix_, la_sys_ );
    // init platform for solution and right hand side vectors
    init_platform_vec<CVector>( la_sys_.Platform, vector_impl_, comm_, couplings_,
            &sol_, &dev_sol_, la_sys_ );
    init_platform_vec<CVector>( la_sys_.Platform, vector_impl_, comm_, couplings_,
            &rhs_, &dev_rhs_, la_sys_ );
    init_platform_ = false;
}

void ConvDiff::prepare_linear_solver ( )
{
    TimingScope tscope ( "Prepare linear solver." );
    LinearSolverFactory<LAD> LinSolFact;
    linear_solver_ = LinSolFact.Get (
                                      param_["LinearSolver"]["Name"].get<std::string>( ) )->
            params ( param_["LinearSolver"] );

    use_ilupp_ = false;
#ifdef WITH_ILUPP
    // prepare ILU++ preconditioner
    use_ilupp_ = param_["LinearSolver"]["UseILUPP"].get<bool>( );
    if ( use_ilupp_ )
    {
        TimingScope tscope ( "Set up ilu++ preconditioner" );
        ilupp_.InitParameter ( param_["ILUPP"]["PreprocessingType"].get<int>( ),
                               param_["ILUPP"]["PreconditionerNumber"].get<int>( ),
                               param_["ILUPP"]["MaxMultilevels"].get<int>( ),
                               param_["ILUPP"]["MemFactor"].get<double>( ),
                               param_["ILUPP"]["PivotThreshold"].get<double>( ),
                               param_["ILUPP"]["MinPivot"].get<double>( ) );
        ilupp_.SetupOperator ( matrix_ );
        linear_solver_->SetupPreconditioner ( ilupp_ );
    }
#endif

    if ( matrix_precond_ != NOPRECOND && !use_ilupp_ )
    {
        set_up_preconditioner ( );
    }
}

void ConvDiff::set_up_preconditioner ( )
{
    TimingScope tscope ( "Set up preconditioner" );
    if ( matrix_precond_ == JACOBI )
    {
        PreconditionerBlockJacobiStand<LAD> *precond;
        precond = new PreconditionerBlockJacobiStand<LAD>;

        precond->Init_Jacobi ( *dev_sol_ );

        precond_ = precond;
    }
    else
    {
        PreconditionerMultiColoring<LAD> *precond;
        precond = new PreconditionerMultiColoring<LAD>;
        if ( matrix_precond_ == GAUSS_SEIDEL )
        {
            precond->Init_GaussSeidel ( );
        }
        else if ( matrix_precond_ == SGAUSS_SEIDEL )
        {
            precond->Init_SymmetricGaussSeidel ( );
        }
        else if ( matrix_precond_ == ILU )
        {
            precond->Init_ILU ( 0 );
        }
        else if ( matrix_precond_ == ILU2 )
        {
            int param1 = param_["LinearAlgebra"]["ILU2Param1"].get<int>( );
            int param2 = param_["LinearAlgebra"]["ILU2Param2"].get<int>( );
            precond->Init_ILU ( param1, param2 );
        }
        precond->Preprocess ( *dev_matrix_, *dev_sol_, &( space_.dof ( ) ) );
        // sync
        MPI_Barrier ( comm_ );

        prepare_lin_alg_structures ( );

        prepare_bc ( );

        assemble_system ( );

        precond_ = precond;
    }
    precond_->SetupOperator ( *dev_matrix_ );
    precond_->Build ( );
    precond_->Print ( );

    linear_solver_->SetupPreconditioner ( *precond_ );
}

void ConvDiff::assemble_system ( )
{
    TimingScope tscope ( "Assemble system." );
    ConvDiffAssembler local_asm ( beta_, nu_ );
    StandardGlobalAssembler<double> global_asm;

    global_asm.assemble_matrix ( space_, local_asm, matrix_ );
    global_asm.assemble_vector ( space_, local_asm, rhs_ );

    // treatment of dirichlet boudary dofs
    if ( !dirichlet_dofs_.empty ( ) )
    {
        matrix_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                                   static_cast < Scalar > ( 1. ) );
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                         vec2ptr ( dirichlet_values_ ) );
        sol_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                         vec2ptr ( dirichlet_values_ ) );
    }

    rhs_.UpdateCouplings ( );
    sol_.UpdateCouplings ( );

    dev_matrix_->CopyStructureFrom ( matrix_ );
    dev_sol_->CopyStructureFrom ( sol_ );
    dev_rhs_->CopyStructureFrom ( rhs_ );

    dev_sol_->CopyFrom ( sol_ );
    dev_rhs_->CopyFrom ( rhs_ );
    dev_matrix_->CopyFrom ( matrix_ );
}

// solve linear system

void ConvDiff::solve_system ( )
{
    linear_solver_->SetupOperator ( *dev_matrix_ );
    // sync
    MPI_Barrier ( comm_ );
    {
        TimingScope tscope ( "Solve system." );

        linear_solver_->Solve ( *dev_rhs_, dev_sol_ );
    }
    // sync
    MPI_Barrier ( comm_ );

    sol_.CopyFrom ( *dev_sol_ );
}

void ConvDiff::visualize_solution ( )
{
    TimingScope tscope ( "Visualization." );
    // Setup visualization object.
    int num_intervals = 2;
    ParallelCellVisualization<double> visu ( space_, num_intervals, comm_, MASTER_RANK );

    std::stringstream input;
    input << num_partitions_ << "_solution_" << refinement_level_;

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

    sol_.UpdateCouplings ( );
    visu.visualize ( EvalFeFunction<LAD>( space_, sol_ ), "u" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );
    visu.write ( input.str ( ) );
}
