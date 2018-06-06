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

// Yet another Laplace solver
// author: Staffan Ronnas

#include "yalaplace.h"

#include <fenv.h>
#include <fstream>
#include <signal.h>

using mesh::MeshPtr;
using mesh::SharedVertexTable;

const int FE_MASK = FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID;

const int MASTER_RANK = 0;

void FPExceptionHandler ( int signalNumber )
{
    std::cerr << "FP exception caught:\n";
    int result = fetestexcept ( FE_MASK );
    std::cerr << "mask = " << FE_MASK << "\n";
    std::cerr << "type = " << result << "\n";
    if ( result & FE_DIVBYZERO ) std::cerr << "divide by zero\n";
    if ( result & FE_UNDERFLOW ) std::cerr << "underflow\n";
    if ( result & FE_OVERFLOW ) std::cerr << "overflow\n";
    if ( result & FE_INVALID ) std::cerr << "invalid\n";
    if ( result & FE_INEXACT ) std::cerr << "inexact\n";
    std::cerr << "\n";
    if ( result & FE_MASK )
    {
        signal ( SIGFPE, SIG_DFL );
        raise ( SIGFPE );
    }
}

void facet_quadrature_selection ( const Element<double>& elem, Quadrature<double>& quadrature )
{
    const FEType<double>::FiniteElement fe_id = elem.get_fe_type ( 0 )->get_my_id ( );

    int fe_deg = 0;

    // compute maxmimum FE degree for all variables
    for ( int v = 0, end_v = elem.get_num_variables ( ); v != end_v; ++v )
    {
        fe_deg = std::max ( fe_deg, elem.get_fe_type ( v )->get_fe_deg ( ) );
    }

    switch ( fe_id )
    {
        case FEType<double>::LAGRANGE_TRI:
        case FEType<double>::LAGRANGE_QUAD:
            quadrature.set_quadrature ( "GaussLine", 3 );
            break;
        case FEType<double>::LAGRANGE_TET:
            quadrature.set_quadrature ( "GaussTriangle", 4 );
            break;
        case FEType<double>::LAGRANGE_HEX:
            quadrature.set_quadrature ( "GaussQuadrilateral", 16 );
            break;
        default:
            assert ( 0 );
    };
}

int main ( int argc, char *argv[] )
{
    MPI_Init ( &argc, &argv );

    std::ofstream info_log ( "yalaplace_info_log" );
    LogKeeper::get_log ( "info" ).set_target ( &info_log );
    std::ofstream debug_log ( "yalaplace_debug_log" );
    LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

    //#pragma STDC FENV_ACCESS ON
    // clear FP exceptions
    feclearexcept ( FE_ALL_EXCEPT );
    // enable FP exceptions
    //    feenableexcept(FE_MASK);
    // set signal handler for FP exceptions
    //    signal(SIGFPE, &FPExceptionHandler);

    if ( argc < 2 )
    {
        std::cerr << "Usage: yalaplace <degree>\n";
        return -1;
    }

    int deg = atoi ( argv[1] );
    if ( deg < 1 )
    {
        std::cerr << "FE degree must be greater than 1.\n";
        return -1;
    }

    LaplaceApp app ( deg );

    // read initial mesh
#if DIMENSION == 3
    app.read_mesh ( std::string ( DATADIR ) + std::string ( "unit_cube.inp" ) );
#elif DIMENSION ==2
    app.read_mesh ( std::string ( DATADIR ) + std::string ( "unit_square.inp" ) );
#endif
    app.distribute_mesh ( );

    // application main loop
    while ( !app.is_done ( ) )
    {
        // prepare system
        app.prepare_system ( );

        // assemble system
        app.assemble_system ( );

        // solve system
        app.solve_system ( );

        // compute error
        const double L2_err = app.compute_error ( LaplaceApp::L2error );
        if ( app.get_rank ( ) == MASTER_RANK )
        {
            std::cout << "L2 error at refinement level " << app.refinement_level ( )
                    << " = " << L2_err << "\n";
        }
        const double H1_err = app.compute_error ( LaplaceApp::H1error );

        if ( app.get_rank ( ) == MASTER_RANK )
        {
            std::cout << "H1 error at refinement level " << app.refinement_level ( )
                    << " = " << H1_err << "\n";

        }

        app.write_visualization ( );

        // refine mesh
        app.refine_mesh ( );
    }

    MPI_Finalize ( );
    return 0;
}

//////////////// Laplace application implementation ////////////////

LaplaceApp::LaplaceApp ( int degree )
: comm_ ( MPI_COMM_WORLD ),
refinement_level_ ( 3 ),
fe_degree_ ( std::vector<int>( 1, degree ) ),
rank_ ( -1 ),
num_partitions_ ( -1 )
{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    // setup linear algebra platform
    la_sys_.Platform = APP_PLATFORM; // NB This decision could be made at runtime
    init_platform ( la_sys_ );

    solver_.InitControl ( 5000, 1.e-14, 1.e-12, 1.e6 );
    solver_.InitParameter ( "NoPreconditioning" ); // for CG
    //    solver_.InitParameter(50, "NoPreconditioning"); // for GMRES
}

LaplaceApp::~LaplaceApp ( )
{
    stop_platform ( la_sys_ );
}

void LaplaceApp::read_mesh ( const std::string& filename )
{
    if ( rank_ == MASTER_RANK )
    {
        mesh::MeshDbViewBuilder builder ( DIM, DIM );

        // TODO: check filename extension to choose reader (factory function)
        ScopedPtr<mesh::Reader>::Type reader ( new mesh::UcdReader ( &builder ) );
        reader->read ( filename.c_str ( ), master_mesh_ );

        // initial refinement -- make sure we have one cell on each partition
        //        while (master_mesh_->num_entities(DIM) < num_partitions_) {
        for ( int r = 0; r < refinement_level_; ++r )
        {
            master_mesh_ = master_mesh_->refine ( );
            //            ++refinement_level_;
        }
    }
}

void LaplaceApp::distribute_mesh ( )
{
    MeshPtr local_mesh = partition_and_distribute ( master_mesh_,
                                                    MASTER_RANK, comm_ );
    assert ( local_mesh != 0 );
    SharedVertexTable shared_verts;
    mesh_ = compute_ghost_cells ( *local_mesh, comm_, shared_verts );
}

void LaplaceApp::prepare_system ( )
{
    space_.Init ( fe_degree_, *mesh_ );

    // setup couplings object
    couplings_.Init ( comm_, space_.dof ( ) );

    prepare_bc ( );
    // Do we really need to do this after every refinement --
    // none of the parameters seem to change...
    // Maybe changing dof_partition_ modifies couplings_?
    prepare_linear_algebra ( );
}

void LaplaceApp::prepare_linear_algebra ( )
{
    // only CPU support for the moment
    assert ( la_sys_.Platform == CPU );

    // if another matrix format is used, setup looks different
    assert ( APP_MATRIX_FORMAT == CSR );

    matrix_.Init ( comm_, couplings_, la_sys_.Platform,
                   APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    x_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    rhs_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );

    // compute matrix graph and use it to build the structure of the matrix
    std::vector<int> diagonal_rows, diagonal_cols, off_diagonal_rows, off_diagonal_cols;

    std::vector < std::vector<bool> > coupling_vars ( 1, std::vector<bool>( 1, true ) );

    InitStructure ( space_, &diagonal_rows, &diagonal_cols, &off_diagonal_rows, &off_diagonal_cols, &coupling_vars );
    couplings_.InitializeCouplings ( off_diagonal_rows, off_diagonal_cols );

    matrix_.InitStructure ( vec2ptr ( diagonal_rows ), vec2ptr ( diagonal_cols ), diagonal_rows.size ( ),
                            vec2ptr ( off_diagonal_rows ), vec2ptr ( off_diagonal_cols ), off_diagonal_rows.size ( ) );
    x_.InitStructure ( );
    rhs_.InitStructure ( );

    x_.Zeros ( );
    matrix_.Zeros ( );
    rhs_.Zeros ( );

#ifndef WEAK_BOUNDARY_CONDITIONS
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // correct solution with dirichlet BC
        x_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }
#endif //NOT WEAK_BOUNDARY_CONDITIONS

}

void LaplaceApp::assemble_system ( )
{

    LocalLaplaceAssembler local_asm ( dirichlet_dofs_, dirichlet_values_ );

    StandardGlobalAssembler<double> global_asm;
    global_asm.assemble_matrix ( space_, local_asm, matrix_ );
    //    matrix_.diagonal().WriteFile("matrix.mat");
    global_asm.assemble_vector ( space_, local_asm, rhs_ );

#ifdef WEAK_BOUNDARY_CONDITIONS
    global_asm.set_quadrature_selection_function ( &facet_quadrature_selection );
    global_asm.should_reset_assembly_target ( false );
    global_asm.assemble_matrix_boundary ( space_, local_asm, matrix_ );
    global_asm.assemble_vector_boundary ( space_, local_asm, rhs_ );
#else
    if ( !dirichlet_dofs_.empty ( ) )
    {
        // set rows corresponding to boundary dofs to 0
        // and the diagonal entry to 1 (coupled_matrix.h)
        matrix_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ),
                                   static_cast < LAD::DataType > ( 1. ) );
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }

#endif //NOT WEAK_BOUNDARY_CONDITIONS

    //    matrix_.diagonal().WriteFile("laplace.mtx");

}

void LaplaceApp::solve_system ( )
{
    // solve linear system
    solver_.SetupOperator ( matrix_ );
    solver_.Solve ( rhs_, &x_ );
}

// Homogenous Dirichlet boundary conditions

struct Dirichlet_0_BC
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        //homogenous Dirichlet boundary
        return std::vector<double>( coords_on_face.size ( ), 0. );
    }
};

void LaplaceApp::prepare_bc ( )
{
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

    Dirichlet_0_BC bc;
    compute_dirichlet_dofs_and_values ( bc, space_, 0,
                                        dirichlet_dofs_, dirichlet_values_ );
}

double LaplaceApp::compute_error ( ErrorNorm norm )
{

    // create error integrator with solution as parameter

    // explicit update of the ghost DoFs, since these are needed in
    // integral evaluation
    x_.UpdateCouplings ( );

    std::vector<double> err ( 0, 0.0 );

    if ( norm == L2error )
    {
        L2ErrorIntegrator<LaplaceApp::ExactSol> error_integrator ( x_ );
        StandardGlobalAssembler<double> global_asm;
        global_asm.assemble_scalar ( space_, error_integrator, err );
    }
    else if ( norm == H1error )
    {
        H1ErrorIntegrator<LaplaceApp::ExactSol> error_integrator ( x_ );
        StandardGlobalAssembler<double> global_asm;
        global_asm.assemble_scalar ( space_, error_integrator, err );
    }
    else
    {
        std::cerr << "Unknow error norm!\n";
        exit ( -1 );
    }
    double total_err = std::accumulate ( err.begin ( ), err.end ( ), 0. );
    double global_err;
    MPI_Reduce ( &total_err, &global_err, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK, comm_ );

    std::vector<double>& err_vec = ( norm == L2error ? L2_err_ : H1_err_ );

    if ( rank_ == MASTER_RANK )
    {
        if ( !err_vec.empty ( ) )
        {
            std::cout << ( norm == L2error ? "L2" : "H1" )
                    << " error quotient = " << err_vec.back ( ) / std::sqrt ( global_err ) << "\n";
        }
        err_vec.push_back ( std::sqrt ( global_err ) );
        return err_vec.back ( );
    }
    else
    {
        return 0.;
    }
}

void LaplaceApp::refine_mesh ( )
{
    // refine the mesh on the master
    if ( rank_ == MASTER_RANK )
    {
        master_mesh_ = master_mesh_->refine ( );
        ++refinement_level_;
        std::cout << "Refinement level " << refinement_level_ << "\n";
    }
    // redistribute mesh
    distribute_mesh ( );
}

void LaplaceApp::write_visualization ( )
{

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

    // Visualize the solution vector and all attributes
    x_.UpdateCouplings ( );
    visu.visualize ( EvalFeFunction<LAD>( space_, x_ ), "u" );

    visu.visualize_cell_data ( material_number, "Material Id" );
    visu.visualize_cell_data ( remote_index, "_remote_index_" );
    visu.visualize_cell_data ( sub_domain, "_sub_domain_" );

    visu.write ( input.str ( ) );
}

bool LaplaceApp::is_done ( ) const
{
    MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
    std::cout << "Ref level = " << refinement_level_ << ", max ref level = " << MAX_REFINEMENT_LEVEL << "\n";

    return refinement_level_ >= MAX_REFINEMENT_LEVEL;
}

void LaplaceApp::debug_output ( )
{
    // NB: this will probably not work in parallel!
    // extract solution vector
    std::vector<int> ind ( x_.size_global ( ) );
    std::vector<LAD::DataType> visu_vec ( x_.size_global ( ) );
    for ( int i = 0; i < static_cast < int > ( ind.size ( ) ); ++i )
    {
        ind[i] = i;
    }

    const int num_dofs = ind.size ( );
    const int NVAR = 1;
    const int dof_per_var = num_dofs / NVAR;
    for ( int v = 0; v < NVAR; ++v )
    {
        const int dof_begin = v * dof_per_var;
        const int dof_end = ( v + 1 ) * dof_per_var;
        x_.GetValues ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( visu_vec ) );
        std::cerr << "x" << v << " = " <<
                string_from_range ( visu_vec.begin ( ) + dof_begin, visu_vec.begin ( ) + dof_end ) << "\n";

        rhs_.GetValues ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( visu_vec ) );
        std::cerr << "b" << v << " = " <<
                string_from_range ( visu_vec.begin ( ) + dof_begin, visu_vec.begin ( ) + dof_end ) << "\n";
    }

    matrix_.diagonal ( ).WriteFile ( "yalaplace_matrix.mat" );
}
