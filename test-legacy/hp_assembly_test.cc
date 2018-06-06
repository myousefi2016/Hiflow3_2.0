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

#include "hp_assembly_test.h"

using mesh::MeshPtr;
using mesh::SharedVertexTable;

const int MASTER_RANK = 0;

const int DEBUG_LEVEL = 2;

int deg0, deg1;

typedef std::vector<double> Coord;

int main ( int argc, char *argv[] )
{

    MPI_Init ( &argc, &argv );
    std::ofstream info_log ( "info_log" );

    LogKeeper::get_log ( "info" ).set_target ( &info_log );
    //    LogKeeper::get_log("info").set_target(&std::cerr);

    std::ofstream debug_log ( "debug_log" );

    LogKeeper::get_log ( "debug" ).set_target ( &debug_log );

    assert ( argc >= 3 );
    deg0 = atoi ( argv[1] );
    deg1 = atoi ( argv[2] );

    LaplaceApp app;

    // read initial mesh
    app.read_mesh ( std::string ( DATADIR ) + std::string ( "two_quads_2d_fix.inp" ) );
    //    app.distribute_mesh();

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
        const double L2_err = app.compute_error ( );
        if ( app.get_rank ( ) == MASTER_RANK )
        {
            std::cout << "Error at refinement level " << app.refinement_level ( )
                    << " = " << L2_err << "\n";

        }

        app.write_visualization ( );

        //        app.debug_output();

        // refine mesh
        app.refine_mesh ( );
    }

    MPI_Finalize ( );
    return 0;
}

//////////////// Laplace application implementation ////////////////

LaplaceApp::LaplaceApp ( )
: comm_ ( MPI_COMM_WORLD ),
refinement_level_ ( 1 ),
fe_degree_ ( std::vector<int>( 1, 1 ) ),
rank_ ( -1 ),
num_partitions_ ( -1 )
{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    // setup linear algebra platform
    la_sys_.Platform = APP_PLATFORM; // NB This decision could be made at runtime
    init_platform ( la_sys_ );

    solver_.InitControl ( 1000, 1.e-14, 1.e-12, 1.e6 );
    solver_.InitParameter ( "NoPreconditioning" ); // for CG
    //    solver_.InitParameter(30, "NoPreconditioning"); // for GMRES

}

LaplaceApp::~LaplaceApp ( )
{
    stop_platform ( la_sys_ );
}

void LaplaceApp::read_mesh ( const std::string& filename )
{
    mesh_ = read_mesh_from_file ( filename, DIM, DIM, 0 );

    // initial refinement -- make sure we have one cell on each partition
    int l = 0;
    while ( l < refinement_level_ )
    {
        mesh_ = mesh_->refine ( );
        ++l;
    }
    old_mesh_ = mesh_;
}

void LaplaceApp::prepare_system ( )
{

    //    space_.Init(fe_degree_, *mesh_);
    std::vector< std::vector<int> > degrees ( 1 );

    //    degrees[0].push_back(deg0);
    //    degrees[0].push_back(deg1);

    //    space_.Init(2, *mesh_);

    degrees[0].resize ( mesh_->num_entities ( DIM ), deg0 );
    for ( int i = 0, end = mesh_->num_entities ( DIM ); i < end; ++i )
    {
        if ( i % 2 == 1 )
        {
            degrees[0][i] = deg1;
        }
    }

    space_.Init_p ( degrees, *mesh_ );

    mesh::AttributePtr degree_attr ( new mesh::IntAttribute ( degrees[0] ) );
    mesh_->add_attribute ( "degree", DIM, degree_attr );

    //    LOG_INFO("dof", "Interface patterns = ");
    //    space_.dof().print_interface_patterns();
    //    LOG_INFO("dof", "Dof numbering = ");
    //    space_.dof().print_numer();
    //    LOG_INFO("dof", "Dof interpolation = " << space_.dof().dof_interpolation());
    permute_constrained_dofs_to_end ( space_.dof ( ) );
#if 1
    LOG_INFO ( "dof", "Reducing ..." );
    space_.dof ( ).dof_interpolation ( ).reduce ( );

    //    LOG_INFO("dof", "Interface patterns = ");
    //    space_.dof().print_interface_patterns();
    //    LOG_INFO("dof", "Dof numbering = ");
    //    space_.dof().print_numer();
    //    LOG_INFO("dof", "Dof interpolation = " << space_.dof().dof_interpolation());
#endif
    // setup couplings object
    couplings_.Init ( comm_, space_.dof ( ) );

    prepare_linear_algebra ( );
}

struct DirichletZero
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        return std::vector<double>( coords_on_face.size ( ), 0. );
    }
};

struct ExactDirichlet
{

    std::vector<double> evaluate ( const mesh::Entity& face,
                                   const std::vector<Coord>& coords_on_face ) const
    {
        std::vector<double> values ( coords_on_face.size ( ), -25. );
        LaplaceApp::ExactSol sol;
        for ( int i = 0, end = values.size ( ); i < end; ++i )
        {
            Vec<LaplaceApp::DIM, double> pt;
            pt[0] = coords_on_face.at ( i )[0];
            pt[1] = coords_on_face.at ( i )[1];
            values.at ( i ) = sol ( pt );
        }
        return values;
    }
};

void LaplaceApp::prepare_bc ( )
{
    dirichlet_dofs_.clear ( );
    dirichlet_values_.clear ( );

#if 1
    DirichletZero dir_zero;

    compute_dirichlet_dofs_and_values ( dir_zero, space_, 0,
                                        dirichlet_dofs_, dirichlet_values_ );
#endif
#if 0
    ExactDirichlet dir_exact;

    compute_dirichlet_dofs_and_values ( dir_exact, space_, 0,
                                        dirichlet_dofs_, dirichlet_values_ );
#endif
    std::cout << "Dirichlet values = " << string_from_range ( dirichlet_values_.begin ( ), dirichlet_values_.end ( ) ) << "\n";
}

void LaplaceApp::prepare_local_refinements ( )
{

    // save pointer to mesh before local refinement
    old_mesh_ = mesh_;

    // refine even-numbered cells
    const int num_cells = mesh_->num_entities ( DIM );
    std::vector<int> refinements ( num_cells, -10 );

    for ( int i = 0; i < num_cells; i += 2 )
    {
        refinements.at ( i ) = 0;
    }

    mesh_ = mesh_->refine ( refinements );
}

void LaplaceApp::prepare_linear_algebra ( )
{
    prepare_bc ( );

    // only CPU support for the moment
    assert ( la_sys_.Platform == CPU );

    // if another matrix format is used, setup looks different
    assert ( APP_MATRIX_FORMAT == CSR );

    // compute matrix graph and use it to build the structure of the matrix
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure ( space_, sparsity );
    couplings_.InitializeCouplings ( sparsity.off_diagonal_rows, sparsity.off_diagonal_cols );

    matrix_.Init ( comm_, couplings_, la_sys_.Platform,
                   APP_LINALG_IMPLEMENTATION, APP_MATRIX_FORMAT );
    x_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );
    rhs_.Init ( comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION );

    matrix_.InitStructure ( vec2ptr ( sparsity.diagonal_rows ),
                            vec2ptr ( sparsity.diagonal_cols ),
                            sparsity.diagonal_rows.size ( ),
                            vec2ptr ( sparsity.off_diagonal_rows ),
                            vec2ptr ( sparsity.off_diagonal_cols ),
                            sparsity.off_diagonal_rows.size ( ) );
    x_.InitStructure ( );
    rhs_.InitStructure ( );

    x_.Zeros ( );
    matrix_.Zeros ( );
    rhs_.Zeros ( );
}

void LaplaceApp::assemble_system ( )
{

    {
        // clear assembly checker file
        std::ofstream octave ( "check_assembly.m" );
        octave.close ( );
    }

    LocalLaplaceAssembler local_asm;
    global_asm_.assemble_matrix ( space_, local_asm, matrix_ );
    global_asm_.assemble_vector ( space_, local_asm, rhs_ );

    // extract RHS vector
    std::vector<int> ind ( rhs_.size_global ( ) );
    std::vector<LAD::DataType> rhs ( rhs_.size_global ( ) );
    for ( int i = 0; i < static_cast < int > ( ind.size ( ) ); ++i )
    {
        ind[i] = i;
    }
    rhs_.GetValues ( vec2ptr ( ind ), ind.size ( ), vec2ptr ( rhs ) );
    LOG_INFO ( "rhs", "Global rhs vector = " << string_from_range ( rhs.begin ( ), rhs.end ( ) ) );

    // write global matrix to file
    matrix_.diagonal ( ).WriteFile ( "matrix.mat" );
    LOG_INFO ( "matrix", "Global matrix written to matrix.mat" );

    std::ofstream octave ( "check_assembly.m", std::ios_base::app );
    octave.precision ( 16 );

    // projection operator into octave
    octave << "% ==== Projection operator for constraints ====\n"
            << "P = eye(" << rhs_.size_global ( ) << ");\n";

    typedef std::vector< std::pair<int, double> >::const_iterator DependencyIterator;

    // modify P for all constrained dofs
    const DofInterpolation& interp = space_.dof ( ).dof_interpolation ( );
    for ( DofInterpolation::const_iterator it = interp.begin ( ); it != interp.end ( ); ++it )
    {

        // dof = index of constrained dof (+1 for octave numbering) .
        // Remove entry corresponding to dof
        const int dof = it->first + 1;
        octave << "P(" << dof << ", " << dof << ") = 0; ";

        // Modify entries corresponding to ddof -> dofs that the constrained dofs depends on.
        for ( DependencyIterator d_it = it->second.begin ( ); d_it != it->second.end ( ); ++d_it )
        {
            const int ddof = d_it->first + 1;
            octave << "P(" << dof << ", " << ddof << ") = " << d_it->second << "; ";
        }
        octave << "\n";
    }
    octave << "\n";

    // Import computed results into octave
    octave << "% ==== Computed results ====\n";
    octave << "A_computed = mmread('matrix.mat');\n"
            << "A_computed = full(A_computed);\n"
            << "b_computed = [" << precise_string_from_range ( rhs.begin ( ), rhs.end ( ) ) << "]';\n\n";

    // Check computed results
    octave << "% === Check ==== \n";

    // create constrained correction to P^T*A*P
    octave << "c = [";
    for ( DofInterpolation::const_iterator it = interp.begin ( ); it != interp.end ( ); ++it )
    {
        octave << it->first + 1 << " ";
    }
    octave << "]; % constrained dofs\n"
            << "C = zeros(" << rhs_.size_global ( ) << ", 1); C(c) = 1; C = diag(C);\n";

    // output results
    octave << "A_check = P'*A*P + C;\n"
            << "b_check = P'*b;\n"
            << "A_computed;\n"
            << "b_computed;\n"
            << "max_err_A = max(max(A_computed - A_check))\n"
            << "max_err_b = max(b_computed - b_check)\n";

    // create automatic check
    octave << "if and(abs(max_err_A) < 1e-12, abs(max_err_b) < 1e-12)\n"
            << "\tdisp \'Passed\'\n"
            << "\texit(0)\n"
            << "else\n"
            << "\tdisp 'Failed'\n"
            << "\tmax_err_A\n"
            << "\tmax_err_b\n"
            << "\texit(-1)\n"
            << "end\n";

    octave.close ( );

    if ( !dirichlet_dofs_.empty ( ) )
    {
        matrix_.diagonalize_rows ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), 1. );
        x_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
        rhs_.SetValues ( vec2ptr ( dirichlet_dofs_ ), dirichlet_dofs_.size ( ), vec2ptr ( dirichlet_values_ ) );
    }
}

void LaplaceApp::solve_system ( )
{
    // solve linear system
    solver_.SetupOperator ( matrix_ );
    solver_.Solve ( rhs_, &x_ );
    interpolate_constrained_vector<double>( space_, x_ );
}

double LaplaceApp::compute_error ( )
{

    // create error integrator with solution as parameter
    x_.UpdateCouplings ( );
    L2ErrorIntegrator<LaplaceApp::ExactSol> error_integrator ( x_ );

    // integrate
    double L2err = 0.;
    global_asm_.integrate_scalar ( space_, error_integrator, L2err );

    double global_L2err;
    MPI_Reduce ( &L2err, &global_L2err, 1, MPI_DOUBLE, MPI_SUM, MASTER_RANK, comm_ );

    std::vector<double> err_vec;
    global_asm_.assemble_scalar ( space_, error_integrator, err_vec );
    LOG_INFO ( "error", "Errors = " << string_from_range ( err_vec.begin ( ), err_vec.end ( ) ) );
    mesh::AttributePtr err_attr ( new mesh::DoubleAttribute ( err_vec ) );
    mesh_->add_attribute ( "error", DIM, err_attr );

    if ( rank_ == MASTER_RANK )
    {
        if ( !err_.empty ( ) )
        {
            std::cout << "Error quotient = " << err_.back ( ) / std::sqrt ( global_L2err ) << "\n";
        }
        err_.push_back ( std::sqrt ( global_L2err ) );
        return err_.back ( );
    }
    else
    {
        return 0.;
    }

}

void LaplaceApp::refine_mesh ( )
{
    mesh_ = old_mesh_->refine ( );
    ++refinement_level_;
    prepare_local_refinements ( );
}

void LaplaceApp::write_visualization ( )
{

    x_.UpdateCouplings ( );

    std::stringstream input;
    input << "hp_assembly_" << refinement_level_;

    ParallelCellVisualization<double> visu ( space_, 1, comm_, MASTER_RANK );
    visu.visualize ( EvalFeFunction<LAD>( space_, x_ ), "u" );
    visu.write ( input.str ( ) );

}

bool LaplaceApp::is_done ( ) const
{
    MPI_Bcast ( &refinement_level_, 1, MPI_INT, MASTER_RANK, comm_ );
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
