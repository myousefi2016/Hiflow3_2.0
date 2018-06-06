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

/// \author Volker Lange, Teresa Beck

#include "vorticity.h"

/// Constructor

template<int DIM>
Vorticity<DIM>::Vorticity ( )
: comm_ ( MPI_COMM_WORLD )
{
    initialized = false;
    was_setup = false;
    vort_solved = false;

    //everything is a pointer (and set to NULL) to save memory
    mesh_ = NULL;
    space_ = NULL;

    couplings_ = NULL;
    sparsity_ = NULL;

    global_asm_ = NULL;
    vort_asm_ = NULL;
    stream_asm_ = NULL;

    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );
}

/// Constructor with initialization and setup for immediat use of class    

template<int DIM>
Vorticity<DIM>::Vorticity ( mesh::MeshPtr mesh, int degrees_vars )
: comm_ ( MPI_COMM_WORLD )
{
    MPI_Comm_rank ( comm_, &rank_ );
    MPI_Comm_size ( comm_, &num_partitions_ );

    initialize ( );
    setup ( mesh, degrees_vars );
}

template<int DIM>
Vorticity<DIM>::~Vorticity ( )
{
    delete global_asm_;
    delete vort_asm_;
    delete space_;
    mesh_.reset ( );
    delete sparsity_;
    delete couplings_;

}

/// calculate the vorticity from a given velocity vector

template<int DIM>
void Vorticity<DIM>::calc_vort ( LAD::VectorType& velo, VectorSpace<double> const& velo_space )
{
    assert ( initialized );
    assert ( was_setup );

    // initializiation
    LAD::VectorType u_;

#ifdef USE_HYPRE
    u_.Init ( comm_, *couplings_ );
#else  
    u_.Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ) );
    u_.InitStructure ( );
#endif  
    u_.Zeros ( );

    LAD::VectorType v_;

#ifdef USE_HYPRE
    v_.Init ( comm_, *couplings_ );
#else
    v_.Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ) );
    v_.InitStructure ( );
#endif  
    v_.Zeros ( );

    // HIFLOW BUG
    //   here();
    //   MPI_Barrier(MPI_COMM_WORLD);
    //   for (int i=0; i<num_partitions(); ++i)
    //   {
    //     if (rank() == i)
    //     {
    //       std::cout << std::endl;
    //       std::cout << "domain " << rank() << " dof.num_partitions():           "           << num_partitions()                 << std::endl;
    //       std::cout << "domain " << rank() << " space.get_nb_var():             "           << space_->get_nb_var()             << std::endl;
    //       std::cout << "domain " << rank() << " dof.ndofs_on_sd(  " << i << "):           " << space_->dof().ndofs_on_sd(sd   ) << std::endl;
    //       std::cout << "domain " << rank() << " dof.ndofs_on_sd(0," << i << "):           " << space_->dof().ndofs_on_sd(0, sd) << std::endl;
    //       std::cout << "domain " << rank() << " dof.ndofs global():             "           << space_->dof().ndofs_global()     << std::endl;
    //       std::cout << std::endl;
    //     } 
    //     MPI_Barrier(MPI_COMM_WORLD);
    //   }
    //   MPI_Barrier(MPI_COMM_WORLD);

    // Determine solution vector "u_"

    // -> calculate indices based on velo vector space
    const int sd = rank ( );
    const int u_var = 0; // assuming u is the first quadrant in given space
    const int dofs_on_sd = space_->dof ( ).ndofs_on_sd ( sd ); // for one vel-comp
    const int u_offset = u_var*dofs_on_sd; // on current process, first dof for u-comp in scalar space has index 'u_offset'
    int domain_offset = velo_space.dof ( ).get_my_dof_offset ( ); // on current process, first dof            in velo   space has index 'domain_offset'
    int *indices = new int[dofs_on_sd]; // indices related to u-comp of velo space
    for ( int i = 0; i < dofs_on_sd; ++i )
        indices[i] = domain_offset + u_offset + i;

    // -> get local data from given solution vector "velo"
    double *values;
    values = new double[dofs_on_sd];
    velo.GetValues ( indices, dofs_on_sd, values );

    // -> calculate indices based on scalar vector space
    domain_offset = space_->dof ( ).get_my_dof_offset ( ); // on current process, first dof            in scalar space has index 'domain_offset'
    for ( int i = 0; i < dofs_on_sd; ++i )
        indices[i] = domain_offset + i;

    // -> set local data to new solution vector "u_"
    u_.SetValues ( indices, dofs_on_sd, values );
#ifdef USE_HYPRE
    u_.Update ( );
#else  
    u_.UpdateCouplings ( );
#endif  

    // Determine solution vector "v_"

    // -> calculate indices based on velo vector space
    const int v_var = 1; // assuming u is the first quadrant in given space
    const int v_offset = v_var*dofs_on_sd; // on current process, first dof for u-comp in scalar space has index 'u_offset'
    domain_offset = velo_space.dof ( ).get_my_dof_offset ( ); // on current process, first dof            in velo   space has index 'domain_offset'
    for ( int i = 0; i < dofs_on_sd; ++i )
        indices[i] = domain_offset + v_offset + i;

    // -> get local data from given solution vector "velo"
    velo.GetValues ( indices, dofs_on_sd, values );

    // -> calculate indices based on scalar vector space
    domain_offset = space_->dof ( ).get_my_dof_offset ( ); // on current process, first dof            in scalar space has index 'domain_offset'
    for ( int i = 0; i < dofs_on_sd; ++i )
        indices[i] = domain_offset + i;

    // -> set local data to new solution vector "u_"
    v_.SetValues ( indices, dofs_on_sd, values );
#ifdef USE_HYPRE
    v_.Update ( );
#else  
    v_.UpdateCouplings ( );
#endif 

    // Clean up
    delete values;
    delete indices;

    // Solve linear system to determine vorticity

    vort_asm_->set_from_function ( false );
    vort_asm_->set_velocity_solution ( &u_, &v_ );

    // Assemble matrix and right-hand-side vector for stream problem.    
    global_asm_->assemble_vector ( *space_, *vort_asm_, rhs_vort_ );

    solver_.SetupOperator ( matrix_vort_ );
#ifndef USE_HYPRE  
    ilupp_.SetupOperator ( matrix_vort_ );
#endif    
    solver_.Solve ( rhs_vort_, &sol_vort_ );
    interpolate_constrained_vector ( *space_, sol_vort_ );
#ifdef USE_HYPRE
    sol_vort_.Update ( );
#else  
    sol_vort_.UpdateCouplings ( );
#endif     

    vort_solved = true;
}

/// calculate the vorticity from a given function

template<int DIM>
void Vorticity<DIM>::calc_vort ( double (*f_ )( Vec<DIM, double> ) )
{
    // f_ is a pointer to a function which returns a double and takes a Vec<2>
    assert ( initialized );
    assert ( was_setup );

    vort_asm_->set_from_function ( true );
    vort_asm_->set_function ( f_ );

    // Assemble matrix and right-hand-side vector for stream problem.    
    global_asm_->assemble_vector ( *space_, *vort_asm_, rhs_vort_ );

    solver_.SetupOperator ( matrix_vort_ );
#ifndef USE_HYPRE  
    ilupp_.SetupOperator ( matrix_vort_ );
#endif
    solver_.Solve ( rhs_vort_, &sol_vort_ );
    interpolate_constrained_vector ( *space_, sol_vort_ );
    // Communicate new solution values.
#ifdef USE_HYPRE
    sol_vort_.Update ( );
#else  
    sol_vort_.UpdateCouplings ( );
#endif 

    vort_solved = true;
}

template<int DIM>
void Vorticity<DIM>::periodify_space ( )
{
    assert ( mesh_ != 0 );
    const int tdim = space ( ).mesh ( ).tdim ( );
    // loop over all cells and init the cell transformations
    for ( mesh::EntityIterator it = space ( ).mesh ( ).begin ( tdim );
          it != space ( ).mesh ( ).end ( tdim );
          ++it )
    {
        std::vector<double> coord_vtx;
        it->get_coordinates ( coord_vtx );
        std::vector<MasterSlave> period = mesh_->get_period ( );
        coord_vtx = unperiodify ( coord_vtx, DIM, period );
        space ( ).fe_manager ( ).get_cell_transformation ( it->index ( ) )->reinit ( coord_vtx );
    }
}

/// Setup mesh, solvers, LA

template<int DIM>
void Vorticity<DIM>::setup ( mesh::MeshPtr mesh, int velo_degree )
{
    assert ( initialized );

    setup_solver ( );
    setup_linear_algebra ( );

    mesh_ = mesh;

    space_->Init ( velo_degree, *mesh_ );

    periodify_space ( );

    couplings_->Clear ( );
    couplings_->Init ( comm_, space ( ).dof ( ) );

    global_asm_->compute_sparsity_structure ( *space_, *sparsity_ );

    couplings_->InitializeCouplings ( sparsity_->off_diagonal_rows, sparsity_->off_diagonal_cols );

#ifdef USE_HYPRE
    matrix_vort_ .Init ( comm_, *couplings_ );
    matrix_stream_.Init ( comm_, *couplings_ );
#else
    matrix_vort_ .Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ), la_matrix_format ( ) );
    matrix_stream_.Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ), la_matrix_format ( ) );
#endif  

    matrix_vort_.InitStructure ( vec2ptr ( sparsity_->diagonal_rows ),
                                 vec2ptr ( sparsity_->diagonal_cols ),
                                 sparsity_->diagonal_rows.size ( ),
                                 vec2ptr ( sparsity_->off_diagonal_rows ),
                                 vec2ptr ( sparsity_->off_diagonal_cols ),
                                 sparsity_->off_diagonal_rows.size ( ) );

    matrix_stream_.InitStructure ( vec2ptr ( sparsity_->diagonal_rows ),
                                   vec2ptr ( sparsity_->diagonal_cols ),
                                   sparsity_->diagonal_rows.size ( ),
                                   vec2ptr ( sparsity_->off_diagonal_rows ),
                                   vec2ptr ( sparsity_->off_diagonal_cols ),
                                   sparsity_->off_diagonal_rows.size ( ) );

#ifdef USE_HYPRE
    sol_ .Init ( comm_, *couplings_ );
    sol_vort_ .Init ( comm_, *couplings_ );
    sol_stream_.Init ( comm_, *couplings_ );
    rhs_vort_ .Init ( comm_, *couplings_ );
    rhs_stream_.Init ( comm_, *couplings_ );
#else
    sol_ .Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ) );
    sol_vort_ .Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ) );
    sol_stream_.Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ) );
    rhs_vort_ .Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ) );
    rhs_stream_.Init ( communicator ( ), *couplings_, la_platform ( ), la_implementation ( ) );

    sol_ .InitStructure ( );
    sol_vort_ .InitStructure ( );
    sol_stream_.InitStructure ( );
    rhs_vort_ .InitStructure ( );
    rhs_stream_.InitStructure ( );
#endif

    rhs_vort_ .Zeros ( );
    rhs_stream_.Zeros ( );
    sol_ .Zeros ( );
    sol_vort_ .Zeros ( );
    sol_stream_.Zeros ( );

    global_asm_->assemble_matrix ( *space_, *vort_asm_, matrix_vort_ );
    global_asm_->assemble_matrix ( *space_, *stream_asm_, matrix_stream_ );
#ifndef USE_HYPRE  
    ilupp_.SetupOperator ( matrix_vort_ );
    ilupp_.SetupOperator ( matrix_stream_ );
#endif  
    was_setup = true;
}

/// initialize the pointers with *new*

template<int DIM>
void Vorticity<DIM>::initialize ( )
{
    space_ = new VectorSpace<double>;

    couplings_ = new Couplings<double>;
    sparsity_ = new SparsityStructure;

    global_asm_ = new HpFemAssembler<double>; //StandardGlobalAssembler;
    vort_asm_ = new VortAssembler<DIM>;
    stream_asm_ = new StreamAssembler<DIM>;

    initialized = true;
}

/// delete all the objections needed for calculation

template<int DIM>
void Vorticity<DIM>::free_mem ( )
{
    // delete all objects
    delete space_;

    delete global_asm_;
    delete vort_asm_;
    delete stream_asm_;

    //  keep couplings and sparsity for further use with solution vectors

    // set pointers to NULL  
    space_ = NULL;

    global_asm_ = NULL;
    vort_asm_ = NULL;
    stream_asm_ = NULL;

    matrix_vort_.Zeros ( );
    matrix_stream_.Zeros ( );

    initialized = false;
}

/// setup the LA with given values

template<int DIM>
void Vorticity<DIM>::setup_linear_algebra ( const std::string platform_str, const std::string impl_str, const std::string matrix_str )
{
    assert ( initialized );

    if ( platform_str == "CPU" )
    {
        la_sys_.Platform = CPU;
    }
    else if ( platform_str == "GPU" )
    {
        la_sys_.Platform = GPU;
    }
    init_platform ( la_sys_ );

    if ( impl_str == "Naive" )
    {
        la_impl_ = NAIVE;
    }
    else if ( impl_str == "BLAS" )
    {
        la_impl_ = BLAS;
    }
    else if ( impl_str == "MKL" )
    {
        la_impl_ = MKL;
    }
    else if ( impl_str == "OPENMP" )
    {
        la_impl_ = OPENMP;
    }
    else if ( impl_str == "SCALAR" )
    {
        la_impl_ = SCALAR;
    }
    else if ( impl_str == "SCALAR_TEX" )
    {
        la_impl_ = SCALAR_TEX;
    }

    if ( matrix_str == "CSR" )
    {
        la_matrix_format_ = CSR;
    }
    else if ( matrix_str == "COO" )
    {
        la_matrix_format_ = COO;
    }
}

/// setup the LA with default values

template<int DIM>
void Vorticity<DIM>::setup_linear_algebra ( )
{
    assert ( initialized );

    la_sys_.Platform = CPU;
    init_platform ( la_sys_ );

    la_impl_ = NAIVE;

    la_matrix_format_ = CSR;
}

/// setup the linear solver and preconditioner ILU++ with given values 

template<int DIM>
void Vorticity<DIM>::setup_solver ( const int max_iter, const double abs_tol, const double rel_tol, const double div_tol )
{
    assert ( initialized );

    solver_.InitControl ( max_iter, abs_tol, rel_tol, div_tol );
    solver_.InitParameter ( 70, "RightPreconditioning" );
#ifndef USE_HYPRE  
    ilupp_.InitParameter ( 0, 1010, 20, 0.8, 1.75, 0.01 );
    solver_.SetupPreconditioner ( ilupp_ );
#endif  
}

/// setup the linear solver and preconditioner ILU++ with default values

template<int DIM>
void Vorticity<DIM>::setup_solver ( )
{
    assert ( initialized );

    solver_.InitControl ( 200, 1.e-9, 1.e-9, 1.e6 ); //max_iter, abs_tol, rel_tol, div_tol
    solver_.InitParameter ( 70, "RightPreconditioning" );
#ifndef USE_HYPRE
    ilupp_.InitParameter ( 0, 1010, 20, 0.8, 1.75, 0.01 );
    solver_.SetupPreconditioner ( ilupp_ );
#endif  
}

/// calculate the velocity from a given vorticity-function

template<int DIM>
void Vorticity<DIM>::calc_velo ( double (*f_ )( Vec<DIM, double> ) )
{
    // f_ is a pointer to a function which returns a double and takes a Vec<2>
    assert ( initialized );
    assert ( was_setup );

    stream_asm_->set_from_function ( true );
    stream_asm_->set_function ( f_ );

    // Assemble matrix and right-hand-side vector for stream problem.    
    global_asm_->assemble_vector ( *space_, *stream_asm_, rhs_stream_ );

    solver_.SetupOperator ( matrix_stream_ );

    solver_.Solve ( rhs_stream_, &sol_stream_ );
    interpolate_constrained_vector ( *space_, sol_stream_ );
    // Communicate new solution values.
#ifdef USE_HYPRE
    sol_stream_.Update ( );
#else  
    sol_stream_.UpdateCouplings ( );
#endif 
    streamsolved = true;
}

/// calculate the stream from a given vorticity-vector

template<int DIM>
void Vorticity<DIM>::calc_stream ( LAD::VectorType& vort )
{
    // f_ is a pointer to a function which returns a double and takes a Vec<2>
    assert ( initialized );
    assert ( was_setup );

    stream_asm_->set_from_function ( false );
    stream_asm_->set_vorticity_solution ( &vort );

    // Assemble matrix and right-hand-side vector for stream problem.    
    global_asm_->assemble_vector ( *space_, *stream_asm_, rhs_stream_ );

    solver_.SetupOperator ( matrix_stream_ );

    solver_.Solve ( rhs_stream_, &sol_stream_ );
    interpolate_constrained_vector ( *space_, sol_stream_ );
    // Communicate new solution values.
#ifdef USE_HYPRE
    sol_stream_.Update ( );
#else  
    sol_stream_.UpdateCouplings ( );
#endif   
    streamsolved = true;
}

// template instanciation
template class Vorticity<2>;
template class Vorticity<3>;
