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

#include "conv_diff.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>

/// ***************************************************************************
/// POSTPROCESSING
/// ***************************************************************************

void MetFlowConvDiffApp::compute_quant_vector ( std::string filename, double time )
{

}

/// Compute some flow characteristics

void MetFlowConvDiffApp::compute_flow_quantities ( std::string filename, LAD::VectorType* sol, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Compute characteristics of TEHD flow " << std::endl;
    Timer timer;
    timer.start ( );

    // Set soultion vector w.r.t. which quantities should be computed
    T_local_asm_.set_solP ( *sol );

    //    compute_quant_scalar();
    compute_quant_vector ( filename, time );

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) )
        std::cout << "  took " << timer.get_duration ( ) << " sec" << std::endl;

    // Reset solution vector
    T_local_asm_.set_solP ( *this->solP_ );
}

/// Compute some flow characteristics

void MetFlowConvDiffApp::compute_solution_norm ( std::string filename, LAD::VectorType* sol, LAD::VectorType* sol_prev, double time )
{
    if ( rank ( ) == master_rank ( ) ) std::cout << "> Compute norm of solution vector " << std::endl;

    Timer timer;
    timer.start ( );

    // Set soultion vector w.r.t. which quantities should be computed
    T_local_asm_prod_.set_left_vector ( *sol );
    T_local_asm_prod_.set_right_vector ( *sol );
    std::vector<bool> no_L2, L2_temp;
    std::vector<bool> no_H1, H1_temp;
    std::vector<bool> no_H2;

    no_L2.resize ( this->num_vars_, false );
    no_H1.resize ( this->num_vars_, false );
    no_H2.resize ( this->num_vars_, false );

    // temp
    L2_temp.resize ( this->num_vars_, false );
    H1_temp.resize ( this->num_vars_, false );
    L2_temp[t_var_] = true;
    H1_temp[t_var_] = true;

    // ********************************************************************************
    // norms
    // temperature
    // L2
    std::vector<double> cell_vals;
    T_local_asm_prod_.set_mode ( L2_temp, no_H1, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( T_local_asm_prod_ ), cell_vals );

    double local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    double global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double temp_l2 = std::sqrt ( global_val );

    // H1
    T_local_asm_prod_.set_mode ( no_L2, H1_temp, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( T_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double temp_h1 = std::sqrt ( global_val );
    double temp_w12 = std::sqrt ( temp_l2 * temp_l2 + temp_h1 * temp_h1 );

    // **************************************************************
    // difference
    LAD::VectorType diff;
    diff.CloneFrom ( *sol );
    diff.Axpy ( *sol_prev, -1. );

#ifdef USE_HYPRE
    diff.Update ( );
#else
    diff.UpdateCouplings ( );
#endif

    T_local_asm_prod_.set_left_vector ( diff );
    T_local_asm_prod_.set_right_vector ( diff );

    // temperature
    // L2
    T_local_asm_prod_.set_mode ( L2_temp, no_H1, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( T_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double diff_temp_l2 = std::sqrt ( global_val );

    // H1
    T_local_asm_prod_.set_mode ( no_L2, H1_temp, no_H2 );
    cell_vals.clear ( );
    global_asm_.assemble_scalar ( *this->space_, boost::ref ( T_local_asm_prod_ ), cell_vals );

    local_val = 0.;
    local_val = std::accumulate ( cell_vals.begin ( ), cell_vals.end ( ), 0. );
    global_val = 0.;
    MPI_Reduce ( &local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, master_rank ( ), comm_ );
    double diff_temp_h1 = std::sqrt ( global_val );
    double diff_temp_w12 = std::sqrt ( diff_temp_l2 * diff_temp_l2 + diff_temp_h1 * diff_temp_h1 );

    // file output
    if ( rank ( ) == master_rank ( ) )
    {
        std::string path = this->root_ + filename + "." + static_cast < ostringstream* > ( &( ostringstream ( ) << this->adapt_counter_ ) )->str ( ) + ".txt";
        ofstream out;
        out.open ( path.c_str ( ), ios::out | ios::app );

        out.precision ( 10 );
        out << std::scientific;
        out << time << " " << temp_l2 << " " << temp_h1 << " " << diff_temp_l2 << " " << diff_temp_h1 << "\n";
        out.close ( );

        std::cout << "  L2 (temp_c):            " << temp_l2 << std::endl
                << "  H1 (temp_c):            " << temp_h1 << std::endl
                << "  L2 (temp_c - temp_p):   " << diff_temp_l2 << std::endl
                << "  H1 (temp_c - temp_p):   " << diff_temp_h1 << std::endl
                << "  rel W_12 diff:          " << diff_temp_w12 / temp_w12 << std::endl
                << " ---------------------    " << std::endl;
    }

    timer.stop ( );
    if ( rank ( ) == master_rank ( ) )
        std::cout << "  took " << timer.get_duration ( ) << " sec" << std::endl;
}

/// Compute characterisitc quantities

void MetFlowConvDiffApp::compute_char_quant ( )
{
}
