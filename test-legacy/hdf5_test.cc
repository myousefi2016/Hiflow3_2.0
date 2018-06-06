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

/// \author Teresa Beck

#include "hiflow.h"
#include<map>
#include<fstream>
#include<iostream>

#include "hdf5.h"

#include "test.h"
#include "common/hdf5_tools.h"

#ifndef WITH_HDF5
#    error "This test only works when HDF5 support is enabled";
#endif

using namespace hiflow;

void write_dataset ( )
{
    std::string filename ( "test_hdf5.h5" );
    H5FilePtr file_ptr_w ( new H5File ( filename, "w", MPI_COMM_WORLD ) );
    H5GroupPtr group_ptr_w ( new H5Group ( file_ptr_w, "indices", "w" ) );

    // Write map<string, string>
    std::map<std::string, std::string> w_map;
    w_map["Ceterum"] = "censeo";
    w_map["Carthaginem"] = "esse";
    w_map["delendam"] = "!";
    write_map ( group_ptr_w, "map<string, string>", w_map );

    // Write single value integer
    int w_sv_i = 1234;
    write_value ( group_ptr_w, "single value integer", w_sv_i );

    // Write array of integers
    int size_ar_i = 10;
    int* w_ar_i;
    w_ar_i = ( int * ) malloc ( sizeof (int )*size_ar_i );
    for ( int i = 0; i < size_ar_i; i++ ) w_ar_i[i] = i;
    write_array ( group_ptr_w, "array of integers", w_ar_i, size_ar_i );
    free ( w_ar_i );

    // Write vector of integers
    int size_v_i = 10;
    std::vector<int> w_v_i ( size_v_i );
    for ( int i = 0; i < size_v_i; i++ ) w_v_i[i] = i;
    write_array ( group_ptr_w, "vector of integers", w_v_i );

    // Write single value string
    std::string w_sv_s = "Discite moniti!";
    write_value ( group_ptr_w, "single value string", w_sv_s );

    // Write vector of strings
    std::vector<std::string> w_v_s ( 3 );
    w_v_s[0] = "Alea";
    w_v_s[1] = "iacta";
    w_v_s[2] = "est!";
    write_array ( group_ptr_w, "vector of strings", w_v_s );

    // Write array of strings
    int size_ar_s = 2;
    std::string *w_ar_s = new std::string[size_ar_s];
    w_ar_s[0] = std::string ( "Cui" );
    w_ar_s[1] = std::string ( "bono?" );
    write_array ( group_ptr_w, "array of strings", w_ar_s, size_ar_s );
    delete[] w_ar_s;
}

void read_dataset ( )
{
    std::string filename ( "test_hdf5.h5" );
    int output = 30;
    H5FilePtr file_ptr_r ( new H5File ( filename, "r", MPI_COMM_WORLD ) );
    H5GroupPtr group_ptr_r ( new H5Group ( file_ptr_r, "indices", "r" ) );

    // Read map<string, string>
    std::map<std::string, std::string> r_map;
    std::map<std::string, std::string>::iterator it;
    read_map ( group_ptr_r, "map<string, string>", r_map );
    std::cout.width ( output );
    std::cout << std::left << "  Map is ";
    for ( it = r_map.begin ( ); it != r_map.end ( ); it++ )
        std::cout << "(" << ( *it ).first << ", " << ( *it ).second << ")  ";
    std::cout << "\n";
    TEST_EQUAL ( r_map["Ceterum"], "censeo" );
    TEST_EQUAL ( r_map["Carthaginem"], "esse" );
    TEST_EQUAL ( r_map["delendam"], "!" );

    // Read single value string
    std::string r_sv_s;
    read_value ( group_ptr_r, "single value string", r_sv_s );
    std::cout.width ( output );
    std::cout << std::left << "  Single string is " << r_sv_s << "\n";
    TEST_EQUAL ( r_sv_s, "Discite moniti!" );

    // Read vector of strings
    std::vector<std::string> r_v_s;
    read_array ( group_ptr_r, "vector of strings", r_v_s );
    std::cout.width ( output );
    std::cout << std::left << "  Vector of strings is ";
    for ( int i = 0; i < static_cast < int > ( r_v_s.size ( ) ); i++ ) std::cout << r_v_s[i] << "  ";
    std::cout << "\n";
    TEST_EQUAL ( r_v_s[0], "Alea" );
    TEST_EQUAL ( r_v_s[1], "iacta" );
    TEST_EQUAL ( r_v_s[2], "est!" );

    // Read array of strings
    int size_r_ar_s = get_dataset_size ( group_ptr_r, "array of strings" );
    std::string *r_ar_s = new std::string[size_r_ar_s];
    read_array ( group_ptr_r, "array of strings", r_ar_s );
    std::cout.width ( output );
    std::cout << std::left << "  Array of strings is ";
    for ( int i = 0; i < size_r_ar_s; i++ ) std::cout << r_ar_s[i] << "  ";
    std::cout << "\n";
    TEST_EQUAL ( r_ar_s[0], "Cui" );
    TEST_EQUAL ( r_ar_s[1], "bono?" );

    // Read single value integer
    int r_sv_i;
    read_value ( group_ptr_r, "single value integer", &r_sv_i );
    std::cout.width ( output );
    std::cout << std::left << "  Single integer is ";
    std::cout << r_sv_i << "\n";
    TEST_EQUAL ( r_sv_i, 1234 );

    // Read array of integers
    int size_r_ar_i = get_dataset_size ( group_ptr_r, "array of integers" );
    int* r_ar_i;
    r_ar_i = ( int * ) malloc ( sizeof (int )*size_r_ar_i );
    read_array ( group_ptr_r, "array of integers", r_ar_i );
    std::cout.width ( output );
    std::cout << std::left << "  Array of int is ";
    for ( int i = 0; i < size_r_ar_i; i++ ) std::cout << r_ar_i[i] << "  ";
    std::cout << "\n";
    for ( int i = 0; i < size_r_ar_i; i++ ) TEST_EQUAL ( r_ar_i[i], i );
    free ( r_ar_i );

    // Read Vector of integers
    std::vector<int> r_v_i;
    read_array ( group_ptr_r, "vector of integers", r_v_i );
    std::cout.width ( output );
    std::cout << std::left << "  Vector of int is ";
    for ( int i = 0; i < static_cast < int > ( r_v_i.size ( ) ); i++ ) std::cout << r_v_i[i] << "  ";
    std::cout << "\n";
    for ( int i = 0; i < static_cast < int > ( r_v_i.size ( ) ); i++ ) TEST_EQUAL ( r_v_i[i], i );
}

void get_time ( void(*application )( ) )
{
    double start, end;
    start = MPI_Wtime ( );
    application ( );
    end = MPI_Wtime ( );
    printf ( "  Time: %4.10f [sec] \n \n", end - start );
}

int main ( int argc, char *argv[] )
{

    MPI_Init ( &argc, &argv );

    printf ( "Writing datasets \n" );
    get_time ( write_dataset );

    printf ( "Reading datasets \n" );
    get_time ( read_dataset );

    MPI_Finalize ( );
    return 0;
}
