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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-07-13

#include "hiflow.h"
#include "linear_algebra/lmp/solvers/cg.h"
#include "linear_algebra/lmp/solvers/gmres.h"
#include "linear_algebra/lmp/lmatrix_dense_cpu.h"
#include "linear_algebra/lmp/init_vec_mat.h"
#include "linear_algebra/lmp/cuda/GPUcublas2_CSR_lMatrix.h"
#include "linear_algebra/lmp/cuda/GPUcublas2_lVector.h"
#include "linear_algebra/lmp/cuda/lmatrix_csr_gpu.h"
#include "tools/ScopedTimer.h"
#include <iostream>
#include <string>

using namespace hiflow::la;

template <class T>
void set_num_threads ( T& t, int nbThreads )
{
}

template <>
void set_num_threads ( CPUopenmp_CSR_lMatrix<double>& t, int nbThreads )
{
    t.set_num_threads ( nbThreads );
}

template <>
void set_num_threads ( CPUmkl_CSR_lMatrix<double>& t, int nbThreads )
{
    t.set_num_threads ( nbThreads );
}

template <>
void set_num_threads ( CPUopenmp_lVector<double>& t, int nbThreads )
{
    t.set_num_threads ( nbThreads );
}

template <>
void set_num_threads ( CPUmkl_lVector<double>& t, int nbThreads )
{
    t.set_num_threads ( nbThreads );
}

template <class VectorType>
void read_vector ( VectorType& v, std::string const& filename )
{
    v.ReadFile ( filename.c_str ( ) );
}

template <>
void read_vector ( GPUblas_lVector<double>& v, std::string const& filename )
{
    CPUsimple_lVector<double> simplev;
    read_vector ( simplev, filename );
    v.Init ( simplev.get_size ( ), "v" );
    v.CopyFrom ( simplev );
}

template <>
void read_vector ( GPUcublas2_lVector<double>& v, std::string const& filename )
{
    CPUsimple_lVector<double> simplev;
    read_vector ( simplev, filename );
    v.Init ( simplev.get_size ( ), "v" );
    v.CopyFrom ( simplev );
}

template <class MatrixType>
void read_matrix ( MatrixType& A, std::string const& filename )
{
    A.ReadFile ( filename.c_str ( ) );
}

template <>
void read_matrix ( GPUscalartex_CSR_lMatrix<double>& A,
                   std::string const& filename )
{
    CPUsimple_CSR_lMatrix<double> simpleA;
    read_matrix ( simpleA, filename );
    A.Init ( simpleA.get_nnz ( ), simpleA.get_num_row ( ), simpleA.get_num_col ( ), "A" );
    A.CopyStructureFrom ( simpleA );
    A.CopyFrom ( simpleA );
}

template <>
void read_matrix ( GPUcublas2_CSR_lMatrix<double>& A,
                   std::string const& filename )
{
    CPUsimple_CSR_lMatrix<double> simpleA;
    read_matrix ( simpleA, filename );
    A.Init ( simpleA.get_nnz ( ), simpleA.get_num_row ( ), simpleA.get_num_col ( ), "A" );
    A.CopyStructureFrom ( simpleA );
    A.CopyFrom ( simpleA );
}

template <class VectorType, class MatrixType>
struct Benchmark
{

    template <class SolverType>
    void operator() ( std::string const& file_A, std::string const& file_b,
            std::string const& file_x_result,
            std::string const& file_x_start, int nbThreads,
            int nbRepetitions, SolverType const& solver )
    {
        VectorType b;
        {
            ScopedTimer timer ( "Read vector b" );
            read_vector ( b, file_b );
            set_num_threads ( b, nbThreads );
        }

        VectorType x;
        {
            ScopedTimer timer ( "Read vector x" );
            read_vector ( x, file_x_start );
            set_num_threads ( x, nbThreads );
        }

        MatrixType A;
        {
            ScopedTimer timer ( "Read matrix A" );
            read_matrix ( A, file_A );
            set_num_threads ( A, nbThreads );
        }

        double rel_eps = 1e-8;
        double abs_eps = 1e-14;
        int max_iter = 1000;
        int print_level = -1;

        for ( int i = 0; i < nbRepetitions - 1; ++i )
        {
            // Overwrite the result with the initial x vector
            read_vector ( x, file_x_start );
            set_num_threads ( x, nbThreads );

            if ( solver ( &x, &b, &A, rel_eps, abs_eps, max_iter, print_level ) != 0 )
                throw std::runtime_error ( "Error in solver." );
        }

        // Do it a last one for printing residue
        print_level = 0;

        // Overwrite the result with the initial x vector
        read_vector ( x, file_x_start );
        set_num_threads ( x, nbThreads );

        x.print ( );
        b.print ( );
        A.print ( );

        if ( solver ( &x, &b, &A, rel_eps, abs_eps, max_iter, print_level ) != 0 )
            throw std::runtime_error ( "Error in solver." );

        // Check result
        VectorType x_result;
        read_vector ( x_result, file_x_result );

        VectorType r ( x.get_size ( ), "r" );
        r.CopyFrom ( x_result );
        r.ScaleAdd ( -1.0, x );

        double dot = r.Dot ( r );
        std::cout << "dot(x - x_result) = " << dot << std::endl;
        if ( std::isnan ( dot ) )
            std::cout << "WARNING: Resulting vector is not a number." << std::endl;
        if ( std::abs ( dot ) > 1e-6 )
            std::cout << "WARNING: Deviation between resulting and expected vector x "
                "larger than 1e-6."
                << std::endl;
    }
};

int main ( int argc, char** argv )
{
    if ( argc != 9 )
    {
        std::string programName = argv[0];
        std::cout << "Usage: " << programName
                << " <A> <b> <x_result> <x_start> <solver> <interface> <number "
                "of threads> <number of repetitions>"
                << std::endl;
        exit ( 1 );
    }

    std::string file_A = argv[1];
    std::string file_b = argv[2];
    std::string file_x_result = argv[3];
    std::string file_x_start = argv[4];
    std::string solver = argv[5];
    std::string interface = argv[6];
    int nbThreads = atoi ( argv[7] );
    int nbRepetitions = atoi ( argv[8] );
    int basisSize = 5;

    if ( solver == "cg" )
    {
        CGFunctor<double, true> solver;
        if ( interface == "simple" )
        {
            Benchmark<CPUsimple_lVector<double>, CPUsimple_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
        }
        else if ( interface == "openmp" )
        {
            Benchmark<CPUopenmp_lVector<double>, CPUopenmp_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
#ifdef WITH_MKL
        }
        else if ( interface == "mkl" )
        {
            Benchmark<CPUmkl_lVector<double>, CPUmkl_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
#endif
#ifdef WITH_CUDA
        }
        else if ( interface == "cuda" )
        {
            Benchmark<GPUblas_lVector<double>, GPUscalartex_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
#endif
#ifdef WITH_CUDA
        }
        else if ( interface == "cublas2" )
        {
            Benchmark<GPUcublas2_lVector<double>, GPUcublas2_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
#endif
        }
        else
        {
            std::cout << "Algo " << interface << " not supported." << std::endl;
            exit ( 1 );
        }
    }
    else if ( solver == "gmres" )
    {
        GMRESFunctor<double, true> solver ( basisSize );
        if ( interface == "simple" )
        {
            Benchmark<CPUsimple_lVector<double>, CPUsimple_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
        }
        else if ( interface == "openmp" )
        {
            Benchmark<CPUopenmp_lVector<double>, CPUopenmp_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
#ifdef WITH_MKL
        }
        else if ( interface == "mkl" )
        {
            Benchmark<CPUmkl_lVector<double>, CPUmkl_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
#endif
#ifdef WITH_CUDA
        }
        else if ( interface == "cuda" )
        {
            Benchmark<GPUblas_lVector<double>, GPUscalartex_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
#endif
#ifdef WITH_CUDA
        }
        else if ( interface == "cublas2" )
        {
            Benchmark<GPUcublas2_lVector<double>, GPUcublas2_CSR_lMatrix<double> >( )(
                    file_A, file_b, file_x_result, file_x_start, nbThreads, nbRepetitions,
                    solver );
#endif
        }
        else
        {
            std::cout << "Algo " << interface << " not supported." << std::endl;
            exit ( 1 );
        }
    }
    else
    {
        std::cout << "Solver " << solver << " not supported." << std::endl;
        exit ( 1 );
    }

    std::cout << "Done successfully. Have a nice day." << std::endl;
}
