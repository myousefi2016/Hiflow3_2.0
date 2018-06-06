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
/// @date 2015-06-30

#include "hiflow.h"
#include "linear_algebra/lmp/cuda/GPUcublas2_CSR_lMatrix.h"
#include "linear_algebra/lmp/cuda/GPUcublas2_lVector.h"
#include "linear_algebra/lmp/cuda/lvector_gpu.h"
#include "linear_algebra/lmp/cuda/lmatrix_csr_gpu.h"
#include <unittestpp.h>

using namespace hiflow::la;

template <class MatrixType, class VectorType>
struct TestFunctor
{

    void operator() ( ) const
    {
        here ( );
        int *indices = new int[3];
        indices[0] = 0;
        indices[1] = 1;
        indices[2] = 2;
        here ( );
        double *values = new double[3];
        values[0] = 1.2;
        values[1] = -2.3;
        values[2] = 0.7;
        here ( );
        VectorType x ( 3, "x" );
        here ( );
        x.SetValues ( indices, 3, values );
        here ( );
        CPUsimple_CSR_lMatrix<double> A ( 9, 3, 3, "A" );
        here ( );
        // Prepare dense CSR structure
        for ( int i = 0; i < 4; ++i ) A.matrix.row[i] = i * 3;
        for ( int i = 0; i < 3; ++i )
            for ( int j = 0; j < 3; ++j ) A.matrix.col[i * 3 + j] = j;
        here ( );
        A.add_value ( 0, 0, 2.0 );
        A.add_value ( 0, 1, -1.0 );
        A.add_value ( 0, 2, 0.0 );
        A.add_value ( 1, 0, -1.0 );
        A.add_value ( 1, 1, 2.0 );
        A.add_value ( 1, 2, -1.0 );
        A.add_value ( 2, 0, 0.0 );
        A.add_value ( 2, 1, -1.0 );
        A.add_value ( 2, 2, 2.0 );
        here ( );
        MatrixType gpuA ( 9, 3, 3, "gpuA" );
        gpuA.CopyStructureFrom ( A );
        gpuA.CopyFrom ( A );
        here ( );
        VectorType y ( 3, "y" );
        y.Zeros ( );
        gpuA.VectorMult ( x, &y );
        here ( );
        y.GetValues ( indices, 3, values );
        here ( );
        CHECK_CLOSE ( 4.7, values[0], 1e-12 );
        CHECK_CLOSE ( -6.5, values[1], 1e-12 );
        CHECK_CLOSE ( 3.7, values[2], 1e-12 );
    }
};

TEST ( GPUscalar_CSR_lMatrix_VectorMult )
{
    std::cout << "\nscalar...\n";
    TestFunctor<GPUscalar_CSR_lMatrix<double>, GPUblas_lVector<double> >( )( );
    std::cout << "... passed\n";
}

TEST ( GPUscalartex_CSR_lMatrix_VectorMult )
{
    std::cout << "\nscalartex\n";
    TestFunctor<GPUscalartex_CSR_lMatrix<double>, GPUblas_lVector<double> >( )( );
    std::cout << "... passed\n";
}

TEST ( GPUcublas2_CSR_lMatrix_VectorMult )
{
    std::cout << "\ncublas2\n";
    TestFunctor<GPUcublas2_CSR_lMatrix<double>, GPUcublas2_lVector<double> >( )( );
    std::cout << "... passed\n";
}

int main ( int argc, char** argv )
{
    return UnitTest::RunAllTests ( );
}
