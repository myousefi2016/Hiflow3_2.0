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
#include "linear_algebra/lmp/solvers/cg.h"
#include "linear_algebra/lmp/lmatrix_dense_cpu.h"
#include "linear_algebra/lmp/cuda/GPUcublas2_CSR_lMatrix.h"
#include "linear_algebra/lmp/cuda/GPUcublas2_lVector.h"
#include "linear_algebra/lmp/init_vec_mat.h"
#include "gtest/gtest.h"

using namespace hiflow::la;

class NoPreconditioner
{
};

/// Workaround for missing C++11 feature std::tuple

template <class MT, class VT, class DT, class PT>
struct TypeList
{
    typedef MT MatrixType;
    typedef VT VectorType;
    typedef DT DataType;
    typedef PT PreconType;
};

/// Test fixture

template <typename T>
class LocalSolverTest : public ::testing::Test
{
};

TYPED_TEST_CASE_P ( LocalSolverTest );

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

/// Test function, doing the real work

template <class TypeParam>
void test ( std::string const& filename )
{
    typedef typename TypeParam::MatrixType MatrixType;
    typedef typename TypeParam::VectorType VectorType;
    typedef typename TypeParam::PreconType PreconType;
    typedef typename TypeParam::DataType DataType;

    VectorType b;
    read_vector ( b, filename + ".rhs.txt" );

    VectorType x;
    read_vector ( x, filename + ".x0.txt" );

    MatrixType A;
    read_matrix ( A, filename + ".A.txt" );

    double rel_eps = 0.0;
    double abs_eps = 0.0;
    int max_iter = 100;
    int print_level = -1;

    EXPECT_EQ ( 0, cg ( &x, &b, &A, rel_eps, abs_eps, max_iter, print_level ) );

    // Check result
    VectorType x_result;
    read_vector ( x_result, filename + ".sol.txt" );

    VectorType r ( x.get_size ( ), "r" );
    r.CopyFrom ( x_result );
    r.ScaleAdd ( -1.0, x );

    double dot = r.Dot ( r );
    EXPECT_LT ( dot, 1.e6 );
}

/// Test body
/// Workaround, because GTest do not support type- and value-parameterized tests together

TYPED_TEST_P ( LocalSolverTest, Equation1 )
{
    test<TypeParam>( "../test/LinearAlgebraTest/data/Equation1" );
}

REGISTER_TYPED_TEST_CASE_P ( LocalSolverTest, Equation1 );

typedef ::testing::Types<
TypeList<CPUsimple_DENSE_lMatrix<double>, CPUsimple_lVector<double>, double, NoPreconditioner>
> MyTypes;

INSTANTIATE_TYPED_TEST_CASE_P ( My, LocalSolverTest, MyTypes );
