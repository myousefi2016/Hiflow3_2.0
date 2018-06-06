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
/// @date 2015-10-30

#include "hiflow.h"
#include "linear_solver/amg/LevelGenerator.h"
#include "linear_solver/amg/StandardCoarsening.h"
#include "linear_solver/amg/StandardInterpolation.h"
#include "linear_solver/amg/TripleProduct.h"
#include "gtest/gtest.h"

using namespace hiflow::la;
using namespace hiflow::AMG;

TEST ( StandardInterpolationTest, A1 )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            StandardInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A1.txt" ).c_str ( ) );

    LevelGeneratorType::ResultType result;
    StandardCoarsening<LevelGeneratorType>::ConnectionType connection = StandardCoarsening<LevelGeneratorType>( )( result, A );
    StandardInterpolation<LevelGeneratorType>( )( result, connection, A );

    EXPECT_EQ ( 3, result.ptr_interpolation_matrix->get_num_row ( ) );
    EXPECT_EQ ( 3, result.ptr_interpolation_matrix->get_num_col ( ) );
    EXPECT_EQ ( 3, result.ptr_interpolation_matrix->get_nnz ( ) );
    EXPECT_DOUBLE_EQ ( 1.0, result.ptr_interpolation_matrix->matrix.val[0] );
    EXPECT_DOUBLE_EQ ( 1.0, result.ptr_interpolation_matrix->matrix.val[1] );
    EXPECT_DOUBLE_EQ ( 1.0, result.ptr_interpolation_matrix->matrix.val[2] );
}

TEST ( StandardInterpolationTest, A2 )
{
    typedef CPUsimple_CSR_lMatrix<double> MatrixType;
    typedef LevelGenerator<
            MatrixType,
            StandardCoarsening,
            StandardInterpolation,
            TripleProduct
            > LevelGeneratorType;

    MatrixType A;
    A.ReadFile ( ( std::string ( DATA_DIR ) + "A2.txt" ).c_str ( ) );
    EXPECT_EQ ( 191, A.get_num_row ( ) );

    LevelGeneratorType::ResultType result;
    StandardCoarsening<LevelGeneratorType>::ConnectionType connection = StandardCoarsening<LevelGeneratorType>( )( result, A );

    EXPECT_EQ ( 135, connection.fine_variables.size ( ) );
    EXPECT_EQ ( 56, connection.coarse_variables.size ( ) );

    StandardInterpolation<LevelGeneratorType>( )( result, connection, A );

    EXPECT_EQ ( 191, result.ptr_interpolation_matrix->get_num_row ( ) );
    EXPECT_EQ ( 56, result.ptr_interpolation_matrix->get_num_col ( ) );
    EXPECT_EQ ( 345, result.ptr_interpolation_matrix->get_nnz ( ) );
}
