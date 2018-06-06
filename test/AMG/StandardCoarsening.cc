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
/// @date 2015-10-21

#include "hiflow.h"
#include "linear_solver/amg/LevelGenerator.h"
#include "linear_solver/amg/StandardCoarsening.h"
#include "linear_solver/amg/StandardInterpolation.h"
#include "linear_solver/amg/TripleProduct.h"
#include "gtest/gtest.h"

using namespace hiflow::la;
using namespace hiflow::AMG;

TEST ( StandardCoarseningTest, A1 )
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

    StandardCoarseningSettings coarsening_settings;
    coarsening_settings.strength_threshold = 0.25;

    LevelGeneratorType::ResultType result;
    StandardCoarsening<LevelGeneratorType>::ConnectionType connection = StandardCoarsening<LevelGeneratorType>( coarsening_settings ) ( result, A );

    EXPECT_EQ ( 3, connection.coarse_variables.size ( ) );
    EXPECT_TRUE ( connection.coarse_variables.count ( 0 ) );
    EXPECT_TRUE ( connection.coarse_variables.count ( 1 ) );
    EXPECT_TRUE ( connection.coarse_variables.count ( 2 ) );
    EXPECT_EQ ( 0, connection.fine_variables.size ( ) );
}

TEST ( StandardCoarseningTest, A2 )
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

    StandardCoarseningSettings coarsening_settings;
    coarsening_settings.strength_threshold = 0.25;

    LevelGeneratorType::ResultType result;
    StandardCoarsening<LevelGeneratorType>::ConnectionType connection = StandardCoarsening<LevelGeneratorType>( coarsening_settings ) ( result, A );

    EXPECT_EQ ( 135, connection.fine_variables.size ( ) );
    EXPECT_EQ ( 56, connection.coarse_variables.size ( ) );

    IndexVector isec;
    std::set_intersection ( connection.coarse_variables.begin ( ), connection.coarse_variables.end ( ),
                            connection.fine_variables.begin ( ), connection.fine_variables.end ( ),
                            std::back_inserter ( isec ) );

    EXPECT_TRUE ( isec.empty ( ) );

    // Only for visualisation of splitting using pyamg-examples
    //for (size_t i(0); i != A.get_num_row(); ++i) {
    //	if (std::find(connection.coarse_variables.begin(), connection.coarse_variables.end(), i) != connection.coarse_variables.end()) std::cout << 1 << ", ";
    //	else std::cout << 0 << ", ";
    //}
}
