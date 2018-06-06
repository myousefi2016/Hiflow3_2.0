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

#include "linear_solver/amg/intersection_size.h"
#include "gtest/gtest.h"
#include <map>
#include <set>
#include <vector>

using namespace hiflow::AMG;

TEST ( intersection_size, vector )
{
    std::vector<int> s1;
    int i1[] = { 2, 3 };
    s1.insert ( s1.begin ( ), i1, i1 + ( sizeof (i1 ) / sizeof (i1[0] ) ) );

    std::vector<int> s2;
    int i2[] = { 3, 4 };
    s2.insert ( s2.begin ( ), i2, i2 + ( sizeof (i2 ) / sizeof (i2[0] ) ) );

    // Check dimension
    EXPECT_EQ ( 1, intersection_size ( s1, s2 ) );
}

TEST ( intersection_size, set )
{
    std::set<int> s1;
    int i1[] = { 2, 3 };
    s1.insert ( i1, i1 + ( sizeof (i1 ) / sizeof (i1[0] ) ) );

    std::set<int> s2;
    int i2[] = { 3, 4 };
    s2.insert ( i2, i2 + ( sizeof (i2 ) / sizeof (i2[0] ) ) );

    // Check dimension
    EXPECT_EQ ( 1, intersection_size ( s1, s2 ) );
}

TEST ( intersection_size, set_with_map )
{
    std::set<int> c1;
    c1.insert ( 1 );
    c1.insert ( 3 );
    c1.insert ( 2 );

    std::map<int, int> c2;
    c2.insert ( std::make_pair ( 1, 1 ) );
    c2.insert ( std::make_pair ( 2, 3 ) );

    // Check dimension
    EXPECT_EQ ( 2, intersection_size ( c1, c2 ) );
}
