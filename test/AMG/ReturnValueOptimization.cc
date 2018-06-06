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

#include "gtest/gtest.h"

struct S1
{
    // Default constructor

    S1 ( ) : i ( 0 )
    {
        ++num_default_constructor_calls;
    };

    // Parameter constructor

    S1 ( int i ) : i ( i )
    {
        ++num_parameter_constructor_calls;
    };

    // Copy constructor

    S1 ( S1 const& other ) : i ( other.i )
    {
        ++num_copy_constructor_calls;
    }

    // Copy assignment operator

    S1& operator= ( S1 const& other )
    {
        ++num_copy_assignment_calls;
        if ( &other != this ) i = other.i;
        return *this;
    }

    int i;

    static int num_parameter_constructor_calls;
    static int num_default_constructor_calls;
    static int num_copy_constructor_calls;
    static int num_copy_assignment_calls;
};

int S1::num_parameter_constructor_calls = 0;
int S1::num_default_constructor_calls = 0;
int S1::num_copy_constructor_calls = 0;
int S1::num_copy_assignment_calls = 0;

struct S2
{

    S1 operator() ( ) const
    {
        return S1 ( );
    }
};

struct S3
{

    S3 ( ) : o1 ( 2 ), o2 ( 3 )
    {
    };

    S1 o1;
    S1 o2;
};

struct S4
{

    S3 operator() ( ) const
    {
        S3 s3;
        return s3;
    }
};

/// Test if copy constructor/assignment will be optimized out for a simple struct by compiler RVO.

TEST ( ReturnValueOptimization, ReturnStruct )
{
    S1 s1 = S2 ( )( );
    EXPECT_EQ ( 1, S1::num_default_constructor_calls );
    EXPECT_EQ ( 0, S1::num_parameter_constructor_calls );
    EXPECT_EQ ( 0, S1::num_copy_constructor_calls );
    EXPECT_EQ ( 0, S1::num_copy_assignment_calls );
}

/// Test if copy constructor/assignment will be optimized out for a struct of structs by compiler RVO.

TEST ( ReturnValueOptimization, ReturnStructOfStruct )
{
    S3 s3 = S4 ( )( );
    EXPECT_EQ ( 1, S1::num_default_constructor_calls );
    EXPECT_EQ ( 2, S1::num_parameter_constructor_calls );
    EXPECT_EQ ( 0, S1::num_copy_constructor_calls );
    EXPECT_EQ ( 0, S1::num_copy_assignment_calls );

    EXPECT_EQ ( 2, s3.o1.i );
    EXPECT_EQ ( 3, s3.o2.i );
}
