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

#include "test.h"
#include "hiflow.h"

using namespace hiflow;

int main ( int argc, char** argv )
{

    TensorIndex<2> tind = make_tensor_index ( 4, 6 );

    std::cout << "tind s = " << tind.size ( ) << ", " << tind.stride ( 0 ) << ", " << tind.stride ( 1 ) << "\n";
    std::cout << "tind d0 = " << tind.dim ( 0 ) << ", " << tind.dim ( 1 ) << "\n";
    Tensor<2> A ( tind );

    for ( int i = 0; i < tind.dim ( 0 ); ++i )
    {
        for ( int j = 0; j < tind.dim ( 1 ); ++j )
        {
            A ( i, j ) = i + j;
        }
    }

    std::vector<double> Av ( A.as_array ( ), A.as_array ( ) + tind.size ( ) );

    std::cout << "A = " << string_from_range ( Av.begin ( ), Av.end ( ) ) << "\n";

    return 0;
}
