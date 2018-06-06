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

#include "hiflow.h"

#include <iostream>
#include <fstream>
#include "test.h"

using namespace hiflow;

int main ( int argc, char** argv )
{
    std::ofstream file ( "test.bin", std::ios::binary );

    //int v = 4;
    //long k = 15;
    //float f = 1.0f;

    std::vector<float> vec ( 4, 3. );

    //  write_binary(file, v);
    //  write_binary(file, k);
    //  write_binary(file, f);
    write_binary ( file, vec );

    return 0;
}
