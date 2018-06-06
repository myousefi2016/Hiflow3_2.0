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

#include "write_csv.h"

void write_csv ( std::vector<std::string> names, std::vector<std::vector<double> > data, std::string filename )
{
    fstream f;
    const char *file = filename.c_str ( );
    f.open ( file, ios::out );
    assert ( names.size ( ) == data.size ( ) );

    // write column names
    for ( int i = 0; i < names.size ( ); ++i )
    {
        f << names[i];
        if ( i != names.size ( ) - 1 )
        {
            f << ", ";
        }
    }
    f << endl;

    // write data
    for ( int i = 0; i < data[0].size ( ); ++i )
    {
        for ( int j = 0; j < names.size ( ); ++j )
        {
            f << data[j][i];
            if ( j != names.size ( ) - 1 )
            {
                f << ", ";
            }
        }
        f << endl;
    }
    f.close ( );
}
