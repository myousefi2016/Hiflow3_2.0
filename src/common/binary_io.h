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

#ifndef HIFLOW_COMMON_BINARY_IO_H
#    define HIFLOW_COMMON_BINARY_IO_H

#    include <iostream>
#    include <cstring>
#    include <climits>
#    include <string>
#    include <vector>
#    include <unistd.h>

/// \brief Functions to read and write binary strings.
/// \author Staffan Ronnas

typedef unsigned array_size_t;

namespace hiflow
{

    template <class T>
    inline
    void read_binary ( std::istream& is, T& val )
    {
        char bytes[sizeof (T )];

        // read separate bytes
        is.read ( static_cast < char* > ( &bytes[0] ), sizeof (T ) );

        // construct value bytewise
        std::memset ( &val, 0x0, sizeof (T ) ); // initialize with 0
        for ( int i = 0; i < sizeof (T ); ++i )
        {
            // set n:th byte (little-endian)
            *( reinterpret_cast < char* > ( &val ) + i ) |= bytes[i];
        }
    }

    template <class T>
    inline
    void read_binary ( std::istream& is, std::vector<T>& val )
    {
        array_size_t sz;
        read_binary ( is, sz );

        std::vector<char> buf ( sz * sizeof (T ) );
        if ( sz > 0 )
        {
            is.read ( reinterpret_cast < char* > ( &buf[0] ), sz * sizeof (T ) );
        }
    }

    inline
    void read_binary ( std::istream& is, std::string& str )
    {
        array_size_t sz;
        read_binary ( is, sz );
        std::vector<char> buf ( sz * sizeof (char ) );
        if ( sz > 0 )
        {
            is.read ( &buf[0], sz * sizeof (char ) );
            str = std::string ( &buf[0], sz );
        }
    }

    template <class T>
    inline
    void write_binary ( std::ostream& os, T val )
    {
        char bytes[sizeof (T )];

        // extract bytes from val
        for ( int i = 0; i < sizeof (T ); ++i )
        {
            bytes[i] = ( *( reinterpret_cast < char* > ( &val ) + i ) ) & 0xFF;
        }

        os.write ( static_cast < char* > ( &bytes[0] ), sizeof (T ) );
    }

    template <class T>
    inline
    void write_binary ( std::ostream& os, const std::vector<T>& vec )
    {
        write_binary ( os, static_cast < array_size_t > ( vec.size ( ) ) );

        if ( vec.size ( ) > 0 )
        {
            os.write ( reinterpret_cast < char* > ( const_cast < T* > ( &vec[0] ) ), vec.size ( ) * sizeof (T ) );
        }
    }

    inline
    void write_binary ( std::ostream& os, const std::string& str )
    {
        write_binary ( os, static_cast < array_size_t > ( str.size ( ) ) );

        if ( str.size ( ) > 0 )
        {
            os.write ( str.data ( ), str.size ( ) * sizeof (char ) );
        }
    }
}
#endif
