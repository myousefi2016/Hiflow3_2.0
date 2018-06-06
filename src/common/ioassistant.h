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

#ifndef HIFLOW_IOASSISTANT_H_
#    define HIFLOW_IOASSISTANT_H_

/// \file ioassistant.h
///
/// \author Michael Schick 

#    include <iostream>
#    include <algorithm>
#    include <fstream>
#    include <vector>
#    include <string>
#    include <sstream>
#    include <cassert>

namespace hiflow
{

    class IOAssistant
    {
      public:

        IOAssistant ( )
        {
        }

        template<class DataType>
        inline void read ( const std::string& filename, std::vector<DataType>& result );
        template<class DataType>
        inline void read_binary ( const std::string& filename, std::vector<DataType>& result );
        template<class DataType>
        inline void write ( const std::string& filename, const std::vector<DataType>& input ) const;
        template<class DataType>
        inline void write_binary ( const std::string& filename, const std::vector<DataType>& input ) const;

    };

    /// INLINE FUNCTIONS

    template<class DataType>
    void IOAssistant::read ( const std::string& filename, std::vector<DataType>& result )
    {
        std::string strline; // Current line

        // Open file

        std::ifstream file ( ( filename + ".txt" ).c_str ( ) );
        if ( !file )
        {
            std::cerr << "\n Can't open file : " << filename + ".txt" << std::endl;
            exit ( 1 );
        }

        // First line

        std::getline ( file, strline );
        std::istringstream istl1 ( strline.c_str ( ) );

        int n;
        istl1 >> n;

        assert ( n > 0 );

        result.resize ( n, 0.0 );

        // Full vector

        for ( int k = 0; k < n; ++k )
        {
            std::getline ( file, strline );
            std::istringstream istlc ( strline.c_str ( ) );

            istlc >> result[k];

        } // for(int k=0; ...

    }

    template<class DataType>
    void IOAssistant::read_binary ( const std::string& filename, std::vector<DataType>& result )
    {
        // Open file

        std::ifstream file;

        file.open ( ( filename + ".bin" ).c_str ( ), std::ios::in | std::ios::binary );

        if ( !file )
        {
            std::cerr << "\n Can't open file : " << filename + ".bin" << std::endl;
            exit ( 1 );
        }

        int n;

        file.read ( ( char* ) &n, sizeof (int ) );

        assert ( n > 0 );

        // Reinit

        result.resize ( n, 0.0 );

        // Full vector

        file.read ( ( char* ) &( result.front ( ) ), sizeof (DataType ) * n );

        file.close ( );
    }

    /// write

    template<class DataType>
    void IOAssistant::write ( const std::string& filename, const std::vector<DataType>& input ) const
    {
        // Open the file

        std::ofstream file ( ( filename + ".txt" ).c_str ( ), std::ios::out );

        if ( !file )
        {
            std::cerr << "\n Can't open file : " << filename + ".txt" << std::endl;
            exit ( 1 );
        }

        // Set output format

        file.setf ( std::ios::scientific, std::ios::floatfield );
        file.precision ( 18 );

        // First line

        file << input.size ( ) << std::endl;

        // Loop on the coordinates

        for ( int k = 0; k < input.size ( ); ++k ) file << input[k] << std::endl;
    }

    template<class DataType>
    void IOAssistant::write_binary ( const std::string& filename, const std::vector<DataType>& input ) const
    {
        // Open the file

        std::ofstream file;
        file.open ( ( filename + ".bin" ).c_str ( ), std::ios::out | std::ios::binary );

        if ( !file )
        {
            std::cerr << "\n Can't open file : " << filename + ".bin" << std::endl;
            exit ( 1 );
        }

        // First line

        int mysize = input.size ( );
        file.write ( ( char* ) &mysize, sizeof (int ) );

        // all the coordinates

        file.write ( ( char* ) &( input.front ( ) ), sizeof (DataType ) * input.size ( ) );

        file.close ( );
    }

}
#endif
