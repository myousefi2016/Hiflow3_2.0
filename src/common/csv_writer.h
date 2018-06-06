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

#ifndef HIFLOW_COMMON_CSV_WRITER_H
#    define HIFLOW_COMMON_CSV_WRITER_H

#    include <string>
#    include <vector>
#    include <iostream>
#    include <fstream>
#    include <sstream>

/// @author Simon Gawlok
/// @brief Class for handling CSV files. Provides functionality to create, write and 
/// read CSV files. Templatization is used such that I/O with data of different types 
/// can be done.

template<class T>
class CSVWriter
{
  public:
    /// Constructor
    /// \param[in] filename Filename of CSV file

    CSVWriter ( std::string filename )
    {
        filename_ = filename;
    }

    /// Standard constructor

    CSVWriter ( )
    {
    }

    /// Initialize/set filename
    /// \param[in] filename Filename of CSV file

    void InitFilename ( std::string filename )
    {
        filename_ = filename;
    }

    /// Initialize CSV file with headings of columns.
    /// WARNING: If file already exists, it will be overwritten!!!
    /// \param[in] col_names Vector of strings containing the headings of the columns
    /// in the created CSV file

    void Init ( std::vector<std::string> col_names )
    {
        // Open CSV file
        std::ofstream out_file ( filename_.c_str ( ), std::ios::out );

        // If file is really open...
        if ( out_file.is_open ( ) )
        {
            for ( int i = 0; i < col_names.size ( ); ++i )
            {
                out_file << col_names[i];

                if ( i == col_names.size ( ) - 1 )
                {
                    out_file << "\n";
                }
                else
                {
                    out_file << ", ";
                }
            }
            // Write output to file
            out_file.flush ( );
        }
        else
        {
            std::cerr << "Could not create CSV file " << filename_ << "!" << std::endl;
            exit ( -1 );
        }
    }

    /// Write row of data to CSV file. Data are appended as last line in
    /// CSV file.
    /// \param[in] row_of_values Vector containing the entries for the next line/dataset
    // in the CSV file.

    void write ( std::vector<T> &row_of_values )
    {
        // Open CSV file
        std::ofstream out_file ( filename_.c_str ( ), std::ios::app );

        // If file is really open...
        if ( out_file.is_open ( ) )
        {
            for ( int i = 0; i < row_of_values.size ( ); ++i )
            {
                out_file << row_of_values[i];

                if ( i == row_of_values.size ( ) - 1 )
                {
                    out_file << "\n";
                }
                else
                {
                    out_file << ", ";
                }
            }
            // Write output to file
            out_file.flush ( );
        }
    }

    /// Read in the values (without headings!) of a CSV file. Data are stored
    /// as a vector (each element is a dataset) of vectors (the datasets themselves).
    /// \param[out] values Read-in content of CSV file

    void read ( std::vector<std::vector<T> > &values )
    {
        // Clear possible old values
        values.clear ( );
        // Open CSV file
        std::ifstream in_file ( filename_.c_str ( ), std::ios::in );

        // If file is really open...
        if ( in_file.is_open ( ) )
        {
            // read first line ("headers")
            std::string line;

            std::getline ( in_file, line );

            // index of current row in values
            int row = 0;

            // start reading real values
            while ( std::getline ( in_file, line ) )
            {
                // create new row in values
                values.push_back ( std::vector<T>( ) );

                line += ",";
                size_t pos = 0;
                while ( ( pos = line.find ( "," ) ) != std::string::npos )
                {
                    std::string val = line.substr ( 0, pos );
                    line = line.substr ( pos + 1 );

                    std::stringstream str_tmp;
                    str_tmp << val;

                    T tmp = T ( );
                    str_tmp >> tmp;
                    values[row].push_back ( tmp );
                }
                ++row;
            }
        }
    }

  private:
    /// Filename of CSV file
    std::string filename_;
};

#endif
