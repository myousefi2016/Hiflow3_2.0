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

#ifndef HIFLOW_COMMON_TABLE_H
#    define HIFLOW_COMMON_TABLE_H

#    include <iomanip>
#    include <iostream>
#    include <map>
#    include <string>
#    include <sstream>
#    include <vector>

/// \author Staffan Ronnas

namespace hiflow
{

    class TableColumn
    {
      public:

        void insert ( const std::string& str );

        int num_rows ( ) const;

        const std::string& value ( int row ) const;

      private:
        static const std::string MISSING_STRING;

        std::vector< std::string > values_;
    };

    class Table
    {
      public:

        Table ( );

        template< class ValueType >
        inline void insert ( const std::string& column_name, const ValueType& value );

        // Print formatted table with given columns.
        void print ( std::ostream& os, const std::vector< std::string >& columns ) const;

        // Print CSV table with given columns.
        void print_csv ( std::ostream& os, const std::vector< std::string >& columns, const std::string& sep = ";" ) const;

        std::ios_base::fmtflags format ( ) const;
        void set_format ( std::ios_base::fmtflags format );

        void set_precision ( int precision );

      private:
        typedef std::map< std::string, TableColumn >::const_iterator ColumnIterator;

        static const std::ios_base::fmtflags DEFAULT_FORMAT;

        // Returns iterators to the requested columns in the same order.
        void find_columns ( const std::vector< std::string >& columns,
                            std::vector< ColumnIterator >& col_iters,
                            int& num_rows ) const;

        // Returns maximum width of any column, but at least min_width.
        int compute_column_width ( const std::vector< ColumnIterator >& col_iters, int num_rows, int min_width ) const;
        void print_divider ( std::ostream& os, int width ) const;

        std::map< std::string, TableColumn > columns_;
        std::ios_base::fmtflags format_;
        int precision_;
    };

    template< class ValueType >
    void Table::insert ( const std::string& column_name, const ValueType& value )
    {
        std::ostringstream sstr;
        sstr << std::setiosflags ( format_ ) << std::setprecision ( precision_ ) << value;
        columns_[column_name].insert ( sstr.str ( ) );
    }

}

#endif
