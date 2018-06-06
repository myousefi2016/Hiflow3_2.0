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

#include "table.h"

#include <iomanip>
#include <numeric>

namespace hiflow
{

    const std::string TableColumn::MISSING_STRING = "";

    void TableColumn::insert ( const std::string& str )
    {
        values_.push_back ( str );
    }

    int TableColumn::num_rows ( ) const
    {
        return values_.size ( );
    }

    const std::string& TableColumn::value ( int row ) const
    {
        if ( row >= 0 && row < num_rows ( ) )
        {
            return values_[row];
        }
        else
        {
            return MISSING_STRING;
        }
    }

    const std::ios_base::fmtflags Table::DEFAULT_FORMAT =
            std::ios_base::left |
            std::ios_base::boolalpha |
            std::ios_base::showpoint |
            std::ios_base::dec |
            std::ios_base::fixed;

    Table::Table ( )
    : format_ ( Table::DEFAULT_FORMAT ), precision_ ( 6 )
    {
    }

    void Table::print ( std::ostream& os, const std::vector< std::string >& columns ) const
    {

        std::vector< ColumnIterator > col_iters;
        int num_rows = 0;
        find_columns ( columns, col_iters, num_rows );

        const int num_cols = col_iters.size ( );
        const int col_width = compute_column_width ( col_iters, num_rows, 10 );
        const int table_width = num_cols * ( col_width + 1 );

        print_divider ( os, table_width );

        for ( int i = 0; i != num_cols - 1; ++i )
        {
            os << std::setw ( col_width ) << col_iters[i]->first << " ";
        }
        os << std::setw ( col_width ) << col_iters[col_iters.size ( ) - 1]->first;
        os << "\n";

        print_divider ( os, table_width );

        for ( int r = 0; r != num_rows; ++r )
        {
            for ( int c = 0; c != num_cols - 1; ++c )
            {
                os << std::setw ( col_width ) << col_iters[c]->second.value ( r ) << " ";
            }
            os << std::setw ( col_width ) << col_iters[num_cols - 1]->second.value ( r );
            os << "\n";
        }

        print_divider ( os, table_width );
    }

    void Table::print_csv ( std::ostream& os, const std::vector< std::string >& columns, const std::string& sep ) const
    {
        std::vector< ColumnIterator > col_iters;
        int num_rows = 0;
        find_columns ( columns, col_iters, num_rows );

        const int num_cols = col_iters.size ( );

        // Print titles.
        for ( int i = 0; i != num_cols - 1; ++i )
        {
            os << col_iters[i]->first << sep;
        }
        os << col_iters[num_cols - 1]->first << "\n";

        // Print rows.
        for ( int r = 0; r != num_rows; ++r )
        {
            for ( int c = 0; c != num_cols - 1; ++c )
            {
                os << col_iters[c]->second.value ( r ) << sep;
            }
            os << col_iters[num_cols - 1]->second.value ( r ) << "\n";
        }
    }

    std::ios_base::fmtflags Table::format ( ) const
    {
        return format_;
    }

    void Table::set_format ( std::ios_base::fmtflags format )
    {
        format_ = format;
    }

    void Table::set_precision ( int precision )
    {
        precision_ = precision;
    }

    void Table::find_columns ( const std::vector< std::string >& columns,
                               std::vector< ColumnIterator >& col_iters,
                               int& num_rows ) const
    {
        col_iters.reserve ( columns.size ( ) );

        num_rows = 0;
        for ( std::vector< std::string >::const_iterator it = columns.begin ( );
              it != columns.end ( ); ++it )
        {
            ColumnIterator col_it = columns_.find ( *it );
            if ( col_it != columns_.end ( ) )
            {
                col_iters.push_back ( col_it );
                num_rows = std::max ( num_rows, col_it->second.num_rows ( ) );
            }
        }
    }

    int Table::compute_column_width ( const std::vector< ColumnIterator >& col_iters, int num_rows, int min_width ) const
    {
        // determine column width
        int col_width = min_width; // minimum column width
        const int num_cols = col_iters.size ( );

        for ( int i = 0; i != num_cols; ++i )
        {
            // All columns will be same width as maximum length of names.
            col_width = std::max ( col_width, int(col_iters[i]->first.length ( ) ) );
        }

        for ( int r = 0; r != num_rows; ++r )
        {
            for ( int c = 0; c != num_cols; ++c )
            {
                col_width = std::max ( col_width, int(col_iters[c]->second.value ( r ).length ( ) ) );
            }
        }
        return col_width;
    }

    void Table::print_divider ( std::ostream& os, int width ) const
    {
        os << std::setfill ( '-' ) << std::setw ( width ) << "\n"
                << std::setfill ( ' ' );
    }
}
