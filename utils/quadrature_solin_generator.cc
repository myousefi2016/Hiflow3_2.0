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

/// Generate Quadrature classes for HiFlow based on Database in book
/// by Pavel Solin.
///

/// \author Staffan Ronnas

/// This utility reads in quadrature rules from a database in the
/// format used in the file values.txt that accompanies the above
/// mentioned book, and generates a quadrature class for HiFlow3 based
/// on the values in the file. The file has to be split up into
/// separate files, one for each type of quadrature that one wants to
/// extract.
///
/// Since the quadrature rules in the Solin book use different
/// reference elements than those in HiFlow3, the values have to be
/// transformed.
///
/// Transformation used: eta_coord_hiflow = (eta_coord_solin - 1) / 2, eta = x, y, z
///                      weight_hiflow = weight_solin / (2 ^ dimension)
///
/// the triangle, where x and y are transformed as
/// x_hiflow = (x_solin + 1) / 2 and the weights as
/// w_hiflow = w_solin / 4.
///
/// The same transformation can be used with all cell types, except the prism.
///
/// Since the text database that comes with book only has 15 decimals,
/// the values here are only accurate to about that number of
/// decimals, and hence no more than so are printed. The extended type
/// long double (128-bit on linux with gcc) is used for reading and
/// manipulating values, so this should not cause any significant loss
/// of precision.
///
/// Note: the extraction of the order information will not work
/// correctly for the tensor-product cells, which have two orders in
/// the original file. To fix this, remove the smaller order from the
/// input, so that it says for instance 'p = 5'.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

struct Quadrature
{
    std::vector< std::vector<long double> > coords;
    std::vector<long double> weights;
    int order;
};

void print_header_code ( std::ostream& os );
void print_copyright ( std::ostream& os );
void print_source_code ( std::ostream& os, const std::vector<Quadrature>& quadratures );
std::vector<int> compute_order_size_map ( const std::vector<Quadrature>& quadratures );

// Name of the quadrature (input argument 1)
const char* g_name;

// Dimension of the quadrature (input argument 2)
int g_dimension;

std::string filename_base;

int main ( int argc, char** argv )
{

    if ( argc < 3 )
    {
        std::cerr << "Usage: quadrature_solin_generator name dimension < quadrature_file\n";
        return 1;
    }

    // Read in arguments
    g_name = argv[1];

    g_dimension = atoi ( argv[2] );
    std::cerr << "Dimension = " << g_dimension << "\n";
    std::cerr << "sizeof(long double) = " << sizeof (long double ) << "\n"
            << "sizeof(double) = " << sizeof (double ) << "\n";

    int num_quadratures = 0;

    std::vector<Quadrature> quadratures;

    long double val;
    std::string comment;

    // Read quadrature rules
    while ( std::cin )
    {
        int num_points = -1;

        std::cin >> num_points;

        std::cerr << "Reading quadrature with " << num_points << " points.\n";

        quadratures.push_back ( Quadrature ( ) );
        Quadrature& curr = quadratures.back ( );
        curr.coords.resize ( g_dimension );

        for ( int p = 0; p < num_points; ++p )
        {

            for ( int i = 0; i < g_dimension; ++i )
            {
                std::cin >> val;
                curr.coords[i].push_back ( ( val + 1. ) / 2. );
            }

            std::cin >> val;
            curr.weights.push_back ( val / std::pow ( 2., g_dimension ) );
            getline ( std::cin, comment ); // ignore end of line
        }

        getline ( std::cin, comment ); // get comment line

        // extract order: search for '=', and take the character directly after
        int pos = comment.find ( "=" );
        curr.order = atoi ( &comment[pos + 1] );

        getline ( std::cin, comment ); // ignore blank line

        ++num_quadratures;
    }

    std::string class_name_lower = std::string ( g_name );
    std::transform ( class_name_lower.begin ( ),
                     class_name_lower.end ( ),
                     class_name_lower.begin ( ), ::tolower );

    // Open output files
    filename_base = std::string ( "q" ) + class_name_lower;

    std::ofstream header_file ( ( filename_base + std::string ( ".h" ) ).c_str ( ) );
    std::ofstream source_file ( ( filename_base + std::string ( ".cc" ) ).c_str ( ) );

    // Print header file.
    print_header_code ( header_file );

    // Print source file.
    print_source_code ( source_file, quadratures );

    return 0;
}

void print_header_code ( std::ostream& os )
{
    print_copyright ( os );

    std::string class_name = std::string ( "Quadrature" ) + std::string ( g_name );
    std::string class_name_upper = class_name;
    std::transform ( class_name_upper.begin ( ),
                     class_name_upper.end ( ),
                     class_name_upper.begin ( ), ::toupper );

    os << "#ifndef HIFLOW_QUADRATURE_" << class_name_upper << "_H\n"
            << "#define HIFLOW_QUADRATURE_" << class_name_upper << "_H\n"
            << "\n"
            << "#include \"quadraturetype.h\"\n"
            << "\n"
            << "namespace hiflow {\n"
            << "\n"
            << "template<class DataType>\n"
            << "class " << class_name << " : public QuadratureType<DataType>\n"
            << "{\n"
            << "public:\n"
            << "  " << class_name << "();\n"
            << "  " << class_name << "<DataType>* clone() const;\n"
            << "\n"
            << "  using QuadratureType<DataType>::x_;\n"
            << "  using QuadratureType<DataType>::y_;\n"
            << "  using QuadratureType<DataType>::z_;\n"
            << "  using QuadratureType<DataType>::w_;\n"
            << "  using QuadratureType<DataType>::index_field_;\n"
            << "  using QuadratureType<DataType>::order_size_map_;\n"
            << "};\n"
            << "} // namespace hiflow\n"
            << "#endif\n";
}

void print_source_code ( std::ostream& os, const std::vector<Quadrature>& quadratures )
{

    // Only 15 decimals are printed, since the original data set does
    // not have more precision than this.
    os << std::setprecision ( 15 );
    print_copyright ( os );

    std::string class_name = std::string ( "Quadrature" ) + std::string ( g_name );

    const char coord_char[3] = { 'x', 'y', 'z' };

    os << "#include \"" << filename_base << ".h\"\n"
            << "\n"
            << "namespace hiflow {\n"
            << "\n"
            << "template<class DataType>\n"
            << class_name << "<DataType>::" << class_name << "()\n"
            << "{\n"
            << "  ////////////////////////////////////////////////\n"
            << "  // Implemented quadratures:\n";

    const int num_quadratures = quadratures.size ( );

    for ( int i = 0; i < num_quadratures; ++i )
    {
        os << "  // Size = " << quadratures.at ( i ).weights.size ( )
                << ", order = " << quadratures.at ( i ).order << "\n";
    }
    os << "  ////////////////////////////////////////////////\n"
            << "  const int total_number_of_quadratures = " << num_quadratures << ";\n"
            << "\n"
            << "  int cntr = 0;\n"
            << "  int size = 0;\n"
            << "\n"
            << "  // Resize vector fields.\n"
            << "  x_.resize(total_number_of_quadratures);\n"
            << "  y_.resize(total_number_of_quadratures);\n"
            << "  z_.resize(total_number_of_quadratures);\n"
            << "  w_.resize(total_number_of_quadratures);\n"
            << "\n";

    std::vector<int> order_size_map = compute_order_size_map ( quadratures );
    os << "  // Create order->size map.\n"
            << "  order_size_map_.clear();\n"
            << "  order_size_map_.resize(" << order_size_map.size ( ) + 1 << ");\n";
    for ( int i = 0; i < static_cast < int > ( order_size_map.size ( ) ); ++i )
    {
        os << "  order_size_map_.at(" << i << ") = " << order_size_map.at ( i ) << ";\n";
    }
    os << "\n";

    os << "  // Fill vector fields with quadrature data.\n"
            << "\n";
    for ( int i = 0; i < num_quadratures; ++i )
    {
        const int quad_size = quadratures.at ( i ).weights.size ( );

        os << "  // -- Size = " << quad_size
                << ", Order = " << quadratures.at ( i ).order << " -----------------\n"
                << "\n"
                << "  size = " << quadratures.at ( i ).weights.size ( ) << ";\n"
                << "\n";

        for ( int d = 0; d < g_dimension; ++d )
        {
            os << "  DataType " << coord_char[d] << i << "[] = {\n";
            for ( int p = 0; p < quad_size; ++p )
            {
                os << "    " << quadratures.at ( i ).coords.at ( d ).at ( p );
                if ( p < ( quad_size - 1 ) )
                {
                    os << ", \n";
                }
                else
                {
                    os << "\n";
                }
            }
            os << "  };\n";
        }

        os << "  DataType w" << i << "[] = {\n";
        for ( int p = 0; p < quad_size; ++p )
        {
            os << "    " << quadratures.at ( i ).weights.at ( p );
            if ( p < ( quad_size - 1 ) )
            {
                os << ", \n";
            }
            else
            {
                os << "\n";
            }
        }
        os << "  };\n"
                << "\n";

        for ( int d = 0; d < g_dimension; ++d )
        {
            os << "  " << coord_char[d] << "_[cntr].reserve(size); "
                    << coord_char[d] << "_[cntr].insert(" << coord_char[d]
                    << "_[cntr].begin(), " << coord_char[d] << i << ", "
                    << coord_char[d] << i << " + size);\n";
        }

        for ( int d = g_dimension; d < 3; ++d )
        {
            os << "  " << coord_char[d] << "_[cntr].resize(size, 0.);\n";
        }

        os << "  " << "w_[cntr].reserve(size); "
                << "w_[cntr].insert(w_[cntr].begin(), w" << i
                << ", w" << i << " + size);\n";

        os << "\n"
                << "  index_field_[size] = cntr;\n"
                << "  ++cntr;\n"
                << "\n";
        os << "  // ----------------------------------------------------\n"
                << "\n";
    }

    os << "}\n"
            << "\n"
            << "template<class DataType>\n"
            << class_name << "<DataType>* " << class_name << "<DataType>::clone() const\n"
            << "{\n"
            << "  return new " << class_name << "<DataType>(*this);\n"
            << "}\n"
            << "\n"
            << "\n"
            << "// Template instantiation.\n"
            << "template class " << class_name << "<double>;\n"
            << "template class " << class_name << "<float>;\n"
            << "\n"
            << "} // namespace hiflow\n";
}

void print_copyright ( std::ostream& os )
{
    os << "// Copyright (C) 2011-2015 Vincent Heuveline\n"
            << "//\n"
            << "// HiFlow3 is free software: you can redistribute it and/or modify it under the\n"
            << "// terms of the GNU Lesser General Public License as published by the Free\n"
            << "// Software Foundation, either version 3 of the License, or (at your option) any\n"
            << "// later version.\n"
            << "// \n"
            << "// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY\n"
            << "// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR\n"
            << "// A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more\n"
            << "// details.\n"
            << "//\n"
            << "// You should have received a copy of the GNU Lesser General Public License\n"
            << "// along with HiFlow3.  If not, see <http://www.gnu.org/licenses/>.\n"
            << "\n"
            << "\n"
            << "/// \\author Code generated with program quadrature_solin_generator.\n"
            << "\n"
            << "// This code was generated with quadrature_solin_generator by Staffan Ronnas.\n"
            << "// The quadrature data comes from the book \n"
            << "// P. Solin, et. al., \"Higher-Order Finite Element Methods\"\n"
            << "// Michael Schick developed the original class on which the generated code is based.\n"
            << "\n"
            << "\n";
}

void print_debug ( const std::vector<Quadrature>& quadratures )
{
    // Print out
    for ( int i = 0; i < static_cast < int > ( quadratures.size ( ) ); ++i )
    {
        const int num_p = quadratures[i].weights.size ( );
        std::cout << std::setprecision ( 20 );
        std::cout << "Quadrature of order " << quadratures[i].order << " with " << num_p << " points\n";

        for ( int d = 0; d < g_dimension; ++d )
        {
            std::copy ( quadratures[i].coords[d].begin ( ),
                        quadratures[i].coords[d].end ( ), std::ostream_iterator<long double>( std::cout, " " ) );
            std::cout << "\n";
        }
        std::copy ( quadratures[i].weights.begin ( ),
                    quadratures[i].weights.end ( ), std::ostream_iterator<long double>( std::cout, " " ) );
        std::cout << "\n\n";
    }
}

std::vector<int> compute_order_size_map ( const std::vector<Quadrature>& quadratures )
{
    const int max_order = quadratures.back ( ).order;
    std::vector<int> order_size_map ( max_order + 1, -1 );
    int i = 0;
    for ( int p = 0; p < max_order + 1; ++p )
    {
        while ( p > quadratures.at ( i ).order )
        {
            ++i;
        }
        order_size_map.at ( p ) = quadratures.at ( i ).weights.size ( );
    }

    return order_size_map;
}
