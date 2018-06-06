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
/// @date 2015-09-29

#ifndef SRC_LINEAR_SOLVER_AMG_DIRECTINTERPOLATION_H_
#    define SRC_LINEAR_SOLVER_AMG_DIRECTINTERPOLATION_H_

#    include "cxx98_prettyprint.h"
#    include "InitFunctor.h"
#    include "linear_algebra/lmp/lmatrix_coo_cpu.h"
#    include "LevelGenerator.h"
#    include "StandardCoarseningInterpolationConnection.h"
#    include "typedefs.h"
#    include "TransposeFunctor.h"
#    include <algorithm>
#    include <map>

//#define PRINT_DEBUG

namespace hiflow
{
    namespace AMG
    {

        /**
         * \brief DirectInterpolation
         *
         * Reference: K. Stueben, Algebraic Multigrid (AMG): An Introduction with Applications
         *            GMD-Report 70, November 1999
         */
        template <class LevelGenerator>
        class DirectInterpolation
        {
          public:

            typedef typename LevelGenerator::MatrixType MatrixType;
            typedef typename LevelGenerator::PtrMatrixType PtrMatrixType;
            typedef typename LevelGenerator::ResultType ResultType;
            typedef typename MatrixType::value_type value_type;
            typedef StandardCoarseningInterpolationConnection ConnectionType;
            typedef NoSettings Settings;
            typedef NoOutput Output;

            DirectInterpolation ( Settings const& settings = Settings ( ) )
            : settings_ ( settings )
            {
            }

            void operator() ( ResultType &result, ConnectionType const& connection, MatrixType const& Af ) const
            {
                size_t coarse_dim = connection.coarse_variables.size ( );
                size_t fine_dim = connection.fine_variables.size ( );
                size_t omega_dim = coarse_dim + fine_dim;

                // Map coarse indices of Af to Ac
                std::map<size_t, size_t> map_coarse_indices;
                {
                    // Workaround, since in C++ it is not possible to initialize different types within a for-loop
                    size_t pos ( 0 );
                    for ( IndexSet::const_iterator iter_cur ( connection.coarse_variables.begin ( ) ),
                          iter_end ( connection.coarse_variables.end ( ) ); iter_cur != iter_end; ++iter_cur, ++pos )
                    {
                        map_coarse_indices[*iter_cur] = pos;
                    }
                }

                // Determine sparse structure of interpolation matrix
                IndexVector rows, cols;
                std::vector<value_type> values;

                // ... fine variables
                for ( IndexSet::const_iterator iter_fine_cur ( connection.fine_variables.begin ( ) ),
                      iter_fine_end ( connection.fine_variables.end ( ) ); iter_fine_cur != iter_fine_end; ++iter_fine_cur )
                {
                    // Find diag element aii and calculate sum_Ei
                    value_type aii ( 0 );
                    value_type sum_Ei ( 0 );
                    for ( int p ( Af.matrix.row[*iter_fine_cur] ); p != Af.matrix.row[*iter_fine_cur + 1]; ++p )
                    {
                        if ( Af.matrix.col[p] == *iter_fine_cur )
                        {
                            aii = Af.matrix.val[p];
                        }
                        else
                        {
                            sum_Ei += Af.matrix.val[p];
                        }
                    }

                    // Strongly connected coarse grid points
                    IndexSet Ci;
                    std::set_intersection ( connection.strong_connections[*iter_fine_cur].begin ( ), connection.strong_connections[*iter_fine_cur].end ( ),
                                            connection.coarse_variables.begin ( ), connection.coarse_variables.end ( ),
                                            std::inserter ( Ci, Ci.end ( ) ) );

                    // Calculate sum_Ci
                    value_type sum_Ci ( 0 );
                    for ( int p ( Af.matrix.row[*iter_fine_cur] ); p != Af.matrix.row[*iter_fine_cur + 1]; ++p )
                    {
                        if ( Ci.find ( Af.matrix.col[p] ) != Ci.end ( ) )
                        {
                            sum_Ci += Af.matrix.val[p];
                        }
                    }

                    // Embrace j-independent factor by dividing the correction factor alpha = sum_Ei(aik) / sum_Ci(aik)
                    // by diagonal element aii and change the sign
                    value_type factor = -sum_Ei / ( aii * sum_Ci );

                    // Find possible coarse variables
                    IndexVector active_coarse;
                    std::set_intersection ( Af.matrix.col + Af.matrix.row[*iter_fine_cur], Af.matrix.col + Af.matrix.row[*iter_fine_cur + 1],
                                            connection.coarse_variables.begin ( ), connection.coarse_variables.end ( ),
                                            std::back_inserter ( active_coarse ) );

#    ifdef PRINT_DEBUG
                    std::cout << "Af.matrix.row = ";
                    for ( size_t pos = Af.matrix.row[*iter_fine_cur]; pos != Af.matrix.row[*iter_fine_cur + 1]; ++pos ) std::cout << Af.matrix.col[pos] << " ";
                    std::cout << std::endl;
                    std::cout << "connection.coarse_variables = ";
                    for ( IndexSet::const_iterator iter_cur = connection.coarse_variables.begin ( ); iter_cur != connection.coarse_variables.end ( ); ++iter_cur ) std::cout << *iter_cur << " ";
                    std::cout << std::endl;
                    std::cout << "active_coarse = ";
                    for ( IndexVector::const_iterator iter_cur = active_coarse.begin ( ); iter_cur != active_coarse.end ( ); ++iter_cur ) std::cout << *iter_cur << " ";
                    std::cout << std::endl;
#    endif

                    for ( IndexVector::const_iterator iter_coarse_cur ( active_coarse.begin ( ) ),
                          iter_coarse_end ( active_coarse.end ( ) ); iter_coarse_cur != iter_coarse_end; ++iter_coarse_cur )
                    {
                        // Find matrix element aij
                        value_type aij ( 0 );
                        for ( int p ( Af.matrix.row[*iter_fine_cur] ); p != Af.matrix.row[*iter_fine_cur + 1]; ++p )
                        {
                            if ( Af.matrix.col[p] == *iter_coarse_cur )
                            {
                                aij = Af.matrix.val[p];
                                break;
                            }
                        }

                        rows.push_back ( *iter_fine_cur );
                        cols.push_back ( map_coarse_indices[*iter_coarse_cur] );
                        values.push_back ( factor * aij );
                    }
                }

                // ... coarse variables
                for ( IndexSet::const_iterator iter_coarse_cur ( connection.coarse_variables.begin ( ) ),
                      iter_coarse_end ( connection.coarse_variables.end ( ) ); iter_coarse_cur != iter_coarse_end; ++iter_coarse_cur )
                {
                    rows.push_back ( *iter_coarse_cur );
                    cols.push_back ( map_coarse_indices[*iter_coarse_cur] );
                    values.push_back ( 1.0 );
                }

#    ifdef PRINT_DEBUG
                std::cout << "rows = ";
                for ( size_t i ( 0 ); i != rows.size ( ); ++i ) std::cout << rows[i] << " ";
                std::cout << std::endl;
                std::cout << "cols = ";
                for ( size_t i ( 0 ); i != cols.size ( ); ++i ) std::cout << cols[i] << " ";
                std::cout << std::endl;
                std::cout << "values = ";
                for ( size_t i ( 0 ); i != values.size ( ); ++i ) std::cout << values[i] << " ";
                std::cout << std::endl;
#    endif

                // Copy sparse structure to temporary coo matrix
                CPUsimple_COO_lMatrix<value_type> coo;
                coo.Init ( values.size ( ), omega_dim, coarse_dim, "Pocc" );
                std::copy ( rows.begin ( ), rows.end ( ), coo.matrix.row );
                std::copy ( cols.begin ( ), cols.end ( ), coo.matrix.col );
                std::copy ( values.begin ( ), values.end ( ), coo.matrix.val );

                // Copy coo matrix to final interpolation matrix
                result.ptr_interpolation_matrix = PtrMatrixType ( new MatrixType );
                result.ptr_interpolation_matrix->ConvertFrom ( coo );

                // Calculate restriction matrix as transposed interpolation matrix
                result.ptr_restriction_matrix = PtrMatrixType ( new MatrixType );
                InitFunctor<MatrixType>( )( *result.ptr_restriction_matrix, omega_dim * coarse_dim, omega_dim, coarse_dim, "R", Af );
                result.ptr_restriction_matrix->CopyStructureFrom ( *result.ptr_interpolation_matrix );
                result.ptr_restriction_matrix->CopyFrom ( *result.ptr_interpolation_matrix );
                TransposeFunctor<MatrixType>( )( *result.ptr_restriction_matrix );
            }

          private:

            Settings settings_;

        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_DIRECTINTERPOLATION_H_ */
