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

#include "refinement_strategies.h"
#include "common/log.h"
#include "common/sort_permutation.h"
#include <cmath>
#include <utility> 

const int DEBUG_LEVEL = 0;

namespace hiflow
{

    void local_fixed_fraction_strategy ( double refine_frac,
                                         double coarsen_frac,
                                         int threshold,
                                         int coarsen_marker,
                                         const std::vector<double>& indicators,
                                         std::vector<int>& adapt_markers )
    {
        adapt_markers.resize ( indicators.size ( ), 0 );
        std::vector<int> sort_ind ( indicators.size ( ), 0 );

        for ( int i = 0; i < indicators.size ( ); ++i )
        {
            sort_ind[i] = i;
        }
        sortingPermutation ( indicators, sort_ind );

        // 1. Mark cells for refinement
        int first_cell = std::floor ( ( 1. - refine_frac ) * sort_ind.size ( ) );
        int num_global_cells = indicators.size ( );

        for ( int l = first_cell; l < sort_ind.size ( ); ++l )
        {
            int cell_index = sort_ind[l];
            adapt_markers[cell_index] = 1;
        }

        // 2.Mark cells for coarsening
        if ( num_global_cells >= threshold )
        {
            int last_cell = std::ceil ( coarsen_frac * sort_ind.size ( ) );
            for ( int l = 0; l < last_cell; ++l )
            {
                int cell_index = sort_ind[l];
                adapt_markers[cell_index] = coarsen_marker;
            }
        }
    }

    void fixed_error_strategy ( double tol,
                                int num_global_cells,
                                double conv_order,
                                int threshold,
                                int coarsen_marker,
                                const std::vector<double>& indicators,
                                std::vector<int>& adapt_markers )
    {
        adapt_markers.resize ( indicators.size ( ), 0 );
        double av_max_error = tol / num_global_cells;
        int num_cells = indicators.size ( );

        for ( int c = 0; c < num_cells; ++c )
        {
            // mark cells for refinement
            if ( indicators[c] > av_max_error )
            {
                adapt_markers[c] = 1;
            }

            // mark cells for coarsening
            if ( num_global_cells >= threshold )
            {
                if ( indicators[c] * std::pow ( 2.0, conv_order ) < av_max_error )
                {
                    adapt_markers[c] = coarsen_marker;
                }
            }
        }
    }
}
