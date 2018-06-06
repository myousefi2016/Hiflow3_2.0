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

#ifndef HIFLOW_ADAPTIVITY_REFINEMENT_STRATEGIES
#    define HIFLOW_ADAPTIVITY_REFINEMENT_STRATEGIES

/// \author Philipp Gerstner

#    include <map>
#    include <string>
#    include <vector>

namespace hiflow
{
    ///
    /// \brief set of functions for building mesh adaption flags out of local error indicators
    /// 

    /// \brief refine / coarsen a fixed fraction of cells on local processor
    /// @param[in] refine_frac refine the refince_frac% cells with largest indicators
    /// @param[in] coarsen_frac coarsen the rcoarsen_frac% cells with smallest indicators
    /// @param[in] threshold start coarsening when more than #threshold cells are present in mesh 
    /// @param[in] coarsen_marker coarsen family l of cells k if \sum_{k \in Family(l)} coarsen_marker(k) <= #Family(l)
    /// @param[in] indicators error indicators
    /// @param[out] adapt_markers resulting adaption markers to be passed to refine routine of mesh 
    void local_fixed_fraction_strategy ( double refine_frac,
                                         double coarsen_frac,
                                         int threshold,
                                         int coarsen_marker,
                                         const std::vector<double>& indicators,
                                         std::vector<int>& adapt_markers );

    /// \brief adapt mesh to obtain a desired accuracy: refine cell k if indicator(k) > tol / num_global_cells \br
    /// coarsen cell k if 2^conv_order * indicator(k) < tol / num_global_cells  
    /// @param[in] tol desired accuracy
    /// @param[in] num_global_cells number of global cells in mesh
    /// @param[in] conv_order assumed convergence order of applied numerical scheme
    /// @param[in] threshold start coarsening when more than #threshold cells are present in mesh 
    /// @param[in] coarsen_marker coarsen family l of cells k if \sum_{k \in Family(l)} coarsen_marker(k) <= #Family(l)
    /// @param[in] indicators error indicators
    /// @param[out] adapt_markers resulting adaption markers to be passed to refine routine of mesh 
    void fixed_error_strategy ( double tol,
                                int num_global_cells,
                                double conv_order,
                                int threshold,
                                int coarsen_marker,
                                const std::vector<double>& indicators,
                                std::vector<int>& adapt_markers );
}
#endif
