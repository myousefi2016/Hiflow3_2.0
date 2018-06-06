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

#include "forcing_eisenstat_walker.h"
#include <cmath>

namespace hiflow
{

    template<class LAD>
    void EWForcing<LAD>::ComputeForcingTerm ( DataType new_residual, DataType lin_solve_accuracy )
    {
        if ( Choice_ == 1 )
        {
            ComputeForcingChoice1 ( new_residual, lin_solve_accuracy );
        }
        else if ( Choice_ == 2 )
        {
            ComputeForcingChoice2 ( new_residual );
        }
        else
        {
            ComputeForcingChoice1 ( new_residual, lin_solve_accuracy );
        }
    }

    template<class LAD>
    void EWForcing<LAD>::ComputeForcingChoice1 ( DataType new_residual, DataType lin_solve_accuracy )
    {
        DataType old_residual_norm = residuals_.back ( );

        DataType eta = std::fabs ( new_residual - lin_solve_accuracy ) / old_residual_norm;

        if ( eta > Maximal_ )
            eta = Maximal_;

        forcing_terms_.push_back ( eta );
        residuals_.push_back ( new_residual );
    }

    template<class LAD>
    void EWForcing<LAD>::ComputeForcingChoice2 ( DataType new_residual )
    {
        DataType old_residual = residuals_.back ( );

        DataType eta = gamma_ * pow ( ( new_residual / old_residual ), alpha_ );

        if ( eta > Maximal_ )
            eta = Maximal_;

        forcing_terms_.push_back ( eta );
        residuals_.push_back ( new_residual );
    }

    template<class LAD>
    EWForcing<LAD>::EWForcing ( )
    {
        Initial_ = 0.5;
        Maximal_ = 0.9;
        Choice_ = 1;
        name_ = "EisenstatWalker";
        gamma_ = 1.0;
        alpha_ = 0.5 * ( 1. + sqrt ( 5.0 ) );
        forcing_terms_.push_back ( Initial_ );
    }

    template<class LAD>
    EWForcing<LAD>::EWForcing ( int choice )
    {
        Initial_ = 0.5;
        Maximal_ = 0.9;
        Choice_ = choice;
        name_ = "EisenstatWalker";
        forcing_terms_.push_back ( Initial_ );
        gamma_ = 1.0;
        alpha_ = 0.5 * ( 1. + sqrt ( 5.0 ) );
    }

    template<class LAD>
    EWForcing<LAD>::EWForcing ( double initial, double max, int choice ) :
    Initial_ ( initial ), Maximal_ ( max ), Choice_ ( choice )
    {
        name_ = "EisenstatWalker";
        forcing_terms_.push_back ( Initial_ );
    }

    // constructor

    template<class LAD>
    EWForcing<LAD>::EWForcing ( double initial, double max, int choice, double gamma, double alpha ) :
    Initial_ ( initial ), Maximal_ ( max ), Choice_ ( choice ), alpha_ ( alpha ), gamma_ ( gamma )
    {
        name_ = "EisenstatWalker";
        forcing_terms_.push_back ( Initial_ );
    }

    template<class LAD>
    EWForcing<LAD>::~EWForcing ( )
    {
    }

    /// template instantiation
    template class EWForcing<la::LADescriptorCoupledD>;
    template class EWForcing<la::LADescriptorCoupledS>;
    template class EWForcing<la::LADescriptorHypreD>;
    template class EWForcing<la::LADescriptorPolynomialChaosD>;
    template class EWForcing<la::LADescriptorPolynomialChaosExpansionD>;

} // namespace hiflow
