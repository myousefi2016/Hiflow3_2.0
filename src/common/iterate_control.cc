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

#include <iostream>
#include <cassert>
#include "iterate_control.h"

/// @author Hendryk Bockelmann, Chandramowli Subramanian

namespace hiflow
{

    /// standard constructor

    IterateControl::IterateControl ( )
    {
        this->Init ( 1000, 1.e-12, 1.e-8, 1.e8 );
        this->set_first ( 0, 0. );
    }

    /// constructor with parameters
    /// @param maxits maximum number of iteration steps
    /// @param atol absolute tolerance of residual to converge
    /// @param rtol relative tolerance of residual to converge
    /// @param dtol relative tolerance of residual to diverge

    IterateControl::IterateControl ( int maxits, double atol,
                                     double rtol, double dtol )
    {
        this->Init ( maxits, atol, rtol, dtol );
    }

    /// Init with parameters.
    /// @param maxits maximum number of iteration steps
    /// @param atol absolute tolerance of residual to converge
    /// @param rtol relative tolerance of residual to converge
    /// @param dtol relative tolerance of residual to diverge

    void IterateControl::Init ( int maxits, double atol, double rtol, double dtol )
    {
        assert ( maxits > 0 );
        assert ( atol >= 0. );
        assert ( rtol >= 0. && rtol < 1. );
        assert ( dtol > 1. );

        this->maxits_ = maxits;
        this->absolute_tol_ = atol;
        this->relative_tol_ = rtol;
        this->divergence_tol_ = dtol;
    }

    /// Set first_step and first_value.
    /// @param first_step first iteration step (usually zero)
    /// @param first_value first residuum

    void IterateControl::set_first ( int first_step, double first_value )
    {
        this->first_step_ = first_step;
        this->first_value_ = first_value;
    }

    /// Check new status, i.e. if iterative method has converged, diverged
    /// or needs to be iterated.
    /// If value is less than absolute tolerance or value/first_value is less,
    /// than relative tolerance, the iterative method succeeded.
    /// If value/first_value is larger than divergence tolerance or
    /// current step is larger than maximum number of iterations, the iterative
    /// method failed. Else it has to be iterated once more.
    /// @param current_step current step
    /// @param current_value current residuum
    /// @return status of iterative method

    IterateControl::State IterateControl::Check ( int current_step,
                                                  double current_value )
    {
        assert ( current_step >= this->first_step ( ) );

        // save first residuum
        if ( current_step == this->first_step ( ) )
            this->set_first_value ( current_value );

        // check
        if ( current_value < this->absolute_tol ( ) )
            this->set_status ( kSuccessAbsoluteTol );

        else if ( current_value < this->relative_tol ( ) * this->first_value ( ) )
            this->set_status ( kSuccessRelativeTol );

        else if ( current_value > this->divergence_tol ( ) * this->first_value ( ) )
            this->set_status ( kFailureDivergenceTol );

        else if ( current_step >= this->maxits ( ) )
            this->set_status ( kFailureMaxitsExceeded );

        else
            this->set_status ( kIterate );

        return this->status ( );
    }

} // namespace hiflow
