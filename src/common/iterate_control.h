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

#ifndef HIFLOW_LINEARSOLVER_ITERATE_CONTROL_H_
#    define HIFLOW_LINEARSOLVER_ITERATE_CONTROL_H_

namespace hiflow
{

    /// @brief Class for test on convergence of iterative methods
    /// @author Hendryk Bockelmann, Chandramowli Subramanian
    ///
    /// Test is done by calling the function @em Check.

    class IterateControl
    {
      public:

        enum State
        {
            kIterate = 0,
            kSuccessAbsoluteTol,
            kSuccessRelativeTol,
            kFailureDivergenceTol,
            kFailureMaxitsExceeded
        };

        IterateControl ( );
        IterateControl ( int maxits, double atol, double rtol, double dtol );

        virtual ~IterateControl ( )
        {
        }

        virtual void Init ( int maxits, double atol, double rtol, double dtol );
        virtual State Check ( int current_step, double current_value );

        void set_first ( int first_step, double first_value );

        double first_value ( ) const
        {
            return first_value_;
        }

        void set_first_value ( double first_value )
        {
            this->first_value_ = first_value;
        }

        int first_step ( ) const
        {
            return first_step_;
        }

        int maxits ( ) const
        {
            return maxits_;
        }

        double absolute_tol ( ) const
        {
            return absolute_tol_;
        }

        double relative_tol ( ) const
        {
            return relative_tol_;
        }

        double divergence_tol ( ) const
        {
            return divergence_tol_;
        }

        State status ( ) const
        {
            return status_;
        }

        void set_status ( State status )
        {
            this->status_ = status;
        }

      protected:
        int maxits_;
        double absolute_tol_;
        double relative_tol_;
        double divergence_tol_;

        double first_value_;
        int first_step_;

        State status_;
    };

} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_ITERATE_CONTROL_H_
