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

#ifndef HIFLOW_NONLINEAR_FORCING_STRATEGY_H_
#    define HIFLOW_NONLINEAR_FORCING_STRATEGY_H_

#    include <string>

namespace hiflow
{

    /// @brief Base class for forcing strategies
    /// @author Michael Schick
    ///
    /// A forcing strategy determines the relative tolerance
    /// for convergence within a linear system of a
    /// Newton solver (Inexact Newton Method)

    template<class LAD>
    class ForcingStrategy
    {
      public:
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;
        typedef typename LAD::DataType DataType;

        ForcingStrategy ( );
        virtual ~ForcingStrategy ( );

        /// Get the current forcing term in stored vector

        DataType GetCurrentForcingTerm ( ) const
        {
            return forcing_terms_.back ( );
        }

        void SetForcingTerm ( DataType forcing )
        {
            forcing_terms_.push_back ( forcing );
        }

        void SetResidualNorm ( DataType norm )
        {
            residuals_.push_back ( norm );
        }

        virtual void ComputeForcingTerm ( DataType new_residual, DataType lin_solve_accuracy ) = 0;

        /// @return name of strategy

        std::string name ( ) const
        {
            return this->name_;
        }

      protected:
        std::string name_;
        std::vector<DataType> forcing_terms_;
        std::vector<DataType> residuals_;

    };

    /// standard constructor

    template<class DataType>
    ForcingStrategy<DataType>::ForcingStrategy ( )
    {
    }

    /// destructor

    template<class DataType>
    ForcingStrategy<DataType>::~ForcingStrategy ( )
    {
    }

} // namespace hiflow

#endif
