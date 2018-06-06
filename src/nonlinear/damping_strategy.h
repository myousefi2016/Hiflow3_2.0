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

#ifndef HIFLOW_NONLINEAR_DAMPING_STRATEGY_H_
#    define HIFLOW_NONLINEAR_DAMPING_STRATEGY_H_

#    include <string>

namespace hiflow
{

    /// Enumerator @em DampingState as return value for the update function

    enum DampingState
    {
        kDampingSuccess = 0,
        kDampingExceeded,
        kDampingError
    };

    /// @brief Base class for damping strategies
    /// @author Tobias Hahn
    ///
    /// A damping strategy implements the update routine within a
    /// Newton solver

    template<class LAD>
    class DampingStrategy
    {
      public:
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;
        typedef typename LAD::DataType DataType;

        DampingStrategy ( );
        virtual ~DampingStrategy ( );

        virtual void Init ( )
        {
        }

        // With RHS
        virtual DampingState Update ( const VectorType& cor, const VectorType& rhs, VectorType* res, VectorType* sol, void* myNewton ) = 0;

        virtual DataType GetResidual ( )
        {
            return this->residual_;
        }

        /// @return name of strategy

        std::string name ( ) const
        {
            return this->name_;
        }

      protected:
        std::string name_;
        DataType residual_;

    };

    /// standard constructor

    template<class DataType>
    DampingStrategy<DataType>::DampingStrategy ( )
    {
    }

    /// destructor

    template<class DataType>
    DampingStrategy<DataType>::~DampingStrategy ( )
    {
    }

} // namespace hiflow

#endif  // HIFLOW_NONLINEARSOLVER_NONLINEAR_SOLVER_H_
