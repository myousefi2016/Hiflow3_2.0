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

#ifndef HIFLOW_NONLINEAR_DAMPING_ARMIJO_H_
#    define HIFLOW_NONLINEAR_DAMPING_ARMIJO_H_

#    include <string>
#    include <vector>
#    include "linear_algebra/la_descriptor.h"
#    include "damping_strategy.h"

namespace hiflow
{

    /// Enumerator @em ArmijoParam User changable

    enum ArmijoParam
    {
        ArmijoInitial,
        ArmijoMinimal,
        ArmijoMaxLoop,
        ArmijoDecrease,
        ArmijoSuffDec
    };

    /// @brief Base class for damping strategies
    /// @author Tobias Hahn
    ///
    /// Solves for x in F(x)=y with nonlinear F

    template<class LAD>
    class ArmijoDamping : public DampingStrategy<LAD>
    {
      public:
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;
        typedef typename LAD::DataType DataType;

        ArmijoDamping ( );
        ArmijoDamping ( DataType init, DataType mini, DataType dec, DataType suffdec, int maxloop );
        ~ArmijoDamping ( );

        DampingState Init ( ArmijoParam param );
        DampingState Init ( ArmijoParam param, int data );
        DampingState Init ( ArmijoParam param, DataType data );

        DampingState Update ( const VectorType& cor, const VectorType& rhs, VectorType* res, VectorType* sol, void* myNewton );

      private:
        std::string name_;
        DataType Initial_;
        DataType Minimal_;
        DataType Decrease_;
        DataType SuffDec_;
        int MaxLoop_;
        std::vector<DataType>* ForcingTerm_;
        void Average ( VectorType* x );
    };

} // namespace hiflow

#endif  // HIFLOW_NONLINEAR_DAMPING_ARMIJO_H_
