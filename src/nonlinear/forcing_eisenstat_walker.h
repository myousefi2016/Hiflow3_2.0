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

#ifndef HIFLOW_NONLINEAR_EISENSTAT_WALKER_H_
#    define HIFLOW_NONLINEAR_EISENSTAT_WALKER_H_

#    include <string>
#    include <vector>
#    include "linear_algebra/la_descriptor.h"
#    include "forcing_strategy.h"

namespace hiflow
{

    /// @brief Eisenstat and Walker computation of forcing terms
    /// @author Michael Schick

    template<class LAD>
    class EWForcing : public ForcingStrategy<LAD>
    {
      public:
        typedef typename LAD::MatrixType MatrixType;
        typedef typename LAD::VectorType VectorType;
        typedef typename LAD::DataType DataType;

        EWForcing ( );
        EWForcing ( int choice );
        EWForcing ( double initial, double max, int choice );
        ~EWForcing ( );

        // constructor
        EWForcing ( double initial, double max, int choice, double gamma, double alpha );

        void ComputeForcingTerm ( DataType new_residual, DataType lin_solve_accuracy );

      private:

        void ComputeForcingChoice1 ( DataType new_residual, DataType lin_solve_accuracy );
        void ComputeForcingChoice2 ( DataType new_residual );

        using ForcingStrategy<LAD>::name_;
        using ForcingStrategy<LAD>::forcing_terms_;
        using ForcingStrategy<LAD>::residuals_;

        DataType Initial_;
        DataType Maximal_;
        int Choice_;

        DataType alpha_, gamma_;
    };

} // namespace hiflow

#endif
