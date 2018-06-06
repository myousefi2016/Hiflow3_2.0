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

/// @author Dimitar Lukarski, Chandramowli Subramanian

#ifndef HIFLOW_LINEARSOLVER_PRECONDITIONER_B_JACOBI_H_
#    define HIFLOW_LINEARSOLVER_PRECONDITIONER_B_JACOBI_H_

#    include "preconditioner.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Base class for all block Jacobi preconditioners.

        template<class LAD>
        class PreconditionerBlockJacobi : public Preconditioner<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            PreconditionerBlockJacobi ( )
            : Preconditioner<LAD>( ), maxits_ ( 1000 )
            {
            }

            virtual ~PreconditionerBlockJacobi ( )
            {
            }

          protected:
            int maxits_;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PRECONDITIONER_B_JACOBI_H_
