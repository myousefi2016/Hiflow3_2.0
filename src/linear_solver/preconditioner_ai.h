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

/// @author Dimitar Lukarski, Nico Trost, Niels Wegh

#ifndef HIFLOW_LINEARSOLVER_PRECONDITIONER_AI_H_
#    define HIFLOW_LINEARSOLVER_PRECONDITIONER_AI_H_

#    include <iostream>

#    include "preconditioner.h"
#    include "linear_algebra/lmp/lpreconditioner_ai.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        class PreconditionerApproximateInverse : public Preconditioner<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            PreconditionerApproximateInverse ( );
            virtual ~PreconditionerApproximateInverse ( );

            /// Setup the local operator for the local preconditioner
            virtual void SetupOperator ( OperatorType& op );
            virtual void SetupVector ( const VectorType& vec );

            /// Sets up paramaters
            virtual void Init_FSAI ( const int power );

            /// Build the preconditioner
            virtual void Build ( void );

            /// Applies the preconditioner on the diagonal block.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if preconditioning succeeded
            virtual LinearSolverState ApplyPreconditioner ( const VectorType& b,
                                                            VectorType* x );

            /// Erase possibly allocated data.
            virtual void Clear ( );

            /// Print possibly info
            virtual void Print ( std::ostream &out = std::cout ) const;

          protected:

            lPreconditioner_ApproximateInverse_FSAI<typename LAD::DataType> *local_ai_;

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PRECONDITIONER_AI_H_
