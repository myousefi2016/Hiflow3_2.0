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

#ifndef HIFLOW_POLYNOMIALCHAOS_PCRESIDUALASSEMBLER_H_
#    define HIFLOW_POLYNOMIALCHAOS_PCRESIDUALASSEMBLER_H_

/// \file pc_residual_assembler.h
/// This abstract class is only needed for matrix free computations. It provides
/// an interface to assembly routines, which replace the matrix vector product.
/// An inherited instance of this class must by provided by the user and defined
/// in the corresponding application.
/// \author Michael Schick

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class LAD>
        class ResidualAssembler
        {
          public:
            typedef typename LAD::VectorType VectorType;

            /// Default constructor

            ResidualAssembler ( )
            {
            }
            /// Default destructor

            virtual ~ResidualAssembler ( )
            {
            }

            /// Abstract inferface, which should provide a functionality for assembling the residual, i.e.,
            /// res = b - Ax, where A is the Galerkin matrix, x the solution vector, and b the Galerkin rhs
            virtual void assemble_residual ( VectorType const& b, VectorType* x, VectorType& res ) const = 0;

            /// If a nonlinear PDE is solved, e.g. by Newton, then this function gives access to the
            /// linearization point and stores its values in the functions argument.

            virtual void get_linearization_point ( VectorType* u ) const
            {
            }
            /// If a nonlinear PDE is solved, e.g. by Newton, then this function gives access to the
            /// linearization point and sets its values provided by this functions argument.

            virtual void set_linearization_point ( VectorType const& u )
            {
            }
            /// Check if a nonlinear PDE is solved, default value is false, i.e., a linear PDE is assumed

            virtual bool nonlinear_problem ( ) const
            {
                return false;
            }
        };

    }
}

#endif
