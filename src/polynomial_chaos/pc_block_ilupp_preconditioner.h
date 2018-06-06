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

#ifndef HIFLOW_POLYNOMIALCHAOS_PCBLOCKILUPP_PRECONDITIONER_H_
#    define HIFLOW_POLYNOMIALCHAOS_PCBLOCKILUPP_PRECONDITIONER_H_

/// \file preconditioner_ilupp_galerkin.h
/// \brief Blockwise ILU++ Galerkin solver for Polynomial Chaos.
/// Applies ILU preconditioning (by external library ILU++) to
/// each Polynomial Chaos mode (blockwise approach)
/// \author Michael Schick

#    include <vector>

#    include "linear_solver/preconditioner_ilupp.h"
#    include "iluplusplus_interface.h"
#    include "linear_algebra/la_descriptor.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class LAD>
        class PCBlockIluppPreconditioner : public la::Preconditioner<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// Default constructor
            PCBlockIluppPreconditioner ( );
            /// Destructor
            virtual ~PCBlockIluppPreconditioner ( );

            /// Inits parameters for ILU++ preconditioner.
            /// \param prepro_type type of preprocessing
            /// \param precond_no number of preconditioner
            /// \param max_levels maximum number of multilevels
            /// \param mem_factor see ILU++ manual
            /// \param threshold see ILU++ manual
            /// \param min_pivot see ILU++ manual

            virtual void InitParameter ( int prepro_type, int precond_no,
                                         int max_levels, double mem_factor,
                                         double threshold, double min_pivot )
            {
                ilupp_.InitParameter ( prepro_type, precond_no, max_levels, mem_factor, threshold, min_pivot );
            }

            virtual void SetupOperator ( OperatorType& op )
            {
                ilupp_.SetupOperator ( *op.Mode ( 0 ) );
            }

            /// Applies the ILU++ preconditioner (blockwise).
            virtual la::LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );

          private:
            /// Associated linear operator
            la::PreconditionerIlupp<la::LADescriptorCoupledD> ilupp_;
        };

    } // namespace polynomialchaos
} // namespace hiflow

#endif  // HIFLOW_POLYNOMIALCHAOS_PCBLOCKILUPP_PRECONDITIONER_H_
