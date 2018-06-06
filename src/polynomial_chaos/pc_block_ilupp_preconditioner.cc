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

/// \author Michael Schick

#include "polynomial_chaos/pc_block_ilupp_preconditioner.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class LAD>
        PCBlockIluppPreconditioner<LAD>::PCBlockIluppPreconditioner ( ) : la::Preconditioner<LAD>( )
        {
        }

        template<class LAD>
        PCBlockIluppPreconditioner<LAD>::~PCBlockIluppPreconditioner ( )
        {
        }

        template<class LAD>
        la::LinearSolverState PCBlockIluppPreconditioner<LAD>::ApplyPreconditioner ( const VectorType& b,
                                                                                     VectorType* x )
        {
            for ( int mode = 0; mode < b.NModes ( ); ++mode )
            {
                this->ilupp_.ApplyPreconditioner ( *b.Mode ( mode ), x->Mode ( mode ) );
            }
            return la::kSolverSuccess;
        }

        /// template instantiation
        template class PCBlockIluppPreconditioner<la::LADescriptorPolynomialChaosD>;

    } // namespace polynomialchaos
} // namespace hiflow
