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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-07-23

#include <assert.h>
#include <sys/time.h>
#include <fstream>
#include "../la_global.h"
#include "../lmp_log.h"
#include "../lvector.h"
#include "../lmatrix.h"
#include "../lpreconditioner.h"
#include "../init_vec_mat.h"
#include "../platform_management.h"
#include "solvers/cg.h"

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        int cg ( lVector<ValueType> *x, const lVector<ValueType> *b,
                 const lMatrix<ValueType> *matrix, const ValueType rel_eps,
                 const ValueType abs_eps, int &max_iter, const int print,
                 lPreconditioner<ValueType> *lPrecond )
        {
            return CGFunctor<ValueType, false>( )( x, b, matrix, rel_eps, abs_eps, max_iter,
                    print, lPrecond );
        }

        /// Instantiate template for double
        template int cg ( lVector<double> *x, const lVector<double> *b,
                          const lMatrix<double> *matrix, const double rel_eps,
                          const double abs_eps, int &max_iter, const int print,
                          lPreconditioner<double> *lPrecond = NULL );

        /// Instantiate template for float
        template int cg ( lVector<float> *x, const lVector<float> *b,
                          const lMatrix<float> *matrix, const float rel_eps,
                          const float abs_eps, int &max_iter, const int print,
                          lPreconditioner<float> *lPrecond = NULL );

    } // namespace hiflow
} // namespace la
