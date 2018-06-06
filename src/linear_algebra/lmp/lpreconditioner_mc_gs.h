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

/// @author Dimitar Lukarski

#ifndef __LPRECONDITIONER_MC_GS_H
#    define __LPRECONDITIONER_MC_GS_H

#    include <iostream>
#    include <stdlib.h>

#    include "lvector.h"
#    include "lmatrix.h"
#    include "lpreconditioner_mc.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Multi-Coloring  Gauss-Seidel Preconditioners for LU solvers
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_MultiColoring_GaussSeidel : public hiflow::la::lPreconditioner_MultiColoring<ValueType>
        {
          public:

            lPreconditioner_MultiColoring_GaussSeidel ( );
            virtual ~lPreconditioner_MultiColoring_GaussSeidel ( );

            virtual void Init ( void );

            virtual void build_aux_matrix ( void );
            virtual void allocate_LU ( void );
            virtual void factorize ( void );
            virtual void scale_D ( void );

            virtual void backwardstep ( void );
            virtual void diagonalstep ( void );
            virtual void forwardstep ( void );

        };

    } // namespace la
} // namespace hiflow

#endif
