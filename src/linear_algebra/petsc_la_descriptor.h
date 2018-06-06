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
/// @date 2015-12-02

#ifndef HIFLOW_LINEARALGEBRA_PETSC_LA_DESCRIPTOR_H_
#    define HIFLOW_LINEARALGEBRA_PETSC_LA_DESCRIPTOR_H_

#    include "linear_algebra/petsc_matrix.h"
#    include "linear_algebra/petsc_vector.h"
#    include <complex>

namespace hiflow
{
    namespace la
    {

        /// @brief Linear algebra descriptors defining PETSc matrix/vector pair. TODO implement complex version

        class LADescriptorPETSc
        {
          public:
            typedef PETScMatrix<double> MatrixType;
            typedef PETScVector<double> VectorType;
            typedef double DataType;
        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_PETSC_LA_DESCRIPTOR_H_
