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

#ifndef __INIT_VEC_MAT_H
#    define __INIT_VEC_MAT_H

#    include <string.h>

#    include "la_global.h"
#    include "lvector.h"
#    include "lmatrix.h"

/// @brief Matrix and Vector Initialization function
/// @author Dimitar Lukarski
///
/// Initialize vectors and matrices on a specific platform and implementation.

namespace hiflow
{
    namespace la
    {

        /// Create a (local) vector on a specific platform and implementation
        /// @return lVector<ValueType>
        /// @param size - size of the vector
        /// @param name - name of the vector
        /// @param platform - the platform
        /// @param implemenation - the implementation of the routines
        template <typename ValueType>
        lVector<ValueType> *init_vector ( const int size,
                                          const std::string name,
                                          const enum PLATFORM &platform,
                                          const enum IMPLEMENTATION &implementation );

        template <typename ValueType>
        lVector<ValueType> *init_vector ( const int size,
                                          const std::string name,
                                          const enum PLATFORM &platform,
                                          const enum IMPLEMENTATION &implementation,
                                          const struct SYSTEM &my_system );

        /// Create a (local) matrix on a specific platform, implementation and matrix format
        /// @return lMatrix<ValueType>
        /// @param init_nnz - number of non-zeros
        /// @param init_num_row - row size
        /// @param init_num_col - colum size
        /// @param name - name of the matrix
        /// @param platform - the platform
        /// @param implemenation - the implementation of the routines
        /// @param matrix_format - the format of the matrix
        template <typename ValueType>
        lMatrix<ValueType> *init_matrix ( const int init_nnz,
                                          const int init_num_row,
                                          const int init_num_col,
                                          const std::string init_name,
                                          const enum PLATFORM &platform,
                                          const enum IMPLEMENTATION &implementation,
                                          const enum MATRIX_FORMAT &matrix_format );

        template <typename ValueType>
        lMatrix<ValueType> *init_matrix ( const int init_nnz,
                                          const int init_num_row,
                                          const int init_num_col,
                                          const std::string init_name,
                                          const enum PLATFORM &platform,
                                          const enum IMPLEMENTATION &implementation,
                                          const enum MATRIX_FORMAT &matrix_format,
                                          const struct SYSTEM &my_system );

    } // namespace la
} // namespace hiflow

#endif
