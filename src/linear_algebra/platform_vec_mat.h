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

/// @author Dimitar Lukarski, Nico Trost, Martin Wlotzka

#ifndef __PLATFORM_VEC_MAT_H
#    define __PLATFORM_VEC_MAT_H

#    include "linear_algebra/la_descriptor.h"

namespace hiflow
{
    namespace la
    {

        class LaCouplings;

        template<typename MatrixType>
        void init_platform_mat ( const PLATFORM platform, const IMPLEMENTATION matrix_impl,
                                 const enum MATRIX_FORMAT matrix_format, enum MATRIX_FREE_PRECOND matrix_precond,
                                 const MPI_Comm& comm_hf, const LaCouplings& cp,
                                 MatrixType *mat, MatrixType **dev_mat );

        template<typename MatrixType>
        void init_platform_mat ( const PLATFORM platform, const IMPLEMENTATION matrix_impl,
                                 const enum MATRIX_FORMAT matrix_format, enum MATRIX_FREE_PRECOND matrix_precond,
                                 const MPI_Comm& comm_hf, const LaCouplings& cp,
                                 MatrixType *mat, MatrixType **dev_mat, const SYSTEM& my_system );

        template<typename VectorType>
        void init_platform_vec ( const PLATFORM platform, const IMPLEMENTATION vector_impl,
                                 const MPI_Comm& comm_hf, const LaCouplings& cp,
                                 VectorType *vec, VectorType **dev_vec );

        template<typename VectorType>
        void init_platform_vec ( const PLATFORM platform, const IMPLEMENTATION vector_impl,
                                 const MPI_Comm& comm_hf, const LaCouplings& cp,
                                 VectorType *vec, VectorType **dev_vec,
                                 const SYSTEM& my_system );

    } // namespace la
} // namespace hiflow

#endif
