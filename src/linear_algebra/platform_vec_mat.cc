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

#include "platform_vec_mat.h"

namespace hiflow
{
    namespace la
    {

        template<typename MatrixType>
        void init_platform_mat ( const PLATFORM platform, const IMPLEMENTATION matrix_impl,
                                 const enum MATRIX_FORMAT matrix_format, enum MATRIX_FREE_PRECOND matrix_precond,
                                 const MPI_Comm& comm, const LaCouplings& la_couplings,
                                 MatrixType *mat, MatrixType **dev_mat )
        {
            MATRIX_FORMAT default_format = CSR;

            // Create new Platform Mat/Vec
            if ( platform != CPU )
            {
                // cpu
                // init empty matrix
                mat->Init ( comm, la_couplings, CPU, NAIVE, default_format );

                // matrix
                ( *dev_mat ) = new MatrixType;
                ( *dev_mat )->Init ( comm, la_couplings, platform, matrix_impl, matrix_format );

            }
            else
            {
                // no accelerator dev (CPU only)
                // init empty matrix
                if ( matrix_format == default_format )
                {
                    mat->Init ( comm, la_couplings, platform, matrix_impl, matrix_format );
                    ( *dev_mat ) = mat;
                }
                else
                {
                    mat->Init ( comm, la_couplings, CPU, NAIVE, default_format );
                    ( *dev_mat ) = new MatrixType;
                    ( *dev_mat )->Init ( comm, la_couplings, platform, matrix_impl, matrix_format );
                }

            }

        }

        template<typename MatrixType>
        void init_platform_mat ( const PLATFORM platform, const IMPLEMENTATION matrix_impl,
                                 const enum MATRIX_FORMAT matrix_format, enum MATRIX_FREE_PRECOND matrix_precond,
                                 const MPI_Comm& comm, const LaCouplings& la_couplings,
                                 MatrixType *mat, MatrixType **dev_mat, const SYSTEM& my_system )
        {
            MATRIX_FORMAT default_format = CSR;

            // Create new Platform Mat/Vec
            if ( platform != CPU )
            {
                // cpu
                // init empty matrix
                mat->Init ( comm, la_couplings, CPU, NAIVE, default_format );
                // matrix
                ( *dev_mat ) = new MatrixType;
                if ( platform == OPENCL )
                    ( *dev_mat )->Init ( comm, la_couplings, platform, matrix_impl, matrix_format, my_system );
                else
                    (*dev_mat )->Init ( comm, la_couplings, platform, matrix_impl, matrix_format );
            }
            else
            {
                // no accelerator dev (CPU only)
                // init empty matrix
                if ( matrix_format == default_format )
                {
                    mat->Init ( comm, la_couplings, platform, matrix_impl, matrix_format );
                    ( *dev_mat ) = mat;
                }
                else
                {
                    mat->Init ( comm, la_couplings, CPU, NAIVE, default_format );
                    ( *dev_mat ) = new MatrixType;
                    ( *dev_mat )->Init ( comm, la_couplings, platform, matrix_impl, matrix_format );
                }

            }

        }

        template<typename VectorType>
        void init_platform_vec ( const PLATFORM platform, const IMPLEMENTATION vector_impl,
                                 const MPI_Comm& comm, const LaCouplings& la_couplings,
                                 VectorType *vec, VectorType **dev_vec )
        {

            // Create new Platform Mat/Vec
            if ( platform != CPU )
            {
                // cpu
                // vectors
                vec->Init ( comm, la_couplings, CPU, NAIVE );

                ( *dev_vec ) = new VectorType;

                ( *dev_vec )->Init ( comm, la_couplings, platform, vector_impl );

            }
            else
            {
                // no accelerator dev (CPU only)

                // vectors
                vec->Init ( comm, la_couplings, platform, vector_impl );
                ( *dev_vec ) = vec;
            }

        }

        template<typename VectorType>
        void init_platform_vec ( const PLATFORM platform, const IMPLEMENTATION vector_impl,
                                 const MPI_Comm& comm, const LaCouplings& la_couplings,
                                 VectorType *vec, VectorType **dev_vec,
                                 const SYSTEM& my_system )
        {

            // Create new Platform Mat/Vec
            if ( platform != CPU )
            {
                // cpu
                // vectors
                vec->Init ( comm, la_couplings, CPU, NAIVE );

                ( *dev_vec ) = new VectorType;
                if ( platform == OPENCL )
                    ( *dev_vec )->Init ( comm, la_couplings, platform, vector_impl, my_system );
                else
                    (*dev_vec )->Init ( comm, la_couplings, platform, vector_impl );
            }
            else
            {
                // no accelerator dev (CPU only)

                // vectors
                vec->Init ( comm, la_couplings, platform, vector_impl );
                ( *dev_vec ) = vec;
            }

        }

        template void init_platform_mat<LADescriptorCoupledD::MatrixType>( const PLATFORM platform, const IMPLEMENTATION matrix_impl,
                const enum MATRIX_FORMAT matrix_format, enum MATRIX_FREE_PRECOND matrix_precond,
                const MPI_Comm& comm, const LaCouplings& la_couplings,
                LADescriptorCoupledD::MatrixType *mat, LADescriptorCoupledD::MatrixType **dev_mat );

        template void init_platform_mat<LADescriptorCoupledD::MatrixType>( const PLATFORM platform, const IMPLEMENTATION matrix_impl,
                const enum MATRIX_FORMAT matrix_format, enum MATRIX_FREE_PRECOND matrix_precond,
                const MPI_Comm& comm, const LaCouplings& la_couplings,
                LADescriptorCoupledD::MatrixType *mat, LADescriptorCoupledD::MatrixType **dev_mat,
                const SYSTEM& my_system );

        template void init_platform_vec<LADescriptorCoupledD::VectorType>( const PLATFORM platform, const IMPLEMENTATION vector_impl,
                const MPI_Comm& comm, const LaCouplings& la_couplings,
                LADescriptorCoupledD::VectorType *vec, LADescriptorCoupledD::VectorType **dev_vec );

        template void init_platform_vec<LADescriptorCoupledD::VectorType>( const PLATFORM platform, const IMPLEMENTATION vector_impl,
                const MPI_Comm& comm, const LaCouplings& la_couplings,
                LADescriptorCoupledD::VectorType *vec, LADescriptorCoupledD::VectorType **dev_vec,
                const SYSTEM& my_system );

    } // namespace la
} // namespace hiflow
