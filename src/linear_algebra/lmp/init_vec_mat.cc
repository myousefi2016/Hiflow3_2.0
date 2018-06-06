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

#include "config.h"

#include "init_vec_mat.h"

#include "lvector_cpu.h"

#include "lmatrix_dense_cpu.h"
#include "lmatrix_csr_cpu.h"
#include "lmatrix_coo_cpu.h"

#include "cuda/lmatrix_csr_gpu.h"
#include "cuda/lmatrix_coo_gpu.h"
#include "cuda/GPUcublas2_lVector.h"
#include "cuda/GPUcublas2_CSR_lMatrix.h"
#include "cuda/lvector_gpu.h"

#include "opencl/lvector_opencl.h"
#include "opencl/lmatrix_csr_opencl.h"

#include "lmp_log.h"

#include <iostream>
#include <stdlib.h>

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        lVector<ValueType> *init_vector ( const int size,
                                          const std::string name,
                                          const enum PLATFORM &platform,
                                          const enum IMPLEMENTATION &implementation )
        {
            if ( platform == CPU )
            {

                switch ( implementation )
                {

                    case NAIVE:
                        return new CPUsimple_lVector<ValueType>( size, name );
                        break;

                    case OPENMP:
                        return new CPUopenmp_lVector<ValueType>( size, name );
                        break;

                    case MKL:
                        return new CPUmkl_lVector<ValueType>( size, name );
                        break;

                    case BLAS:
                        return new CPUcblas_lVector<ValueType>( size, name );
                        break;

                    default:
                        ;
                }

            }

#ifdef WITH_CUDA
            if ( platform == GPU )
            {

                switch ( implementation )
                {

                    case BLAS:
                        return new GPUblas_lVector<ValueType>( size, name );
                        break;

                    case CUBLAS2:
                        return new GPUcublas2_lVector<ValueType>( size, name );
                        break;

                    default:
                        ;
                }

            }
#endif

            LOG_ERROR ( "init_vector() incompatibility PLATFORM/IMPLEMENTATION" );
            LOG_ERROR ( " Platform ID=" << platform << " Implementation ID=" << implementation );
            exit ( -1 );
            return NULL;
        }

        template <typename ValueType>
        lVector<ValueType> *init_vector ( const int size,
                                          const std::string name,
                                          const enum PLATFORM &platform,
                                          const enum IMPLEMENTATION &implementation,
                                          const struct SYSTEM &my_system )
        {
            if ( platform == CPU )
            {

                switch ( implementation )
                {

                    case NAIVE:
                        return new CPUsimple_lVector<ValueType>( size, name );
                        break;

                    case OPENMP:
                        return new CPUopenmp_lVector<ValueType>( size, name );
                        break;

                    case MKL:
                        return new CPUmkl_lVector<ValueType>( size, name );
                        break;

                    case BLAS:
                        return new CPUcblas_lVector<ValueType>( size, name );
                        break;

                    default:
                        ;
                }

            }

#ifdef WITH_CUDA
            if ( platform == GPU )
            {

                switch ( implementation )
                {

                    case BLAS:
                        return new GPUblas_lVector<ValueType>( size, name );
                        break;

                    default:
                        ;
                }

            }
#endif

#ifdef WITH_OPENCL
            if ( platform == OPENCL )
            {
                lVector<ValueType> *new_vector = new OPENCL_lVector<ValueType>( my_system.my_manager );
                new_vector->Init ( size, name );
                return new_vector;
            }
#endif

            LOG_ERROR ( "init_vector() incompatibility PLATFORM/IMPLEMENTATION" );
            LOG_ERROR ( " Platform ID=" << platform << " Implementation ID=" << implementation );
            exit ( -1 );
            return NULL;
        }

        template <typename ValueType>
        lMatrix<ValueType> *init_matrix ( const int init_nnz,
                                          const int init_num_row,
                                          const int init_num_col,
                                          const std::string init_name,
                                          const enum PLATFORM &platform,
                                          const enum IMPLEMENTATION &implementation,
                                          const enum MATRIX_FORMAT &matrix_format )
        {
            // DENSE
            if ( matrix_format == DENSE )
            {
                if ( platform == CPU )
                {
                    switch ( implementation )
                    {
                        case NAIVE:
                            return new CPUsimple_DENSE_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;
                    }
                }
                // CSR
            }
            else if ( matrix_format == CSR )
            {
                if ( platform == CPU )
                {

                    switch ( implementation )
                    {

                        case NAIVE:
                            return new CPUsimple_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case OPENMP:
                            return new CPUopenmp_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case MKL:
                            return new CPUmkl_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        default:
                            ;
                    }
                }

#ifdef WITH_CUDA
                if ( platform == GPU )
                {

                    switch ( implementation )
                    {

                        case SCALAR:
                            return new GPUscalar_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case SCALAR_TEX:
                            return new GPUscalartex_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case CUBLAS2:
                            return new GPUcublas2_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        default:
                            ;
                    }
                }
#endif

                // COO
            }
            else if ( matrix_format == COO )
            {
                if ( platform == CPU )
                {

                    switch ( implementation )
                    {

                        case NAIVE:
                            return new CPUsimple_COO_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case OPENMP:
                            return new CPUopenmp_COO_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case MKL:
                            return new CPUmkl_COO_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        default:
                            ;
                    }
                }

#ifdef WITH_CUDA
                if ( platform == GPU )
                {
                    // no COO on GPU
                }
#endif
            }

            LOG_ERROR ( "init_matrix() incompatibility PLATFORM/IMPLEMENTATION/FORMAT" );
            LOG_ERROR ( "Platform ID=" << platform << " Implementation ID=" << implementation
                        << " Matrix Format ID=" << matrix_format );
            exit ( -1 );
            return NULL;

        }

        template <typename ValueType>
        lMatrix<ValueType> *init_matrix ( const int init_nnz,
                                          const int init_num_row,
                                          const int init_num_col,
                                          const std::string init_name,
                                          const enum PLATFORM &platform,
                                          const enum IMPLEMENTATION &implementation,
                                          const enum MATRIX_FORMAT &matrix_format,
                                          const struct SYSTEM &my_system )
        {
            // DENSE
            if ( matrix_format == DENSE )
            {
                if ( platform == CPU )
                {
                    switch ( implementation )
                    {
                        case NAIVE:
                            return new CPUsimple_DENSE_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;
                    }
                }
                // CSR
            }
            else if ( matrix_format == CSR )
            {
                if ( platform == CPU )
                {

                    switch ( implementation )
                    {

                        case NAIVE:
                            return new CPUsimple_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case OPENMP:
                            return new CPUopenmp_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case MKL:
                            return new CPUmkl_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        default:
                            ;
                    }
                }

#ifdef WITH_CUDA
                if ( platform == GPU )
                {

                    switch ( implementation )
                    {

                        case SCALAR:
                            return new GPUscalar_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case SCALAR_TEX:
                            return new GPUscalartex_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case CUBLAS2:
                            return new GPUcublas2_CSR_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        default:
                            ;
                    }
                }
#endif
            }

            // COO
            if ( matrix_format == COO )
            {
                if ( platform == CPU )
                {

                    switch ( implementation )
                    {

                        case NAIVE:
                            return new CPUsimple_COO_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case OPENMP:
                            return new CPUopenmp_COO_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        case MKL:
                            return new CPUmkl_COO_lMatrix<ValueType>( init_nnz, init_num_row, init_num_col, init_name );
                            break;

                        default:
                            ;
                    }

                }

#ifdef WITH_CUDA
                if ( platform == GPU )
                {
                    // no COO on GPU
                }
#endif
            }

#ifdef WITH_OPENCL
            if ( matrix_format == CSR && platform == OPENCL )
            {
                lMatrix<ValueType> *new_matrix = new OPENCL_CSR_lMatrix<ValueType>( my_system.my_manager );
                new_matrix->Init ( init_nnz, init_num_row, init_num_col, init_name );
                return new_matrix;
            }
#endif

            LOG_ERROR ( "init_matrix() incompatibility PLATFORM/IMPLEMENTATION/FORMAT" );
            LOG_ERROR ( "Platform ID=" << platform << " Implementation ID=" << implementation
                        << " Matrix Format ID=" << matrix_format );
            exit ( -1 );
            return NULL;

        }

        template lVector<double> *init_vector ( const int size,
                                                const std::string name,
                                                const enum PLATFORM &platform,
                                                const enum IMPLEMENTATION &implementation );

        template lVector<double> *init_vector ( const int size,
                                                const std::string name,
                                                const enum PLATFORM &platform,
                                                const enum IMPLEMENTATION &implementation,
                                                const struct SYSTEM &my_system );

        template lVector<float> *init_vector ( const int size,
                                               const std::string name,
                                               const enum PLATFORM &platform,
                                               const enum IMPLEMENTATION &implementation );

        template lVector<float> *init_vector ( const int size,
                                               const std::string name,
                                               const enum PLATFORM &platform,
                                               const enum IMPLEMENTATION &implementation,
                                               const struct SYSTEM &my_system );

        template lMatrix<double> *init_matrix ( const int init_nnz,
                                                const int init_num_row,
                                                const int init_num_col,
                                                const std::string init_name,
                                                const enum PLATFORM &platform,
                                                const enum IMPLEMENTATION &implementation,
                                                const enum MATRIX_FORMAT &matrix_format );

        template lMatrix<double> *init_matrix ( const int init_nnz,
                                                const int init_num_row,
                                                const int init_num_col,
                                                const std::string init_name,
                                                const enum PLATFORM &platform,
                                                const enum IMPLEMENTATION &implementation,
                                                const enum MATRIX_FORMAT &matrix_format,
                                                const struct SYSTEM &my_system );

        template lMatrix<float> *init_matrix ( const int init_nnz,
                                               const int init_num_row,
                                               const int init_num_col,
                                               const std::string init_name,
                                               const enum PLATFORM &platform,
                                               const enum IMPLEMENTATION &implementation,
                                               const enum MATRIX_FORMAT &matrix_format );

        template lMatrix<float> *init_matrix ( const int init_nnz,
                                               const int init_num_row,
                                               const int init_num_col,
                                               const std::string init_name,
                                               const enum PLATFORM &platform,
                                               const enum IMPLEMENTATION &implementation,
                                               const enum MATRIX_FORMAT &matrix_format,
                                               const struct SYSTEM &my_system );

    } // namespace la
} // namespace hiflow
