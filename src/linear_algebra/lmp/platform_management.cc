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

#include <iostream>
#include <stdlib.h>
#include <typeinfo>

#include "config.h"

#include "la_global.h"
#include "platform_management.h"

#include "cuda/platform_management_gpu.h"
#include "opencl/platform_management_opencl.h"

#include "lmp_log.h"

namespace hiflow
{
    namespace la
    {

        // select platform

        void select_platform ( struct SYSTEM &my_system,
                               enum PLATFORM &v_platform, enum IMPLEMENTATION &v_implementation,
                               enum PLATFORM &m_platform, enum IMPLEMENTATION &m_implementation,
                               enum MATRIX_FORMAT &m_format, enum MATRIX_FREE_PRECOND &m_precond )
        {

            int input;

            std::cout << "Hardware Platform: CPU=0, GPU=1, OPENCL=2 >";
            std::cin >> input;
            my_system.Platform = PLATFORM ( input );

            // to ensure the compatibility
            v_platform = my_system.Platform;
            m_platform = my_system.Platform;

            std::cout << "lVector Implementaton: (CPU)NAIVE=0, (CPU/CBLAS,GPU/CUBLAS)BLAS=1, (CPU)OpenMP=2, (CPU)MKL=3, (CPU/GPU)OPENCL=6 >";
            std::cin >> input;
            v_implementation = IMPLEMENTATION ( input );

            std::cout << "lMatrix Implementaton:"
                    << " (CPU)NAIVE=0, (none)BLAS=1 (CPU)OpenMP=2"
                    << " (CPU)MKL=3, (GPU)SCALAR=4, (GPU)SCALAR+TEX_CACHE=5, (CPU/GPU)OPENCL=6 >";
            std::cin >> input;
            m_implementation = IMPLEMENTATION ( input );

            my_system.GPU_CUBLAS = ( v_implementation == BLAS );

            std::cout << "lMatrix Format: (CPU,GPU)CSR=1, (CPU)COO=2 >";
            std::cin >> input;
            m_format = MATRIX_FORMAT ( input );

            //  std::cout << "lMatrix-free-Preconditioner (CPU-only, CSR, COO):"
            //            << " No Precond=0, Jacobi=1, Gauss Seidel=2, Symmetric Gauss Seidel=3, SOR=4, SSOR=5 >";
            //  std::cin >> input ;
            //  m_precond = MATRIX_FREE_PRECOND (input) ;
            m_precond = MATRIX_FREE_PRECOND ( 0 );

            /*
            #ifndef OPENCL
              std::cout << "Precision: FLOAT=0, DOUBLE=1 >";
              std::cin >> input ;
              my_system.is_double = PRECISION (input) ;
            #endif
             */

            // manual selecting ? -> no rank
            my_system.rank = -1;
        }

        template <typename ValueType>
        void select_platform ( struct SYSTEM &my_system,
                               enum PLATFORM &v_platform, enum IMPLEMENTATION &v_implementation,
                               enum PLATFORM &m_platform, enum IMPLEMENTATION &m_implementation,
                               enum MATRIX_FORMAT &m_format, enum MATRIX_FREE_PRECOND &m_precond )
        {

            int input;

            std::cout << "Hardware Platform: CPU=0, GPU=1, OPENCL=2 >";
            std::cin >> input;
            my_system.Platform = PLATFORM ( input );

            // to ensure the compatibility
            v_platform = my_system.Platform;
            m_platform = my_system.Platform;

            std::cout << "lVector Implementaton: (CPU)NAIVE=0, (CPU/CBLAS,GPU/CUBLAS)BLAS=1, (CPU)OpenMP=2, (CPU)MKL=3, (CPU/GPU)OPENCL=6 >";
            std::cin >> input;
            v_implementation = IMPLEMENTATION ( input );

            std::cout << "lMatrix Implementaton:"
                    << " (CPU)NAIVE=0, (none)BLAS=1 (CPU)OpenMP=2"
                    << " (CPU)MKL=3, (GPU)SCALAR=4, (GPU)SCALAR+TEX_CACHE=5, (CPU/GPU)OPENCL=6 >";
            std::cin >> input;
            m_implementation = IMPLEMENTATION ( input );

            my_system.GPU_CUBLAS = ( v_implementation == BLAS );

            std::cout << "lMatrix Format: (CPU,GPU)CSR=1, (CPU)COO=2 >";
            std::cin >> input;
            m_format = MATRIX_FORMAT ( input );

            //  std::cout << "lMatrix-free-Preconditioner (CPU-only, CSR, COO):"
            //            << " No Precond=0, Jacobi=1, Gauss Seidel=2, Symmetric Gauss Seidel=3, SOR=4, SSOR=5 >";
            //  std::cin >> input ;
            //  m_precond = MATRIX_FREE_PRECOND (input) ;
            m_precond = MATRIX_FREE_PRECOND ( 0 );

            if ( typeid (ValueType ) == typeid (float ) )
                my_system.is_double = false;
            else
                my_system.is_double = true;
        }

        // init platform CPU/GPU

        void init_platform ( struct SYSTEM &my_system )
        {
            //  my_system.rank = -1;

            if ( my_system.Platform == CPU )
            {
                // do nothing :)
                //    std::cout << "init_platform using CPU ...done" << std::endl ;
            }
#ifdef WITH_CUDA
            else if ( my_system.Platform == GPU )
            {
                if ( my_system.rank < 0 )
                {
                    // select the next available device
                    init_platform_gpu ( my_system.GPU_CUBLAS );
                }
                else
                {
                    // select rank specific device
                    int num_gpu_per_node = 1;
                    init_platform_gpu ( my_system.GPU_CUBLAS, my_system.rank % num_gpu_per_node );
                }
            }
#endif
#ifdef WITH_OPENCL
            else if ( my_system.Platform == OPENCL )
            {

                if ( my_system.rank < 0 )
                {
                    // ask what to select
                    init_platform_opencl ( my_system );
                }
                else
                {
                    // select rank specific device
                    int opencl_platform = 0;
                    int opencl_devices = 2; // devince per node
                    init_platform_opencl ( my_system, opencl_platform, my_system.rank % opencl_devices );
                }

            }
#endif
            else
            {
                LOG_ERROR ( "init_platform() unknown platform" );
                exit ( -1 );
            }
            my_system.Initialized = true;
        }

        // stop platform CPU/GPU

        void stop_platform ( struct SYSTEM &my_system )
        {
            if ( my_system.Platform != CPU )
            {
#ifdef WITH_CUDA
                if ( my_system.Platform == GPU )
                {
                    stop_platform_gpu ( my_system.GPU_CUBLAS );
                }
#endif

#ifdef WITH_OPENCL
                if ( my_system.Platform == OPENCL )
                {
                    stop_platform_opencl ( my_system );
                }
#endif
                LOG_ERROR ( "stop_platform() unknown platform" );
                exit ( -1 );
            }
            my_system.Initialized = false;
        }

        template void select_platform<float>( struct SYSTEM &my_system,
                enum PLATFORM &v_platform, enum IMPLEMENTATION &v_implementation,
                enum PLATFORM &m_platform, enum IMPLEMENTATION &m_implementation,
                enum MATRIX_FORMAT &m_format, enum MATRIX_FREE_PRECOND &m_precond );

        template void select_platform<double>( struct SYSTEM &my_system,
                enum PLATFORM &v_platform, enum IMPLEMENTATION &v_implementation,
                enum PLATFORM &m_platform, enum IMPLEMENTATION &m_implementation,
                enum MATRIX_FORMAT &m_format, enum MATRIX_FREE_PRECOND &m_precond );

    } // namespace la
} // namespace hiflow
