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

/// @author Nico Trost, Benedikt Galler, Dimitar Lukarski

#ifndef __LMATRIX_CSR_OPENCL_H
#    define __LMATRIX_CSR_OPENCL_H

#    include "config.h"

#    include <iostream>
#    include <stdlib.h>

#    include "opencl/lvector_opencl.h"
#    include "../lmatrix_csr.h"

#    ifdef WITH_OPENCL

#        ifdef __APPLE__
#            include <cl.h>
#        else
#            include <CL/cl.h>
#        endif

#        include "opencl/mem_opencl.h"
#        include "opencl/opencl_global.h"

#    endif

/// @brief Provides the base matrix (CSR) class for OPENCL
/// @author Nico Trost, Benedikt Galler, Dimitar Lukarski
///
/// CPU_CSR_lMatrix maintains the matrix;
/// provides accessing functions and the
/// assignment operator (assign matrix from cpu and opencl)

template <typename ValueType>
class OPENCL_CSR_lMatrix : public CSR_lMatrix<ValueType>
{
  public:

    OPENCL_CSR_lMatrix ( );

#    ifdef WITH_OPENCL
    OPENCL_CSR_lMatrix ( opencl_manager *man );
#    endif

    /// The virtual destructor free the buffers of the matrix
    virtual ~OPENCL_CSR_lMatrix ( );

#    ifdef WITH_OPENCL
    virtual void set_opencl_manager ( opencl_manager *man );
    opencl_manager* get_opencl_manager ( );
#    endif

    virtual void CopyFrom ( const lMatrix<ValueType> &mat2 );
    virtual void CopyTo ( lMatrix<ValueType> &mat2 ) const;

    virtual lMatrix<ValueType> &operator= ( const lMatrix<ValueType> &mat2 );

    virtual void CopyStructureFrom ( const lMatrix<ValueType> &mat2 );
    virtual void CopyStructureTo ( lMatrix<ValueType> &mat2 ) const;

    virtual lMatrix<ValueType> *CloneWithoutContent ( ) const;

    virtual void ConvertFrom ( const lMatrix<ValueType> &mat2 );

    virtual void Init ( const int init_nnz,
                        const int init_num_row,
                        const int init_num_col,
                        const std::string init_name );
    virtual void Clear ( void );

    virtual void Zeros ( void );

    virtual lMatrix<ValueType> *extract_submatrix ( const int start_row, const int start_col, const int end_row, const int end_col ) const;

#    ifdef WITH_OPENCL
    virtual void set_thread_blocks ( cl_device_id my_device, const int size );
#    endif

    virtual void VectorMultAdd ( const lVector<ValueType> &invec,
                                 lVector<ValueType> *outvec ) const;
    virtual void VectorMultAdd ( const OPENCL_lVector<ValueType> &invec,
                                 OPENCL_lVector<ValueType> *outvec ) const;

    virtual void VectorMult ( const lVector<ValueType> &invec,
                              lVector<ValueType> *outvec ) const;
    virtual void VectorMult ( const OPENCL_lVector<ValueType> &invec,
                              OPENCL_lVector<ValueType> *outvec ) const;

    virtual void CastFrom ( const lMatrix<double>& other );
    virtual void CastFrom ( const lMatrix<float>& other );

    virtual void CastTo ( lMatrix<double>& other ) const;
    virtual void CastTo ( lMatrix<float>& other ) const;

    /// Stores the diagonal elements as first element per row.
    virtual void SwapDiagElementsToRowFront ( void );

    virtual void VectorMultNoDiag ( const lVector<ValueType> &in,
                                    lVector<ValueType> *out ) const;

  protected:
    size_t local_threads;
    size_t global_threads;
    bool manager_initialized;

#    ifdef WITH_OPENCL
    cl_mem val;
    cl_mem col;
    cl_mem row;
    opencl_manager *my_manager;
#    endif

};

#endif
