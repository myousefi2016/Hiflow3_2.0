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

/// @author Dimitar Lukarski, Martin Wlotzka

#ifndef __LMATRIX_CSR_GPU_H
#    define __LMATRIX_CSR_GPU_H

#    include <iostream>
#    include <stdlib.h>

#    include "lvector_gpu.h"
#    include "../lmatrix_csr_cpu.h"
#    include "../lmatrix_formats.h"
#    include "../lvector.h"
#    include "../lvector_cpu.h"
#    include "../lmatrix.h"
#    include "../lmatrix_csr.h"

/// @brief Provides the base matrix (CSR) class for GPU
/// @author Dimitar Lukarski
///
/// CPU_CSR_lMatrix maintains the matrix;
/// access functions and copy operators (to and from a cpu matrix)

template <typename ValueType>
class GPU_CSR_lMatrix : public CSR_lMatrix<ValueType>
{
  public:

    GPU_CSR_lMatrix ( );
    virtual ~GPU_CSR_lMatrix ( );

    hiflow::la::CSR_lMatrixType<ValueType> matrix; // it will be allocated on the device
    // in side:
    //  ValueType *val; // value
    //  int *col; // col index
    //  int *row; // row pointer

    // only copy
    virtual hiflow::la::lMatrix<ValueType> &operator= ( const hiflow::la::lMatrix<ValueType> &mat2 );
    virtual void CopyFrom ( const hiflow::la::lMatrix<ValueType> &mat2 );
    virtual void CopyTo ( hiflow::la::lMatrix<ValueType> &mat2 ) const;
    virtual void CopyStructureFrom ( const hiflow::la::lMatrix<ValueType> &mat2 );
    virtual void CopyStructureTo ( hiflow::la::lMatrix<ValueType> &mat2 ) const;
    virtual void ConvertFrom ( const hiflow::la::lMatrix<ValueType> &mat2 );

    virtual void Init ( const int init_nnz,
                        const int init_num_row,
                        const int init_num_col,
                        const std::string init_name );
    virtual void Clear ( void );

    virtual void Zeros ( void );

    virtual void set_thread_blocks ( void );
    virtual void set_thread_blocks ( int thread_block, int thread_block_size );

    virtual void Sync ( void ) const;

    virtual void CastFrom ( const hiflow::la::lMatrix<double>& other );
    virtual void CastFrom ( const hiflow::la::lMatrix<float>& other );

    virtual void CastFrom ( const GPU_CSR_lMatrix<double>& other );
    virtual void CastFrom ( const GPU_CSR_lMatrix<float>& other );

    virtual void CastFrom ( const CPU_CSR_lMatrix<double>& other );
    virtual void CastFrom ( const CPU_CSR_lMatrix<float>& other );

    virtual void CastTo ( hiflow::la::lMatrix<double>& other ) const;
    virtual void CastTo ( hiflow::la::lMatrix<float>& other ) const;

    virtual void CastTo ( GPU_CSR_lMatrix<double>& other ) const;
    virtual void CastTo ( GPU_CSR_lMatrix<float>& other ) const;

    virtual void CastTo ( CPU_CSR_lMatrix<double>& other ) const;
    virtual void CastTo ( CPU_CSR_lMatrix<float>& other ) const;

    virtual void SwapDiagElementsToRowFront ( void );

  protected:
    int thread_block_size_;
    int thread_block_;

};

/// @brief Provides GPU/Cuda CSR Scalar
/// Matrix-Vector Multiplication
/// @author Dimitar Lukarski

template <typename ValueType>
class GPUscalar_CSR_lMatrix : public GPU_CSR_lMatrix<ValueType>
{
  public:

    GPUscalar_CSR_lMatrix ( );
    /// The constructor call the Init() function to the do allocation
    /// of the matrix
    GPUscalar_CSR_lMatrix ( int init_nnz,
                            int init_num_row,
                            int init_num_col,
                            std::string init_name );
    virtual ~GPUscalar_CSR_lMatrix ( );

    virtual void VectorMultAdd ( const hiflow::la::lVector<ValueType> &invec,
                                 hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMultAdd ( const GPU_lVector<ValueType> &invec,
                                 GPU_lVector<ValueType> *outvec ) const;

    virtual void VectorMult ( const hiflow::la::lVector<ValueType> &invec,
                              hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMult ( const GPU_lVector<ValueType> &invec,
                              GPU_lVector<ValueType> *outvec ) const;

    virtual void VectorMultNoDiag ( const GPU_lVector<ValueType> &in,
                                    GPU_lVector<ValueType> *out ) const;

    virtual void VectorMultNoDiag ( const hiflow::la::lVector<ValueType> &in,
                                    hiflow::la::lVector<ValueType> *out ) const;

    virtual void BlocksPsgauss_seidel ( const hiflow::la::lVector<ValueType> &invec,
                                        hiflow::la::lVector<ValueType> *outvec,
                                        const int num_blocks ) const;

};

/// @brief Provides GPU/Cuda CSR Scalar + Texture Caching
/// Matrix-Vector Multiplication
/// @author Dimitar Lukarski

template <typename ValueType>
class GPUscalartex_CSR_lMatrix : public GPU_CSR_lMatrix<ValueType>
{
  public:

    GPUscalartex_CSR_lMatrix ( );
    /// The constructor call the Init() function to the do allocation
    /// of the matrix
    GPUscalartex_CSR_lMatrix ( int init_nnz,
                               int init_num_row,
                               int init_num_col,
                               std::string init_name );
    virtual ~GPUscalartex_CSR_lMatrix ( );

    virtual void VectorMultAdd ( const hiflow::la::lVector<ValueType> &invec,
                                 hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMultAdd ( const GPU_lVector<ValueType> &invec,
                                 GPU_lVector<ValueType> *outvec ) const;

    virtual void VectorMult ( const hiflow::la::lVector<ValueType> &invec,
                              hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMult ( const GPU_lVector<ValueType> &invec,
                              GPU_lVector<ValueType> *outvec ) const;

    virtual void VectorMultNoDiag ( const GPU_lVector<ValueType> &invec,
                                    GPU_lVector<ValueType> *outvec ) const;

    virtual void VectorMultNoDiag ( const hiflow::la::lVector<ValueType> &in,
                                    hiflow::la::lVector<ValueType> *out ) const;

};

#endif
