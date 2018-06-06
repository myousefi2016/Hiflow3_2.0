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

#ifndef __LMATRIX_COO_GPU_H
#    define __LMATRIX_COO_GPU_H

#    include <iostream>
#    include <stdlib.h>

#    include "../lmatrix_coo_cpu.h"

/// @brief Provides the base matrix (COO) class for GPU
/// @author Dimitar Lukarski
///
/// CPU_COO_lMatrix maintains the matrix;
/// access functions and copy operators (to and from a cpu matrix)

template <typename ValueType>
class GPU_COO_lMatrix : public COO_lMatrix<ValueType>
{
  public:

    GPU_COO_lMatrix ( );
    virtual ~GPU_COO_lMatrix ( );

    hiflow::la::COO_lMatrixType<ValueType> matrix;
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

    virtual void CastFrom ( const GPU_COO_lMatrix<double>& other );
    virtual void CastFrom ( const GPU_COO_lMatrix<float>& other );

    virtual void CastFrom ( const CPU_COO_lMatrix<double>& other );
    virtual void CastFrom ( const CPU_COO_lMatrix<float>& other );

    virtual void CastTo ( hiflow::la::lMatrix<double>& other ) const;
    virtual void CastTo ( hiflow::la::lMatrix<float>& other ) const;

    virtual void CastTo ( GPU_COO_lMatrix<double>& other ) const;
    virtual void CastTo ( GPU_COO_lMatrix<float>& other ) const;

    virtual void CastTo ( CPU_COO_lMatrix<double>& other ) const;
    virtual void CastTo ( CPU_COO_lMatrix<float>& other ) const;

  protected:
    int thread_block_size_;
    int thread_block_;

};

#endif
