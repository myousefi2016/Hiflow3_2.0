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

#ifndef __LMATRIX_COO_CPU_H
#    define __LMATRIX_COO_CPU_H

#    include <iostream>
#    include <stdlib.h>

#    include "lmatrix_formats.h"
#    include "lvector_cpu.h"
#    include "lmatrix_coo.h"

/// @brief Provides the base matrix (COO) class for CPU
/// @author Dimitar Lukarski
///
/// CPU_COO_lMatrix maintains the matrix;
/// provides access functions, copy operators,
/// transformation from different cpu matrix types (CSR),
/// basic implementation of matrix-free-preconditioner

template <typename ValueType>
class CPU_COO_lMatrix : public COO_lMatrix<ValueType>
{
  public:

    CPU_COO_lMatrix ( );
    virtual ~CPU_COO_lMatrix ( );

    hiflow::la::COO_lMatrixType<ValueType> matrix;
    // in side:
    //  ValueType *val; // value
    //  int *col; // col index
    //  int *row; // row index

    virtual void get_as_coo ( std::vector<ValueType>& vals,
                              std::vector<int>& rows,
                              std::vector<int>& cols ) const;

    // copy and converting
    virtual hiflow::la::lMatrix<ValueType> &operator= ( const hiflow::la::lMatrix<ValueType> &mat2 );
    virtual void CopyFrom ( const hiflow::la::lMatrix<ValueType> &mat2 );
    virtual void CopyTo ( hiflow::la::lMatrix<ValueType> &mat2 ) const;
    virtual void CopyStructureFrom ( const hiflow::la::lMatrix<ValueType> &mat2 );
    virtual void CopyStructureTo ( hiflow::la::lMatrix<ValueType> &mat2 ) const;
    virtual void ConvertFrom ( const hiflow::la::lMatrix<ValueType> &mat2 );

    virtual void CastFrom ( const hiflow::la::lMatrix<double>& other );
    virtual void CastFrom ( const hiflow::la::lMatrix<float>& other );

    virtual void CastFrom ( const CPU_COO_lMatrix<double>& other ) = 0;
    virtual void CastFrom ( const CPU_COO_lMatrix<float>& other ) = 0;

    virtual void CastTo ( hiflow::la::lMatrix<double>& other ) const;
    virtual void CastTo ( hiflow::la::lMatrix<float>& other ) const;

    virtual void CastTo ( CPU_COO_lMatrix<double>& other ) const = 0;
    virtual void CastTo ( CPU_COO_lMatrix<float>& other ) const = 0;

    virtual void VectorMultNoDiag ( const hiflow::la::lVector<ValueType>& in,
                                    hiflow::la::lVector<ValueType>* out ) const;
    virtual void VectorMultNoDiag ( const CPU_lVector<ValueType>& in,
                                    CPU_lVector<ValueType>* out ) const;

    virtual void Init ( const int init_nnz,
                        const int init_num_row,
                        const int init_num_col,
                        const std::string init_name );
    virtual void Clear ( void );

    virtual void Zeros ( void );

    virtual void Sort ( void );

    virtual void Reorder ( const int *index );

    virtual void VectorMultAdd ( const hiflow::la::lVector<ValueType> &invec,
                                 hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMultAdd ( const CPU_lVector<ValueType> &invec,
                                 CPU_lVector<ValueType> *outvec ) const;

    virtual void ReadFile ( const char* filename );
    virtual void WriteFile ( const char* filename ) const;

    virtual void Pjacobi ( const hiflow::la::lVector<ValueType> &invec,
                           hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void Pgauss_seidel ( const hiflow::la::lVector<ValueType> &invec,
                                 hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void Psgauss_seidel ( const hiflow::la::lVector<ValueType> &invec,
                                  hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void Psor ( const ValueType omega,
                        const hiflow::la::lVector<ValueType> &invec,
                        hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void Pssor ( const ValueType omega,
                         const hiflow::la::lVector<ValueType> &invec,
                         hiflow::la::lVector<ValueType> *outvec ) const;

    /// Transform a local CSR matrix to a local COO matrix;
    /// the local matrix should be initialized
    /// @param rows - the row pointers set
    /// @param cols - the column index set
    /// @param data - the values of the coo matrix
    /// @param num_rows - the number of rows in the coo matrix
    /// @param num_cols - the number of cols in the coo matrix
    /// @param num_nnz - the number of nnz in the coo matrix
    virtual void TransformFromCSR ( const int * Ap,
                                    const int * Aj,
                                    const ValueType * Ax,
                                    const int num_rows,
                                    const int num_cols,
                                    const int num_nonzeros );
};

/// @brief Provides CPU naive/simple sequential implementation
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUsimple_COO_lMatrix : public CPU_COO_lMatrix<ValueType>
{
  public:
    CPUsimple_COO_lMatrix ( );
    CPUsimple_COO_lMatrix ( int init_nnz,
                            int init_num_col,
                            int init_num_row,
                            std::string init_name );
    virtual ~CPUsimple_COO_lMatrix ( );

    virtual void VectorMult ( const hiflow::la::lVector<ValueType> &invec,
                              hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMult ( const CPU_lVector<ValueType> &invec,
                              CPU_lVector<ValueType> *outvec ) const;

    virtual void CastFrom ( const CPU_COO_lMatrix<double>& other );
    virtual void CastFrom ( const CPU_COO_lMatrix<float>& other );

    virtual void CastTo ( CPU_COO_lMatrix<double>& other ) const;
    virtual void CastTo ( CPU_COO_lMatrix<float>& other ) const;
};

/// @brief Provides CPU OpenMP parallel implementation
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUopenmp_COO_lMatrix : public CPU_COO_lMatrix<ValueType>
{
  public:
    CPUopenmp_COO_lMatrix ( );
    CPUopenmp_COO_lMatrix ( int init_nnz,
                            int init_num_row,
                            int init_num_col,
                            std::string init_name );
    virtual ~CPUopenmp_COO_lMatrix ( );

    virtual void CloneFrom ( const hiflow::la::lMatrix<ValueType> &other );

    virtual void VectorMultAdd ( const hiflow::la::lVector<ValueType> &invec,
                                 hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMultAdd ( const CPU_lVector<ValueType> &invec,
                                 CPU_lVector<ValueType> *outvec ) const;

    virtual void VectorMult ( const hiflow::la::lVector<ValueType> &invec,
                              hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMult ( const CPU_lVector<ValueType> &invec,
                              CPU_lVector<ValueType> *outvec ) const;

    virtual void set_num_threads ( void );
    virtual void set_num_threads ( int num_thread );

    virtual int num_threads ( void ) const
    {
        return this->num_threads_;
    }

    virtual void CastFrom ( const CPU_COO_lMatrix<double>& other );
    virtual void CastFrom ( const CPU_COO_lMatrix<float>& other );

    virtual void CastTo ( CPU_COO_lMatrix<double>& other ) const;
    virtual void CastTo ( CPU_COO_lMatrix<float>& other ) const;

    virtual void VectorMultNoDiag ( const CPU_lVector<ValueType>& in,
                                    CPU_lVector<ValueType>* out ) const;

  protected:
    int num_threads_;

};

// @brief Provides wrapper to CPU Intel/MKL implementation
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUmkl_COO_lMatrix : public CPU_COO_lMatrix<ValueType>
{
  public:
    CPUmkl_COO_lMatrix ( );
    CPUmkl_COO_lMatrix ( int init_nnz,
                         int init_num_row,
                         int init_num_col,
                         std::string init_name );
    virtual ~CPUmkl_COO_lMatrix ( );

    virtual void CloneFrom ( const hiflow::la::lMatrix<ValueType> &other );

    virtual void VectorMult ( const hiflow::la::lVector<ValueType> &invec,
                              hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMult ( const CPU_lVector<ValueType> &invec,
                              CPU_lVector<ValueType> *outvec ) const;

    virtual void set_num_threads ( void );
    virtual void set_num_threads ( int num_thread );

    virtual int num_threads ( void ) const
    {
        return this->num_threads_;
    }

    virtual void CastFrom ( const CPU_COO_lMatrix<double>& other );
    virtual void CastFrom ( const CPU_COO_lMatrix<float>& other );

    virtual void CastTo ( CPU_COO_lMatrix<double>& other ) const;
    virtual void CastTo ( CPU_COO_lMatrix<float>& other ) const;

  protected:
    int num_threads_;

};

#endif
