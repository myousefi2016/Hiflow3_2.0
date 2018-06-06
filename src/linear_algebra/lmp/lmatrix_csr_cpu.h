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

#ifndef __LMATRIX_CSR_CPU_H
#    define __LMATRIX_CSR_CPU_H

#    include <iostream>
#    include <stdlib.h>

#    include "lmatrix_formats.h"
#    include "lvector_cpu.h"
#    include "lmatrix_csr.h"

/// @brief Provides the base matrix (CSR) class for CPU
/// @author Dimitar Lukarski
///
/// CPU_CSR_lMatrix maintains the matrix;
/// provides access functions, copy operators,
/// transformation from different cpu matrix types (COO),
/// basic implementation of matrix-free-preconditioner,
/// sequential implementation of the VectorMultAdd fct

template <typename ValueType>
class CPU_CSR_lMatrix : public CSR_lMatrix<ValueType>
{
  public:

    /// Similar to standard library
    typedef ValueType value_type;

    CPU_CSR_lMatrix ( );
    virtual ~CPU_CSR_lMatrix ( );

    // @@ Bernd Doser, 2015-10-16: Data should be private.
    hiflow::la::CSR_lMatrixType<ValueType> matrix;
    // in side:
    //  ValueType *val; // value
    //  int *col; // col index
    //  int *row; // row pointer

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

    virtual void CastFrom ( const CPU_CSR_lMatrix<double>& other ) = 0;
    virtual void CastFrom ( const CPU_CSR_lMatrix<float>& other ) = 0;

    virtual void CastTo ( hiflow::la::lMatrix<double>& other ) const;
    virtual void CastTo ( hiflow::la::lMatrix<float>& other ) const;

    virtual void CastTo ( CPU_CSR_lMatrix<double>& other ) const = 0;
    virtual void CastTo ( CPU_CSR_lMatrix<float>& other ) const = 0;

    virtual void SwapDiagElementsToRowFront ( void );

    virtual void VectorMultNoDiag ( const hiflow::la::lVector<ValueType>& in,
                                    hiflow::la::lVector<ValueType>* out ) const;
    virtual void VectorMultNoDiag ( const CPU_lVector<ValueType>& in,
                                    CPU_lVector<ValueType>* out ) const;

    virtual void Init ( const int init_nnz,
                        const int init_num_row,
                        const int init_num_col,
                        const std::string init_name );
    virtual void Clear ( void );

    virtual void ZeroRows ( const int *index_set,
                            const int size,
                            const ValueType alpha );

    virtual void ZeroCols ( const int *index_set,
                            const int size,
                            const ValueType alpha );

    virtual void Zeros ( void );

    /// Compress matrix by removing zero elements
    virtual void Compress ( void );

    virtual void Reorder ( const int *index );

    virtual void Multicoloring ( int &ncolors, int **color_sizes, int **permut_index ) const;

    virtual void Levelscheduling ( int &nlevels, int **level_sizes, int **permut_index ) const;

    virtual void ilu0 ( void );

    virtual void ilup ( const int p );
    virtual void ilup ( const hiflow::la::lMatrix<ValueType> &mat, const int p );

    virtual void ilu_solve ( const hiflow::la::lVector<ValueType> &invec,
                             hiflow::la::lVector<ValueType> *outvec ) const;

    virtual void Scale ( const ValueType alpha );

    virtual void ScaleOffdiag ( const ValueType alpha );

    virtual void init_structure ( const int *rows,
                                  const int *cols );
    virtual void add_value ( const int row,
                             const int col,
                             const ValueType val );
    virtual void get_value ( const int row,
                             const int col,
                             ValueType *val ) const;

    virtual void add_values ( const int* rows,
                              int num_rows,
                              const int* cols,
                              int num_cols,
                              const ValueType* values );

    virtual void get_add_values ( const int* rows,
                                  int num_rows,
                                  const int* cols,
                                  int num_cols,
                                  const int* cols_target,
                                  int num_cols_target,
                                  ValueType* values ) const;

    virtual void VectorMultAdd_submatrix ( const int* rows,
                                           int num_rows,
                                           const int* cols,
                                           int num_cols,
                                           const int* cols_input,
                                           const ValueType* in_values,
                                           ValueType* out_values ) const;

    virtual void VectorMultAdd_submatrix_vanka ( const int* rows,
                                                 int num_rows,
                                                 const CPU_lVector< ValueType > &invec,
                                                 ValueType* out_values ) const;

    virtual void VectorMultAdd_submatrix_vanka ( const int* rows,
                                                 int num_rows,
                                                 const hiflow::la::lVector< ValueType > &invec,
                                                 ValueType* out_values ) const;

    virtual void transpose_me ( void );
    virtual void compress_me ( void );
    virtual bool issymmetric ( void );

    virtual void delete_diagonal ( void );
    virtual void delete_offdiagonal ( void );
    virtual void delete_lower_triangular ( void );
    virtual void delete_strictly_lower_triangular ( void );
    virtual void delete_upper_triangular ( void );
    virtual void delete_strictly_upper_triangular ( void );

    virtual hiflow::la::lMatrix<ValueType> *extract_submatrix ( const int start_row, const int start_col,
                                                                const int end_row, const int end_col ) const;

    virtual void extract_diagelements ( const int start_i, const int end_i, hiflow::la::lVector<ValueType> *vec ) const;
    virtual void extract_invdiagelements ( const int start_i, const int end_i, hiflow::la::lVector<ValueType> *vec ) const;

    virtual void VectorMultAdd ( const hiflow::la::lVector<ValueType> &invec,
                                 hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMultAdd ( const CPU_lVector<ValueType> &invec,
                                 CPU_lVector<ValueType> *outvec ) const;

    virtual hiflow::la::lMatrix<ValueType> *MatrixMult ( const hiflow::la::lMatrix<ValueType> &inmat ) const;
    virtual hiflow::la::lMatrix<ValueType> *MatrixMult ( const CPU_CSR_lMatrix<ValueType> &inmat ) const;

    virtual hiflow::la::lMatrix<ValueType> *MatrixMultSupStructure ( const hiflow::la::lMatrix<ValueType> &inmat ) const;
    virtual hiflow::la::lMatrix<ValueType> *MatrixMultSupStructure ( const CPU_CSR_lMatrix<ValueType> &inmat ) const;

    virtual hiflow::la::lMatrix<ValueType> *MatrixSupSPower ( const int p ) const;

    virtual void MatrixAdd ( const hiflow::la::lMatrix<ValueType> &inmat );
    virtual void MatrixAdd ( const CPU_CSR_lMatrix<ValueType> &inmat );

    virtual void GershgorinSpectrum ( ValueType *lambda_min, ValueType *lambda_max ) const;

    virtual void ReadFile ( const char* filename );
    virtual void WriteFile ( const char* filename ) const;

    virtual void Pjacobi ( const hiflow::la::lVector<ValueType> &invec,
                           hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void Pgauss_seidel ( const hiflow::la::lVector<ValueType> &invec,
                                 hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void Psgauss_seidel ( const hiflow::la::lVector<ValueType> &invec,
                                  hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void BlockPsgauss_seidel ( const hiflow::la::lVector<ValueType> &invec,
                                       hiflow::la::lVector<ValueType> *outvec,
                                       const int start_i,
                                       const int end_i ) const;

    virtual void BlocksPsgauss_seidel ( const hiflow::la::lVector<ValueType> &invec,
                                        hiflow::la::lVector<ValueType> *outvec,
                                        const int num_blocks ) const;

    virtual void Psor ( const ValueType omega,
                        const hiflow::la::lVector<ValueType> &invec,
                        hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void Pssor ( const ValueType omega,
                         const hiflow::la::lVector<ValueType> &invec,
                         hiflow::la::lVector<ValueType> *outvec ) const;

    /// Transform a local COO matrix to a local CSR matrix;
    /// the local matrix should be initialized;
    /// The elements afterward are sorted
    /// @param rows - the row index set
    /// @param cols - the column index set
    /// @param data - the values of the coo matrix
    /// @param num_rows - the number of rows in the coo matrix
    /// @param num_cols - the number of cols in the coo matrix
    /// @param num_nnz - the number of nnz in the coo matrix
    virtual void TransformFromCOO ( const int * rows,
                                    const int * cols,
                                    const ValueType * data,
                                    const int num_rows,
                                    const int num_cols,
                                    const int num_nonzeros );
};

/// @brief Provides CPU naive/simple sequential implementation
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUsimple_CSR_lMatrix : public CPU_CSR_lMatrix<ValueType>
{
  public:

    CPUsimple_CSR_lMatrix ( );
    /// The constructor call the Init() function to do the allocation
    /// of the matrix
    CPUsimple_CSR_lMatrix ( int init_nnz,
                            int init_num_row,
                            int init_num_col,
                            std::string init_name );
    virtual ~CPUsimple_CSR_lMatrix ( );

    virtual void VectorMult ( const hiflow::la::lVector<ValueType> &invec,
                              hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMult ( const CPU_lVector<ValueType> &invec,
                              CPU_lVector<ValueType> *outvec ) const;

    virtual void CastFrom ( const CPU_CSR_lMatrix<double>& other );
    virtual void CastFrom ( const CPU_CSR_lMatrix<float>& other );

    virtual void CastTo ( CPU_CSR_lMatrix<double>& other ) const;
    virtual void CastTo ( CPU_CSR_lMatrix<float>& other ) const;
};

/// @brief Provides CPU OpenMP parallel implementation
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUopenmp_CSR_lMatrix : public CPU_CSR_lMatrix<ValueType>
{
  public:

    CPUopenmp_CSR_lMatrix ( );
    CPUopenmp_CSR_lMatrix ( int init_nnz,
                            int init_num_row,
                            int init_num_col,
                            std::string init_name );
    virtual ~CPUopenmp_CSR_lMatrix ( );

    virtual void CloneFrom ( const hiflow::la::lMatrix<ValueType> &other );

    virtual void VectorMultAdd ( const hiflow::la::lVector<ValueType> &invec,
                                 hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMultAdd ( const CPU_lVector<ValueType> &invec,
                                 CPU_lVector<ValueType> *outvec ) const;

    virtual void VectorMult ( const hiflow::la::lVector<ValueType> &invec,
                              hiflow::la::lVector<ValueType> *outvec ) const;
    virtual void VectorMult ( const CPU_lVector<ValueType> &invec,
                              CPU_lVector<ValueType> *outvec ) const;

    virtual void BlocksPsgauss_seidel ( const hiflow::la::lVector<ValueType> &invec,
                                        hiflow::la::lVector<ValueType> *outvec,
                                        const int num_blocks ) const;

    virtual void add_values ( const int* rows,
                              int num_rows,
                              const int* cols,
                              int num_cols,
                              const ValueType* values );

    virtual void get_add_values ( const int* rows,
                                  int num_rows,
                                  const int* cols,
                                  int num_cols,
                                  const int* cols_target,
                                  int num_cols_target,
                                  ValueType* values ) const;

    virtual void VectorMultAdd_submatrix ( const int* rows,
                                           int num_rows,
                                           const int* cols,
                                           int num_cols,
                                           const int* cols_input,
                                           const ValueType* in_values,
                                           ValueType* out_values ) const;

    virtual void VectorMultAdd_submatrix_vanka ( const int* rows,
                                                 int num_rows,
                                                 const CPU_lVector< ValueType > &invec,
                                                 ValueType* out_values ) const;

    virtual void VectorMultAdd_submatrix_vanka ( const int* rows,
                                                 int num_rows,
                                                 const hiflow::la::lVector< ValueType > &invec,
                                                 ValueType* out_values ) const;

    virtual void set_num_threads ( void );
    virtual void set_num_threads ( int num_thread );

    virtual int num_threads ( void ) const
    {
        return this->num_threads_;
    }

    virtual void CastFrom ( const CPU_CSR_lMatrix<double>& other );
    virtual void CastFrom ( const CPU_CSR_lMatrix<float>& other );

    virtual void CastTo ( CPU_CSR_lMatrix<double>& other ) const;
    virtual void CastTo ( CPU_CSR_lMatrix<float>& other ) const;

    virtual void VectorMultNoDiag ( const CPU_lVector<ValueType>& in,
                                    CPU_lVector<ValueType>* out ) const;

  protected:
    int num_threads_;

};

/// @brief Provides wrapper to CPU Intel/MKL implementation
/// @author Dimitar Lukarski

template <typename ValueType>
class CPUmkl_CSR_lMatrix : public CPU_CSR_lMatrix<ValueType>
{
  public:

    CPUmkl_CSR_lMatrix ( );
    /// The constructor call the Init() function to the do allocation
    /// of the matrix
    CPUmkl_CSR_lMatrix ( int init_nnz,
                         int init_num_row,
                         int init_num_col,
                         std::string init_name );
    virtual ~CPUmkl_CSR_lMatrix ( );

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

    virtual void CastFrom ( const CPU_CSR_lMatrix<double>& other );
    virtual void CastFrom ( const CPU_CSR_lMatrix<float>& other );

    virtual void CastTo ( CPU_CSR_lMatrix<double>& other ) const;
    virtual void CastTo ( CPU_CSR_lMatrix<float>& other ) const;

  protected:
    int num_threads_;

};

#endif
