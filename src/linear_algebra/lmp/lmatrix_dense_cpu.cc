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

/// @author Niels Wegh, Simon Gawlok, Martin Wlotzka

#include "lmatrix_dense_cpu.h"

#include "lmp_mem.h"
#include "init_vec_mat.h"
#include "lmp_log.h"

#include <iostream>
#include <cassert>
#include <cstdlib>

using namespace hiflow::la;

// class CPU_DenselMatrix

/// standard constructor

template <typename ValueType>
CPU_DENSE_lMatrix<ValueType>::CPU_DENSE_lMatrix ( )
{
    this->platform_name_ = "CPU";
    this->platform_id_ = CPU;

    this->matrix.val = NULL;
}

/// destructor

template <typename ValueType>
CPU_DENSE_lMatrix<ValueType>::~CPU_DENSE_lMatrix ( )
{
    this->Clear ( );
}

/// Clears allocated memory.

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Clear ( )
{
    if ( this->matrix.val != NULL )
        delete [] this->matrix.val;
    this->matrix.val = NULL;

    this->nnz_ = 0;
    this->num_row_ = 0;
    this->num_col_ = 0;
    this->name_ = "";
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Zeros ( void )
{
    memsethost ( this->matrix.val, 0, this->get_nnz ( ), sizeof (ValueType ) );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Init ( const int init_nnz,
                                          const int init_num_row,
                                          const int init_num_col,
                                          const std::string init_name )
{

    assert ( init_nnz == init_num_row * init_num_col );
    assert ( init_num_row >= 0 );
    assert ( init_num_col >= 0 );

    this->Clear ( );

    // allocate
    this->matrix.val = new ValueType[init_nnz];
    assert ( this->matrix.val != NULL );

    this->name_ = init_name;
    this->nnz_ = init_nnz;
    this->num_row_ = init_num_row;
    this->num_col_ = init_num_col;

    this->Zeros ( );

}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Scale ( const ValueType alpha )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::Scale() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ScaleOffdiag ( const ValueType alpha )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ScaleOffdiag() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
lMatrix<ValueType> &CPU_DENSE_lMatrix<ValueType>::operator= ( const lMatrix<ValueType> &mat2 )
{

    if ( this == &mat2 )
        return *this;

    this->CopyFrom ( mat2 );
    return *this;

}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::CopyFrom ( const lMatrix<ValueType> &mat2 )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::CopyFrom() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::CopyTo ( lMatrix<ValueType> &mat2 ) const
{
    mat2.CopyFrom ( *this );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::CopyStructureFrom ( const lMatrix<ValueType> &mat2 )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::CopyStructureFrom() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::CopyStructureTo ( lMatrix<ValueType> &mat2 ) const
{
    mat2.CopyStructureFrom ( *this );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ConvertFrom ( const hiflow::la::lMatrix<ValueType> &mat2 )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ConvertFrom() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Psgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::Psgauss_seidel() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::BlocksPsgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec,
                                                          const int num_blocks ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::BlocksPsgauss_seidel() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::BlockPsgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec,
                                                         const int start_i, const int end_i ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::BlockPsgauss_seidel() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Pssor ( const ValueType omega, const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::Pssor() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Psor ( const ValueType omega, const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::Psor() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Pgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::Pgauss_seidel() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Pjacobi ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::Pjacobi() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::TransformFromCOO ( const int * rows,
                                                      const int * cols,
                                                      const ValueType * data,
                                                      const int num_rows,
                                                      const int num_cols,
                                                      const int num_nonzeros )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::TransformFromCOO() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ZeroRows ( const int *index_set,
                                              const int size,
                                              const ValueType alpha )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ZeroRows() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ZeroCols ( const int *index_set,
                                              const int size,
                                              const ValueType alpha )
{

    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ZeroCols() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::init_structure ( const int *rows, const int *cols )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::init_structure() only interface - not implemented yet" );
    exit ( -1 );
}

// Adds a value to the (i,j)-th element of the matrix.
// @param row - row index (i)
// @param col - col index (j)
// @param val - value added to the (i,j)-th element
// Note
// the function does not check if the pair (row,col) exists in the structure!

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::add_value ( const int row, const int col, ValueType val )
{

    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    this->matrix.val[row * this->get_num_col ( ) + col] += val;

}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::add_values ( const int* rows, int num_rows, const int* cols, int num_cols, const ValueType* values )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );
    for ( int i = 0; i < num_rows; ++i )
    {
        const int row_offset = rows[i] * this->get_num_col ( );
        const int val_offset = i*num_cols;
        for ( int j = 0; j < num_cols; ++j )
        {
            this->matrix.val[row_offset + cols[j]] += values[val_offset + j];
        }
    }
}

// Returns the (i,j)-th element of the matrix.
// param row - row index (i)
// param col - col index (j)
// return val - the (i,j)-th element
// Note
// the function does not check if the pair (row,col) exists in the structure!

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::get_value ( const int row, const int col, ValueType *val ) const
{

    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    *val = this->matrix.val[row * this->get_num_col ( ) + col];

}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::get_add_values ( const int* rows,
                                                    int num_rows,
                                                    const int* cols,
                                                    int num_cols,
                                                    const int* cols_target,
                                                    int num_cols_target,
                                                    ValueType* values ) const
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int i = 0; i < num_rows; ++i )
    {
        const int offset1 = i*num_cols_target;
        const int offset2 = rows[i] * this->get_num_col ( );
        for ( int j = 0; j < num_cols; ++j )
        {
            values[offset1 + cols_target[j]] += this->matrix.val[offset2 + cols[j]];
        }
    }
}

// Matrix vector multiplication.
// param invec - the matrix is multiplied with this vector
// return outvec - the resulting vector
// Note
// the function does not check if the pair (row,col) exists in the structure!

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::VectorMultAdd ( const lVector<ValueType> &invec,
                                                   lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::VectorMultAdd unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMultAdd ( *casted_invec, casted_outvec );

}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::VectorMultAdd ( const CPU_lVector<ValueType> &invec,
                                                   CPU_lVector<ValueType> *outvec ) const
{

    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( ( outvec->get_size ( ) == this->get_num_row ( ) ) ||
             ( invec .get_size ( ) == 0 ) );

    if ( this->get_nnz ( ) > 0 )
        for ( int i = 0; i<this->get_num_row ( ); ++i )
        {
            const int offset = i * this->get_num_col ( );
            for ( int j = 0; j<this->get_num_col ( ); ++j )
            {
                outvec->buffer[i] += this->matrix.val[offset + j] * invec.buffer[ j ];
            }
        }

}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::CastFrom ( const hiflow::la::lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( const CPU_DENSE_lMatrix<double>* other_cpu_dense
         = dynamic_cast < const CPU_DENSE_lMatrix<double>* > ( &other ) )
    {
        this->CastFrom ( *other_cpu_dense );
    }
    else
    {
        LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::CastFrom<double> called with non-CPU_DENSE matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::CastFrom ( const hiflow::la::lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( const CPU_DENSE_lMatrix<float>* other_cpu_dense
         = dynamic_cast < const CPU_DENSE_lMatrix<float>* > ( &other ) )
    {
        this->CastFrom ( *other_cpu_dense );
    }
    else
    {
        LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::CastFrom<float> called with non-CPU_DENSE matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::CastTo ( hiflow::la::lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( CPU_DENSE_lMatrix<double>* other_cpu_dense
         = dynamic_cast < CPU_DENSE_lMatrix<double>* > ( &other ) )
    {
        this->CastTo ( *other_cpu_dense );
    }
    else
    {
        LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::CastTo<double> called with non-CPU_DENSE matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::CastTo ( hiflow::la::lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( CPU_DENSE_lMatrix<float>* other_cpu_dense
         = dynamic_cast < CPU_DENSE_lMatrix<float>* > ( &other ) )
    {
        this->CastTo ( *other_cpu_dense );
    }
    else
    {
        LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::CastTo<double> called with non-CPU_DENSE matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::VectorMultNoDiag ( const hiflow::la::lVector<ValueType>& in,
                                                      hiflow::la::lVector<ValueType>* out ) const
{
    assert ( in.get_size ( ) == this->get_num_col ( ) );
    assert ( out->get_size ( ) == this->get_num_row ( ) );

    const CPU_lVector<ValueType>* cpu_in = dynamic_cast < const CPU_lVector<ValueType>* > ( &in );
    CPU_lVector<ValueType>* cpu_out = dynamic_cast < CPU_lVector<ValueType>* > ( out );

    if ( cpu_in && cpu_out )
    {
        this->VectorMultNoDiag ( *cpu_in, cpu_out );
    }
    else
    {
        LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::VectorMultNoDiag called with non-CPU vector argument." );
        this->print ( );
        in.print ( );
        out->print ( );
    }
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::VectorMultNoDiag ( const CPU_lVector<ValueType>& in,
                                                      CPU_lVector<ValueType>* out ) const
{
    assert ( in.get_size ( ) == this->get_num_col ( ) );
    assert ( out->get_size ( ) == this->get_num_row ( ) );

    for ( int i = 0; i < this->get_num_row ( ); ++i )
    {
        out->buffer[i] = 0.0;
        const int offset = i * this->get_num_col ( );

        for ( int j = 0; j < this->get_num_col ( ); ++j )
        {
            if ( i != j )
            {
                out->buffer[i] += this->matrix.val[offset + j] * in.buffer[j];
            }
        }
    }
}

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_DENSE_lMatrix<ValueType>::MatrixMult ( const lMatrix<ValueType> &inmat ) const
{

    const CPU_DENSE_lMatrix<ValueType> *casted_inmat = dynamic_cast < const CPU_DENSE_lMatrix<ValueType>* > ( &inmat );

    if ( casted_inmat == NULL )
    {
        LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::MatrixMult unsupported inmatrix" );
        this->print ( );
        inmat.print ( );
        exit ( -1 );
    }

    return (this->MatrixMult ( *casted_inmat ) );

}

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_DENSE_lMatrix<ValueType>::MatrixMult ( const CPU_DENSE_lMatrix<ValueType> &inmat ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::MatrixMult() only interface - not implemented yet" );
    exit ( -1 );
    return NULL;
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::MatrixAdd ( const lMatrix<ValueType> &inmat )
{

    const CPU_DENSE_lMatrix<ValueType> *casted_inmat = dynamic_cast < const CPU_DENSE_lMatrix<ValueType>* > ( &inmat );

    if ( casted_inmat == NULL )
    {
        LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::MatrixAdd unsupported inmatrix" );
        this->print ( );
        inmat.print ( );
        exit ( -1 );
    }

    this->MatrixAdd ( *casted_inmat );

}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::MatrixAdd ( const CPU_DENSE_lMatrix<ValueType> &inmat )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::MatrixAdd() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::GershgorinSpectrum ( ValueType *lambda_min, ValueType *lambda_max ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::GershgorinSpectrum() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ReadFile ( const char * filename )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ReadFile() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::WriteFile ( const char * filename ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::WriteFile() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_DENSE_lMatrix<ValueType>::extract_submatrix ( const int start_row, const int start_col,
                                                                                  const int end_row, const int end_col ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::extract_submatrix() only interface - not implemented yet" );
    exit ( -1 );
    return NULL;
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Reorder ( const int *index )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::Reorder() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ilu0 ( void )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ilu0() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ilup ( const int p )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ilup() not yet implemented" );
    this->print ( );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ilusp ( const int p, const int ncolors, const int *color_sizes, const int *permut_index )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ilusp() not yet implemented" );
    this->print ( );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::ilu_solve ( const hiflow::la::lVector<ValueType> &invec,
                                               hiflow::la::lVector<ValueType> *outvec ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::ilu_solve() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::Multicoloring ( int &ncolors, int **color_sizes, int **permut_index ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::Multicoloring() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::extract_diagelements ( const int start_i, const int end_i, lVector<ValueType> *vec ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::extract_diagelements() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::extract_invdiagelements ( const int start_i, const int end_i, lVector<ValueType> *vec ) const
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::extract_invdiagelements() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::compress_me ( )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::compress_me() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::copyLtoU ( )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::copyLtoU() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::delete_diagonal ( )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::delete_diagonal() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::delete_offdiagonal ( )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::delete_offdiagonal() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::delete_lower_triangular ( )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::delete_lower_triangular() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::delete_strictly_lower_triangular ( )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::delete_strictly_lower_triangular() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::delete_upper_triangular ( )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::delete_upper_triangular() only interface - not implemented yet" );
    exit ( -1 );
}

template <typename ValueType>
void CPU_DENSE_lMatrix<ValueType>::delete_strictly_upper_triangular ( )
{
    LOG_ERROR ( "CPU_DENSE_lMatrix<ValueType>::delete_strictly_upper_triangular() only interface - not implemented yet" );
    exit ( -1 );
}

template class CPU_DENSE_lMatrix<float>;
template class CPU_DENSE_lMatrix<double>;
