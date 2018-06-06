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

#include "lvector_cpu.h"
#include "lmatrix_coo_cpu.h"
#include "lmatrix_csr_cpu.h"
#include "lmp_mem.h"
#include "lmp_log.h"

extern "C"
{
#include "mmio.h"
}

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

using namespace hiflow::la;

const int DEBUG_LEVEL = 1;

// class CPU_COO_lMatrix

template <typename ValueType>
CPU_COO_lMatrix<ValueType>::CPU_COO_lMatrix ( )
{
    this->platform_name_ = "CPU (x86)";
    this->platform_id_ = CPU;

    this->matrix.val = NULL;
    this->matrix.col = NULL;
    this->matrix.row = NULL;
}

template <typename ValueType>
CPU_COO_lMatrix<ValueType>::~CPU_COO_lMatrix ( )
{
    this->Clear ( );
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Clear ( )
{
    if ( this->matrix.val != NULL )
        delete [] this->matrix.val;
    this->matrix.val = NULL;

    if ( this->matrix.col != NULL )
        delete [] this->matrix.col;
    this->matrix.col = NULL;

    if ( this->matrix.row != NULL )
        delete [] this->matrix.row;
    this->matrix.row = NULL;

    this->nnz_ = 0;
    this->num_row_ = 0;
    this->num_col_ = 0;
    this->name_ = "";
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::get_as_coo ( std::vector<ValueType>& val,
                                              std::vector<int>& col,
                                              std::vector<int>& row ) const
{
    val.resize ( this->nnz_ );
    row.resize ( this->nnz_ );
    col.resize ( this->nnz_ );

    memcpy ( &( val[0] ), this->matrix.val, this->nnz_ * sizeof (ValueType ) );
    memcpy ( &( col[0] ), this->matrix.col, this->nnz_ * sizeof (int ) );
    memcpy ( &( row[0] ), this->matrix.row, this->nnz_ * sizeof (int ) );
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Zeros ( void )
{
    memsethost ( this->matrix.val, 0, this->get_nnz ( ), sizeof (ValueType ) );
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Init ( const int init_nnz,
                                        const int init_num_row,
                                        const int init_num_col,
                                        const std::string init_name )
{

    assert ( init_nnz >= 0 );
    assert ( init_num_row >= 0 );
    assert ( init_num_col >= 0 );

    this->Clear ( );

    // allocate
    this->matrix.val = new ValueType[init_nnz];
    assert ( this->matrix.val != NULL );

    this->matrix.col = new int[init_nnz];
    assert ( this->matrix.col != NULL );

    this->matrix.row = new int[init_nnz];
    assert ( this->matrix.row != NULL );

    this->name_ = init_name;
    this->nnz_ = init_nnz;
    this->num_row_ = init_num_row;
    this->num_col_ = init_num_col;

    this->Zeros ( );

}

// This function is adapted from Matrix Market
// http://math.nist.gov/MatrixMarket/mmio-c.html

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::ReadFile ( const char * filename )
{
    this->Clear ( );

    this->name_.append ( filename );

    FILE *fid;
    MM_typecode matcode;

    fid = fopen ( filename, "r" );

    if ( fid == NULL )
    {
        printf ( "Unable to open file %s\n", filename );
        exit ( 1 );
    }

    if ( mm_read_banner ( fid, &matcode ) != 0 )
    {
        printf ( "Could not process lMatrix Market banner.\n" );
        exit ( 1 );
    }

    if ( !mm_is_valid ( matcode ) )
    {
        printf ( "Invalid lMatrix Market file.\n" );
        exit ( 1 );
    }

    if ( !( ( mm_is_real ( matcode ) || mm_is_integer ( matcode ) || mm_is_pattern ( matcode ) ) && mm_is_coordinate ( matcode ) && mm_is_sparse ( matcode ) ) )
    {
        printf ( "Sorry, this application does not support " );
        printf ( "Market Market type: [%s]\n", mm_typecode_to_str ( matcode ) );
        printf ( "Only sparse real-valued or pattern coordinate matrices are supported\n" );
        exit ( 1 );
    }

    int num_rows, num_cols, num_nonzeros;
    if ( mm_read_mtx_crd_size ( fid, &num_rows, &num_cols, &num_nonzeros ) != 0 )
        exit ( 1 );

    this->num_row_ = num_rows;
    this->num_col_ = num_cols;
    this->nnz_ = num_nonzeros;

    this->matrix.col = new int[this->get_nnz ( )];
    assert ( this->matrix.col != NULL );

    this->matrix.row = new int[this->get_nnz ( )];
    assert ( this->matrix.row != NULL );

    this->matrix.val = new ValueType[this->get_nnz ( )];
    assert ( this->matrix.val != NULL );

    LOG_INFO ( "fileio", "Reading sparse matrix from file " << filename );

    if ( mm_is_real ( matcode ) || mm_is_integer ( matcode ) )
    {
        for ( int i = 0; i < this->get_nnz ( ); ++i )
        {
            int ROW, COL;
            double VAL; // always read in a double and convert later if necessary

            if ( fscanf ( fid, " %d %d %lf \n", &ROW, &COL, &VAL ) != EOF )
            {

                this->matrix.row[i] = ROW - 1;
                this->matrix.col[i] = COL - 1;
                this->matrix.val[i] = static_cast < ValueType > ( VAL );
            }
            else
            {
                break;
            }
        }
    }
    else
    {
        printf ( "Unrecognized data type\n" );
        exit ( 1 );
    }

    LOG_INFO ( "fileio", "Reading sparse matrix from file done." );

    if ( mm_is_symmetric ( matcode ) )
    { //duplicate off diagonal entries
        int off_diagonals = 0;
        for ( int i = 0; i < this->get_nnz ( ); ++i )
        {
            if ( matrix.row[i] != matrix.col[i] )
                ++off_diagonals;
        }

        int true_nonzeros = 2 * off_diagonals + ( this->get_nnz ( ) - off_diagonals );

        int* new_row = new int[true_nonzeros];
        int* new_col = new int[true_nonzeros];
        ValueType* new_val = new ValueType[true_nonzeros];

        int ptr = 0;
        for ( int i = 0; i < this->get_nnz ( ); ++i )
        {
            if ( matrix.row[i] != matrix.col[i] )
            {
                new_row[ptr] = matrix.row[i];
                new_col[ptr] = matrix.col[i];
                new_val[ptr] = matrix.val[i];
                ptr++;
                new_col[ptr] = matrix.row[i];
                new_row[ptr] = matrix.col[i];
                new_val[ptr] = matrix.val[i];
                ptr++;
            }
            else
            {
                new_row[ptr] = matrix.row[i];
                new_col[ptr] = matrix.col[i];
                new_val[ptr] = matrix.val[i];
                ptr++;
            }
        }

        delete [] this->matrix.row;
        delete [] this->matrix.col;
        delete [] this->matrix.val;

        matrix.row = new_row;
        matrix.col = new_col;
        matrix.val = new_val;

        this->nnz_ = true_nonzeros;
    } //end symmetric case

}

// This function is adapted from Matrix Market
// http://math.nist.gov/MatrixMarket/mmio-c.html

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::WriteFile ( const char * filename ) const
{

    FILE *fid;
    MM_typecode matcode;

    printf ( "Writing sparse matrix to file (%s):", filename );
    fflush ( stdout );

    fid = fopen ( filename, "w" );

    if ( fid == NULL )
    {
        printf ( "Unable to open file %s\n", filename );
        exit ( 1 );
    }

    int i;

    mm_initialize_typecode ( &matcode );
    mm_set_matrix ( &matcode );
    mm_set_coordinate ( &matcode );
    mm_set_real ( &matcode );

    mm_write_banner ( fid, matcode );
    mm_write_mtx_crd_size ( fid, this->get_num_row ( ), this->get_num_col ( ), this->get_nnz ( ) ); // row, col, nnz

    /* NOTE: matrix market files use 1-based indices, i.e. first element
      of a vector has index 1, not 0.  */

    for ( i = 0; i<this->get_nnz ( ); ++i )
        fprintf ( fid, "%d %d %24.14e\n", this->matrix.col[i] + 1, this->matrix.row[i] + 1, this->matrix.val[i] );

    fclose ( fid );

    printf ( " done\n" );

}

template <typename ValueType>
lMatrix<ValueType> &CPU_COO_lMatrix<ValueType>::operator= ( const lMatrix<ValueType> &mat2 )
{

    if ( this == &mat2 )
        return *this;

    this->CopyFrom ( mat2 );
    return *this;

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::CopyFrom ( const lMatrix<ValueType> &mat2 )
{

    if ( this != &mat2 )
    {

        // CPU COO = CPU COO
        if ( const CPU_COO_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const CPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            if ( this->get_nnz ( ) > 0 )
            {
                memcpyhost ( this->matrix.val, casted_mat->matrix.val, this->get_nnz ( ) );
            }

        }
        else
        {

            // if mat is not a CPU matrix
            mat2.CopyTo ( *this );

        }
    }

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::CopyTo ( lMatrix<ValueType> &mat2 ) const
{
    mat2.CopyFrom ( *this );
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::CastFrom ( const hiflow::la::lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( const COO_lMatrix<double>* other_coo
         = dynamic_cast < const COO_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in COO format
        if ( const CPU_COO_lMatrix<double>* other_cpu_coo
             = dynamic_cast < const CPU_COO_lMatrix<double>* > ( other_coo ) )
        {
            // CPU from CPU
            this->CastFrom ( *other_cpu_coo );
        }
        else
        {
            // CPU from non-CPU via non-CPU to CPU
            other.CastTo ( *this );
        }
    }
    else
    {
        LOG_ERROR ( "CPU_COO_lMatrix<ValueType>::CastFrom<double> called with non-COO matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::CastFrom ( const hiflow::la::lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( const COO_lMatrix<double>* other_coo
         = dynamic_cast < const COO_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in COO format
        if ( const CPU_COO_lMatrix<float>* other_cpu_coo
             = dynamic_cast < const CPU_COO_lMatrix<float>* > ( other_coo ) )
        {
            // CPU from CPU
            this->CastFrom ( *other_cpu_coo );
        }
        else
        {
            // CPU from non-CPU via non-CPU to CPU
            other.CastTo ( *this );
        }
    }
    else
    {
        LOG_ERROR ( "CPU_COO_lMatrix<ValueType>::CastFrom<float> called with non-COO matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::CastTo ( hiflow::la::lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( COO_lMatrix<double>* other_coo
         = dynamic_cast < COO_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in COO format
        if ( CPU_COO_lMatrix<double>* other_cpu_coo
             = dynamic_cast < CPU_COO_lMatrix<double>* > ( other_coo ) )
        {
            // CPU to CPU
            this->CastTo ( *other_cpu_coo );
        }
        else
        {
            // CPU to non-CPU via non-CPU from CPU
            other.CastFrom ( *this );
        }
    }
    else
    {
        LOG_ERROR ( "CPU_COO_lMatrix<ValueType>::CastTo<double> called with non-COO matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::CastTo ( hiflow::la::lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( COO_lMatrix<double>* other_coo
         = dynamic_cast < COO_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in COO format
        if ( CPU_COO_lMatrix<float>* other_cpu_coo
             = dynamic_cast < CPU_COO_lMatrix<float>* > ( &other ) )
        {
            // CPU to CPU
            this->CastTo ( *other_cpu_coo );
        }
        else
        {
            // CPU to non-CPU via non-CPU from CPU
            other.CastFrom ( *this );
        }
    }
    else
    {
        LOG_ERROR ( "CPU_COO_lMatrix<ValueType>::CastTo<float> called with non-COO matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::CopyStructureFrom ( const lMatrix<ValueType> &mat2 )
{

    if ( this != &mat2 )
    {

        this->Init ( mat2.get_nnz ( ), mat2.get_num_row ( ), mat2.get_num_col ( ), mat2.get_name ( ) );

        // CPU COO = CPU COO
        if ( const CPU_COO_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const CPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            if ( this->get_nnz ( ) > 0 )
            {
                memcpyhost ( this->matrix.col, casted_mat->matrix.col, this->get_nnz ( ) );
                memcpyhost ( this->matrix.row, casted_mat->matrix.row, this->get_nnz ( ) );
            }

        }
        else
        {

            // if mat is not a CPU matrix
            mat2.CopyStructureTo ( *this );
        }
    }
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::CopyStructureTo ( lMatrix<ValueType> &mat2 ) const
{
    mat2.CopyStructureFrom ( *this );
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::ConvertFrom ( const lMatrix<ValueType> &mat2 )
{

    // CPU COO = CSR CPU
    if ( const CPU_CSR_lMatrix<ValueType> *casted_mat =
         dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &mat2 ) )
    {

        this->Init ( mat2.get_nnz ( ), mat2.get_num_row ( ), mat2.get_num_col ( ), mat2.get_name ( ) );

        assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
        assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
        assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

        if ( this->get_nnz ( ) > 0 )
            this->TransformFromCSR ( casted_mat->matrix.row, casted_mat->matrix.col, casted_mat->matrix.val,
                                     this->get_num_row ( ), this->get_num_col ( ), this->get_nnz ( ) );

    }
    else
    {
        LOG_ERROR ( "CPU_COO_lMatrix<ValueType>::ConvertFrom; unsupported matrix type" );
        this->print ( );
        mat2.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Psgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_coo_matrix<ValueType>::Psgauss_seidel unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    assert ( casted_invec ->get_size ( ) > 0 );
    assert ( casted_outvec->get_size ( ) > 0 );
    assert ( casted_invec ->get_size ( ) == this->get_num_col ( ) );
    assert ( casted_outvec->get_size ( ) == this->get_num_row ( ) );

    // M = (D+L) D^-1 (D+R)
    // Mz=r
    // (D+L) y = r
    // forward step

    int last_j = 0;

    for ( int i = 0; i<this->get_num_row ( ); ++i )
    {

        casted_outvec->buffer[i] = casted_invec->buffer[i];

        for ( int j = 0; j<this->get_nnz ( ); ++j )
        {

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] == this->matrix.row[j] ) )
                last_j = j;

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] < this->matrix.row[j] ) )
                casted_outvec->buffer[i] -= this->matrix.val[j] * casted_outvec->buffer[ this->matrix.col[j] ];

        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

    // Dy
    for ( int i = 0; i<this->get_nnz ( ); ++i )
        if ( this->matrix.row[i] == this->matrix.col[i] )
        {
            assert ( this->matrix.val[i] != 0.0 );
            casted_outvec->buffer[this->matrix.row[i]] *= this->matrix.val[i];
        }

    // backward
    // (D+R)r = Dy

    for ( int i = this->get_num_row ( ) - 1; i >= 0; --i )
    {

        for ( int j = 0; j<this->get_nnz ( ); ++j )
        {

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] == this->matrix.row[j] ) )
                last_j = j;

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] > this->matrix.row[j] ) )
                casted_outvec->buffer[i] -= this->matrix.val[j] * casted_outvec->buffer[ this->matrix.col[j] ];

        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Pssor ( const ValueType omega, const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_coo_matrix<ValueType>::Pssor unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    assert ( casted_invec ->get_size ( ) > 0 );
    assert ( casted_outvec->get_size ( ) > 0 );
    assert ( casted_invec ->get_size ( ) == this->get_num_col ( ) );
    assert ( casted_outvec->get_size ( ) == this->get_num_row ( ) );

    // M = (D+wL) D^-1 (D+wR)
    // Mz=r

    // (D+wL) y = r
    // forward step

    int last_j = 0;
    for ( int i = 0; i<this->get_num_row ( ); ++i )
    {

        casted_outvec->buffer[i] = casted_invec->buffer[i];

        for ( int j = 0; j<this->get_nnz ( ); ++j )
        {

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] == this->matrix.row[j] ) )
                last_j = j;

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] < this->matrix.row[j] ) )
                casted_outvec->buffer[i] -= omega * this->matrix.val[j] * casted_outvec->buffer[ this->matrix.col[j] ];

        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

    // Dy
    for ( int i = 0; i<this->get_nnz ( ); ++i )
        if ( this->matrix.row[i] == this->matrix.col[i] )
        {
            assert ( this->matrix.val[i] != 0.0 );
            casted_outvec->buffer[this->matrix.row[i]] *= this->matrix.val[i];
        }

    // backward
    // (D+Rw)r = Dy

    for ( int i = this->get_num_row ( ) - 1; i >= 0; --i )
    {

        for ( int j = 0; j<this->get_nnz ( ); ++j )
        {

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] == this->matrix.row[j] ) )
                last_j = j;

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] > this->matrix.row[j] ) )
                casted_outvec->buffer[i] -= omega * this->matrix.val[j] * casted_outvec->buffer[ this->matrix.col[j] ];

        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Psor ( const ValueType omega, const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_coo_matrix<ValueType>::Psor unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    assert ( casted_invec ->get_size ( ) > 0 );
    assert ( casted_outvec->get_size ( ) > 0 );
    assert ( casted_invec ->get_size ( ) == this->get_num_col ( ) );
    assert ( casted_outvec->get_size ( ) == this->get_num_row ( ) );

    // M = (D/omega + L)
    int diag_j = 0;
    for ( int i = 0; i<this->get_num_row ( ); ++i )
    {

        casted_outvec->buffer[i] = casted_invec->buffer[i];

        for ( int j = 0; j<this->get_nnz ( ); ++j )
        {

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] == this->matrix.row[j] ) )
                diag_j = j;

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] < this->matrix.row[j] ) )
                casted_outvec->buffer[i] -= this->matrix.val[j] * casted_outvec->buffer[ this->matrix.col[j] ];

        }

        assert ( this->matrix.val[diag_j] != 0.0 );
        casted_outvec->buffer[i] = casted_outvec->buffer[i] * omega / ( this->matrix.val[diag_j] );

    }

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Pgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_coo_matrix<ValueType>::Pgauss_seidel unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    assert ( casted_invec ->get_size ( ) > 0 );
    assert ( casted_outvec->get_size ( ) > 0 );
    assert ( casted_invec ->get_size ( ) == this->get_num_col ( ) );
    assert ( casted_outvec->get_size ( ) == this->get_num_row ( ) );

    // M = D+L
    // Mz=r

    //  Mz=b ;
    //  for (i = 0; i < n; ++i)
    //  {
    //    z[i] = r[i];
    //
    //    for (j = 0; j < i; ++j)
    //      z[i] -= M[i, j] * z[j];
    //
    //    z[i] /= M[i, i];
    //  }

    int last_j = 0;
    for ( int i = 0; i<this->get_num_row ( ); ++i )
    {

        casted_outvec->buffer[i] = casted_invec->buffer[i];

        for ( int j = 0; j<this->get_nnz ( ); ++j )
        {

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] == this->matrix.row[j] ) )
                last_j = j;

            if ( ( this->matrix.row[j] == i ) && ( this->matrix.col[j] < this->matrix.row[j] ) )
                casted_outvec->buffer[i] -= this->matrix.val[j] * casted_outvec->buffer[ this->matrix.col[j] ];

        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Pjacobi ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_coo_matrix<ValueType>::Pjacobi unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    assert ( casted_invec ->get_size ( ) > 0 );
    assert ( casted_outvec->get_size ( ) > 0 );
    assert ( casted_invec ->get_size ( ) == this->get_num_col ( ) );
    assert ( casted_outvec->get_size ( ) == this->get_num_row ( ) );

    for ( int i = 0; i<this->get_nnz ( ); ++i )
        if ( this->matrix.row[i] == this->matrix.col[i] )
        {
            assert ( this->matrix.val[i] != 0.0 );
            casted_outvec->buffer[this->matrix.row[i]] = casted_invec->buffer[this->matrix.row[i]] / this->matrix.val[i];
        }

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::TransformFromCSR ( const int * Ap,
                                                    const int * Aj,
                                                    const ValueType * Ax,
                                                    const int num_rows,
                                                    const int num_cols,
                                                    const int num_nonzeros )
{

    assert ( this->get_nnz ( ) == num_nonzeros );
    assert ( this->get_num_row ( ) == num_rows );
    assert ( this->get_num_col ( ) == num_cols );

    for ( int i = 0; i < num_rows; ++i )
    {
        int row_start = Ap[i];
        int row_end = Ap[i + 1];
        for ( int jj = row_start; jj < row_end; ++jj )
        {
            this->matrix.row[jj] = i;
        }
    }

    for ( int i = 0; i < num_nonzeros; ++i )
    {
        this->matrix.col[i] = Aj[i];
        this->matrix.val[i] = Ax[i];
    }

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::VectorMultAdd ( const lVector<ValueType> &invec,
                                                 lVector<ValueType> *outvec ) const
{

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPU_COO_lMatrix<ValueType>::VectorMultAdd unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMultAdd ( *casted_invec, casted_outvec );

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::VectorMultAdd ( const CPU_lVector<ValueType> &invec,
                                                 CPU_lVector<ValueType> *outvec ) const
{

    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( ( outvec->get_size ( ) == this->get_num_row ( ) ) ||
             ( invec .get_size ( ) == 0 ) );

    if ( this->get_nnz ( ) > 0 )
        for ( int i = 0; i<this->get_nnz ( ); ++i )
            outvec->buffer[this->matrix.row[i]] += this->matrix.val[i] * invec.buffer[this->matrix.col[i]];

}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::VectorMultNoDiag ( const hiflow::la::lVector<ValueType>& in,
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
        LOG_ERROR ( "CPU_COO_lMatrix<ValueType>::VectorMultNoDiag called with non-CPU vector argument." );
        this->print ( );
        in.print ( );
        out->print ( );
    }
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::VectorMultNoDiag ( const CPU_lVector<ValueType>& in,
                                                    CPU_lVector<ValueType>* out ) const
{
    assert ( in.get_size ( ) == this->get_num_col ( ) );
    assert ( out->get_size ( ) == this->get_num_row ( ) );

    for ( int i = 0; i < this->get_num_row ( ); ++i )
    {
        out->buffer[i] = 0.0;
    }

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        if ( this->matrix.row[i] != this->matrix.col[i] )
        {
            out->buffer[this->matrix.row[i]] += this->matrix.val[i] * in.buffer[this->matrix.col[i]];
        }
    }
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Sort ( void )
{

    ValueType tv;
    int ti;

    LOG_DEBUG ( 3, "CPU_COO_lMatrix<ValueType>::Sort() bubble sort - Sorting the cols...." );

    // bubble sort (by cols)
    for ( int i = 0; i<this->get_nnz ( ) - 1; ++i )
        for ( int j = 0; j<this->get_nnz ( ) - i - 1; ++j )
            if ( this->matrix.col[j] > this->matrix.col[j + 1] )
            {

                ti = this->matrix.col[j];
                this->matrix.col[j] = this->matrix.col[j + 1];
                this->matrix.col[j + 1] = ti;

                ti = this->matrix.row[j];
                this->matrix.row[j] = this->matrix.row[j + 1];
                this->matrix.row[j + 1] = ti;

                tv = this->matrix.val[j];
                this->matrix.val[j] = this->matrix.val[j + 1];
                this->matrix.val[j + 1] = tv;

            }

    LOG_DEBUG ( 3, "CPU_COO_lMatrix<ValueType>::Sort() bubble sort - Sorting the  rows...." );
    // bubble sort (by rows)
    for ( int i = 0; i<this->get_nnz ( ) - 1; ++i )
        for ( int j = 0; j<this->get_nnz ( ) - i - 1; ++j )
            if ( this->matrix.row[j] > this->matrix.row[j + 1] )
            {

                ti = this->matrix.col[j];
                this->matrix.col[j] = this->matrix.col[j + 1];
                this->matrix.col[j + 1] = ti;

                ti = this->matrix.row[j];
                this->matrix.row[j] = this->matrix.row[j + 1];
                this->matrix.row[j + 1] = ti;

                tv = this->matrix.val[j];
                this->matrix.val[j] = this->matrix.val[j + 1];
                this->matrix.val[j + 1] = tv;

            }

    LOG_DEBUG ( 3, "CPU_COO_lMatrix<ValueType>::Sort() bubble sort - Sorting - done..." );
}

template <typename ValueType>
void CPU_COO_lMatrix<ValueType>::Reorder ( const int *index )
{

    int *col, *row, *ind;

    col = new int[this->get_nnz ( )];
    assert ( col != NULL );

    row = new int[this->get_nnz ( )];
    assert ( row != NULL );

    ind = new int[this->get_num_row ( )];
    assert ( ind != NULL );

    // permute the index
    // from index dst -> index src
    for ( int i = 0; i<this->get_num_row ( ); ++i )
        ind[ index[i] ] = i;

    for ( int i = 0; i<this->get_nnz ( ); ++i )
    {
        col[i] = this->matrix.col[i];
        row[i] = this->matrix.row[i];
    }

    for ( int i = 0; i<this->get_nnz ( ); ++i )
    {
        this->matrix.col[ i ] = ind[col[i]];
        this->matrix.row[ i ] = ind[row[i]];
    }

    delete [] col;
    delete [] row;
    delete [] ind;

}

template class CPU_COO_lMatrix<float>;
template class CPU_COO_lMatrix<double>;
