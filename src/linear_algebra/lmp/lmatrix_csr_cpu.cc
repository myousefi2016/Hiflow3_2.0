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

/// @author Dimitar Lukarski, Simon Gawlok, Martin Wlotzka

#include "lvector_cpu.h"
#include "lmatrix_csr_cpu.h"
#include "lmatrix_coo_cpu.h"
#include "lmp_mem.h"
#include "init_vec_mat.h"
#include "lmp_log.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include "common/log.h"

using namespace hiflow::la;

// class CPU_CSR_lMatrix

const int DEBUG_LEVEL = 0;

template <typename ValueType>
CPU_CSR_lMatrix<ValueType>::CPU_CSR_lMatrix ( )
{
    this->platform_name_ = "CPU";
    this->platform_id_ = CPU;

    this->matrix.val = NULL;
    this->matrix.col = NULL;
    this->matrix.row = NULL;
}

template <typename ValueType>
CPU_CSR_lMatrix<ValueType>::~CPU_CSR_lMatrix ( )
{
    this->Clear ( );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Clear ( )
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
void CPU_CSR_lMatrix<ValueType>::get_as_coo ( std::vector<ValueType>& val,
                                              std::vector<int>& col,
                                              std::vector<int>& row ) const
{
    val.resize ( this->nnz_ );
    row.resize ( this->nnz_ );
    col.resize ( this->nnz_ );

    memcpy ( &( val[0] ), this->matrix.val, this->nnz_ * sizeof (ValueType ) );
    memcpy ( &( col[0] ), this->matrix.col, this->nnz_ * sizeof (int ) );

    for ( int i = 0; i < this->num_row_; ++i )
    {
        for ( int j = this->matrix.row[i]; j < this->matrix.row[i + 1]; ++j )
        {
            row[j] = i;
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Zeros ( void )
{
    memsethost ( this->matrix.val, 0, this->get_nnz ( ), sizeof (ValueType ) );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Init ( const int init_nnz,
                                        const int init_num_row,
                                        const int init_num_col,
                                        const std::string init_name )
{

    assert ( init_nnz >= 0 );
    assert ( init_num_row >= 0 );
    assert ( init_num_col >= 0 );

    this->Clear ( );

    LOG_DEBUG ( 2, "init_nnz = " << init_nnz << " init_num_row = " << init_num_row << " init_num_col = " << init_num_col );

    // allocate
    this->matrix.val = new ValueType[init_nnz];
    assert ( this->matrix.val != NULL );

    this->matrix.col = new int[init_nnz];
    assert ( this->matrix.col != NULL );

    this->matrix.row = new int[init_num_row + 1];
    assert ( this->matrix.row != NULL );

    this->matrix.row[init_num_row] = init_nnz;

    this->name_ = init_name;
    this->nnz_ = init_nnz;
    this->num_row_ = init_num_row;
    this->num_col_ = init_num_col;

    this->Zeros ( );

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Compress ( )
{
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<ValueType> values;

    rows.reserve ( this->get_num_row ( ) );
    cols.reserve ( this->get_nnz ( ) );
    values.reserve ( this->get_nnz ( ) );
    if ( this->get_nnz ( ) > 0 )
    {
        // determine new matrix structure and values
        rows.push_back ( 0 );
        int current_elem = 0;
        for ( int i = 0, e_i = this->get_num_row ( ); i != e_i; ++i )
        {
            for ( int j = this->matrix.row[i], e_j = this->matrix.row[i + 1]; j != e_j; ++j )
            {
                if ( this->matrix.val[j] != static_cast < ValueType > ( 0.0 ) )
                {
                    values.push_back ( this->matrix.val[j] );
                    cols.push_back ( this->matrix.col[j] );
                    ++current_elem;
                }
            }
            rows.push_back ( current_elem );
        }

        // set new matrix structure and values
        this->nnz_ = current_elem;

        // row pointer
        delete [] this->matrix.row;
        this->matrix.row = new int[rows.size ( )];
        assert ( this->matrix.row != NULL );

        for ( size_t i = 0, i_e = rows.size ( ); i != i_e; ++i )
        {
            this->matrix.row[i] = rows[i];
        }

        // col pointer
        delete [] this->matrix.col;
        this->matrix.col = new int[cols.size ( )];
        assert ( this->matrix.col != NULL );

        for ( size_t i = 0, i_e = cols.size ( ); i != i_e; ++i )
        {
            this->matrix.col[i] = cols[i];
        }

        // values
        delete [] this->matrix.val;
        this->matrix.val = new ValueType[values.size ( )];
        assert ( this->matrix.val != NULL );

        for ( size_t i = 0, i_e = values.size ( ); i != i_e; ++i )
        {
            this->matrix.val[i] = values[i];
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Scale ( const ValueType alpha )
{

    for ( size_t j = 0, j_e = this->get_nnz ( ); j != j_e; ++j )
        this->matrix.val[j] *= alpha;

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ScaleOffdiag ( const ValueType alpha )
{

    if ( this->get_nnz ( ) > 0 )
        for ( size_t i = 0, i_e = this->get_num_row ( ); i != i_e; ++i )
            for ( size_t j = this->matrix.row[i], j_e = this->matrix.row[i + 1]; j != j_e; ++j )
                if ( i != this->matrix.col[j] )
                    this->matrix.val[j] *= alpha;

}

template <typename ValueType>
lMatrix<ValueType> &CPU_CSR_lMatrix<ValueType>::operator= ( const lMatrix<ValueType> &mat2 )
{

    if ( this == &mat2 )
        return *this;

    this->CopyFrom ( mat2 );
    return *this;

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::CopyFrom ( const lMatrix<ValueType> &mat2 )
{

    if ( this != &mat2 )
    {

        // CPU CSR = CPU CSR
        if ( const CPU_CSR_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &mat2 ) )
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
void CPU_CSR_lMatrix<ValueType>::CopyTo ( lMatrix<ValueType> &mat2 ) const
{
    mat2.CopyFrom ( *this );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::CastFrom ( const hiflow::la::lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( const CSR_lMatrix<double>* other_csr
         = dynamic_cast < const CSR_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in CSR format
        if ( const CPU_CSR_lMatrix<double>* other_cpu_csr
             = dynamic_cast < const CPU_CSR_lMatrix<double>* > ( other_csr ) )
        {
            // CPU from CPU
            this->CastFrom ( *other_cpu_csr );
        }
        else
        {
            // CPU from non-CPU via non-CPU to CPU
            other.CastTo ( *this );
        }
    }
    else
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::CastFrom<double> called with non-CSR matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::CastFrom ( const hiflow::la::lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( const CSR_lMatrix<double>* other_csr
         = dynamic_cast < const CSR_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in CSR format
        if ( const CPU_CSR_lMatrix<float>* other_cpu_csr
             = dynamic_cast < const CPU_CSR_lMatrix<float>* > ( other_csr ) )
        {
            // CPU from CPU
            this->CastFrom ( *other_cpu_csr );
        }
        else
        {
            // CPU from non-CPU via non-CPU to CPU
            other.CastTo ( *this );
        }
    }
    else
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::CastFrom<float> called with non-CSR matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::CastTo ( hiflow::la::lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( CSR_lMatrix<double>* other_csr
         = dynamic_cast < CSR_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in CSR format
        if ( CPU_CSR_lMatrix<double>* other_cpu_csr
             = dynamic_cast < CPU_CSR_lMatrix<double>* > ( other_csr ) )
        {
            // CPU to CPU
            this->CastTo ( *other_cpu_csr );
        }
        else
        {
            // CPU to non-CPU via non-CPU from CPU
            other.CastFrom ( *this );
        }
    }
    else
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::CastTo<double> called with non-CSR matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::CastTo ( hiflow::la::lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( CSR_lMatrix<double>* other_csr
         = dynamic_cast < CSR_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in CSR format
        if ( CPU_CSR_lMatrix<float>* other_cpu_csr
             = dynamic_cast < CPU_CSR_lMatrix<float>* > ( &other ) )
        {
            // CPU to CPU
            this->CastTo ( *other_cpu_csr );
        }
        else
        {
            // CPU to non-CPU via non-CPU from CPU
            other.CastFrom ( *this );
        }
    }
    else
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::CastTo<float> called with non-CSR matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::CopyStructureFrom ( const lMatrix<ValueType> &mat2 )
{

    if ( this != &mat2 )
    {

        this->Init ( mat2.get_nnz ( ), mat2.get_num_row ( ), mat2.get_num_col ( ), mat2.get_name ( ) );

        // CPU CSR = CPU CSR
        if ( const CPU_CSR_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            if ( this->get_nnz ( ) > 0 )
            {
                memcpyhost ( this->matrix.col, casted_mat->matrix.col, this->get_nnz ( ) );
                memcpyhost ( this->matrix.row, casted_mat->matrix.row, this->get_num_row ( ) + 1 );
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
void CPU_CSR_lMatrix<ValueType>::CopyStructureTo ( lMatrix<ValueType> &mat2 ) const
{
    mat2.CopyStructureFrom ( *this );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ConvertFrom ( const hiflow::la::lMatrix<ValueType> &mat2 )
{

    // CPU CSR = COO CPU
    if ( const CPU_COO_lMatrix<ValueType> *casted_mat =
         dynamic_cast < const CPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
    {

        this->Init ( casted_mat->get_nnz ( ), casted_mat->get_num_row ( ), casted_mat->get_num_col ( ), casted_mat->get_name ( ) );

        assert ( this->get_nnz ( ) == casted_mat->get_nnz ( ) );
        assert ( this->get_num_row ( ) == casted_mat->get_num_row ( ) );
        assert ( this->get_num_col ( ) == casted_mat->get_num_col ( ) );

        if ( this->get_nnz ( ) > 0 )
            this->TransformFromCOO ( casted_mat->matrix.row, casted_mat->matrix.col, casted_mat->matrix.val,
                                     this->get_num_row ( ), this->get_num_col ( ), this->get_nnz ( ) );

    }
    else
    {
        // unsupported type
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::ConverFrom unsupported matrix type" );
        this->print ( );
        mat2.print ( );
        exit ( -1 );
    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::SwapDiagElementsToRowFront ( void )
{
    for ( int i = 0; i < this->get_num_row ( ); ++i )
    {
        for ( int j = this->matrix.row[i]; j < this->matrix.row[i + 1]; ++j )
        {
            if ( i == this->matrix.col[j] )
            {
                int j0 = this->matrix.row[i];

                std::swap ( this->matrix.val[j0], this->matrix.val[j] );

                std::swap ( this->matrix.col[j0], this->matrix.col[j] );

                break;
            }
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::VectorMultNoDiag ( const hiflow::la::lVector<ValueType>& in,
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
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::VectorMultNoDiag called with non-CPU vector argument." );
        this->print ( );
        in.print ( );
        out->print ( );
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::VectorMultNoDiag ( const CPU_lVector<ValueType>& in,
                                                    CPU_lVector<ValueType>* out ) const
{
    assert ( in.get_size ( ) == this->get_num_col ( ) );
    assert ( out->get_size ( ) == this->get_num_row ( ) );

    for ( int i = 0; i < this->get_num_row ( ); ++i )
    {
        ValueType res = 0.0;
        for ( int j = this->matrix.row[i]; j < this->matrix.row[i + 1]; ++j )
        {
            if ( i != this->matrix.col[j] )
            {
                res += this->matrix.val[j] * in.buffer[this->matrix.col[j]];
            }
        }
        out->buffer[i] = res;
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Psgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_csr_matrix<ValueType>::Psgauss_seidel unsupported in or out vector" );
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

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {
            if ( i > this->matrix.col[j] )
                casted_outvec->buffer[i] -= casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
            if ( i == this->matrix.col[j] )
            {
                last_j = j;
                break;
            }
        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

    // Dy
    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( i == this->matrix.col[j] )
            {
                casted_outvec->buffer[i] *= this->matrix.val[j];
                break;
            }

    // backward
    // (D+R)z = Dy
    for ( int i = this->get_num_row ( ) - 1; i >= 0; --i )
    {

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {

            if ( i == this->matrix.col[j] )
                last_j = j;
            if ( i < this->matrix.col[j] )
                casted_outvec->buffer[i] -= casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::BlocksPsgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec,
                                                        const int num_blocks ) const
{

    int start_i = 0;
    int end_i = 0;
    int step_i = invec.get_size ( ) / num_blocks;

    for ( int np = 0; np < num_blocks; ++np )
    {
        start_i = np*step_i;
        end_i = ( np + 1 ) * step_i;

        if ( np == num_blocks - 1 ) end_i = invec.get_size ( );

        this->BlockPsgauss_seidel ( invec, outvec, start_i, end_i );

    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::BlockPsgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec,
                                                       const int start_i, const int end_i ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_csr_matrix<ValueType>::BlockPsgauss_seidel unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    assert ( start_i >= 0 );
    assert ( end_i > 0 );
    assert ( end_i <= this->get_num_row ( ) );
    assert ( end_i > start_i );

    assert ( casted_invec ->get_size ( ) > 0 );
    assert ( casted_outvec->get_size ( ) > 0 );
    assert ( casted_invec ->get_size ( ) == this->get_num_col ( ) );
    assert ( casted_outvec->get_size ( ) == this->get_num_row ( ) );

    // M = (D+L) D^-1 (D+R)
    // Mz=r
    // (D+L) y = r
    // forward step

    int last_j = 0;
    for ( int i = start_i; i < end_i; ++i )
    {

        casted_outvec->buffer[i] = casted_invec->buffer[i];

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {

            if ( ( this->matrix.col[j] >= start_i ) &&
                 ( this->matrix.col[j] < end_i ) )
            {

                if ( i > this->matrix.col[j] )
                    casted_outvec->buffer[i] -= casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
                if ( i == this->matrix.col[j] )
                {
                    last_j = j;
                    break;
                }
            }

        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

    // Dy
    for ( int i = start_i; i < end_i; ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( i == this->matrix.col[j] )
            {
                casted_outvec->buffer[i] *= this->matrix.val[j];
                break;
            }

    // backward
    // (D+R)z = Dy
    for ( int i = end_i - 1; i >= start_i; --i )
    {

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {

            if ( ( this->matrix.col[j] >= start_i ) &&
                 ( this->matrix.col[j] < end_i ) )
            {

                if ( i == this->matrix.col[j] )
                    last_j = j;
                if ( i < this->matrix.col[j] )
                    casted_outvec->buffer[i] -= casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
            }

        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Pssor ( const ValueType omega, const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_csr_matrix<ValueType>::Pssor unsupported in or out vector" );
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

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {
            if ( i > this->matrix.col[j] )
                casted_outvec->buffer[i] -= omega * casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
            if ( i == this->matrix.col[j] )
            {
                last_j = j;
                break;
            }
        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

    // Dy
    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( i == this->matrix.col[j] )
            {
                casted_outvec->buffer[i] *= this->matrix.val[j];
                break;
            }

    // backward
    // (D+wR)z = Dy
    for ( int i = this->get_num_row ( ) - 1; i >= 0; --i )
    {

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {

            if ( i == this->matrix.col[j] )
                last_j = j;
            if ( i < this->matrix.col[j] )
                casted_outvec->buffer[i] -= omega * casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Psor ( const ValueType omega, const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_csr_matrix<ValueType>::Psor unsupported in or out vector" );
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
    int last_j = 0;
    for ( int i = 0; i<this->get_num_row ( ); ++i )
    {

        casted_outvec->buffer[i] = casted_invec->buffer[i];

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {
            if ( i > this->matrix.col[j] )
                casted_outvec->buffer[i] -= casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
            if ( i == this->matrix.col[j] )
            {
                last_j = j;
                break;
            }
        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] = casted_outvec->buffer[i] * omega / this->matrix.val[last_j];

    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Pgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{

    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_csr_matrix<ValueType>::Pgauss_seidel unsupported in or out vector" );
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

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {
            if ( i > this->matrix.col[j] )
                casted_outvec->buffer[i] -= casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
            if ( i == this->matrix.col[j] )
            {
                last_j = j;
                break;
            }
        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Pjacobi ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPUsimple_csr_matrix<ValueType>::Pjacobi unsupported matrix or in/out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );

    }

    assert ( casted_invec ->get_size ( ) > 0 );
    assert ( casted_outvec->get_size ( ) > 0 );
    assert ( casted_invec ->get_size ( ) == this->get_num_col ( ) );
    assert ( casted_outvec->get_size ( ) == this->get_num_row ( ) );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( i == this->matrix.col[j] )
            {
                assert ( this->matrix.val[j] );
                casted_outvec->buffer[i] = casted_invec->buffer[ this->matrix.col[j] ] / this->matrix.val[j];
            }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::TransformFromCOO ( const int * rows,
                                                    const int * cols,
                                                    const ValueType * data,
                                                    const int num_rows,
                                                    const int num_cols,
                                                    const int num_nonzeros )
{
    assert ( this->get_nnz ( ) == num_nonzeros );
    assert ( this->get_num_row ( ) == num_rows );
    assert ( this->get_num_col ( ) == num_cols );

    // seto to zero
    memsethost ( this->matrix.row, 0, num_rows + 1, sizeof (int ) );
    memsethost ( this->matrix.col, 0, num_nonzeros, sizeof (int ) );

    if ( data != NULL )
        memsethost ( this->matrix.val, 0, num_nonzeros, sizeof (ValueType ) );

    // build the row offsets
    for ( int i = 0; i != num_nonzeros; ++i )
        ++( this->matrix.row[rows[i]] );
    this->matrix.row[num_rows] = num_nonzeros;

    int *index_row = new int[num_rows]; // auxiliary index for accessing the CSR cols, vals

    // accumulate the row offsets
    int acum = 0;
    for ( int i = 0; i != num_rows; ++i )
    {
        int current_row = this->matrix.row[i];
        this->matrix.row[i] = acum;
        index_row[i] = acum; // = this->matrix.row[i]
        acum = acum + current_row;
    }

    // for init structure only - no data
    if ( data == NULL )
    {

        for ( int i = 0; i != num_nonzeros; ++i )
        {
            // inverse mapping
            this->matrix.col[ index_row[ rows[i] ] ] = cols[i];

            // keep the track of the current row index
            ++( index_row[ rows[i] ] ); // = index_row[ rows[i] ] + 1 ;
        }

    }
    else
    {

        for ( int i = 0; i != num_nonzeros; ++i )
        {
            const int temp = index_row[ rows[i] ];
            // inverse mapping
            this->matrix.col[ temp ] = cols[i];
            this->matrix.val[ temp ] = data[i];

            // keep the track of the current row index
            ++( index_row[ rows[i] ] ); // = index_row[ rows[i] ] + 1 ;
        }
    }

    delete[] index_row;

    // Sorting the Cols (per Row)
    ValueType tv;
    int ti;

    for ( int k = 0, k_e = this->get_num_row ( ); k != k_e; ++k )
    {
        for ( int i = this->matrix.row[k], i_e = this->matrix.row[k + 1]; i != i_e; ++i )
        {
            int current_min_index = i;
            int current_min = this->matrix.col[i];
            for ( int j = i + 1; j < i_e; ++j )
            {
                if ( this->matrix.col[j] < current_min )
                {
                    current_min_index = j;
                    current_min = this->matrix.col[j];
                }
            }
            ti = this->matrix.col[i];
            this->matrix.col[i] = this->matrix.col[current_min_index];
            this->matrix.col[current_min_index] = ti;

            tv = this->matrix.val[i];
            this->matrix.val[i] = this->matrix.val[current_min_index];
            this->matrix.val[current_min_index] = tv;
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ZeroRows ( const int *index_set,
                                            const int size,
                                            const ValueType alpha )
{

    // the index_set is not checked is it in the range
    for ( int i = 0; i != size; ++i )
    {
        const int current_index = index_set[i];
        for ( int j = this->matrix.row[current_index], j_e = this->matrix.row[current_index + 1]; j != j_e; ++j )
        {
            if ( current_index == this->matrix.col[j] )
            {
                this->matrix.val[j] = alpha;
            }
            else
            {
                this->matrix.val[j] = 0.0;
            }
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ZeroCols ( const int *index_set,
                                            const int size,
                                            const ValueType alpha )
{

    assert ( size > 0 );

    for ( int k = 0; k<this->get_num_row ( ); ++k )
    {
        for ( int j = this->matrix.row[k]; j<this->matrix.row[k + 1]; ++j )
        {
            for ( int i = 0; i < size; ++i )
            {
                if ( index_set[i] == this->matrix.col[j] )
                {
                    if ( k == index_set[i] )
                    {
                        this->matrix.val[j] = alpha;
                    }
                    else
                    {
                        this->matrix.val[j] = 0.0;
                    }
                }
            }
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::init_structure ( const int *rows, const int *cols )
{

    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    if ( DEBUG_LEVEL >= 2 )
    {
        for ( int l = 0; l<this->get_nnz ( ); ++l )
        {
            std::cout << "i: " << rows[l] << ", j: " << cols[l] << std::endl;
        }
    }

    this->TransformFromCOO ( rows, cols, NULL,
                             this->get_num_row ( ),
                             this->get_num_col ( ),
                             this->get_nnz ( ) );

}

// Note
// the function does not check is the pair (row,col) exsit in the structure!

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::add_value ( const int row, const int col, ValueType val )
{

    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( size_t j = this->matrix.row[row], e_j = this->matrix.row[row + 1]; j != e_j; ++j )
    {
        if ( col == this->matrix.col[j] )
        {
            this->matrix.val[j] += val;
            return;
        }
    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::add_values ( const int* rows, int num_rows, const int* cols, int num_cols, const ValueType* values )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int i = 0; i != num_rows; ++i )
    {
        int k = 0;
        const int row_offset = i*num_cols;
        for ( size_t j = this->matrix.row[rows[i]], j_e = this->matrix.row[rows[i] + 1]; k < num_cols && j != j_e; ++j )
        {
            if ( cols[k] == this->matrix.col[j] )
            {
                this->matrix.val[j] += values[row_offset + k];
                ++k;
            }
            else if ( cols[k] < this->matrix.col[j] )
            {
                while ( cols[k] < this->matrix.col[j] && k < num_cols )
                {
                    ++k;
                }
                if ( k >= num_cols )
                {
                    break;
                }
                if ( cols[k] == this->matrix.col[j] )
                {
                    this->matrix.val[j] += values[row_offset + k];
                    ++k;
                }
            }
        }
    }
}

// Note
// the function does not check is the pair (row,col) exists in the structure!

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::get_value ( const int row, const int col, ValueType *val ) const
{

    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int j = this->matrix.row[row]; j<this->matrix.row[row + 1]; ++j )
    {
        if ( col == this->matrix.col[j] )
        {
            *val = this->matrix.val[j];
            return;
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::get_add_values ( const int* rows,
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

    for ( int i = 0; i != num_rows; ++i )
    {
        int k = 0;
        const int row_offset = i*num_cols_target;
        for ( size_t j = this->matrix.row[rows[i]], j_e = this->matrix.row[rows[i] + 1]; k < num_cols && j != j_e; ++j )
        {
            if ( cols[k] == this->matrix.col[j] )
            {
                values[row_offset + cols_target[k]] += this->matrix.val[j];
                ++k;
            }
            else if ( cols[k] < this->matrix.col[j] )
            {
                while ( cols[k] < this->matrix.col[j] )
                {
                    ++k;
                }
                if ( k >= num_cols )
                {
                    break;
                }
                if ( cols[k] == this->matrix.col[j] )
                {
                    values[row_offset + cols_target[k]] += this->matrix.val[j];
                    ++k;
                }
            }
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::VectorMultAdd_submatrix ( const int* rows,
                                                           int num_rows,
                                                           const int* cols,
                                                           int num_cols,
                                                           const int* cols_input,
                                                           const ValueType* in_values,
                                                           ValueType* out_values ) const
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int i = 0; i < num_rows; ++i )
    {
        int k = 0;
        for ( size_t j = this->matrix.row[rows[i]]; k < num_cols && j<this->matrix.row[rows[i] + 1]; ++j )
        {
            if ( cols[k] == this->matrix.col[j] )
            {
                out_values[i] += this->matrix.val[j] * in_values[cols_input[k]];
                ++k;
            }
            else if ( cols[k] < this->matrix.col[j] )
            {
                while ( cols[k] < this->matrix.col[j] )
                {
                    ++k;
                }
                if ( k >= num_cols )
                {
                    break;
                }
                if ( cols[k] == this->matrix.col[j] )
                {
                    out_values[i] += this->matrix.val[j] * in_values[cols_input[k]];
                    ++k;
                }
            }
        }
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::VectorMultAdd_submatrix_vanka ( const int* rows,
                                                                 int num_rows,
                                                                 const CPU_lVector< ValueType > &invec,
                                                                 ValueType* out_values ) const
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int i = 0; i != num_rows; ++i )
    {
        ValueType res = 0.;
        for ( size_t j = this->matrix.row[rows[i]], j_e = this->matrix.row[rows[i] + 1]; j != j_e; ++j )
        {
            res += this->matrix.val[j] * invec.buffer[this->matrix.col[j]];
        }
        out_values[i] += res;
    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::VectorMultAdd_submatrix_vanka ( const int* rows,
                                                                 int num_rows,
                                                                 const hiflow::la::lVector< ValueType > &invec,
                                                                 ValueType* out_values ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );

    if ( casted_invec == NULL )
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::VectorMultAdd_submatrix_vanka unsupported in vector" );
        this->print ( );
        invec.print ( );
        exit ( -1 );
    }

    this->VectorMultAdd_submatrix_vanka ( rows, num_rows, *casted_invec, out_values );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::VectorMultAdd ( const lVector<ValueType> &invec,
                                                 lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::VectorMultAdd unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMultAdd ( *casted_invec, casted_outvec );

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::VectorMultAdd ( const CPU_lVector<ValueType> &invec,
                                                 CPU_lVector<ValueType> *outvec ) const
{

    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( ( outvec->get_size ( ) == this->get_num_row ( ) ) ||
             ( invec .get_size ( ) == 0 ) );

    if ( this->get_nnz ( ) > 0 )
    {
        for ( size_t i = 0, e_i = this->get_num_row ( ); i != e_i; ++i )
        {
            ValueType res = 0.;
            for ( size_t j = this->matrix.row[i], e_j = this->matrix.row[i + 1]; j != e_j; ++j )
            {
                res += this->matrix.val[j] * invec.buffer[ this->matrix.col[j] ];
            }
            outvec->buffer[i] += res;
        }
    }

}

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_CSR_lMatrix<ValueType>::MatrixMult ( const hiflow::la::lMatrix<ValueType> &inmat ) const
{

    const CPU_CSR_lMatrix<ValueType> *casted_inmat = dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &inmat );

    if ( casted_inmat == NULL )
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::MatrixMult unsupported inmatrix" );
        this->print ( );
        inmat.print ( );
        exit ( -1 );
    }

    return (this->MatrixMult ( *casted_inmat ) );

}

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_CSR_lMatrix<ValueType>::MatrixMult ( const CPU_CSR_lMatrix<ValueType> &inmat ) const
{
    assert ( inmat.get_num_row ( ) == this->get_num_col ( ) );

    std::string mat_name;
    mat_name = "Matrix Mult = ";
    mat_name.append ( this->get_name ( ) );
    mat_name.append ( " x " );
    mat_name.append ( inmat.get_name ( ) );

    LOG_DEBUG ( 0, mat_name );

    std::vector<int> rows, cols;

    int mat_nnz = 0;
    ValueType val = 0.0;
    int c = 0;

    LOG_DEBUG ( 2, "MatrixMatrix Multiplication row = " << this->get_num_row ( ) << " nnz = " << this->get_nnz ( ) );

    // determine structure of resulting matrix
    if ( inmat.get_num_row ( ) > 0 )
    {
        for ( int row = 0; row < this->get_num_row ( ); ++row )
        {
            for ( int col = 0; col < inmat.get_num_col ( ); ++col )
            {
                val = 0.0;

                for ( int i = this->matrix.row[row]; i < this->matrix.row[row + 1]; ++i )
                {
                    c = this->matrix.col[i];

                    for ( int j = inmat.matrix.row[c]; j < inmat.matrix.row[c + 1]; ++j )
                    {
                        if ( col == inmat.matrix.col[j] )
                        {
                            val += this->matrix.val[i] * inmat.matrix.val[j];
                            break;
                        }
                    }
                }

                if ( val != 0.0 )
                {
                    ++mat_nnz;
                    rows.push_back ( row );
                    cols.push_back ( col );
                }
            }

        }
    }

    hiflow::la::lMatrix<ValueType> *outmat;
    outmat = init_matrix<ValueType>( mat_nnz,
            this->get_num_row ( ),
            inmat.get_num_col ( ),
            mat_name,
            this->get_platform ( ),
            this->get_implementation ( ),
            this->get_matrix_format ( ) ); // CSR

    outmat->init_structure ( &( rows.front ( ) ), &( cols.front ( ) ) );

    outmat->Zeros ( );

    // compute matrix
    if ( inmat.get_num_row ( ) > 0 )
    {
        for ( int row = 0; row < this->get_num_row ( ); ++row )
        {
            for ( int col = 0; col < inmat.get_num_col ( ); ++col )
            {
                val = 0.0;

                for ( int i = this->matrix.row[row]; i < this->matrix.row[row + 1]; ++i )
                {
                    c = this->matrix.col[i];

                    for ( int j = inmat.matrix.row[c]; j < inmat.matrix.row[c + 1]; ++j )
                    {
                        if ( col == inmat.matrix.col[j] )
                        {
                            val += this->matrix.val[i] * inmat.matrix.val[j];
                            break;
                        }
                    }
                }

                if ( val != 0.0 )
                {
                    outmat->add_value ( row, col, val );
                }
            }
        }
    }

    outmat->compress_me ( );

    return outmat;
}

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_CSR_lMatrix<ValueType>::MatrixMultSupStructure ( const hiflow::la::lMatrix<ValueType> &inmat ) const
{

    const CPU_CSR_lMatrix<ValueType> *casted_inmat = dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &inmat );

    if ( casted_inmat == NULL )
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::MatrixMultSupStructure unsupported inmatrix" );
        this->print ( );
        inmat.print ( );
        exit ( -1 );
    }

    return (this->MatrixMultSupStructure ( *casted_inmat ) );

}

// Thanks to Felix Riehn

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_CSR_lMatrix<ValueType>::MatrixMultSupStructure ( const CPU_CSR_lMatrix<ValueType> &inmat ) const
{

    assert ( inmat.get_num_row ( ) == this->get_num_col ( ) );
    assert ( inmat.get_num_col ( ) == this->get_num_row ( ) );

    CPUsimple_CSR_lMatrix<ValueType> *outmat = new CPUsimple_CSR_lMatrix<ValueType>;

    std::string mat_name;
    mat_name = "Matrix Mult=";
    mat_name.append ( this->get_name ( ) );
    mat_name.append ( " x " );
    mat_name.append ( inmat.get_name ( ) );

    //------------------------------
    //multiply matrix pattern(this)
    //with external matrix(inmat)
    //according to R.E.Bank and C.C.Douglas
    // A(this)*B(inmat)=C(outmat)
    //------------------------------

    int nrow = this->get_num_row ( );

    LOG_DEBUG ( 1, "Building a MatrixMultSupStructure A" );

    LOG_DEBUG ( 2, "MatrixMultSupStructure number of nodes: " << nrow
                << " number of nonzeros: " << this->get_nnz ( ) );

    //vector container to store resulting matrix pattern
    std::vector<int> row, col;
    std::vector<int> index;

    index.resize ( nrow );
    row.resize ( nrow + 1 );
    col.resize ( this->get_nnz ( )*2.5 );

    //initial counters
    row[0] = 0;
    int k = 0, jlast = -1, length = 0, j = 0;

    for ( int h = 0; h < nrow; ++h )
    {
        index[h] = -2;
    }

    int new_nz = 0;
    //row loop
    for ( int i = 0; i < nrow; ++i )
    {
        if ( i % ( ( nrow - 1 - ( nrow - 1 ) % 10 ) / 10 ) == 0 )
        {
            LOG_DEBUG ( 3, "Building " << float(100 * float(( i + 1 ) / float(nrow - 1 ) ) ) << "%" );
        }
        jlast = -1;
        length = 0;
        //loop on nonzeros of A in row i
        for ( int nz_ii = this->matrix.row[i]; nz_ii<this->matrix.row[i + 1]; ++nz_ii )
        {
            k = this->matrix.col[nz_ii];
            //loop on nonzeros of B in row k
            for ( int nz_k = inmat.matrix.row[k]; nz_k < inmat.matrix.row[k + 1]; ++nz_k )
            {
                j = inmat.matrix.col[nz_k];
                //store new column positions
                if ( index[j]<-1 )
                {
                    index[j] = jlast;
                    jlast = j;
                    length++;
                }
            }
        }

        //fill findings to csr structure
        row[i + 1] = row[i] + length;

        if ( static_cast < int > ( col.size ( ) ) < new_nz + length )
        {
            col.resize ( col.size ( ) + this->get_nnz ( ) );
            col.resize ( col.size ( ) + this->get_nnz ( ) );
        }

        while ( jlast >= 0 )
        {
            length--;
            col[row[i] + length] = jlast;
            jlast = index[jlast];
            index[col[row[i] + length]] = -2;
            new_nz++;
        }

    }

    LOG_DEBUG ( 2, "New matrix pattern has " << new_nz << " non-zero entries" );
    //write to matrix
    outmat->Init ( new_nz, nrow, nrow, "C=AB" );
    outmat->matrix.row[nrow] = new_nz;
    for ( int i = 0; i < nrow; ++i )
        outmat->matrix.row[i] = row[i];
    outmat->matrix.row[nrow] = new_nz;

    for ( int h = 0; h < new_nz; ++h )
    {
        outmat->matrix.col[h] = col[h];
        outmat->matrix.val[h] = 1;
    }

    col.clear ( );
    row.clear ( );
    index.clear ( );

    CPUsimple_COO_lMatrix<ValueType> matrix_coo;
    matrix_coo.ConvertFrom ( *outmat );
    outmat->ConvertFrom ( matrix_coo );

    //  outmat->WriteFile("output.mtx");
    return outmat;

}

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_CSR_lMatrix<ValueType>::MatrixSupSPower ( const int p ) const
{

    assert ( p > 0 );

    lMatrix<ValueType> *outmat = new CPUsimple_CSR_lMatrix<ValueType>;
    outmat->CloneFrom ( *this );

    lMatrix<ValueType> *tmp;

    for ( int pp = 0; pp < p - 1; ++pp )
    {
        tmp = outmat->MatrixMultSupStructure ( *this );

        outmat->Clear ( );
        delete outmat;

        outmat = tmp;
        tmp = NULL;
    }

    return outmat;

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::MatrixAdd ( const hiflow::la::lMatrix<ValueType> &inmat )
{

    const CPU_CSR_lMatrix<ValueType> *casted_inmat = dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &inmat );

    if ( casted_inmat == NULL )
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::MatrixAdd unsupported inmatrix" );
        this->print ( );
        inmat.print ( );
        exit ( -1 );
    }

    this->MatrixAdd ( *casted_inmat );

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::MatrixAdd ( const CPU_CSR_lMatrix<ValueType> &inmat )
{

    assert ( inmat.get_num_row ( ) == this->get_num_col ( ) );
    assert ( inmat.get_num_col ( ) == this->get_num_row ( ) );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            for ( int k = inmat.matrix.row[i]; k < inmat.matrix.row[i + 1]; ++k )
                if ( this->matrix.col[j] == inmat.matrix.col[k] )
                    this->matrix.val[j] += inmat.matrix.val[k];

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::GershgorinSpectrum ( ValueType *lambda_min, ValueType *lambda_max ) const
{

    assert ( this->get_num_row ( ) == this->get_num_col ( ) );

    ValueType offdiag = 0.0;
    ValueType diag = 0.0;

    *lambda_min = 0.0;
    *lambda_max = 0.0;

    for ( int i = 0; i<this->get_num_row ( ); ++i )
    {
        offdiag = 0.0;

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] != i )
            {
                offdiag += std::abs ( this->matrix.val[j] );
            }
            else
            {
                diag = this->matrix.val[j];
            }

        *lambda_min = ( ( *lambda_min )<( diag - offdiag ) ) ? ( *lambda_min ) : ( diag - offdiag );
        *lambda_max = ( ( *lambda_max )>( diag + offdiag ) ) ? ( *lambda_max ) : ( diag + offdiag );

    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ReadFile ( const char * filename )
{
    this->Clear ( );

    CPUsimple_COO_lMatrix<ValueType> matrix_coo;

    matrix_coo.ReadFile ( filename );

    this->Init ( matrix_coo.get_nnz ( ),
                 matrix_coo.get_num_row ( ),
                 matrix_coo.get_num_col ( ),
                 matrix_coo.get_name ( ) );
    this->ConvertFrom ( matrix_coo ); // the transformation is only there

    matrix_coo.Clear ( );
    this->issymmetric ( );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::WriteFile ( const char * filename ) const
{
    CPUsimple_COO_lMatrix<ValueType> matrix_coo;

    matrix_coo.Init ( this->get_nnz ( ),
                      this->get_num_row ( ),
                      this->get_num_col ( ),
                      this->get_name ( ) );
    matrix_coo.ConvertFrom ( *this ); // the transformation is only there

    matrix_coo.WriteFile ( filename );

    matrix_coo.Clear ( );
}

template <typename ValueType>
hiflow::la::lMatrix<ValueType> *CPU_CSR_lMatrix<ValueType>::extract_submatrix ( const int start_row, const int start_col,
                                                                                const int end_row, const int end_col ) const
{
    assert ( start_row >= 0 );
    assert ( start_col >= 0 );
    assert ( this->get_nnz ( ) > 0 );
    assert ( end_row <= this->get_num_row ( ) );
    assert ( end_col <= this->get_num_col ( ) );
    assert ( end_row > start_row );
    assert ( end_col > start_col );

    hiflow::la::lMatrix<ValueType> *sub_matrix;
    std::string sub_mat_name;
    sub_mat_name = "sub matrix from ";
    sub_mat_name.append ( this->get_name ( ) );

    int sub_mat_nnz = 0;

    std::vector<int> rows, cols;

    for ( int i = start_row; i < end_row; ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( ( this->matrix.col[j] >= start_col ) &&
                 ( this->matrix.col[j] < end_col ) )
            {
                ++sub_mat_nnz;
                rows.push_back ( i - start_row );
                cols.push_back ( this->matrix.col[j] - start_col );
            }

    sub_matrix = this->CloneWithoutContent ( );
    sub_matrix->Init ( sub_mat_nnz,
                       end_row - start_row,
                       end_col - start_col,
                       sub_mat_name );

    //  sub_matrix = init_matrix<ValueType>(sub_mat_nnz,
    //                                      end_row - start_row,
    //                                      end_col - start_col,
    //                                      sub_mat_name,
    //                                      this->get_platform(),
    //                                      this->get_implementation(),
    //                                      this->get_matrix_format()); // CSR

    if ( sub_mat_nnz > 0 )
    {
        sub_matrix->init_structure ( &( rows.front ( ) ), &( cols.front ( ) ) );

        sub_matrix->Zeros ( );
    }

    int sub_mat_ind = 0;

    // extract the sub matrix
    for ( int i = start_row; i < end_row; ++i )
    {

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {

            if ( ( this->matrix.col[j] >= start_col ) &&
                 ( this->matrix.col[j] < end_col ) )
            {

                sub_matrix->add_value ( rows[sub_mat_ind],
                                        cols[sub_mat_ind],
                                        this->matrix.val[j] );

                ++sub_mat_ind;
            }
        }

    }

    return sub_matrix;

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Reorder ( const int *index )
{

    CPU_COO_lMatrix<ValueType> *mat_coo = new CPUsimple_COO_lMatrix<ValueType>;

    mat_coo->ConvertFrom ( *this );
    //  mat_coo->print();
    mat_coo->Reorder ( index );

    this->Clear ( );
    this->ConvertFrom ( *mat_coo );
    //  this->print();
    delete mat_coo;

}

// Thanks to Felix Riehn

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::transpose_me ( void )
{
    CPU_COO_lMatrix<ValueType> *mat_coo = new CPUsimple_COO_lMatrix<ValueType>;
    mat_coo->ConvertFrom ( *this );

    std::swap ( mat_coo->matrix.row, mat_coo->matrix.col );
    mat_coo->swap_dimensions ( );

    this->ConvertFrom ( *mat_coo );

    delete mat_coo;
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ilu0 ( void )
{
    int nrow = this->get_num_row ( );

    LOG_DEBUG ( 1, "lMatrix ILU(0) number of nodes: " << this->get_num_row ( ) << " number of nonzeros: " << this->get_nnz ( ) );

    //position of diagonal element in nnz sequence
    int diag_k = 0;
    double tmp2 = 0.0;

    for ( int row_i = 1; row_i < nrow; ++row_i )
    {
        if ( row_i % ( ( nrow - 1 - ( nrow - 1 ) % 20 ) / 20 ) == 0 )
        {
            LOG_DEBUG ( 3, "Building ILU0 " << float(100 * float(( row_i + 1 ) / float(nrow - 1 ) ) ) << "%" );
        }

        for ( int nz_ik = this->matrix.row[row_i]; this->matrix.col[nz_ik] < row_i; ++nz_ik )
        {

            for ( int p = this->matrix.row[this->matrix.col[nz_ik]];
                  p<this->matrix.row[this->matrix.col[nz_ik] + 1]; ++p )
                if ( this->matrix.col[p] == this->matrix.col[nz_ik] )
                    diag_k = p;

            if ( this->matrix.col[diag_k] != this->matrix.col[nz_ik] )
            {
                LOG_DEBUG ( 0, "ILU0 facotization - the matrix has vanishing diagonal elements!!!" );
            }

            this->matrix.val[nz_ik] = this->matrix.val[nz_ik] / ( 1.0 * this->matrix.val[diag_k] );

            for ( int nz_ij = nz_ik + 1; nz_ij<this->matrix.row[row_i + 1]; ++nz_ij )
            {

                //ordering independent query for a_kj ;-)
                tmp2 = 0.0;
                for ( int p = this->matrix.row[this->matrix.col[nz_ik]];
                      p<this->matrix.row[this->matrix.col[nz_ik] + 1]; ++p )
                    if ( this->matrix.col[p] == this->matrix.col[nz_ij] )
                        tmp2 = this->matrix.val[p];

                this->matrix.val[nz_ij] =
                        this->matrix.val[nz_ij] - this->matrix.val[nz_ik] * tmp2;

            }

        }
    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ilup ( const int p )
{
    //----------------
    // powering matrix
    //----------------
    // to get the full powered matrix pattern
    // one has to make sure there is no destructive
    // intereference between the elements while powering
    // -> hence the matrix to be powered is filled with ones first

    lMatrix<ValueType> *outmat;

    outmat = this->MatrixSupSPower ( p + 1 );

    outmat->Zeros ( ); // set the whole matrix to zero
    outmat->MatrixAdd ( *this ); // outmat = outmat + inmat

    this->ilup ( *outmat, p );

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ilup ( const lMatrix<ValueType> &mat, const int p )
{

    CPU_CSR_lMatrix<ValueType> *outmat = new CPUsimple_CSR_lMatrix<ValueType>;

    outmat->CloneFrom ( mat );
    outmat->Zeros ( );
    outmat->MatrixAdd ( *this );

    //-------------------
    // factorizing matrix
    //-------------------
    CPU_CSR_lMatrix<ValueType> *lev = new CPUsimple_CSR_lMatrix<ValueType>;
    lev->CloneFrom ( *outmat );

    int nrow = outmat->get_num_row ( );
    int new_nnz = 0;

    LOG_DEBUG ( 1, " ILU(p) size of powered pattern: " << outmat->get_nnz ( ) );

    //set levels
    for ( int i = 0; i < nrow; ++i )
    {
        for ( int nzj = outmat->matrix.row[i]; nzj < outmat->matrix.row[i + 1]; ++nzj )
        {
            lev->matrix.val[nzj] = 9999.0;
            if ( outmat->matrix.val[nzj] * outmat->matrix.val[nzj] > 0.0 )
            {
                lev->matrix.val[nzj] = 0.0;
            }
        }
    }

    int nz_kk = 0;
    double pp = ( double ) p; // mat_power-1.0;
    double lev_kj = 0;
    double val_kj = 0;
    double lev_tmp = 0;

    LOG_DEBUG ( 1, " ILU(p) factorizing order: " << p );

    for ( int i = 1; i < nrow; ++i )
    {

        if ( i % ( ( nrow - 1 - ( nrow - 1 ) % 10 ) / 10 ) == 0 )
        {
            LOG_DEBUG ( 3, "ILU(" << p << ") " << float(100 * float(( i + 1 ) / float(nrow - 1 ) ) ) << "%" );
        }

        for ( int nz_ik = outmat->matrix.row[i]; outmat->matrix.col[nz_ik] < i; ++nz_ik )
        {

            if ( lev->matrix.val[nz_ik] <= pp )
            {
                //        LOG_INFO("ilup", " k=" << outmat->matrix.col[nz_ik] );

                //find diagonal element in line k
                for ( int nzh = outmat->matrix.row[outmat->matrix.col[nz_ik]];
                      nzh < outmat->matrix.row[outmat->matrix.col[nz_ik] + 1]; ++nzh )
                    if ( outmat->matrix.col[nzh] == outmat->matrix.col[nz_ik] )
                    {
                        nz_kk = nzh;
                        break;
                    }
                //norm elements
                outmat->matrix.val[nz_ik] = outmat->matrix.val[nz_ik] / ( 1.0 * outmat->matrix.val[nz_kk] );

                // LOG_INFO("ilup", "lower nnz=" << nz_ik+1 << " upper=" << outmat->matrix.row[i+1] );
                for ( int nz_ij = nz_ik + 1; nz_ij < outmat->matrix.row[i + 1]; ++nz_ij )
                {
                    //find element kj if it doesn't exist it is zero
                    val_kj = 0.0;
                    lev_kj = 9999;
                    for ( int nzh = outmat->matrix.row[outmat->matrix.col[nz_ik]];
                          nzh < outmat->matrix.row[outmat->matrix.col[nz_ik] + 1]; ++nzh )
                    {
                        if ( outmat->matrix.col[nzh] == outmat->matrix.col[nz_ij] )
                        {
                            val_kj = outmat->matrix.val[nzh];
                            lev_kj = lev->matrix.val[nzh];
                            break;
                        }
                    }
                    lev_tmp = 1 + lev_kj + lev->matrix.val[nz_ik];

                    outmat->matrix.val[nz_ij] = outmat->matrix.val[nz_ij] - val_kj * outmat->matrix.val[nz_ik];
                    if ( lev->matrix.val[nz_ij] > lev_tmp )
                        lev->matrix.val[nz_ij] = lev_tmp;
                }

            }

        }
        for ( int nzh = outmat->matrix.row[i]; nzh < outmat->matrix.row[i + 1]; ++nzh )
            if ( lev->matrix.val[nzh] > pp )
            {
                outmat->matrix.val[nzh] = 0.0;
                lev->matrix.val[nzh] = 9999.0;
            }
    }

    int testcount = 0;
    for ( int i = 0; i < nrow; ++i )
    {
        for ( int nzj = outmat->matrix.row[i]; nzj < outmat->matrix.row[i + 1]; ++nzj )
        {
            if ( lev->matrix.val[nzj] <= pp )
                testcount++;
        }
    }

    LOG_DEBUG ( 1, " ILU(p) pattern size after factorization: " << testcount );
    new_nnz = testcount;

    // compressing the matrix by levels

    // init the memory for the matrix
    this->Clear ( ); //
    this->Init ( new_nnz, nrow, nrow, "ILU" );

    char data_info[255];
    sprintf ( data_info, "(%d)", p );
    this->name_.append ( data_info );

    int nzcount = -1;
    int nz_row = 0;
    for ( int i = 0; i < nrow; ++i )
    {
        nz_row = -1;
        for ( int nzj = outmat->matrix.row[i]; nzj < outmat->matrix.row[i + 1]; ++nzj )
        {
            if ( lev->matrix.val[nzj] <= pp )
            {
                nzcount++;
                nz_row++;
                //this->add_value(i,outmat->matrix.col[nzj],outmat->matrix.val[nzj]);
                if ( nz_row == 0 )
                    this->matrix.row[i] = nzcount;

                this->matrix.val[nzcount] = outmat->matrix.val[nzj];
                this->matrix.col[nzcount] = outmat->matrix.col[nzj];

            }
        }
    }

    this->matrix.row[this->get_num_row ( )] = new_nnz;

    delete outmat;

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::ilu_solve ( const hiflow::la::lVector<ValueType> &invec,
                                             hiflow::la::lVector<ValueType> *outvec ) const
{
    const CPU_lVector<ValueType> *casted_invec = dynamic_cast < const CPU_lVector<ValueType>* > ( &invec );
    CPU_lVector<ValueType> *casted_outvec = dynamic_cast < CPU_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::ilu_solve unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    assert ( casted_invec ->get_size ( ) > 0 );
    assert ( casted_outvec->get_size ( ) > 0 );
    assert ( casted_invec ->get_size ( ) == this->get_num_col ( ) );
    assert ( casted_outvec->get_size ( ) == this->get_num_row ( ) );

    // M = (1+L) (D+U)
    // Mz=r
    // forward step

    int last_j = 0;
    for ( int i = 0; i<this->get_num_row ( ); ++i )
    {

        casted_outvec->buffer[i] = casted_invec->buffer[i];

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {
            if ( i > this->matrix.col[j] )
                casted_outvec->buffer[i] -= casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
            if ( i == this->matrix.col[j] )
            {
                last_j = j;
                break;
            }
        }

    }

    // backward
    // (D+U)z = y
    for ( int i = this->get_num_row ( ) - 1; i >= 0; --i )
    {

        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
        {

            if ( i == this->matrix.col[j] )
                last_j = j;
            if ( i < this->matrix.col[j] )
                casted_outvec->buffer[i] -= casted_outvec->buffer[ this->matrix.col[j] ] * this->matrix.val[j];
        }

        assert ( this->matrix.val[last_j] != 0.0 );
        casted_outvec->buffer[i] /= this->matrix.val[last_j];

    }

}

// Thanks to Felix Riehn

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Multicoloring ( int &ncolors, int **color_sizes, int **permut_index ) const
{
    //variables to store row's first nnz and dependencies of the nodes
    int m, m1, dep;

    //highest colour aka number of colours
    int hicol = 0;

    //node colour array (square matrix)
    int *nc = new int[this->get_num_row ( )];
    for ( int i = 0; i<this->get_num_row ( ); ++i )
        nc[i] = 0;

    //loop over rows of matrix/nodes in grid
    for ( int j = 0; j<this->get_num_row ( ); ++j )
    {
        //set colour of current node to 1
        nc[j]++;

        m = this->matrix.row[j];
        m1 = this->matrix.row[j + 1];

        LOG_DEBUG ( 9, "nnz intervall: " << m << ".." << m1 - 1 );
        LOG_DEBUG ( 9, " dependences of node " << j << " : colour " );

        //loop over (nnz's/nodes) (that are in row j/that node j is linked to)
        for ( int k = m; k < m1; ++k )
            for ( int i = m; i < m1; ++i )
            {

                //current node/column that node/element j depends of
                dep = this->matrix.col[i];

                //exclude diagonal elements
                if ( dep != j )
                {

                    //check whether any linked node has same colour as node j if so change colour
                    if ( nc[dep] == nc[j] )
                    {
                        ++nc[j];
                    }

                    if ( nc[j] > hicol )
                        hicol = nc[j];

                }
            }

    }

    //adapt to interface
    ncolors = hicol;

    //check for colour zero (shouldn't exist)
    for ( int i = 0; i<this->get_num_row ( ); ++i )
        if ( nc[i] == 0 )
        {

            LOG_ERROR ( "multi-coloring - ZERO COLOR! (shouldn't exist)" );
            this->print ( );
            LOG_ERROR ( "multi-coloring - not exiting !!!" );
        }

    // colour multiplicity
    ( *color_sizes ) = new int[ncolors];

    //permutation index
    ( *permut_index ) = new int[this->get_num_row ( )];

    //accounting array for the ordering by colour
    int* colcount = new int[ncolors];

    for ( int k = 0; k < ncolors; ++k )
    {
        ( *color_sizes )[k] = 0;
        colcount[k] = 0;
    }

    //get colour multiplicity (here there is a shift in colour notation by -1!!)
    //which should be changed cause it's not good style
    for ( int i = 0; i<this->get_num_row ( ); ++i )
    {
        ++( *color_sizes )[nc[i] - 1];
    }

    //check for consistency between number of nodes == ncolours*multiplicity

    int tmp = 0;
    for ( int k = 0; k < ncolors; ++k )
        tmp += ( *color_sizes )[k];

    if ( tmp != this->get_num_row ( ) )
    {
        LOG_ERROR ( "multi-coloring - check for consistency between number of nodes - FAILED !!!" );
        this->print ( );
        LOG_ERROR ( "multi-coloring - not exiting !!!" );

    }

    //reorder nodes by color
    int position;
    for ( int j = 0; j<this->get_num_row ( ); ++j )
    {
        position = 0;

        for ( int k = 0; k < nc[j] - 1; ++k )
            position += ( *color_sizes )[k];

        position += colcount[nc[j] - 1];
        ( *permut_index )[position] = j;

        ++colcount[nc[j] - 1];
    }

}

// Thanks to Niels Wegh

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::Levelscheduling ( int &nlevels, int **level_sizes, int **permut_index ) const
{
    // src code
    // Niels Wegh

    int nrow = this->get_num_row ( );
    ( *permut_index ) = new int[this->get_num_row ( )];

    std::vector<int> Level;
    std::vector<int> Knoten ( nrow );
    std::vector<int> newKnoten ( nrow );
    std::vector<int>::iterator iter ( Knoten.begin ( ) );
    std::vector<int> depth ( nrow );
    std::vector<bool> besuch ( nrow );

    //presettings
    for ( int i = 0; iter != Knoten.end ( ); ++i, ++iter )
    {
        *iter = i;
        besuch[i] = false;
    }

    //newKnoten=Knoten;
    int done = 0;
    bool Stop = false;
    int nlev = 0;
    //int k = 0;

    // depth
    for ( int L = 0; L < nrow; ++L )
    {
        if ( done == nrow )
        {
            nlev = Level.size ( );
            break;
        }
        for ( int k = 0; k < nrow; ++k )
        {
            if ( besuch[k] == true )
                continue;
            if ( depth[k] == L )
            {

                // this loops scans the L part if there is any node that has to be calculated before this node
                for ( int j = this->matrix.row[k]; j<this->matrix.row[k + 1]; ++j )
                {
                    if ( this->matrix.col[j] < k )
                    {
                        if ( besuch[this->matrix.col[j]] == false )
                        {
                            Stop = true;
                            break;
                        }
                    }
                }
                if ( Stop == true )
                {
                    Stop = false;
                    continue;
                }
                ( *permut_index )[done] = k;
                done += 1;
                besuch[k] = true;

                for ( int j = this->matrix.row[k]; j<this->matrix.row[k + 1]; ++j )
                { // makes the depth for the nodes that are behind this node.
                    if ( this->matrix.col[j] > k )
                    {
                        if ( besuch[this->matrix.col[j]] == true )
                        {
                            LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::Levelscheduling - internal builder error" )
                        }
                        if ( depth[this->matrix.col[j]] <= depth[k] )
                        {
                            depth[this->matrix.col[j]] = depth[k] + 1;
                        }
                    }
                }
            }
        }
        Level.push_back ( done );
        LOG_DEBUG ( 5, "CPU_CSR_lMatrix<ValueType>::Levelscheduling - Anzahl der Knoten im Level(" << L << ") :" << Level[L] );
    }

    LOG_DEBUG ( 2, "CPU_CSR_lMatrix<ValueType>::Levelscheduling  - Anzahl der Level: " << nlev << " und done: " << done );

    nlevels = nlev;
    ( *level_sizes ) = new int[nlevels];

    int sum = 0;
    for ( int i = 0; i < nlevels; ++i )
    {
        ( *level_sizes )[i] = Level[i] - sum;
        LOG_DEBUG ( 3, "CPU_CSR_lMatrix<ValueType>::Levelscheduling - L=" << Level[i] - sum );
        sum = Level[i];
    }

}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::extract_diagelements ( const int start_i, const int end_i, lVector<ValueType> *vec ) const
{

    ValueType *tmp_buff = new ValueType[end_i - start_i];
    assert ( tmp_buff != NULL );

    // extract the diag elements
    for ( int i = start_i; i < end_i; ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] == i )
                tmp_buff[i - start_i] = this->matrix.val[j];

    vec->SetBlockValues ( 0, end_i - start_i,
                          tmp_buff );

    delete [] tmp_buff;
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::extract_invdiagelements ( const int start_i, const int end_i, lVector<ValueType> *vec ) const
{

    ValueType *tmp_buff = new ValueType[end_i - start_i];
    assert ( tmp_buff != NULL );

    // extract the diag elements
    for ( int i = start_i; i < end_i; ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] == i )
                tmp_buff[i - start_i] = 1 / this->matrix.val[j];

    vec->SetBlockValues ( 0, end_i - start_i,
                          tmp_buff );

    delete [] tmp_buff;

}

// Thanks to Felix Riehn

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::compress_me ( )
{
    if ( this->get_nnz ( ) > 0 )
    {
        int nrow = this->get_num_row ( );

        int new_nnz = 0;
        for ( int i = 0; i < nrow; ++i )
        {
            for ( int nzj = this->matrix.row[i]; nzj<this->matrix.row[i + 1]; ++nzj )
            {
                if ( this->matrix.val[nzj] != 0.0 ) ++new_nnz;
            }
        }

        LOG_DEBUG ( 1, "lMatrix->compress_me() pattern size after compressing: "
                    << new_nnz << " before:" << this->get_nnz ( )
                    << " difference:" << this->get_nnz ( ) - new_nnz );

        if ( new_nnz < this->get_nnz ( ) )
        {
            CPU_CSR_lMatrix<ValueType> *mat = new CPUsimple_CSR_lMatrix<ValueType>;
            mat->CloneFrom ( *this );

            char data_info[255];
            sprintf ( data_info, "compressed with %d", this->get_nnz ( ) - new_nnz );

            // compressing

            // init the memory for the matrix
            this->Clear ( ); //
            this->Init ( new_nnz, mat->get_num_row ( ), mat->get_num_col ( ), mat->get_name ( ) );
            this->name_.append ( data_info );

            int nzcount = -1;
            int nz_row = 0;
            for ( int i = 0; i < nrow; ++i )
            {
                nz_row = -1;
                for ( int nzj = mat->matrix.row[i]; nzj < mat->matrix.row[i + 1]; ++nzj )
                {
                    if ( mat->matrix.val[nzj] != 0.0 )
                    {
                        ++nzcount;
                        ++nz_row;

                        if ( nz_row == 0 ) this->matrix.row[i] = nzcount;

                        this->matrix.val[nzcount] = mat->matrix.val[nzj];
                        this->matrix.col[nzcount] = mat->matrix.col[nzj];
                    }
                }
            }

            assert ( new_nnz == nzcount + 1 );

            this->matrix.row[this->get_num_row ( )] = new_nnz;

            delete mat;
        }
    }
}

template <typename ValueType>
bool CPU_CSR_lMatrix<ValueType>::issymmetric ( )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    bool half1 = false;
    bool half2 = false;

    CPUsimple_CSR_lMatrix<ValueType> tmp;
    tmp.CloneFrom ( *this );
    tmp.transpose_me ( );

    tmp.Scale ( -1.0 );
    tmp.MatrixAdd ( *this );

    tmp.compress_me ( );
    if ( tmp.get_nnz ( ) == 0 )
        half1 = true;

    tmp.Clear ( );
    tmp.CloneFrom ( *this );
    CPUsimple_CSR_lMatrix<ValueType> tmp2;
    tmp2.CloneFrom ( *this );
    tmp2.transpose_me ( );

    tmp.Scale ( -1.0 );
    tmp.MatrixAdd ( tmp2 );
    tmp.compress_me ( );

    if ( tmp.get_nnz ( ) == 0 )
        half2 = true;

    if ( ( half1 == true ) || ( half2 == true ) )
    {

        this->symmetric_ = true;
        return true;

    }
    else
    {

        this->symmetric_ = false;
        return false;

    }
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::delete_diagonal ( )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] == i )
                this->matrix.val[j] = 0.0;

    // keep the diagonal entries
    // this->compress_me();
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::delete_offdiagonal ( )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] != i )
                this->matrix.val[j] = 0.0;

    this->compress_me ( );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::delete_lower_triangular ( )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] > i )
                this->matrix.val[j] = 0.0;

    this->compress_me ( );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::delete_strictly_lower_triangular ( )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    this->delete_lower_triangular ( );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] == i )
                this->matrix.val[j] = 0.0;

    // the first row is missing - do not compress!!!
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::delete_upper_triangular ( )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] < i )
                this->matrix.val[j] = 0.0;

    this->compress_me ( );
}

template <typename ValueType>
void CPU_CSR_lMatrix<ValueType>::delete_strictly_upper_triangular ( )
{
    assert ( this->get_nnz ( ) > 0 );
    assert ( this->get_num_row ( ) > 0 );
    assert ( this->get_num_col ( ) > 0 );

    this->delete_upper_triangular ( );

    for ( int i = 0; i<this->get_num_row ( ); ++i )
        for ( int j = this->matrix.row[i]; j<this->matrix.row[i + 1]; ++j )
            if ( this->matrix.col[j] == i )
                this->matrix.val[j] = 0.0;

    // the last row is missing - do not compress!!!
}

template class CPU_CSR_lMatrix<float>;
template class CPU_CSR_lMatrix<double>;
